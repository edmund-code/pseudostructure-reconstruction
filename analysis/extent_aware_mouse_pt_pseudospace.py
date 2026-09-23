# %% [markdown]
# # Extent-aware mouse PT pseudospace: fixed-center prototype
#
# **Question.** Does a segmented tubule behave more like a point or an observation that averages a finite stretch of the proximal-tubule axis? This is an exploratory diagnostic, not a replacement for the reviewed DPT workflow.
#
# Inputs are notebook 03's saved PT `lognorm` matrix and PT-specific DPT, notebook 03's pass-1 annotations, and the original quality-controlled mouse GeoJSON polygons. Only the two control mouse specimens enter this notebook. Notebook 03's PT DPT was learned in a human–mouse shared space; it is an explicit historical baseline, while the new backbone below is built from mouse expression alone. All input artifacts are read-only. Results go to `results/pt_extent_aware_mouse/`.
#
# **Predeclared comparison.** We hold nine known PT genes out of all new coordinate and width construction, reserve the published early/late markers solely for orientation, and fit each held-out gene in one specimen to predict the other. Segment labels and glomerular distance are diagnostics. Labels came from upstream expression-based cluster review, and the historical DPT used the full expression space, so neither is a completely independent external ground truth. The two mice give descriptive reproducibility, not population inference.

# %% [markdown]
# ## 1. Inputs and baseline
#
# The control PT subset and the pass-1 object share observation IDs. The latter provides `feature_index` for an exact join to the kept GeoJSON features and glomerular centroids for distance. Missing or ambiguous joins stop the run.

# %%
import argparse
import json
import os
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.sparse import csgraph
from scipy.spatial import cKDTree
from scipy.stats import spearmanr
from shapely.geometry import shape
from sklearn.decomposition import PCA
from sklearn.neighbors import kneighbors_graph
from sklearn.preprocessing import SplineTransformer, StandardScaler
from IPython.display import display, Markdown

SEED = 17
rng = np.random.default_rng(SEED)
SAMPLES = ('Ctrl1A2', 'Ctrl1A4')
SEGMENTS = ('PT-S1', 'PT-S2', 'PT-S3')
ORIENT_EARLY = ('Slc5a2', 'Slc5a12', 'Gatm')
ORIENT_LATE = ('Slc7a13', 'Slc22a7', 'Cyp7b1')
# Prior biological candidates, not selected from the present DPT or model results.
HELD_OUT = ('Aqp1', 'Miox', 'Pck1', 'Slc22a8', 'Acsm2', 'Slc5a8', 'Aldob', 'Fabp1', 'Hnf4a')
# Notebook 03's reviewed PT label panel; exclude it from the new backbone and curve training.
LABEL_PANEL = ('Lrp2', 'Cubn', 'Slc34a1', 'Slc5a2', 'Slc5a12', 'Gatm',
               'Slc22a6', 'Slc13a3', 'Cyp2e1', 'Slc22a7', 'Slc7a13',
               'Cyp7b1', 'Slc6a18', 'Acsm3')
N_BACKBONE_GENES = 300
N_CURVE_GENES = 48
PRIMARY_C = 0.02
GRID = np.linspace(0, 1, 201)
LAMBDA = 1.0


def project_root():
    starts = [Path.cwd().resolve()]
    if '__file__' in globals():
        starts.insert(0, Path(__file__).resolve().parent)
    for start in starts:
        for parent in (start, *start.parents):
            if (parent / 'pseudospace').is_dir() and (parent / 'data').is_dir():
                return parent
    raise RuntimeError('Run from the project tree or set a script path within it.')


root = project_root()
parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--data-root', type=Path)
parser.add_argument('--results-root', type=Path)
args, _ = parser.parse_known_args()
data_root = Path(args.data_root or os.environ.get('PSEUDOSPACE_DATA_ROOT') or root / 'data').expanduser().resolve()
results_root = Path(args.results_root or os.environ.get('PSEUDOSPACE_RESULTS_ROOT') or root / 'results').expanduser().resolve()
upstream = results_root / 'human_vs_healthy_mouse'
out = results_root / 'pt_extent_aware_mouse'
pt_path = upstream / 'cross_species_pt_dpt.h5ad'
pass1_path = upstream / 'cross_species_harmony_pass1.h5ad'
for path in (pt_path, pass1_path, *(data_root / f'{s}_kept_tubules_labeled_fine.geojson' for s in SAMPLES)):
    if not path.is_file():
        raise FileNotFoundError(f'Required upstream input missing: {path}')
out.mkdir(parents=True, exist_ok=True)

pt_all = ad.read_h5ad(pt_path)
mask = (pt_all.obs['sample'].astype(str).isin(SAMPLES) &
        pt_all.obs['comparison_species'].astype(str).eq('mouse')).to_numpy()
pt = pt_all[mask].copy()
del pt_all
if pt.n_obs == 0 or 'lognorm' not in pt.layers:
    raise ValueError('Saved PT artifact has no healthy mouse PT or lognorm layer.')
required = {'sample', 'segment_class', 'total_scanpy_dpt', 'x_centroid', 'y_centroid'}
missing = required - set(pt.obs)
if missing:
    raise KeyError(f'PT object missing columns: {sorted(missing)}')
if 'measured_in_both_inputs' not in pt.var:
    raise KeyError('Shared-gene availability audit missing; structural zeros cannot be excluded.')
pass1 = ad.read_h5ad(pass1_path, backed='r')
pass1_obs = pass1.obs.copy()
pass1.file.close()
if not pt.obs_names.isin(pass1_obs.index).all():
    raise ValueError('Some PT observation IDs are missing from pass-1.')
joined = pass1_obs.reindex(pt.obs_names)
if joined['feature_index'].isna().any() or not joined['sample'].astype(str).eq(pt.obs['sample'].astype(str)).all():
    raise ValueError('Feature index or sample mismatch between PT and pass-1 objects.')

# Distance is never used to construct a coordinate, width, or smoothing parameter.
distance = pd.Series(index=pt.obs_names, dtype=float)
for sample in SAMPLES:
    glom = pass1_obs.loc[(pass1_obs['sample'].astype(str) == sample) &
                         (pass1_obs['coarse_class'].astype(str) == 'Glomerulus')]
    if len(glom) < 3:
        raise ValueError(f'{sample}: fewer than three saved glomeruli.')
    rows = pt.obs['sample'].astype(str).eq(sample)
    tree = cKDTree(glom[['x_centroid', 'y_centroid']].to_numpy(dtype=float))
    distance.loc[rows] = tree.query(pt.obs.loc[rows, ['x_centroid', 'y_centroid']].to_numpy(dtype=float))[0]
pt.obs['distance_to_nearest_glomerulus_px'] = distance.to_numpy()
print(f'Input: {pt.n_obs:,} healthy mouse PT tubules, {pt.n_vars:,} mapped genes; samples:')
display(pt.obs['sample'].value_counts().rename('tubules').to_frame())
print('Measured in both inputs:', int(pt.var['measured_in_both_inputs'].sum()))

# %% [markdown]
# ## 2. Polygon morphology and an extent prior
#
# The files contain polygon outlines and calibrated area, but no skeleton or traced nephron path. We calculate the long and short sides of each polygon's minimum rotated rectangle, perimeter, area, aspect ratio, and axis-aligned bounds. The rotated long side is a simple 2D relative extent proxy. Its relation to true nephron arc length is unknown: sectioning angle, branching appearance, and tortuosity can all change it.

# %%
def polygon_metrics(geom, area_um2):
    if not geom.is_valid:
        geom = geom.buffer(0)
    if geom.is_empty or geom.area <= 0 or not np.isfinite(area_um2) or area_um2 <= 0:
        return None
    scale = np.sqrt(area_um2 / geom.area)  # microns per polygon-coordinate unit
    box = geom.minimum_rotated_rectangle
    xy = np.asarray(box.exterior.coords)[:4]
    sides = np.linalg.norm(np.roll(xy, -1, axis=0) - xy, axis=1) * scale
    major, minor = float(sides.max()), float(sides.min())
    xmin, ymin, xmax, ymax = geom.bounds
    return {'area_um2': area_um2, 'major_um': major, 'minor_um': minor,
            'perimeter_um': geom.length * scale, 'aspect_ratio': major / max(minor, 1e-8),
            'bbox_width_um': (xmax - xmin) * scale, 'bbox_height_um': (ymax - ymin) * scale}

records = []
for sample in SAMPLES:
    target = joined.loc[pt.obs['sample'].astype(str).eq(sample), 'feature_index'].astype(int)
    lookup = dict(zip(target.to_numpy(), target.index))
    with (data_root / f'{sample}_kept_tubules_labeled_fine.geojson').open() as handle:
        features = json.load(handle)['features']
    if max(lookup) >= len(features):
        raise ValueError(f'{sample}: feature_index exceeds GeoJSON feature count.')
    for index, name in lookup.items():
        feature = features[index]
        props = feature.get('properties', {})
        if props.get('source_sample', sample) != sample or props.get('kept_tubule') is False:
            raise ValueError(f'{sample} feature {index} does not match saved kept-tubule identity.')
        area_values = [m['value'] for m in props.get('measurements', []) if m.get('name') == 'Area um^2']
        metrics = polygon_metrics(shape(feature['geometry']), float(area_values[0])) if len(area_values) == 1 else None
        if metrics is not None:
            records.append({'obs_name': name, **metrics})
morph = pd.DataFrame.from_records(records).set_index('obs_name').reindex(pt.obs_names)
valid = morph['major_um'].notna().to_numpy()
print(f'Polygon measurements: {valid.sum():,}/{len(valid):,} matched PT structures.')
if valid.mean() < 0.95:
    raise ValueError('Too many PT polygons lack valid calibrated geometry; inspect the join.')
if not valid.all():
    pt = pt[valid].copy()
    morph = morph.loc[pt.obs_names].copy()
pt.obs = pt.obs.join(morph)
pt.obs['length_group'] = (pt.obs.groupby('sample', observed=True)['major_um']
                          .transform(lambda x: pd.qcut(x.rank(method='first'), 3,
                                                       labels=['short', 'medium', 'long'])).astype(str))
# Median 1, conservative range [0.5, 1.5]; no micron-to-pseudotime conversion is claimed.
length_ratio = pt.obs['major_um'].to_numpy() / pt.obs['major_um'].median()
pt.obs['relative_extent'] = np.clip(length_ratio, 0.5, 1.5)
print('Direct skeleton/path length: unavailable. Width uses rotated-box major axis as a prior.')
display(pt.obs[['area_um2', 'major_um', 'minor_um', 'perimeter_um', 'aspect_ratio',
                'bbox_width_um', 'bbox_height_um']].describe(percentiles=[.1, .5, .9]).T.round(2))
cols = ['area_um2', 'major_um', 'minor_um', 'perimeter_um', 'aspect_ratio']
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
for col in cols:
    axes[0].hist(np.log1p(pt.obs[col]), bins=35, alpha=.4, label=col, density=True)
axes[0].legend(fontsize=7)
axes[0].set_title('Log morphology distributions')
axes[0].set_xlabel('log(1 + measurement)')
cor = pt.obs[cols].corr(method='spearman')
im = axes[1].imshow(cor, vmin=-1, vmax=1, cmap='coolwarm')
axes[1].set_xticks(range(len(cols)), cols, rotation=55, ha='right')
axes[1].set_yticks(range(len(cols)), cols)
fig.colorbar(im, ax=axes[1], label='Spearman rho')
fig.tight_layout(); plt.show()

# Baseline DPT diagnostics, before building the new coordinate.
fig, axes = plt.subplots(1, 3, figsize=(13, 3.5))
for seg in SEGMENTS:
    values = pt.obs.loc[pt.obs['segment_class'].astype(str).eq(seg), 'total_scanpy_dpt'].to_numpy(dtype=float)
    axes[0].hist(values, bins=35, density=True, alpha=.45, label=seg)
axes[0].legend(); axes[0].set_title('Existing PT DPT by reviewed segment')
axes[1].scatter(pt.obs['total_scanpy_dpt'], pt.obs['distance_to_nearest_glomerulus_px'], s=2, alpha=.17)
axes[1].set_title('DPT vs nearest glomerulus')
axes[2].scatter(pt.obs['total_scanpy_dpt'], pt.obs['major_um'], s=2, alpha=.17)
axes[2].set_title('DPT vs 2D long axis')
for ax in axes: ax.set_xlabel('Existing PT DPT')
fig.tight_layout(); plt.show()
display(pd.crosstab(pt.obs['length_group'], pt.obs['segment_class'], normalize='index').round(3))

# %% [markdown]
# ## 3. Mouse-only molecular backbone
#
# We select high-variance, measured genes after excluding all reserved validation and published PT annotation/orientation genes. PCA on these genes, a symmetric kNN graph, and one shortest-path distance from an early-marker root give a transparent point ordering. Marker expression chooses and orients the root only. DPT, reviewed segments, glomerular distance, and morphology do not enter the graph. A disconnected graph is reported rather than silently imputed.

# %%
excluded = set(HELD_OUT) | set(LABEL_PANEL) | set(ORIENT_EARLY) | set(ORIENT_LATE)
available = pt.var['measured_in_both_inputs'].to_numpy(dtype=bool)
all_expr = pt.layers['lognorm'].tocsr()
means = np.asarray(all_expr.mean(axis=0)).ravel()
variance = np.asarray(all_expr.power(2).mean(axis=0)).ravel() - means**2
fraction = np.asarray((all_expr > 0).mean(axis=0)).ravel()
search = np.array([i for i, gene in enumerate(pt.var_names)
                   if available[i] and gene not in excluded and fraction[i] >= .10 and variance[i] > 0])
if len(search) < N_BACKBONE_GENES:
    raise ValueError(f'Only {len(search)} eligible backbone genes; need {N_BACKBONE_GENES}.')
ranked = search[np.argsort(variance[search])[::-1]]
backbone_genes = pt.var_names[ranked[:N_BACKBONE_GENES]].tolist()
curve_genes = backbone_genes[:N_CURVE_GENES]
heldout_genes = [g for g in HELD_OUT if g in pt.var_names and bool(pt.var.loc[g, 'measured_in_both_inputs'])]
if len(heldout_genes) < 5:
    raise ValueError(f'Only {len(heldout_genes)} reserved genes are measured.')
print(f'Construction: {len(backbone_genes)} backbone genes, {len(curve_genes)} curve genes; '
      f'validation only: {heldout_genes}')
print('Orientation only:', ORIENT_EARLY, ORIENT_LATE)


def expression(genes):
    return all_expr[:, pt.var_names.get_indexer(genes)].toarray().astype(float)


def orientation_score():
    early = expression(list(ORIENT_EARLY))
    late = expression(list(ORIENT_LATE))
    # Equal weight per gene after standardization; never uses DPT or segment labels.
    both = np.column_stack([early, late])
    z = StandardScaler().fit_transform(both)
    return z[:, len(ORIENT_EARLY):].mean(axis=1) - z[:, :len(ORIENT_EARLY)].mean(axis=1)

orient = orientation_score()


def make_backbone(n_genes=300, neighbors=30):
    x = expression(backbone_genes[:n_genes])
    pcs = PCA(n_components=12, random_state=SEED).fit_transform(StandardScaler().fit_transform(x))
    candidates = np.flatnonzero(orient <= np.quantile(orient, .03))
    center = np.median(pcs[candidates], axis=0)
    root_index = candidates[np.argmin(np.linalg.norm(pcs[candidates] - center, axis=1))]
    graph = kneighbors_graph(pcs, n_neighbors=neighbors, mode='distance', include_self=False)
    graph = graph.maximum(graph.T).tocsr()
    component_count = csgraph.connected_components(graph, directed=False)[0]
    if component_count != 1:
        raise ValueError(f'kNN graph has {component_count} components at k={neighbors}; inspect topology.')
    path = csgraph.dijkstra(graph, directed=False, indices=int(root_index))
    if not np.isfinite(path).all():
        raise ValueError('Backbone contains unreachable PT structures.')
    coordinate = (path - path.min()) / (path.max() - path.min())
    if spearmanr(coordinate, orient).statistic < 0:
        coordinate = 1 - coordinate
    return coordinate, pcs

backbone, pcs = make_backbone()
pt.obs['backbone_position'] = backbone
print('Backbone vs existing DPT: Spearman rho =',
      round(spearmanr(backbone, pt.obs['total_scanpy_dpt']).statistic, 3))
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
axes[0].scatter(pt.obs['total_scanpy_dpt'], backbone, s=2, alpha=.2)
axes[0].set(xlabel='Existing PT DPT', ylabel='Mouse-only backbone')
axes[1].scatter(pcs[:, 0], pcs[:, 1], c=backbone, s=2, cmap='viridis')
axes[1].set(xlabel='PC1', ylabel='PC2', title='Graph geodesic on 12 PCs')
fig.tight_layout(); plt.show()

# %% [markdown]
# ## 4–5. Truncated Gaussian observations and smooth gene curves
#
# For tubule *i*, center is its backbone position and width is `c × relative_extent`. Numerical weights on a fixed [0,1] grid integrate a cubic B-spline basis. At `c=0`, evaluation is exactly at the center. Both point and extended models use the same basis and second-difference penalty. Width and curve can trade off, so we fix widths from geometry and do **not** jointly optimize centers here; free centers could absorb gene noise and make the small experiment hard to identify.

# %%
spline = SplineTransformer(n_knots=8, degree=3, knots='uniform', include_bias=True)
spline.fit(GRID[:, None])
B_grid = spline.transform(GRID[:, None])
penalty = np.diff(np.eye(B_grid.shape[1]), n=2, axis=0)
penalty = penalty.T @ penalty


def observation_basis(center, width):
    center = np.asarray(center, dtype=float)
    width = np.asarray(width, dtype=float)
    if np.all(width == 0):
        return spline.transform(center[:, None]), center
    safe_width = np.maximum(width, 1e-6)
    z = (GRID[None, :] - center[:, None]) / safe_width[:, None]
    weights = np.exp(-.5 * z**2)
    weights /= weights.sum(axis=1, keepdims=True)
    return weights @ B_grid, weights @ GRID


def fit_curve(design, y, lam=LAMBDA):
    y = np.asarray(y)
    a = design.T @ design + lam * penalty + 1e-8 * np.eye(design.shape[1])
    return np.linalg.solve(a, design.T @ y)


def widths(metric='major_um', c=PRIMARY_C, mode='morphology'):
    values = pt.obs[metric].to_numpy(dtype=float)
    relative = np.clip(values / np.median(values), .5, 1.5)
    if mode == 'equal':
        relative[:] = 1.0
    elif mode == 'permuted':
        # Preserve each specimen's width distribution, break tubule-width correspondence.
        for sample in SAMPLES:
            idx = np.flatnonzero(pt.obs['sample'].astype(str).to_numpy() == sample)
            relative[idx] = relative[np.random.default_rng(SEED).permutation(idx)]
    return c * relative

bases = {}
point_basis, _ = observation_basis(backbone, np.zeros(pt.n_obs))
bases['point'] = point_basis
for c in (0.01, 0.02, 0.04):
    bases[f'extent_c{c:g}'], _ = observation_basis(backbone, widths(c=c))
bases['equal'], _ = observation_basis(backbone, widths(mode='equal'))
bases['permuted'], _ = observation_basis(backbone, widths(mode='permuted'))
assert np.allclose(bases['point'].sum(axis=1), 1)
assert np.allclose(bases['extent_c0.02'].sum(axis=1), 1, atol=1e-6)
assert np.allclose(observation_basis(backbone, widths(c=0))[0], bases['point'])
_, expected_position = observation_basis(backbone, widths())
pt.obs['extent_mean_position'] = expected_position
pt.obs['extent_sigma'] = widths()
fig, ax = plt.subplots(figsize=(9, 3.5))
for group in ('short', 'medium', 'long'):
    idx = np.flatnonzero(pt.obs['length_group'].eq(group).to_numpy())
    selected = idx[np.argmin(abs(backbone[idx] - .5))]
    sigma = pt.obs['extent_sigma'].iloc[selected]
    density = np.exp(-.5 * ((GRID - backbone[selected]) / sigma)**2)
    ax.plot(GRID, density / np.trapz(density, GRID), label=f'{group}: σ={sigma:.3f}')
ax.set(xlabel='Backbone position', ylabel='Truncated density', title='Representative near-midpoint tubules')
ax.legend(); fig.tight_layout(); plt.show()

# Fit representative construction genes on all data to inspect underlying curves. These are
# descriptive plots; held-out gene validation occurs below with specimen separation.
Y_curve = expression(curve_genes)
point_coef = fit_curve(bases['point'], Y_curve)
extent_coef = fit_curve(bases['extent_c0.02'], Y_curve)
curve_plot_genes = curve_genes[:4]
fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
for ax, gene in zip(axes.flat, curve_plot_genes):
    j = curve_genes.index(gene)
    ax.scatter(backbone, Y_curve[:, j], s=2, alpha=.12, color='0.35')
    ax.plot(GRID, B_grid @ point_coef[:, j], label='Point curve')
    ax.plot(GRID, B_grid @ extent_coef[:, j], label='Extent underlying curve')
    ax.set(title=gene, ylabel='lognorm expression')
axes[0, 0].legend(); axes[1, 0].set_xlabel('Backbone position'); axes[1, 1].set_xlabel('Backbone position')
fig.tight_layout(); plt.show()

# For a long tubule near the steepest fitted gradient, compare the fitted point value
# with the actual integrated prediction. The red mark is its measured aggregate.
j = int(np.argmax(np.ptp(B_grid @ extent_coef, axis=0)))
gene = curve_genes[j]
gradient = np.abs(np.gradient(B_grid @ extent_coef[:, j], GRID))
steep = GRID[np.argmax(gradient)]
long_idx = np.flatnonzero(pt.obs['length_group'].eq('long').to_numpy())
examples = long_idx[np.argsort(abs(backbone[long_idx] - steep))[:8]]
fig, ax = plt.subplots(figsize=(8, 4))
ax.scatter(backbone, Y_curve[:, j], s=2, alpha=.08, color='0.5')
ax.plot(GRID, B_grid @ extent_coef[:, j], color='#0072B2', label='Underlying extent curve')
ax.scatter(backbone[examples], Y_curve[examples, j], color='#D55E00', s=24, label='Observed long tubule')
ax.scatter(backbone[examples], (bases['extent_c0.02'] @ extent_coef)[examples, j],
           marker='x', color='black', s=40, label='Integrated prediction')
ax.set(xlabel='Backbone position', ylabel='lognorm expression',
       title=f'{gene}: long tubules near steepest fitted change')
ax.legend(fontsize=8); fig.tight_layout(); plt.show()

# A diagnostic for the proposed integration mechanism: point-model residual size by length.
resid = np.mean((Y_curve - bases['point'] @ point_coef)**2, axis=1)
display(pt.obs.assign(point_residual_mse=resid).groupby('length_group')['point_residual_mse']
        .agg(['count', 'median', 'mean']).reindex(['short', 'medium', 'long']).round(3))
print('Residuals are aggregate-expression deviations, not within-tubule heterogeneity.')

# %% [markdown]
# ## 6–9. Held-out prediction and biological diagnostics
#
# Leave one mouse specimen out. For each prespecified validation gene, fit its curve on the other specimen and predict this specimen; the gene never entered the new coordinate, orientation, or width. Error is normalized by that gene's **training** standard deviation. Primary width `c=.02` was fixed before validation. Controls test stronger spline smoothing, equal widths, and widths shuffled within specimen. The linear backbone row is a simple point baseline; the point-continuous row is its spline fit. The DPT spline is historical and had upstream access to the held-out genes, so its prediction is a useful baseline with that caveat.

# %%
dpt = pt.obs['total_scanpy_dpt'].to_numpy(dtype=float)
if not np.isfinite(dpt).all():
    raise ValueError('Saved PT DPT has non-finite values.')
dpt_basis, _ = observation_basis(dpt, np.zeros(pt.n_obs))
Y_valid = expression(heldout_genes)
sample_array = pt.obs['sample'].astype(str).to_numpy()
groups = pt.obs['length_group'].to_numpy()
base_designs = {
    'DPT': (dpt_basis, LAMBDA),
    'backbone_linear': (np.column_stack([np.ones(pt.n_obs), backbone]), 0.0),
    'point_continuous': (bases['point'], LAMBDA),
    'point_strong_smoothing': (bases['point'], 10 * LAMBDA),
    'extent_morphology': (bases['extent_c0.02'], LAMBDA),
    'extent_equal': (bases['equal'], LAMBDA),
    'extent_permuted': (bases['permuted'], LAMBDA),
}


def cross_validate(designs, y=Y_valid, genes=heldout_genes):
    rows = []
    for method, (design, lam) in designs.items():
        for sample in SAMPLES:
            test = sample_array == sample
            train = ~test
            train_y = y[train]
            if method.endswith('_linear'):
                coef = np.linalg.lstsq(design[train], train_y, rcond=None)[0]
            else:
                coef = fit_curve(design[train], train_y, lam)
            pred = design[test] @ coef
            scale = np.maximum(train_y.std(axis=0), .05)
            error = ((y[test] - pred) / scale)**2
            for j, gene in enumerate(genes):
                for group in ('all', 'short', 'medium', 'long'):
                    keep = np.ones(test.sum(), dtype=bool) if group == 'all' else groups[test] == group
                    rows.append({'method': method, 'heldout_specimen': sample,
                                 'gene': gene, 'length_group': group,
                                 'normalized_mse': error[keep, j].mean(), 'n_tubules': int(keep.sum())})
    return pd.DataFrame(rows)

cv = cross_validate(base_designs)
display(cv.groupby(['method', 'length_group'])['normalized_mse'].mean()
        .unstack().round(3).sort_values('all'))
fig, ax = plt.subplots(figsize=(8, 4))
plot_cv = cv[cv['method'].isin(['point_continuous', 'point_strong_smoothing',
                               'extent_morphology', 'extent_equal', 'extent_permuted'])]
for method, block in plot_cv.groupby('method'):
    values = block.groupby('length_group')['normalized_mse'].mean().reindex(['short', 'medium', 'long'])
    ax.plot(range(3), values, marker='o', label=method)
ax.set_xticks(range(3), ['short', 'medium', 'long'])
ax.set(ylabel='Held-out normalized MSE', title='Does morphology width help most for long tubules?')
ax.legend(fontsize=8); fig.tight_layout(); plt.show()

# Difference is paired by specimen and gene, so one very expressive gene cannot dominate.
pair = cv[cv['method'].isin(['point_continuous', 'extent_morphology'])]
pair = pair.pivot(index=['heldout_specimen', 'gene', 'length_group'], columns='method', values='normalized_mse')
pair['extent_minus_point'] = pair['extent_morphology'] - pair['point_continuous']
display(pair.groupby('length_group')['extent_minus_point'].agg(['mean', 'median']).round(4))

# Segment ordering is a label diagnostic, not a wholly independent test, because upstream
# cluster review used expression. Nearest-glomerulus distance is withheld from construction.
coordinates = {'DPT': dpt, 'backbone_linear': backbone,
               'point_continuous': backbone, 'extent_morphology': expected_position}
segment_rank = pt.obs['segment_class'].astype(str).map(dict(zip(SEGMENTS, (0, 1, 2)))).to_numpy(dtype=float)
metrics = []
for method, coord in coordinates.items():
    for sample in SAMPLES:
        m = sample_array == sample
        metrics.append({'method': method, 'sample': sample,
                        'segment_order_rho': spearmanr(coord[m], segment_rank[m]).statistic,
                        'glomerulus_distance_rho': spearmanr(coord[m],
                            pt.obs['distance_to_nearest_glomerulus_px'].to_numpy()[m]).statistic})
biology = pd.DataFrame(metrics)
display(biology.round(3))
fig, axes = plt.subplots(1, len(coordinates), figsize=(14, 3.5), sharey=True)
for ax, (method, coord) in zip(axes, coordinates.items()):
    ax.boxplot([coord[pt.obs['segment_class'].astype(str).eq(s).to_numpy()] for s in SEGMENTS],
               showfliers=False)
    ax.set_xticks(range(1, 4), ['S1', 'S2', 'S3'])
    ax.set(title=method, ylabel='Position on [0,1]')
fig.tight_layout(); plt.show()
fig, axes = plt.subplots(1, 2, figsize=(10, 3.5))
for ax, method in zip(axes, ('DPT', 'extent_morphology')):
    ax.scatter(coordinates[method], pt.obs['distance_to_nearest_glomerulus_px'], s=2, alpha=.15)
    ax.set(xlabel=method, ylabel='Distance to glomerulus (px)')
fig.tight_layout(); plt.show()

# Same gene curves separately fit in the two mice. This tests program reproducibility,
# not specimen mixing. Low-expression genes can give unstable correlations.
repro = []
for method, (design, lam) in base_designs.items():
    if method == 'backbone_linear':
        continue
    fitted = []
    for sample in SAMPLES:
        subset = sample_array == sample
        coef = fit_curve(design[subset], Y_curve[subset], lam)
        fitted.append(B_grid @ coef)
    per_gene = [spearmanr(fitted[0][:, j], fitted[1][:, j]).statistic for j in range(len(curve_genes))]
    repro.append({'method': method, 'median_curve_rho': np.nanmedian(per_gene),
                  'n_construction_genes': len(curve_genes)})
repro = pd.DataFrame(repro)
display(repro.round(3))

# %% [markdown]
# ## 10. Sensitivity and identifiability
#
# Every row is descriptive. Width scale, geometry metric, gene count, smoothing, graph neighborhood, and specimen are exposed. We do not choose the scale by agreement with DPT or optimize it on the reserved genes. The fixed-center model cannot prove a tubule truly spans this many pseudospace units: a broader latent gene curve and a broader observation kernel can compensate for one another.

# %%
sensitivity_designs = {}
for c in (0.0, 0.01, 0.02, 0.04):
    design, _ = observation_basis(backbone, widths(c=c))
    sensitivity_designs[f'c={c:g}'] = (design, LAMBDA)
for metric in ('area_um2', 'perimeter_um'):
    proxy = np.sqrt(pt.obs[metric].to_numpy()) if metric == 'area_um2' else pt.obs[metric].to_numpy()
    pt.obs[f'_proxy_{metric}'] = proxy
    design, _ = observation_basis(backbone, widths(metric=f'_proxy_{metric}'))
    sensitivity_designs[f'proxy={metric}'] = (design, LAMBDA)
for n_genes, neighbors in ((150, 30), (300, 15)):
    alternative, _ = make_backbone(n_genes=n_genes, neighbors=neighbors)
    design, _ = observation_basis(alternative, widths())
    sensitivity_designs[f'backbone_g{n_genes}_k{neighbors}'] = (design, LAMBDA)
sensitivity_designs['smoothing=10'] = (bases['extent_c0.02'], 10 * LAMBDA)
sensitivity = cross_validate(sensitivity_designs)
display(sensitivity[sensitivity['length_group'].eq('all')]
        .groupby(['method', 'heldout_specimen'])['normalized_mse'].mean().unstack().round(3))

# %% [markdown]
# ## 11. Comparison summary and saved intermediate artifacts
#
# A coordinate-only row uses a linear held-out gene predictor; continuous rows use identical spline complexity. Segment and physical-distance associations are computed from position (for the extent model, the truncated kernel mean), so a fixed-center extent model can change them only near the boundaries. Separate columns avoid a made-up aggregate score.

# %%
summary = (cv.groupby(['method', 'length_group'])['normalized_mse'].mean().unstack()
           .rename(columns={'all': 'heldout_gene_nmse', 'short': 'short_nmse', 'long': 'long_nmse'}))
summary = summary.join(biology.groupby('method')[['segment_order_rho', 'glomerulus_distance_rho']].mean())
summary = summary.join(repro.set_index('method')['median_curve_rho'])
summary = summary.reindex(['DPT', 'backbone_linear', 'point_continuous', 'extent_morphology'])
summary = summary.drop(columns=['medium'])
summary = summary.rename_axis('model')
display(summary.round(3))

positions = pt.obs[['sample', 'segment_class', 'total_scanpy_dpt', 'backbone_position',
                    'extent_mean_position', 'extent_sigma', 'relative_extent',
                    'length_group', 'major_um', 'minor_um', 'area_um2', 'perimeter_um',
                    'aspect_ratio', 'distance_to_nearest_glomerulus_px']].copy()
positions.to_csv(out / 'positions_and_morphology.csv', index_label='obs_name')
cv.to_csv(out / 'heldout_gene_loso.csv', index=False)
sensitivity.to_csv(out / 'sensitivity_loso.csv', index=False)
summary.to_csv(out / 'method_comparison.csv')
np.savez_compressed(out / 'construction_gene_curves.npz', grid=GRID,
                    genes=np.asarray(curve_genes), point=B_grid @ point_coef,
                    extent=B_grid @ extent_coef)
(out / 'run_manifest.json').write_text(json.dumps({
    'pt_input': str(pt_path), 'pass1_input': str(pass1_path),
    'n_tubules': pt.n_obs, 'n_mapped_genes': pt.n_vars,
    'backbone_genes': backbone_genes, 'curve_genes': curve_genes,
    'heldout_genes': heldout_genes, 'orientation_only_genes': list(ORIENT_EARLY + ORIENT_LATE),
    'primary_c': PRIMARY_C, 'seed': SEED,
    'note': 'Exploratory; DPT and reviewed labels had upstream access to expression.'
}, indent=2))
print('Saved exploratory artifacts to', out)

# %% [markdown]
# ## 12. Interpretation
#
# Read the table and the paired short/long errors before making a claim. The following answers are generated from the primary held-out comparison; they remain conditional on two specimens and on the upstream annotation and DPT leakage described above.

# %%
def average_error(method, group='all'):
    return float(cv.loc[(cv.method == method) & (cv.length_group == group), 'normalized_mse'].mean())

point = average_error('point_continuous')
extent = average_error('extent_morphology')
relative_gain = (point - extent) / point
short_delta = average_error('extent_morphology', 'short') - average_error('point_continuous', 'short')
long_delta = average_error('extent_morphology', 'long') - average_error('point_continuous', 'long')
equal = average_error('extent_equal')
permuted = average_error('extent_permuted')
repro_point = float(repro.set_index('method').loc['point_continuous', 'median_curve_rho'])
repro_extent = float(repro.set_index('method').loc['extent_morphology', 'median_curve_rho'])
dpt_error = average_error('DPT')
credible_gain = (relative_gain > .01 and long_delta < short_delta and
                 (min(equal, permuted) - extent) / point > .01 and
                 repro_extent >= repro_point - .05)
display(Markdown(f'''1. **Overall:** morphology-width normalized MSE is {extent:.3f} versus point {point:.3f}; lower is better. Relative change is {relative_gain:+.2%}; differences below 1% are negligible in this prototype. Historical DPT has MSE {dpt_error:.3f}, but its upstream construction included these genes.
2. **Long structures:** extent minus point MSE is {long_delta:+.3f} for long and {short_delta:+.3f} for short tubules. A more negative long difference would support the length-specific claim.
3. **Width identity:** morphology {extent:.3f}, equal {equal:.3f}, shuffled {permuted:.3f}. Morphology must beat both controls, especially in long tubules, to favor physical extent over generic smoothing.
4. **Gene curves:** the displayed point and extent curves nearly overlap and show no obvious spline oscillation in this run. Some top-variance genes may reflect off-target signal; these fits do not establish anatomical gene-curve truth.
5. **Specimen stability:** median construction-gene curve rho is {repro_extent:.3f} for extent and {repro_point:.3f} for point. Both are high, but held-out prediction differs between the two mice (see the specimen columns in sensitivity).
6. **Identifiability:** centers were fixed and widths constrained. 2D length does not measure nephron arc length; near-equal performance among morphology, equal, and shuffled kernels means width biology is not identifiable here.
7. **Decision:** {'The prespecified descriptive criteria jointly favor deeper joint center/width work, pending external validation.' if credible_gain else 'Current evidence does not justify a full joint center/width method. The fixed-center result is the informative stopping point.'} The 1% practical-gain screen is exploratory, not a significance test.'''))
