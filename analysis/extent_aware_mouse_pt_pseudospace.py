# %% [markdown]
# # Extent-aware human–mouse PT pseudospace: fixed-center prototype
#
# **Question.** Does a segmented tubule behave more like a point or an observation that averages a finite stretch of the proximal-tubule axis? This is an exploratory diagnostic, not a replacement for the reviewed DPT workflow.
#
# Inputs are notebook 03's shared PT `lognorm` matrix and PT-specific DPT, pass-1 annotations, and the exact source GeoJSON polygons for two control mice and two human sections. The section named `HUK1_MED1` is healthy **cortex** from the same donor as `HUK1_COR1`. The new backbone is learned jointly from genes measured in both species, independently of DPT; all input artifacts are read-only. Results go to `results/pt_extent_aware_human_mouse/`.
#
# **Predeclared comparison.** Nine PT genes are excluded from new coordinate and width construction. Published early/late markers orient the backbone only. For each held-out gene, one specimen predicts the other specimen of the same species; a separate cross-species transfer check is descriptive. Segment labels and glomerular distance are diagnostics. The reviewed labels and historical DPT had upstream access to expression, so they are not fully external validation. The two human sections are one donor, not two biological replicates.

# %% [markdown]
# ## 1. Inputs and baseline
#
# Notebook 03's complete PT object contains both species on its shared PT DPT. Pass-1 provides `feature_index` for exact joins to each sample's source GeoJSON and glomerular centroids for distance. Missing or ambiguous joins stop the run. Distances are computed within specimen and never compared in absolute pixels across specimens.

# %%
import argparse
import json
import os
import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csgraph
from scipy.spatial import cKDTree
from scipy.stats import spearmanr
from shapely.geometry import shape
from sklearn.decomposition import PCA
from sklearn.neighbors import kneighbors_graph
from sklearn.preprocessing import SplineTransformer, StandardScaler
from IPython.display import display, Image, Markdown

SEED = 17
rng = np.random.default_rng(SEED)
SPECIES_SAMPLES = {'mouse': ('Ctrl1A2', 'Ctrl1A4'),
                   'human': ('HUK1_COR1', 'HUK1_MED1')}
SAMPLES = tuple(sample for pair in SPECIES_SAMPLES.values() for sample in pair)
GEOJSON_FILES = {**{s: f'{s}_kept_tubules_labeled_fine.geojson' for s in SPECIES_SAMPLES['mouse']},
                 **{s: f'{s}_v2.geojson' for s in SPECIES_SAMPLES['human']}}
SPECIES_COLORS = {'mouse': '#0072B2', 'human': '#D55E00'}
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
if str(root) not in sys.path:
    sys.path.insert(0, str(root))
from pseudospace.heatmaps import plot_marker_heatmap
parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--data-root', type=Path)
parser.add_argument('--results-root', type=Path)
args, _ = parser.parse_known_args()
data_root = Path(args.data_root or os.environ.get('PSEUDOSPACE_DATA_ROOT') or root / 'data').expanduser().resolve()
results_root = Path(args.results_root or os.environ.get('PSEUDOSPACE_RESULTS_ROOT') or root / 'results').expanduser().resolve()
upstream = results_root / 'human_vs_healthy_mouse'
out = results_root / 'pt_extent_aware_human_mouse'
pt_path = upstream / 'cross_species_pt_dpt.h5ad'
pass1_path = upstream / 'cross_species_harmony_pass1.h5ad'
for path in (pt_path, pass1_path, *(data_root / GEOJSON_FILES[s] for s in SAMPLES)):
    if not path.is_file():
        raise FileNotFoundError(f'Required upstream input missing: {path}')
out.mkdir(parents=True, exist_ok=True)

pt = ad.read_h5ad(pt_path)
if pt.n_obs == 0 or 'lognorm' not in pt.layers:
    raise ValueError('Saved shared PT artifact has no structures or lognorm layer.')
for species, samples in SPECIES_SAMPLES.items():
    mask = pt.obs['sample'].astype(str).isin(samples)
    if not mask.any() or not pt.obs.loc[mask, 'comparison_species'].astype(str).eq(species).all():
        raise ValueError(f'{species} sample/species mismatch in saved PT object.')
if not pt.obs['sample'].astype(str).isin(SAMPLES).all():
    raise ValueError('Unexpected PT specimen in saved object.')
if not pt.obs.loc[pt.obs['comparison_species'].astype(str).eq('human'), 'region'].astype(str).eq('cortex').all():
    raise ValueError('Human PT sections are expected to be healthy cortex, including HUK1_MED1.')
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
print(f'Input: {pt.n_obs:,} shared human–mouse PT tubules, {pt.n_vars:,} mapped genes; samples:')
display(pt.obs.groupby(['comparison_species', 'sample'], observed=True).size().rename('tubules').to_frame())
print('Both human sections are healthy cortex from one donor.')
print('Measured in both inputs:', int(pt.var['measured_in_both_inputs'].sum()))

# %% [markdown]
# ## 2. Polygon morphology and an extent prior
#
# All four source GeoJSON files contain polygon outlines and calibrated area, but no skeleton or traced nephron path. We calculate rotated-box sides, perimeter, area, aspect ratio, and axis-aligned bounds. The long side is a 2D relative extent proxy, normalized **within specimen** before assigning widths. Sectioning angle and tortuosity prevent conversion to true nephron arc length.

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
    with (data_root / GEOJSON_FILES[sample]).open() as handle:
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
# Specimen median 1, conservative range [0.5, 1.5]; no micron-to-pseudotime conversion.
length_ratio = (pt.obs['major_um'] / pt.obs.groupby('sample', observed=True)['major_um'].transform('median')).to_numpy()
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
fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharex='col')
for row, species in enumerate(SPECIES_SAMPLES):
    sub = pt.obs['comparison_species'].astype(str).eq(species)
    for seg in SEGMENTS:
        m = sub & pt.obs['segment_class'].astype(str).eq(seg)
        axes[row, 0].hist(pt.obs.loc[m, 'total_scanpy_dpt'], bins=35, density=True, alpha=.45, label=seg)
    axes[row, 0].legend(fontsize=8)
    axes[row, 1].scatter(pt.obs.loc[sub, 'total_scanpy_dpt'],
                         pt.obs.loc[sub, 'distance_to_nearest_glomerulus_px'], s=2, alpha=.17)
    axes[row, 2].scatter(pt.obs.loc[sub, 'total_scanpy_dpt'], pt.obs.loc[sub, 'major_um'], s=2, alpha=.17)
    axes[row, 0].set_ylabel(species.title())
for ax, title in zip(axes[0], ('DPT by reviewed segment', 'DPT vs glomerulus distance', 'DPT vs 2D long axis')):
    ax.set_title(title)
for ax in axes[1]: ax.set_xlabel('Existing shared PT DPT')
fig.tight_layout(); plt.show()
display(pd.crosstab([pt.obs['comparison_species'], pt.obs['length_group']],
                    pt.obs['segment_class'], normalize='index').round(3))

# %% [markdown]
# ## Notebook 03 marker dot plots and reference heatmaps
#
# Notebook 03's cluster-review dot plot and DE heatmap are shown as **upstream context**: they helped assign the saved labels and cannot validate this new method. The new PT dot plots use that notebook's S1/S2/S3 marker panel, split by species and by morphology-length group. These same markers are excluded from backbone-gene selection.

# %%
for filename in ('celltyping/cluster_reference_gene_dotplot.png',
                 'celltyping/coarse_cluster_de_heatmap.png',
                 'curves/paired_top_gene_heatmap_lognorm.png',
                 'curves/paired_top_gene_heatmap_shape.png'):
    figure = upstream / filename
    if figure.is_file():
        display(Markdown(f'**Notebook 03 reference:** `{filename}`'))
        display(Image(filename=str(figure), width=900))
    else:
        print('Notebook 03 reference figure unavailable:', figure)

PT_MARKERS = {'PT-S1': ['Slc5a2', 'Slc5a12', 'Gatm'],
              'PT-S2': ['Slc22a6', 'Slc13a3', 'Cyp2e1'],
              'PT-S3': ['Slc22a7', 'Cyp7b1']}
PT_COLORS = {'PT-S1': '#4C9BD3', 'PT-S2': '#5B8E55', 'PT-S3': '#D6B48A'}
for suffix, labels in {
    'species_segment': [f'{sp} {seg}' for sp in SPECIES_SAMPLES for seg in SEGMENTS],
    'species_length': [f'{sp} {length}' for sp in SPECIES_SAMPLES for length in ('short', 'medium', 'long')],
}.items():
    second = pt.obs['segment_class'] if suffix == 'species_segment' else pt.obs['length_group']
    pt.obs[suffix] = pd.Categorical(pt.obs['comparison_species'].astype(str) + ' ' + second.astype(str),
                                    categories=labels, ordered=True)
    dot = sc.pl.dotplot(pt, PT_MARKERS, groupby=suffix, layer='lognorm',
                         standard_scale='var', show=False, return_fig=True)
    dot.savefig(out / f'pt_marker_dotplot_{suffix}.png')
    dot.show()

# %% [markdown]
# ## 3. Joint molecular backbone
#
# We select genes measured in **both** inputs by average within-species variance, excluding reserved validation and published PT annotation/orientation genes. Each gene is standardized within species before joint PCA so species-level abundance offsets do not dominate the graph; this is a transparent alignment assumption and can hide genuine species shifts. A symmetric kNN graph and shortest path from an early-marker root yield one shared point coordinate. DPT, segment labels, glomerular distance, and morphology do not enter it.

# %%
excluded = set(HELD_OUT) | set(LABEL_PANEL) | set(ORIENT_EARLY) | set(ORIENT_LATE)
available = pt.var['measured_in_both_inputs'].to_numpy(dtype=bool)
all_expr = pt.layers['lognorm'].tocsr()
species_array = pt.obs['comparison_species'].astype(str).to_numpy()
within_variance = []
within_fraction = []
for species in SPECIES_SAMPLES:
    matrix = all_expr[species_array == species]
    mean = np.asarray(matrix.mean(axis=0)).ravel()
    within_variance.append(np.asarray(matrix.power(2).mean(axis=0)).ravel() - mean**2)
    within_fraction.append(np.asarray((matrix > 0).mean(axis=0)).ravel())
variance = np.mean(within_variance, axis=0)
fraction = np.min(within_fraction, axis=0)
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
    z = np.empty_like(both)
    for species in SPECIES_SAMPLES:
        m = species_array == species
        z[m] = StandardScaler().fit_transform(both[m])
    return z[:, len(ORIENT_EARLY):].mean(axis=1) - z[:, :len(ORIENT_EARLY)].mean(axis=1)

orient = orientation_score()


def make_backbone(n_genes=300, neighbors=30):
    x = expression(backbone_genes[:n_genes])
    standardized = np.empty_like(x)
    candidates = []
    for species in SPECIES_SAMPLES:
        m = species_array == species
        standardized[m] = StandardScaler().fit_transform(x[m])
        candidates.extend(np.flatnonzero(m & (orient <= np.quantile(orient[m], .03))))
    pcs = PCA(n_components=12, random_state=SEED).fit_transform(standardized)
    candidates = np.asarray(candidates, dtype=int)
    center = np.median(pcs[candidates], axis=0)
    root_index = candidates[np.argmin(np.linalg.norm(pcs[candidates] - center, axis=1))]
    graph = kneighbors_graph(pcs, n_neighbors=neighbors, mode='distance', include_self=False)
    graph = graph.maximum(graph.T).tocsr()
    edge_row, edge_col = graph.nonzero()
    n_cross = int(np.sum(species_array[edge_row] != species_array[edge_col]))
    cross_fraction = n_cross / len(edge_row)
    print(f'Joint graph: {n_genes} genes, k={neighbors}, {n_cross:,}/{len(edge_row):,} '
          f'({cross_fraction:.1%}) cross-species edges.')
    component_count = csgraph.connected_components(graph, directed=False)[0]
    if component_count != 1:
        raise ValueError(f'kNN graph has {component_count} components at k={neighbors}; inspect topology.')
    path = csgraph.dijkstra(graph, directed=False, indices=int(root_index))
    if not np.isfinite(path).all():
        raise ValueError('Backbone contains unreachable PT structures.')
    coordinate = (path - path.min()) / (path.max() - path.min())
    if spearmanr(coordinate, orient).statistic < 0:
        coordinate = 1 - coordinate
    return coordinate, pcs, cross_fraction

backbone, pcs, cross_edge_fraction = make_backbone()
pt.obs['backbone_position'] = backbone
print('Backbone vs existing DPT: Spearman rho =',
      round(spearmanr(backbone, pt.obs['total_scanpy_dpt']).statistic, 3))
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
for species in SPECIES_SAMPLES:
    m = species_array == species
    axes[0].scatter(pt.obs['total_scanpy_dpt'].to_numpy()[m], backbone[m], s=2, alpha=.15,
                    color=SPECIES_COLORS[species], label=species)
axes[0].legend(markerscale=4)
axes[0].set(xlabel='Existing shared PT DPT', ylabel='Joint backbone')
axes[1].scatter(pcs[:, 0], pcs[:, 1], c=backbone, s=2, cmap='viridis')
axes[1].set(xlabel='PC1', ylabel='PC2', title='Graph geodesic on 12 PCs')
fig.tight_layout(); plt.show()

# Restrict paired species figures to support observed in both species; no extrapolated
# tail is shown as an inferred conserved program.
backbone_support = {sp: np.quantile(backbone[species_array == sp], [.01, .99])
                   for sp in SPECIES_SAMPLES}
common_low = max(bounds[0] for bounds in backbone_support.values())
common_high = min(bounds[1] for bounds in backbone_support.values())
if common_high <= common_low:
    raise ValueError('Joint backbone has no mouse–human 1–99% common support.')
print('Backbone 1–99% support:', backbone_support, 'shared:', (common_low, common_high))
fig, ax = plt.subplots(figsize=(9, 3))
for species in SPECIES_SAMPLES:
    ax.hist(backbone[species_array == species], bins=60, density=True,
            alpha=.4, color=SPECIES_COLORS[species], label=species)
ax.axvspan(common_low, common_high, color='0.85', alpha=.3, label='shared 1–99% support')
ax.set(xlabel='Joint backbone position', ylabel='Density', title='Species occupancy of the joint backbone')
ax.legend(); fig.tight_layout(); plt.show()

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
    medians = pt.obs.groupby('sample', observed=True)[metric].transform('median').to_numpy(dtype=float)
    relative = np.clip(values / medians, .5, 1.5)
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

# Notebook 03's same marker-gradient panel, now shown both on its shared PT DPT and
# on the new distribution's mean coordinate. Each species is plotted separately;
# the label strip is contextual, not independent evidence.
marker_panel = {group: [{'label': gene, 'candidates': [gene]} for gene in genes]
                for group, genes in PT_MARKERS.items()}
for position in ('total_scanpy_dpt', 'extent_mean_position'):
    for species in SPECIES_SAMPLES:
        plot_marker_heatmap(
            pt, marker_panel, list(PT_MARKERS), PT_COLORS,
            pseudotime_col=position, cluster_col=None,
            mask=species_array == species,
            title=f'{species.title()} PT markers on {position}',
            output_name=f'{species}_{position}_marker_heatmap.png',
            strip_col='segment_class', strip_order=list(SEGMENTS), strip_colors=PT_COLORS,
            n_bins=120, output_dir=out, project_dir=root,
        )
fig, ax = plt.subplots(figsize=(9, 3.5))
for group in ('short', 'medium', 'long'):
    idx = np.flatnonzero((pt.obs['length_group'].eq(group) &
                         pt.obs['comparison_species'].astype(str).eq('mouse')).to_numpy())
    selected = idx[np.argmin(abs(backbone[idx] - .5))]
    sigma = pt.obs['extent_sigma'].iloc[selected]
    density = np.exp(-.5 * ((GRID - backbone[selected]) / sigma)**2)
    ax.plot(GRID, density / np.trapz(density, GRID), label=f'{group}: σ={sigma:.3f}')
ax.set(xlabel='Joint backbone position', ylabel='Truncated density',
       title='Representative mouse tubules; width normalized within specimen')
ax.legend(); fig.tight_layout(); plt.show()

# Species-specific curves preserve expression-level differences on the common backbone.
# Held-out gene validation occurs below with specimen separation.
Y_curve = expression(curve_genes)
curve_coefficients = {}
for species in SPECIES_SAMPLES:
    m = species_array == species
    curve_coefficients[species] = {
        'point': fit_curve(bases['point'][m], Y_curve[m]),
        'extent': fit_curve(bases['extent_c0.02'][m], Y_curve[m]),
    }
point_coef = curve_coefficients['mouse']['point']
extent_coef = curve_coefficients['mouse']['extent']
curve_plot_genes = curve_genes[:4]
fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
for ax, gene in zip(axes.flat, curve_plot_genes):
    j = curve_genes.index(gene)
    mouse = species_array == 'mouse'
    ax.scatter(backbone[mouse], Y_curve[mouse, j], s=2, alpha=.12, color='0.35')
    mouse_grid = (GRID >= backbone_support['mouse'][0]) & (GRID <= backbone_support['mouse'][1])
    ax.plot(GRID[mouse_grid], (B_grid @ point_coef[:, j])[mouse_grid], label='Point curve')
    ax.plot(GRID[mouse_grid], (B_grid @ extent_coef[:, j])[mouse_grid], label='Extent underlying curve')
    ax.set(title=gene, ylabel='lognorm expression')
axes[0, 0].legend(); axes[1, 0].set_xlabel('Backbone position'); axes[1, 1].set_xlabel('Backbone position')
fig.tight_layout(); plt.show()

# Notebook 03's paired-heatmap convention: common expression scale preserves level,
# and within-curve z-score isolates shape. Genes are the predeclared top 24
# construction-variance genes, not selected by human–mouse divergence or DPT.
heat_genes = curve_genes[:24]
heat_grid = GRID[(GRID >= common_low) & (GRID <= common_high)]
heat_basis = spline.transform(heat_grid[:, None])
paired = {sp: heat_basis @ curve_coefficients[sp]['extent'][:, :len(heat_genes)]
          for sp in SPECIES_SAMPLES}
for view in ('common_lognorm', 'shape_only'):
    values = paired if view == 'common_lognorm' else {
        sp: np.clip((matrix - matrix.mean(axis=0)) /
                    np.maximum(matrix.std(axis=0), 1e-8), -2.5, 2.5)
        for sp, matrix in paired.items()
    }
    low, high = (min(matrix.min() for matrix in values.values()),
                 max(matrix.max() for matrix in values.values())) if view == 'common_lognorm' else (-2.5, 2.5)
    fig, axes = plt.subplots(1, 2, figsize=(10, 7), sharey=True)
    for ax, species in zip(axes, SPECIES_SAMPLES):
        im = ax.imshow(values[species].T, aspect='auto', interpolation='nearest',
                       cmap='viridis' if view == 'common_lognorm' else 'bwr',
                       vmin=low, vmax=high, extent=[heat_grid[0], heat_grid[-1], len(heat_genes)-.5, -.5])
        ax.set(title=species.title(), xlabel='Joint backbone position')
    axes[0].set_yticks(np.arange(len(heat_genes)), heat_genes, fontsize=8)
    fig.colorbar(im, ax=axes, shrink=.7, label='lognorm' if view == 'common_lognorm' else 'Within-curve z-score')
    fig.suptitle(f'Species-specific extent gene curves on shared support: {view}')
    fig.savefig(out / f'paired_extent_gene_heatmap_{view}.png', dpi=180, bbox_inches='tight')
    plt.show()


# For a long tubule near the steepest fitted gradient, compare the fitted point value
# with the actual integrated prediction. The red mark is its measured aggregate.
j = int(np.argmax(np.ptp(B_grid @ extent_coef, axis=0)))
gene = curve_genes[j]
gradient = np.abs(np.gradient(B_grid @ extent_coef[:, j], GRID))
steep = GRID[np.argmax(gradient)]
long_idx = np.flatnonzero((pt.obs['length_group'].eq('long') &
                               pt.obs['comparison_species'].astype(str).eq('mouse')).to_numpy())
examples = long_idx[np.argsort(abs(backbone[long_idx] - steep))[:8]]
fig, ax = plt.subplots(figsize=(8, 4))
ax.scatter(backbone[mouse], Y_curve[mouse, j], s=2, alpha=.08, color='0.5')
ax.plot(GRID[mouse_grid], (B_grid @ extent_coef[:, j])[mouse_grid],
        color='#0072B2', label='Underlying extent curve')
ax.scatter(backbone[examples], Y_curve[examples, j], color='#D55E00', s=24, label='Observed long tubule')
ax.scatter(backbone[examples], (bases['extent_c0.02'] @ extent_coef)[examples, j],
           marker='x', color='black', s=40, label='Integrated prediction')
ax.set(xlabel='Backbone position', ylabel='lognorm expression',
       title=f'{gene}: long tubules near steepest fitted change')
ax.legend(fontsize=8); fig.tight_layout(); plt.show()

# A diagnostic for the proposed integration mechanism: point-model residual size by length.
point_pred = np.empty_like(Y_curve)
for species in SPECIES_SAMPLES:
    m = species_array == species
    point_pred[m] = bases['point'][m] @ curve_coefficients[species]['point']
resid = np.mean((Y_curve - point_pred)**2, axis=1)
display(pt.obs.assign(point_residual_mse=resid).groupby(['comparison_species', 'length_group'], observed=True)['point_residual_mse']
        .agg(['count', 'median', 'mean']).round(3))
print('Residuals are aggregate-expression deviations, not within-tubule heterogeneity.')

# %% [markdown]
# ## 6–9. Held-out prediction and biological diagnostics
#
# Leave one specimen out **within each species**. For each prespecified validation gene, fit on the other specimen/section of that species and predict the held-out one; the gene never entered new coordinate, orientation, or width. Error is normalized by its training standard deviation. Human section-to-section agreement is within one donor. A separate mouse→human and human→mouse transfer check asks whether raw expression curves transport across species; species-level abundance differences can dominate that check. DPT is historical and had upstream access to held-out genes.

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


def cross_validate(designs, y=Y_valid, genes=heldout_genes, shared_support=True):
    rows = []
    for method, (design, lam) in designs.items():
        for species, samples in SPECIES_SAMPLES.items():
            for sample in samples:
                test = sample_array == sample
                train = (species_array == species) & ~test
                train_y = y[train]
                if method.endswith('_linear'):
                    coef = np.linalg.lstsq(design[train], train_y, rcond=None)[0]
                else:
                    coef = fit_curve(design[train], train_y, lam)
                pred = design[test] @ coef
                scale = np.maximum(train_y.std(axis=0), .05)
                error = ((y[test] - pred) / scale)**2
                in_support = np.ones(test.sum(), dtype=bool)
                if shared_support:
                    for coordinate in (dpt, backbone):
                        low, high = np.quantile(coordinate[train], [.01, .99])
                        in_support &= (coordinate[test] >= low) & (coordinate[test] <= high)
                for j, gene in enumerate(genes):
                    for group in ('all', 'short', 'medium', 'long'):
                        keep = in_support if group == 'all' else in_support & (groups[test] == group)
                        if not keep.any():
                            raise ValueError(f'No common-support {group} tubules for {sample}.')
                        rows.append({'method': method, 'species': species, 'heldout_specimen': sample,
                                     'gene': gene, 'length_group': group,
                                     'normalized_mse': error[keep, j].mean(), 'n_tubules': int(keep.sum())})
    return pd.DataFrame(rows)

# Primary scores use the same test tubules for every method: intersection of the
# training specimen's 1–99% DPT and backbone support. Full-range errors remain as an audit.
support_rows = []
for sample in SAMPLES:
    species = species_array[sample_array == sample][0]
    train = (species_array == species) & (sample_array != sample)
    test = sample_array == sample
    for name, coordinate in (('DPT', dpt), ('joint_backbone', backbone)):
        low, high = np.quantile(coordinate[train], [.01, .99])
        support_rows.append({'sample': sample, 'coordinate': name,
                             'fraction_test_outside_training_support':
                             float(np.mean((coordinate[test] < low) | (coordinate[test] > high)))})
support_audit = pd.DataFrame(support_rows)
display(support_audit.round(3))
cv = cross_validate(base_designs)
cv_all_positions = cross_validate(base_designs, shared_support=False)
print('Primary CV uses common training support; full-range results are saved separately.')
display(cv.groupby(['species', 'method', 'length_group'], observed=True)['normalized_mse']
        .mean().unstack().round(3))
fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharey=False)
plot_cv = cv[cv['method'].isin(['point_continuous', 'point_strong_smoothing',
                               'extent_morphology', 'extent_equal', 'extent_permuted'])]
for ax, species in zip(axes, SPECIES_SAMPLES):
    for method, block in plot_cv[plot_cv.species.eq(species)].groupby('method'):
        values = block.groupby('length_group')['normalized_mse'].mean().reindex(['short', 'medium', 'long'])
        ax.plot(range(3), values, marker='o', label=method)
    ax.set_xticks(range(3), ['short', 'medium', 'long'])
    ax.set(ylabel='Held-out normalized MSE', title=f'{species.title()}: does length matter?')
axes[1].legend(fontsize=7, bbox_to_anchor=(1.02, 1), loc='upper left')
fig.tight_layout(); plt.show()

# Difference is paired by specimen and gene, so one very expressive gene cannot dominate.
pair = cv[cv['method'].isin(['point_continuous', 'extent_morphology'])]
pair = pair.pivot(index=['species', 'heldout_specimen', 'gene', 'length_group'], columns='method', values='normalized_mse')
pair['extent_minus_point'] = pair['extent_morphology'] - pair['point_continuous']
display(pair.groupby(['species', 'length_group'])['extent_minus_point'].agg(['mean', 'median']).round(4))

# Segment ordering is a label diagnostic, not a wholly independent test, because upstream
# cluster review used expression. Nearest-glomerulus distance is withheld from construction.
coordinates = {'DPT': dpt, 'backbone_linear': backbone,
               'point_continuous': backbone, 'extent_morphology': expected_position}
segment_rank = pt.obs['segment_class'].astype(str).map(dict(zip(SEGMENTS, (0, 1, 2)))).to_numpy(dtype=float)
metrics = []
for method, coord in coordinates.items():
    for sample in SAMPLES:
        m = sample_array == sample
        metrics.append({'method': method, 'species': species_array[m][0], 'sample': sample,
                        'segment_order_rho': spearmanr(coord[m], segment_rank[m]).statistic,
                        'glomerulus_distance_rho': spearmanr(coord[m],
                            pt.obs['distance_to_nearest_glomerulus_px'].to_numpy()[m]).statistic})
biology = pd.DataFrame(metrics)
display(biology.round(3))
fig, axes = plt.subplots(2, len(coordinates), figsize=(14, 7), sharey=True)
for row, species in enumerate(SPECIES_SAMPLES):
    m = species_array == species
    for ax, (method, coord) in zip(axes[row], coordinates.items()):
        ax.boxplot([coord[m & pt.obs['segment_class'].astype(str).eq(seg).to_numpy()]
                    for seg in SEGMENTS], showfliers=False)
        ax.set_xticks(range(1, 4), ['S1', 'S2', 'S3'])
        ax.set(title=f'{species}: {method}', ylabel='Position on [0,1]')
fig.tight_layout(); plt.show()
fig, axes = plt.subplots(2, 2, figsize=(10, 7))
for row, species in enumerate(SPECIES_SAMPLES):
    m = species_array == species
    for ax, method in zip(axes[row], ('DPT', 'extent_morphology')):
        ax.scatter(coordinates[method][m],
                   pt.obs['distance_to_nearest_glomerulus_px'].to_numpy()[m], s=2, alpha=.15)
        ax.set(xlabel=f'{species}: {method}', ylabel='Nearest glomerulus (px)')
fig.tight_layout(); plt.show()

# Fit the same construction-gene curves separately to each specimen/section.
# Human agreement reflects repeat sections of one donor, not donor replication.
repro = []
for method, (design, lam) in base_designs.items():
    if method == 'backbone_linear':
        continue
    for species, samples in SPECIES_SAMPLES.items():
        fitted = []
        for sample in samples:
            subset = sample_array == sample
            coef = fit_curve(design[subset], Y_curve[subset], lam)
            fitted.append(B_grid @ coef)
        per_gene = [spearmanr(fitted[0][:, j], fitted[1][:, j]).statistic
                    for j in range(len(curve_genes))]
        repro.append({'species': species, 'method': method,
                      'median_curve_rho': np.nanmedian(per_gene),
                      'n_construction_genes': len(curve_genes)})
repro = pd.DataFrame(repro)
display(repro.round(3))

# Cross-species transfer is a secondary diagnostic on raw lognorm expression.
# It is sensitive to species-level expression shifts as well as coordinate quality.
transfer_rows = []
for method in ('DPT', 'point_continuous', 'extent_morphology'):
    design, lam = base_designs[method]
    for source in SPECIES_SAMPLES:
        target = next(sp for sp in SPECIES_SAMPLES if sp != source)
        train, test = species_array == source, species_array == target
        coef = fit_curve(design[train], Y_valid[train], lam)
        scale = np.maximum(Y_valid[train].std(axis=0), .05)
        error = ((Y_valid[test] - design[test] @ coef) / scale)**2
        transfer_rows.append({'method': method, 'source_species': source,
                              'target_species': target,
                              'normalized_mse': float(error.mean())})
transfer = pd.DataFrame(transfer_rows)
display(transfer.round(3))

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
    alternative, _, _ = make_backbone(n_genes=n_genes, neighbors=neighbors)
    design, _ = observation_basis(alternative, widths())
    sensitivity_designs[f'backbone_g{n_genes}_k{neighbors}'] = (design, LAMBDA)
sensitivity_designs['smoothing=10'] = (bases['extent_c0.02'], 10 * LAMBDA)
sensitivity = cross_validate(sensitivity_designs)
display(sensitivity[sensitivity['length_group'].eq('all')]
        .groupby(['method', 'heldout_specimen'])['normalized_mse'].mean().unstack().round(3))

# %% [markdown]
# ## 11. Comparison summary and saved intermediate artifacts
#
# A coordinate-only row uses a linear held-out gene predictor; continuous rows use identical spline complexity. Segment and physical-distance associations come from position (for extent, the truncated-kernel mean), so fixed centers change them only near boundaries. Primary prediction scores use the same shared training support for every model; full-range scores and the excluded fraction are saved for audit. Separate columns avoid an aggregate winner score.

# %%
summary = (cv.groupby(['species', 'method', 'length_group'])['normalized_mse'].mean().unstack()
           .rename(columns={'all': 'heldout_gene_nmse', 'short': 'short_nmse', 'long': 'long_nmse'}))
summary = summary.join(biology.groupby(['species', 'method'])[['segment_order_rho', 'glomerulus_distance_rho']].mean())
summary = summary.join(repro.set_index(['species', 'method'])['median_curve_rho'])
summary = summary.reindex(pd.MultiIndex.from_product(
    [SPECIES_SAMPLES, ('DPT', 'backbone_linear', 'point_continuous', 'extent_morphology')],
    names=['species', 'model']))
summary = summary.drop(columns=['medium'])
display(summary.round(3))

positions = pt.obs[['comparison_species', 'sample', 'region', 'segment_class', 'total_scanpy_dpt', 'backbone_position',
                    'extent_mean_position', 'extent_sigma', 'relative_extent',
                    'length_group', 'major_um', 'minor_um', 'area_um2', 'perimeter_um',
                    'aspect_ratio', 'distance_to_nearest_glomerulus_px']].copy()
positions.to_csv(out / 'positions_and_morphology.csv', index_label='obs_name')
cv.to_csv(out / 'heldout_gene_loso.csv', index=False)
cv_all_positions.to_csv(out / 'heldout_gene_loso_all_positions.csv', index=False)
sensitivity.to_csv(out / 'sensitivity_loso.csv', index=False)
summary.to_csv(out / 'method_comparison.csv')
support_audit.to_csv(out / 'cv_support_audit.csv', index=False)
transfer.to_csv(out / 'cross_species_transfer.csv', index=False)
np.savez_compressed(out / 'construction_gene_curves.npz', grid=GRID,
                    genes=np.asarray(curve_genes),
                    **{f'{species}_{method}': B_grid @ coefficient
                       for species, fits in curve_coefficients.items()
                       for method, coefficient in fits.items()})
(out / 'run_manifest.json').write_text(json.dumps({
    'pt_input': str(pt_path), 'pass1_input': str(pass1_path),
    'n_tubules': pt.n_obs, 'n_mapped_genes': pt.n_vars,
    'species_samples': SPECIES_SAMPLES, 'geometry_sources': GEOJSON_FILES,
    'backbone_genes': backbone_genes, 'curve_genes': curve_genes,
    'heldout_genes': heldout_genes, 'orientation_only_genes': list(ORIENT_EARLY + ORIENT_LATE),
    'primary_c': PRIMARY_C, 'seed': SEED,
    'cross_species_knn_edge_fraction': cross_edge_fraction,
    'note': 'Exploratory; DPT and reviewed labels had upstream access to expression.'
}, indent=2))
print('Saved exploratory artifacts to', out)

# %% [markdown]
# ## 12. Interpretation
#
# Read each species row and the paired short/long errors before making a claim. The following answers use the primary within-species held-out comparisons. Human sections share one donor; DPT and reviewed labels retain the upstream expression-dependence caveat.

# %%
def average_error(species, method, group='all'):
    selected = cv.species.eq(species) & cv.method.eq(method) & cv.length_group.eq(group)
    return float(cv.loc[selected, 'normalized_mse'].mean())

interpretation = []
species_support = []
for species in SPECIES_SAMPLES:
    point = average_error(species, 'point_continuous')
    extent = average_error(species, 'extent_morphology')
    dpt_error = average_error(species, 'DPT')
    equal = average_error(species, 'extent_equal')
    permuted = average_error(species, 'extent_permuted')
    short_delta = average_error(species, 'extent_morphology', 'short') - average_error(species, 'point_continuous', 'short')
    long_delta = average_error(species, 'extent_morphology', 'long') - average_error(species, 'point_continuous', 'long')
    curve_rho = repro.set_index(['species', 'method'])
    point_rho = float(curve_rho.loc[(species, 'point_continuous'), 'median_curve_rho'])
    extent_rho = float(curve_rho.loc[(species, 'extent_morphology'), 'median_curve_rho'])
    relative_gain = (point - extent) / point
    b = biology.set_index(['species', 'method'])
    point_order = float(b.loc[(species, 'point_continuous'), 'segment_order_rho'].mean())
    dpt_order = float(b.loc[(species, 'DPT'), 'segment_order_rho'].mean())
    supported = (relative_gain > .01 and long_delta < short_delta and
                 (min(equal, permuted) - extent) / point > .01 and extent_rho >= point_rho - .05)
    species_support.append(supported)
    interpretation.append(
        f'**{species.title()}:** held-out NMSE extent {extent:.3f}, point {point:.3f}, '
        f'DPT {dpt_error:.3f} (historical gene leakage). Extent change {relative_gain:+.2%}; '
        f'long minus point {long_delta:+.4f}, short minus point {short_delta:+.4f}. '
        f'Equal-width {equal:.3f}, shuffled-width {permuted:.3f}. '
        f'Within-species curve rho extent {extent_rho:.3f}, point {point_rho:.3f}. '
        f'Segment-order rho backbone {point_order:.3f}, DPT {dpt_order:.3f}. '
        f'The 1% descriptive width-specific screen is {"met" if supported else "not met"}.')
display(Markdown(
    '\n\n'.join(interpretation) + '\n\n' +
    '**Interpretation.** The paired heatmaps and marker plots show whether fitted curves follow '
    'expected PT programs. Marker panels and reviewed labels are contextual because they informed '
    'upstream annotation. Similar point and extent curves, or similar equal/permuted-width errors, '
    'leave physical width non-identifiable. Human section agreement is not donor replication. '
    f'Only {cross_edge_fraction:.1%} of joint kNN edges link species; shared-axis alignment remains provisional. '
    'Cross-species transfer on raw lognorm expression also measures abundance shifts, so it cannot '
    'by itself select a coordinate. ' +
    ('Both species pass the descriptive screen; external donors and direct length measurements are needed before a joint model.'
     if all(species_support) else
     'The combined evidence does not justify a full joint center/width model yet.')
))
