# %% [markdown]
# # 13 · Reviewed nephron and PT scFates pseudospace
#
# Reproduce notebook 03's reviewed two-pass Harmony path from the four original
# tubule-by-gene matrices and HCOP table: cluster the full cohort, retain the reviewed
# tubular nephron, reintegrate it, and then take the PT subset. No notebook 03 output is read.
#
# Compare nonbranching scFates and Scanpy DPT on both the full nephron and PT subset.
# The PT scFates curve supplies notebook 12's primary coordinate. The full-nephron
# nonbranching curve is a diagnostic forced through potentially branching biology.
# See https://scfates.readthedocs.io/en/latest/Basic_Curved_trajectory_analysis.html .
#
# Two human sections come from one healthy cortex donor; `HUK1_MED1` is a source name, not
# medullary tissue. The coordinates describe structures, not independently replicated donors.

# %% [markdown]
# ## 1 · Inputs and fixed choices
#
# Use `--data-root` / `--results-root` in the script mirror or the `PSEUDOSPACE_*` environment
# variables in a notebook. The cohort, reviewed labels, and clustering fingerprint are pinned
# to the healthy cohort. If clustering changes, review every cluster and update the fingerprint
# and labels together before making a coordinate.

# %%
import argparse
import hashlib
import os
import sys
from pathlib import Path

start = Path(__file__).resolve().parent if '__file__' in globals() else Path.cwd()
project = next(p for p in (start, *start.parents) if (p / 'pseudospace').is_dir())
sys.path.insert(0, str(project))
parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--data-root', type=Path)
parser.add_argument('--results-root', type=Path)
args, _ = parser.parse_known_args()
data_root = Path(args.data_root or os.environ.get('PSEUDOSPACE_DATA_ROOT') or project / 'data').expanduser()
results_root = Path(args.results_root or os.environ.get('PSEUDOSPACE_RESULTS_ROOT') or project / 'results').expanduser()
output = results_root / 'minimal_pt_scfates'
output.mkdir(parents=True, exist_ok=True)
samples = {'mouse': ('Ctrl1A2', 'Ctrl1A4'), 'human': ('HUK1_COR1', 'HUK1_MED1')}
source_paths = {name: data_root / 'tubule_by_gene' / f'{name}_tubule_by_gene_caleb.h5ad'
                for names in samples.values() for name in names}
ortholog_path = data_root / 'human_mouse_hcop_fifteen_column.txt.gz'
for path in (*source_paths.values(), ortholog_path):
    if not path.is_file():
        raise FileNotFoundError(path)

SEED = 0
NOTEBOOK_LOGIC_VERSION = '13-nephron-scfates-v2'
HARMONY_VERSION = '2.0.5'
SCFATES_NODES = 30
TRAJECTORY_DIMS = 5
PT_NEIGHBORS = 30
PASS2_LOGIC_VERSION = '13-nephron-pass2-v1'
REVIEWED_CLUSTER_LABELS = {
    # Confirmed against this partition's marker detection and the reviewed notebook 03
    # cluster memberships. The PT rows have Slc5a2 0.91 (S1), Slc22a6 0.94 (S2),
    # and Slc22a7/Slc7a13 0.74/0.78 (S3) detection, respectively.
    '0': 'PT-S1', '1': 'PT-S2', '2': 'CNT_CD', '3': 'AL', '4': 'DCT',
    '5': 'Glomerulus', '6': 'Unassigned', '7': 'PT-S3', '8': 'CCD',
    '9': 'AL', '10': 'SmoothMuscle',
}
REVIEWED_FINGERPRINT = {
    'n_structures': 26839, 'n_reporting_genes': 2000, 'resolution': 0.7,
    'n_neighbors': 30, 'random_state': 0, 'n_clusters': 11,
    'membership_sha1': '41fd77bb98da',
}

# %% [markdown]
# ## 2 · Rebuild and review the full-cohort clusters
#
# Keep notebook 03's input/QC/Harmony/Leiden settings that define its reviewed cluster IDs;
# then retain its nephron classes for the second integration. The accepted
# map creates target columns even for unmeasured human genes. As in notebook 03, keep those
# columns for the pinned clustering, but notebook 12 later excludes every unmeasured gene
# from analysis. Changing this gene set would change the Leiden partition.
#
# Mouse structures came from the quality-controlled fine GeoJSON exports: no extra mouse
# gene-count threshold. Human structures require 100 detected mapped genes. The lowest 5% of
# `n_spots` is removed separately in each species, among gene-rule survivors.

# %%
import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scFates as scf
from IPython.display import display
from scipy import sparse
from scipy.stats import spearmanr
import rpy2.robjects as ro

from pseudospace.cross_species import (read_ortholog_table, build_one_to_one_ortholog_map,
    load_cross_species_samples, combine_cross_species)
from pseudospace.harmony import run_harmony_rpy2, select_harmony_hvgs_by_condition
from pseudospace.io_qc import annotate_mito_ribo_mouse_symbols
from pseudospace.stage_cache import cached_anndata, digest
from pseudospace.vocabulary import KEEP_TUBULE_CLASSES, coarse_for

# Validation stays outside the cache. A changed R installation must never silently reuse it.
installed_harmony = str(ro.r('as.character(packageVersion("harmony"))')[0])
if installed_harmony != HARMONY_VERSION:
    raise RuntimeError(f'R harmony {HARMONY_VERSION} required; found {installed_harmony}.')
if scf.__version__ != '1.2.5':
    raise RuntimeError(f'scFates 1.2.5 required; found {scf.__version__}.')
ortholog_table = read_ortholog_table(ortholog_path)
ortholog_map = build_one_to_one_ortholog_map(ortholog_table, min_support=3)
ortholog_map.to_csv(output / 'ortholog_map_used.csv', index=False)
files = list(source_paths.values())

def build_reviewed_embedding():
    adatas, _ = load_cross_species_samples(
        data_root / 'tubule_by_gene', files, ortholog_map,
        mouse_samples=samples['mouse'], human_samples=samples['human'],
        human_regions={name: 'cortex' for name in samples['human']})
    obj = annotate_mito_ribo_mouse_symbols(
        combine_cross_species(adatas, require_measured_in_both=False))
    species = obj.obs.comparison_species.astype(str).to_numpy()
    detected = np.asarray((obj.X > 0).sum(axis=1)).ravel()
    gene_rule = (species == 'mouse') | (detected >= 100)
    spots = obj.obs.n_spots.to_numpy(float)
    keep = gene_rule.copy()
    for group in samples:
        group_mask = species == group
        floor = np.quantile(spots[group_mask & gene_rule], .05)
        keep &= ~group_mask | (spots >= floor)
    obj = obj[keep].copy()
    obj.layers['counts'] = obj.X.copy()
    detected_per_gene = np.asarray((obj.X > 0).sum(axis=0)).ravel()
    total_per_gene = np.asarray(obj.X.sum(axis=0)).ravel()
    keep_gene = ((detected_per_gene >= np.ceil(.05 * obj.n_obs))
                 & (total_per_gene >= 20))
    obj = obj[:, keep_gene].copy()
    sc.pp.normalize_total(obj, target_sum=1e4)
    sc.pp.log1p(obj)
    obj.layers['lognorm'] = obj.X.copy()
    obj = select_harmony_hvgs_by_condition(
        obj, group_key='comparison_species', groups=('mouse', 'human'),
        mode='intersection', min_mean=.0125, max_mean=3, min_disp=.5)
    hvg = obj[:, obj.var.highly_variable_for_harmony.to_numpy(bool)].copy()
    hvg.X = hvg.layers['lognorm'].copy()
    sc.tl.pca(hvg, n_comps=50, random_state=SEED)
    ro.r(f'set.seed({SEED})')
    hvg = run_harmony_rpy2(hvg, batch_key='sample', n_pcs=50, theta=6,
                           lambda_val=1, max_iter=30, tau=0)
    obj.obsm['X_harmony'] = hvg.obsm['X_harmony'].copy()
    sc.pp.highly_variable_genes(obj, n_top_genes=min(2000, obj.n_vars), flavor='seurat')
    sc.pp.neighbors(obj, n_neighbors=30, use_rep='X_harmony', random_state=SEED)
    sc.tl.leiden(obj, resolution=.7, key_added='leiden_coarse', flavor='igraph',
                 n_iterations=2, directed=False, random_state=SEED)
    sc.tl.umap(obj, random_state=SEED)
    return obj

cache_code = digest([NOTEBOOK_LOGIC_VERSION,
    Path(run_harmony_rpy2.__code__.co_filename),
    Path(combine_cross_species.__code__.co_filename),
    Path(annotate_mito_ribo_mouse_symbols.__code__.co_filename)])
cohort = cached_anndata('reviewed_embedding', build_reviewed_embedding,
    root=output / 'stage_cache',
    params={'seed': SEED, 'harmony': HARMONY_VERSION, 'scanpy': sc.__version__,
            'logic': NOTEBOOK_LOGIC_VERSION},
    inputs={'orthologs': digest(ortholog_path),
            'matrices': {name: digest(path) for name, path in source_paths.items()}},
    code=cache_code)
fingerprint = {'n_structures': cohort.n_obs,
    'n_reporting_genes': int(cohort.var.highly_variable.sum()),
    'resolution': .7, 'n_neighbors': 30, 'random_state': SEED,
    'n_clusters': cohort.obs.leiden_coarse.nunique(),
    'membership_sha1': hashlib.sha1(','.join(
        cohort.obs.leiden_coarse.astype(str)).encode()).hexdigest()[:12]}
display(pd.Series(fingerprint, name='this_run'))
review_genes = ('Slc5a2', 'Slc22a6', 'Slc22a7', 'Slc7a13', 'Umod',
                'Slc12a3', 'Aqp2', 'Nphs2', 'Acta2')
review_matrix = cohort[:, list(review_genes)].X
review_matrix = review_matrix.toarray() if sparse.issparse(review_matrix) else np.asarray(review_matrix)
cluster_ids = cohort.obs.leiden_coarse.astype(str).to_numpy()
cluster_review = pd.DataFrame([
    {'cluster': cluster, 'n_structures': int((cluster_ids == cluster).sum()),
     **{gene: float((review_matrix[cluster_ids == cluster, j] > 0).mean())
        for j, gene in enumerate(review_genes)}}
    for cluster in sorted(set(cluster_ids), key=int)])
cluster_review.to_csv(output / 'reviewed_cluster_marker_detection.csv', index=False)
display(cluster_review.round(2))
if fingerprint != REVIEWED_FINGERPRINT:
    raise ValueError('Reviewed cluster fingerprint changed. Reinspect each cluster in notebook 03 '
                     f'before updating labels here: {fingerprint}')
if set(cohort.obs.leiden_coarse.astype(str)) != set(REVIEWED_CLUSTER_LABELS):
    raise ValueError('Reviewed cluster map does not cover the current partition.')
cohort.obs['segment_class'] = cohort.obs.leiden_coarse.astype(str).map(REVIEWED_CLUSTER_LABELS)
cohort.obs['coarse_class'] = cohort.obs.segment_class.map(coarse_for)
cohort.obs['broad_tubule_marker_call'] = cohort.obs.coarse_class
display(cohort.obs.groupby(['sample', 'segment_class'], observed=True).size().unstack(fill_value=0))

# Centroids are input metadata, not inferred from the expression graph. Inspect PT spatial
# coverage in every specimen before relying on a trajectory from segmented structures.
fig, axes = plt.subplots(2, 2, figsize=(10, 10), layout='constrained')
for ax, name in zip(axes.flat, source_paths):
    rows = cohort.obs.loc[cohort.obs['sample'].astype(str).eq(name)]
    if not np.isfinite(rows[['x_centroid', 'y_centroid']].to_numpy(float)).all():
        raise ValueError(f'{name}: non-finite source centroids.')
    ax.scatter(rows.x_centroid, rows.y_centroid, s=1, color='0.85', rasterized=True)
    for label, color in (('PT-S1', '#4C9BD3'), ('PT-S2', '#5B8E55'), ('PT-S3', '#D6B48A')):
        selected = rows.loc[rows.segment_class.eq(label)]
        ax.scatter(selected.x_centroid, selected.y_centroid, s=2, alpha=.45,
                   color=color, label=label, rasterized=True)
    ax.set(title=name, xlabel='x centroid', ylabel='y centroid', aspect='equal')
axes[0, 0].legend(frameon=False, markerscale=3)
fig.savefig(output / 'pt_centroid_coverage.pdf', bbox_inches='tight')
plt.show()
plt.close(fig)

# %% [markdown]
# ## 3 · Retain reviewed nephron classes and run pass-2 Harmony
#
# The first-pass Leiden labels stay fixed. Remove glomerular, smooth-muscle and unresolved
# clusters, then repeat species-balanced HVG selection, PCA, Harmony, neighbors and UMAP on
# the retained tubular nephron. PT is selected only after this reintegration. No cluster
# has a reviewed DTL label in this cohort, so full-nephron coverage means all reviewed
# tubular structures, not a complete set of anatomical segments.

# %%
nephron_mask = cohort.obs.coarse_class.isin(KEEP_TUBULE_CLASSES).to_numpy()
filter_counts = (cohort.obs.assign(pass2_decision=np.where(nephron_mask, 'retained', 'removed'))
    .groupby(['coarse_class', 'pass2_decision'], observed=True).size()
    .rename('n_structures').reset_index())
filter_counts.to_csv(output / 'pass2_nephron_filter_counts.csv', index=False)
display(filter_counts)
if nephron_mask.sum() < 100 or cohort.obs.loc[nephron_mask, 'comparison_species'].nunique() != 2:
    raise ValueError('Pass-2 nephron needs at least 100 structures from both species.')

def build_nephron_embedding():
    obj = cohort[nephron_mask].copy()
    obj = select_harmony_hvgs_by_condition(
        obj, group_key='comparison_species', groups=('mouse', 'human'),
        mode='intersection', min_mean=.0125, max_mean=3, min_disp=.5)
    hvg = obj[:, obj.var.highly_variable_for_harmony.to_numpy(bool)].copy()
    hvg.X = hvg.layers['lognorm'].copy()
    sc.tl.pca(hvg, n_comps=50, random_state=SEED)
    ro.r(f'set.seed({SEED})')
    hvg = run_harmony_rpy2(hvg, batch_key='sample', n_pcs=50, theta=6,
                           lambda_val=1, max_iter=30, tau=0)
    obj.obsm['X_harmony'] = hvg.obsm['X_harmony'].copy()
    sc.pp.neighbors(obj, n_neighbors=30, use_rep='X_harmony', random_state=SEED)
    sc.tl.umap(obj, random_state=SEED)
    return obj

nephron = cached_anndata('pass2_nephron_embedding', build_nephron_embedding,
    root=output / 'stage_cache',
    params={'seed': SEED, 'harmony': HARMONY_VERSION, 'scanpy': sc.__version__,
            'logic': PASS2_LOGIC_VERSION},
    inputs={'pass1_embedding': digest(np.asarray(cohort.obsm['X_harmony'])),
            'retained_ids': digest(cohort.obs_names[nephron_mask].tolist()),
            'reviewed_labels': digest(cohort.obs.segment_class.astype(str).tolist())},
    code=digest([PASS2_LOGIC_VERSION, Path(run_harmony_rpy2.__code__.co_filename)]))
if not nephron.obs_names.equals(cohort.obs_names[nephron_mask]):
    raise ValueError('Pass-2 cache does not match the reviewed nephron membership.')
pt = nephron[nephron.obs.coarse_class.eq('PT')].copy()
if set(pt.obs.segment_class) != {'PT-S1', 'PT-S2', 'PT-S3'}:
    raise ValueError('The reviewed PT subset must contain S1, S2, and S3.')
pt.obs.groupby(['comparison_species', 'sample', 'segment_class'], observed=True).size().rename(
    'structures').to_csv(output / 'pt_cluster_counts.csv')
display(pt.obs.groupby(['sample', 'segment_class'], observed=True).size().unstack(fill_value=0))

fig, axes = plt.subplots(1, 3, figsize=(15, 4), layout='constrained')
for ax, key, title in zip(axes, ('sample', 'coarse_class', 'segment_class'),
                          ('Specimen', 'Nephron family', 'Reviewed segment')):
    codes = nephron.obs[key].astype('category')
    for category in codes.cat.categories:
        mask = codes.eq(category).to_numpy()
        ax.scatter(*nephron.obsm['X_umap'][mask, :2].T, s=1.5, alpha=.55,
                   label=str(category), rasterized=True)
    ax.set(title=title, xlabel='UMAP 1', ylabel='UMAP 2')
    ax.legend(frameon=False, markerscale=3, fontsize=7, ncol=2)
fig.savefig(output / 'pass2_nephron_umap.pdf', bbox_inches='tight')
plt.show()
plt.close(fig)

# %% [markdown]
# ## 4 · Compare nonbranching scFates and DPT on nephron and PT
#
# Each comparison fits both methods to the same first five pass-2 Harmony dimensions.
# The marker panels orient the paths after cluster review. DPT starts near the scFates
# early root. The full-nephron single curve is a topology diagnostic, not a claim that
# every nephron segment lies along one anatomical path. The PT scFates curve is the
# coordinate saved for notebook 12.

# %%
early_genes = ('Slc5a2', 'Slc5a12', 'Gatm')
late_genes = ('Slc22a7', 'Slc7a13', 'Cyp7b1')

def fit_coordinates(obj, name, n_neighbors):
    axis_genes = [g for g in (*early_genes, *late_genes) if g in obj.var_names]
    for panel in (early_genes, late_genes):
        if len(set(panel) & set(axis_genes)) < 2:
            raise ValueError(f'{name}: too few measured axis markers in {panel}.')
    values = obj[:, axis_genes].X
    values = values.toarray() if sparse.issparse(values) else np.asarray(values)
    z = np.empty_like(values, dtype=float)
    for species in samples:
        mask = obj.obs.comparison_species.astype(str).eq(species).to_numpy()
        block = values[mask]
        z[mask] = (block - block.mean(axis=0)) / np.maximum(block.std(axis=0), 1e-6)
    obj.obs['pt_marker_axis'] = (
        z[:, [axis_genes.index(g) for g in late_genes if g in axis_genes]].mean(axis=1)
        - z[:, [axis_genes.index(g) for g in early_genes if g in axis_genes]].mean(axis=1))
    obj.obs['early_tip_score'] = -obj.obs.pt_marker_axis
    obj.obsm['X_trajectory'] = np.asarray(obj.obsm['X_harmony'][:, :TRAJECTORY_DIMS], dtype=float)
    scf.tl.curve(obj, Nodes=SCFATES_NODES, use_rep='X_trajectory',
                 ndims_rep=TRAJECTORY_DIMS, seed=SEED)
    graph = obj.uns['graph']
    degree = np.asarray(graph['B']).astype(bool).sum(axis=0)
    if not (len(graph['tips']) == 2 and len(graph['forks']) == 0
            and np.count_nonzero(degree == 1) == 2 and np.all(degree <= 2)):
        raise ValueError(f'{name}: scFates graph is not one nonbranching path.')
    scf.tl.root(obj, 'early_tip_score', tips_only=True)
    root_node = int(graph['root'])
    membership = np.asarray(obj.obsm['X_R'][:, root_node]).ravel()
    near_root = np.flatnonzero(membership >= np.quantile(membership, .99))
    root_cell = near_root[np.argmax(obj.obs.early_tip_score.to_numpy()[near_root])]
    # One seeded mapping keeps the producer compact; more maps would quantify projection uncertainty.
    scf.tl.pseudotime(obj, n_jobs=2, n_map=1, seed=SEED)
    scfates_time = obj.obs.t.to_numpy(float)
    if not np.isfinite(scfates_time).all() or np.ptp(scfates_time) <= 0:
        raise ValueError(f'{name}: scFates gave non-finite or constant pseudotime.')
    obj.obs['shared_pseudospace'] = (scfates_time - scfates_time.min()) / np.ptp(scfates_time)

    sc.pp.neighbors(obj, n_neighbors=n_neighbors, use_rep='X_trajectory',
                    random_state=SEED, key_added='trajectory_neighbors')
    sc.tl.diffmap(obj, neighbors_key='trajectory_neighbors', random_state=SEED)
    obj.uns['iroot'] = int(root_cell)
    sc.tl.dpt(obj, neighbors_key='trajectory_neighbors')
    dpt_time = obj.obs.dpt_pseudotime.to_numpy(float)
    finite = np.isfinite(dpt_time)
    if finite.sum() < .95 * obj.n_obs or np.ptp(dpt_time[finite]) <= 0:
        raise ValueError(f'{name}: DPT has too many disconnected or constant structures.')
    normalized = np.full(obj.n_obs, np.nan)
    normalized[finite] = (dpt_time[finite] - dpt_time[finite].min()) / np.ptp(dpt_time[finite])
    obj.obs['dpt_pseudospace'] = normalized
    print(f'{name}: {obj.n_obs:,} structures; {int((~finite).sum()):,} DPT-disconnected; '
          f'scFates root node {root_node}, DPT root {obj.obs_names[root_cell]}')
    return graph, root_node, root_cell

nephron_graph, nephron_root_node, nephron_root_cell = fit_coordinates(
    nephron, 'nephron', max(PT_NEIGHBORS, int(np.sqrt(nephron.n_obs))))
pt_graph, pt_root_node, pt_root_cell = fit_coordinates(pt, 'PT', PT_NEIGHBORS)
if not np.isfinite(pt.obs.dpt_pseudospace).all():
    raise ValueError('PT DPT must cover every saved structure.')

def compare_methods(obj, name):
    rows = []
    for specimen, obs in obj.obs.groupby('sample', observed=True):
        finite = obs.dpt_pseudospace.notna()
        valid = obs.loc[finite]
        rows.append({'subset': name, 'specimen': specimen,
            'species': obs.comparison_species.iloc[0], 'n_structures': len(obs),
            'n_dpt_connected': int(finite.sum()),
            'spearman_scfates_dpt': spearmanr(valid.shared_pseudospace, valid.dpt_pseudospace).statistic,
            'spearman_scfates_marker_axis': spearmanr(obs.shared_pseudospace, obs.pt_marker_axis).statistic,
            'spearman_dpt_marker_axis': spearmanr(valid.dpt_pseudospace, valid.pt_marker_axis).statistic})
    obj.obs[['comparison_species', 'sample', 'segment_class', 'coarse_class',
             'pt_marker_axis', 'shared_pseudospace', 'dpt_pseudospace']].to_csv(
        output / f'{name}_coordinate_comparison.csv')
    return pd.DataFrame(rows)

nephron_comparison = compare_methods(nephron, 'nephron')
comparison = compare_methods(pt, 'pt')
pd.concat([nephron_comparison, comparison]).to_csv(
    output / 'coordinate_agreement_by_specimen_and_subset.csv', index=False)
display(pd.concat([nephron_comparison, comparison]).round(3))
segment_order = ('PT-S1', 'PT-S2', 'PT-S3')
order_rows = []
for name, obs in pt.obs.groupby('sample', observed=True):
    medians = obs.groupby('segment_class', observed=True).shared_pseudospace.median().reindex(segment_order)
    order_rows.append({'specimen': name, **medians.to_dict(),
                       'S2_minus_S1': medians['PT-S2'] - medians['PT-S1'],
                       'S3_minus_S2': medians['PT-S3'] - medians['PT-S2']})
    if medians.isna().any() or medians['PT-S3'] <= max(medians['PT-S1'], medians['PT-S2']):
        raise ValueError(f'{name}: scFates did not place S3 after S1/S2: {medians.to_dict()}')
order_table = pd.DataFrame(order_rows)
order_table.to_csv(output / 'pt_segment_order_by_specimen.csv', index=False)
display(order_table.round(3))
comparison.to_csv(output / 'coordinate_agreement_by_specimen.csv', index=False)

# %% [markdown]
# ## 5 · Inspect full-nephron and PT paths and method agreement
#
# The full-nephron panels show whether a forced single scFates curve preserves the reviewed
# families, and whether it agrees with DPT within each specimen. The PT panels repeat those
# checks on the subset that notebook 12 actually uses. Discordance is a coordinate diagnostic.

# %%
fig, axes = plt.subplots(1, 3, figsize=(15, 4), layout='constrained')
xy = nephron.obsm['X_umap'][:, :2]
for family in KEEP_TUBULE_CLASSES:
    mask = nephron.obs.coarse_class.eq(family).to_numpy()
    if mask.any():
        axes[0].scatter(xy[mask, 0], xy[mask, 1], s=2, alpha=.4,
                        label=family, rasterized=True)
axes[0].legend(frameon=False, markerscale=4)
axes[0].set(title='Reviewed nephron families', xlabel='UMAP 1', ylabel='UMAP 2')
for ax, column, title in zip(axes[1:], ('shared_pseudospace', 'dpt_pseudospace'),
                             ('scFates · full nephron', 'DPT · full nephron')):
    points = ax.scatter(xy[:, 0], xy[:, 1], c=nephron.obs[column], s=2,
                        cmap='viridis', rasterized=True)
    fig.colorbar(points, ax=ax, label='Pseudospace')
    ax.set(title=title, xlabel='UMAP 1', ylabel='UMAP 2')
fig.savefig(output / 'nephron_scfates_dpt_umap.pdf', bbox_inches='tight')
plt.show()
plt.close(fig)

# %%
fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True, sharey=True, layout='constrained')
for ax, (specimen, obs) in zip(axes.flat, nephron.obs.groupby('sample', observed=True)):
    for family in KEEP_TUBULE_CLASSES:
        rows = obs.loc[obs.coarse_class.eq(family)]
        ax.scatter(rows.shared_pseudospace, rows.dpt_pseudospace, s=3, alpha=.2,
                   label=family)
    rho = nephron_comparison.set_index('specimen').loc[specimen, 'spearman_scfates_dpt']
    ax.plot([0, 1], [0, 1], color='0.3', lw=.8, ls='--')
    ax.set(title=f'{specimen} · Spearman {rho:.2f}', xlim=(0, 1), ylim=(0, 1),
           xlabel='scFates pseudospace', ylabel='DPT pseudospace')
axes[0, 0].legend(frameon=False, markerscale=3)
fig.savefig(output / 'nephron_scfates_vs_dpt_by_specimen.pdf', bbox_inches='tight')
plt.show()
plt.close(fig)

# %%
from pseudospace.vocabulary import SEGMENT_DISPLAY_ORDER

reviewed_segments = [segment for segment in SEGMENT_DISPLAY_ORDER
                     if segment in set(nephron.obs.segment_class.astype(str))]
fig, axes = plt.subplots(2, 2, figsize=(16, 9), sharey=True, layout='constrained')
for row, (column, title) in enumerate((('shared_pseudospace', 'scFates'),
                                        ('dpt_pseudospace', 'DPT'))):
    for col, species in enumerate(samples):
        obs = nephron.obs.loc[nephron.obs.comparison_species.eq(species)]
        values = [obs.loc[obs.segment_class.astype(str).eq(segment), column]
                  .dropna().to_numpy(float) for segment in reviewed_segments]
        axes[row, col].boxplot([v if len(v) else np.array([np.nan]) for v in values],
                               showfliers=False)
        axes[row, col].set_xticks(range(1, len(reviewed_segments) + 1))
        axes[row, col].set_xticklabels([f'{segment}\n(n={len(v):,})'
                                       for segment, v in zip(reviewed_segments, values)])
        plt.setp(axes[row, col].get_xticklabels(), rotation=45, ha='right')
        axes[row, col].set(title=f'{title} · {species.title()} (n={len(obs):,})',
                           ylim=(0, 1), ylabel='Pseudospace' if col == 0 else '')
fig.suptitle('Full-nephron pseudospace by reviewed segment', y=1.02)
fig.savefig(output / 'nephron_segment_order_scfates_dpt.pdf', bbox_inches='tight')
fig.savefig(output / 'nephron_segment_order_scfates_dpt.png', dpi=150, bbox_inches='tight')
plt.show()
plt.close(fig)

# %%
colors = {'PT-S1': '#4C9BD3', 'PT-S2': '#5B8E55', 'PT-S3': '#D6B48A'}
fig, axes = plt.subplots(1, 3, figsize=(14, 4), layout='constrained')
xy = pt.obsm['X_trajectory'][:, :2]
for label in segment_order:
    mask = pt.obs.segment_class.eq(label).to_numpy()
    axes[0].scatter(xy[mask, 0], xy[mask, 1], s=2, alpha=.25, color=colors[label], label=label)
nodes = np.asarray(pt_graph['F'])[:2].T
for a, b in np.argwhere(np.triu(np.asarray(pt_graph['B']) > 0)):
    axes[0].plot(nodes[[a, b], 0], nodes[[a, b], 1], color='black', lw=1)
axes[0].scatter(nodes[:, 0], nodes[:, 1], s=9, color='black')
axes[0].legend(frameon=False, markerscale=4)
axes[0].set(xlabel='Harmony dimension 1', ylabel='Harmony dimension 2',
            title='Reviewed PT clusters and scFates curve')
for ax, column, title in zip(axes[1:], ('shared_pseudospace', 'dpt_pseudospace'),
                             ('scFates (primary)', 'DPT (comparator)')):
    points = ax.scatter(xy[:, 0], xy[:, 1], c=pt.obs[column], s=2, cmap='viridis', rasterized=True)
    fig.colorbar(points, ax=ax, label='Pseudospace')
    ax.set(xlabel='Harmony dimension 1', ylabel='Harmony dimension 2', title=title)
fig.savefig(output / 'pt_scfates_dpt_embedding.pdf', bbox_inches='tight')
plt.show()
plt.close(fig)

# %%
fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True, sharey=True, layout='constrained')
for ax, (name, obs) in zip(axes.flat, pt.obs.groupby('sample', observed=True)):
    for label in segment_order:
        rows = obs.loc[obs.segment_class.eq(label)]
        ax.scatter(rows.shared_pseudospace, rows.dpt_pseudospace, s=3, alpha=.18,
                   color=colors[label], label=label)
    ax.plot([0, 1], [0, 1], color='0.3', lw=.8, ls='--')
    rho = comparison.set_index('specimen').loc[name, 'spearman_scfates_dpt']
    ax.set(title=f'{name} · Spearman {rho:.2f}', xlim=(0, 1), ylim=(0, 1),
           xlabel='scFates pseudospace', ylabel='DPT pseudospace')
axes[0, 0].legend(frameon=False, markerscale=3)
fig.savefig(output / 'pt_scfates_vs_dpt_by_specimen.pdf', bbox_inches='tight')
plt.show()
plt.close(fig)

# %%
fig, axes = plt.subplots(2, 2, figsize=(16, 9), sharey=True, layout='constrained')
for row, (column, title) in enumerate((('shared_pseudospace', 'scFates'),
                                        ('dpt_pseudospace', 'DPT'))):
    for col, species in enumerate(samples):
        obs = pt.obs.loc[pt.obs.comparison_species.eq(species)]
        values = [obs.loc[obs.segment_class.eq(label), column].dropna().to_numpy(float)
                  for label in segment_order]
        axes[row, col].boxplot(values, showfliers=False)
        axes[row, col].set_xticks(range(1, len(segment_order) + 1))
        axes[row, col].set_xticklabels([f'{label}\n(n={len(v):,})'
                                       for label, v in zip(segment_order, values)])
        plt.setp(axes[row, col].get_xticklabels(), rotation=45, ha='right')
        axes[row, col].set(title=f'{title} · {species.title()} (n={len(obs):,})',
                           ylim=(0, 1), ylabel='Pseudospace' if col == 0 else '')
fig.suptitle('PT-specific pseudospace by reviewed segment', y=1.02)
fig.savefig(output / 'pt_segment_order_scfates_dpt.pdf', bbox_inches='tight')
fig.savefig(output / 'pt_segment_order_scfates_dpt.png', dpi=150, bbox_inches='tight')
plt.show()
plt.close(fig)

# %% [markdown]
# ### Marker-gradient heatmaps for both coordinates
#
# Use notebook 03's fine reference panels as visual checks, without changing the reviewed
# labels. Each figure holds mouse and human rows with scFates and DPT columns in the same
# gene order and colour range. As in notebook 03, each panel trims 5% from both coordinate
# tails, bins 120 positions, smooths across bins, and z-scores each gene along its own
# coordinate. Interpolated gaps and within-gene scaling make these diagnostic plots, not
# measurements of absolute expression or independent significance tests.

# %%
from pseudospace.heatmaps import binned_interpolated_expression

NEPHRON_HEATMAP_MARKERS = {
    'PT-S1': ('Lrp2', 'Cubn', 'Slc34a1', 'Slc5a2', 'Slc5a12', 'Gatm'),
    'PT-S2': ('Slc22a6', 'Slc13a3', 'Cyp2e1'),
    'PT-S3': ('Slc22a7', 'Slc7a13', 'Cyp7b1', 'Slc6a18', 'Acsm3'),
    'DTL': ('Aqp1', 'Slc14a2', 'Corin', 'Fst'),
    'ATL': ('Clcnka', 'Sptssb', 'Akr1b3'),
    'TAL': ('Slc12a1', 'Umod', 'Kcnj1', 'Cldn16'),
    'DCT': ('Slc12a3', 'Pvalb', 'Trpm6', 'Egf'),
    'CNT': ('Calb1', 'Hsd11b2', 'Slc8a1'),
    'CCD': ('Aqp2', 'Aqp3', 'Fxyd4'),
    'OMCD': ('Atp6v0d2', 'Rhcg', 'Foxi1'),
    'IMCD': ('Aqp4', 'Slc14a2', 'Wnt7b'),
}
PT_HEATMAP_MARKERS = {
    'PT-S1': ('Slc5a2', 'Slc5a12', 'Gatm'),
    'PT-S2': ('Slc22a6', 'Slc13a3', 'Cyp2e1'),
    'PT-S3': ('Slc22a7', 'Cyp7b1'),
}

def compare_marker_heatmaps(obj, marker_groups, scope):
    measured = obj.var.measured_in_both_inputs.astype(bool)
    panels = {group: [gene for gene in genes if gene in obj.var_names and measured.loc[gene]]
              for group, genes in marker_groups.items()}
    coverage = pd.DataFrame([{'program': group, 'requested_genes': len(genes),
                              'measured_genes': len(panels[group]),
                              'genes_used': '; '.join(panels[group])}
                             for group, genes in marker_groups.items()])
    coverage.to_csv(output / f'{scope}_marker_heatmap_coverage.csv', index=False)
    if any(not genes for genes in panels.values()):
        raise ValueError(f'{scope}: a reference program has no jointly measured genes.')
    rows = [(group, gene) for group, genes in panels.items() for gene in genes]
    unique_genes = list(dict.fromkeys(gene for _, gene in rows))
    gene_idx = [unique_genes.index(gene) for _, gene in rows]
    labels = [f'{group} · {gene}' if gene == panels[group][0] else gene
              for group, gene in rows]
    fig, axes = plt.subplots(2, 2, figsize=(15, max(6, .34 * len(rows) + 2)),
                             sharex=True, sharey=True, layout='constrained')
    image = None
    for row, species in enumerate(samples):
        species_mask = obj.obs.comparison_species.astype(str).eq(species).to_numpy()
        for col, (column, method) in enumerate((('shared_pseudospace', 'scFates'),
                                                 ('dpt_pseudospace', 'DPT'))):
            subset = obj[species_mask & np.isfinite(obj.obs[column].to_numpy(float)),
                         unique_genes].copy()
            order = np.argsort(subset.obs[column].to_numpy(float), kind='stable')
            n_trim = int(np.floor(.05 * len(order)))
            subset = subset[order[n_trim:len(order) - n_trim]].copy()
            _, heatmap_z, _, edges, _ = binned_interpolated_expression(
                subset, gene_idx, column, n_bins=120, smooth_sigma=2.5)
            if heatmap_z.shape[0] != len(rows) or not np.isfinite(heatmap_z).all():
                raise ValueError(f'{scope} {species} {method}: invalid marker heatmap.')
            ax = axes[row, col]
            image = ax.imshow(np.clip(heatmap_z, -2, 2), aspect='auto', cmap='bwr',
                              vmin=-2, vmax=2, interpolation='nearest',
                              extent=(edges[0], edges[-1], len(rows) - .5, -.5))
            ax.set(xlim=(0, 1), title=f'{species.title()} · {method} (n={subset.n_obs:,})',
                   xlabel='Pseudospace' if row == 1 else '')
            ax.set_yticks(np.arange(len(rows)))
            if col == 0:
                ax.set_yticklabels(labels, fontsize=7 if scope == 'nephron' else 9)
            else:
                ax.tick_params(labelleft=False)
            boundary = 0
            for genes in panels.values():
                boundary += len(genes)
                if boundary < len(rows):
                    ax.axhline(boundary - .5, color='black', lw=.5)
    fig.colorbar(image, ax=axes, label='Within-gene z-score', shrink=.65)
    fig.suptitle(f'{scope.upper()} marker gradients · scFates versus DPT', fontsize=15)
    fig.savefig(output / f'{scope}_marker_heatmap_scfates_dpt.pdf', bbox_inches='tight')
    fig.savefig(output / f'{scope}_marker_heatmap_scfates_dpt.png', dpi=180,
                bbox_inches='tight')
    plt.show()
    plt.close(fig)
    return coverage

display(compare_marker_heatmaps(nephron, NEPHRON_HEATMAP_MARKERS, 'nephron'))
display(compare_marker_heatmaps(pt, PT_HEATMAP_MARKERS, 'pt'))

# %% [markdown]
# ### Current coordinate check
#
# Inspect `coordinate_agreement_by_specimen_and_subset.csv`, the reviewed-segment plots,
# and the paired marker heatmaps together. A high pooled correlation can hide reversed or flat
# ordering in one specimen. Two mouse specimens and two cortex sections of one human donor
# give structure-level diagnostics, not independent donor replication. The producer refits
# clustering and pass-2 Harmony, so changes from notebook 03 are not coordinate-only effects.

# %% [markdown]
# ## 6 · Save the nephron diagnostic and notebook 12 PT contract
#
# Save metadata-only coordinates for both scopes. Notebook 12 reads the PT file and reloads
# expression from raw counts. Both human slices retain the `cortex` annotation.

# %%
contract = ('comparison_species', 'sample', 'region', 'segment_class', 'coarse_class',
            'broad_tubule_marker_call', 'shared_pseudospace', 'dpt_pseudospace')
for name, obj, root_node, root_cell in (
    ('nephron', nephron, nephron_root_node, nephron_root_cell),
    ('pt', pt, pt_root_node, pt_root_cell)):
    saved = ad.AnnData(obs=obj.obs.loc[:, contract].copy())
    saved.uns['coordinate_method'] = (
        'scFates 1.2.5 tl.curve; 30 nodes; 5 pass-2 Harmony dimensions; no forks')
    saved.uns['comparison_method'] = 'Scanpy DPT; same pass-2 embedding and early endpoint'
    saved.uns['reviewed_cluster_fingerprint'] = fingerprint
    saved.uns['scfates_root_node'] = root_node
    saved.uns['dpt_root_structure'] = str(obj.obs_names[root_cell])
    saved.write_h5ad(output / f'cross_species_{name}_scfates.h5ad')
    assert saved.obs_names.equals(obj.obs_names)
    assert saved.obs.shared_pseudospace.between(0, 1).all()
print('Notebook 12 inputs:', output / 'cross_species_pt_scfates.h5ad',
      output / 'ortholog_map_used.csv', sep='\n  ')
