# %% [markdown]
# # 13 · Minimal PT scFates pseudospace
#
# Build only the two files notebook 12 needs: a reviewed PT structure table with a shared
# coordinate, and the accepted human–mouse ortholog map. The input is the four original
# tubule-by-gene matrices plus the HCOP table. No notebook 03 output is read.
#
# **Primary coordinate:** scFates' nonbranching principal curve on the shared Harmony PT
# embedding. **Comparator:** Scanpy DPT on the same PT embedding and rooted at the same end.
# scFates' curve is a constrained single path, not a principal tree or branch analysis.
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
NOTEBOOK_LOGIC_VERSION = '13-minimal-pt-scfates-v1'
HARMONY_VERSION = '2.0.5'
SCFATES_NODES = 30
PT_DIMS = 5
PT_NEIGHBORS = 30
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
# ## 2 · Rebuild the reviewed PT subset
#
# Keep notebook 03's input/QC/Harmony/Leiden settings that define its reviewed cluster IDs;
# omit its global-nephron trajectory and all downstream expression analyses. The accepted
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
cohort.obs['coarse_class'] = np.where(cohort.obs.segment_class.str.startswith('PT-'), 'PT',
    cohort.obs.segment_class)
cohort.obs['broad_tubule_marker_call'] = cohort.obs.coarse_class
pt = cohort[cohort.obs.coarse_class.eq('PT')].copy()
if set(pt.obs.segment_class) != {'PT-S1', 'PT-S2', 'PT-S3'}:
    raise ValueError('The reviewed PT subset must contain S1, S2, and S3.')
pt.obs.groupby(['comparison_species', 'sample', 'segment_class'], observed=True).size().rename(
    'structures').to_csv(output / 'pt_cluster_counts.csv')
display(pt.obs.groupby(['sample', 'segment_class'], observed=True).size().unstack(fill_value=0))

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
# ## 3 · Fit one nonbranching scFates curve and a DPT comparator
#
# Fit both methods in the first five Harmony dimensions of the reviewed PT subset. Predeclared
# early and late marker panels orient the path; they do not assign clusters or pick pathways.
# scFates chooses the early endpoint of its **single curve**. DPT starts at a structure near
# that same endpoint. Every saved coordinate is scaled to [0, 1] on all PT structures before
# notebook 12 takes its four-specimen common-support interval.

# %%
early_genes = ('Slc5a2', 'Slc5a12', 'Gatm')
late_genes = ('Slc22a7', 'Slc7a13', 'Cyp7b1')
for panel in (early_genes, late_genes):
    if len(set(panel) & set(pt.var_names)) < 2:
        raise ValueError(f'Too few measured axis markers in {panel}.')
axis_genes = [g for g in (*early_genes, *late_genes) if g in pt.var_names]
axis_values = pt[:, axis_genes].X
axis_values = axis_values.toarray() if sparse.issparse(axis_values) else np.asarray(axis_values)
axis_z = np.empty_like(axis_values, dtype=float)
for group in samples:
    mask = pt.obs.comparison_species.astype(str).eq(group).to_numpy()
    block = axis_values[mask]
    axis_z[mask] = (block - block.mean(axis=0)) / np.maximum(block.std(axis=0), 1e-6)
marker_axis = axis_z[:, [axis_genes.index(g) for g in late_genes if g in axis_genes]].mean(axis=1) - axis_z[:, [axis_genes.index(g) for g in early_genes if g in axis_genes]].mean(axis=1)
pt.obs['pt_marker_axis'] = marker_axis
pt.obs['early_tip_score'] = -marker_axis
pt.obsm['X_pt'] = np.asarray(pt.obsm['X_harmony'][:, :PT_DIMS], dtype=float)
scf.tl.curve(pt, Nodes=SCFATES_NODES, use_rep='X_pt', ndims_rep=PT_DIMS, seed=SEED)
graph = pt.uns['graph']
degree = np.asarray(graph['B']).astype(bool).sum(axis=0)
if not (len(graph['tips']) == 2 and len(graph['forks']) == 0
        and np.count_nonzero(degree == 1) == 2 and np.all(degree <= 2)):
    raise ValueError('scFates principal graph is not a single nonbranching path.')
scf.tl.root(pt, 'early_tip_score', tips_only=True)
root_node = int(pt.uns['graph']['root'])
membership = np.asarray(pt.obsm['X_R'][:, root_node]).ravel()
near_root = np.flatnonzero(membership >= np.quantile(membership, .99))
root_cell = near_root[np.argmax(pt.obs.early_tip_score.to_numpy()[near_root])]
# One seeded mapping keeps this producer small; increase n_map to assess projection uncertainty.
scf.tl.pseudotime(pt, n_jobs=2, n_map=1, seed=SEED)
scfates_time = pt.obs.t.to_numpy(float)
if not np.isfinite(scfates_time).all() or np.ptp(scfates_time) <= 0:
    raise ValueError('scFates did not give every PT structure finite pseudotime.')
pt.obs['shared_pseudospace'] = (scfates_time - scfates_time.min()) / np.ptp(scfates_time)

sc.pp.neighbors(pt, n_neighbors=PT_NEIGHBORS, use_rep='X_pt', random_state=SEED,
                key_added='pt_neighbors')
sc.tl.diffmap(pt, neighbors_key='pt_neighbors', random_state=SEED)
pt.uns['iroot'] = int(root_cell)
sc.tl.dpt(pt, neighbors_key='pt_neighbors')
dpt_time = pt.obs.dpt_pseudotime.to_numpy(float)
if not np.isfinite(dpt_time).all() or np.ptp(dpt_time) <= 0:
    raise ValueError('DPT is disconnected from the shared PT root.')
pt.obs['dpt_pseudospace'] = (dpt_time - dpt_time.min()) / np.ptp(dpt_time)
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
pt.obs[['comparison_species', 'sample', 'segment_class', 'pt_marker_axis',
        'shared_pseudospace', 'dpt_pseudospace']].to_csv(output / 'pt_coordinate_comparison.csv')
comparison = pd.DataFrame([{'specimen': name, 'species': obs.comparison_species.iloc[0],
    'n_structures': len(obs),
    'spearman_scfates_dpt': spearmanr(obs.shared_pseudospace, obs.dpt_pseudospace).statistic,
    'spearman_scfates_marker_axis': spearmanr(obs.shared_pseudospace, obs.pt_marker_axis).statistic,
    'spearman_dpt_marker_axis': spearmanr(obs.dpt_pseudospace, obs.pt_marker_axis).statistic}
    for name, obs in pt.obs.groupby('sample', observed=True)])
comparison.to_csv(output / 'coordinate_agreement_by_specimen.csv', index=False)
display(comparison.round(3))

# %% [markdown]
# ## 4 · Inspect the curve, ordering, and method agreement
#
# The embedding uses the same five dimensions for scFates and DPT. A high overall correlation
# can hide a reversed or flattened segment in one specimen; inspect each specimen and the
# S1/S2/S3 distributions before interpreting notebook 12. Disagreement is a coordinate
# sensitivity signal, not an extra test of species biology.

# %%
colors = {'PT-S1': '#4C9BD3', 'PT-S2': '#5B8E55', 'PT-S3': '#D6B48A'}
fig, axes = plt.subplots(1, 3, figsize=(14, 4), layout='constrained')
xy = pt.obsm['X_pt'][:, :2]
for label in segment_order:
    mask = pt.obs.segment_class.eq(label).to_numpy()
    axes[0].scatter(xy[mask, 0], xy[mask, 1], s=2, alpha=.25, color=colors[label], label=label)
nodes = np.asarray(graph['F'])[:2].T
for a, b in np.argwhere(np.triu(np.asarray(graph['B']) > 0)):
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
fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True, layout='constrained')
for row, (column, title) in enumerate((('shared_pseudospace', 'scFates'),
                                        ('dpt_pseudospace', 'DPT'))):
    for col, species in enumerate(samples):
        obs = pt.obs.loc[pt.obs.comparison_species.eq(species)]
        values = [obs.loc[obs.segment_class.eq(label), column].to_numpy() for label in segment_order]
        box = axes[row, col].boxplot(values, tick_labels=segment_order, showfliers=False,
                                      patch_artist=True)
        for patch, label in zip(box['boxes'], segment_order):
            patch.set_facecolor(colors[label])
        axes[row, col].set(title=f'{title} · {species}', ylim=(0, 1), ylabel='Pseudospace')
fig.savefig(output / 'pt_segment_order_scfates_dpt.pdf', bbox_inches='tight')
plt.show()
plt.close(fig)

# %% [markdown]
# ### Current coordinate check
#
# In the current four-specimen run, scFates places the PT-S3 median after both PT-S1 and PT-S2 in every specimen. Human PT-S1 and PT-S2 overlap: their median difference is +0.006 in `HUK1_COR1` and −0.001 in `HUK1_MED1`. The S1/S2 cluster labels therefore should not be read as sharply separated positions in human cortex.
#
# Within-specimen scFates–DPT Spearman correlations are 0.49–0.72. The two coordinates are related but not interchangeable; notebook 12 uses **scFates only** as its primary coordinate. These are structure-level diagnostics from two mouse specimens and two cortex sections of one human donor, not biological replication. The producer refits Harmony and clustering, so its PT membership can differ from notebook 03; a comparison with the old notebook 12 run is not a coordinate-only sensitivity analysis. Inspect the saved per-specimen tables and plots if the input data or integration changes.

# %% [markdown]
# ## 5 · Save the small notebook 12 contract
#
# Notebook 12 reloads raw counts itself. Save metadata only: specimen identity, reviewed PT
# label, the primary scFates coordinate, and DPT for sensitivity inspection. Both human slices
# retain the `cortex` annotation. The ortholog CSV was written before embedding, with the same
# accepted pairs used to build the cross-species object.

# %%
contract = ('comparison_species', 'sample', 'region', 'segment_class', 'coarse_class',
            'broad_tubule_marker_call', 'shared_pseudospace', 'dpt_pseudospace')
saved = ad.AnnData(obs=pt.obs.loc[:, contract].copy())
saved.uns['coordinate_method'] = 'scFates 1.2.5 tl.curve; 30 nodes; 5 Harmony dimensions; no forks'
saved.uns['comparison_method'] = 'Scanpy DPT; same PT embedding and early endpoint'
saved.uns['reviewed_cluster_fingerprint'] = fingerprint
saved.uns['scfates_root_node'] = root_node
saved.uns['dpt_root_structure'] = str(pt.obs_names[root_cell])
saved.write_h5ad(output / 'cross_species_pt_scfates.h5ad')
assert saved.obs_names.equals(pt.obs_names)
assert saved.obs.shared_pseudospace.between(0, 1).all()
print('Notebook 12 inputs:', output / 'cross_species_pt_scfates.h5ad',
      output / 'ortholog_map_used.csv', sep='\n  ')
