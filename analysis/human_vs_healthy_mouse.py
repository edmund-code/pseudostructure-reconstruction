"""Human versus healthy-mouse nephron pseudospace reconstruction.

The reconstruction follows the validated mouse-only workflow: reviewed marker panels, deterministic
R Harmony 2.0.5, one pass-1 clustering checkpoint, a label-preserving pass-2 integration, and
Scanpy DPT rooted and oriented by the same nephron marker logic. Human counts are first converted
to an auditable mouse-symbol ortholog space.

The comparison contains one human donor sampled in cortex and medulla and two mouse control
specimens. Species differences are therefore descriptive effect sizes, not confirmatory tests.
"""

# %% [markdown]
# # Human versus healthy-mouse pseudospace
#
# This notebook rebuilds a shared nephron trajectory from the two healthy mouse controls and the
# available human cortex/medulla samples. It deliberately mirrors the finalized mouse-only
# reconstruction wherever the data permit. Cross-species mapping and the one-human-donor
# inference boundary are the two explicit differences.

# %% [markdown]
# # Section 0 - configuration
#
# Paths, deterministic integration parameters, support-based HCOP mapping, cell-typing panels,
# label vocabulary, and heatmap panels. All downstream code is defined in this notebook.

# %%
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import warnings
from pathlib import Path

# These must be set before importing matplotlib, numba, Scanpy, or modules that import them.
CACHE_ROOT = Path(os.environ.get('PSEUDOSPACE_CACHE_ROOT', '/tmp/pseudospace_human_mouse_cache'))
for _subdir in ('matplotlib', 'numba', 'xdg'):
    (CACHE_ROOT / _subdir).mkdir(parents=True, exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR', str(CACHE_ROOT / 'matplotlib'))
os.environ.setdefault('NUMBA_CACHE_DIR', str(CACHE_ROOT / 'numba'))
os.environ.setdefault('XDG_CACHE_HOME', str(CACHE_ROOT / 'xdg'))

import numpy as np
import pandas as pd

warnings.filterwarnings('ignore', category=FutureWarning)


def _find_project_dir() -> Path:
    starts = [Path.cwd().resolve()]
    if '__file__' in globals():
        starts.insert(0, Path(__file__).resolve().parent)
    for start in starts:
        for candidate in (start, *start.parents):
            if (candidate / 'pseudospace').is_dir() and (candidate / 'data').is_dir():
                return candidate
    raise RuntimeError('Could not locate the repository root containing pseudospace/ and data/')


def _workflow_roots(project_dir: Path) -> tuple[Path, Path]:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--data-root', type=Path)
    parser.add_argument('--results-root', type=Path)
    args, _ = parser.parse_known_args()
    data_root = args.data_root or os.environ.get('PSEUDOSPACE_DATA_ROOT') or project_dir / 'data'
    results_root = (args.results_root or os.environ.get('PSEUDOSPACE_RESULTS_ROOT')
                    or project_dir / 'results')
    return Path(data_root).expanduser().resolve(), Path(results_root).expanduser().resolve()


PROJECT_DIR = _find_project_dir()
if str(PROJECT_DIR) not in sys.path:
    sys.path.insert(0, str(PROJECT_DIR))
DATA_ROOT, RESULTS_ROOT = _workflow_roots(PROJECT_DIR)

TUBULE_BY_GENE_DIR = DATA_ROOT / 'tubule_by_gene'
ORTHOLOG_TABLE_PATH = DATA_ROOT / 'human_mouse_hcop_fifteen_column.txt.gz'
PATHWAY_LIBRARY_DIR = DATA_ROOT / 'mouse_vs_human' / 'pathway_gene_sets'

RESULTS_DIR = RESULTS_ROOT / 'human_vs_healthy_mouse'
HARMONY_OUTPUT_PATH = RESULTS_DIR / 'cross_species_harmony_pass1.h5ad'
PASS2_HARMONY_OUTPUT_PATH = RESULTS_DIR / 'cross_species_harmony_pass2.h5ad'
DPT_OUTPUT_PATH = RESULTS_DIR / 'cross_species_dpt.h5ad'
CELLTYPING_DIR = RESULTS_DIR / 'celltyping'
HEATMAP_OUTPUT_DIR = RESULTS_DIR / 'heatmaps'
CURVE_OUTPUT_DIR = RESULTS_DIR / 'curves'
DIAGNOSTIC_DIR = RESULTS_DIR / 'diagnostics'
for _directory in (RESULTS_DIR, CELLTYPING_DIR, HEATMAP_OUTPUT_DIR, CURVE_OUTPUT_DIR,
                   DIAGNOSTIC_DIR):
    _directory.mkdir(parents=True, exist_ok=True)

MOUSE_SAMPLES = ('Ctrl1A2', 'Ctrl1A4')
HUMAN_SAMPLES = ('HUK1_COR1', 'HUK1_MED1')
SAMPLE_ORDER = (*MOUSE_SAMPLES, *HUMAN_SAMPLES)
# Note: HUK1_MED1 is named 'medulla' in raw Visium HD metadata for catalog traceability,
# but nephron segment composition (1,897 PT and 292 DCT1 tubules) demonstrates it is Cortex 2.
HUMAN_REGIONS = {'HUK1_COR1': 'cortex', 'HUK1_MED1': 'medulla'}
SPECIES_GROUPS = ('mouse', 'human')
BATCH_KEY = 'sample'  # Both Pass 1 and Pass 2 integrate on sample (Ctrl1A2, Ctrl1A4, HUK1_COR1, HUK1_MED1)
RANDOM_STATE = 0

MIN_GENES_PER_TUBULE = 100
MIN_GENE_TUBULE_FRACTION = 0.05
MIN_GENE_TOTAL_COUNTS = 20
NORMALIZE_TARGET_SUM = 1e4
N_HVGS = 2000

HARMONY_EXPECTED_VERSION = '2.0.5'
HARMONY_HVG_MIN_MEAN = 0.0125
HARMONY_HVG_MAX_MEAN = 3
HARMONY_HVG_MIN_DISP = 0.5
HARMONY_N_PCS = 50
HARMONY_THETA = 6
HARMONY_LAMBDA = 1
HARMONY_MAX_ITER = 30
HARMONY_TAU = 0
HARMONY_PCA_N_COMPS = 50
HARMONY_NEIGHBORS_N = 25
HARMONY_UMAP_MIN_DIST = 0.3
HARMONY_UMAP_SPREAD = 1.0

N_NEIGHBORS = 30
COARSE_RESOLUTION = 0.7
INTEGRATION_QC_MAX_CELLS = 10000
PAGA_CONNECTIVITY_THRESHOLD = 0.01
N_DIFFMAP_COMPONENTS_TO_TEST = 10
EIGENVALUE_FLOOR = 1e-8
MIN_VALID_FOR_SPEARMAN = 20
N_BINS = 120
GAM_LAMBDA_GRID = np.logspace(-3, 3, 13)

# Mouse symbols are used throughout after human-to-mouse conversion.
NEPHRON_AXIS_MARKERS = {
    'Podocyte':     ['Nphs1', 'Nphs2', 'Podxl'],
    'PT-S1':        ['Slc5a2', 'Slc5a12', 'Gatm'],
    'PT-S2':        ['Slc22a6', 'Slc13a3', 'Cyp2e1'],
    'PT-S3':        ['Slc7a13', 'Slc22a7', 'Cyp7b1'],
    'DTL1':         ['Corin', 'Uncx', 'Slc14a2'],
    'DTL2':         ['Fst', 'Cdh13', 'Stk32a'],
    'DTL3':         ['Nr2e3', 'Dmkn'],
    'ATL':          ['Clcnka', 'Sptssb', 'Akr1b3'],
    'mTAL':         ['Slc12a1', 'Umod', 'Cldn10', 'Ptger3'],
    'cTAL':         ['Slc12a1', 'Umod', 'Cldn16', 'Kcnj10'],
    'Macula-densa': ['Nos1', 'Ptgs2', 'Pappa2'],
    'DCT1':         ['Slc12a3', 'Pvalb', 'Trpm6', 'Egf'],
    'DCT2':         ['Slc12a3', 'Trpv5', 'Slc8a1', 'S100g'],
    'CNT':          ['Calb1', 'Slc8a1', 'Aqp2', 'Rhcg'],
    'CCD':          ['Aqp2', 'Aqp3', 'Kit'],
    'OMCD':         ['Aqp2', 'Atp6v0d2', 'Rhcg'],
    'IMCD':         ['Aqp2', 'Aqp4', 'Slc14a2', 'Wnt7b'],
}
NEPHRON_SEGMENT_ORDER = list(NEPHRON_AXIS_MARKERS)

NON_TUBULE_MARKERS = {
    'Vessel': ['Pecam1', 'Cdh5', 'Kdr', 'Emcn', 'Flt1', 'Egfl7'],
    'Stroma': ['Col1a1', 'Col3a1', 'Col1a2', 'Dcn', 'Pdgfrb', 'Mgp'],
    'SmoothMuscle': ['Acta2', 'Myh11', 'Tagln', 'Rgs5', 'Pdgfra'],
    'Immune': ['Ptprc', 'C1qa', 'C1qb', 'Cd52', 'Lyz2', 'Cd74'],
}

SEGMENT_TO_COARSE = {
    'Podocyte': 'Glomerulus',
    'PT-S1': 'PT', 'PT-S2': 'PT', 'PT-S3': 'PT',
    'DTL1': 'DTL', 'DTL2': 'DTL', 'DTL3': 'DTL',
    'ATL': 'AL', 'mTAL': 'AL', 'cTAL': 'AL', 'Macula-densa': 'AL',
    'DCT1': 'DCT', 'DCT2': 'DCT',
    'CNT': 'CNT_CD', 'CCD': 'CNT_CD', 'OMCD': 'CNT_CD', 'IMCD': 'CNT_CD',
}
COARSE_ORDER = ['PT', 'DTL', 'AL', 'DCT', 'CNT_CD']
KEEP_TUBULE_CLASSES = list(COARSE_ORDER)
REMOVE_CLASSES = ['Glomerulus', 'Vessel', 'Stroma', 'SmoothMuscle', 'Immune', 'Unassigned']
SEGMENT_TO_COARSE.update({family: family for family in COARSE_ORDER})
SEGMENT_TO_COARSE['Glomerulus'] = 'Glomerulus'
LABEL_VOCABULARY = (NEPHRON_SEGMENT_ORDER + COARSE_ORDER + ['Glomerulus']
                    + [value for value in REMOVE_CLASSES if value != 'Glomerulus'])

SEGMENT_DISPLAY_ORDER = []
for _segment in NEPHRON_SEGMENT_ORDER:
    _family = SEGMENT_TO_COARSE.get(_segment)
    if _family in COARSE_ORDER and _family not in SEGMENT_DISPLAY_ORDER:
        SEGMENT_DISPLAY_ORDER.append(_family)
    SEGMENT_DISPLAY_ORDER.append(_segment)
SEGMENT_DISPLAY_ORDER += [value for value in REMOVE_CLASSES if value not in SEGMENT_DISPLAY_ORDER]

TOTAL_POSITION_MARKERS = {
    'early': sorted({gene for group in ('PT-S1', 'PT-S2')
                     for gene in NEPHRON_AXIS_MARKERS[group]}),
    'late': sorted({gene for group in ('OMCD', 'IMCD')
                    for gene in NEPHRON_AXIS_MARKERS[group]}),
}

HEATMAP_MARKERS = {
    'PT': {
        'PT-S1': ['Slc5a2', 'Slc5a12', 'Gatm'],
        'PT-S2': ['Slc22a6', 'Slc13a3', 'Cyp2e1'],
        'PT-S3': ['Slc22a7', 'Cyp7b1'],
    },
    'AL': {
        'ATL': ['Sptssb', 'Tacstd2'],
        'mTAL': ['Slc12a1', 'Umod', 'Cldn10', 'Ptger3'],
        'cTAL': ['Kcnj10', 'Bcl6'],
    },
    'DCT': {
        'DCT1': ['Pvalb', 'Egf'],
        'DCT2': ['Slc12a3', 'Calb1', 'S100g'],
    },
    'CNT_CD': {
        'CNT': ['Calb1', 'Hsd11b2'],
        'CCD': ['Aqp2', 'Aqp3'],
        'OMCD': ['Atp6v0d2', 'Rhcg'],
        'IMCD': ['Aqp4', 'Slc14a2'],
    },
}
HEATMAP_COLORS = {
    'PT-S1': '#4C9BD3', 'PT-S2': '#5B8E55', 'PT-S3': '#D6B48A',
    'ATL': '#26A69A', 'mTAL': '#F58518', 'cTAL': '#E8A33D',
    'DCT1': '#D81B60', 'DCT2': '#EC6EA5',
    'CNT': '#76B7B2', 'CCD': '#4E79A7', 'OMCD': '#9C755F', 'IMCD': '#593C8F',
}
SEGMENT_STRIP_COLORS = {
    **HEATMAP_COLORS,
    'PT': '#4C9BD3', 'DTL': '#7E57C2', 'AL': '#F58518',
    'DCT': '#D81B60', 'CNT_CD': '#8D6E63', 'Glomerulus': '#9E9E9E',
}

CURATED_MARKER_WHITELIST = sorted({
    gene for genes in [*NEPHRON_AXIS_MARKERS.values(), *NON_TUBULE_MARKERS.values()]
    for gene in genes
} | {
    gene for family in HEATMAP_MARKERS.values() for genes in family.values() for gene in genes
})


def _marker_dict(markers: dict[str, list[str]]) -> dict[str, list[dict[str, object]]]:
    return {
        group: [{'label': gene, 'candidates': [gene]} for gene in genes]
        for group, genes in markers.items()
    }


def _flatten_marker_panel(groups: list[str]) -> dict[str, list[dict[str, object]]]:
    return {
        group: _marker_dict(markers)[group]
        for family in HEATMAP_MARKERS.values()
        for group, markers in family.items()
        if group in groups
    }


print(f'Results directory: {RESULTS_DIR.relative_to(PROJECT_DIR)}')
print(f'Mouse controls: {MOUSE_SAMPLES}; human samples: {HUMAN_SAMPLES}')
print('Ortholog mapping: HCOP reciprocal best matches with support from at least 3 databases')
print(f'Cell-typing/heatmap genes protected from filtering: {len(CURATED_MARKER_WHITELIST)}')

# %% [markdown]
# # Section 1 - pass-1 Harmony integration
#
# Human genes are converted to HCOP support-based reciprocal best mouse matches before concatenation.
# HVGs are then selected independently in mouse and human and intersected. As in the mouse-only
# workflow, sample is the Harmony batch; species is retained as biology and is never the batch key.

# %%
import scanpy as sc
import matplotlib.pyplot as plt
import rpy2.robjects as ro
from IPython.display import Image, display
from scipy import sparse
from scipy.ndimage import gaussian_filter1d
from scipy.spatial import cKDTree
from scipy.stats import spearmanr
from sklearn.metrics import silhouette_score

from pseudospace.cross_species import (
    build_one_to_one_ortholog_map,
    combine_cross_species,
    load_cross_species_samples,
    plot_cross_species_marker_alignment,
    plot_pre_filtering_tubule_qc,
    read_ortholog_table,
)
from pseudospace.io_qc import annotate_mito_ribo_mouse_symbols, sample_name_from_path
from pseudospace.harmony import run_harmony_rpy2, select_harmony_hvgs_by_condition
from pseudospace.markers import (
    assign_cluster_labels,
    auto_discover_cluster_identities,
    available_markers_by_group,
    build_gene_lookup,
    plot_coarse_cluster_visualizations,
)
from pseudospace.heatmaps import plot_marker_heatmap
from pseudospace.trajectory import (
    choose_diffusion_component,
    choose_root_global,
    choose_root_pt_cluster,
    orient_and_normalize,
    recompute_subset_dpt,
    save_total_pseudotime_anndata,
)

import matplotlib as _mpl
import matplotlib.patheffects as PathEffects
_mpl_cfg = Path(_mpl.get_configdir()).resolve()
if not _mpl_cfg.is_relative_to(CACHE_ROOT.resolve()):
    warnings.warn(
        f'matplotlib is caching to {_mpl_cfg}, not under {CACHE_ROOT.resolve()}. '
        'Run the configuration cell before importing Scanpy/matplotlib.', stacklevel=2)
else:
    print(f'Cache redirects active ({CACHE_ROOT}).')

_installed_harmony = str(ro.r('as.character(packageVersion("harmony"))')[0])
if _installed_harmony != HARMONY_EXPECTED_VERSION:
    raise RuntimeError(
        f'R harmony {HARMONY_EXPECTED_VERSION} is required for the pinned clustering; '
        f'found {_installed_harmony}. Activate the kidney-pseudospace environment.'
    )
print(f'R harmony version: {_installed_harmony}')

# %%
if not ORTHOLOG_TABLE_PATH.exists():
    raise FileNotFoundError(f'Missing required HCOP table: {ORTHOLOG_TABLE_PATH}')

ortholog_table = read_ortholog_table(ORTHOLOG_TABLE_PATH)
ambiguous_human_symbols = ortholog_table.loc[
    ortholog_table.groupby(ortholog_table['human_symbol'].str.upper())['mouse_symbol']
    .transform('nunique').gt(1),
    'human_symbol',
].astype(str).unique().tolist()
ortholog_map = build_one_to_one_ortholog_map(
    ortholog_table,
    min_support=3,
)
ortholog_map.to_csv(RESULTS_DIR / 'ortholog_map_used.csv', index=False)

files = sorted(TUBULE_BY_GENE_DIR.glob('*_tubule_by_gene_caleb.h5ad'))
found_samples = sorted(sample_name_from_path(path) for path in files
                       if sample_name_from_path(path) in SAMPLE_ORDER)
if found_samples != sorted(SAMPLE_ORDER):
    raise FileNotFoundError(
        f'Expected exactly {sorted(SAMPLE_ORDER)}, found {found_samples}. '
        'Check data/tubule_by_gene/.'
    )

adatas, human_mapping_report = load_cross_species_samples(
    TUBULE_BY_GENE_DIR,
    files,
    ortholog_map,
    mouse_samples=MOUSE_SAMPLES,
    human_samples=HUMAN_SAMPLES,
    human_regions=HUMAN_REGIONS,
    ambiguous_human_symbols=ambiguous_human_symbols,
)
human_mapping_report.to_csv(RESULTS_DIR / 'human_gene_mapping_report.csv', index=False)

adata_combined = combine_cross_species(adatas)
adata_combined = annotate_mito_ribo_mouse_symbols(adata_combined)
print(f'Combined shared space: {adata_combined.n_obs:,} structures x '
      f'{adata_combined.n_vars:,} ortholog genes')

# %%
print(f'Before cell filtering: {adata_combined.n_obs:,}')
sc.pp.filter_cells(adata_combined, min_genes=MIN_GENES_PER_TUBULE)
print(f'After cell filtering:  {adata_combined.n_obs:,}')

adata_combined.layers['counts'] = adata_combined.X.copy()
sc.pp.normalize_total(adata_combined, target_sum=NORMALIZE_TARGET_SUM)
sc.pp.log1p(adata_combined)
adata_combined.layers['lognorm'] = adata_combined.X.copy()

adata_combined = select_harmony_hvgs_by_condition(
    adata_combined,
    group_key='comparison_species',
    groups=SPECIES_GROUPS,
    mode='intersection',
    min_mean=HARMONY_HVG_MIN_MEAN,
    max_mean=HARMONY_HVG_MAX_MEAN,
    min_disp=HARMONY_HVG_MIN_DISP,
)
hvg_col = 'highly_variable_for_harmony'
if int(adata_combined.var[hvg_col].sum()) < 50:
    raise ValueError('Fewer than 50 mouse/human-conserved HVGs were retained.')

cohort_summary = (
    adata_combined.obs.groupby(['comparison_species', 'sample', 'region'], observed=True)
    .size().rename('n_structures').reset_index()
)
cohort_summary.to_csv(RESULTS_DIR / 'cohort_summary.csv', index=False)
display(cohort_summary)

# %%
# Uncorrected PCA/UMAP is retained as the before-Harmony reference.
adata_hvg = adata_combined[:, adata_combined.var[hvg_col]].copy()
adata_hvg.X = adata_hvg.layers['lognorm'].copy()
sc.tl.pca(adata_hvg, n_comps=HARMONY_PCA_N_COMPS, random_state=RANDOM_STATE)
adata_combined.obsm['X_pca'] = adata_hvg.obsm['X_pca'].copy()

sc.pp.neighbors(
    adata_combined, use_rep='X_pca', n_neighbors=HARMONY_NEIGHBORS_N,
    key_added='pre_harmony', random_state=RANDOM_STATE,
)
sc.tl.umap(
    adata_combined, neighbors_key='pre_harmony', min_dist=HARMONY_UMAP_MIN_DIST,
    spread=HARMONY_UMAP_SPREAD, random_state=RANDOM_STATE,
)
adata_combined.obsm['X_umap_pre_harmony'] = adata_combined.obsm['X_umap'].copy()
sc.pl.umap(
    adata_combined, color=['sample', 'comparison_species'], ncols=2, frameon=False,
    title=['Pre-Harmony: sample', 'Pre-Harmony: species'],
)

# %%
print(f'Running deterministic R Harmony on {adata_hvg.n_obs:,} structures x '
      f'{adata_hvg.n_vars:,} shared HVGs.')
ro.r(f'set.seed({RANDOM_STATE})')
adata_hvg = run_harmony_rpy2(
    adata_hvg,
    batch_key=BATCH_KEY,
    n_pcs=HARMONY_N_PCS,
    theta=HARMONY_THETA,
    lambda_val=HARMONY_LAMBDA,
    max_iter=HARMONY_MAX_ITER,
    tau=HARMONY_TAU,
)
adata_combined.obsm['X_harmony'] = adata_hvg.obsm['X_harmony'].copy()
adata_combined.obsm['X_pca'] = adata_hvg.obsm['X_pca'].copy()

sc.pp.neighbors(
    adata_combined, use_rep='X_harmony', n_neighbors=HARMONY_NEIGHBORS_N,
    random_state=RANDOM_STATE,
)
sc.tl.umap(
    adata_combined, min_dist=HARMONY_UMAP_MIN_DIST, spread=HARMONY_UMAP_SPREAD,
    random_state=RANDOM_STATE,
)
sc.pl.umap(
    adata_combined,
    color=['sample', 'comparison_species', 'region'],
    ncols=3,
    frameon=False,
)
adata_combined.write(HARMONY_OUTPUT_PATH)
print(f'Saved pass-1 object: {HARMONY_OUTPUT_PATH.relative_to(PROJECT_DIR)}')

# %% [markdown]
# # Section 2 - cell typing, pass-2 Harmony, and shared Scanpy DPT
#
# The sequence is identical to the finalized mouse notebook: whitelist-protected gene filtering,
# one coarse Leiden run, marker dotplot and reviewed manual labels, removal of non-nephron
# structures, label-preserving pass-2 Harmony, integration checks, and a marker-oriented DPT.
# Labels are never recomputed after pass 2.

# %%
# Whitelist-protected expression filter.
adata_all = sc.read_h5ad(HARMONY_OUTPUT_PATH)
counts = adata_all.layers['counts'].copy() if 'counts' in adata_all.layers else adata_all.X.copy()
adata_all.layers['counts'] = counts.copy()
expressed = np.asarray((counts > 0).sum(axis=0)).ravel()
total_counts = np.asarray(counts.sum(axis=0)).ravel()
min_cells = max(1, int(np.ceil(MIN_GENE_TUBULE_FRACTION * adata_all.n_obs)))
gene_keep = (expressed >= min_cells) & (total_counts >= MIN_GENE_TOTAL_COUNTS)
whitelist_hits = adata_all.var_names.str.upper().isin(
    {gene.upper() for gene in CURATED_MARKER_WHITELIST}
)
gene_keep |= np.asarray(whitelist_hits)

adata_all.var['n_structures_expressed'] = expressed
adata_all.var['total_counts'] = total_counts
adata_all.var['passes_expression_count_filter'] = gene_keep
adata_all = adata_all[:, gene_keep].copy()

adata_all.X = adata_all.layers['counts'].copy()
adata_all.uns.pop('log1p', None)
sc.pp.normalize_total(adata_all, target_sum=NORMALIZE_TARGET_SUM)
sc.pp.log1p(adata_all)
adata_all.layers['lognorm'] = adata_all.X.copy()
adata_all_expression_filtered = adata_all.copy()
print(f'Expression-filtered matrix: {adata_all.n_obs:,} x {adata_all.n_vars:,}')

# %%
# The one clustering run, on the pass-1 Harmony representation.
gene_lookup = build_gene_lookup(adata_all)
sc.pp.highly_variable_genes(
    adata_all, n_top_genes=min(N_HVGS, adata_all.n_vars), flavor='seurat'
)
whitelist_present = {
    gene_lookup[gene.upper()]
    for gene in CURATED_MARKER_WHITELIST
    if gene.upper() in gene_lookup
}
feature_mask = (
    adata_all.var['highly_variable'].to_numpy(dtype=bool)
    | adata_all.var_names.isin(whitelist_present)
)
adata_all.var['selected_for_clustering'] = feature_mask
adata_cluster = adata_all[:, feature_mask].copy()
gene_lookup = build_gene_lookup(adata_cluster)

sc.pp.neighbors(
    adata_cluster, n_neighbors=N_NEIGHBORS, use_rep='X_harmony',
    random_state=RANDOM_STATE,
)
sc.tl.leiden(
    adata_cluster,
    resolution=COARSE_RESOLUTION,
    key_added='leiden_coarse',
    flavor='igraph',
    n_iterations=2,
    directed=False,
    random_state=RANDOM_STATE,
)
sc.tl.umap(
    adata_cluster, min_dist=HARMONY_UMAP_MIN_DIST, spread=HARMONY_UMAP_SPREAD,
    random_state=RANDOM_STATE,
)
sc.tl.rank_genes_groups(
    adata_cluster, 'leiden_coarse', method='wilcoxon', pts=True,
    key_added='rank_coarse',
)

CLUSTERING_FINGERPRINT = {
    'n_cells': int(adata_cluster.n_obs),
    'n_features': int(adata_cluster.n_vars),
    'resolution': COARSE_RESOLUTION,
    'n_neighbors': N_NEIGHBORS,
    'random_state': RANDOM_STATE,
    'n_clusters': int(adata_cluster.obs['leiden_coarse'].nunique()),
    'membership_sha1': hashlib.sha1(
        ','.join(adata_cluster.obs['leiden_coarse'].astype(str)).encode()
    ).hexdigest()[:12],
}
print('Clustering fingerprint:')
for _key, _value in CLUSTERING_FINGERPRINT.items():
    print(f'  {_key}: {_value}')

# %%
coarse_panel = {**NEPHRON_AXIS_MARKERS, **NON_TUBULE_MARKERS}
resolved_panel, missing_panel = available_markers_by_group(
    adata_cluster, coarse_panel, gene_lookup
)
print('Markers absent from the shared ortholog space:', missing_panel)

sc.pl.dotplot(
    adata_cluster, resolved_panel, groupby='leiden_coarse',
    standard_scale='var', show=False,
)
plt.savefig(CELLTYPING_DIR / 'coarse_marker_dotplot.png', dpi=180, bbox_inches='tight')
plt.show()

fig, axes = plt.subplots(1, 3, figsize=(18, 5))
sc.pl.embedding(adata_cluster, basis='umap', color='leiden_coarse', ax=axes[0], frameon=False, show=False)
coords_coarse = adata_cluster.obsm['X_umap']
for c_id in sorted(adata_cluster.obs['leiden_coarse'].unique(), key=lambda x: int(x) if str(x).isdigit() else str(x)):
    mask = (adata_cluster.obs['leiden_coarse'] == c_id).to_numpy()
    if mask.any():
        x_c, y_c = np.median(coords_coarse[mask], axis=0)
        txt = axes[0].text(x_c, y_c, str(c_id), ha='center', va='center', fontsize=11, fontweight='bold', color='black')
        txt.set_path_effects([PathEffects.withStroke(linewidth=3.5, foreground='white')])

sc.pl.embedding(adata_cluster, basis='umap', color='comparison_species', ax=axes[1], frameon=False, show=False)
sc.pl.embedding(adata_cluster, basis='umap', color='sample', ax=axes[2], frameon=False, show=False)
plt.tight_layout()
plt.savefig(CELLTYPING_DIR / 'coarse_leiden_umap.png', dpi=180, bbox_inches='tight')
plt.show()

top_de = sc.get.rank_genes_groups_df(adata_cluster, group=None, key='rank_coarse')
top20 = (
    top_de.sort_values(['group', 'scores'], ascending=[True, False])
    .groupby('group', observed=True).head(20)
)
top20.to_csv(CELLTYPING_DIR / 'coarse_cluster_top_markers_full.csv', index=False)
top8 = top20.groupby('group', observed=True)['names'].apply(
    lambda values: '; '.join(values.head(8))
)
top8.to_csv(CELLTYPING_DIR / 'coarse_cluster_top_markers.csv')
print(top8.to_string())

# %%
# Panel-score suggestions and per-cell second opinion. The manual map below remains authoritative.
score_columns = {}
for panel_name, genes in resolved_panel.items():
    score_name = f'_score_{panel_name}'
    sc.tl.score_genes(
        adata_cluster, genes, score_name=score_name, use_raw=False,
        random_state=RANDOM_STATE,
    )
    score_columns[panel_name] = score_name

panel_scores = (
    adata_cluster.obs.groupby('leiden_coarse', observed=True)[list(score_columns.values())]
    .mean().rename(columns={value: key for key, value in score_columns.items()})
)
panel_scores.to_csv(CELLTYPING_DIR / 'coarse_cluster_panel_scores.csv')

score_matrix = adata_cluster.obs[list(score_columns.values())].to_numpy(dtype=float)
panel_names = np.array(list(score_columns))
ranked = np.argsort(score_matrix, axis=1)
primary_segment = panel_names[ranked[:, -1]]
secondary_segment = panel_names[ranked[:, -2]]
adata_cluster.obs['celltype_primary_segment'] = primary_segment
adata_cluster.obs['celltype_secondary_segment'] = secondary_segment
adata_cluster.obs['celltype_primary'] = [
    SEGMENT_TO_COARSE.get(value, value) for value in primary_segment
]
adata_cluster.obs['celltype_secondary'] = [
    SEGMENT_TO_COARSE.get(value, value) for value in secondary_segment
]
adata_cluster.obs['celltype_margin'] = (
    np.take_along_axis(score_matrix, ranked[:, -1:], 1).ravel()
    - np.take_along_axis(score_matrix, ranked[:, -2:-1], 1).ravel()
)
PER_CELL_CALL_COLUMNS = [
    'celltype_primary', 'celltype_secondary',
    'celltype_primary_segment', 'celltype_secondary_segment', 'celltype_margin',
]
adata_cluster.obs.drop(columns=list(score_columns.values()), inplace=True)

display(panel_scores.round(3))
display(pd.crosstab(
    adata_cluster.obs['leiden_coarse'],
    adata_cluster.obs['celltype_primary_segment'],
))
print('Suggested labels require visual review:')
for cluster, row in panel_scores.iterrows():
    print(f'  {str(cluster)!r}: {str(row.idxmax())!r}')

# %% [markdown]
# ## Dynamic Pass-1 Cluster Discovery and Biological Annotations
#
# Rather than enforcing a fragile, hardcoded manual checkpoint with hash checks that fails closed
# on any upstream parameter variation, cluster segment identities and species breakdowns are
# dynamically rediscovered after each run using marker panel scoring, Wilcoxon differential
# expression, and biological nephron segmentation heuristics. Diagnostic visualizations
# (species breakdown, panel scores heatmap, diagnostic DE markers) are automatically generated.

# %%
auto_labels, cluster_summary = auto_discover_cluster_identities(
    adata_cluster,
    panel_scores,
    top_de,
    segment_to_coarse=SEGMENT_TO_COARSE,
    remove_classes=REMOVE_CLASSES,
    coarse_order=COARSE_ORDER,
)

cluster_summary.to_csv(CELLTYPING_DIR / 'coarse_cluster_summary.csv', index=False)
print('\n=== Dynamically Discovered Pass-1 Cluster Identities & Breakdowns ===')
print(cluster_summary.to_string(index=False))

# Generate diagnostic visualizations explaining what each cluster is
generated_plots = plot_coarse_cluster_visualizations(
    adata_cluster,
    panel_scores,
    top_de,
    cluster_summary,
    CELLTYPING_DIR,
)
print('Generated cluster diagnostic visualizations:', generated_plots)
for _fname in (
    'coarse_cluster_tubule_likelihood_matrix.png',
    'coarse_cluster_tubule_composition_stacked.png',
    'unknown_cluster_diagnostic_profile.png',
):
    _plot_file = CELLTYPING_DIR / _fname
    if _plot_file.is_file():
        display(Image(filename=str(_plot_file)))

assign_cluster_labels(
    adata_cluster, 'leiden_coarse', auto_labels, 'segment_class',
    unknown='Unassigned',
)
for column in ['leiden_coarse', 'segment_class', *PER_CELL_CALL_COLUMNS]:
    values = adata_cluster.obs[column].reindex(adata_all.obs_names)
    if values.isna().any():
        raise RuntimeError(f'Failed to transfer {column} to the expression object.')
    adata_all.obs[column] = (
        values.to_numpy(dtype=float)
        if column == 'celltype_margin'
        else values.astype(str).to_numpy()
    )

adata_all.obs['coarse_class'] = (
    adata_all.obs['segment_class']
    .map(lambda label: SEGMENT_TO_COARSE.get(
        label, label if label in REMOVE_CLASSES else 'Unassigned'
    ))
    .astype(str)
)
adata_all.obs['broad_tubule_marker_call'] = adata_all.obs['coarse_class']
display(pd.crosstab(adata_all.obs['segment_class'], adata_all.obs['coarse_class']))

# Coarse Leiden annotations UMAP with text annotations at cluster centroids
fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
sc.pl.embedding(
    adata_cluster, basis='umap', color='segment_class',
    ax=axes[0], show=False, frameon=False, title='Coarse Leiden annotations',
)
coords_annot = adata_cluster.obsm['X_umap']
for c_id in sorted(adata_cluster.obs['leiden_coarse'].unique(), key=lambda x: int(x) if str(x).isdigit() else str(x)):
    mask = (adata_cluster.obs['leiden_coarse'] == c_id).to_numpy()
    if mask.any():
        x_c, y_c = np.median(coords_annot[mask], axis=0)
        seg = auto_labels.get(str(c_id), str(c_id))
        txt = axes[0].text(x_c, y_c, f'{c_id}: {seg}', ha='center', va='center', fontsize=10, fontweight='bold', color='black')
        txt.set_path_effects([PathEffects.withStroke(linewidth=3.5, foreground='white')])

sc.pl.embedding(
    adata_cluster, basis='umap', color='celltype_margin',
    ax=axes[1], show=False, frameon=False, color_map='viridis',
    title='Panel-score margin (second-opinion confidence)',
)
_annot_qc_path = CELLTYPING_DIR / 'coarse_leiden_annotations_umap.png'
plt.savefig(_annot_qc_path, dpi=180, bbox_inches='tight')
plt.show()

# %%
# Preserve every structure and label on pass 1, then remove non-nephron classes for pass 2.
adata_pass1 = sc.read_h5ad(HARMONY_OUTPUT_PATH)
for column in [
    'segment_class', 'coarse_class', 'broad_tubule_marker_call', 'leiden_coarse',
    *PER_CELL_CALL_COLUMNS,
]:
    values = adata_all.obs.loc[adata_pass1.obs_names, column]
    adata_pass1.obs[column] = (
        values.to_numpy(dtype=float)
        if column == 'celltype_margin'
        else values.astype(str).to_numpy()
    )
adata_pass1.write(HARMONY_OUTPUT_PATH)

# Pre-filtering tubule QC visualization to audit retained vs removed structures
plot_pre_filtering_tubule_qc(
    adata_pass1,
    output_path=CELLTYPING_DIR / 'pre_filtering_tubule_qc.png',
)

keep_mask = adata_all.obs['coarse_class'].isin(KEEP_TUBULE_CLASSES).to_numpy()
for class_name in REMOVE_CLASSES:
    n_removed = int(adata_all.obs['coarse_class'].eq(class_name).sum())
    if n_removed:
        print(f'Removing {n_removed:,} {class_name} structures')
adata_tubule = adata_all[keep_mask].copy()
present_segments = set(adata_tubule.obs['segment_class'].astype(str))
adata_tubule.obs['coarse_class'] = pd.Categorical(
    adata_tubule.obs['coarse_class'], categories=COARSE_ORDER, ordered=True
)
adata_tubule.obs['segment_class'] = pd.Categorical(
    adata_tubule.obs['segment_class'],
    categories=[value for value in SEGMENT_DISPLAY_ORDER if value in present_segments],
    ordered=True,
)
print(f'Tubule-only matrix: {adata_tubule.n_obs:,} x {adata_tubule.n_vars:,}')
display(pd.crosstab(
    adata_tubule.obs['segment_class'],
    adata_tubule.obs['comparison_species'],
))

# %%
# Pass-2 Harmony on retained nephron structures; labels are carried, never recomputed.
adata_tubule = select_harmony_hvgs_by_condition(
    adata_tubule,
    group_key='comparison_species',
    groups=SPECIES_GROUPS,
    mode='intersection',
    min_mean=HARMONY_HVG_MIN_MEAN,
    max_mean=HARMONY_HVG_MAX_MEAN,
    min_disp=HARMONY_HVG_MIN_DISP,
)
tub_hvg = 'highly_variable_for_harmony'
adata_tubule_hvg = adata_tubule[:, adata_tubule.var[tub_hvg]].copy()
adata_tubule_hvg.X = adata_tubule_hvg.layers['lognorm'].copy()
sc.tl.pca(
    adata_tubule_hvg, n_comps=HARMONY_PCA_N_COMPS,
    random_state=RANDOM_STATE,
)
adata_tubule.obsm['X_pca'] = adata_tubule_hvg.obsm['X_pca'].copy()

ro.r(f'set.seed({RANDOM_STATE})')
adata_tubule_hvg = run_harmony_rpy2(
    adata_tubule_hvg,
    batch_key=BATCH_KEY,
    n_pcs=HARMONY_N_PCS,
    theta=HARMONY_THETA,
    lambda_val=HARMONY_LAMBDA,
    max_iter=HARMONY_MAX_ITER,
    tau=HARMONY_TAU,
)
adata_tubule.obsm['X_harmony'] = adata_tubule_hvg.obsm['X_harmony'].copy()
sc.pp.neighbors(
    adata_tubule, use_rep='X_harmony', n_neighbors=HARMONY_NEIGHBORS_N,
    random_state=RANDOM_STATE,
)
sc.tl.umap(
    adata_tubule, min_dist=HARMONY_UMAP_MIN_DIST, spread=HARMONY_UMAP_SPREAD,
    random_state=RANDOM_STATE,
)
adata_tubule.write(PASS2_HARMONY_OUTPUT_PATH)

sc.pl.embedding(
    adata_tubule,
    basis='umap',
    color=['sample', 'comparison_species', 'coarse_class', 'segment_class'],
    ncols=2,
    frameon=False,
    show=False,
)
plt.savefig(CELLTYPING_DIR / 'pass2_umap.png', dpi=180, bbox_inches='tight')
plt.show()

# Side-by-side (Human vs Mouse) canonical marker alignment on the Pass-2 coordinates
plot_cross_species_marker_alignment(
    adata_tubule,
    output_path=CELLTYPING_DIR / 'cross_species_marker_alignment_pass2.png',
)
plt.show()

# %%
# Integration audit: batch mixing must be interpreted beside biological conservation.
_rng = np.random.default_rng(RANDOM_STATE)
_qc_idx = (
    _rng.choice(adata_tubule.n_obs, INTEGRATION_QC_MAX_CELLS, replace=False)
    if adata_tubule.n_obs > INTEGRATION_QC_MAX_CELLS
    else np.arange(adata_tubule.n_obs)
)
_batch_labels = adata_tubule.obs[BATCH_KEY].astype(str).to_numpy()[_qc_idx]
_segment_labels = adata_tubule.obs['coarse_class'].astype(str).to_numpy()[_qc_idx]
_species_labels = adata_tubule.obs['comparison_species'].astype(str).to_numpy()[_qc_idx]
integration_qc = pd.DataFrame([
    {
        'embedding': embedding,
        'sample_asw_lower_is_better': float(silhouette_score(
            adata_tubule.obsm[embedding][_qc_idx], _batch_labels
        )),
        'segment_asw_higher_is_better': float(silhouette_score(
            adata_tubule.obsm[embedding][_qc_idx], _segment_labels
        )),
        'species_asw_descriptive': float(silhouette_score(
            adata_tubule.obsm[embedding][_qc_idx], _species_labels
        )),
    }
    for embedding in ('X_pca', 'X_harmony')
])
integration_qc['theta'] = HARMONY_THETA
integration_qc['n_structures_scored'] = len(_qc_idx)
integration_qc.to_csv(CELLTYPING_DIR / 'integration_qc_silhouette.csv', index=False)
display(integration_qc.round(4))

# %%
# Shared early-to-late marker axis.
CONTINUUM_PREFIX = 'total'
tubule_lookup = build_gene_lookup(adata_tubule)
adata_total = adata_tubule.copy()
for role in ('early', 'late'):
    genes = [
        tubule_lookup[gene.upper()]
        for gene in TOTAL_POSITION_MARKERS[role]
        if gene.upper() in tubule_lookup
    ]
    if len(genes) < 2:
        raise ValueError(f'Insufficient {role} axis markers: {genes}')
    sc.tl.score_genes(
        adata_total, genes,
        score_name=f'{CONTINUUM_PREFIX}_{role}_trajectory_score',
        use_raw=False,
        random_state=RANDOM_STATE,
    )
    print(f'{role} axis genes: {genes}')

adata_total.obs['total_marker_axis'] = (
    adata_total.obs['total_late_trajectory_score'].astype(float)
    - adata_total.obs['total_early_trajectory_score'].astype(float)
)
adata_total.obs['total_segment_marker_call'] = (
    adata_total.obs['segment_class'].astype(str).to_numpy()
)

# %%
# Graph topology must be inspected before interpreting a single continuous DPT.
marker_axis = adata_total.obs['total_marker_axis'].to_numpy(dtype=float)
trajectory_neighbors = 'trajectory_neighbors'
n_neighbors = max(N_NEIGHBORS, int(np.sqrt(adata_total.n_obs)))
sc.pp.neighbors(
    adata_total,
    n_neighbors=n_neighbors,
    use_rep='X_harmony',
    key_added=trajectory_neighbors,
    random_state=RANDOM_STATE,
)
sc.tl.paga(
    adata_total, groups='segment_class', neighbors_key=trajectory_neighbors
)
sc.pl.paga(
    adata_total, threshold=PAGA_CONNECTIVITY_THRESHOLD,
    color='segment_class', frameon=False, show=False,
)
plt.savefig(CELLTYPING_DIR / 'paga_segment_topology.png', dpi=180, bbox_inches='tight')
plt.show()

paga_connectivity = pd.DataFrame(
    np.asarray(adata_total.uns['paga']['connectivities'].todense()),
    index=adata_total.obs['segment_class'].cat.categories,
    columns=adata_total.obs['segment_class'].cat.categories,
)
paga_connectivity.to_csv(CELLTYPING_DIR / 'paga_connectivity_matrix.csv')
present_order = [
    segment for segment in SEGMENT_DISPLAY_ORDER
    if segment in paga_connectivity.index
]
consecutive_paga = pd.DataFrame([
    {
        'from': left,
        'to': right,
        'paga_connectivity': float(paga_connectivity.loc[left, right]),
        'connected': bool(
            paga_connectivity.loc[left, right] >= PAGA_CONNECTIVITY_THRESHOLD
        ),
    }
    for left, right in zip(present_order[:-1], present_order[1:])
])
consecutive_paga.to_csv(
    CELLTYPING_DIR / 'paga_consecutive_segment_connectivity.csv', index=False
)
display(consecutive_paga.round(4))

# %%
sc.tl.diffmap(
    adata_total, neighbors_key=trajectory_neighbors, random_state=RANDOM_STATE
)
diagnostic_component, component_df = choose_diffusion_component(
    adata_total,
    marker_axis,
    n_components_to_test=N_DIFFMAP_COMPONENTS_TO_TEST,
    eigenvalue_floor=EIGENVALUE_FLOOR,
    min_valid=MIN_VALID_FOR_SPEARMAN,
)
print(f'Best marker-correlated diffusion component: {diagnostic_component} (diagnostic only)')
display(component_df)

iroot, root_cluster = choose_root_pt_cluster(
    adata_total,
    marker_axis,
    leiden_col='leiden_coarse',
    rep_key='X_harmony',
    family_col='broad_tubule_marker_call',
    pt_label='PT',
)
if iroot is None:
    iroot = choose_root_global(
        marker_axis, adata_total, 'X_harmony', bottom_quantile=0.01
    )
    root_cluster = 'global bottom-1% fallback'
adata_total.uns['iroot'] = int(iroot)
sc.tl.dpt(adata_total, neighbors_key=trajectory_neighbors)
total_scanpy_dpt = orient_and_normalize(
    adata_total.obs['dpt_pseudotime'].to_numpy(dtype=float),
    marker_axis,
    min_valid=MIN_VALID_FOR_SPEARMAN,
)
adata_total.obs['total_scanpy_dpt'] = total_scanpy_dpt
adata_total.obs['shared_pseudospace'] = total_scanpy_dpt

nonfinite = ~np.isfinite(total_scanpy_dpt)
if nonfinite.any():
    print(f'Dropping {int(nonfinite.sum()):,} root-disconnected structures.')
    adata_total = adata_total[~nonfinite].copy()
    marker_axis = marker_axis[~nonfinite]
    total_scanpy_dpt = total_scanpy_dpt[~nonfinite]

print('Spearman(DPT, marker axis) =',
      round(spearmanr(total_scanpy_dpt, marker_axis).correlation, 3))
sc.pl.embedding(
    adata_total,
    basis='umap',
    color=['total_scanpy_dpt', 'segment_class', 'comparison_species'],
    cmap='viridis',
    frameon=False,
    ncols=3,
    show=False,
)
plt.savefig(CELLTYPING_DIR / 'total_dpt_umap.png', dpi=180, bbox_inches='tight')
plt.show()

# %%
segment_dpt = (
    adata_total.obs.groupby('segment_class', observed=True)['total_scanpy_dpt']
    .agg(['size', 'mean', 'median'])
    .reindex([
        value for value in SEGMENT_DISPLAY_ORDER
        if value in set(adata_total.obs['segment_class'].astype(str))
    ])
)
segment_dpt.to_csv(CELLTYPING_DIR / 'dpt_by_segment.csv')
display(segment_dpt.round(3))

species_segment_dpt = (
    adata_total.obs.groupby(
        ['comparison_species', 'segment_class'], observed=True
    )['total_scanpy_dpt']
    .agg(['size', 'median'])
    .reset_index()
)
species_segment_dpt.to_csv(
    CELLTYPING_DIR / 'dpt_by_species_and_segment.csv', index=False
)
display(species_segment_dpt.round(3))

fig, axes = plt.subplots(1, 2, figsize=(16, 4.5), sharey=True)
order = [s for s in SEGMENT_DISPLAY_ORDER if s in set(adata_total.obs['segment_class'])]

for ax, species in zip(axes, SPECIES_GROUPS):
    sp_mask = adata_total.obs['comparison_species'].astype(str).eq(species)
    data = [
        adata_total.obs.loc[sp_mask & (adata_total.obs['segment_class'].astype(str) == seg),
                            'total_scanpy_dpt'].dropna().to_numpy(dtype=float)
        for seg in order
    ]
    ax.boxplot([d if len(d) > 0 else np.array([np.nan]) for d in data], showfliers=False)
    ax.set_xticks(range(1, len(order) + 1))
    ax.set_xticklabels([f'{seg}\n(n={len(d):,})' for seg, d in zip(order, data)])
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
    ax.set_title(f'{species.title()} (n={int(sp_mask.sum()):,})')
    ax.set_ylabel('total_scanpy_dpt' if ax == axes[0] else '')

fig.suptitle('Pseudospace by segment (anatomical order left to right)', y=1.02)
fig.savefig(CELLTYPING_DIR / 'dpt_by_segment.png', dpi=150, bbox_inches='tight')
plt.show()

# %%
# Save the canonical shared DPT file with filtered expression and pass-2 coordinates.
DPT_EXTRA_OBS_COLS = [
    'comparison_species', 'region', 'x_centroid', 'y_centroid',
    'segment_class', 'coarse_class',
    'leiden_coarse', 'celltype_primary_segment', 'celltype_secondary_segment',
    'total_marker_axis', 'total_early_trajectory_score',
    'total_late_trajectory_score', 'shared_pseudospace',
]
adata_dpt_saved = save_total_pseudotime_anndata(
    adata_total,
    'total_scanpy_dpt',
    DPT_OUTPUT_PATH,
    expression_filtered_adata=adata_all_expression_filtered,
    annotation_adata=adata_total,
    project_dir=PROJECT_DIR,
    forbidden_labels=REMOVE_CLASSES + ['nan', 'none'],
    assert_no_nan=True,
    extra_obs_cols=DPT_EXTRA_OBS_COLS,
)
positions = adata_total.obs_names.get_indexer(adata_dpt_saved.obs_names)
if (positions < 0).any():
    raise RuntimeError('DPT output contains structures absent from adata_total.')

for column in DPT_EXTRA_OBS_COLS:
    if column in adata_total.obs.columns:
        values = adata_total.obs[column].to_numpy()[positions]
        adata_dpt_saved.obs[column] = (
            values.astype(float)
            if pd.api.types.is_numeric_dtype(adata_total.obs[column])
            else values.astype(str)
        )
for key in ('X_harmony', 'X_pca', 'X_umap'):
    if key in adata_total.obsm:
        adata_dpt_saved.obsm[key] = np.asarray(adata_total.obsm[key])[positions]
adata_dpt_saved.write(DPT_OUTPUT_PATH)

print(f'Saved canonical shared DPT: {DPT_OUTPUT_PATH.relative_to(PROJECT_DIR)}')
print(f'Root cluster: {root_cluster}; retained structures: {adata_dpt_saved.n_obs:,}')

# %% [markdown]
# # Section 3 - marker heatmaps and family-specific DPT
#
# Marker panels and display logic match the mouse-only workflow, but each plot is shown separately
# for mouse and human. Family-specific DPT is reconstructed on the combined shared embedding and
# then visualized by species; coordinates are never independently percentile-ranked by species.

# %%
adata_heatmap = sc.read_h5ad(DPT_OUTPUT_PATH)
if 'lognorm' in adata_heatmap.layers:
    adata_heatmap.X = adata_heatmap.layers['lognorm'].copy()

global_group_order = [
    group for family in ('PT', 'AL', 'DCT', 'CNT_CD')
    for group in HEATMAP_MARKERS[family]
]
global_panel = {
    group: values
    for family in ('PT', 'AL', 'DCT', 'CNT_CD')
    for group, values in _marker_dict(HEATMAP_MARKERS[family]).items()
}
for species in SPECIES_GROUPS:
    species_mask = (
        adata_heatmap.obs['comparison_species'].astype(str).eq(species).to_numpy()
    )
    plot_marker_heatmap(
        adata_heatmap,
        global_panel,
        global_group_order,
        HEATMAP_COLORS,
        pseudotime_col='total_scanpy_dpt',
        cluster_col=None,
        mask=species_mask,
        title=f'{species.title()} nephron markers on the shared DPT',
        output_name=f'{species}_total_marker_heatmap.png',
        strip_col='segment_class',
        strip_order=SEGMENT_DISPLAY_ORDER,
        strip_colors=SEGMENT_STRIP_COLORS,
        n_bins=N_BINS,
        output_dir=HEATMAP_OUTPUT_DIR,
        project_dir=PROJECT_DIR,
    )

# %%
subset_dpt_rows = []
for family, family_markers in HEATMAP_MARKERS.items():
    if family not in set(adata_heatmap.obs['broad_tubule_marker_call'].astype(str)):
        subset_dpt_rows.append({
            'segment': family,
            'status': 'skipped: family absent from reviewed labels',
        })
        continue
    output_col = f'{family.lower()}_subset_dpt'
    try:
        report = recompute_subset_dpt(
            adata_heatmap,
            family,
            _marker_dict(family_markers),
            list(family_markers),
            output_col=output_col,
            n_neighbors=N_NEIGHBORS,
            random_state=RANDOM_STATE,
        )
        report['status'] = 'completed'
        subset_dpt_rows.append(report)
    except (KeyError, ValueError) as exc:
        subset_dpt_rows.append({'segment': family, 'status': f'skipped: {exc}'})
        print(f'WARNING: {family} subset DPT skipped: {exc}')

subset_dpt_report = pd.DataFrame(subset_dpt_rows)
subset_dpt_report.to_csv(HEATMAP_OUTPUT_DIR / 'subset_dpt_report.csv', index=False)
display(subset_dpt_report)

# %%
for family, family_markers in HEATMAP_MARKERS.items():
    pseudotime_col = f'{family.lower()}_subset_dpt'
    if pseudotime_col not in adata_heatmap.obs:
        continue
    for species in SPECIES_GROUPS:
        mask = (
            adata_heatmap.obs['comparison_species'].astype(str).eq(species)
            & adata_heatmap.obs['broad_tubule_marker_call'].astype(str).eq(family)
        ).to_numpy()
        if int(mask.sum()) < 30:
            print(f'Skipping {species} {family} heatmap: only {int(mask.sum())} structures.')
            continue
        plot_marker_heatmap(
            adata_heatmap,
            _marker_dict(family_markers),
            list(family_markers),
            HEATMAP_COLORS,
            pseudotime_col=pseudotime_col,
            cluster_col=None,
            mask=mask,
            title=f'{species.title()} {family} markers on shared {family} DPT',
            output_name=f'{species}_{family.lower()}_subset_heatmap.png',
            strip_col='sample',
            strip_order=SAMPLE_ORDER,
            n_bins=min(N_BINS, 80),
            output_dir=HEATMAP_OUTPUT_DIR,
            project_dir=PROJECT_DIR,
        )

# Persist family coordinates on the canonical DPT object.
adata_heatmap.write(DPT_OUTPUT_PATH)

# %% [markdown]
# # Section 4 - human versus healthy-mouse PT trajectories
#
# The same nested level/shape GAM used in the mouse notebook is fitted on PT structures over the
# mouse/human common DPT support. Here c=0 is mouse and c=1 is human. The human cortex and medulla
# samples come from one donor, so the fitted differences are descriptive effect sizes only.
# Human regions are shown separately as sensitivity curves but are not treated as independent
# biological replicates.

# %%
from pseudospace.levelshape import fit_single_condition_curves, run_level_shape
from pseudospace.stats_gam import as_csr, gam_internal_knots, resolve_present

SECTION4_CONFIG = {
    'n_internal_knots': 9,
    'lambda_grid': GAM_LAMBDA_GRID,
    'common_support_pct': (1, 99),
    'grid_points': 101,
    'min_detected_fraction': 0.02,
    'min_mean_expression': 0.0,
    'pathway_min_genes': 10,
    'pathway_max_genes': 150,
    'pathway_libraries': [
        'Reactome_2022', 'MSigDB_Hallmark_2020', 'KEGG_2019_Mouse',
    ],
    'top_gene_plots': 12,
    'top_pathway_plots': 12,
    'heatmap_genes': 40,
}
SPECIES_COLORS = {'mouse': '#0072B2', 'human': '#D55E00'}
figure_index = []


def _grid_axes(n_panels: int, ncols: int = 3, panel_height: float = 2.8):
    nrows = max(1, int(np.ceil(n_panels / ncols)))
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(4.2 * ncols, panel_height * nrows),
        sharex=True, squeeze=False,
    )
    return fig, axes.ravel()


def _save_figure(fig, filename: str):
    path = CURVE_OUTPUT_DIR / filename
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.show()
    plt.close(fig)
    figure_index.append(filename)


# %%
adata_full = sc.read_h5ad(DPT_OUTPUT_PATH)
if 'lognorm' not in adata_full.layers:
    raise KeyError(f'{DPT_OUTPUT_PATH.name} is missing layers["lognorm"].')
adata_full.X = adata_full.layers['lognorm'].copy()

unexpected_classes = sorted(
    set(adata_full.obs['broad_tubule_marker_call'].astype(str))
    - set(KEEP_TUBULE_CLASSES)
)
if unexpected_classes:
    raise ValueError(f'DPT file contains excluded classes: {unexpected_classes}')

pt_mask = (
    adata_full.obs['broad_tubule_marker_call'].astype(str).eq('PT')
    & adata_full.obs['comparison_species'].astype(str).isin(SPECIES_GROUPS)
    & np.isfinite(adata_full.obs['total_scanpy_dpt'].to_numpy(dtype=float))
).to_numpy()
adata_pt = adata_full[pt_mask].copy()
if adata_pt.n_obs < 100 or adata_pt.obs['comparison_species'].nunique() != 2:
    raise ValueError('Insufficient PT structures from both species.')

s = adata_pt.obs['total_scanpy_dpt'].to_numpy(dtype=float)
species = adata_pt.obs['comparison_species'].astype(str).to_numpy()
c = (species == 'human').astype(float)
samples = adata_pt.obs['sample'].astype(str).to_numpy()
adata_pt.obs['inference_unit'] = np.where(
    species == 'human', 'HUK1', samples
)

p_lo, p_hi = SECTION4_CONFIG['common_support_pct']
lo = max(np.percentile(s[species == group], p_lo) for group in SPECIES_GROUPS)
hi = min(np.percentile(s[species == group], p_hi) for group in SPECIES_GROUPS)
if not lo < hi:
    raise ValueError(f'Mouse and human PT have no common DPT support: {lo:.4f} >= {hi:.4f}.')
grid = np.linspace(lo, hi, SECTION4_CONFIG['grid_points'])
knots = gam_internal_knots(
    s, basis_df=3 + SECTION4_CONFIG['n_internal_knots']
)

Y_all = as_csr(adata_pt.layers['lognorm'])
detected = np.asarray((Y_all > 0).sum(axis=0)).ravel()
gene_mean_all = np.asarray(Y_all.mean(axis=0)).ravel()
min_detected = int(np.ceil(
    SECTION4_CONFIG['min_detected_fraction'] * adata_pt.n_obs
))
tested = (
    (detected >= min_detected)
    & (gene_mean_all >= SECTION4_CONFIG['min_mean_expression'])
)
gene_names = adata_pt.var_names.to_numpy()[tested]
Y_genes = Y_all[:, tested].tocsr().astype(np.float64)
gene_lookup_tested = {str(gene).upper(): index for index, gene in enumerate(gene_names)}

pt_cohort_summary = (
    adata_pt.obs.groupby(
        ['comparison_species', 'sample', 'region', 'inference_unit'], observed=True
    )
    .agg(
        n_pt_structures=('sample', 'size'),
        dpt_min=('total_scanpy_dpt', 'min'),
        dpt_max=('total_scanpy_dpt', 'max'),
    )
    .reset_index()
)
pt_cohort_summary.to_csv(CURVE_OUTPUT_DIR / 'pt_cohort_summary.csv', index=False)
display(pt_cohort_summary.round(3))
print(f'PT structures: {adata_pt.n_obs:,}; tested genes: {len(gene_names):,}')
print(f'Common mouse/human DPT support: [{lo:.3f}, {hi:.3f}]')
print('Inference units: two mouse specimens versus one human donor; no species p-values.')

# %% [markdown]
# ## 4.1 - gene-level level/shape decomposition
#
# M0 is one shared smooth, M1 adds a constant species level shift, and M2 adds a
# species-by-pseudospace interaction. Species-effect RMS summarizes total fitted separation;
# shape RMS summarizes non-parallel remodeling after removing the mean level gap.

# %%
gene_fit = run_level_shape(
    Y_genes, s, c, knots, grid, SECTION4_CONFIG['lambda_grid']
)
axis_basis = {
    gene.upper()
    for role in TOTAL_POSITION_MARKERS.values()
    for gene in role
}
gene_results = pd.DataFrame({
    'gene': gene_names,
    'species_effect_rms': gene_fit['condition_effect_rms'],
    'species_effect_max_abs': gene_fit['condition_effect_max_abs'],
    'level_effect_human_minus_mouse': gene_fit['level_effect'],
    'shape_rms': gene_fit['shape_rms'],
    'curve_spearman': gene_fit['curve_spearman'],
    'mouse_amplitude': gene_fit['amplitude_healthy'],
    'human_amplitude': gene_fit['amplitude_aki'],
    'level_F_cellwise_uncalibrated': gene_fit['level_F'],
    'shape_F_cellwise_uncalibrated': gene_fit['shape_F'],
})
gene_results['axis_basis_gene'] = gene_results['gene'].str.upper().isin(axis_basis)
gene_results['technical_gene'] = gene_results['gene'].str.match(
    r'^(mt-|Rpl|Rps|Mrpl|Mrps)', case=False
)
gene_results['trajectory_summary'] = np.select(
    [
        gene_results['curve_spearman'] >= 0.70,
        gene_results['curve_spearman'] < 0,
    ],
    ['similar fitted shape', 'opposite fitted directions'],
    default='dissimilar fitted shape',
)
gene_results['primary_eligible'] = (
    ~gene_results['axis_basis_gene'] & ~gene_results['technical_gene']
)
gene_results = gene_results.sort_values(
    ['primary_eligible', 'species_effect_rms'],
    ascending=[False, False],
).reset_index(drop=True)
gene_results.to_csv(
    CURVE_OUTPUT_DIR / 'gene_trajectory_comparison.csv', index=False
)
display(gene_results.head(20)[[
    'gene', 'species_effect_rms', 'level_effect_human_minus_mouse',
    'shape_rms', 'curve_spearman', 'axis_basis_gene',
]])

# %%
top_gene_rows = gene_results[gene_results['primary_eligible']].head(
    SECTION4_CONFIG['top_gene_plots']
)
top_genes = top_gene_rows['gene'].tolist()
top_gene_idx = np.array(
    [gene_lookup_tested[gene.upper()] for gene in top_genes], dtype=int
)

sample_gene_curves = {}
for sample_name in sorted(set(samples)):
    sample_mask = samples == sample_name
    curves, support_mask = fit_single_condition_curves(
        Y_genes[sample_mask][:, top_gene_idx],
        s[sample_mask],
        knots,
        grid,
        SECTION4_CONFIG['lambda_grid'],
        gene_fit['lam_idx'][top_gene_idx],
        support_pct=SECTION4_CONFIG['common_support_pct'],
    )
    sample_gene_curves[sample_name] = curves

curve_rows = []
for gene_position, (gene, gene_index) in enumerate(zip(top_genes, top_gene_idx)):
    for group, curves in (
        ('mouse', gene_fit['curve_healthy']),
        ('human', gene_fit['curve_aki']),
    ):
        for x_value, expression in zip(grid, curves[gene_index]):
            curve_rows.append({
                'gene': gene, 'species': group,
                'pseudospace': x_value, 'fitted_lognorm': expression,
            })
pd.DataFrame(curve_rows).to_csv(
    CURVE_OUTPUT_DIR / 'top_gene_fitted_curves.csv', index=False
)

sample_curve_rows = []
for sample_name, curves in sample_gene_curves.items():
    sample_species = str(
        adata_pt.obs.loc[adata_pt.obs['sample'].astype(str).eq(sample_name),
                         'comparison_species'].iloc[0]
    )
    for gene_position, gene in enumerate(top_genes):
        for x_value, expression in zip(grid, curves[gene_position]):
            sample_curve_rows.append({
                'gene': gene, 'sample': sample_name, 'species': sample_species,
                'pseudospace': x_value, 'fitted_lognorm': expression,
            })
pd.DataFrame(sample_curve_rows).to_csv(
    CURVE_OUTPUT_DIR / 'top_gene_per_sample_curves.csv', index=False
)

# %%
fig, axes = _grid_axes(len(top_genes), ncols=3)
for axis, gene, gene_index, gene_position in zip(
    axes, top_genes, top_gene_idx, range(len(top_genes))
):
    for sample_name, curves in sample_gene_curves.items():
        sample_species = (
            'human' if sample_name in HUMAN_SAMPLES else 'mouse'
        )
        axis.plot(
            grid, curves[gene_position], ls='--', lw=0.9, alpha=0.55,
            color=SPECIES_COLORS[sample_species],
        )
    axis.plot(
        grid, gene_fit['curve_healthy'][gene_index], lw=2.2,
        color=SPECIES_COLORS['mouse'], label='mouse',
    )
    axis.plot(
        grid, gene_fit['curve_aki'][gene_index], lw=2.2,
        color=SPECIES_COLORS['human'], label='human',
    )
    axis.set_title(gene)
    axis.set_xlabel('Shared PT DPT')
    axis.set_ylabel('Fitted log-normalized expression')
for axis in axes[len(top_genes):]:
    axis.set_visible(False)
axes[0].legend(frameon=False)
fig.suptitle(
    'Top descriptive mouse-human PT trajectory differences\n'
    'Dashed lines are specimens/regions; human regions share one donor',
    fontsize=12,
)
_save_figure(fig, 'top_gene_trajectory_curves.png')

# %%
# Paired standardized fitted-curve heatmap for the broader collaborator shortlist.
heatmap_rows = gene_results[gene_results['primary_eligible']].head(
    SECTION4_CONFIG['heatmap_genes']
)
heatmap_idx = np.array(
    [gene_lookup_tested[gene.upper()] for gene in heatmap_rows['gene']], dtype=int
)
mouse_heat = gene_fit['curve_healthy'][heatmap_idx]
human_heat = gene_fit['curve_aki'][heatmap_idx]
pooled = np.concatenate([mouse_heat, human_heat], axis=1)
row_mean = pooled.mean(axis=1, keepdims=True)
row_std = pooled.std(axis=1, keepdims=True)
row_std[row_std < 1e-8] = 1.0
mouse_z = np.clip((mouse_heat - row_mean) / row_std, -2.5, 2.5)
human_z = np.clip((human_heat - row_mean) / row_std, -2.5, 2.5)

fig, axes = plt.subplots(
    1, 2, figsize=(9.5, max(6, 0.19 * len(heatmap_rows))),
    sharey=True,
)
for axis, values, title in zip(
    axes, (mouse_z, human_z), ('Healthy mouse', 'Human')
):
    image = axis.imshow(
        values, aspect='auto', interpolation='nearest', cmap='bwr',
        vmin=-2.5, vmax=2.5,
        extent=[grid[0], grid[-1], len(heatmap_rows) - 0.5, -0.5],
    )
    axis.set_title(title)
    axis.set_xlabel('Shared PT DPT')
axes[0].set_yticks(np.arange(len(heatmap_rows)))
axes[0].set_yticklabels(heatmap_rows['gene'])
fig.colorbar(image, ax=axes, label='Pooled gene-wise z-score', shrink=0.65)
fig.suptitle('Paired fitted PT trajectories ranked by mouse-human effect size')
_save_figure(fig, 'paired_top_gene_heatmap.png')

# %% [markdown]
# ## 4.2 - canonical PT marker-gradient check
#
# These axis-adjacent genes are displayed separately from the discovery ranking. Their purpose is
# to show whether S1, S2, and S3 positional programs occupy comparable locations in the shared
# coordinate, not to claim species differences.

# %%
pt_marker_groups = {
    key.replace('PT-', ''): genes
    for key, genes in HEATMAP_MARKERS['PT'].items()
}
fig, axes = plt.subplots(1, 3, figsize=(12, 3.4), sharex=True)
marker_curve_rows = []
for axis, (group, requested) in zip(axes, pt_marker_groups.items()):
    present = [
        gene for gene in resolve_present(requested, gene_names)
        if gene.upper() in gene_lookup_tested
    ]
    for gene in present:
        gene_index = gene_lookup_tested[gene.upper()]
        axis.plot(
            grid, gene_fit['curve_healthy'][gene_index],
            color=SPECIES_COLORS['mouse'], lw=1.8,
            label='mouse' if gene == present[0] else None,
        )
        axis.plot(
            grid, gene_fit['curve_aki'][gene_index],
            color=SPECIES_COLORS['human'], lw=1.8, ls='--',
            label='human' if gene == present[0] else None,
        )
        for species_name, values in (
            ('mouse', gene_fit['curve_healthy'][gene_index]),
            ('human', gene_fit['curve_aki'][gene_index]),
        ):
            for x_value, expression in zip(grid, values):
                marker_curve_rows.append({
                    'marker_group': group, 'gene': gene,
                    'species': species_name, 'pseudospace': x_value,
                    'fitted_lognorm': expression,
                })
    axis.set_title(f'{group}: {", ".join(present)}')
    axis.set_xlabel('Shared PT DPT')
    axis.set_ylabel('Fitted log-normalized expression')
axes[0].legend(frameon=False)
fig.suptitle('Canonical PT marker gradients on the shared coordinate')
_save_figure(fig, 'canonical_pt_marker_curves.png')
pd.DataFrame(marker_curve_rows).to_csv(
    CURVE_OUTPUT_DIR / 'canonical_pt_marker_curves.csv', index=False
)

# %% [markdown]
# ## 4.3 - predefined pathway/module trajectories
#
# Pathway scores are the mean of PT-standardized member-gene expression. All eligible predefined
# Hallmark, KEGG Mouse, and Reactome sets are modeled; pathways are ranked by descriptive
# mouse-human fitted-curve separation.

# %%
gene_mean = np.asarray(Y_genes.mean(axis=0)).ravel()
gene_sq = np.asarray(Y_genes.multiply(Y_genes).mean(axis=0)).ravel()
gene_std = np.sqrt(np.maximum(gene_sq - gene_mean ** 2, 1e-12))
background_lookup = {gene.upper(): gene for gene in gene_names}
pathway_rows = []
pathway_mapping_rows = []

for library in SECTION4_CONFIG['pathway_libraries']:
    library_path = PATHWAY_LIBRARY_DIR / f'{library}.json'
    if not library_path.exists():
        print(f'WARNING: missing pathway library {library_path.name}')
        continue
    library_sets = json.loads(library_path.read_text())
    requested_symbols = set()
    matched_symbols = set()
    kept_count = 0
    for pathway, members in library_sets.items():
        if isinstance(members, dict):
            members = members.get('genes', [])
        requested = {str(gene).upper() for gene in members}
        present = sorted({
            background_lookup[gene]
            for gene in requested
            if gene in background_lookup
        })
        requested_symbols |= requested
        matched_symbols |= {gene for gene in requested if gene in background_lookup}
        if (
            SECTION4_CONFIG['pathway_min_genes']
            <= len(present)
            <= SECTION4_CONFIG['pathway_max_genes']
        ):
            pathway_rows.append({
                'library': library,
                'pathway': pathway,
                'genes_present': present,
                'n_genes_present': len(present),
            })
            kept_count += 1
    pathway_mapping_rows.append({
        'library': library,
        'n_pathways_retained': kept_count,
        'unique_symbols_requested': len(requested_symbols),
        'unique_symbols_matched': len(matched_symbols),
        'symbol_match_rate': (
            len(matched_symbols) / max(len(requested_symbols), 1)
        ),
    })

retained_pathways = pd.DataFrame(pathway_rows)
pathway_mapping = pd.DataFrame(pathway_mapping_rows)
pathway_mapping.to_csv(
    CURVE_OUTPUT_DIR / 'pathway_symbol_mapping_report.csv', index=False
)
display(pathway_mapping.round(3))

if retained_pathways.empty:
    pathway_results = pd.DataFrame()
    pathway_module_scores = np.empty((adata_pt.n_obs, 0))
else:
    pathway_module_scores = np.empty(
        (adata_pt.n_obs, len(retained_pathways)), dtype=np.float64
    )
    Y_csc = Y_genes.tocsc()
    for pathway_index, members in enumerate(retained_pathways['genes_present']):
        indices = [gene_lookup_tested[gene.upper()] for gene in members]
        values = Y_csc[:, indices].toarray()
        pathway_module_scores[:, pathway_index] = (
            (values - gene_mean[indices]) / gene_std[indices]
        ).mean(axis=1)

    pathway_fit = run_level_shape(
        pathway_module_scores,
        s,
        c,
        knots,
        grid,
        SECTION4_CONFIG['lambda_grid'],
    )
    pathway_results = retained_pathways.drop(columns='genes_present').copy()
    pathway_results['genes_present'] = retained_pathways['genes_present'].map(
        lambda values: '; '.join(values)
    )
    pathway_results['species_effect_rms'] = pathway_fit['condition_effect_rms']
    pathway_results['level_effect_human_minus_mouse'] = pathway_fit['level_effect']
    pathway_results['shape_rms'] = pathway_fit['shape_rms']
    pathway_results['curve_spearman'] = pathway_fit['curve_spearman']
    pathway_results = pathway_results.sort_values(
        'species_effect_rms', ascending=False
    ).reset_index(drop=True)
pathway_results.to_csv(
    CURVE_OUTPUT_DIR / 'pathway_trajectory_comparison.csv', index=False
)
display(pathway_results.head(20) if len(pathway_results) else pathway_results)

# %%
if len(pathway_results):
    top_pathway_rows = pathway_results.head(
        SECTION4_CONFIG['top_pathway_plots']
    )
    pathway_key_to_index = {
        (row.library, row.pathway): index
        for index, row in retained_pathways.iterrows()
    }
    top_pathway_idx = np.array([
        pathway_key_to_index[(row.library, row.pathway)]
        for row in top_pathway_rows.itertuples()
    ], dtype=int)

    sample_pathway_curves = {}
    for sample_name in sorted(set(samples)):
        sample_mask = samples == sample_name
        curves, _ = fit_single_condition_curves(
            pathway_module_scores[sample_mask][:, top_pathway_idx],
            s[sample_mask],
            knots,
            grid,
            SECTION4_CONFIG['lambda_grid'],
            pathway_fit['lam_idx'][top_pathway_idx],
            support_pct=SECTION4_CONFIG['common_support_pct'],
        )
        sample_pathway_curves[sample_name] = curves

    fig, axes = _grid_axes(len(top_pathway_idx), ncols=3)
    pathway_source_rows = []
    for plot_index, (axis, pathway_index, row) in enumerate(zip(
        axes, top_pathway_idx, top_pathway_rows.itertuples()
    )):
        for sample_name, curves in sample_pathway_curves.items():
            sample_species = 'human' if sample_name in HUMAN_SAMPLES else 'mouse'
            axis.plot(
                grid, curves[plot_index], ls='--', lw=0.8, alpha=0.5,
                color=SPECIES_COLORS[sample_species],
            )
        for species_name, values in (
            ('mouse', pathway_fit['curve_healthy'][pathway_index]),
            ('human', pathway_fit['curve_aki'][pathway_index]),
        ):
            axis.plot(
                grid, values, lw=2.1, color=SPECIES_COLORS[species_name],
                label=species_name,
            )
            for x_value, score in zip(grid, values):
                pathway_source_rows.append({
                    'library': row.library,
                    'pathway': row.pathway,
                    'species': species_name,
                    'pseudospace': x_value,
                    'fitted_module_z_score': score,
                })
        axis.set_title(row.pathway, fontsize=8)
        axis.set_xlabel('Shared PT DPT')
        axis.set_ylabel('Fitted module z-score')
    for axis in axes[len(top_pathway_idx):]:
        axis.set_visible(False)
    axes[0].legend(frameon=False)
    fig.suptitle(
        'Top descriptive mouse-human pathway trajectory differences\n'
        'Dashed lines are specimens/regions',
        fontsize=12,
    )
    _save_figure(fig, 'top_pathway_trajectory_curves.png')
    pd.DataFrame(pathway_source_rows).to_csv(
        CURVE_OUTPUT_DIR / 'top_pathway_fitted_curves.csv', index=False
    )

# %% [markdown]
# # Section 5 - three-axis concordance
#
# The marker axis and shared DPT are compared with distance to the nearest reviewed glomerular
# structure within each sample. This spatial proxy is diagnostic only. It is computed separately
# per sample because slide coordinates are not comparable across specimens. In the current
# reviewed clustering no cluster has sufficiently specific glomerular identity, so this section
# records a transparent skip rather than using the mixed cluster 7 as a false spatial anchor.

# %%
adata_pass1 = sc.read_h5ad(HARMONY_OUTPUT_PATH)
dpt_for_spatial = sc.read_h5ad(DPT_OUTPUT_PATH)
required_coordinates = {'x_centroid', 'y_centroid'}
concordance_rows = []

if required_coordinates.issubset(adata_pass1.obs.columns):
    for sample_name in SAMPLE_ORDER:
        pass1_sample = adata_pass1.obs['sample'].astype(str).eq(sample_name).to_numpy()
        glomerulus = (
            pass1_sample
            & adata_pass1.obs['coarse_class'].astype(str).eq('Glomerulus').to_numpy()
        )
        dpt_sample = dpt_for_spatial.obs['sample'].astype(str).eq(sample_name).to_numpy()
        if glomerulus.sum() < 2 or dpt_sample.sum() < 10:
            concordance_rows.append({
                'sample': sample_name,
                'status': 'skipped: insufficient glomerular anchors or tubules',
            })
            continue

        anchor_xy = adata_pass1.obs.loc[
            glomerulus, ['x_centroid', 'y_centroid']
        ].to_numpy(dtype=float)
        query_xy = dpt_for_spatial.obs.loc[
            dpt_sample, ['x_centroid', 'y_centroid']
        ].to_numpy(dtype=float)
        valid_anchor = np.isfinite(anchor_xy).all(axis=1)
        valid_query = np.isfinite(query_xy).all(axis=1)
        distances = np.full(query_xy.shape[0], np.nan)
        distances[valid_query] = cKDTree(anchor_xy[valid_anchor]).query(
            query_xy[valid_query], k=1
        )[0]

        for axis_name in ('total_marker_axis', 'total_scanpy_dpt'):
            values = dpt_for_spatial.obs.loc[
                dpt_sample, axis_name
            ].to_numpy(dtype=float)
            finite = np.isfinite(values) & np.isfinite(distances)
            rho = (
                spearmanr(values[finite], distances[finite]).correlation
                if finite.sum() >= 10 else np.nan
            )
            concordance_rows.append({
                'sample': sample_name,
                'species': (
                    'human' if sample_name in HUMAN_SAMPLES else 'mouse'
                ),
                'axis': axis_name,
                'physical_proxy': 'distance_to_nearest_glomerulus',
                'n_structures': int(finite.sum()),
                'n_glomerular_anchors': int(valid_anchor.sum()),
                'spearman': float(rho),
                'status': 'completed',
            })
else:
    concordance_rows.append({
        'status': 'skipped: centroid coordinates absent from pass-1 object',
    })

concordance = pd.DataFrame(concordance_rows)
concordance.to_csv(DIAGNOSTIC_DIR / 'three_axis_concordance.csv', index=False)
display(concordance.round(3))

# %% [markdown]
# # Section 6 - validation and summary
#
# The final cell writes a compact methods/inference note and lists the key outputs. Cross-species
# conclusions must remain donor-descriptive until independent human donors are added.

# %%
analysis_notes = f"""# Human versus healthy-mouse pseudospace

- Shared matrix before cell typing: {adata_combined.n_obs:,} structures x {adata_combined.n_vars:,} ortholog genes.
- Canonical DPT object: {adata_dpt_saved.n_obs:,} retained nephron structures.
- PT comparison: {adata_pt.n_obs:,} structures over shared support [{lo:.4f}, {hi:.4f}].
- Mouse controls: {', '.join(MOUSE_SAMPLES)}.
- Human samples: {', '.join(HUMAN_SAMPLES)}; both come from one donor.
- Ortholog mapping: HCOP reciprocal best matches with support from at least 3 databases; identical symbols are preferred when support ties. No manual ortholog overrides.
- Integration: deterministic R Harmony {HARMONY_EXPECTED_VERSION}; sample is the batch key and species is never regressed.
- Cell labels: one reviewed pass-1 Leiden map pinned to fingerprint {CLUSTERING_FINGERPRINT['membership_sha1']}; no relabeling after pass 2.
- Species comparisons: descriptive level, shape, and fitted-curve effect sizes only.
- Human cortex/medulla sampling and donor identity are inseparable from species in this cohort.
- Human cortex and medulla curves are sensitivity views, not independent biological replicates.
- No cluster has sufficiently specific glomerular identity; distance-to-glomerulus concordance is skipped.
- Columns ending cellwise_uncalibrated are diagnostics and must not be interpreted as species tests.
"""
(RESULTS_DIR / 'analysis_notes.md').write_text(analysis_notes)

print('Cross-species workflow complete.')
print(f'Results: {RESULTS_DIR.relative_to(PROJECT_DIR)}')
print('Figures:')
for filename in figure_index:
    print(f'  - {filename}')
print('Top descriptive gene effects:')
display(gene_results.head(10)[[
    'gene', 'species_effect_rms', 'shape_rms',
    'level_effect_human_minus_mouse', 'curve_spearman',
]])
