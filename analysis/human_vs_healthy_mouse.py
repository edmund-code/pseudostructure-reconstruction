# %% [markdown]
# # Human versus healthy-mouse exploratory classification and PT pseudospace
#
# Broad identities, structure QC, species-specific validation, and integration diagnostics precede an exploratory PT trajectory. Fine segment labels and a shared whole-nephron axis are not validated. Both human slices are healthy cortex from one donor.
#

# %% [markdown]
# # Section 0 - configuration
#
# Paths, deterministic integration parameters, support-based HCOP mapping, cell-typing panels,
# label vocabulary, and heatmap panels. All downstream code is defined in this notebook.
#
#

# %%
# Purpose: from __future__ import annotations
# %load_ext autoreload
# %autoreload 2

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

# Define only locations needed to start the workflow.
TUBULE_BY_GENE_DIR = DATA_ROOT / 'tubule_by_gene'
ORTHOLOG_TABLE_PATH = DATA_ROOT / 'human_mouse_hcop_fifteen_column.txt.gz'

RESULTS_DIR = RESULTS_ROOT / 'human_vs_healthy_mouse'
HARMONY_OUTPUT_PATH = RESULTS_DIR / 'cross_species_harmony_pass1.h5ad'
PASS2_HARMONY_OUTPUT_PATH = RESULTS_DIR / 'cross_species_harmony_pass2.h5ad'
GLOBAL_DPT_OUTPUT_PATH = RESULTS_DIR / 'cross_species_nephron_global_dpt.h5ad'
# The PT file below is a post-global-DPT subset for later PT analyses.
DPT_OUTPUT_PATH = RESULTS_DIR / 'cross_species_pt_dpt.h5ad'
CELLTYPING_DIR = RESULTS_DIR / 'celltyping'
HEATMAP_OUTPUT_DIR = RESULTS_DIR / 'heatmaps'
CURVE_OUTPUT_DIR = RESULTS_DIR / 'curves'
DIAGNOSTIC_DIR = RESULTS_DIR / 'diagnostics'
for _directory in (RESULTS_DIR, CELLTYPING_DIR, HEATMAP_OUTPUT_DIR, CURVE_OUTPUT_DIR,
                   DIAGNOSTIC_DIR):
    _directory.mkdir(parents=True, exist_ok=True)



print(f'Project: {PROJECT_DIR.name}')
print(f'Results will be written to: {RESULTS_DIR.relative_to(PROJECT_DIR)}')


# %% [markdown]
# # Section 1 - pass-1 integration and independent diagnostics
#
# Human genes are mapped to reviewed mouse orthologs. The embedding uses intersected species HVGs. Sample correction can remove species and region biology because sample, species, and human region are confounded. Mixing is evaluated within broad compartments and is not optimized in isolation.
#
# **Navigation:** settings → data/QC → uncorrected reference → species-first validation → Harmony. Each code cell produces one table, figure, saved object, or clearly named state transition.
#

# %% [markdown]
# ## 1.1 - pass-1 settings and dependencies
#
# **Purpose:** define only the parameters required for the following integration pass.
#

# %%
# Purpose: # Settings used only for pass-1 loading, QC, and integration.
# Settings used only for pass-1 loading, QC, and integration.
# Keep these adjacent to the pass that consumes them so reruns are easy to audit.
MOUSE_SAMPLES = ('Ctrl1A2', 'Ctrl1A4')
HUMAN_SAMPLES = ('HUK1_COR1', 'HUK1_MED1')
SAMPLE_ORDER = (*MOUSE_SAMPLES, *HUMAN_SAMPLES)
# Both human slices are healthy CORTEX. `HUK1_MED1` is only *named* medulla at source (the Visium
# matrix is HUK1_MED); the processed segmentation is cortex tissue, so it must never be treated as
# medullary in an analysis or a write-up. `sample` keeps the source name for provenance.
HUMAN_REGIONS = {'HUK1_COR1': 'cortex', 'HUK1_MED1': 'cortex'}
SPECIES_GROUPS = ('mouse', 'human')
BATCH_KEY = 'sample'
RANDOM_STATE = 0

MIN_GENES_PER_TUBULE = 100
# This gene filter is applied immediately after structure QC, before HVG selection
# and Harmony, so low-support genes cannot influence the integration embedding.
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
COARSE_RESOLUTION = 0.7

print('Pass-1 settings:', {'min_genes': MIN_GENES_PER_TUBULE,
                            'min_gene_fraction': MIN_GENE_TUBULE_FRACTION,
                            'min_gene_total_counts': MIN_GENE_TOTAL_COUNTS, 'HVGs': N_HVGS,
                            'batch_key': BATCH_KEY, 'seed': RANDOM_STATE})


# %%
# Purpose: import scanpy as sc
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
import rpy2.robjects as ro
from IPython.display import display
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
from pseudospace.heatmaps import plot_marker_heatmap
from pseudospace.markers import build_gene_lookup
from pseudospace.trajectory import (
    choose_diffusion_component,
    choose_root_global,
    choose_root_pt_cluster,
    orient_and_normalize,
    recompute_subset_dpt,
    save_total_pseudotime_anndata,
)

import matplotlib as _mpl
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
# Purpose: if not ORTHOLOG_TABLE_PATH.exists():
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
# Purpose: audit QC retention before structures are filtered.
# `shared_n_genes` is the number of detected genes in the mapped, shared human–mouse
# ortholog matrix for each structure; it is not the number of genes in the native input.
# A structure passes only when `shared_n_genes >= MIN_GENES_PER_TUBULE` (currently 100).
# Within each `(sample, region)` group, `fraction_removed = 1 - mean(passes_min_genes)`.
# Equivalently: `(number before filtering - number passing) / number before filtering`.

qc_before = adata_combined.obs.copy()
qc_before['shared_n_genes'] = np.asarray((adata_combined.X > 0).sum(axis=1)).ravel()
qc_before['passes_min_genes'] = qc_before['shared_n_genes'] >= MIN_GENES_PER_TUBULE
qc_before['shared_total_counts'] = np.asarray(adata_combined.X.sum(axis=1)).ravel()
qc_before.to_csv(DIAGNOSTIC_DIR / 'structure_qc_before_filter.csv')
qc_retention = qc_before.groupby(['sample', 'region'], observed=True)['passes_min_genes'].agg(['size', 'sum', 'mean'])
qc_retention['fraction_removed'] = 1 - qc_retention['mean']
qc_retention.to_csv(DIAGNOSTIC_DIR / 'qc_retention_by_sample.csv')
# Start a structure-level audit. `sum` is the number retained because True counts as 1.
tubule_filter_audit = pd.DataFrame([{
    'filter_stage': 'shared-ortholog minimum-gene QC (>= 100 detected genes)',
    'filter_type': 'tubule filter',
    'n_input_tubules': int(len(qc_before)),
    'n_retained_tubules': int(qc_before['passes_min_genes'].sum()),
    'n_removed_tubules': int((~qc_before['passes_min_genes']).sum()),
}])
tubule_filter_audit['fraction_removed'] = (
    tubule_filter_audit['n_removed_tubules'] / tubule_filter_audit['n_input_tubules']
)
tubule_filter_audit.to_csv(DIAGNOSTIC_DIR / 'tubule_filter_audit.csv', index=False)
display(qc_retention)
display(tubule_filter_audit)
qc_metrics = [c for c in qc_before if any(k in c.lower() for k in ('count', 'gene', 'mito', 'pct_counts_mt', 'area')) and pd.api.types.is_numeric_dtype(qc_before[c]) and c != 'passes_min_genes']
qc_before.groupby('sample', observed=True)[qc_metrics].describe().to_csv(DIAGNOSTIC_DIR / 'qc_distribution_by_sample.csv')
print('Area fields:', [c for c in qc_metrics if 'area' in c.lower()])
print('Boundary fields available:', [c for c in qc_before if any(k in c.lower() for k in ('boundary', 'polygon', 'contour'))])
print('Centroid plots below locate structures; segmentation boundaries require source-image review.')


# %%
# Purpose: apply structure QC, then low-support gene QC before Harmony HVG selection.
print(f'Before structure filtering: {adata_combined.n_obs:,}')
sc.pp.filter_cells(adata_combined, min_genes=MIN_GENES_PER_TUBULE)
print(f'After structure filtering:  {adata_combined.n_obs:,}')

# The same measurement once the tubule filter has run, for the before/after plot below. Taken
# here, before the gene filter drops columns, so it reflects the tubule filter alone.
_qc_after = pd.DataFrame({
    'sample': adata_combined.obs['sample'].astype(str).to_numpy(),
    'shared_n_genes': np.asarray((adata_combined.X > 0).sum(axis=1)).ravel(),
    'total_counts': np.asarray(adata_combined.X.sum(axis=1)).ravel(),
}, index=adata_combined.obs_names)

adata_combined.layers['counts'] = adata_combined.X.copy()
expressed = np.asarray((adata_combined.layers['counts'] > 0).sum(axis=0)).ravel()
total_counts = np.asarray(adata_combined.layers['counts'].sum(axis=0)).ravel()
min_cells = max(1, int(np.ceil(MIN_GENE_TUBULE_FRACTION * adata_combined.n_obs)))
gene_keep = (expressed >= min_cells) & (total_counts >= MIN_GENE_TOTAL_COUNTS)
adata_combined.var['n_structures_expressed'] = expressed
adata_combined.var['total_counts'] = total_counts
adata_combined.var['passes_expression_count_filter'] = gene_keep
n_genes_before_expression_filter = int(adata_combined.n_vars)
# The normalisation below divides by the RETAINED panel's total, so a structure that lost genes to
# the expression filter is scaled up relative to one that did not. Record the native (pre-filter)
# total here so Section 4 can repeat the relative-expression comparison with that denominator.
pre_filter_total_counts = np.asarray(adata_combined.layers['counts'].sum(axis=1)).ravel()
adata_combined = adata_combined[:, gene_keep].copy()
tubule_filter_audit.loc[len(tubule_filter_audit)] = {
    'filter_stage': 'expression-count filter before Harmony (genes only)',
    'filter_type': 'gene filter; no structures removed',
    'n_input_tubules': int(adata_combined.n_obs),
    'n_retained_tubules': int(adata_combined.n_obs),
    'n_removed_tubules': 0,
    'fraction_removed': 0.0,
}
tubule_filter_audit.to_csv(DIAGNOSTIC_DIR / 'tubule_filter_audit.csv', index=False)
print(f'Gene filter before Harmony: {n_genes_before_expression_filter:,} → {adata_combined.n_vars:,} genes; '
      f'minimum support = {min_cells:,} structures and {MIN_GENE_TOTAL_COUNTS:,} total counts.')
sc.pp.normalize_total(adata_combined, target_sum=NORMALIZE_TARGET_SUM)
sc.pp.log1p(adata_combined)
adata_combined.layers['lognorm'] = adata_combined.X.copy()

# Denominator audit: how much of each structure's native UMI survives ortholog mapping and the gene
# expression filter, i.e. what the normalisation actually divides by.
adata_combined.obs['pre_filter_total_counts'] = pre_filter_total_counts
adata_combined.obs['retained_panel_total_counts'] = np.asarray(
    adata_combined.layers['counts'].sum(axis=1)
).ravel()
adata_combined.obs['retained_panel_fraction'] = (
    adata_combined.obs['retained_panel_total_counts']
    / adata_combined.obs['pre_filter_total_counts'].clip(lower=1)
)
normalisation_audit = (
    adata_combined.obs.groupby(['comparison_species', 'sample'], observed=True)['retained_panel_fraction']
    .agg(['size', 'mean', 'median', 'min']).reset_index()
)
normalisation_audit.to_csv(DIAGNOSTIC_DIR / 'normalisation_denominator_audit.csv', index=False)
display(normalisation_audit.round(3))
print('The retained panel is the normalisation denominator; Section 4 repeats the comparison with '
      'the native total (normalisation_sensitivity.csv).')

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
OUTPUT_DIR = DIAGNOSTIC_DIR

# --- Per-structure size around the tubule filter: before and after, per sample ------------------
# This is the distribution `sc.pp.filter_cells(min_genes=MIN_GENES_PER_TUBULE)` acts on, and what it
# left behind. `shared_n_genes` is the ortholog-space count, not the native input count, so a
# sample cut disproportionately by the human->mouse mapping is visible here rather than only as a
# smaller n later. The pre-filter frame comes from the QC audit two cells up; the post-filter frame
# was captured in the filter cell, before the gene filter dropped columns.
_before = qc_before[['sample', 'shared_n_genes']].copy()
_before['sample'] = _before['sample'].astype(str)
_before['total_counts'] = np.asarray(qc_before['shared_total_counts']).ravel()
_before['passes_min_genes'] = _before['shared_n_genes'] >= MIN_GENES_PER_TUBULE
_after = _qc_after
_before.to_csv(OUTPUT_DIR / 'tubule_counts_before_filter.csv')
_after.to_csv(OUTPUT_DIR / 'tubule_counts_after_filter.csv')
_retention = _before.groupby('sample', observed=True)['passes_min_genes'].agg(['size', 'sum', 'mean'])
_retention['fraction_removed'] = 1 - _retention['mean']
print(f"Tubule filter (>= {MIN_GENES_PER_TUBULE} detected genes): "
      f"{int(_before['passes_min_genes'].sum()):,}/{len(_before):,} structures pass, "
      f"{1 - _before['passes_min_genes'].mean():.2%} removed")
print(_retention.to_string())

_samples = sorted(_before['sample'].unique())
_colors = dict(zip(_samples, plt.get_cmap('tab10').colors))
_GENE_COLUMN = 'shared_n_genes' if 'shared_n_genes' in _before.columns else 'n_genes_by_counts'
_NOUN = 'structure' if _GENE_COLUMN == 'shared_n_genes' else 'tubule'
_GENE_AXIS = ('detected genes per structure (shared ortholog space)' if _GENE_COLUMN == 'shared_n_genes'
              else 'detected genes per tubule')

fig, axes = plt.subplots(2, 3, figsize=(21, 10))

for column, (frame, state) in enumerate(((_before, 'before'), (_after, 'after'))):
    for sample in _samples:
        sub = frame[frame['sample'] == sample]
        axes[0, column].hist(sub[_GENE_COLUMN], bins=100, histtype='step', lw=1.2,
                             color=_colors[sample], label=f'{sample} (n={len(sub):,})')
    if state == 'before':
        axes[0, column].axvline(MIN_GENES_PER_TUBULE, color='crimson', ls='--', lw=1.6,
                                label=f'threshold = {MIN_GENES_PER_TUBULE} genes')
    axes[0, column].set_yscale('log')
    axes[0, column].set_xlabel(_GENE_AXIS)
    axes[0, column].set_ylabel(f'{_NOUN}s (log scale)')
    axes[0, column].set_title(f'Detected genes per {_NOUN}, {state} filtering')
    axes[0, column].legend(fontsize=7)

    for sample in _samples:
        sub = frame[frame['sample'] == sample]
        axes[1, column].hist(np.log10(sub['total_counts'].to_numpy() + 1), bins=100,
                             histtype='step', lw=1.2, color=_colors[sample], label=sample)
    axes[1, column].set_xlabel(f'log10(total counts + 1) per {_NOUN}')
    axes[1, column].set_ylabel(f'{_NOUN}s')
    axes[1, column].set_title(f'Total counts per {_NOUN}, {state} filtering')
    axes[1, column].legend(fontsize=7)

# How much each sample lost, and what that did to its size distribution.
_removed = (1 - _retention['mean'].reindex(_samples).to_numpy()) * 100
axes[0, 2].bar(_samples, _removed, color=[_colors[sample] for sample in _samples])
for index, value in enumerate(_removed):
    axes[0, 2].text(index, value, f'{value:.0f}%', ha='center', va='bottom', fontsize=9)
axes[0, 2].set_ylim(0, max(_removed) * 1.25)
axes[0, 2].set_ylabel(f'% of {_NOUN}s removed')
axes[0, 2].set_title('Removed by the tubule filter, per sample')
axes[0, 2].tick_params(axis='x', rotation=30)

_positions = np.arange(len(_samples)) * 3.0
axes[1, 2].boxplot([np.log10(_before.loc[_before['sample'] == sample, 'total_counts'].to_numpy() + 1)
                    for sample in _samples], positions=_positions - 0.6, widths=0.9,
                   showfliers=False, patch_artist=True,
                   boxprops=dict(facecolor='#CCCCCC', edgecolor='#888888'),
                   medianprops=dict(color='black'))
axes[1, 2].boxplot([np.log10(_after.loc[_after['sample'] == sample, 'total_counts'].to_numpy() + 1)
                    for sample in _samples], positions=_positions + 0.6, widths=0.9,
                   showfliers=False, patch_artist=True,
                   boxprops=dict(facecolor='#7FB3D5', edgecolor='#2E6DA4'),
                   medianprops=dict(color='black'))
axes[1, 2].set_xticks(_positions)
axes[1, 2].set_xticklabels(_samples, rotation=30, ha='right')
axes[1, 2].set_ylabel('log10(total counts + 1)')
axes[1, 2].set_title('Size per sample, before (grey) vs after (blue)')

fig.suptitle(f'{_NOUN.capitalize()} size around the tubule filter '
             f'(>= {MIN_GENES_PER_TUBULE} detected genes)', fontsize=14)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(OUTPUT_DIR / 'tubule_count_distribution.png', dpi=200, bbox_inches='tight')
plt.show()


# %%
# Purpose: # Uncorrected PCA/UMAP is retained as the before-Harmony reference.
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


print('Completed:', '# Uncorrected PCA/UMAP is retained as the before-Harmony reference.')


# %%
# Purpose: print(f'Running deterministic R Harmony on {adata_hvg.n_obs:,} structures x '
print(f'Running deterministic R Harmony on {adata_hvg.n_obs:,} structures x '
      f'{adata_hvg.n_vars:,} shared HVGs.')
ro.r(f'set.seed({RANDOM_STATE})')
adata_hvg = run_harmony_rpy2(
    adata_hvg,
    batch_key=BATCH_KEY,  # 'sample': integrates cross-species & intra-species slide batches simultaneously
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
# # Section 2 - cluster review, nephron reintegration, and global DPT
#
# Retain every QC-passing structure in pass 1. Broad identities are manual review decisions based on cluster differential expression. After review, remove only structures outside the tubular nephron continuum, re-harmonize the retained nephron structures, and inspect global DPT. Subset PT only afterward for the downstream comparison.
#
# **Navigation:** filter genes → cluster → inspect DE and spatial context → assign broad labels → retain tubular-nephron classes → re-harmonize → inspect global DPT → subset PT.
#

# %% [markdown]
# ## 2.1 - annotation settings
#
# **Purpose:** set filtering and conservative-label rules before applying them.
#

# %%
# Cluster review settings.
# Labels below are manual review decisions based on cluster-specific differential
# expression tables produced in the next cells. They are not computed from marker scores.
# Re-review this map whenever clustering parameters or the input cohort changes.
# Update these reviewed calls if clustering parameters, data, or cluster IDs change.
# These calls were reviewed from the current run's DE heatmap, reference-gene
# dotplot/evidence panel, and sample composition. They are broad labels only.
REVIEWED_CLUSTER_LABELS = {
    '0': 'PT-S1',        # PT-S1 +2.67, next PT-S2 +1.38
    '1': 'PT-S2',        # PT-S2 +2.39, next PT-S1 +0.94
    '2': 'CNT_CD',       # OMCD +2.13 vs CNT +2.09 -- panels do not separate; family call
    '3': 'CNT_CD',       # IMCD +2.61 vs CCD +1.97, but DE is Aqp2+Aqp3 (CCD) -- family call
    '4': 'DCT1',         # DCT1 +2.42, next DCT2 +1.92; DE Slc12a3/Trpm6
    '5': 'AL',           # cTAL +1.76 vs mTAL +1.46 and 921/955 of the cluster is mouse -- family call
    '6': 'Glomerulus',   # Podocyte +3.14 with Vessel +2.87: mixed glomerular structure
    '7': 'PT-S3',        # PT-S3 +2.73, next +0.06
    '8': 'Unassigned',   # best panel only +0.65 -- no segment identity
    '9': 'SmoothMuscle', # SmoothMuscle +3.13
    '10': 'mTAL',        # mTAL +2.07 at full panel coverage; DE Slc12a1/Umod/Kcnj1/Cldn10
}
ASSIGNMENT_INTERPRETATION = {
    '0': 'PT S1: S1 solute transport and proximal metabolic program (PT-S1 panel best, PT-S2 second).',
    '1': 'PT S2: proximal solute transport (Slc22a6, Slc13a3); the S1 call that preceded the v2 re-segmentation no longer matched.',
    '2': 'Collecting duct: CNT and OMCD panels tie and the DE mixes Calb1/Hsd11b2/Slc8a1 (CNT) with Rhcg (OMCD), so the family is reported.',
    '3': 'Collecting duct: the IMCD panel scores highest but the DE is Aqp2+Aqp3 (CCD-like), so the family is reported rather than a contested fine call.',
    '4': 'DCT: distal convoluted tubule, DCT1 panel best (Slc12a3, Trpm6).',
    '5': 'Ascending limb; 921 of 955 structures are mouse, so this is a mouse-driven cluster and gets the family rather than a thin/thick call.',
    '6': 'Mixed glomerular structure: podocyte genes plus endothelial/perivascular signal (Vessel panel second).',
    '7': 'PT S3: late proximal program (Slc22a7, Slc7a13); four of five PT-S3 panel genes are top DE.',
    '8': 'No segment identity: the best panel reaches only +0.65 and the top DE is a stress/matrix program, so it stays unresolved.',
    '9': 'Smooth-muscle/perivascular program (Acta2, Myh11, Tagln); non-nephron.',
    '10': 'TAL, medullary: thick ascending limb, full panel coverage with Cldn10 (the medullary marker) in the DE.',
}

# Reference genes are visualization aids only. They are never used to cluster, score,
# or assign labels. A reviewer interprets their expression alongside cluster DE.
SEGMENT_REFERENCE_PANEL = {
    'Podocyte': ['Nphs1', 'Nphs2', 'Podxl'],
    'PT-S1': ['Lrp2', 'Cubn', 'Slc34a1', 'Slc5a2', 'Slc5a12', 'Gatm'],
    'PT-S2': ['Slc22a6', 'Slc13a3', 'Cyp2e1'],
    'PT-S3': ['Slc22a7', 'Slc7a13', 'Cyp7b1', 'Slc6a18', 'Acsm3'],
    'DTL': ['Aqp1', 'Slc14a2', 'Corin', 'Fst'],
    'ATL': ['Clcnka', 'Sptssb', 'Akr1b3'],
    'TAL': ['Slc12a1', 'Umod', 'Kcnj1', 'Cldn16'],
    'DCT': ['Slc12a3', 'Pvalb', 'Trpm6', 'Egf'],
    'CNT': ['Calb1', 'Hsd11b2', 'Slc8a1'],
    'CCD': ['Aqp2', 'Aqp3', 'Fxyd4'],
    'OMCD': ['Atp6v0d2', 'Rhcg', 'Foxi1'],
    'IMCD': ['Aqp4', 'Slc14a2', 'Wnt7b'],
    'Vessel': ['Pecam1', 'Cdh5', 'Kdr'],
    'Stroma': ['Col1a1', 'Col3a1', 'Dcn'],
    'SmoothMuscle': ['Acta2', 'Myh11', 'Tagln'],
    'Immune': ['Ptprc', 'Lyz2', 'C1qa'],
}

# Segment vocabulary, families and rollup come from the shared package definition so this
# workflow and the mouse-only one cannot drift apart. See pseudospace/vocabulary.py.
from pseudospace import vocabulary as segment_vocabulary

COARSE_ORDER = list(segment_vocabulary.COARSE_FAMILIES)
KEEP_TUBULE_CLASSES = list(segment_vocabulary.KEEP_TUBULE_CLASSES)
REMOVE_CLASSES = list(segment_vocabulary.REMOVE_CLASSES)
LABEL_VOCABULARY = list(segment_vocabulary.LABEL_VOCABULARY)
# Defined here rather than down in the global-DPT cell where it is also used: the cluster
# spatial review above already needs it, and a top-to-bottom run failed on the ordering.
SEGMENT_DISPLAY_ORDER = list(segment_vocabulary.SEGMENT_DISPLAY_ORDER)

# Family palette. Fine labels inherit their family colour, so a map holding a mix of fine and
# coarse labels stays legible -- and the family is what the coarse continuum reports anyway.
FAMILY_COLORS = {
    'PT': '#4C9BD3', 'DTL': '#9575CD', 'AL': '#F58518', 'DCT': '#D81B60', 'CNT_CD': '#76B7B2',
    'Glomerulus': '#F2C94C', 'Vessel': '#E15759', 'Stroma': '#B07AA1',
    'SmoothMuscle': '#7E57C2', 'Immune': '#59A14F', 'Unassigned': '#9E9E9E',
}


def label_color(label: str) -> str:
    return FAMILY_COLORS.get(segment_vocabulary.coarse_for(str(label)), '#BDBDBD')
TRAJECTORY_COMPARTMENT = 'PT'
N_NEIGHBORS = 30

print('Manual review map:', REVIEWED_CLUSTER_LABELS)
print('All clusters remain available for exploration; PT trajectory requires an explicit current-run PT label.')


# %%
# Purpose: load the already expression-filtered pass-1 object for cluster review.
# The low-support gene filter was deliberately applied before Harmony/HVG selection.
adata_all = sc.read_h5ad(HARMONY_OUTPUT_PATH)
counts = adata_all.layers['counts'].copy() if 'counts' in adata_all.layers else adata_all.X.copy()
adata_all.layers['counts'] = counts.copy()
required_gene_qc = {'n_structures_expressed', 'total_counts', 'passes_expression_count_filter'}
missing_gene_qc = required_gene_qc - set(adata_all.var.columns)
if missing_gene_qc:
    raise KeyError('Pass-1 object lacks pre-Harmony gene-QC annotations: ' + ', '.join(sorted(missing_gene_qc)))
if not adata_all.var['passes_expression_count_filter'].astype(bool).all():
    raise ValueError('Pass-1 object contains genes that failed the pre-Harmony expression filter.')
if 'lognorm' not in adata_all.layers:
    adata_all.X = adata_all.layers['counts'].copy()
    sc.pp.normalize_total(adata_all, target_sum=NORMALIZE_TARGET_SUM)
    sc.pp.log1p(adata_all)
    adata_all.layers['lognorm'] = adata_all.X.copy()
else:
    adata_all.X = adata_all.layers['lognorm'].copy()
adata_all_expression_filtered = adata_all.copy()
print(f'Pre-Harmony expression-filtered matrix: {adata_all.n_obs:,} x {adata_all.n_vars:,}')
display(tubule_filter_audit)



# %%
# Purpose: create a clustering-ready object and discover cluster-specific genes.
# Use highly variable genes only for the reporting object. The Harmony coordinates,
# calculated earlier from shared HVGs, remain the representation used for neighbors.
sc.pp.highly_variable_genes(
    adata_all, n_top_genes=min(N_HVGS, adata_all.n_vars), flavor='seurat'
)
feature_mask = adata_all.var['highly_variable'].to_numpy(dtype=bool)
adata_all.var['selected_for_clustering'] = feature_mask
adata_cluster = adata_all[:, feature_mask].copy()

# Cluster in the integrated representation without supplying any cell-type marker panel.
sc.pp.neighbors(
    adata_cluster, n_neighbors=N_NEIGHBORS, use_rep='X_harmony',
    random_state=RANDOM_STATE,
)
sc.tl.leiden(
    adata_cluster, resolution=COARSE_RESOLUTION, key_added='leiden_coarse',
    flavor='igraph', n_iterations=2, directed=False, random_state=RANDOM_STATE,
)
sc.tl.umap(
    adata_cluster, min_dist=HARMONY_UMAP_MIN_DIST, spread=HARMONY_UMAP_SPREAD,
    random_state=RANDOM_STATE,
)

# Transfer cluster membership to the full retained-gene matrix before testing.
# Clustering uses HVGs; interpretation tests every retained gene, not only HVGs.
adata_all.obs['leiden_coarse'] = adata_cluster.obs['leiden_coarse'].reindex(adata_all.obs_names)
# Differential expression is the evidence used to review cluster identity.
sc.tl.rank_genes_groups(
    adata_all, 'leiden_coarse', method='wilcoxon', pts=True,
    key_added='rank_coarse',
)
import hashlib

CLUSTERING_FINGERPRINT = {
    'n_structures': int(adata_cluster.n_obs),
    'n_reporting_genes': int(adata_cluster.n_vars),
    'resolution': COARSE_RESOLUTION,
    'n_neighbors': N_NEIGHBORS,
    'random_state': RANDOM_STATE,
    'n_clusters': int(adata_cluster.obs['leiden_coarse'].nunique()),
    # Recorded for provenance only -- nothing compares against it. It ties the reviewed labels
    # to the exact membership they were written against.
    'membership_sha1': hashlib.sha1(
        ','.join(adata_cluster.obs['leiden_coarse'].astype(str)).encode()).hexdigest()[:12],
}
display(pd.Series(CLUSTERING_FINGERPRINT, name='value'))


# %% [markdown]
# The 2,000 HVGs plus markers select genes for marker reporting only. Neighbors use the existing Harmony representation built from intersected HVGs; changing this reporting list does not rebuild that representation.
#

# %%
# Purpose: visually review cluster gene programs before assigning any labels.
# Cluster membership was calculated without this reference list. The figures below
# combine data-driven DE with reference-gene context; neither figure assigns a label.

# 1. Label each cluster directly on the UMAP, then show species and sample composition.
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
for axis, color, title in zip(
    axes,
    ['leiden_coarse', 'comparison_species', 'sample'],
    ['Cluster IDs', 'Species', 'Sample'],
):
    sc.pl.embedding(
        adata_cluster, basis='umap', color=color, ax=axis, frameon=False,
        show=False, legend_loc='on data' if color == 'leiden_coarse' else 'right margin',
        title=title,
    )
fig.tight_layout()
fig.savefig(CELLTYPING_DIR / 'coarse_leiden_umap.png', dpi=180, bbox_inches='tight')
plt.show()

# 2. Show the highest-ranked genes discovered from all retained genes for every cluster.
top_de = sc.get.rank_genes_groups_df(adata_all, group=None, key='rank_coarse')
top_markers = (
    top_de.sort_values(['group', 'scores'], ascending=[True, False])
    .groupby('group', observed=True).head(20)
)
top_markers.to_csv(CELLTYPING_DIR / 'coarse_cluster_top_markers_full.csv', index=False)
top_marker_summary = top_markers.groupby('group', observed=True)['names'].apply(
    lambda values: '; '.join(values.head(12))
)
top_marker_summary.to_csv(CELLTYPING_DIR / 'coarse_cluster_top_markers.csv')
display(top_marker_summary.to_frame('top_12_differential_genes'))

sc.pl.rank_genes_groups_heatmap(
    adata_all, key='rank_coarse', n_genes=8, groupby='leiden_coarse',
    standard_scale='var', show=False,
)
plt.savefig(CELLTYPING_DIR / 'coarse_cluster_de_heatmap.png', dpi=180, bbox_inches='tight')
plt.show()

# 3. Overlay reference genes as a visual aid. They are not a scoring or labeling rule.
present_reference_panel = {
    segment: [gene for gene in genes if gene in adata_all.var_names]
    for segment, genes in SEGMENT_REFERENCE_PANEL.items()
}
present_reference_panel = {
    segment: genes for segment, genes in present_reference_panel.items() if genes
}
sc.pl.dotplot(
    adata_all, present_reference_panel, groupby='leiden_coarse',
    standard_scale='var', show=False,
)
plt.savefig(CELLTYPING_DIR / 'cluster_reference_gene_dotplot.png', dpi=180, bbox_inches='tight')
plt.show()

# 4. Convert DE/reference overlap into a compact, non-binding review cue.
# A cue is only the count of reference genes appearing among the top 20 DE genes.
top_gene_sets = top_markers.groupby('group', observed=True)['names'].agg(
    lambda genes: {str(gene).upper() for gene in genes}
)
evidence_rows = []
for cluster, genes in top_gene_sets.items():
    for segment, reference_genes in SEGMENT_REFERENCE_PANEL.items():
        hits = sorted(genes & {gene.upper() for gene in reference_genes})
        evidence_rows.append({
            'cluster': str(cluster), 'reference_segment': segment,
            'n_top_DE_reference_hits': len(hits),
            'top_DE_reference_hits': '; '.join(hits),
        })
reference_evidence = pd.DataFrame(evidence_rows)
reference_evidence.to_csv(
    CELLTYPING_DIR / 'cluster_reference_gene_evidence.csv', index=False
)
evidence_matrix = reference_evidence.pivot(
    index='cluster', columns='reference_segment', values='n_top_DE_reference_hits'
).fillna(0).astype(int)
fig, axis = plt.subplots(figsize=(14, max(3, 0.45 * len(evidence_matrix))))
image = axis.imshow(evidence_matrix.to_numpy(), cmap='Blues', aspect='auto')
axis.set_xticks(range(len(evidence_matrix.columns)), evidence_matrix.columns, rotation=45, ha='right')
axis.set_yticks(range(len(evidence_matrix.index)), evidence_matrix.index)
axis.set_xlabel('Reference program (review aid only)')
axis.set_ylabel('Leiden cluster')
axis.set_title('Reference genes recovered among each cluster’s top 20 DE genes')
for row in range(evidence_matrix.shape[0]):
    for col in range(evidence_matrix.shape[1]):
        value = evidence_matrix.iat[row, col]
        if value:
            axis.text(col, row, str(value), ha='center', va='center', fontsize=8)
fig.colorbar(image, ax=axis, label='Number of top-DE reference genes')
fig.tight_layout()
fig.savefig(CELLTYPING_DIR / 'cluster_reference_evidence_heatmap.png', dpi=180, bbox_inches='tight')
plt.show()

# 5. Show whether a cluster is dominated by one specimen or shared across samples.
cluster_by_sample = pd.crosstab(
    adata_cluster.obs['leiden_coarse'], adata_cluster.obs['sample'], normalize='index'
)
cluster_by_sample.to_csv(CELLTYPING_DIR / 'cluster_sample_fractions.csv')
cluster_by_sample.plot(kind='bar', stacked=True, figsize=(10, 4), colormap='tab20')
plt.ylabel('Fraction of structures in cluster')
plt.xlabel('Leiden cluster')
plt.title('Sample composition within each cluster')
plt.legend(title='Sample', bbox_to_anchor=(1.02, 1), loc='upper left')
plt.tight_layout()
plt.savefig(CELLTYPING_DIR / 'cluster_sample_composition.png', dpi=180, bbox_inches='tight')
plt.show()

cluster_composition = pd.crosstab(
    adata_cluster.obs['leiden_coarse'],
    [adata_cluster.obs['comparison_species'], adata_cluster.obs['sample']],
)
cluster_composition.to_csv(CELLTYPING_DIR / 'cluster_species_sample_counts.csv')
display(cluster_composition)


# %% [markdown]
# ## Review cluster differential expression and assign broad labels
#
# Cluster identity is determined by reviewing cluster-specific differential-expression results,
# the integrated UMAP, and sample composition. No predefined marker panel assigns labels.
# The manual review map is intentionally explicit and must be revisited after reclustering.
#

# %%
# Purpose: apply reviewed broad labels after inspecting differential-expression output.
# No label is inferred here. Update REVIEWED_CLUSTER_LABELS only after reviewing
# coarse_cluster_top_markers.csv, the UMAP, and the cluster composition table above.
cluster_ids = sorted(adata_cluster.obs['leiden_coarse'].astype(str).unique())
unsupported_labels = segment_vocabulary.invalid_labels(REVIEWED_CLUSTER_LABELS.values())
if unsupported_labels:
    raise ValueError(
        f'Unsupported manually reviewed labels: {unsupported_labels}. '
        f'Valid labels are {list(segment_vocabulary.LABEL_VOCABULARY)}.'
    )

# A previous-resolution map is allowed but never applied to IDs absent from this run.
# This is a warning, not an error, so cluster exploration is always available.
obsolete_review = sorted(set(REVIEWED_CLUSTER_LABELS) - set(cluster_ids))
if obsolete_review:
    print(
        'NOTE: review labels for absent cluster IDs were ignored: '
        + ', '.join(obsolete_review)
    )
current_labels = {
    cluster: REVIEWED_CLUSTER_LABELS.get(cluster, 'Unassigned')
    for cluster in cluster_ids
}

cluster_summary = pd.DataFrame({
    'cluster': cluster_ids,
    'broad_label': [current_labels[cluster] for cluster in cluster_ids],
})
cluster_sizes = adata_cluster.obs['leiden_coarse'].astype(str).value_counts()
cluster_summary['n_structures'] = cluster_summary['cluster'].map(cluster_sizes).astype(int)
cluster_summary['review_status'] = np.where(
    cluster_summary['cluster'].isin(REVIEWED_CLUSTER_LABELS),
    'manually reviewed from cluster differential expression',
    'unreviewed; available for exploration',
)
cluster_summary.to_csv(CELLTYPING_DIR / 'coarse_cluster_summary.csv', index=False)
display(cluster_summary)

for obj in (adata_cluster, adata_all):
    if obj is adata_all:
        obj.obs['leiden_coarse'] = adata_cluster.obs['leiden_coarse'].reindex(obj.obs_names)
    # Current clusters without an explicit review stay unassigned but remain available
    # in the saved pass-1 object and all cluster-level exploration outputs.
    obj.obs['segment_class'] = obj.obs['leiden_coarse'].astype(str).map(current_labels)
    # Roll the fine label up to its family rather than copying it, so a cluster labelled PT-S1
    # still reports PT on the coarse continuum.
    obj.obs['coarse_class'] = obj.obs['segment_class'].astype(str).map(segment_vocabulary.coarse_for)
    obj.obs['broad_tubule_marker_call'] = obj.obs['coarse_class']
    obj.obs['review_status'] = np.where(
        obj.obs['leiden_coarse'].astype(str).isin(REVIEWED_CLUSTER_LABELS),
        'manually reviewed from cluster differential expression',
        'unreviewed; excluded from PT trajectory until explicitly selected',
    )
    # Pass 2 retains the tubular nephron; PT is selected only after global DPT review.
    obj.obs['nephron_eligible'] = obj.obs['coarse_class'].isin(KEEP_TUBULE_CLASSES)
    obj.obs['trajectory_eligible'] = obj.obs['coarse_class'].eq(TRAJECTORY_COMPARTMENT)
    obj.obs['trajectory_exclusion_reason'] = np.where(
        obj.obs['nephron_eligible'], '',
        np.where(obj.obs['coarse_class'].eq('Glomerulus'),
                 'renal corpuscle; excluded from the tubular-nephron continuum',
                 np.where(obj.obs['coarse_class'].eq('Unassigned'),
                          'ambiguous or mixed cluster; excluded pending review',
                          'non-nephron structure')),
    )

for column in adata_all.obs:
    adata_combined.obs[column] = adata_all.obs[column].reindex(adata_combined.obs_names)
adata_combined.write(HARMONY_OUTPUT_PATH)
adata_all.obs.to_csv(CELLTYPING_DIR / 'structure_annotations_and_eligibility.csv')

# Visualize the reviewed calls directly over their original Leiden clusters.
fig, axis = plt.subplots(figsize=(9, 7))
sc.pl.embedding(
    adata_cluster, basis='umap', color='segment_class', ax=axis,
    frameon=False, show=False, title='Reviewed broad labels on Leiden clusters',
)
coordinates = np.asarray(adata_cluster.obsm['X_umap'])
cluster_values = adata_cluster.obs['leiden_coarse'].astype(str).to_numpy()
for cluster in cluster_ids:
    mask = cluster_values == cluster
    if mask.any():
        x_pos, y_pos = np.median(coordinates[mask], axis=0)
        label = f'{cluster}: {current_labels[cluster]}'
        text = axis.text(
            x_pos, y_pos, label, ha='center', va='center', fontsize=9,
            fontweight='bold', color='black',
        )
        text.set_path_effects([
            PathEffects.withStroke(linewidth=3, foreground='white')
        ])
fig.tight_layout()
fig.savefig(
    CELLTYPING_DIR / 'reviewed_cluster_labels_on_umap.png',
    dpi=180, bbox_inches='tight',
)
plt.show()

# Write and display the evidence used for each manual call. This report does not
# recompute labels; it makes the current review decision auditable.
reference_hits = (
    reference_evidence[reference_evidence['n_top_DE_reference_hits'].gt(0)]
    .groupby('cluster', observed=True)
    .apply(
        lambda rows: '; '.join(
            f"{row.reference_segment}: {row.top_DE_reference_hits}"
            for row in rows.itertuples()
        )
    )
)
assignment_evidence = cluster_summary.copy()
assignment_evidence['interpretation'] = assignment_evidence['cluster'].map(
    ASSIGNMENT_INTERPRETATION
)
assignment_evidence['top_DE_genes'] = assignment_evidence['cluster'].map(
    top_marker_summary
)
assignment_evidence['reference_gene_evidence'] = assignment_evidence['cluster'].map(
    reference_hits
).fillna('No reference genes among the top 20 DE genes')
assignment_evidence.to_csv(
    CELLTYPING_DIR / 'reviewed_cluster_assignment_evidence.csv', index=False
)
display(assignment_evidence[[
    'cluster', 'broad_label', 'interpretation', 'top_DE_genes',
    'reference_gene_evidence',
]])

# Export the same evidence as a compact figure for handoff and manuscript records.
import textwrap
figure_rows = []
for row in assignment_evidence.itertuples():
    figure_rows.append([
        row.cluster,
        row.broad_label,
        textwrap.fill(row.interpretation, width=34),
        textwrap.fill(row.top_DE_genes, width=46),
        textwrap.fill(row.reference_gene_evidence, width=42),
    ])
fig, axis = plt.subplots(figsize=(22, max(6, 0.85 * len(figure_rows) + 1)))
axis.axis('off')
table = axis.table(
    cellText=figure_rows,
    colLabels=[
        'Leiden cluster', 'Reviewed label', 'Interpretation',
        'Top differential genes', 'Reference-gene evidence',
    ],
    cellLoc='left', colLoc='left', loc='center',
    colWidths=[0.07, 0.10, 0.20, 0.31, 0.32],
)
table.auto_set_font_size(False)
table.set_fontsize(7.5)
table.scale(1, 3.1)
for (row_index, column_index), cell in table.get_celld().items():
    if row_index == 0:
        cell.set_text_props(weight='bold')
        cell.set_facecolor('#D9EAF7')
fig.suptitle(
    'Evidence supporting reviewed broad cluster labels\nReference matches are review aids, not automated assignments',
    y=0.98, fontsize=13,
)
fig.tight_layout(rect=(0, 0, 1, 0.95))
fig.savefig(
    CELLTYPING_DIR / 'reviewed_cluster_assignment_evidence.png',
    dpi=180, bbox_inches='tight',
)
plt.show()


# %% [markdown]
# ## Spatial check: reviewed clusters in kidney coordinates
#
# Each point is a segmented structure at its centroid. Colors indicate the reviewed broad
# label (top row) or the original Leiden cluster (bottom row). These are coordinate-only
# sanity checks: they do not display tissue images or segmentation boundaries.
#
#

# %%
# Purpose: check whether reviewed labels occupy plausible spatial territories.
# This is a sanity check only; expression-based labels remain the primary evidence.
required_spatial_columns = {'x_centroid', 'y_centroid', 'sample', 'segment_class', 'leiden_coarse'}
missing_spatial_columns = required_spatial_columns - set(adata_all.obs.columns)
if missing_spatial_columns:
    print('Spatial cluster plot skipped; missing:', sorted(missing_spatial_columns))
else:
    spatial = adata_all.obs.copy()
    spatial['x_centroid'] = pd.to_numeric(spatial['x_centroid'], errors='coerce')
    spatial['y_centroid'] = pd.to_numeric(spatial['y_centroid'], errors='coerce')
    spatial = spatial.loc[
        np.isfinite(spatial['x_centroid']) & np.isfinite(spatial['y_centroid'])
    ].copy()

    observed_labels = list(dict.fromkeys(spatial['segment_class'].astype(str)))
    label_order = [
        label for label in SEGMENT_DISPLAY_ORDER if label in observed_labels
    ] + [label for label in observed_labels if label not in SEGMENT_DISPLAY_ORDER]
    label_colors = {label: label_color(label) for label in label_order}

    cluster_order = sorted(spatial['leiden_coarse'].astype(str).unique(), key=int)
    cluster_cmap = plt.get_cmap('tab20', max(len(cluster_order), 1))
    cluster_colors = {
        cluster: cluster_cmap(index)
        for index, cluster in enumerate(cluster_order)
    }

    fig, axes = plt.subplots(2, len(SAMPLE_ORDER), figsize=(4.2 * len(SAMPLE_ORDER), 8))
    for column, sample_name in enumerate(SAMPLE_ORDER):
        sample_data = spatial.loc[spatial['sample'].astype(str).eq(sample_name)]
        for label in label_order:
            subset = sample_data.loc[sample_data['segment_class'].astype(str).eq(label)]
            axes[0, column].scatter(
                subset['x_centroid'], subset['y_centroid'], s=2.5,
                color=label_colors[label], alpha=0.72, linewidths=0,
                label=label if column == 0 else None,
            )
        for cluster in cluster_order:
            subset = sample_data.loc[sample_data['leiden_coarse'].astype(str).eq(cluster)]
            axes[1, column].scatter(
                subset['x_centroid'], subset['y_centroid'], s=2.5,
                color=cluster_colors[cluster], alpha=0.72, linewidths=0,
                label=f'Cluster {cluster}' if column == 0 else None,
            )
        for row in range(2):
            axes[row, column].set_aspect('equal')
            axes[row, column].invert_yaxis()
            axes[row, column].set_xticks([])
            axes[row, column].set_yticks([])
        axes[0, column].set_title(f'{sample_name}\nreviewed labels')
        axes[1, column].set_title(f'{sample_name}\nLeiden clusters')

    label_handles, label_names = axes[0, 0].get_legend_handles_labels()
    cluster_handles, cluster_names = axes[1, 0].get_legend_handles_labels()
    fig.legend(
        label_handles, label_names, title='Reviewed broad label',
        loc='center left', bbox_to_anchor=(1.0, 0.72), frameon=False,
    )
    fig.legend(
        cluster_handles, cluster_names, title='Leiden cluster',
        loc='center left', bbox_to_anchor=(1.0, 0.24), frameon=False,
    )
    fig.suptitle(
        'Spatial distribution of reviewed labels and original clusters\nEach point is a structure centroid; y-axis is inverted to match image coordinates',
        y=0.98, fontsize=13,
    )
    fig.tight_layout(rect=(0, 0, 0.88, 0.93))
    fig.savefig(
        CELLTYPING_DIR / 'reviewed_labels_and_clusters_in_kidney_coordinates.png',
        dpi=220, bbox_inches='tight',
    )
    plt.show()


# %% [markdown]
# ## Segmentation QC: representative structures from each cluster
#
# Each crop is drawn from the original segmentation GeoJSON. The selected structure is colored by
# its reviewed broad label; nearby segmentation polygons are gray. This checks whether the cluster
# corresponds to plausible structure morphology and whether segmentation boundaries look sensible.
#
# The data package does not contain matched mouse tissue images, so these are geometry-based
# segmentation crops rather than histology overlays.
#
#

# %%
# Purpose: inspect several original segmentation polygons from every Leiden cluster.
# Representative structures are selected deterministically across samples so each rerun is comparable.
SEGMENTATION_GEOJSON_BY_SAMPLE = {
    'Ctrl1A2': DATA_ROOT / 'Ctrl_1A2_v4.geojson',
    'Ctrl1A4': DATA_ROOT / 'Ctrl_1A4_v4.geojson',
    'HUK1_COR1': DATA_ROOT / 'HUK1_COR1_v2.geojson',
    'HUK1_MED1': DATA_ROOT / 'HUK1_MED1_v2.geojson',
}
N_SEGMENTATION_EXAMPLES_PER_CLUSTER = 3
# Full-resolution TIFFs use the same pixel coordinate system as the segmentation GeoJSON.
# They are uncompressed and several gigabytes each, so crops are read with np.memmap rather
# than loading an entire slide into memory.
TIFF_BACKGROUND_BY_SAMPLE = {
    'Ctrl1A2': DATA_ROOT / 'histology_images' / 'Ctrl_1A2.tif',
    'Ctrl1A4': DATA_ROOT / 'histology_images' / 'Ctrl_1A4.tif',
    'HUK1_COR1': DATA_ROOT / 'histology_images' / 'HUK1_COR_1.tif',
    'HUK1_MED1': DATA_ROOT / 'histology_images' / 'HUK1_MED_1.tif',
}
def _outer_rings(feature):
    # Return the exterior ring(s) without requiring a geometry library.
    geometry = feature.get('geometry') or {}
    if geometry.get('type') == 'Polygon':
        return [np.asarray(geometry['coordinates'][0], dtype=float)]
    if geometry.get('type') == 'MultiPolygon':
        return [
            np.asarray(polygon[0], dtype=float)
            for polygon in geometry['coordinates']
        ]
    return []

def _draw_feature(axis, feature, *, color, linewidth=0.7, alpha=1.0, fill=False):
    for ring in _outer_rings(feature):
        if len(ring) < 3:
            continue
        axis.plot(ring[:, 0], ring[:, 1], color=color, lw=linewidth, alpha=alpha)
        if fill:
            axis.fill(ring[:, 0], ring[:, 1], color=color, alpha=0.16)

missing_geojson = {
    sample: str(path) for sample, path in SEGMENTATION_GEOJSON_BY_SAMPLE.items()
    if not path.exists()
}
required_crop_columns = {'sample', 'feature_index', 'x_centroid', 'y_centroid', 'leiden_coarse', 'segment_class'}
missing_crop_columns = required_crop_columns - set(adata_all.obs.columns)

if missing_geojson or missing_crop_columns:
    print('Segmentation crop panel skipped.')
    if missing_geojson:
        print('Missing GeoJSON:', missing_geojson)
    if missing_crop_columns:
        print('Missing observation columns:', sorted(missing_crop_columns))
else:
    # Load each source segmentation once, then index its polygons by feature_index.
    segmentation_features = {
        sample: json.loads(path.read_text()).get('features', [])
        for sample, path in SEGMENTATION_GEOJSON_BY_SAMPLE.items()
    }
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None
    histology_memmaps = {}
    for sample, image_path in TIFF_BACKGROUND_BY_SAMPLE.items():
        if not image_path.exists():
            continue
        image = Image.open(image_path)
        if image.mode != 'RGB' or not image.tile or image.tile[0].codec_name != 'raw':
            raise ValueError(
                f'{image_path.name} must be a single-plane uncompressed RGB TIFF for crop QC.'
            )
        width, height = image.size
        offset = int(image.tile[0].offset)
        histology_memmaps[sample] = np.memmap(
            image_path, dtype=np.uint8, mode='r', offset=offset,
            shape=(height, width, 3),
        )
        image.close()
    crop_obs = adata_all.obs.copy()
    crop_obs['feature_index'] = pd.to_numeric(crop_obs['feature_index'], errors='coerce')
    crop_obs['x_centroid'] = pd.to_numeric(crop_obs['x_centroid'], errors='coerce')
    crop_obs['y_centroid'] = pd.to_numeric(crop_obs['y_centroid'], errors='coerce')
    crop_obs = crop_obs.dropna(subset=['feature_index', 'x_centroid', 'y_centroid']).copy()
    crop_obs['feature_index'] = crop_obs['feature_index'].astype(int)
    crop_obs['leiden_coarse'] = crop_obs['leiden_coarse'].astype(str)

    cluster_order = sorted(crop_obs['leiden_coarse'].unique(), key=int)
    label_colors = {
        label: label_color(label)
        for label in dict.fromkeys(crop_obs['segment_class'].astype(str))
    }
    selected_rows = []
    for cluster in cluster_order:
        candidates = crop_obs.loc[crop_obs['leiden_coarse'].eq(cluster)].sort_values(
            ['sample', 'feature_index']
        )
        # Spread examples across the ordered structures rather than selecting one local patch.
        positions = np.linspace(
            0, len(candidates) - 1,
            min(N_SEGMENTATION_EXAMPLES_PER_CLUSTER, len(candidates)),
            dtype=int,
        )
        selected_rows.extend(candidates.iloc[np.unique(positions)].to_dict('records'))

    fig, axes = plt.subplots(
        len(cluster_order), N_SEGMENTATION_EXAMPLES_PER_CLUSTER,
        figsize=(4.2 * N_SEGMENTATION_EXAMPLES_PER_CLUSTER, 3.8 * len(cluster_order)),
        squeeze=False,
    )
    for axis in axes.ravel():
        axis.set_axis_off()

    for row_index, cluster in enumerate(cluster_order):
        examples = [row for row in selected_rows if row['leiden_coarse'] == cluster]
        for column_index, row in enumerate(examples):
            axis = axes[row_index, column_index]
            sample_name = str(row['sample'])
            feature_index = int(row['feature_index'])
            features = segmentation_features[sample_name]
            if not 0 <= feature_index < len(features):
                axis.text(0.5, 0.5, 'feature index absent', ha='center', va='center')
                continue

            target = features[feature_index]
            rings = _outer_rings(target)
            if not rings:
                axis.text(0.5, 0.5, 'polygon absent', ha='center', va='center')
                continue
            points = np.concatenate(rings, axis=0)
            x_min, y_min = points.min(axis=0)
            x_max, y_max = points.max(axis=0)
            pad = max(150.0, 1.7 * max(x_max - x_min, y_max - y_min))

            # Read only the requested full-resolution TIFF crop beneath the polygons.
            # The GeoJSON and TIFF share pixel coordinates, so no scale factor is applied.
            if sample_name in histology_memmaps:
                slide = histology_memmaps[sample_name]
                crop_left = max(0, int(np.floor(x_min - pad)))
                crop_right = min(slide.shape[1], int(np.ceil(x_max + pad)))
                crop_top = max(0, int(np.floor(y_min - pad)))
                crop_bottom = min(slide.shape[0], int(np.ceil(y_max + pad)))
                axis.imshow(
                    slide[crop_top:crop_bottom, crop_left:crop_right],
                    extent=[crop_left, crop_right, crop_bottom, crop_top],
                    origin='upper', alpha=0.94, zorder=0,
                )

            # Find nearby QC-kept structures by centroid and draw their source polygons.
            sample_obs = crop_obs.loc[crop_obs['sample'].astype(str).eq(sample_name)]
            nearby = sample_obs.loc[
                sample_obs['x_centroid'].between(x_min - pad, x_max + pad)
                & sample_obs['y_centroid'].between(y_min - pad, y_max + pad)
            ]
            for neighbor_index in nearby['feature_index'].astype(int):
                if neighbor_index == feature_index or not 0 <= neighbor_index < len(features):
                    continue
                _draw_feature(axis, features[neighbor_index], color='#BDBDBD', linewidth=0.45, alpha=0.7)

            label = str(row['segment_class'])
            _draw_feature(
                axis, target, color=label_colors.get(label, '#333333'),
                linewidth=1.8, alpha=1.0, fill=True,
            )
            axis.set_xlim(x_min - pad, x_max + pad)
            axis.set_ylim(y_max + pad, y_min - pad)
            axis.set_aspect('equal')
            background_note = (
                'full-resolution TIFF background'
                if sample_name in histology_memmaps
                else 'segmentation only: TIFF unavailable'
            )
            axis.set_title(
                f'Cluster {cluster} | {label}\n{sample_name}, feature {feature_index} | {background_note}',
                fontsize=8,
            )
    fig.suptitle(
        'Representative segmentation crops by Leiden cluster\nColor: selected structure; gray: nearby segmentations; crops use full-resolution TIFF background',
        y=0.995, fontsize=14,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    fig.savefig(
        CELLTYPING_DIR / 'representative_segmentations_by_cluster.png',
        dpi=220, bbox_inches='tight',
    )
    plt.show()


# %% [markdown]
# ## Nephron reintegration and global DPT review
#
# The pass-1 object remains complete. Pass 2 retains only manually reviewed tubular-nephron classes (`PT`, loop, TAL, DCT, and CNT/CD). The glomerular cluster remains available in pass 1 and spatial QC, but it is excluded here because the renal corpuscle is not part of the tubular continuum; unassigned/mixed clusters are also excluded pending review. Harmony and DPT are then rebuilt across the retained nephron cohort. PT is selected only after this global inspection.
#

# %%
# Purpose: retain reviewed tubular-nephron classes for pass-2 reintegration.
# This post-review filter removes structures that cannot lie on the tubular continuum;
# it does not pre-select PT before integration.
n_before_nephron_filter = int(adata_all.n_obs)
nephron_mask = adata_all.obs['nephron_eligible'].to_numpy(dtype=bool)
adata_tubule = adata_all[nephron_mask].copy()

# `fraction_removed` is the count removed at this step divided by all pass-1 QC-passing structures.
nephron_filter_counts = (adata_all.obs.assign(pass2_decision=np.where(nephron_mask, 'retained', 'removed'))
    .groupby(['coarse_class', 'pass2_decision'], observed=True).size()
    .rename('n_structures').reset_index())
nephron_filter_counts['fraction_of_pass1'] = nephron_filter_counts['n_structures'] / n_before_nephron_filter
nephron_filter_counts.to_csv(CELLTYPING_DIR / 'pass2_nephron_filter_counts.csv', index=False)
display(nephron_filter_counts.sort_values(['pass2_decision', 'coarse_class']))

tubule_filter_audit.loc[len(tubule_filter_audit)] = {
    'filter_stage': 'reviewed non-nephron and unresolved-cluster removal before pass 2',
    'filter_type': 'nephron filter', 'n_input_tubules': n_before_nephron_filter,
    'n_retained_tubules': int(adata_tubule.n_obs),
    'n_removed_tubules': n_before_nephron_filter - int(adata_tubule.n_obs),
    'fraction_removed': 1 - int(adata_tubule.n_obs) / n_before_nephron_filter,
}
tubule_filter_audit.to_csv(DIAGNOSTIC_DIR / 'tubule_filter_audit.csv', index=False)
if adata_tubule.n_obs < 100 or adata_tubule.obs['comparison_species'].nunique() < 2:
    raise ValueError('Too few reviewed nephron structures from both species for pass-2 Harmony.')
for column in ('coarse_class', 'segment_class'):
    adata_tubule.obs[column] = pd.Categorical(adata_tubule.obs[column])
print(f'Pass-2 nephron cohort: {adata_tubule.n_obs:,} retained; {n_before_nephron_filter - adata_tubule.n_obs:,} removed.')
display(tubule_filter_audit)


# %%
# Purpose: # Pass-2 Harmony on the reviewed tubular-nephron cohort; labels are carried, never recomputed.
# Pass-2 Harmony on the reviewed tubular-nephron cohort; labels are carried, never recomputed.
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
    batch_key=BATCH_KEY,  # 'sample': aligns across 4 biological samples (Ctrl1A2, Ctrl1A4, HUK1_COR1, HUK1_MED1)
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


# %%
# Purpose: # Pass-2 nephron mixing is descriptive; it cannot validate segment conservation.
# Pass-2 nephron mixing is descriptive; it cannot validate segment conservation.
def safe_asw(embedding, labels):
    return float(silhouette_score(embedding, labels)) if 1 < len(np.unique(labels)) < len(labels) else np.nan
rng = np.random.default_rng(RANDOM_STATE)
idx = rng.choice(adata_tubule.n_obs, min(2000, adata_tubule.n_obs), replace=False)
integration_qc = pd.DataFrame([dict(stage='pass2_nephron', embedding=rep, sample_asw=safe_asw(adata_tubule.obsm[rep][idx], adata_tubule.obs['sample'].to_numpy()[idx]), species_asw=safe_asw(adata_tubule.obsm[rep][idx], adata_tubule.obs['comparison_species'].to_numpy()[idx])) for rep in ('X_pca', 'X_harmony')])
integration_qc.to_csv(CELLTYPING_DIR / 'integration_qc_silhouette.csv', index=False)
display(integration_qc)


# %% [markdown]
# ## Reference only: PT positional genes for global-DPT orientation
#
# These genes do not cluster or label structures. They only orient the global nephron graph: the root should lie at the low end of this proximal marker axis. Inspect PAGA and the segment distributions before treating a branched nephron graph as one line.
#

# %%
# Purpose: define reference genes used only after manual PT cluster review.
# These lists do not participate in clustering, differential expression, or broad labels.
NEPHRON_AXIS_MARKERS = {
    'PT-S1': ['Slc5a2', 'Slc5a12', 'Gatm'],
    'PT-S2': ['Slc22a6', 'Slc13a3', 'Cyp2e1'],
    'PT-S3': ['Slc7a13', 'Slc22a7', 'Cyp7b1'],
}
print('Reference-only PT groups:', ', '.join(NEPHRON_AXIS_MARKERS))


# %% [markdown]
# ### 2.4 - PT trajectory settings
#
# **Purpose:** declare the PT marker axis and graph diagnostics just before constructing pseudotime.
#

# %%
# Purpose: # Global-DPT orientation settings: defined where the trajectory is constructed.
# Global-DPT orientation settings: defined where the trajectory is constructed.
TOTAL_POSITION_MARKERS = {
    'early': sorted({gene for group in ('PT-S1', 'PT-S2') for gene in NEPHRON_AXIS_MARKERS[group]}),
    'late': sorted(NEPHRON_AXIS_MARKERS['PT-S3']),
}
PAGA_CONNECTIVITY_THRESHOLD = 0.01
N_DIFFMAP_COMPONENTS_TO_TEST, EIGENVALUE_FLOOR, MIN_VALID_FOR_SPEARMAN = 10, 1e-8, 20
print('PT reference axis for global-DPT orientation:', {key: len(value) for key, value in TOTAL_POSITION_MARKERS.items()})


# %%
# Purpose: # Score a PT reference axis on the re-harmonized whole-nephron cohort.
# Legacy total_* column names now refer to the global tubular-nephron coordinate.
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
# Purpose: # Inspect global-nephron graph topology before interpreting a single continuous DPT.
# Inspect global-nephron graph topology before interpreting a single continuous DPT.
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
adata_total.obs['leiden_coarse'] = adata_total.obs['leiden_coarse'].astype('category').cat.remove_unused_categories()
if adata_total.obs['leiden_coarse'].nunique() > 1:
    sc.tl.paga(
        adata_total, groups='leiden_coarse', neighbors_key=trajectory_neighbors
    )
    sc.pl.paga(
        adata_total, threshold=PAGA_CONNECTIVITY_THRESHOLD,
        color='leiden_coarse', frameon=False, show=False,
    )
    plt.savefig(CELLTYPING_DIR / 'paga_segment_topology.png', dpi=180, bbox_inches='tight')
    plt.show()
    
    paga_connectivity = pd.DataFrame(
        np.asarray(adata_total.uns['paga']['connectivities'].todense()),
        index=adata_total.obs['leiden_coarse'].cat.categories,
        columns=adata_total.obs['leiden_coarse'].cat.categories,
    )
    paga_connectivity.to_csv(CELLTYPING_DIR / 'paga_connectivity_matrix.csv')
    paga_edges = paga_connectivity.stack().rename('connectivity').reset_index()
    paga_edges.columns = ['from_cluster', 'to_cluster', 'connectivity']
    paga_edges.to_csv(CELLTYPING_DIR / 'global_nephron_paga_cluster_edges.csv', index=False)
    display(paga_edges)
else:
    print('The retained nephron cohort has one cluster; PAGA between-cluster topology is unavailable.')


# %%
# Purpose: sc.tl.diffmap(
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

# DPT can leave graph components disconnected from the selected root; remove only those.
n_before_root_connectivity = int(adata_total.n_obs)
nonfinite = ~np.isfinite(total_scanpy_dpt)
if nonfinite.any():
    print(f'Dropping {int(nonfinite.sum()):,} root-disconnected structures.')
    adata_total = adata_total[~nonfinite].copy()
    marker_axis = marker_axis[~nonfinite]
    total_scanpy_dpt = total_scanpy_dpt[~nonfinite]

tubule_filter_audit.loc[len(tubule_filter_audit)] = {
    'filter_stage': 'connected to selected global-nephron DPT root',
    'filter_type': 'tubule filter',
    'n_input_tubules': n_before_root_connectivity,
    'n_retained_tubules': int(adata_total.n_obs),
    'n_removed_tubules': n_before_root_connectivity - int(adata_total.n_obs),
    'fraction_removed': 1 - int(adata_total.n_obs) / n_before_root_connectivity,
}
tubule_filter_audit.to_csv(DIAGNOSTIC_DIR / 'tubule_filter_audit.csv', index=False)
print('Spearman(DPT, marker axis) =',
      round(spearmanr(total_scanpy_dpt, marker_axis).correlation, 3))
display(tubule_filter_audit)
sc.pl.embedding(
    adata_total,
    basis='umap',
    color=['total_scanpy_dpt', 'segment_class', 'comparison_species'],
    cmap='viridis',
    frameon=False,
    ncols=3,
    show=False,
)
plt.savefig(CELLTYPING_DIR / 'global_nephron_dpt_umap.png', dpi=180, bbox_inches='tight')
plt.show()


# %% [markdown]
# ### Global-DPT summary figures: reviewed segments versus marker programs
#
# Three different quantities share this section, and two of them put the same `PT-S1`/`PT-S2`/`PT-S3`
# labels on their x-axis. Read this table before comparing any two bars.
#
# | Figure | x-axis labels | How it is made | Counts mean |
# | --- | --- | --- | --- |
# | `global_nephron_dpt_umap.png` | — (colour only) | embedding coloured by `total_scanpy_dpt`, `segment_class`, `comparison_species` | — |
# | `global_nephron_dpt_by_segment.png` (next cell) | reviewed `segment_class`: `PT-S1/S2/S3`, `AL`, `mTAL`, `DCT1`, `CNT_CD` | the manual cluster review map `REVIEWED_CLUSTER_LABELS`, keyed by Leiden cluster ID | structures carrying that reviewed label |
# | `global_nephron_dpt_by_fine_reference_program.png` (below) | marker programs: `PT-S1` ... `IMCD` | per-structure argmax of marker-panel scores: mean lognorm per panel over the genes detected in >= 5% of that species' structures, z-scored **within species**, >= 2 usable genes per species, `Unresolved` when the winner is within 0.15 SD of the runner-up | structures whose strongest *marker program* is that program — a diagnostic, not a label |
#
# `Unresolved` is never drawn; it is counted in each panel title and audited in the CSV. The reviewed
# figure is the one that reports segment sizes. The program figure must not be read as segment counts
# and must not be compared bar-for-bar with the reviewed figure: the same label appears on both axes
# with two different meanings.
#
# Two standing caveats: (1) the reviewed map is keyed by Leiden cluster ID and this notebook records
# `CLUSTERING_FINGERPRINT` without pinning the map to it, so re-read `REVIEWED_CLUSTER_LABELS` after
# any change to the clustering inputs; (2) the marker programs are reference panels, not annotations,
# and panels that lost genes to ortholog mapping (for example `ATL`, which keeps only `Clcnka`) cannot
# separate the segments they name.
#

# %%
# Purpose: segment_dpt = (
segment_dpt = (
    adata_total.obs.groupby('segment_class', observed=True)['total_scanpy_dpt']
    .agg(['size', 'mean', 'median'])
    .reindex([
        value for value in SEGMENT_DISPLAY_ORDER
        if value in set(adata_total.obs['segment_class'].astype(str))
    ])
)
segment_dpt.to_csv(CELLTYPING_DIR / 'global_dpt_by_segment.csv')
display(segment_dpt.round(3))

species_segment_dpt = (
    adata_total.obs.groupby(
        ['comparison_species', 'segment_class'], observed=True
    )['total_scanpy_dpt']
    .agg(['size', 'median'])
    .reset_index()
)
species_segment_dpt.to_csv(
    CELLTYPING_DIR / 'global_dpt_by_species_and_segment.csv', index=False
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
    labels = [f'{seg}\n(n={len(d):,})' for seg, d in zip(order, data)]
    ax.set_xticklabels(labels)
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
    ax.set_title(f'{species.title()} (n={int(sp_mask.sum()):,})')
    ax.set_ylabel('total_scanpy_dpt' if ax == axes[0] else '')

fig.suptitle('Global-nephron DPT by reviewed segment (inspect branches; not a forced linear anatomy)', y=1.02)
fig.savefig(CELLTYPING_DIR / 'global_nephron_dpt_by_segment.png', dpi=150, bbox_inches='tight')
plt.show()


# %% [markdown]
# ### Global-DPT visualization: broad labels and fine reference programs
#
# First inspect the global DPT colored by the reviewed broad labels and original Leiden IDs.
# Then inspect expression of fine reference programs along the same coordinate, separately by
# species. These fine markers are visual validation only: they do not generate or overwrite
# any cluster labels.
#

# %%
# Purpose: visualize the global DPT with reviewed broad labels and fine reference programs.
# Broad labels are manual cluster-review calls. Fine reference genes below are plotted only
# to assess whether each part of the global graph expresses the expected nephron program.
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
for axis, color, title in zip(
    axes,
    ['coarse_class', 'leiden_coarse', 'total_scanpy_dpt'],
    ['Reviewed broad nephron labels', 'Original Leiden clusters', 'Global-nephron DPT'],
):
    sc.pl.embedding(
        adata_total, basis='umap', color=color, ax=axis, frameon=False,
        show=False, title=title, cmap='viridis',
        legend_loc='on data' if color == 'leiden_coarse' else 'right margin',
    )
fig.tight_layout()
fig.savefig(CELLTYPING_DIR / 'global_dpt_broad_labels_umap.png', dpi=180, bbox_inches='tight')
plt.show()

# Keep only fine tubular programs represented in the reference list and in this matrix.
# The coverage table makes absent orthologs/genes explicit instead of silently omitting them.
# This explicit order is intentionally limited to tubular-nephron programs.
# Do not add glomerular, vascular, stromal, smooth-muscle, immune, or unresolved programs here.
FINE_NEPHRON_PROGRAM_ORDER = [
    'PT-S1', 'PT-S2', 'PT-S3', 'DTL', 'ATL', 'TAL', 'DCT', 'CNT', 'CCD', 'OMCD', 'IMCD',
]
global_fine_panel = {
    program: [gene for gene in SEGMENT_REFERENCE_PANEL[program] if gene in adata_total.var_names]
    for program in FINE_NEPHRON_PROGRAM_ORDER
    if program in SEGMENT_REFERENCE_PANEL
}
global_fine_panel = {program: genes for program, genes in global_fine_panel.items() if genes}
fine_marker_coverage = pd.DataFrame([
    {'fine_program': program,
     'requested_genes': len(SEGMENT_REFERENCE_PANEL[program]),
     'present_genes': len(genes),
     'genes_plotted': '; '.join(genes)}
    for program, genes in global_fine_panel.items()
])
fine_marker_coverage.to_csv(CELLTYPING_DIR / 'global_dpt_fine_marker_coverage.csv', index=False)
display(fine_marker_coverage)

# Give each fine program a stable color strip. It is a display key only, not an annotation.
global_fine_order = list(global_fine_panel)
global_fine_colors = {
    program: plt.get_cmap('tab20')(index % 20)
    for index, program in enumerate(global_fine_order)
}
for species in SPECIES_GROUPS:
    species_mask = adata_total.obs['comparison_species'].astype(str).eq(species).to_numpy()
    plot_marker_heatmap(
        adata_total,
        {program: [{'label': gene, 'candidates': [gene]} for gene in genes]
         for program, genes in global_fine_panel.items()},
        global_fine_order,
        global_fine_colors,
        pseudotime_col='total_scanpy_dpt',
        cluster_col=None,
        mask=species_mask,
        title=f'{species.title()} fine nephron reference programs along global DPT',
        output_name=f'{species}_global_dpt_fine_marker_heatmap.png',
        strip_col='coarse_class',
        strip_order=COARSE_ORDER,
        strip_colors=FAMILY_COLORS,
        n_bins=120,
        output_dir=CELLTYPING_DIR,
        project_dir=PROJECT_DIR,
    )
print('Saved broad-label UMAP and species-specific fine-marker global-DPT heatmaps.')


# %%
# Purpose: summarize global DPT by the dominant fine reference program.
# HOW THIS FIGURE IS BUILT -- keep this in sync with the markdown above and the notes cell:
#   1. panel score = mean lognorm expression of that program's usable panel genes
#   2. a gene is usable in a species only when it is detected in >= MIN_PANEL_GENE_DETECTED_FRACTION
#      of that species' structures, so a panel that is silent there cannot outrank one that is not
#      (human loses Cyp2e1, Slc14a2, most PT-S3 genes, and the whole IMCD panel)
#   3. a program needs >= MIN_PROGRAM_GENES_FOR_CALL usable genes in that species to be callable
#   4. scores are z-scored WITHIN each species, so the species-wide expression offset cannot decide
#      the call (mouse sits up to 0.7 SD above human on PT-S3/IMCD before standardisation)
#   5. call = strongest program only when its z is positive and beats the runner-up by
#      MARGIN_TO_SECOND SD, otherwise 'Unresolved'
# This is a marker-expression diagnostic, not an automated fine-segment annotation, and it is not a
# reviewed-segment count: the reviewed labels are `segment_class` (figure above).
MIN_PROGRAM_GENES_FOR_CALL = 2
MIN_PANEL_GENE_DETECTED_FRACTION = MIN_GENE_TUBULE_FRACTION
MARGIN_TO_SECOND = 0.15

expression_for_fine_scores = (
    adata_total.layers['lognorm'] if 'lognorm' in adata_total.layers else adata_total.X
)
species_for_fine_scores = adata_total.obs['comparison_species'].astype(str)


def _panel_mean_expression(genes: list[str]) -> np.ndarray:
    """Mean lognorm expression over `genes`, sparse or dense."""
    values = expression_for_fine_scores[:, [adata_total.var_names.get_loc(gene) for gene in genes]]
    return np.asarray(values.mean(axis=1)).ravel()


def _gene_detected_fraction(gene: str) -> pd.Series:
    """Fraction of each species' structures with a non-zero value for `gene`."""
    column = expression_for_fine_scores[:, adata_total.var_names.get_loc(gene)]
    dense = np.asarray(column.todense()).ravel() if hasattr(column, 'todense') else np.asarray(column).ravel()
    return pd.Series(dense > 0, index=adata_total.obs_names).groupby(
        species_for_fine_scores, observed=True
    ).mean()


fine_panel_gene_detected_fraction = {
    gene: _gene_detected_fraction(gene)
    for gene in sorted({gene for genes in global_fine_panel.values() for gene in genes})
}

# Eligibility is per species: the diagnostic must not rank a panel that is silent in that species.
callable_fine_panels = {species: {} for species in SPECIES_GROUPS}
eligibility_rows = []
for program, genes in global_fine_panel.items():
    for species in SPECIES_GROUPS:
        usable = [
            gene for gene in genes
            if float(fine_panel_gene_detected_fraction[gene].get(species, 0.0))
            >= MIN_PANEL_GENE_DETECTED_FRACTION
        ]
        callable_here = len(usable) >= MIN_PROGRAM_GENES_FOR_CALL
        if callable_here:
            callable_fine_panels[species][program] = usable
        eligibility_rows.append({
            'program': program,
            'species': species,
            'panel_genes_present': len(genes),
            'panel_genes_usable': len(usable),
            'usable_genes': '; '.join(usable) if usable else 'none',
            'callable': callable_here,
        })
fine_program_eligibility = pd.DataFrame(eligibility_rows)
fine_program_eligibility.to_csv(CELLTYPING_DIR / 'global_dpt_fine_program_eligibility.csv', index=False)
display(fine_program_eligibility)

for species in SPECIES_GROUPS:
    if len(callable_fine_panels[species]) < 2:
        raise ValueError(
            f'Fewer than two fine programs are callable in {species}; '
            'a dominant-program call is not defined there.'
        )

fine_program_scores = pd.DataFrame(
    np.nan, index=adata_total.obs_names, columns=list(global_fine_panel)
)
for species in SPECIES_GROUPS:
    species_mask = species_for_fine_scores.eq(species).to_numpy()
    for program, usable_genes in callable_fine_panels[species].items():
        fine_program_scores.loc[species_mask, program] = _panel_mean_expression(usable_genes)[species_mask]


def _within_species_zscore(column: pd.Series) -> pd.Series:
    spread = column.std(ddof=0)
    return (column - column.mean()) / spread if spread > 0 else column * np.nan


# Panels that are not callable in a species stay NaN there, so they are excluded from that
# species' argmax instead of being ranked against a baseline they cannot meet.
standardized_fine_scores = fine_program_scores.groupby(
    species_for_fine_scores, observed=True
).transform(_within_species_zscore)
ordered_scores = np.sort(standardized_fine_scores.fillna(-np.inf).to_numpy(), axis=1)
max_program = standardized_fine_scores.idxmax(axis=1)
max_score = ordered_scores[:, -1]
margin_to_second = ordered_scores[:, -1] - ordered_scores[:, -2]
adata_total.obs['dominant_fine_reference_program'] = np.where(
    (max_score > 0) & (margin_to_second >= MARGIN_TO_SECOND), max_program, 'Unresolved'
)

# Audit 1: within-species standardisation must leave no program with a species head start.
program_species_bias = (
    standardized_fine_scores.groupby(species_for_fine_scores, observed=True).mean().T
)
program_species_bias['human_minus_mouse'] = (
    program_species_bias['human'] - program_species_bias['mouse']
)
program_species_bias.round(3).to_csv(CELLTYPING_DIR / 'global_dpt_fine_program_species_bias.csv')
display(program_species_bias.round(3))

# Audit 2: the diagnostic is not allowed to replace the reviewed labels, so it is reported
# against them instead of being presented alone.
diagnostic_vs_reviewed = (
    adata_total.obs.assign(
        diagnostic_program=adata_total.obs['dominant_fine_reference_program'].astype(str)
    )
    .groupby(['comparison_species', 'segment_class', 'diagnostic_program'], observed=True)
    .size().rename('n_structures').reset_index()
)
diagnostic_vs_reviewed.to_csv(
    CELLTYPING_DIR / 'global_dpt_fine_program_vs_reviewed_segment.csv', index=False
)
diagnostic_agreement = (
    adata_total.obs['segment_class'].astype(str)
    .eq(adata_total.obs['dominant_fine_reference_program'].astype(str))
    .groupby(adata_total.obs['comparison_species'].astype(str), observed=True).mean()
)
print('Exact agreement with the reviewed segment_class:',
      ', '.join(f'{species}={value:.1%}' for species, value in diagnostic_agreement.items()))
print('Programs not tested per species:',
      {species: sorted(set(global_fine_panel) - set(callable_fine_panels[species]))
       for species in SPECIES_GROUPS})

# The boxplot is deliberately separate from the reviewed broad-label plot above: it asks whether
# fine marker programs occupy sensible portions of the already constructed DPT. `Unresolved` is
# counted in each panel title and audited in the CSV, never drawn, so it cannot be mistaken for an
# anatomical segment.
called_fine_programs = set(adata_total.obs['dominant_fine_reference_program'].astype(str))
fine_boxplot_order = [
    program for program in FINE_NEPHRON_PROGRAM_ORDER
    if program in global_fine_panel and program in called_fine_programs
]
fine_dpt_summary = (
    adata_total.obs.groupby(['comparison_species', 'dominant_fine_reference_program'], observed=True)['total_scanpy_dpt']
    .agg(['size', 'median', 'mean']).reset_index()
)
fine_dpt_summary.to_csv(CELLTYPING_DIR / 'global_dpt_by_dominant_fine_reference_program.csv', index=False)
display(fine_dpt_summary.round(3))

fig, axes = plt.subplots(1, 2, figsize=(max(15, 1.05 * len(fine_boxplot_order)), 5.4), sharey=True)
for axis, species in zip(axes, SPECIES_GROUPS):
    species_mask = adata_total.obs['comparison_species'].astype(str).eq(species)
    data = [adata_total.obs.loc[
        species_mask & adata_total.obs['dominant_fine_reference_program'].astype(str).eq(program),
        'total_scanpy_dpt'
    ].dropna().to_numpy(dtype=float) for program in fine_boxplot_order]
    axis.boxplot([values if len(values) else np.array([np.nan]) for values in data], showfliers=False)
    axis.set_xticks(range(1, len(fine_boxplot_order) + 1))
    axis.set_xticklabels([
        f'{program}\n(n={len(values):,})' if callable_fine_panels[species].get(program)
        else f'{program}\n(not tested)'
        for program, values in zip(fine_boxplot_order, data)
    ])
    plt.setp(axis.get_xticklabels(), rotation=45, ha='right')
    n_unresolved = int(
        (species_mask & adata_total.obs['dominant_fine_reference_program'].astype(str).eq('Unresolved')).sum()
    )
    axis.set_title(f'{species.title()} (unresolved {n_unresolved:,} of {int(species_mask.sum()):,})')
    axis.set_xlabel('reference program (marker-panel argmax, not a segment label)')
    axis.set_ylabel('Global-nephron DPT' if axis is axes[0] else '')
fig.suptitle('Global DPT by dominant fine reference program (diagnostic, not a fine label)', y=1.07)
fig.text(
    0.5, 1.005,
    'panel = mean lognorm of the genes detected in >= 5% of that species; z-scored within species; '
    '>= 2 usable genes; argmax with a 0.15 SD margin; "not tested" = panel not expressed in that species',
    ha='center', va='bottom', fontsize=8.5,
)
fig.tight_layout()
fig.savefig(CELLTYPING_DIR / 'global_nephron_dpt_by_fine_reference_program.png', dpi=180, bbox_inches='tight')
plt.show()


# %%
# Purpose: save global DPT, then recompute a PT-specific DPT for downstream analyses.
# The global artifact remains the whole-nephron review result. After inspecting it, the
# reviewed PT clusters receive a new PT-only neighbor graph, root, and S1-to-S3 DPT.
DPT_EXTRA_OBS_COLS = ['comparison_species', 'region', 'x_centroid', 'y_centroid',
    'segment_class', 'coarse_class', 'leiden_coarse', 'total_marker_axis',
    'total_early_trajectory_score', 'total_late_trajectory_score', 'global_nephron_dpt',
    'shared_pseudospace',
    # carried so Section 4 can renormalise with the native denominator
    'pre_filter_total_counts']
adata_total.write(GLOBAL_DPT_OUTPUT_PATH)
print(f'Saved global nephron DPT: {GLOBAL_DPT_OUTPUT_PATH.relative_to(PROJECT_DIR)}')

# Compare the ROLLED-UP family, not the raw label: a cluster reviewed as PT-S1/PT-S2/PT-S3 is
# still the PT compartment. Comparing the raw label against 'PT' silently stopped matching once
# the review map started carrying fine segment labels.
current_pt_clusters = sorted(
    cluster for cluster, label in current_labels.items()
    if segment_vocabulary.coarse_for(label) == TRAJECTORY_COMPARTMENT
)
if not current_pt_clusters:
    raise ValueError('No current cluster is manually labeled PT; update the review map after inspecting the cluster evidence.')
adata_pt_global = adata_total[adata_total.obs['coarse_class'].astype(str).eq(TRAJECTORY_COMPARTMENT)].copy()
if adata_pt_global.n_obs < 100 or adata_pt_global.obs['comparison_species'].nunique() < 2:
    raise ValueError('Insufficient PT structures from both species after global-nephron DPT.')
tubule_filter_audit.loc[len(tubule_filter_audit)] = {
    'filter_stage': 'PT subset after global-nephron DPT inspection',
    'filter_type': 'downstream analysis subset', 'n_input_tubules': int(adata_total.n_obs),
    'n_retained_tubules': int(adata_pt_global.n_obs),
    'n_removed_tubules': int(adata_total.n_obs - adata_pt_global.n_obs),
    'fraction_removed': 1 - int(adata_pt_global.n_obs) / int(adata_total.n_obs),
}

# Preserve the global coordinate, then construct a fresh PT-only DPT. Fine PT marker
# modules orient this coordinate but do not assign PT cells or create fine clusters.
adata_pt_global.obs['global_nephron_dpt'] = adata_pt_global.obs['total_scanpy_dpt'].to_numpy(dtype=float)
# `recompute_subset_dpt` accepts marker records so aliases could be supplied;
# wrap these reference gene strings in that format without changing the marker list itself.
pt_dpt_marker_modules = {
    group: [{'label': gene, 'candidates': [gene]} for gene in genes]
    for group, genes in NEPHRON_AXIS_MARKERS.items()
}
pt_dpt_diagnostic = recompute_subset_dpt(
    adata_pt_global, 'PT', pt_dpt_marker_modules,
    ['PT-S1', 'PT-S2', 'PT-S3'], output_col='pt_specific_dpt',
    n_neighbors=N_NEIGHBORS, random_state=RANDOM_STATE,
)
adata_pt_global.obs['total_scanpy_dpt'] = adata_pt_global.obs['pt_specific_dpt'].to_numpy(dtype=float)
adata_pt_global.obs['shared_pseudospace'] = adata_pt_global.obs['total_scanpy_dpt'].to_numpy(dtype=float)
finite_pt_dpt = np.isfinite(adata_pt_global.obs['total_scanpy_dpt'].to_numpy(dtype=float))
if not finite_pt_dpt.all():
    print(f'Removing {int((~finite_pt_dpt).sum()):,} PT structures disconnected from the PT-specific DPT root.')
    adata_pt_global = adata_pt_global[finite_pt_dpt].copy()
tubule_filter_audit.loc[len(tubule_filter_audit)] = {
    'filter_stage': 'connected to PT-specific DPT root', 'filter_type': 'downstream analysis subset',
    'n_input_tubules': int(finite_pt_dpt.size), 'n_retained_tubules': int(adata_pt_global.n_obs),
    'n_removed_tubules': int((~finite_pt_dpt).sum()),
    'fraction_removed': float((~finite_pt_dpt).mean()),
}
tubule_filter_audit.to_csv(DIAGNOSTIC_DIR / 'tubule_filter_audit.csv', index=False)
pt_dpt_diagnostic_frame = pd.DataFrame([pt_dpt_diagnostic])
pt_dpt_diagnostic_frame.to_csv(CELLTYPING_DIR / 'pt_specific_dpt_diagnostics.csv', index=False)
display(pt_dpt_diagnostic_frame)

sc.pl.embedding(adata_pt_global, basis='umap', color=['total_scanpy_dpt', 'comparison_species'],
                cmap='viridis', frameon=False, ncols=2, show=False)
plt.savefig(CELLTYPING_DIR / 'pt_specific_dpt_umap.png', dpi=180, bbox_inches='tight')
plt.show()

adata_dpt_saved = save_total_pseudotime_anndata(adata_pt_global, 'total_scanpy_dpt', DPT_OUTPUT_PATH,
    expression_filtered_adata=adata_all_expression_filtered, annotation_adata=adata_pt_global,
    project_dir=PROJECT_DIR, forbidden_labels=REMOVE_CLASSES + ['nan', 'none'],
    assert_no_nan=True, extra_obs_cols=DPT_EXTRA_OBS_COLS)
positions = adata_pt_global.obs_names.get_indexer(adata_dpt_saved.obs_names)
if (positions < 0).any(): raise RuntimeError('PT DPT output contains structures absent from the PT-specific subset.')
for column in DPT_EXTRA_OBS_COLS:
    if column in adata_pt_global.obs.columns:
        values = adata_pt_global.obs[column].to_numpy()[positions]
        adata_dpt_saved.obs[column] = values.astype(float) if pd.api.types.is_numeric_dtype(adata_pt_global.obs[column]) else values.astype(str)
for key in ('X_harmony', 'X_pca', 'X_umap'):
    if key in adata_pt_global.obsm: adata_dpt_saved.obsm[key] = np.asarray(adata_pt_global.obsm[key])[positions]
adata_dpt_saved.write(DPT_OUTPUT_PATH)
print(f'Saved PT-specific DPT: {DPT_OUTPUT_PATH.relative_to(PROJECT_DIR)}')
print(f'Global root cluster: {root_cluster}; global retained: {adata_total.n_obs:,}; PT-specific DPT: {adata_dpt_saved.n_obs:,}')


# %% [markdown]
# # Section 3 - PT-specific marker-gradient heatmaps
#
# These plots use a DPT recomputed only within the reviewed PT clusters after global-nephron
# DPT review. The fine S1/S2/S3 markers orient that coordinate but do not create fine labels.
#

# %% [markdown]
# ## 3.1 - heatmap settings
#
# **Purpose:** define heatmap-specific marker panels, colors, and binning at point of use.
#

# %%
# Purpose: # Marker panels and display settings used only by the heatmap section.
# Marker panels and display settings used only by the heatmap section.
HEATMAP_MARKERS = {
    'PT': {
        'PT-S1': ['Slc5a2', 'Slc5a12', 'Gatm'],
        'PT-S2': ['Slc22a6', 'Slc13a3', 'Cyp2e1'],
        'PT-S3': ['Slc22a7', 'Cyp7b1'],
    },
}
HEATMAP_COLORS = {'PT-S1': '#4C9BD3', 'PT-S2': '#5B8E55', 'PT-S3': '#D6B48A', 'ATL': '#26A69A', 'mTAL': '#F58518', 'cTAL': '#E8A33D', 'DCT1': '#D81B60', 'DCT2': '#EC6EA5', 'CNT': '#76B7B2', 'CCD': '#4E79A7', 'OMCD': '#9C755F', 'IMCD': '#593C8F'}
SEGMENT_STRIP_COLORS = {**HEATMAP_COLORS, 'PT': '#4C9BD3', 'DTL': '#7E57C2', 'AL': '#F58518', 'DCT': '#D81B60', 'CNT_CD': '#8D6E63', 'Glomerulus': '#9E9E9E'}
N_BINS = 120
HEATMAP_OUTPUT_DIR = RESULTS_DIR / 'heatmaps'; HEATMAP_OUTPUT_DIR.mkdir(exist_ok=True)
def _marker_dict(markers):
    return {group: [{'label': gene, 'candidates': [gene]} for gene in genes] for group, genes in markers.items()}
print(f'Heatmap families: {", ".join(HEATMAP_MARKERS)}')


# %%
# Purpose: adata_heatmap = sc.read_h5ad(DPT_OUTPUT_PATH)
adata_heatmap = sc.read_h5ad(DPT_OUTPUT_PATH)
if 'lognorm' in adata_heatmap.layers:
    adata_heatmap.X = adata_heatmap.layers['lognorm'].copy()

global_group_order = [
    group for family in ('PT',)
    for group in HEATMAP_MARKERS[family]
]
global_panel = {
    group: values
    for family in ('PT',)
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
        title=f'{species.title()} PT markers on the PT-specific DPT',
        output_name=f'{species}_total_marker_heatmap.png',
        strip_col='segment_class',
        strip_order=SEGMENT_DISPLAY_ORDER,
        strip_colors=SEGMENT_STRIP_COLORS,
        n_bins=N_BINS,
        output_dir=HEATMAP_OUTPUT_DIR,
        project_dir=PROJECT_DIR,
    )

print('Completed:', 'adata_heatmap = sc.read_h5ad(DPT_OUTPUT_PATH)')


# %% [markdown]
# # Section 4 - human versus healthy-mouse PT trajectories
#
# The same nested level/shape GAM used in the mouse notebook is fitted on PT structures over the
# mouse/human common DPT support. Here c=0 is mouse and c=1 is human. Both human slices are healthy
# cortex from one donor, so the fitted differences are descriptive effect sizes only. The two human
# slices are shown separately as sensitivity curves but are not treated as independent biological
# replicates.
#
#
# **Output:** descriptive (not inferential) gene and pathway trajectory comparisons, with figure files in `curves/`.
#
# Peak positions are coordinates on THIS notebook's DPT construction (the PT-specific recomputed DPT), so they are not comparable with notebook 02's global-nephron PT coordinate.
#

# %%
# Purpose: from pseudospace.levelshape import fit_single_condition_curves, run_level_shape, summarize_curve_effects
from pseudospace.levelshape import (
    fit_single_condition_curves,
    run_level_shape,
    summarize_curve_effects,
)
from pseudospace.pathways import (
    build_pathway_membership,
    member_gene_evidence,
    summarize_pathway_redundancy,
)
from pseudospace.stats_gam import as_csr, gam_internal_knots, resolve_present

# This input is needed only for the pathway comparison below.
PATHWAY_LIBRARY_DIR = DATA_ROOT / 'mouse_vs_human' / 'pathway_gene_sets'
GAM_LAMBDA_GRID = np.logspace(-3, 3, 13)

SECTION4_CONFIG = {
    'n_internal_knots': 9,
    'lambda_grid': GAM_LAMBDA_GRID,
    'common_support_pct': (1, 99),
    'grid_points': 101,
    'min_detected_fraction': 0.02,
    'min_mean_expression': 0.0,
    'pathway_min_genes': 10,
    # No longer applied to membership: kept so the walkthrough can report which pathways an
    # upper bound of this size would have removed.
    'pathway_max_genes': 150,
    'pathway_overlap_threshold': 0.6,
    'redundancy_prioritized_pathways': 40,
    'pseudobulk_bins': 20,
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
# Purpose: adata_full = sc.read_h5ad(DPT_OUTPUT_PATH)
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
#
#

# %%
# Purpose: gene_fit = run_level_shape(
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
# `species_effect_rms` alone mixes a constant vertical offset with a redistribution along
# pseudospace, and `shape_rms` alone reports a weaker gradient as a shape change. The effect
# summary separates the three questions: `level_fraction` is the share of the squared total effect
# that is vertical offset, `pattern_rms_z` is free of both level and amplitude, and the amplitude
# columns describe gradient strength on its own.
gene_curve_effects = summarize_curve_effects(
    gene_fit['curve_healthy'], gene_fit['curve_aki'], feature_names=gene_names
)
gene_results = gene_results.merge(
    gene_curve_effects.rename(columns={'feature': 'gene'}), on='gene', how='left'
)
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

# Normalisation sensitivity: the primary fit divides by the retained panel's total. Re-fit with the
# native (pre-filter) denominator recorded in the QC cell, so a level-dominated ranking cannot be an
# artifact of the panel denominator alone.
if ('counts' in adata_pt.layers) and ('pre_filter_total_counts' in adata_pt.obs.columns):
    native_scale = NORMALIZE_TARGET_SUM / np.maximum(
        adata_pt.obs['pre_filter_total_counts'].to_numpy(dtype=float), 1.0
    )
    Y_native = sparse.diags(native_scale) @ as_csr(adata_pt.layers['counts'])[:, tested]
    Y_native.data = np.log1p(Y_native.data)
    native_fit = run_level_shape(Y_native, s, c, knots, grid, SECTION4_CONFIG['lambda_grid'])
    gene_results['level_effect_native_denominator'] = native_fit['level_effect']
    gene_results['shape_rms_native_denominator'] = native_fit['shape_rms']
    top_native = gene_results[gene_results['primary_eligible']].head(20)
    normalisation_sensitivity = pd.DataFrame([{
        'level_effect_spearman_panel_vs_native': float(spearmanr(
            gene_results['level_effect_human_minus_mouse'],
            gene_results['level_effect_native_denominator'],
        ).correlation),
        'top20_sign_agreement': float((
            np.sign(top_native['level_effect_human_minus_mouse'])
            == np.sign(top_native['level_effect_native_denominator'])
        ).mean()),
        'top20_median_level_fraction_panel_denominator': float(top_native['level_fraction'].median()),
        'top20_median_level_fraction_native_denominator': float(np.median(
            (top_native['level_effect_native_denominator'] ** 2)
            / np.maximum(top_native['species_effect_rms'] ** 2, 1e-12)
        )),
    }])
    normalisation_sensitivity.to_csv(CURVE_OUTPUT_DIR / 'normalisation_sensitivity.csv', index=False)
    display(normalisation_sensitivity.round(3))
    print('Agreement between the two denominators means the vertical-offset finding survives the '
          'retained-panel denominator; it does not make the species offset biological.')
else:
    print('Normalisation sensitivity skipped: PT object lacks layers["counts"] or '
          'obs["pre_filter_total_counts"].')
gene_results.to_csv(
    CURVE_OUTPUT_DIR / 'gene_trajectory_comparison.csv', index=False
)
display(gene_results.head(20)[[
    'gene', 'species_effect_rms', 'level_effect_human_minus_mouse',
    'shape_rms', 'level_fraction', 'pattern_rms_z', 'difference_type',
    'curve_spearman', 'axis_basis_gene',
]])

# Magnitude audit. Every structure is normalised to the same total, so a broad vertical offset can
# dominate the ranking; say how much of it is level rather than spatial redistribution instead of
# leaving the reader to infer it from the effect size.
top_effect_genes = gene_results[gene_results['primary_eligible']].head(20)
print('Top-20 eligible genes: median level fraction '
      f"{top_effect_genes['level_fraction'].median():.3f}; "
      f"lower in human {int((top_effect_genes['level_effect_human_minus_mouse'] < 0).sum())}/20; "
      f"types {top_effect_genes['difference_type'].value_counts().to_dict()}")



# %%
# Purpose: top_gene_rows = gene_results[gene_results['primary_eligible']].head(
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


print('Completed:', 'top_gene_rows = gene_results[gene_results["primary_eligible"]].head(')


# %%
# Purpose: fig, axes = _grid_axes(len(top_genes), ncols=3)
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


print('Completed:', 'fig, axes = _grid_axes(len(top_genes), ncols=3)')


# %%
# Purpose: # Paired standardized fitted-curve heatmap for the broader collaborator shortlist.
# Paired standardized fitted-curve heatmap for the broader collaborator shortlist.
heatmap_rows = gene_results[gene_results['primary_eligible']].head(
    SECTION4_CONFIG['heatmap_genes']
)
heatmap_idx = np.array(
    [gene_lookup_tested[gene.upper()] for gene in heatmap_rows['gene']], dtype=int
)
mouse_heat = gene_fit['curve_healthy'][heatmap_idx]
human_heat = gene_fit['curve_aki'][heatmap_idx]

# Two coordinated views of the same genes, because one heatmap cannot answer both questions:
#   (1) common scale - fitted curves in lognorm units with ONE colour scale for both species, so a
#       species-wide level or amplitude difference stays visible;
#   (2) shape only - each curve standardised over the grid, so peak position and profile similarity
#       are readable even where one species expresses the gene lower overall.
# A single pooled z-score (the previous version) mixes the two: it removes the between-species level
# difference by construction and then hides it inside the standardisation.
common_vmin = float(np.nanmin([mouse_heat.min(), human_heat.min()]))
common_vmax = float(np.nanmax([mouse_heat.max(), human_heat.max()]))
fig, axes = plt.subplots(1, 2, figsize=(9.5, max(6, 0.19 * len(heatmap_rows))), sharey=True)
for axis, values, title in zip(axes, (mouse_heat, human_heat), ('Healthy mouse', 'Human')):
    image = axis.imshow(
        values, aspect='auto', interpolation='nearest', cmap='viridis',
        vmin=common_vmin, vmax=common_vmax,
        extent=[grid[0], grid[-1], len(heatmap_rows) - 0.5, -0.5],
    )
    axis.set_title(title)
    axis.set_xlabel('Shared PT DPT')
axes[0].set_yticks(np.arange(len(heatmap_rows)))
axes[0].set_yticklabels(heatmap_rows['gene'])
fig.colorbar(image, ax=axes, label='Fitted log-normalised expression (shared scale)', shrink=0.65)
fig.suptitle('Paired fitted PT trajectories on a common expression scale\n'
             'level and amplitude preserved; colour is comparable across species')
_save_figure(fig, 'paired_top_gene_heatmap_lognorm.png')

def _shape_only(curves):
    mean = np.nanmean(curves, axis=1, keepdims=True)
    spread = np.nanstd(curves, axis=1, keepdims=True)
    spread[~np.isfinite(spread) | (spread < 1e-8)] = 1.0
    return np.clip((curves - mean) / spread, -2.5, 2.5)

fig, axes = plt.subplots(1, 2, figsize=(9.5, max(6, 0.19 * len(heatmap_rows))), sharey=True)
for axis, values, title in zip(
    axes, (_shape_only(mouse_heat), _shape_only(human_heat)), ('Healthy mouse', 'Human')
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
fig.colorbar(image, ax=axes, label='Within-curve z-score over pseudospace', shrink=0.65)
fig.suptitle('Same genes, shape only\nlevel and amplitude removed: peak position and width only')
_save_figure(fig, 'paired_top_gene_heatmap_shape.png')


print('Completed:', '# Paired standardized fitted-curve heatmap for the broader collaborato')


# %% [markdown]
# ## 4.2 - canonical PT marker-gradient check
#
# These axis-adjacent genes are displayed separately from the discovery ranking. Their purpose is
# to show whether S1, S2, and S3 positional programs occupy comparable locations in the shared
# coordinate, not to claim species differences.
#
#

# %%
# Purpose: pt_marker_groups = {
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


print('Completed:', 'pt_marker_groups = {')


# %% [markdown]
# ## 4.3 - predefined pathway/module trajectories
#
# Pathway scores are the mean of PT-standardized member-gene expression. All eligible predefined
# Hallmark, KEGG Mouse, and Reactome sets are modeled; pathways are ranked by descriptive
# mouse-human fitted-curve separation.
#
#

# %%
# Purpose: gene_mean = np.asarray(Y_genes.mean(axis=0)).ravel()
gene_mean = np.asarray(Y_genes.mean(axis=0)).ravel()
gene_sq = np.asarray(Y_genes.multiply(Y_genes).mean(axis=0)).ravel()
gene_std = np.sqrt(np.maximum(gene_sq - gene_mean ** 2, 1e-12))

# Pathway membership through the same accepted ortholog map that built the matrix, with every
# coverage stage recorded. `adata_pt.var_names` is the assayed matrix and `gene_names` the genes
# that passed the detection filter, so "never in the library", "no accepted ortholog", "assayed but
# filtered" and "used in the score" stay distinguishable. Matching by uppercase alone cannot.
pathway_coverage = []
pathway_mapping_rows = []
for library in SECTION4_CONFIG['pathway_libraries']:
    library_path = PATHWAY_LIBRARY_DIR / f'{library}.json'
    if not library_path.exists():
        print(f'WARNING: missing pathway library {library_path.name}')
        continue
    library_sets = json.loads(library_path.read_text())
    membership = build_pathway_membership(
        library_sets,
        adata_pt.var_names,
        ortholog_map=ortholog_map,
        library_name=library,
        min_genes=SECTION4_CONFIG['pathway_min_genes'],
        # The old 150-member cap removed whole Hallmark sets before they were ever scored. Pathways
        # are kept and the ones an upper bound would remove are counted and named instead.
        max_genes=None,
        tested=gene_names,
    )
    pathway_coverage.append(membership)
    pathway_mapping_rows.append({
        'library': library,
        'n_pathways_total': int(len(membership)),
        'n_pathways_retained': int(membership['retained'].sum()),
        'n_pathway_members_requested': int(membership['n_requested'].sum()),
        'n_pathway_members_assayed': int(membership['n_assayed'].sum()),
        'n_pathway_members_with_ortholog': int(membership['n_with_ortholog'].sum()),
        'n_pathways_over_previous_cap': int(
            (membership['n_assayed'] > SECTION4_CONFIG['pathway_max_genes']).sum()
        ),
    })

pathway_coverage = pd.concat(pathway_coverage, ignore_index=True)
pathway_coverage.to_csv(CURVE_OUTPUT_DIR / 'pathway_membership_coverage.csv', index=False)
pathway_mapping = pd.DataFrame(pathway_mapping_rows)
pathway_mapping.to_csv(CURVE_OUTPUT_DIR / 'pathway_symbol_mapping_report.csv', index=False)
display(pathway_mapping)

retained_pathways = pathway_coverage[pathway_coverage['retained']].reset_index(drop=True)
print(f'Pathways retained: {len(retained_pathways):,} of {len(pathway_coverage):,}; excluded '
      f'{int((~pathway_coverage["retained"]).sum()):,} '
      f'{pathway_coverage.loc[~pathway_coverage["retained"], "exclusion_reason"].value_counts().to_dict()}')
print(f'Pathways a {SECTION4_CONFIG["pathway_max_genes"]}-member upper bound would remove: '
      f'{int((pathway_coverage["n_assayed"] > SECTION4_CONFIG["pathway_max_genes"]).sum()):,}')

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
    # Same level/amplitude/pattern split as for genes; module scores are already standardised, so
    # the level term here is the species offset in pooled-z units.
    pathway_curve_effects = summarize_curve_effects(
        pathway_fit['curve_healthy'], pathway_fit['curve_aki'],
        feature_names=[f'{library}: {pathway}' for library, pathway in
                       zip(retained_pathways['library'], retained_pathways['pathway'])],
    )
    for column in ('level_fraction', 'shape_fraction', 'pattern_rms_z', 'amplitude_ratio',
                   'amplitude_log2_ratio', 'difference_type'):
        pathway_results[column] = pathway_curve_effects[column].to_numpy()
    pathway_results = pathway_results.sort_values(
        'species_effect_rms', ascending=False
    ).reset_index(drop=True)
    top_effect_pathways = pathway_results.head(20)
    print('Top-20 pathways: median level fraction '
          f"{top_effect_pathways['level_fraction'].median():.3f}; "
          f"lower in human {int((top_effect_pathways['level_effect_human_minus_mouse'] < 0).sum())}/20; "
          f"types {top_effect_pathways['difference_type'].value_counts().to_dict()}")
pathway_results.to_csv(
    CURVE_OUTPUT_DIR / 'pathway_trajectory_comparison.csv', index=False
)
display(pathway_results.head(20)[[
    'library', 'pathway', 'n_genes_present', 'species_effect_rms',
    'level_effect_human_minus_mouse', 'shape_rms', 'level_fraction', 'pattern_rms_z',
    'difference_type',
]] if len(pathway_results) else pathway_results)

# Member-gene evidence beside every prioritized pathway: how many members move with the aggregate
# trend, and whether a single gene carries it. A pathway name is not evidence of a distinct program
# when its member set is another pathway's member set.
if len(pathway_results):
    prioritized_pathways = pathway_results.head(
        SECTION4_CONFIG['redundancy_prioritized_pathways']
    ).copy()
    prioritized_pathways['member_list'] = prioritized_pathways['genes_present'].map(
        lambda value: [gene for gene in str(value).split('; ') if gene]
    )

    # Redundancy is grouped on the prioritized subset, not the whole library: Reactome is a
    # hierarchy of nested parent/child sets, so a library-wide grouping collapses into one group
    # and says nothing about the specific claims on the shortlist.
    redundancy_pairs, redundancy_groups = summarize_pathway_redundancy(
        prioritized_pathways,
        gene_column='member_list',
        label_columns=('library', 'pathway'),
        overlap_threshold=SECTION4_CONFIG['pathway_overlap_threshold'],
    )
    redundancy_pairs.to_csv(CURVE_OUTPUT_DIR / 'pathway_redundancy_pairs.csv', index=False)
    redundancy_groups.to_csv(CURVE_OUTPUT_DIR / 'pathway_redundancy_groups.csv', index=False)
    print(f'Overlapping prioritized pairs (>= '
          f'{SECTION4_CONFIG["pathway_overlap_threshold"]:.0%} of the smaller member set): '
          f'{len(redundancy_pairs):,}; largest group: '
          f'{int(redundancy_groups["group_size"].max()) if len(redundancy_groups) else 0}')
    if len(redundancy_pairs):
        display(redundancy_pairs.sort_values('n_shared', ascending=False).head(10))

    top_pathway_evidence = []
    for row in pathway_results.head(SECTION4_CONFIG['top_pathway_plots']).itertuples():
        members = [gene for gene in str(row.genes_present).split('; ') if gene]
        evidence = member_gene_evidence(
            members, gene_results, effect_column='level_effect_human_minus_mouse'
        )
        top_pathway_evidence.append({
            'library': row.library,
            'pathway': row.pathway,
            'species_effect_rms': row.species_effect_rms,
            'level_effect_human_minus_mouse': row.level_effect_human_minus_mouse,
            **evidence.to_dict(),
        })
    top_pathway_evidence = pd.DataFrame(top_pathway_evidence)
    top_pathway_evidence.to_csv(
        CURVE_OUTPUT_DIR / 'pathway_top_member_evidence.csv', index=False
    )
    display(top_pathway_evidence[[
        'pathway', 'n_members_present', 'fraction_members_agreeing', 'strongest_member',
        'strongest_member_effect', 'member_effect_mean',
        'member_effect_mean_without_strongest', 'sign_flips_without_strongest',
    ]].round(3))


# %%
# Purpose: if len(pathway_results):
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


print('Completed:', 'if len(pathway_results):')


# %% [markdown]
# ## 4.4 - Curve modules: positional programs and response programs
#
# A module is a group of genes that share a *shape*, not a peak. Two catalogs are discovered from the
# fitted curves:
#
# * **positional** - cluster the reference (healthy mouse) curves: which genes share a normal PT pattern.
# * **response** - cluster the mean-centred difference curves `Delta f(s) = f_human(s) - f_mouse(s)`:
#   which genes change in the same way, independent of how big their overall offset is. Broad
#   suppression, selective loss of the late program, peak relocation and gradient flattening then
#   separate instead of all being "large effect".
#
# Clustering is on the correlation distance between grid-standardised curves, so level and amplitude do
# not decide membership; peak position, width, monotonicity and amplitude annotate a module afterwards.
# Peaks alone are deliberately not the grouping key (a narrow spike and a broad plateau would merge),
# and dynamic time warping is not used because aligning peaks would erase the positional difference
# being measured. Stability is reported under specimen omission (the balanced curves recomputed
# without each specimen) and under a different data mixture (pooled fit vs specimen-balanced), and
# module-pathway overlap is tested against the genes **eligible for module discovery**, BH-corrected
# across every (module, gene set) pair.
#
# Every member list and the per-gene assignment are written to `curves/module_*.csv`.
#

# %%
# Purpose: discover positional and response curve modules, with stability and enrichment.
from pseudospace.modules import (
    difference_curves,
    discover_curve_modules,
    enrich_modules,
    module_stability,
)
from pseudospace.specimen import specimen_balanced_curves

# Per-specimen curves for EVERY tested gene (not only the plotted shortlist): the pooled fit cannot
# be used one specimen at a time because the condition indicator is constant within a specimen.
specimen_full_curves = {}
for sample_name in sorted(set(samples)):
    mask = samples == sample_name
    curves, _ = fit_single_condition_curves(
        Y_genes[mask], s[mask], knots, grid, SECTION4_CONFIG['lambda_grid'],
        gene_fit['lam_idx'], support_pct=SECTION4_CONFIG['common_support_pct'],
    )
    specimen_full_curves[sample_name] = curves

mouse_specimens = [name for name in specimen_full_curves if name in MOUSE_SAMPLES]
human_specimens = [name for name in specimen_full_curves if name in HUMAN_SAMPLES]
balanced_reference = specimen_balanced_curves({k: specimen_full_curves[k] for k in mouse_specimens})
balanced_comparison = specimen_balanced_curves({k: specimen_full_curves[k] for k in human_specimens})
specimen_balanced_frame = pd.DataFrame({
    'gene': gene_names,
    'balanced_mouse_amplitude': np.nanmax(balanced_reference, axis=1) - np.nanmin(balanced_reference, axis=1),
    'balanced_human_amplitude': np.nanmax(balanced_comparison, axis=1) - np.nanmin(balanced_comparison, axis=1),
    'balanced_level_effect_human_minus_mouse': np.nanmean(balanced_comparison - balanced_reference, axis=1),
})
specimen_balanced_frame.to_csv(CURVE_OUTPUT_DIR / 'specimen_balanced_gene_amplitudes.csv', index=False)

MODULE_CONFIG = {
    'max_distance': 0.4,
    'min_amplitude': 0.05,
    'min_features': 10,
}
positional_labels, positional_modules = discover_curve_modules(
    gene_fit['curve_healthy'], grid,
    max_distance=MODULE_CONFIG['max_distance'],
    min_amplitude=MODULE_CONFIG['min_amplitude'],
    min_features=MODULE_CONFIG['min_features'],
    feature_names=gene_names,
)
response_delta = difference_curves(gene_fit['curve_healthy'], gene_fit['curve_aki'])
response_centered = difference_curves(gene_fit['curve_healthy'], gene_fit['curve_aki'], center=True)
response_labels, response_modules = discover_curve_modules(
    response_centered, grid,
    max_distance=MODULE_CONFIG['max_distance'],
    min_amplitude=MODULE_CONFIG['min_amplitude'],
    min_features=MODULE_CONFIG['min_features'],
    feature_names=gene_names,
)

gene_module_table = pd.DataFrame({
    'gene': gene_names,
    'positional_module': positional_labels.to_numpy(),
    'response_module': response_labels.to_numpy(),
    'delta_mean_human_minus_mouse': np.nanmean(response_delta, axis=1),
    'delta_pattern_rms': np.sqrt(np.nanmean(
        (response_centered - np.nanmean(response_centered, axis=1, keepdims=True)) ** 2, axis=1)),
})
gene_module_table.to_csv(CURVE_OUTPUT_DIR / 'module_gene_assignment.csv', index=False)
positional_modules.to_csv(CURVE_OUTPUT_DIR / 'module_positional_catalog.csv', index=False)
response_modules.to_csv(CURVE_OUTPUT_DIR / 'module_response_catalog.csv', index=False)
display(positional_modules)
display(response_modules)

# Stability: same discovery on curves that drop one specimen (from the per-specimen fits) and on a
# different data mixture (specimen-balanced instead of pooled).
def _response_modules_from(curve_reference, curve_comparison):
    labels, _ = discover_curve_modules(
        difference_curves(curve_reference, curve_comparison, center=True), grid,
        max_distance=MODULE_CONFIG['max_distance'],
        min_amplitude=MODULE_CONFIG['min_amplitude'],
        min_features=MODULE_CONFIG['min_features'], feature_names=gene_names,
    )
    return labels

stability_runs = {'pooled_fit': response_labels}
for dropped in sorted(specimen_full_curves):
    kept = [name for name in specimen_full_curves if name != dropped]
    kept_mouse = [name for name in kept if name in MOUSE_SAMPLES]
    kept_human = [name for name in kept if name in HUMAN_SAMPLES]
    if not kept_mouse or not kept_human:
        continue                      # dropping this specimen collapses a side of the comparison
    stability_runs[f'without_{dropped}'] = _response_modules_from(
        specimen_balanced_curves({k: specimen_full_curves[k] for k in kept_mouse}),
        specimen_balanced_curves({k: specimen_full_curves[k] for k in kept_human}),
    )
stability_runs['specimen_balanced_mixture'] = _response_modules_from(
    balanced_reference, balanced_comparison
)
module_stability_table = module_stability(stability_runs)
module_stability_table.to_csv(CURVE_OUTPUT_DIR / 'module_response_stability.csv', index=False)
display(module_stability_table.round(3))

# Enrichment against the pathway library, with the discovery-eligible genes as the background and
# correction across every tested (module, gene set) pair.
module_gene_sets = {
    f'{row.library}: {row.pathway}': row.genes_present
    for row in pathway_coverage[pathway_coverage['retained']].itertuples()
}
module_enrichment = enrich_modules(
    {module: members.index for module, members in response_labels.groupby(response_labels, observed=True)
     if module != 'unassigned'},
    module_gene_sets,
    background=gene_names,
)
module_enrichment.to_csv(CURVE_OUTPUT_DIR / 'module_pathway_enrichment.csv', index=False)
print('Response modules:', int((response_labels != 'unassigned').sum()), 'genes in',
      int(response_labels[response_labels != 'unassigned'].nunique()), 'modules;',
      'positional modules:', int(positional_labels[positional_labels != 'unassigned'].nunique()))
print('Module-pathway pairs tested:', len(module_enrichment), '; significant after BH (q<0.05):',
      int((module_enrichment['p_value_adjusted'] < 0.05).sum()))
display(module_enrichment.head(15)[[
    'module', 'gene_set', 'n_module_genes', 'n_overlap', 'expected_overlap',
    'p_value', 'p_value_adjusted',
]].round(4))

# A broad suppression is not a module of its own in the clustered catalog: report those genes
# separately so "everything went down together" cannot be read as a discovered program.
level_only = gene_module_table[
    (gene_results.set_index('gene').reindex(gene_module_table['gene'])['difference_type'].to_numpy()
     == 'level shift')
]
level_only.to_csv(CURVE_OUTPUT_DIR / 'module_level_only_genes.csv', index=False)
print(f'Genes whose difference is a level shift rather than a pattern change: {len(level_only):,}; '
      'module discovery is run on the mean-centred difference, so they are annotated, not clustered.')


# %% [markdown]
# ## 4.5 - Specimen-level summaries and coordinate robustness
#
# The pooled GAMs weight tubules, so a specimen with more structures contributes more everywhere and
# the specimen mixture can drift along pseudospace. The specimen-balanced curve is the primary
# descriptive summary from here on, with the individual specimen curves kept visible in the gene
# figures. Pseudobulk profiles (raw counts summed per specimen and pseudospace bin) are written for
# count-based or specimen-level modelling; bins from one specimen are repeated measurements along a
# single coordinate, not replicates, so the number of independent units stays the number of specimens.
#
# Coordinate robustness: dropping a specimen and recomputing the global DPT on the remainder shows how
# much the coordinate depends on any one specimen. A transcriptome-level feature exclusion would need
# HVG selection and Harmony to be re-run, which is why the mouse workflow's physical-axis sensitivity
# analysis remains the feature-exclusion check.
#

# %%
# Purpose: specimen-level summaries and coordinate robustness.
from pseudospace.specimen import coordinate_agreement, pseudobulk_profiles

# Pseudobulk profiles from raw counts: one row per (specimen, pseudospace bin). Bins from the same
# specimen are repeated measurements along one coordinate, not replicates.
pseudobulk = pseudobulk_profiles(
    adata_pt.layers['counts'], samples, s,
    n_bins=SECTION4_CONFIG['pseudobulk_bins'],
    gene_names=list(adata_pt.var_names),
)
# Keep the wide shape (one row per specimen and bin, one column per gene): the long form is the same
# information at 10,000x the rows and nothing downstream needs it.
pseudobulk.to_csv(CURVE_OUTPUT_DIR / 'pseudobulk_profiles.csv', index=False)
pseudobulk_columns = [column for column in pseudobulk if column.startswith('count_')]
print(f'Pseudobulk: {len(pseudobulk)} specimen-bin rows x {len(pseudobulk_columns)} genes; '
      f'{pseudobulk["specimen"].nunique()} specimens; '
      f'median {int(pseudobulk["n_structures"].median())} structures per row.')
display(pseudobulk[[column for column in pseudobulk.columns
                    if not column.startswith('count_')]].describe().round(3))

# Coordinate robustness: drop each specimen, rebuild neighbours + diffusion map + DPT on the rest with
# the same root rule, and compare the resulting coordinate with the full one.
coordinate_rows = []
for dropped in sorted(set(samples)):
    keep = samples != dropped
    if len(set(samples[keep])) < 2 or np.unique(c[keep]).size < 2:
        continue
    subset = adata_pt[keep].copy()
    trajectory_key = 'robustness_neighbors'
    sc.pp.neighbors(subset, use_rep='X_harmony', n_neighbors=N_NEIGHBORS, key_added=trajectory_key,
                    random_state=RANDOM_STATE)
    sc.tl.diffmap(subset, neighbors_key=trajectory_key, random_state=RANDOM_STATE)
    # Root at the structure that currently sits at the lowest DPT, so the recomputed coordinate is
    # compared from the same starting point rather than from an unrelated one.
    subset.uns['iroot'] = int(np.argmin(subset.obs['total_scanpy_dpt'].to_numpy(dtype=float)))
    sc.tl.dpt(subset, neighbors_key=trajectory_key)
    recomputed = np.asarray(subset.obs['dpt_pseudotime'], dtype=float)
    reference = adata_pt.obs['total_scanpy_dpt'].to_numpy(dtype=float)[keep]
    agreement = coordinate_agreement(reference, recomputed)
    agreement.update({'dropped_sample': dropped, 'n_structures': int(keep.sum())})
    coordinate_rows.append(agreement)
coordinate_robustness = pd.DataFrame(coordinate_rows)
if len(coordinate_robustness):
    coordinate_robustness.to_csv(CURVE_OUTPUT_DIR / 'coordinate_specimen_robustness.csv', index=False)
    display(coordinate_robustness.round(3))
    print('Specimen omission keeps the coordinate if the rank concordance stays high; a low value for '
          'one specimen means that specimen drives the axis.')
else:
    print('Coordinate robustness skipped: no specimen can be dropped without collapsing a condition.')


# %% [markdown]
# ## 4.6 - Signed rankings and competitive enrichment
#
# Enrichment is only meaningful per question. One unsigned `shape_rms` ranking cannot tell "late-PT
# induction" from "loss of an early-PT program", so five signed per-gene statistics are ranked
# separately - level offset, amplitude change, redistribution towards late pseudospace, and the early
# and late contrasts - and each is tested competitively against the pathway library with the
# **discovery-eligible genes as the background**.
#
# Two p-values are reported per (question, gene set):
#
# * a label-permutation p over random gene sets of the same size, which ignores within-set correlation,
# * a correlation-aware z-test that inflates the variance of the set mean by the average squared
#   correlation of the set's genes measured in the expression matrix, which is the CAMERA point: genes
#   in a pathway are not independent, and the naive mean overstates the evidence.
#
# Both are corrected across every tested pair. A pathway that shares most of its members with another
# pathway is already flagged in the redundancy tables, and a pathway claim resting on one gene is
# visible in the member-gene evidence, so neither has to be discovered again from an enrichment list.
# Signed downstream signatures (PROGENy-style) are scored when a gene-to-weight matrix is supplied;
# none is bundled, so that step is skipped unless one is present.
#

# %%
# Purpose: signed rankings and competitive (correlation-aware) enrichment per question.
from pseudospace.enrichment import (
    camera_like_enrichment,
    score_signed_signatures,
    signed_gene_rankings,
)

signed_rankings = signed_gene_rankings(
    gene_fit['curve_healthy'], gene_fit['curve_aki'], grid, gene_names=gene_names
)
signed_rankings = signed_rankings.merge(
    gene_results[['gene', 'level_fraction', 'shape_fraction', 'pattern_rms_z', 'difference_type']],
    on='gene', how='left',
)
signed_rankings.to_csv(CURVE_OUTPUT_DIR / 'gene_signed_rankings.csv', index=False)
display(signed_rankings.head(10).round(3))

ENRICHMENT_QUESTIONS = {
    'level_shift': 'level_effect',
    'amplitude_change': 'amplitude_log2_ratio',
    'redistribution': 'redistribution',
    'early_contrast': 'early_delta',
    'late_contrast': 'late_delta',
}
N_ENRICHMENT_PERMUTATIONS = 500
enrichment_summary = []
for question, column in ENRICHMENT_QUESTIONS.items():
    statistics = signed_rankings[['gene', column]].dropna().set_index('gene')[column]
    table = camera_like_enrichment(
        statistics,
        module_gene_sets,
        expression=Y_genes,
        gene_names=gene_names,
        background=gene_names,
        n_permutations=N_ENRICHMENT_PERMUTATIONS,
        min_set_size=SECTION4_CONFIG['pathway_min_genes'],
    )
    table.insert(0, 'question', question)
    table.insert(2, 'ranking_column', column)
    table.to_csv(CURVE_OUTPUT_DIR / f'signed_enrichment_{question}.csv', index=False)
    significant = table[table['p_value_permutation_adjusted'] < 0.05]
    enrichment_summary.append({
        'question': question,
        'n_genes_ranked': int(len(statistics)),
        'n_sets_tested': int(table['p_value_permutation'].notna().sum()),
        'n_significant_after_BH': int(len(significant)),
        'median_correlation_inflation': float(table['correlation_inflation'].median()),
    })
    if len(significant):
        display(significant.head(5)[[
            'question', 'gene_set', 'n_genes_tested', 'set_mean_statistic',
            'p_value_permutation', 'p_value_permutation_adjusted',
            'p_value_correlation_aware', 'correlation_inflation',
        ]].round(4))
enrichment_summary = pd.DataFrame(enrichment_summary)
enrichment_summary.to_csv(CURVE_OUTPUT_DIR / 'signed_enrichment_summary.csv', index=False)
display(enrichment_summary)

# Signed downstream-response signatures (PROGENy-style). No signature matrix is bundled: nothing
# biological is invented here, and the step is skipped until one is supplied as
# {"signature": {"GENE": weight}} where a negative weight means the gene represses the signature.
SIGNATURE_PATH = DATA_ROOT / 'mouse_vs_human' / 'signed_signature_matrix.json'
if SIGNATURE_PATH.exists():
    signature_scores, signature_coverage = score_signed_signatures(
        Y_genes, json.loads(SIGNATURE_PATH.read_text()), gene_names=gene_names,
    )
    signature_coverage.to_csv(CURVE_OUTPUT_DIR / 'signed_signature_coverage.csv', index=False)
    usable_signatures = list(signature_coverage.loc[signature_coverage['usable'], 'signature'])
    if usable_signatures:
        signature_fit = run_level_shape(
            signature_scores[usable_signatures].to_numpy(dtype=float),
            s, c, knots, grid, SECTION4_CONFIG['lambda_grid'],
        )
        signature_summary = pd.DataFrame({
            'signature': usable_signatures,
            'level_effect_human_minus_mouse': signature_fit['level_effect'],
            'shape_rms': signature_fit['shape_rms'],
            'curve_spearman': signature_fit['curve_spearman'],
        })
        signature_summary.to_csv(CURVE_OUTPUT_DIR / 'signed_signature_trajectories.csv', index=False)
        display(signature_summary.round(3))
    display(signature_coverage)
else:
    print(f'Signed-signature scoring skipped: {SIGNATURE_PATH.name} is not present. Supply a '
          'PROGENy-style {"signature": {"GENE": weight}} matrix to enable it; pathway-component '
          'expression alone does not establish pathway activity.')


# %% [markdown]
# ## 4.7 - Magnitude checks: detection, abundance bias, within-species gradients, matched regions
#
# Before any top-ranked gene is described as a species difference, four magnitude checks make the
# confounders visible:
#
# 1. **Retained orthologs must actually be assayed in both species.** Detection fractions and mean
#    abundance per species, with the share of genes detected on both sides, written to
#    `between_species_detection_abundance.csv`.
# 2. **Abundance-dependent bias.** Expression ratio against abundance (`ratio_vs_abundance.png`): a
#    trend means the ratio carries a technical component, not only biology.
# 3. **Within-species gradients.** Comparing the early-to-late rise *inside* each species removes a
#    constant gene-specific offset between species, which is what dominates the pooled level term. It
#    cannot remove position-dependent capture differences or the donor confound.
# 4. **Matched-region, per-specimen contrasts.** Each human specimen is compared with each mouse
#    specimen within the same pseudospace bin, so the comparison needs no shared absolute scale and the
#    spread across the four specimen pairs says how reproducible the direction is.
#

# %%
# Purpose: magnitude checks before any species claim.
# 1. Detection and abundance per species: a retained ortholog that is not assayed on one side cannot
#    support a cross-species statement.
species_is_mouse = (species == 'mouse')
detection_abundance = pd.DataFrame({
    'gene': gene_names,
    'detected_fraction_mouse': np.asarray((Y_genes[species_is_mouse] > 0).mean(axis=0)).ravel(),
    'detected_fraction_human': np.asarray((Y_genes[~species_is_mouse] > 0).mean(axis=0)).ravel(),
    'mean_lognorm_mouse': np.asarray(Y_genes[species_is_mouse].mean(axis=0)).ravel(),
    'mean_lognorm_human': np.asarray(Y_genes[~species_is_mouse].mean(axis=0)).ravel(),
})
detection_abundance['log2_ratio_human_over_mouse'] = np.log2(
    (detection_abundance['mean_lognorm_human'] + 1e-3)
    / (detection_abundance['mean_lognorm_mouse'] + 1e-3)
)
detection_abundance['mean_abundance'] = 0.5 * (
    detection_abundance['mean_lognorm_mouse'] + detection_abundance['mean_lognorm_human']
)
detection_abundance['detected_in_both'] = (
    (detection_abundance['detected_fraction_mouse'] > 0)
    & (detection_abundance['detected_fraction_human'] > 0)
)
detection_abundance.to_csv(CURVE_OUTPUT_DIR / 'between_species_detection_abundance.csv', index=False)
print(f'Genes tested: {len(detection_abundance):,}; detected in both species: '
      f'{int(detection_abundance["detected_in_both"].sum()):,} '
      f'({detection_abundance["detected_in_both"].mean():.1%}); median log2 human/mouse ratio '
      f'{detection_abundance["log2_ratio_human_over_mouse"].median():.2f}')

# 2. Abundance-dependent bias: bin by abundance and look at the median ratio per bin, so a systematic
#    trend cannot be mistaken for a gene-specific species difference.
abundance_bins = pd.qcut(detection_abundance['mean_abundance'], 10, duplicates='drop')
ratio_by_abundance = detection_abundance.groupby(abundance_bins, observed=True)[
    'log2_ratio_human_over_mouse'
].agg(['size', 'median', 'mean', 'std']).reset_index()
ratio_by_abundance['abundance_midpoint'] = ratio_by_abundance['mean_abundance'].apply(
    lambda interval: float(interval.mid)
)
ratio_by_abundance.to_csv(CURVE_OUTPUT_DIR / 'ratio_vs_abundance_bins.csv', index=False)
display(ratio_by_abundance.round(3))
fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
axes[0].scatter(detection_abundance['mean_abundance'],
                detection_abundance['log2_ratio_human_over_mouse'], s=3, alpha=0.25,
                color='#4C72B0', linewidths=0)
axes[0].axhline(0.0, color='black', lw=0.8)
axes[0].plot(ratio_by_abundance['abundance_midpoint'], ratio_by_abundance['median'],
             color='crimson', lw=2, marker='o', ms=4, label='median per abundance decile')
axes[0].legend(frameon=False, fontsize=8)
axes[0].set_xlabel('Mean fitted lognorm expression (both species)')
axes[0].set_ylabel('log2(human / mouse)')
axes[0].set_title('Ratio versus abundance')
axes[1].hist(detection_abundance['log2_ratio_human_over_mouse'].dropna(), bins=80,
             color='#4C72B0')
axes[1].axvline(0.0, color='black', lw=0.8)
axes[1].set_xlabel('log2(human / mouse)')
axes[1].set_title('Distribution of mean-expression ratios')
fig.suptitle('Magnitude checks: is the ratio abundance-dependent?')
fig.tight_layout()
fig.savefig(CURVE_OUTPUT_DIR / 'ratio_vs_abundance.png', dpi=180, bbox_inches='tight')
plt.show()

# 3. Within-species early-to-late gradients (constant offset cancels).
early_third = grid <= np.quantile(grid, 1 / 3)
late_third = grid >= np.quantile(grid, 2 / 3)
within_species_gradients = pd.DataFrame({
    'gene': gene_names,
    'mouse_early_to_late': np.nanmean(balanced_reference[:, late_third], axis=1)
        - np.nanmean(balanced_reference[:, early_third], axis=1),
    'human_early_to_late': np.nanmean(balanced_comparison[:, late_third], axis=1)
        - np.nanmean(balanced_comparison[:, early_third], axis=1),
})
within_species_gradients['gradient_difference_human_minus_mouse'] = (
    within_species_gradients['human_early_to_late']
    - within_species_gradients['mouse_early_to_late']
)
within_species_gradients = within_species_gradients.merge(
    gene_results[['gene', 'level_effect_human_minus_mouse', 'difference_type']], on='gene', how='left'
)
within_species_gradients.to_csv(CURVE_OUTPUT_DIR / 'within_species_gradients.csv', index=False)
finite_gradients = within_species_gradients.dropna(subset=['mouse_early_to_late', 'human_early_to_late'])
print('Within-species early-to-late gradients: Spearman(mouse, human) =',
      round(spearmanr(finite_gradients['mouse_early_to_late'],
                      finite_gradients['human_early_to_late']).correlation, 3),
      '; opposite-sign gradients:', int((np.sign(finite_gradients['mouse_early_to_late'])
                                         != np.sign(finite_gradients['human_early_to_late'])).sum()),
      'of', len(finite_gradients))
highest_ratio = detection_abundance.nlargest(20, 'log2_ratio_human_over_mouse')['gene']
display(within_species_gradients.set_index('gene').loc[highest_ratio].round(3))

# 4. Matched-region, per-specimen contrast: compare each human specimen with each mouse specimen
#    inside the same pseudospace bin, on within-row fractions so no shared absolute scale is assumed.
pseudobulk_counts = pseudobulk[[c for c in pseudobulk.columns if c.startswith('count_')]].to_numpy(dtype=float)
row_totals = np.maximum(pseudobulk_counts.sum(axis=1), 1.0)
pseudobulk_fraction = pd.DataFrame(
    pseudobulk_counts / row_totals[:, None],
    columns=[c.replace('count_', '') for c in pseudobulk.columns if c.startswith('count_')],
)
pseudobulk_fraction['specimen'] = pseudobulk['specimen'].to_numpy()
pseudobulk_fraction['pseudospace_bin'] = pseudobulk['pseudospace_bin'].to_numpy()
human_pseudobulk = pseudobulk_fraction['specimen'].isin(HUMAN_SAMPLES)
pair_ratios = []
for bin_index, block in pseudobulk_fraction.groupby('pseudospace_bin', observed=True):
    human_rows = block[block['specimen'].isin(HUMAN_SAMPLES)]
    mouse_rows = block[~block['specimen'].isin(HUMAN_SAMPLES)]
    if human_rows.empty or mouse_rows.empty:
        continue
    genes = [column for column in block.columns if column not in ('specimen', 'pseudospace_bin')]
    human_values = human_rows[genes].to_numpy(dtype=float)
    mouse_values = mouse_rows[genes].to_numpy(dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        ratios = np.log2((human_values[:, None, :] + 1e-6) / (mouse_values[None, :, :] + 1e-6))
    for human_index, human_name in enumerate(human_rows['specimen']):
        for mouse_index, mouse_name in enumerate(mouse_rows['specimen']):
            pair_ratios.append(pd.DataFrame({
                'bin': int(bin_index),
                'human_specimen': human_name,
                'mouse_specimen': mouse_name,
                'gene': genes,
                'log2_ratio': ratios[human_index, mouse_index],
            }))
if pair_ratios:
    matched_region_ratios = pd.concat(pair_ratios, ignore_index=True)
    matched_region_ratios.to_csv(CURVE_OUTPUT_DIR / 'matched_region_specimen_ratios.csv', index=False)
    per_gene = matched_region_ratios.groupby('gene')['log2_ratio'].agg(
        ['median', 'std', 'size']).reset_index()
    # Every gene is measured in the same specimen pairs, so the unit count is a single number.
    per_gene['n_specimen_pairs'] = matched_region_ratios[
        ['human_specimen', 'mouse_specimen']].drop_duplicates().shape[0]
    per_gene.to_csv(CURVE_OUTPUT_DIR / 'matched_region_per_gene.csv', index=False)
    display(per_gene.sort_values('median', ascending=False).head(10).round(3))
    print('Matched-region contrast: bins x specimen pairs =',
          len(matched_region_ratios), 'rows; genes with a consistent direction across bins and pairs:',
          int((per_gene['median'].abs() > 0.25).sum()))
else:
    print('Matched-region contrast skipped: no bin holds both a human and a mouse specimen.')

# 5. Capture summary, so an abundance difference can be read against the measurement itself.
capture_columns = [c for c in ('n_genes_by_counts', 'total_counts', 'n_spots', 'pct_counts_mt',
                              'x_centroid')
                   if c in adata_pt.obs.columns]
capture_summary = adata_pt.obs.groupby(['comparison_species', 'sample', 'region'], observed=True)[
    capture_columns].median().reset_index()
capture_summary.to_csv(CURVE_OUTPUT_DIR / 'capture_summary.csv', index=False)
display(capture_summary)


# %% [markdown]
# # Section 5 - exploratory PT spatial concordance
#
# Distance to reviewed glomerular anchors is a diagnostic proxy, computed separately by sample. Broad marker suggestions alone do not establish a reviewed glomerulus; this comparison is skipped without explicit reviewed anchors.

# %%
# Purpose: adata_pass1 = sc.read_h5ad(HARMONY_OUTPUT_PATH)
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
            & adata_pass1.obs['leiden_coarse'].astype(str).isin([k for k,v in REVIEWED_CLUSTER_LABELS.items() if v == 'Glomerulus']).to_numpy()
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
        if valid_anchor.sum() < 2:
            concordance_rows.append({'sample': sample_name, 'status': 'skipped: insufficient finite reviewed anchors'})
            continue
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
#
#

# %%
# Purpose: analysis_notes = f"""# Human versus healthy-mouse pseudospace
analysis_notes = f"""# Human versus healthy-mouse pseudospace

- Shared matrix before cell typing: {adata_combined.n_obs:,} structures x {adata_combined.n_vars:,} ortholog genes.
- Global nephron DPT: {adata_total.n_obs:,} root-connected tubular-nephron structures; downstream PT subset: {adata_dpt_saved.n_obs:,} structures.
- PT comparison: {adata_pt.n_obs:,} structures over shared support [{lo:.4f}, {hi:.4f}].
- Mouse controls: {', '.join(MOUSE_SAMPLES)}.
- Human samples: {', '.join(HUMAN_SAMPLES)}; both come from one donor.
- Ortholog mapping: HCOP reciprocal best matches with support from at least 3 databases; identical symbols are preferred when support ties. No manual ortholog overrides.
- Integration: deterministic R Harmony {HARMONY_EXPECTED_VERSION}; sample is the batch key and sample correction may remove species/region biology.
- Cell labels: provisional broad pass-1 labels; run fingerprint {CLUSTERING_FINGERPRINT['membership_sha1']}; no relabeling after pass 2.
- Species comparisons: descriptive level, shape, and fitted-curve effect sizes only.
- Human sampling (two healthy cortex slices, named HUK1_COR1/HUK1_MED1 at source) and donor identity are inseparable from species in this cohort.
- The two human slice curves are sensitivity views, not independent biological replicates.
- Ambiguous and non-tubular structures are retained in pass 1 but excluded before pass-2 nephron Harmony; the glomerular cluster is excluded because it is not part of the tubular continuum.
- Global DPT is inspected before PT is selected. `cross_species_nephron_global_dpt.h5ad` carries it as `total_scanpy_dpt`; in `cross_species_pt_dpt.h5ad` that column is the PT-specific coordinate and the global one is kept as `global_nephron_dpt`.
- Reviewed labels are keyed by Leiden cluster ID and are not pinned to `CLUSTERING_FINGERPRINT` in this notebook; re-read `REVIEWED_CLUSTER_LABELS` after any change to the clustering inputs.
- Fine-program diagnostic figure (`global_nephron_dpt_by_fine_reference_program.png`): panel = mean lognormalised expression of the panel genes detected in >= 5% of that species' structures (>= 2 usable genes), z-scored within species, argmax with a 0.15 SD margin. It is a marker check, not a label source, and its counts are not reviewed segment counts; see global_dpt_fine_program_vs_reviewed_segment.csv, global_dpt_fine_program_eligibility.csv and global_dpt_fine_program_species_bias.csv.
- Columns ending cellwise_uncalibrated are diagnostics and must not be interpreted as species tests.
- Peak positions are coordinates on this notebook's own DPT construction (PT-specific, recomputed after global-DPT review and not comparable with the mouse-only workflow's global-nephron coordinate).
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

