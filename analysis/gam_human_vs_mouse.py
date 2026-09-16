# %% [markdown]
# # GAM analysis of the human vs healthy-mouse PT pseudospace
#
# Two analyses, in order:
#
# * **Part 1 - gene analysis.** The nested level/shape GAM for every tested gene, the canonical PT
#   marker-gradient check, and the specimen-balanced weighting.
# * **Part 2 - pathway analysis.** The same GAM on predefined pathway scores, with membership
#   coverage, redundancy grouping and member-level evidence.
#
# Everything above Part 1 is setup and shared input. This notebook **consumes** notebook 03's upstream
# artifacts - the PT-specific DPT coordinate (`cross_species_pt_dpt.h5ad`), the reviewed cluster labels
# carried on it, and the accepted ortholog map (`ortholog_map_used.csv`) - and never rebuilds Harmony,
# the diffusion map, DPT or the clustering. Curve-module discovery, the signed-ranking enrichment
# sweeps and the magnitude/detection audits stay in notebook 03.
#
# Run notebook 03 first: this notebook stops with a clear error if the upstream artifacts are missing.
#

# %% [markdown]
# ## Setup
#
# Configuration, the upstream dependencies, and the coordinate and gene set every
# analysis below shares.

# %%
# Purpose: from __future__ import annotations
# %load_ext autoreload
# %autoreload 2


import argparse
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


# Upstream run: notebook 03 owns the pseudospace, the cluster labels and the ortholog
# map. This notebook reads them and never writes into that directory.
UPSTREAM_DIR = RESULTS_ROOT / 'human_vs_healthy_mouse'
RESULTS_DIR = UPSTREAM_DIR
OUTPUT_DIR = RESULTS_ROOT / 'gam_human_vs_mouse'
# The PT file below is a post-global-DPT subset for later PT analyses.
DPT_OUTPUT_PATH = RESULTS_DIR / 'cross_species_pt_dpt.h5ad'
CURVE_OUTPUT_DIR = OUTPUT_DIR      # every table and figure written here
DIAGNOSTIC_DIR = OUTPUT_DIR / 'diagnostics'
for _directory in (OUTPUT_DIR, DIAGNOSTIC_DIR):
    _directory.mkdir(parents=True, exist_ok=True)

STAGE_CACHE_ENABLED = os.environ.get('PSEUDOSPACE_STAGE_CACHE', '1').strip().lower() not in ('0', 'false', 'no', '')
STAGE_CACHE_DIR = UPSTREAM_DIR / "stage_cache"   # reuse 03's cache: same fits hit
STAGE_CACHE_DIR.mkdir(parents=True, exist_ok=True)
print(f'Stage cache: {"on" if STAGE_CACHE_ENABLED else "off"} ({STAGE_CACHE_DIR})')



print(f'Project: {PROJECT_DIR.name}')
print(f'Reading upstream artifacts from: {UPSTREAM_DIR.relative_to(PROJECT_DIR)}')
print(f'Results will be written to: {OUTPUT_DIR.relative_to(PROJECT_DIR)}')


# %%
# Purpose: # Settings used only for pass-1 loading, QC, and integration.
# The specimen groups the fits compare, and the normalisation target the native-denominator
# sensitivity uses. Pass-1 loading, QC and Harmony happen upstream in notebook 03; this notebook
# only reads their result, so their parameters do not belong here.
MOUSE_SAMPLES = ('Ctrl1A2', 'Ctrl1A4')
HUMAN_SAMPLES = ('HUK1_COR1', 'HUK1_MED1')
SPECIES_GROUPS = ('mouse', 'human')

# This gene filter is applied immediately after structure QC, before HVG selection
# and Harmony, so low-support genes cannot influence the integration embedding.
NORMALIZE_TARGET_SUM = 1e4



# %%
# Purpose: import scanpy as sc
import scanpy as sc
import matplotlib.pyplot as plt
from IPython.display import display
from scipy import sparse
from scipy.stats import spearmanr

from pseudospace.stage_cache import (
    cached_run_level_shape,
    digest,
)

import matplotlib as _mpl
_mpl_cfg = Path(_mpl.get_configdir()).resolve()
if not _mpl_cfg.is_relative_to(CACHE_ROOT.resolve()):
    warnings.warn(
        f'matplotlib is caching to {_mpl_cfg}, not under {CACHE_ROOT.resolve()}. '
        'Run the configuration cell before importing Scanpy/matplotlib.', stacklevel=2)
else:
    print(f'Cache redirects active ({CACHE_ROOT}).')



# %% [markdown]
# ### Ortholog map (upstream)
#
# The accepted human-to-mouse map notebook 03 wrote. Pathway membership resolves through it, and
# the analysis gene set is gated on which features each input actually measured
# (`measured_in_both_inputs`).

# %%
# Read the map 03 produced rather than re-deriving it: deriving it needs the HCOP table plus the
# curated overrides, and 03 already wrote the exact mapping this analysis must use.
ORTHOLOG_MAP_PATH = RESULTS_DIR / 'ortholog_map_used.csv'
if not ORTHOLOG_MAP_PATH.exists():
    raise FileNotFoundError(
        f'{ORTHOLOG_MAP_PATH} is missing. Run analysis/notebooks/03_human_vs_healthy_mouse.ipynb '
        'first: this notebook consumes its pseudospace coordinate and ortholog map.'
    )
ortholog_map = pd.read_csv(ORTHOLOG_MAP_PATH)
print(f'Ortholog map: {ORTHOLOG_MAP_PATH.relative_to(PROJECT_DIR)} ({len(ortholog_map):,} pairs, '
      f'statuses {ortholog_map["mapping_status"].value_counts().to_dict()})')

# %%
# The retained tubular-nephron classes come from the shared package definition, so this
# notebook cannot drift from the workflows that produced its input.
from pseudospace.vocabulary import KEEP_TUBULE_CLASSES
print('Retained tubular-nephron classes:', ', '.join(KEEP_TUBULE_CLASSES))


# %% [markdown]
# ### Coordinate and analysis gene set
#
# The PT pseudospace coordinate and the gene set every fit below uses. The coordinate comes
# straight from 03; the gene set is the detectable genes that **both** inputs actually measured,
# so a structural zero from a missing feature can never enter a comparison.

# %%
# Purpose: from pseudospace.levelshape import fit_single_condition_curves, run_level_shape, summarize_curve_effects
from pseudospace.levelshape import (
    fit_single_condition_curves,
    summarize_curve_effects,
)
from pseudospace.pathways import (
    build_pathway_membership,
    member_gene_evidence,
    summarize_pathway_redundancy,
)
from pseudospace.stats_gam import (
    as_csr,
    gam_internal_knots,
    resolve_present,
)

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
# `measured_in_both_inputs` is the availability rule: a gene absent from one input's feature list is
# a structural zero there, so it cannot enter any comparison (ranking, curve fit, pathway membership
# or enrichment background). Everything downstream of `tested`/`gene_names` is gated by this.
measured_in_both = adata_pt.var['measured_in_both_inputs'].to_numpy(dtype=bool)
tested = (
    (detected >= min_detected)
    & (gene_mean_all >= SECTION4_CONFIG['min_mean_expression'])
    & measured_in_both
)
print(f'Analysis gene set: {int(tested.sum()):,} of {adata_pt.n_vars:,} genes are detectable and '
      f'measured in both inputs ({int((~measured_in_both).sum()):,} excluded as unmeasured).')
gene_names = adata_pt.var_names.to_numpy()[tested]
Y_genes = Y_all[:, tested].tocsr().astype(np.float64)
gene_lookup_tested = {str(gene).upper(): index for index, gene in enumerate(gene_names)}
# Fingerprint of the tested-gene matrix, reused by every cached stage below so it is hashed once
# per run rather than once per call site.
Y_GENES_FINGERPRINT = digest(Y_genes)

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
# ### Reference markers
#
# Axis-adjacent genes. They flag the early/late programs in the gene tables and orient the read
# of the fitted curves; nothing here is fitted or selected on them, and because the axis was built
# from the same markers upstream they are not independent validation.

# %%
# Purpose: define reference genes used only after manual PT cluster review.
# These lists do not participate in clustering, differential expression, or broad labels.
NEPHRON_AXIS_MARKERS = {
    'PT-S1': ['Slc5a2', 'Slc5a12', 'Gatm'],
    'PT-S2': ['Slc22a6', 'Slc13a3', 'Cyp2e1'],
    'PT-S3': ['Slc7a13', 'Slc22a7', 'Cyp7b1'],
}
print('Reference-only PT groups:', ', '.join(NEPHRON_AXIS_MARKERS))

# Purpose: # Global-DPT orientation settings: defined where the trajectory is constructed.
# Global-DPT orientation settings: defined where the trajectory is constructed.
TOTAL_POSITION_MARKERS = {
    'early': sorted({gene for group in ('PT-S1', 'PT-S2') for gene in NEPHRON_AXIS_MARKERS[group]}),
    'late': sorted(NEPHRON_AXIS_MARKERS['PT-S3']),
}
print('PT reference axis for global-DPT orientation:', {key: len(value) for key, value in TOTAL_POSITION_MARKERS.items()})


# %% [markdown]
# ## Part 1 - gene analysis

# %% [markdown]
# ### 1.1 - gene-level level/shape decomposition
#
# M0 is one shared smooth, M1 adds a constant species offset and M2 adds a species-by-pseudospace
# interaction; the likelihood-ratio tests choose between them, and the effect split then separates a
# vertical offset (level), a gradient-strength change (amplitude) and a spatial redistribution
# (pattern). The comparison is repeated with the native total as denominator, so a result cannot be
# attributed to the normalisation choice alone.
#

# %%
# Purpose: gene_fit = run_level_shape(
gene_fit = cached_run_level_shape(
    Y_genes, s, c, knots, grid, SECTION4_CONFIG['lambda_grid'],
    stage='section4_gene_fit', root=STAGE_CACHE_DIR, y_fingerprint=Y_GENES_FINGERPRINT,
    enabled=STAGE_CACHE_ENABLED,
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
# Only the columns the fit table does not already carry are merged: `shape_rms`, `level_effect` and
# `condition_effect_rms` come from the fit itself, and a straight merge would suffix them
# (`shape_rms_x`/`shape_rms_y`) and break every table that selects them.
gene_curve_effects = gene_curve_effects.rename(columns={'feature': 'gene'})
gene_results = gene_results.merge(
    gene_curve_effects[['gene', 'level_fraction', 'shape_fraction', 'pattern_rms_z', 'amplitude_reference', 'amplitude_comparison', 'amplitude_ratio', 'amplitude_log2_ratio', 'pattern_status', 'difference_type']],
    on='gene', how='left',
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
    native_fit = cached_run_level_shape(
        Y_native, s, c, knots, grid, SECTION4_CONFIG['lambda_grid'],
        stage='section4_gene_fit_native_denominator', root=STAGE_CACHE_DIR,
        y_fingerprint=digest(Y_native), enabled=STAGE_CACHE_ENABLED,
    )
    # `gene_results` is sorted by effect size at this point, while the fit arrays are in gene order:
    # assigning them positionally compared different genes with one another. Merge on the gene name,
    # and carry the alternative fit's own total effect so the level fraction uses the right
    # denominator (the panel fit's `species_effect_rms` belongs to the other normalisation).
    native_effects = pd.DataFrame({
        'gene': gene_names,
        'level_effect_native_denominator': native_fit['level_effect'],
        'shape_rms_native_denominator': native_fit['shape_rms'],
        'condition_effect_rms_native_denominator': native_fit['condition_effect_rms'],
    })
    gene_results = gene_results.merge(native_effects, on='gene', how='left')
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
            / np.maximum(top_native['condition_effect_rms_native_denominator'] ** 2, 1e-12)
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
    'shape_rms', 'level_fraction', 'pattern_rms_z', 'pattern_status', 'difference_type',
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
# ### 1.2 - canonical PT marker-gradient check
#
# The reference markers plotted against the fitted mouse and human curves: they show where each PT
# program sits along the shared coordinate and whether that placement is consistent between species.
#

# %%
# Reference marker panels for the check below. They do not participate in any fit: they are
# plotted against the fitted curves to show where each PT program sits on the axis.
# Purpose: # Marker panels and display settings used only by the heatmap section.
# Marker panels and display settings used only by the heatmap section.
HEATMAP_MARKERS = {
    'PT': {
        'PT-S1': ['Slc5a2', 'Slc5a12', 'Gatm'],
        'PT-S2': ['Slc22a6', 'Slc13a3', 'Cyp2e1'],
        'PT-S3': ['Slc22a7', 'Cyp7b1'],
    },
}
def _marker_dict(markers):
    return {group: [{'label': gene, 'candidates': [gene]} for gene in genes] for group, genes in markers.items()}
print(f'Heatmap families: {", ".join(HEATMAP_MARKERS)}')

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
# ### 1.3 - specimen-balanced weighting
#
# The pooled fits weight tubules, so a specimen that contributed more structures shapes every curve
# more. Fitting each specimen separately and averaging the curves gives each specimen equal weight -
# the version that cannot be driven by one specimen's share of the data. Both weightings are reported.
#

# %%
# Purpose: per-specimen GAM fits and the specimen-balanced mixture.
# The pooled fit above weights TUBULES, so a specimen that contributed more structures shapes every
# curve more. Fitting each specimen separately and averaging the curves gives each specimen equal
# weight; comparing the two weightings says how much of a result is a weighting artefact.
from pseudospace.specimen import specimen_balanced_curves

specimen_full_curves = {}
for sample_name in sorted(set(samples)):
    mask = samples == sample_name
    specimen_full_curves[sample_name], _ = fit_single_condition_curves(
        Y_genes[mask], s[mask], knots, grid, SECTION4_CONFIG['lambda_grid'],
        gene_fit['lam_idx'], support_pct=SECTION4_CONFIG['common_support_pct'],
    )

mouse_specimens = [name for name in specimen_full_curves if name in MOUSE_SAMPLES]
human_specimens = [name for name in specimen_full_curves if name in HUMAN_SAMPLES]
balanced_reference = specimen_balanced_curves({n: specimen_full_curves[n] for n in mouse_specimens})
balanced_comparison = specimen_balanced_curves({n: specimen_full_curves[n] for n in human_specimens})

specimen_balanced_frame = pd.DataFrame({
    'gene': gene_names,
    'balanced_mouse_amplitude': (np.nanmax(balanced_reference, axis=1)
                                 - np.nanmin(balanced_reference, axis=1)),
    'balanced_human_amplitude': (np.nanmax(balanced_comparison, axis=1)
                                 - np.nanmin(balanced_comparison, axis=1)),
    'balanced_level_effect_human_minus_mouse': np.nanmean(
        balanced_comparison - balanced_reference, axis=1),
})
specimen_balanced_frame.to_csv(CURVE_OUTPUT_DIR / 'specimen_balanced_gene_amplitudes.csv',
                               index=False)
print(f'Specimen fits: {len(specimen_full_curves)} specimens ({len(mouse_specimens)} mouse, '
      f'{len(human_specimens)} human); balanced curve matrix {balanced_reference.shape}.')

# The same estimator on both weightings: a result that survives this is not an artefact of how many
# tubules each specimen contributed.
pooled_effects = summarize_curve_effects(
    gene_fit['curve_healthy'], gene_fit['curve_aki'], feature_names=gene_names
).rename(columns={'feature': 'gene'})
balanced_effects = summarize_curve_effects(
    balanced_reference, balanced_comparison, feature_names=gene_names
).rename(columns={'feature': 'gene'})
weighting = pooled_effects.merge(balanced_effects, on='gene', suffixes=('_pooled', '_balanced'))
weighting_agreement = pd.DataFrame([
    {
        'metric': metric,
        'spearman_pooled_vs_balanced': float(weighting[metric + '_pooled'].corr(
            weighting[metric + '_balanced'], method='spearman')),
        'top20_sign_agreement': float(np.mean(
            np.sign(weighting.reindex(weighting[metric + '_pooled'].abs().nlargest(20).index)[
                metric + '_pooled'])
            == np.sign(weighting.reindex(weighting[metric + '_pooled'].abs().nlargest(20).index)[
                metric + '_balanced']))),
    }
    for metric in ('level_effect', 'amplitude_log2_ratio', 'pattern_rms_z', 'level_fraction')
])
weighting_agreement.to_csv(CURVE_OUTPUT_DIR / 'pooled_vs_specimen_balanced_agreement.csv',
                           index=False)
display(weighting_agreement.round(3))
print('Where the two weightings disagree, the pooled number is the one that depends on how many '
      'tubules each specimen contributed.')

fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
axes[0].scatter(weighting['level_effect_pooled'], weighting['level_effect_balanced'], s=4,
                alpha=0.25, color='#4C9BD3')
axes[0].set_xlabel('level effect, pooled fit')
axes[0].set_ylabel('level effect, balanced fit')
axes[0].set_title(f"level: Spearman {weighting_agreement.loc[0, 'spearman_pooled_vs_balanced']:+.3f}")
axes[1].scatter(weighting['pattern_rms_z_pooled'], weighting['pattern_rms_z_balanced'], s=4,
                alpha=0.25, color='#E15759')
axes[1].set_xlabel('pattern RMS (z), pooled')
axes[1].set_ylabel('pattern RMS (z), balanced')
axes[1].set_title(f"pattern: Spearman {weighting_agreement.loc[2, 'spearman_pooled_vs_balanced']:+.3f}")
fig.tight_layout()
fig.savefig(CURVE_OUTPUT_DIR / 'pooled_vs_specimen_balanced.png', dpi=180)
plt.show()


# %% [markdown]
# ## Part 2 - pathway analysis

# %% [markdown]
# ### 2.1 - predefined pathway trajectories
#
# Pathway scores are the mean of PT-standardised member genes, fitted with the same nested GAM as the
# genes. Membership is resolved through the ortholog map with explicit coverage stages, and redundant
# sets are grouped rather than reported as independent hits.
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
    pathway_fit = None      # bound on the empty path so later cells cannot NameError
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

    pathway_fit = cached_run_level_shape(
        pathway_module_scores, s, c, knots, grid, SECTION4_CONFIG['lambda_grid'],
        stage='section4_pathway_fit', root=STAGE_CACHE_DIR,
        y_fingerprint=digest(pathway_module_scores), enabled=STAGE_CACHE_ENABLED,
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
                   'amplitude_log2_ratio', 'pattern_status', 'difference_type'):
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
    'pattern_status', 'difference_type',
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


# %% [markdown]
# ### 2.2 - member evidence
#
# Each prioritized pathway ships member-level evidence: how many of its members move with the
# aggregate trend, which gene contributes most, and whether the direction survives without it.
# A pathway name is not evidence of a distinct program when its member set is another
# pathway's member set.

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


# %%
# Purpose: close with what was written and what must not be read into it.
written = sorted(p.name for p in CURVE_OUTPUT_DIR.glob('*') if p.is_file())
print('Tables and figures written to', CURVE_OUTPUT_DIR.relative_to(PROJECT_DIR) + ':')
for name in written:
    print('  ', name)
print()
print(f'PT structures analysed: {adata_pt.n_obs:,}; tested genes: {len(gene_names):,} '
      f'({int((~measured_in_both).sum()):,} genes excluded as not measured in both inputs).')
print()
print('Reading rules for these numbers:')
print('  - the coordinate is upstream (notebook 03). Nothing here re-derives Harmony, the diffusion')
print('    map, DPT or the cluster labels, so a coordinate change has to be made in 03.')
print('  - genes absent from either input feature list are excluded from every fit (Part 1 and')
print('    Part 2 both): a structural zero is not a measurement (measured_in_both_inputs).')
print('  - the effect split separates level, amplitude and pattern, and a pattern verdict is only')
print('    offered when both curves carry a gradient (see pattern_status).')
print('  - the inference unit is the specimen: two mouse specimens versus one human donor, so these')
print('    are descriptive summaries, not species-level tests.')
print('  - axis-adjacent marker genes check orientation; they do not validate it independently.')

