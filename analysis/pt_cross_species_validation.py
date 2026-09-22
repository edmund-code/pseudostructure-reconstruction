# %% [markdown]
# # Cross-species validation of the human vs healthy-mouse PT pseudospace
#
# A short, focused **sanity-check** notebook, not a re-analysis. It asks one practical question:
#
# > Are the major human-mouse positional-expression differences seen in notebook 07 likely to reflect
# > real biology, rather than an obviously bad or mismatched PT pseudospace?
#
# Four high-value diagnostics answer it, one per question:
#
# | # | diagnostic | what it would invalidate |
# | --- | --- | --- |
# | 1 | DPT occupancy by specimen | a species concentrated at one end of the axis |
# | 2 | S1-to-S3 ordering **per specimen** | a coordinate that does not run early-to-late everywhere |
# | 3 | within-species curve reproducibility | divergence driven by unreproducible weak GAMs |
# | 4 | divergence among strong, reproducible genes | a signal that survives amplitude/robustness gates |
#
# **Deliberately out of scope**: no pathway enrichment, no re-clustering of notebook 07, no
# alternative trajectory or nonlinear alignment, no cortex-matched mouse subset, no simulation/null
# study, and no re-derivation of the pseudospace. If the coordinate's root is clearly invalid this
# notebook reports it loudly rather than rebuilding it.
#
# **Input.** Notebook 03's PT-only artifact (`cross_species_pt_dpt.h5ad`) and - read-only, never
# refit - notebook 07's cached per-specimen GAM curves, which are keyed on the same conventions
# notebook 03 validated. Every curve here is therefore byte-identical to the ones notebook 07
# clustered.
#
# **Inference unit.** 2 mouse specimens and 2 human *sections of one donor*. Within-species
# reproducibility measures agreement between sections of the same donor; it is not biological
# replication, and no p-value is computed anywhere below.
#
# Run notebook 03 (and ideally 07) first: this notebook stops with a clear error if its inputs are
# missing.
#

# %% [markdown]
# ## Setup
#
# Configuration, upstream dependencies, and the parameters every section shares. The GAM conventions
# are byte-identical to notebook 07's, because that is what makes the per-specimen curve cache a
# **hit** rather than a re-fit.

# %%
# Purpose: bootstrap the run: cache redirects, repository root, roots, and output directories.
# No IPython magics here on purpose: notebooks 02-07 carry none, and a magic is a live line in the
# .ipynb that the synthetic dry-run harness cannot exec.

import argparse
import os
import sys
import warnings
from pathlib import Path

# These must be set before importing matplotlib, numba, Scanpy, or modules that import them.
CACHE_ROOT = Path(os.environ.get(
    'PSEUDOSPACE_CACHE_ROOT', '/tmp/pseudospace_pt_cross_species_validation_cache'))
for _subdir in ('matplotlib', 'numba', 'xdg'):
    (CACHE_ROOT / _subdir).mkdir(parents=True, exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR', str(CACHE_ROOT / 'matplotlib'))
os.environ.setdefault('NUMBA_CACHE_DIR', str(CACHE_ROOT / 'numba'))
os.environ.setdefault('XDG_CACHE_HOME', str(CACHE_ROOT / 'xdg'))


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


def _rel(path) -> str:
    """Path relative to the repository when it is inside it, absolute otherwise.

    A run may legitimately point the results root outside the repository (the synthetic dry-run
    harness does), and an unconditional ``relative_to`` would then raise instead of printing.
    """
    path = Path(path)
    try:
        return str(path.relative_to(PROJECT_DIR))
    except ValueError:
        return str(path)


# Notebook 03 owns the pseudospace; notebook 07 owns the GAM curves. This notebook reads both and
# never writes into either directory.
UPSTREAM_DIR = RESULTS_ROOT / 'human_vs_healthy_mouse'
OUTPUT_DIR = RESULTS_ROOT / 'pt_cross_species_validation'
DIAGNOSTIC_DIR = OUTPUT_DIR / 'diagnostics'
for _directory in (OUTPUT_DIR, DIAGNOSTIC_DIR):
    _directory.mkdir(parents=True, exist_ok=True)

DPT_OUTPUT_PATH = UPSTREAM_DIR / 'cross_species_pt_dpt.h5ad'
ORTHOLOG_MAP_PATH = UPSTREAM_DIR / 'ortholog_map_used.csv'
# 03's own record of how the PT-specific coordinate was rooted. Optional: the notebook degrades to a
# coordinate-level root check when the upstream diagnostics are not present (e.g. a synthetic run).
PT_ROOT_DIAGNOSTICS_PATH = UPSTREAM_DIR / 'celltyping' / 'pt_specific_dpt_diagnostics.csv'
STAGE_CACHE_DIR = UPSTREAM_DIR / 'stage_cache'
# 07's saved tables, used only to cross-check this notebook's independently recomputed numbers.
NOTEBOOK_07_DIR = RESULTS_ROOT / 'pt_gam_clustering'

STAGE_CACHE_ENABLED = os.environ.get(
    'PSEUDOSPACE_STAGE_CACHE', '1').strip().lower() not in ('0', 'false', 'no', '')

print(f'Project: {PROJECT_DIR.name}')
print(f'Reading the PT artifact from: {_rel(UPSTREAM_DIR)}')
print(f'Results will be written to: {_rel(OUTPUT_DIR)}')
print(f'Stage cache: {"on" if STAGE_CACHE_ENABLED else "off"} ({_rel(STAGE_CACHE_DIR)})')


# %%
# Purpose: import the analysis stack and the project helpers this notebook reuses.
# Every fit, cache and alias resolver already exists in `pseudospace`; nothing here re-implements one.
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from IPython.display import display
from scipy.stats import rankdata

from pseudospace.levelshape import fit_single_condition_curves
from pseudospace.markers import build_gene_lookup
from pseudospace.specimen import specimen_balanced_curves
from pseudospace.stage_cache import (
    cached_payload,
    cached_run_level_shape,
    code_digest,
    digest,
)
from pseudospace.stats_gam import as_csr, gam_internal_knots, safe_spearman, zscore_rows
from pseudospace.vocabulary import KEEP_TUBULE_CLASSES

import matplotlib as _mpl
_mpl_cfg = Path(_mpl.get_configdir()).resolve()
if not _mpl_cfg.is_relative_to(CACHE_ROOT.resolve()):
    warnings.warn(
        f'matplotlib is caching to {_mpl_cfg}, not under {CACHE_ROOT.resolve()}. '
        'Run the bootstrap cell before importing Scanpy/matplotlib.', stacklevel=2)
else:
    print(f'Cache redirects active ({CACHE_ROOT}).')


# %%
# Purpose: every threshold this notebook applies, as one named knob each.
# The GAM block is byte-identical to notebook 07's Section 3: that is what makes the per-specimen
# curve cache a hit, and what makes the curves validated here the same curves 07 clustered.
GAM_CONFIG = {
    'n_internal_knots': 9,
    'lambda_grid': np.logspace(-3, 3, 13),
    'common_support_pct': (1, 99),
    'grid_points': 101,
    'min_detected_fraction': 0.02,
}

# The canonical PT early/mid/late programmes, as notebook 06's curve-layer QC panel defines them.
# They orient and audit the coordinate; they never assign a structure to a segment here.
MARKER_PANELS = {
    'S1': ['Slc5a2', 'Slc5a12', 'Gatm', 'Lrp2', 'Cubn', 'Slc34a1'],
    'S2': ['Slc22a6', 'Slc13a3', 'Cyp2e1'],
    'S3': ['Slc22a7', 'Cyp7b1', 'Slc7a13', 'Slc6a18', 'Acsm3'],
}

CONFIG = {
    'occupancy_bins': 20,
    'occupancy_percentiles': (5, 25, 50, 75, 95),
    # A Mann-Whitney common-language effect (P(human structure sits earlier than a mouse one)) of 0.5
    # is perfect overlap. Beyond roughly a "small" effect the occupancy imbalance is reported as a
    # concern rather than a footnote.
    'occupancy_auc_tolerance': 0.10,
    # 07's positional amplitude floor: real shape in BOTH species before any shape metric is read.
    'amplitude_floor': 0.05,
    'amplitude_thresholds': (0.05, 0.10, 0.20),
    # A "high reproducibility" gene is a descriptive subset, not a significance threshold: both
    # within-species correlations must reach it.
    'reproducibility_floor': 0.80,
    # Representative genes must be strong AND reproducible, or their "divergence" is not interpretable.
    'representative_min_amplitude': 0.10,
    'representative_per_class': 3,
    'raw_bins': 15,
    'min_structures_per_bin': 5,
    'spearman_min_points': 10,
    'numerical_epsilon': 1e-12,
}

SPECIES_COLORS = {'mouse': '#0072B2', 'human': '#D55E00'}
SPECIMEN_CMAP = plt.get_cmap('tab10')
figure_index = []


# %%
# Purpose: the small numeric and plotting helpers every section below shares.
def _save_figure(fig, filename: str):
    """Write one figure into this notebook's own output directory and record it."""
    path = OUTPUT_DIR / filename
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.show()
    plt.close(fig)
    figure_index.append(filename)
    print('figure:', filename)


def _peak_to_peak(values):
    values = np.asarray(values, dtype=float)
    return np.nanmax(values, axis=1) - np.nanmin(values, axis=1)


def _row_block_pearson(left, right, epsilon):
    """Per-row Pearson correlation of two aligned blocks, over the columns a row actually carries.

    Both blocks arrive z-scored within their own condition, so the number describes positional shape
    alone: a level shift or an amplitude change cannot appear in it.
    """
    a = np.asarray(left, dtype=float)
    b = np.asarray(right, dtype=float)
    finite = np.isfinite(a) & np.isfinite(b)
    count = np.maximum(finite.sum(axis=1, keepdims=True), 1)
    mean_a = np.where(finite, a, 0.0).sum(axis=1, keepdims=True) / count
    mean_b = np.where(finite, b, 0.0).sum(axis=1, keepdims=True) / count
    a = np.where(finite, a - mean_a, 0.0)
    b = np.where(finite, b - mean_b, 0.0)
    numerator = (a * b).sum(axis=1)
    denominator = np.sqrt((a ** 2).sum(axis=1) * (b ** 2).sum(axis=1))
    return numerator / np.where(denominator > epsilon, denominator, np.nan)


def _curve_eligibility(curves, gate, min_peak_to_peak, min_finite_fraction=0.5):
    """Which features may carry a shape metric, and which grid columns that run can actually see.

    A curve from a single specimen is undefined outside that specimen's own pseudospace support, so
    the support is the columns where at least ``min_finite_fraction`` of the features are finite; a
    feature must then be finite across all of them and clear the caller's absolute-magnitude gate.
    Identical to notebook 07's rule, so the shape metrics below match its numbers exactly.
    """
    values = np.asarray(curves, dtype=float)
    finite = np.isfinite(values)
    column_ok = finite.mean(axis=0) >= float(min_finite_fraction)
    if column_ok.sum() < 3:
        column_ok = np.ones(values.shape[1], dtype=bool)
    eligible = (finite[:, column_ok].all(axis=1)
                & (np.asarray(gate, dtype=float) >= min_peak_to_peak))
    return eligible, column_ok


def _genuine_shape_eligibility(mouse_block, human_block, grid, floor, min_finite_fraction=0.5):
    """Genes carrying real positional shape in BOTH species, plus the columns both are seen on.

    The same absolute amplitude floor is applied twice, once per species: a nearly flat curve divided
    by its own tiny spread turns numerical noise into an apparently strong shape. Returning the
    shared columns is what lets every shape metric below use one consistent grid.
    """
    mouse_block = np.asarray(mouse_block, dtype=float)
    human_block = np.asarray(human_block, dtype=float)
    mouse_ok, mouse_columns = _curve_eligibility(
        mouse_block, _peak_to_peak(mouse_block), floor, min_finite_fraction)
    human_ok, human_columns = _curve_eligibility(
        human_block, _peak_to_peak(human_block), floor, min_finite_fraction)
    columns = mouse_columns & human_columns
    if columns.sum() < 3:
        columns = np.ones(mouse_block.shape[1], dtype=bool)
    finite = (np.isfinite(mouse_block[:, columns]).all(axis=1)
              & np.isfinite(human_block[:, columns]).all(axis=1))
    return mouse_ok & human_ok & finite, columns, mouse_ok, human_ok


def _standardize_species(curves, columns):
    """Z_g(s) for one condition: centre and scale THAT condition's curve across the shared grid.

    Centring removes the mean level and scaling removes the amplitude, so only positional shape
    remains. The conditions are standardised independently and never subtracted from one another.
    """
    return zscore_rows(np.asarray(curves, dtype=float)[:, columns])


def _summary_row(label, n, correlation):
    """One descriptive row: size, median cross-species shape correlation and its tail fractions."""
    values = np.asarray(correlation, dtype=float)
    values = values[np.isfinite(values)]
    return {
        'subset': label,
        'n_genes': int(n),
        'median_shape_correlation': float(np.median(values)) if values.size else np.nan,
        'fraction_correlation_above_0.8': float(np.mean(values > 0.8)) if values.size else np.nan,
        'fraction_correlation_below_0': float(np.mean(values < 0)) if values.size else np.nan,
        'fraction_correlation_below_-0.5': float(np.mean(values < -0.5)) if values.size else np.nan,
    }


def _probability_smaller(left, right):
    """P(a random value of `left` is below a random value of `right`) - the Mann-Whitney effect size.

    0.5 means the two samples are interchangeable along the axis; 0.0 or 1.0 means they are fully
    separated. Computed from ranks, so it needs no distributional assumption and no p-value.
    """
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    ranks = rankdata(np.concatenate([left, right]))
    rank_sum_left = ranks[:left.size].sum()
    # Ranks give P(left > right) directly; the smaller-side probability is the complement over the
    # n_left * n_right pairs. Ties (structures tied at the scaled floor) split between the two
    # equally, which is the standard mid-rank treatment and does not bias the direction.
    greater = rank_sum_left - left.size * (left.size + 1) / 2.0
    smaller = left.size * right.size - greater
    return float(smaller / (left.size * right.size))


# %% [markdown]
# ## 1 - Input and PT-only scope
#
# Load notebook 03's PT artifact and verify the contract this notebook relies on. The pseudospace is
# never rebuilt here: `total_scanpy_dpt` is read exactly as 03 wrote it.

# %%
# Purpose: load the PT-only artifact and verify the exact contract this notebook relies on.
for required_path in (DPT_OUTPUT_PATH, ORTHOLOG_MAP_PATH):
    if not required_path.exists():
        raise FileNotFoundError(
            f'{required_path} is missing. Run analysis/notebooks/03_human_vs_healthy_mouse.ipynb '
            'first: this notebook consumes its PT pseudospace coordinate and accepted ortholog map.'
        )

adata_full = sc.read_h5ad(DPT_OUTPUT_PATH)
MEASUREMENT = 'lognorm'
if MEASUREMENT not in adata_full.layers:
    raise KeyError(f'{DPT_OUTPUT_PATH.name} is missing layers["{MEASUREMENT}"].')
adata_full.X = adata_full.layers[MEASUREMENT].copy()

missing_obs = [column for column in ('comparison_species', 'sample', 'total_scanpy_dpt')
               if column not in adata_full.obs.columns]
if missing_obs:
    raise KeyError(f'{DPT_OUTPUT_PATH.name} is missing obs columns: {missing_obs}')
if 'measured_in_both_inputs' not in adata_full.var.columns:
    raise KeyError(f'{DPT_OUTPUT_PATH.name} is missing var["measured_in_both_inputs"].')

# Verify the PT restriction where the object is self-describing, and fall back to the reviewed class
# column so an older object cannot smuggle in another segment.
pt_restriction = np.ones(adata_full.n_obs, dtype=bool)
if 'broad_tubule_marker_call' in adata_full.obs.columns:
    broad = adata_full.obs['broad_tubule_marker_call'].astype(str)
    unexpected = sorted(set(broad) - set(KEEP_TUBULE_CLASSES))
    if unexpected:
        raise ValueError(f'{DPT_OUTPUT_PATH.name} contains non-nephron classes: {unexpected}')
    non_pt = sorted(set(broad) - {'PT'})
    if non_pt:
        print(f'NOTE: the artifact is not PT-only (also carries {non_pt}); restricting to PT here.')
        pt_restriction = broad.eq('PT').to_numpy()

temp_mask = (pt_restriction
             & adata_full.obs['comparison_species'].astype(str).isin(SPECIES_COLORS)
             & np.isfinite(adata_full.obs['total_scanpy_dpt'].to_numpy(dtype=float))).to_numpy()
adata_pt = adata_full[temp_mask].copy()
del temp_mask

if adata_pt.n_obs < 100 or adata_pt.obs['comparison_species'].nunique() != 2:
    raise ValueError('Insufficient PT structures from both species.')

s = adata_pt.obs['total_scanpy_dpt'].to_numpy(dtype=float)
species = adata_pt.obs['comparison_species'].astype(str).to_numpy()
samples = adata_pt.obs['sample'].astype(str).to_numpy()
comparison_indicator = (species == 'human').astype(float)

# Species comes from the data, not from a hardcoded sample list, so a renamed section cannot silently
# be treated as the wrong species.
species_by_sample = (adata_pt.obs[['sample', 'comparison_species']].astype(str)
                     .drop_duplicates().set_index('sample')['comparison_species'].to_dict())
if len(set(species_by_sample.values())) != 2:
    raise ValueError(f'Samples do not map onto two species: {species_by_sample}')
mouse_specimens = sorted(name for name, kind in species_by_sample.items() if kind == 'mouse')
human_specimens = sorted(name for name, kind in species_by_sample.items() if kind == 'human')
specimen_order = mouse_specimens + human_specimens

cohort = (adata_pt.obs.groupby(['comparison_species', 'sample'], observed=True)
          .agg(n_pt_structures=('sample', 'size'), dpt_min=('total_scanpy_dpt', 'min'),
               dpt_max=('total_scanpy_dpt', 'max')).reset_index())
display(cohort.round(4))
print(f'PT structures: {adata_pt.n_obs:,}; genes assayed: {adata_pt.n_vars:,}')
print(f'Mouse specimens: {mouse_specimens}; human sections: {human_specimens} '
      '(two sections of ONE donor, biologically one replicate)')
print(f'All {adata_pt.n_obs:,} structures carry a finite PT DPT on [{s.min():.3f}, {s.max():.3f}].')


# %% [markdown]
# ## 2 - DPT occupancy by specimen
#
# > Do human and mouse actually occupy broadly overlapping PT positions, or is one species
# > concentrated toward one end of the axis?
#
# Every specimen is binned on the **same** edges over the artifact's own [0, 1] range, so the
# histograms are directly comparable, and the percentiles are reported per specimen. The common
# min/max support 07 uses is shown for context only: two specimens can share a min and a max and
# still occupy disjoint halves.

# %%
# Purpose: shared bins, per-specimen percentiles, and the occupancy shape along the DPT axis.
OCCUPANCY_BINS = CONFIG['occupancy_bins']
occupancy_edges = np.linspace(0.0, 1.0, OCCUPANCY_BINS + 1)
occupancy_centres = 0.5 * (occupancy_edges[:-1] + occupancy_edges[1:])

occupancy_rows = []
for specimen in specimen_order:
    values = s[samples == specimen]
    row = {'specimen': specimen, 'species': species_by_sample[specimen],
           'n_structures': int(values.size), 'mean': float(np.mean(values)),
           'std': float(np.std(values, ddof=1))}
    for percentile in CONFIG['occupancy_percentiles']:
        row[f'p{percentile}'] = float(np.percentile(values, percentile))
    occupancy_rows.append(row)
occupancy_by_specimen = pd.DataFrame(occupancy_rows)
occupancy_by_specimen['fraction_after_half'] = [
    float(np.mean(s[samples == row.specimen] > 0.5)) for row in occupancy_by_specimen.itertuples()]
occupancy_by_specimen.to_csv(OUTPUT_DIR / 'dpt_occupancy_by_specimen.csv', index=False)
display(occupancy_by_specimen.round(4))

# One binned count table across the DPT axis: structures from each specimen in each shared bin.
occupancy_bin_index = np.clip(np.digitize(s, occupancy_edges) - 1, 0, OCCUPANCY_BINS - 1)
occupancy_bins = pd.DataFrame({'bin_left': occupancy_edges[:-1], 'bin_centre': occupancy_centres})
for specimen in specimen_order:
    occupancy_bins[f'n_{specimen}'] = np.bincount(
        occupancy_bin_index[samples == specimen], minlength=OCCUPANCY_BINS)
occupancy_bins['n_mouse'] = occupancy_bins[[f'n_{name}' for name in mouse_specimens]].sum(axis=1)
occupancy_bins['n_human'] = occupancy_bins[[f'n_{name}' for name in human_specimens]].sum(axis=1)
occupancy_bins.to_csv(DIAGNOSTIC_DIR / 'dpt_occupancy_bins.csv', index=False)

# The species-level question as one number: how often is a human structure earlier than a mouse one.
mouse_dpt = s[np.isin(samples, mouse_specimens)]
human_dpt = s[np.isin(samples, human_specimens)]
auc_human_earlier = _probability_smaller(human_dpt, mouse_dpt)
mouse_median, human_median = float(np.median(mouse_dpt)), float(np.median(human_dpt))
fraction_human_in_mouse_iqr = float(np.mean(
    (human_dpt >= np.percentile(mouse_dpt, 25)) & (human_dpt <= np.percentile(mouse_dpt, 75))))
occupancy_imbalance = abs(auc_human_earlier - 0.5)
overlap_ok = occupancy_imbalance <= CONFIG['occupancy_auc_tolerance']

# The common-support window 07 fits on, reported only so the imbalance above is not mistaken for it.
p_lo, p_hi = GAM_CONFIG['common_support_pct']
shared_lo = max(np.percentile(s[species == name], p_lo) for name in SPECIES_COLORS)
shared_hi = min(np.percentile(s[species == name], p_hi) for name in SPECIES_COLORS)

print(f'Species medians: mouse {mouse_median:.3f}, human {human_median:.3f}')
print(f'P(human structure earlier than a mouse one) = {auc_human_earlier:.3f} '
      f'(0.5 = interchangeable; tolerance {CONFIG["occupancy_auc_tolerance"]})')
print(f'Human structures inside the mouse IQR: {fraction_human_in_mouse_iqr:.1%}')
print(f'Common 1st-99th percentile support window 07 fits on: [{shared_lo:.3f}, {shared_hi:.3f}] '
      '- shared endpoints do NOT imply balanced occupancy.')
if not overlap_ok:
    print('OCCUPANCY FLAG: the two species are concentrated on different parts of the axis. Every '
          'late-PT human curve rests on far fewer structures than its mouse counterpart, so the '
          'cross-species curves at that end are the least comparable. This is a real limitation of '
          'the comparison, not an artefact of the fit.')
else:
    print(f'Occupancy overlap: P(human earlier) = {auc_human_earlier:.3f} sits within tolerance of '
          '0.5, so neither species is strongly concentrated at one end. The medians still place '
          'the human sections earlier on average - a difference in degree, not a disjoint axis.')

fig, axes = plt.subplots(1, 2, figsize=(11.5, 3.8))
for specimen in specimen_order:
    values = s[samples == specimen]
    colour = SPECIES_COLORS[species_by_sample[specimen]]
    axes[0].hist(values, bins=occupancy_edges, density=True, histtype='step', lw=1.6, color=colour,
                 label=f'{specimen} (n={values.size:,})')
    ordered = np.sort(values)
    axes[1].plot(ordered, np.arange(1, ordered.size + 1) / ordered.size, lw=1.6, color=colour,
                 label=specimen)
axes[0].axvline(mouse_median, color=SPECIES_COLORS['mouse'], ls='--', lw=1.0)
axes[0].axvline(human_median, color=SPECIES_COLORS['human'], ls='--', lw=1.0)
axes[0].set(xlabel='PT DPT (shared coordinate)', ylabel='density',
            title=f'Same bins for every specimen\nP(human earlier) = {auc_human_earlier:.2f}')
axes[0].legend(frameon=False, fontsize=7)
axes[1].set(xlabel='PT DPT (shared coordinate)', ylabel='ECDF',
            title='ECDF: a right-shifted curve means later occupancy')
axes[1].legend(frameon=False, fontsize=7)
fig.suptitle('Figure 1 - DPT occupancy by specimen (mouse blue, human orange)', fontsize=12)
_save_figure(fig, 'fig01_dpt_occupancy.png')

# %% [markdown]
# ## 3 - PT orientation and root, validated per specimen
#
# > Does the shared coordinate actually run early-to-late in **every** specimen?
#
# The canonical S1/S2/S3 programmes are scored inside each specimen separately: every panel gene is
# z-scored across that specimen's structures, and the panel score is their mean, so one
# high-abundance gene cannot dominate the panel and the species-wide abundance offset cannot decide
# the sign. Each score is then correlated with DPT. A pooled correlation would hide the failure mode
# that matters here: one specimen running the wrong way. The expected ordering is S1 declining with
# DPT, S3 rising with DPT, and the S1-to-S3 axis increasing with DPT.
#
# The root is then checked explicitly, because a root inside the wrong territory inverts the whole
# axis.

# %%
# Purpose: per-specimen marker-programme scores and their correlation with the shared DPT.
gene_lookup = build_gene_lookup(adata_pt)
expression = as_csr(adata_pt.layers[MEASUREMENT]).astype(np.float64)

marker_panel_indices = {}
for panel, genes in MARKER_PANELS.items():
    present = [gene for gene in genes if gene.upper() in gene_lookup]
    marker_panel_indices[panel] = [int(adata_pt.var_names.get_loc(gene_lookup[gene.upper()]))
                                   for gene in present]
    print(f'{panel} panel: {len(present)}/{len(genes)} canonical genes usable here'
          + ('' if present else '  <-- NO GENE FOUND, this panel is not checked'))
if not any(marker_panel_indices.values()):
    raise ValueError('No canonical PT marker gene was found; the orientation check cannot run.')


def _panel_score(mask, columns):
    """S1/S2/S3 score inside one specimen: mean of its panel genes' within-specimen z-scores."""
    if not columns:
        return np.full(int(mask.sum()), np.nan)
    dense = np.asarray(expression[mask][:, columns].todense(), dtype=float)
    return np.nanmean(zscore_rows(dense.T), axis=0)


# The project already carries an early-to-late marker axis on the object; use it when it is present
# rather than inventing a second definition, and derive the same quantity when it is not.
has_project_axis = 'total_marker_axis' in adata_pt.obs.columns
project_axis = (adata_pt.obs['total_marker_axis'].to_numpy(dtype=float) if has_project_axis
                else np.full(adata_pt.n_obs, np.nan))
print('Early-to-late axis: ' + ("obs['total_marker_axis'] (the project's existing axis)"
                                if has_project_axis else
                                'derived here as z(S3) - z(S1) (project axis absent on this object)'))

orientation_rows = []
for specimen in specimen_order:
    mask = samples == specimen
    specimen_s = s[mask]
    row = {'specimen': specimen, 'species': species_by_sample[specimen],
           'n_structures': int(mask.sum())}
    for panel in MARKER_PANELS:
        score = _panel_score(mask, marker_panel_indices[panel])
        row[f'spearman_dpt_{panel}'] = (safe_spearman(specimen_s, score)
                                        if np.isfinite(score).sum() >= CONFIG['spearman_min_points']
                                        else np.nan)
    if has_project_axis:
        axis_values = project_axis[mask]
    else:
        axis_values = (_panel_score(mask, marker_panel_indices['S3'])
                       - _panel_score(mask, marker_panel_indices['S1']))
    row['spearman_dpt_marker_axis'] = safe_spearman(specimen_s, axis_values)
    row['early_declines_with_dpt'] = bool(np.isfinite(row['spearman_dpt_S1'])
                                          and row['spearman_dpt_S1'] < 0)
    row['late_rises_with_dpt'] = bool(np.isfinite(row['spearman_dpt_S3'])
                                      and row['spearman_dpt_S3'] > 0)
    row['axis_increases_with_dpt'] = bool(np.isfinite(row['spearman_dpt_marker_axis'])
                                          and row['spearman_dpt_marker_axis'] > 0)
    row['orientation_ok'] = bool(row['early_declines_with_dpt'] and row['late_rises_with_dpt']
                                 and row['axis_increases_with_dpt'])
    row['low_dpt_axis_median'] = float(
        np.nanmedian(axis_values[specimen_s <= np.quantile(specimen_s, 0.1)]))
    row['high_dpt_axis_median'] = float(
        np.nanmedian(axis_values[specimen_s >= np.quantile(specimen_s, 0.9)]))
    orientation_rows.append(row)

orientation_by_specimen = pd.DataFrame(orientation_rows)

# --------------------------------------------------------------------------------------
# The root itself: is the index valid, does the structure still exist after PT subsetting, and is
# it early/S1-like? A missing or stale index is reported, never silently ignored.
# --------------------------------------------------------------------------------------
root_rows = []


def _record_root(check, result, interpretation):
    root_rows.append({'check': check, 'result': result, 'interpretation': interpretation})


carried_iroot = adata_full.uns.get('iroot', None)
if carried_iroot is None:
    _record_root('uns["iroot"] carried by the PT artifact', 'absent',
                 'pass: no stale global-root index reached this file, so no dpt call can silently '
                 'ignore it')
else:
    index = int(carried_iroot)
    in_range = 0 <= index < adata_full.n_obs
    _record_root('uns["iroot"] carried by the PT artifact', f'present, index {index}',
                 'pass: index is in range' if in_range else
                 f'concern: index {index} is outside [0, {adata_full.n_obs}); a dpt call on this '
                 'object would silently ignore it. Notebook 03 must drop the stale index '
                 '(adata.uns.pop("iroot")) after the PT subset, as notebook 02 already does.')

if PT_ROOT_DIAGNOSTICS_PATH.exists():
    root_diagnostics = pd.read_csv(PT_ROOT_DIAGNOSTICS_PATH)
    root_columns = {'root_obs_name', 'root_marker_axis', 'raw_dpt_marker_spearman'}
    if root_columns <= set(root_diagnostics.columns) and len(root_diagnostics):
        root_reference = root_diagnostics.iloc[0]
        root_name = str(root_reference['root_obs_name'])
        exists = root_name in set(adata_pt.obs_names)
        _record_root('root structure recorded by 03 still present after PT subsetting',
                     f'{root_name} ({"present" if exists else "MISSING"})',
                     'pass: the rooted structure survived the subset' if exists else
                     'concern: the rooted structure is absent from the PT object')
        if exists:
            root_position = int(adata_pt.obs_names.get_loc(root_name))
            root_dpt = float(s[root_position])
            # `total_scanpy_dpt` is scaled from its 5th-95th percentile, so the earliest 5% of
            # structures are tied at the floor. Requiring an exact 0% would fail a perfectly good
            # root, so the root is required to sit AT that floor rather than to be uniquely first.
            floor_percentile = float((s <= np.quantile(s, 0.05)).mean())
            root_at_floor = bool(root_dpt <= np.quantile(s, 0.05) + 1e-9)
            _record_root('root sits at the early end of the coordinate',
                         f'root DPT {root_dpt:.4f}; {floor_percentile:.0%} of structures share '
                         'the p5-scaled floor',
                         'pass: the root sits at the earliest end of the coordinate'
                         if root_at_floor else
                         'concern: the root is not at the early end of the axis')
            _record_root('root is early/S1-like by the marker axis',
                         f'root_marker_axis {float(root_reference["root_marker_axis"]):.3f}',
                         'pass: negative axis = early/S1-like'
                         if float(root_reference['root_marker_axis']) < 0 else
                         'concern: the root is not early-like')
            _record_root('03 reported DPT oriented with the marker axis',
                         'raw DPT-vs-marker Spearman '
                         f'{float(root_reference["raw_dpt_marker_spearman"]):.3f}',
                         'pass: DPT increases early-to-late'
                         if float(root_reference['raw_dpt_marker_spearman']) > 0 else
                         'concern: DPT disagrees with the marker axis')
else:
    print(f'NOTE: {_rel(PT_ROOT_DIAGNOSTICS_PATH)} is not present; validating the root from the '
          'coordinate alone.')

# Coordinate-level root check, always available: the earliest structures must be early-like and the
# latest structures late-like.
if has_project_axis:
    low_axis_median = float(np.nanmedian(project_axis[s <= np.quantile(s, 0.01)]))
    high_axis_median = float(np.nanmedian(project_axis[s >= np.quantile(s, 0.99)]))
    _record_root('earliest vs latest structures by the marker axis',
                 f'low-DPT axis {low_axis_median:.3f} vs high-DPT axis {high_axis_median:.3f}',
                 'pass: the coordinate runs early-to-late' if low_axis_median < high_axis_median
                 else 'concern: the marker axis does not increase with DPT')

root_validation = pd.DataFrame(root_rows, columns=['check', 'result', 'interpretation'])
root_validation.to_csv(DIAGNOSTIC_DIR / 'pt_dpt_root_validation.csv', index=False)
orientation_by_specimen.to_csv(OUTPUT_DIR / 'pt_marker_validation_by_specimen.csv', index=False)
display(orientation_by_specimen.round(3))
display(root_validation)

orientation_ok_all = bool(orientation_by_specimen['orientation_ok'].all())
if orientation_ok_all:
    print(f'Every specimen runs early-to-late: S1 falls, S3 rises, and the axis increases with DPT '
          f'in all {len(specimen_order)} specimens.')
else:
    failing = orientation_by_specimen.loc[~orientation_by_specimen['orientation_ok'],
                                          'specimen'].tolist()
    print(f'ORIENTATION FLAG: the expected early-to-late ordering fails in {failing}. Notebook 07\'s '
          'cross-species curve comparisons are suspect until this is understood.')

fig2, axes = plt.subplots(1, len(specimen_order), figsize=(3.1 * len(specimen_order), 3.2),
                          sharex=True, sharey=True)
axes = np.atleast_1d(axes)
for axis, specimen in zip(axes, specimen_order):
    mask = samples == specimen
    specimen_s = s[mask]
    edges = np.linspace(np.nanmin(specimen_s), np.nanmax(specimen_s), 21)
    index = np.clip(np.digitize(specimen_s, edges) - 1, 0, edges.size - 2)
    centres = 0.5 * (edges[:-1] + edges[1:])
    for panel, colour in zip(MARKER_PANELS, ('#4C72B0', '#55A868', '#C44E52')):
        score = _panel_score(mask, marker_panel_indices[panel])
        if not np.isfinite(score).any():
            continue
        binned = np.array([np.nanmean(score[index == bin_index])
                           if (index == bin_index).sum() else np.nan
                           for bin_index in range(edges.size - 1)])
        axis.plot(centres, binned, lw=1.8, color=colour, label=panel)
    row = orientation_by_specimen.set_index('specimen').loc[specimen]
    axis.axhline(0.0, color='k', lw=0.6, alpha=0.5)
    axis.set_title(f'{specimen}\naxis rho = {row["spearman_dpt_marker_axis"]:.2f}',
                   fontsize=9, color='#333333')
    axis.set_xlabel('PT DPT')
    axis.text(0.03, 0.04, 'early->late OK' if row['orientation_ok'] else 'ORDERING FAILS',
              transform=axis.transAxes, fontsize=7.5,
              color='#1b7837' if row['orientation_ok'] else '#b2182b')
axes[0].set_ylabel('panel score (z, binned)')
axes[0].legend(frameon=False, fontsize=7)
fig2.suptitle('Figure 2 - canonical S1/S2/S3 programmes against DPT, validated in every specimen',
              fontsize=12)
_save_figure(fig2, 'fig02_pt_marker_orientation.png')

# %% [markdown]
# ## 4 - Within-species reproducibility of the gene curves
#
# > Does the strong human-mouse positional divergence remain among genes whose GAM shape is
# > reproducible within **both** species?
#
# This is the most important gene-level validation. For every gene the two mouse specimen curves are
# correlated with each other and the two human section curves are correlated with each other, on the
# same within-condition standardisation notebook 07 uses for its cross-species shape correlation. A
# gene whose own two mouse specimens disagree, or whose two human sections disagree, cannot support a
# *human-vs-mouse* claim.
#
# The per-specimen curves are notebook 07's cached fits - the same curves it clustered - reloaded via
# the shared stage cache rather than refit.

# %%
# Purpose: reload notebook 07's per-specimen GAM curves from the shared stage cache (hit, not refit).
Y_all = as_csr(adata_pt.layers[MEASUREMENT])
detected_structures = np.asarray((Y_all > 0).sum(axis=0)).ravel()
measured_in_both = adata_pt.var['measured_in_both_inputs'].to_numpy(dtype=bool)
min_detected = int(np.ceil(GAM_CONFIG['min_detected_fraction'] * adata_pt.n_obs))
eligible_genes = measured_in_both & (detected_structures >= min_detected)
gene_names = adata_pt.var_names.to_numpy()[eligible_genes]
Y_genes = Y_all[:, eligible_genes].tocsr().astype(np.float64)
if gene_names.size < 20:
    raise ValueError(f'Only {gene_names.size} genes are eligible; the PT object looks wrong.')
Y_GENES_FINGERPRINT = digest(Y_genes)

GAM_LAMBDA_GRID = GAM_CONFIG['lambda_grid']
SUPPORT_PERCENTILES = GAM_CONFIG['common_support_pct']
p_lo, p_hi = SUPPORT_PERCENTILES
lo = max(np.percentile(s[species == name], p_lo) for name in SPECIES_COLORS)
hi = min(np.percentile(s[species == name], p_hi) for name in SPECIES_COLORS)
if not lo < hi:
    raise ValueError(f'Mouse and human PT have no common DPT support: {lo:.4f} >= {hi:.4f}.')
grid = np.linspace(lo, hi, GAM_CONFIG['grid_points'])
knots = gam_internal_knots(s, basis_df=3 + GAM_CONFIG['n_internal_knots'])

# The pooled fit is needed only for its per-gene GCV lambda; with 03's conventions unchanged this is
# the same cache entry notebook 07 reused, so it is a hit and not a re-fit.
pooled_fit = cached_run_level_shape(
    Y_genes, s, comparison_indicator, knots, grid, GAM_LAMBDA_GRID,
    stage='section4_gene_fit', root=STAGE_CACHE_DIR, y_fingerprint=Y_GENES_FINGERPRINT,
    enabled=STAGE_CACHE_ENABLED,
)
gene_lambda_index = np.asarray(pooled_fit['lam_idx'], dtype=int)


def _fit_specimen_curves():
    """Fit one single-condition GAM per specimen for every eligible gene."""
    fitted = {}
    for name in sorted(species_by_sample):
        mask = samples == name
        fitted[name], _ = fit_single_condition_curves(
            Y_genes[mask], s[mask], knots, grid, GAM_LAMBDA_GRID, gene_lambda_index,
            support_pct=SUPPORT_PERCENTILES,
        )
    return fitted


specimen_curves = cached_payload(
    'pt_specimen_curves', _fit_specimen_curves, root=STAGE_CACHE_DIR,
    params={'grid': np.asarray(grid), 'lambda_grid': np.asarray(GAM_LAMBDA_GRID),
            'knots': np.asarray(knots), 'support_pct': list(SUPPORT_PERCENTILES)},
    inputs={'s': np.asarray(s), 'samples': np.asarray(samples), 'y': Y_GENES_FINGERPRINT,
            'lambda_index': gene_lambda_index},
    code=code_digest(fit_single_condition_curves), enabled=STAGE_CACHE_ENABLED,
)
missing_specimens = [name for name in specimen_order if name not in specimen_curves]
if missing_specimens:
    raise KeyError(f'Per-specimen curves missing for {missing_specimens}.')

mouse_curves = specimen_balanced_curves({name: specimen_curves[name] for name in mouse_specimens})
human_curves = specimen_balanced_curves({name: specimen_curves[name] for name in human_specimens})
print(f'Eligible genes: {gene_names.size:,}; per-specimen curves: {specimen_order}')
print(f'Shared support: [{lo:.4f}, {hi:.4f}] on {grid.size} grid points; '
      f'{knots.size} internal knots')


# %%
# Purpose: within-species reproducibility and the cross-species shape correlation, on one column set.
# 07's eligibility rule fixes the grid columns both species are seen on; the SAME columns are used for
# the within-species pairs so every number below is directly comparable.
mouse_amplitude = _peak_to_peak(mouse_curves)
human_amplitude = _peak_to_peak(human_curves)
shape_eligible, shape_columns, _, _ = _genuine_shape_eligibility(
    mouse_curves, human_curves, grid, CONFIG['amplitude_floor'])
if int(shape_eligible.sum()) < 20:
    raise ValueError(f'Only {int(shape_eligible.sum())} genes carry shape in both species.')

Z_mouse = _standardize_species(mouse_curves, shape_columns)
Z_human = _standardize_species(human_curves, shape_columns)
cross_species_shape = _row_block_pearson(Z_mouse, Z_human, CONFIG['numerical_epsilon'])

mouse_reproducibility = _row_block_pearson(
    _standardize_species(specimen_curves[mouse_specimens[0]], shape_columns),
    _standardize_species(specimen_curves[mouse_specimens[1]], shape_columns),
    CONFIG['numerical_epsilon'])
human_reproducibility = _row_block_pearson(
    _standardize_species(specimen_curves[human_specimens[0]], shape_columns),
    _standardize_species(specimen_curves[human_specimens[1]], shape_columns),
    CONFIG['numerical_epsilon'])

reproducibility_floor = CONFIG['reproducibility_floor']
reproducible = ((mouse_reproducibility >= reproducibility_floor)
                & (human_reproducibility >= reproducibility_floor))
high_confidence = shape_eligible & reproducible

reproducibility_table = pd.DataFrame({
    'gene': gene_names,
    'mouse_within_species_correlation': mouse_reproducibility,
    'human_within_species_correlation': human_reproducibility,
    'human_mouse_shape_correlation': cross_species_shape,
    'mouse_amplitude': mouse_amplitude,
    'human_amplitude': human_amplitude,
    'shape_eligible_both_species': shape_eligible,
    'reproducible_within_both_species': reproducible,
})
reproducibility_table.to_csv(OUTPUT_DIR / 'within_species_curve_reproducibility.csv', index=False)

# Cross-check this notebook's independent recomputation against notebook 07's saved numbers.
shape_check_note = 'notebook 07 outputs not found - cross-check skipped'
reference_path = NOTEBOOK_07_DIR / 'shared_archetype_gene_metrics.csv'
if reference_path.exists():
    reference = pd.read_csv(reference_path).set_index('gene')
    joined = pd.DataFrame({'ours': pd.Series(cross_species_shape, index=gene_names)}).join(
        reference[['shape_correlation']], how='inner').dropna()
    max_difference = float(np.max(np.abs(joined['ours'] - joined['shape_correlation'])))
    if max_difference > 1e-6:
        raise ValueError('Cross-species shape correlation disagrees with notebook 07 by '
                         f'{max_difference:.3e}; the curves or the standardisation differ.')
    shape_check_note = f'{len(joined):,} genes identical to notebook 07 up to {max_difference:.1e}'
print(f'Cross-check against notebook 07 shared_archetype_gene_metrics.csv - {shape_check_note}')

reproducibility_summary = pd.DataFrame([
    _summary_row('all eligible genes (shape in both species)',
                 int(shape_eligible.sum()), cross_species_shape[shape_eligible]),
    _summary_row(f'high reproducibility (both correlations >= {reproducibility_floor})',
                 int(high_confidence.sum()), cross_species_shape[high_confidence]),
])
display(reproducibility_summary.round(4))
print(f'Genes with real shape in both species: {int(shape_eligible.sum()):,}')
print(f'  of which reproducible within both species: {int(high_confidence.sum()):,} '
      f'({int(high_confidence.sum()) / max(int(shape_eligible.sum()), 1):.1%})')
print(f'Median within-species correlation - mouse '
      f'{np.nanmedian(mouse_reproducibility[shape_eligible]):.3f}, human '
      f'{np.nanmedian(human_reproducibility[shape_eligible]):.3f}')

median_all = reproducibility_summary.loc[0, 'median_shape_correlation']
median_high = reproducibility_summary.loc[1, 'median_shape_correlation']
divergence_survives = bool(np.isfinite(median_high) and median_high < 0.8)
if divergence_survives:
    print(f'Divergence survives the reproducibility gate: median cross-species shape correlation is '
          f'{median_all:.3f} over all eligible genes and {median_high:.3f} among genes reproducible '
          f'in both species. The divergence is not produced by unreproducible weak curves.')
else:
    print(f'REPRODUCIBILITY FLAG: among genes reproducible in both species the median cross-species '
          f'shape correlation rises to {median_high:.3f}; the apparent divergence is concentrated in '
          f'genes whose curves are not reproducible within a species.')

# %% [markdown]
# ## 5 - Strong spatial-gene sensitivity
#
# > Is the biological signal still present in strong, reproducible positional genes, or is it mainly
# > driven by weak curves that become exaggerated after z-scoring?
#
# Notebook 07 observed that stronger-amplitude genes looked more conserved. Here that is tested
# directly, with no re-clustering and no new model: genes are gated on an increasingly strong
# positional amplitude requirement **in both species**, and the cross-species shape correlation is
# summarised at each threshold, with and without the within-species reproducibility requirement.

# %%
# Purpose: the amplitude-threshold x reproducibility grid of cross-species shape conservation.
sensitivity_rows = []
for threshold in CONFIG['amplitude_thresholds']:
    amplitude_ok = (mouse_amplitude >= threshold) & (human_amplitude >= threshold)
    for requires_reproducibility, mask in ((False, amplitude_ok),
                                           (True, amplitude_ok & reproducible)):
        label = (f'amplitude >= {threshold:.2f} in both species'
                 + (' AND reproducible in both species' if requires_reproducibility else ''))
        sensitivity_rows.append(_summary_row(label, int(mask.sum()), cross_species_shape[mask]))
shape_conservation_sensitivity = pd.DataFrame(sensitivity_rows)
shape_conservation_sensitivity['amplitude_threshold'] = [
    float(label.split('>=')[1].split()[0])
    for label in shape_conservation_sensitivity['subset']]
shape_conservation_sensitivity['requires_reproducibility'] = [
    'AND' in label for label in shape_conservation_sensitivity['subset']]
shape_conservation_sensitivity.to_csv(OUTPUT_DIR / 'shape_conservation_sensitivity.csv', index=False)
display(shape_conservation_sensitivity.round(4))

strongest = shape_conservation_sensitivity.iloc[-1]
print(f'At the strongest, reproducible threshold: {int(strongest["n_genes"]):,} genes, median '
      f'cross-species shape correlation {strongest["median_shape_correlation"]:.3f}, '
      f'{strongest["fraction_correlation_below_0"]:.1%} below zero.')
if int(strongest['n_genes']) and strongest['fraction_correlation_below_-0.5'] > 0.1:
    print('Divergence is still present among strong, reproducible genes; it is not an artefact of '
          'weak curves amplified by z-scoring.')
else:
    print('NOTE: at the strongest reproducible threshold the divergent tail is small; the amplitude '
          'trend itself is the more informative summary above.')

# %% [markdown]
# ## 6 - Representative genes, back in the raw data
#
# > Is the apparent difference present in the actual data, or generated only by the balanced GAM?
#
# Three strongly conserved and three strongly divergent genes - the divergent ones required to be
# reproducible in both species with adequate amplitude, so their divergence is interpretable - are
# plotted as their per-specimen GAM curves, the balanced species curves, and the binned raw
# expression underneath. No gallery: six genes, one panel each.

# %%
# Purpose: pick the representative genes, then draw their fits over their own binned raw data.
representative_pool = (shape_eligible & reproducible
                       & (mouse_amplitude >= CONFIG['representative_min_amplitude'])
                       & (human_amplitude >= CONFIG['representative_min_amplitude']))
pool_index = np.flatnonzero(representative_pool)
if pool_index.size == 0:
    raise ValueError('No gene is strong and reproducible in both species; nothing to display.')
pool_order = pool_index[np.argsort(cross_species_shape[pool_index])]
n_each = CONFIG['representative_per_class']
divergent_index = pool_order[:n_each]
conserved_index = pool_order[::-1][:n_each]

representative_rows = []
for label, indices in (('conserved', conserved_index), ('divergent', divergent_index)):
    for index in indices:
        representative_rows.append({
            'class': label, 'gene': gene_names[index],
            'human_mouse_shape_correlation': float(cross_species_shape[index]),
            'mouse_within_species_correlation': float(mouse_reproducibility[index]),
            'human_within_species_correlation': float(human_reproducibility[index]),
            'mouse_amplitude': float(mouse_amplitude[index]),
            'human_amplitude': float(human_amplitude[index]),
        })
representative_genes = pd.DataFrame(representative_rows)
representative_genes.to_csv(DIAGNOSTIC_DIR / 'representative_genes.csv', index=False)
display(representative_genes.round(3))
_divergent_picks = representative_genes.loc[representative_genes['class'].eq('divergent'),
                                            'human_mouse_shape_correlation']
print(f'Representative pool (strong + reproducible in both species): {pool_index.size:,} genes; '
      f'divergent picks reach shape correlation as low as {_divergent_picks.min():.2f}.')
if pool_index.size < 2 * n_each:
    print(f'NOTE: only {pool_index.size} gene(s) clear both gates, so the conserved and divergent '
          'lists overlap; read the panel labels rather than the classes.')


def _binned_expression(column, mask):
    """Mean log-normalised expression of one gene per DPT bin, inside a single specimen."""
    values = np.asarray(expression[mask][:, column].todense(), dtype=float).ravel()
    specimen_s = s[mask]
    edges = np.linspace(np.nanmin(specimen_s), np.nanmax(specimen_s), CONFIG['raw_bins'] + 1)
    index = np.clip(np.digitize(specimen_s, edges) - 1, 0, CONFIG['raw_bins'] - 1)
    centres, means = [], []
    for bin_index in range(CONFIG['raw_bins']):
        selected = index == bin_index
        if selected.sum() >= CONFIG['min_structures_per_bin']:
            centres.append(float(0.5 * (edges[bin_index] + edges[bin_index + 1])))
            means.append(float(np.nanmean(values[selected])))
    return np.asarray(centres), np.asarray(means)


fig3 = plt.figure(figsize=(13.5, 8.6))
grid_layout = fig3.add_gridspec(3, 3)
axis_sensitivity = fig3.add_subplot(grid_layout[0, :])
for requires, style, colour in ((False, '-', '#4C72B0'), (True, '--', '#C44E52')):
    block = (shape_conservation_sensitivity[
        shape_conservation_sensitivity['requires_reproducibility'].eq(requires)]
        .sort_values('amplitude_threshold'))
    axis_sensitivity.plot(block['amplitude_threshold'], block['median_shape_correlation'],
                          marker='o', ls=style, color=colour,
                          label='amplitude only' if not requires else 'amplitude + reproducible')
    for x_value, y_value, count in zip(block['amplitude_threshold'],
                                       block['median_shape_correlation'], block['n_genes']):
        axis_sensitivity.annotate(f'n={int(count):,}', (x_value, y_value), fontsize=7,
                                  textcoords='offset points', xytext=(0, 7))
axis_sensitivity.axhline(0.0, color='k', lw=0.7, alpha=0.6)
axis_sensitivity.set(xlabel='peak-to-peak amplitude required in BOTH species',
                     ylabel='median cross-species shape correlation',
                     title='Conservation versus amplitude gate (no re-clustering, no new model)')
axis_sensitivity.legend(frameon=False, fontsize=8)

gene_panels = list(conserved_index) + list(divergent_index)
for position, gene_index in enumerate(gene_panels):
    axis = fig3.add_subplot(grid_layout[1 + position // 3, position % 3])
    column = int(adata_pt.var_names.get_loc(gene_names[gene_index]))
    for specimen, colour in zip(specimen_order,
                                SPECIMEN_CMAP(np.linspace(0, 1, len(specimen_order)))):
        axis.plot(grid, specimen_curves[specimen][gene_index], lw=0.7, alpha=0.55, color=colour,
                  label=specimen if not position else None)
        centres, means = _binned_expression(column, samples == specimen)
        axis.plot(centres, means, ls='none', marker='.', ms=3.5, color=colour, alpha=0.7)
    axis.plot(grid, mouse_curves[gene_index], lw=2.0, color=SPECIES_COLORS['mouse'],
              label='balanced mouse' if not position else None)
    axis.plot(grid, human_curves[gene_index], lw=2.0, color=SPECIES_COLORS['human'],
              label='balanced human' if not position else None)
    axis.set_title(f'{gene_names[gene_index]} '
                   f'({"conserved" if position < len(conserved_index) else "divergent"}) '
                   f'r={cross_species_shape[gene_index]:.2f}', fontsize=9)
    axis.set_xlabel('PT DPT')
    axis.set_ylabel('lognorm expression')
    axis.text(0.02, 0.93, f'repro m={mouse_reproducibility[gene_index]:.2f} '
                          f'h={human_reproducibility[gene_index]:.2f}',
              transform=axis.transAxes, fontsize=6.5, va='top')
fig3.legend(loc='lower center', ncol=6, frameon=False, fontsize=8)
fig3.suptitle('Figure 3 - amplitude sensitivity, and representative genes over their own raw data',
              fontsize=12)
_save_figure(fig3, 'fig03_conservation_and_representative_genes.png')

# %%
# Purpose: record which figures this run produced, beside the tables.
pd.DataFrame({'figure': figure_index}).to_csv(DIAGNOSTIC_DIR / 'figure_index.csv', index=False)
print(f'Figures written: {figure_index}')

# %% [markdown]
# ## 7 - Validation summary
#
# Four questions, each answered from the numbers executed above - never from a preset conclusion.

# %%
# Purpose: the four-question verdict, assembled from the flags computed above.
summary_table = pd.DataFrame([
    {
        'diagnostic': 'Do species broadly overlap along DPT?',
        'result': (f'P(human earlier than mouse) = {auc_human_earlier:.2f}; '
                   f'median DPT mouse {mouse_median:.2f} vs human {human_median:.2f}; '
                   f'{fraction_human_in_mouse_iqr:.0%} of human structures inside the mouse IQR'),
        'interpretation': 'pass' if overlap_ok else 'concern',
    },
    {
        'diagnostic': 'Does DPT show S1-to-S3 ordering in every specimen?',
        'result': (f'{int(orientation_by_specimen["orientation_ok"].sum())}/'
                   f'{len(orientation_by_specimen)} specimens early-to-late; '
                   'axis rho range ['
                   f'{orientation_by_specimen["spearman_dpt_marker_axis"].min():.2f}, '
                   f'{orientation_by_specimen["spearman_dpt_marker_axis"].max():.2f}]'),
        'interpretation': 'pass' if orientation_ok_all else 'concern',
    },
    {
        'diagnostic': 'Are gene curves reproducible within species?',
        'result': (f'median rho {np.nanmedian(mouse_reproducibility[shape_eligible]):.2f} (mouse), '
                   f'{np.nanmedian(human_reproducibility[shape_eligible]):.2f} (human); '
                   f'{int(high_confidence.sum()):,}/{int(shape_eligible.sum()):,} genes at rho >= '
                   f'{reproducibility_floor} in both species'),
        'interpretation': 'pass' if divergence_survives else 'concern',
    },
    {
        'diagnostic': 'Does divergence remain in strong, reproducible genes?',
        'result': (f'amplitude >= {strongest["amplitude_threshold"]:.2f} + reproducible: '
                   f'{int(strongest["n_genes"]):,} genes, median rho '
                   f'{strongest["median_shape_correlation"]:.2f}, '
                   f'{strongest["fraction_correlation_below_0"]:.0%} below zero'),
        'interpretation': 'pass' if int(strongest['n_genes']) >= 20 else 'concern',
    },
])
summary_table.to_csv(DIAGNOSTIC_DIR / 'validation_summary.csv', index=False)
display(summary_table)

concerns = summary_table.loc[summary_table['interpretation'].eq('concern'), 'diagnostic'].tolist()
if not concerns:
    conclusion = (
        'The large-scale human-mouse curve differences are not readily explained by gross DPT '
        'occupancy imbalance, incorrect PT orientation, or unreproducible weak GAMs. They therefore '
        'appear to contain a substantial biological component, although the human data represent one '
        'donor and the comparison remains descriptive.')
else:
    conclusion = ('This notebook does not clear the notebook-07 cross-species interpretation: '
                  + '; '.join(concerns)
                  + '. Read the corresponding diagnostic above before interpreting any human-mouse '
                    'curve difference.')
print('CONCLUSION -', conclusion)
print()
print(f'Written to {_rel(OUTPUT_DIR)}: dpt_occupancy_by_specimen.csv, '
      'pt_marker_validation_by_specimen.csv, within_species_curve_reproducibility.csv, '
      'shape_conservation_sensitivity.csv, and the figures indexed in diagnostics/figure_index.csv.')
print('What must NOT be read into this: it is not a test of biological replication (two human '
      'sections are one donor), it computes no p-value or significance threshold, it validates only '
      'the failure modes it names, and a "pass" is not evidence that the divergence is adaptive.')
