# %% [markdown]
# # Human versus healthy-mouse proximal tubule: continuous spatial conservation and rewiring
#
# Notebook 05 asked which genes are higher or lower across the whole proximal tubule (PT). That
# question averages one coordinate away: a gene whose zonation is intact but displaced, or whose
# gradient is inverted, contributes nothing to a whole-PT mean difference. This notebook asks the
# spatial question instead, and decomposes every one-to-one ortholog into
#
# ```
# level difference  +  zonation amplitude difference  +  positional displacement  +  residual shape divergence
# ```
#
# The objects that carry the analysis are **fitted expression curves along PT pseudospace** — a
# balanced mouse curve and a balanced human curve per gene, each the equal-weight mean of its
# specimens' own fitted curves. Absolute abundance and spatial pattern are then analysed separately
# by construction: `level` is the mean vertical offset of the two curves, `amplitude` is peak-to-trough
# of each curve, `shape` compares curves after each is standardised within its species, and `phase` is
# a constrained horizontal displacement that has to earn its interpretation against a fixed maximum.
#
# ## Read this before the numbers
#
# * **The two human slices (`HUK1_COR1`, `HUK1_MED1`) come from one donor.** They are not independent
#   human biological replicates and no column in this notebook is a population-level human-versus-mouse
#   test. Every contrast is a descriptive effect size on a specific balanced-curve pair.
# * Agreement between the two human slices, and between the two mouse specimens, is used as a
#   **robustness** measure — "is this curve reproducible within a species" — never as replication-based
#   inference. Section 5 does leave-one-slice-out as *slice sensitivity*, not as a replication test.
# * **Species is inseparable from sampling and from Harmony batch correction** in this cohort
#   (`sample` is the batch key). A cross-species difference in a fitted curve is a difference between
#   these two mice and this one donor, processed this way.
# * Both human slices are anatomically cortex; the `MED1` label does not make the tissue medullary.
# * Genes that were used to build or orient the PT coordinate are flagged (`axis_basis_gene`) and are
#   excluded from every unbiased discovery list, but they are retained for validation plots.
# * Direction conventions, used everywhere below: **positive level effect = higher in human**;
#   **positive shift = the human positional program occurs later along PT than the mouse program**;
#   positive early-to-late gradient = expression increases from early to late PT.
#
# ## Layout
#
# | Section | Question it answers |
# | --- | --- |
# | 0 | Setup, cohort, gene universe, direction conventions |
# | 1 | Balanced species curves and within-species reproducibility |
# | 2 | Landmark-anchored cross-species registration of PT position |
# | 3 | Level / amplitude / shape / position / phase / inversion metrics per gene |
# | 4 | Mutually exclusive spatial phenotypes |
# | 5 | Sensitivity of every phenotype to the analysis choices |
# | 6 | Primary figures (1A–1E) |
# | 7 | Clustering of the residual "complex shape rewiring" curves only |
# | 8 | Pathway enrichment against the spatial phenotypes |
# | 9 | Pathway-level gene x pseudospace visualisation |
# | 10 | Conventional whole-PT signal versus continuous spatial signal |
# | 11 | Spatially conserved functions |
# | 12 | Literature validation (deferred — placeholder) |
# | 13 | Tables, notes and the auto-generated summary |
#

# %%
# Purpose: configuration, direction conventions, roots and the output layout for this notebook.

from __future__ import annotations

import os
import re
import sys
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import sparse
from scipy.stats import hypergeom, pearsonr, spearmanr
from IPython.display import display


def _find_project_dir() -> Path:
    """Walk up from the working directory (and this file) until the repository root is found."""
    starts = [Path.cwd().resolve()]
    if "__file__" in globals():
        starts.insert(0, Path(__file__).resolve().parent)

    for start in starts:
        for candidate in (start, *start.parents):
            if (candidate / "pseudospace").is_dir() and (candidate / "data").is_dir():
                return candidate

    raise RuntimeError("Could not locate the repository root containing pseudospace/ and data/.")


PROJECT_DIR = _find_project_dir()

if str(PROJECT_DIR) not in sys.path:
    sys.path.insert(0, str(PROJECT_DIR))

# Roots: the PSEUDOSPACE_* environment variables win, otherwise the repository layout is used.
DATA_ROOT = Path(os.environ.get("PSEUDOSPACE_DATA_ROOT", PROJECT_DIR / "data")).resolve()
RESULTS_ROOT = Path(os.environ.get("PSEUDOSPACE_RESULTS_ROOT", PROJECT_DIR / "results")).resolve()

UPSTREAM_DIR = RESULTS_ROOT / "human_vs_healthy_mouse"
DPT_OUTPUT_PATH = UPSTREAM_DIR / "cross_species_pt_dpt.h5ad"
ORTHOLOG_MAP_PATH = UPSTREAM_DIR / "ortholog_map_used.csv"
PATHWAY_MEMBERSHIP_PATH = UPSTREAM_DIR / "curves" / "pathway_membership_coverage.csv"
PATHWAY_REDUNDANCY_GROUPS_PATH = UPSTREAM_DIR / "curves" / "pathway_redundancy_groups.csv"

# 05's saved objects. These are read by the cross-check cell only, never as an analysis input.
PREVIOUS_CURVES_PATH = UPSTREAM_DIR / "gene_rewiring" / "balanced_gene_curves.npz"
PREVIOUS_ATLAS_PATH = UPSTREAM_DIR / "gene_rewiring" / "gene_curve_atlas.csv"

OUT_DIR = UPSTREAM_DIR / "spatial_rewiring"
FIG_DIR = OUT_DIR / "figures"
TABLE_DIR = OUT_DIR / "tables"

# The stage cache is deliberately 03/05's directory rather than a new one: the pooled fit is keyed on
# the fitter's source plus the fit parameters, so pointing here turns that stage into a cache hit
# instead of a silent recomputation.
STAGE_CACHE_DIR = UPSTREAM_DIR / "stage_cache"

for _directory in (OUT_DIR, FIG_DIR, TABLE_DIR, STAGE_CACHE_DIR):
    _directory.mkdir(parents=True, exist_ok=True)

from pseudospace.levelshape import fit_single_condition_curves, summarize_curve_effects
from pseudospace.specimen import specimen_balanced_curves
from pseudospace.pathways import summarize_pathway_redundancy
from pseudospace.stats_gam import as_csr, bh_adjust, gam_internal_knots
from pseudospace.stage_cache import (
    cache_status,
    cached_payload,
    cached_run_level_shape,
    digest,
    stage_cache_enabled,
)

STAGE_CACHE_ENABLED = stage_cache_enabled(None)

# ======================================================================================
# Analysis settings. Every threshold the notebook uses lives here, and the diagnostics cell in
# section 4 reports the empirical distributions these values are read against.
# ======================================================================================

CONFIG = {
    # --- curve fitting (identical to 03/05 so the shared cache applies) ------------------
    "grid_points": 101,
    "common_support_pct": (1, 99),
    "n_internal_knots": 9,
    "lambda_grid": np.logspace(-3, 3, 13),

    # --- gene universe -------------------------------------------------------------------
    "min_detected_fraction_all": 0.02,     # as in 03/05: >= 2% of all PT structures
    "min_detection_each_species": 0.05,    # a landmark or pathway member must be seen in both
    "min_pattern_grid_points": 40,         # grid points both species must support before any metric

    # --- zonation ------------------------------------------------------------------------
    "amplitude_patterned": 0.15,           # peak-to-trough above this = a patterned curve
    "amplitude_flat": 0.075,               # ... and below this = a flat curve
    "amplitude_degenerate": 0.05,          # below this the curve has no usable positional signal

    # --- within-species reproducibility ---------------------------------------------------
    "within_species_corr": 0.50,           # specimen/slice curves must agree at least this well

    # --- shape (standardised, level- and amplitude-free) ----------------------------------
    "conserved_corr": 0.70,                # high shape agreement
    "rewired_max_corr": 0.30,              # low shape agreement

    # --- amplitude change -----------------------------------------------------------------
    "amplitude_change_major": 0.75,        # |log2 human/mouse amplitude|, ~1.68x
    "amplitude_change_small": 0.40,        # ~1.32x: still "no major amplitude change"

    # --- phase ---------------------------------------------------------------------------
    "max_registration_shift": 0.20,        # never search further than 20% of PT (unit axis)
    "shift_step": 0.01,
    "minimum_shift": 0.08,                 # a shift must exceed this to be interpreted
    "registered_corr": 0.80,               # ... and must reach this shape agreement
    "shift_improvement": 0.10,             # ... and improve on the unshifted correlation by this

    # --- gradient inversion ---------------------------------------------------------------
    "gradient_min_abs": 0.10,              # both species need an appreciable early-to-late slope
    "inversion_max_corr": 0.30,

    # --- pathway enrichment ---------------------------------------------------------------
    "pathway_min_members": 10,
    "pathway_max_members": 500,
    "enrichment_fdr": 0.05,
    "enrichment_min_members_flagged": 3,

    # --- clustering of residual rewiring --------------------------------------------------
    "cluster_k_min": 4,
    "cluster_k_max": 8,
    "cluster_seed": 0,
    "cluster_min_genes": 60,               # fewer than this and the section is skipped, not forced

    # --- figures --------------------------------------------------------------------------
    "top_labels_per_panel": 12,
    "representative_genes_per_class": 3,
    "heatmap_genes_per_class": 6,
}

SPECIES_COLORS = {
    "mouse": "#0072B2",
    "human": "#D55E00",
}

# Direction conventions, restated once so the sign of every column is unambiguous.
DIRECTION_CONVENTIONS = {
    "level_effect_human_minus_mouse": "positive = higher in human",
    "amplitude_log2_ratio_human_over_mouse": "positive = stronger zonation in human",
    "shift_human_minus_mouse": "positive = human program occurs later along PT",
    "gradient_early_to_late": "positive = expression rises from early to late PT",
    "pattern_rms_z": "standardised RMS difference; 0 = identical spatial pattern",
}

# Canonical PT landmark programs, taken from 03 (NEPHRON_AXIS_MARKERS and the PT fine panels) with
# the S3 additions used for the mouse-only axis. They anchor the cross-species registration in
# section 2 and are excluded from unbiased discovery.
LANDMARK_PANELS = {
    "S1": ["Slc5a2", "Slc5a12", "Gatm", "Lrp2", "Cubn", "Slc34a1"],
    "S2": ["Slc22a6", "Slc13a3", "Cyp2e1"],
    "S3": ["Slc22a7", "Cyp7b1", "Slc7a13", "Slc6a18", "Acsm3"],
}

# Every gene that contributed to building or orienting the PT coordinate.
AXIS_BASIS_GENES = tuple(sorted({gene for panel in LANDMARK_PANELS.values() for gene in panel}))

# Obvious technical features: mitochondrial and cytoplasmic ribosomal protein genes.
TECHNICAL_GENE_PATTERNS = (
    r"^MT-",            # human mitochondrial
    r"^MT(ND|CO|ATP|CYB|RNR)\d",   # mouse mitochondrial (Mt-Nd1 ...)
    r"^RP[SL]\d+[A-Z]?$",          # cytoplasmic ribosomal proteins
    r"^MRP[SL]\d+$",               # mitochondrial ribosomal proteins
)


def _save_figure(fig, name, dpi=240):
    """Write a figure into `figures/` and report the path relative to the repository root."""
    path = FIG_DIR / name
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    print(f"saved {path.relative_to(PROJECT_DIR)}")
    return path


def _save_table(frame, name, **kwargs):
    """Write a table into `tables/` and report the path and the shape that was written."""
    path = TABLE_DIR / name
    frame.to_csv(path, index=False, **kwargs)
    print(f"saved {path.relative_to(PROJECT_DIR)}  ({len(frame):,} rows x {frame.shape[1]} cols)")
    return path


def _report(label, kept, total=None, why="", indent=2):
    """One consistent line for 'how many features entered this step and what was dropped'."""
    pad = " " * indent
    if total is None:
        print(f"{pad}{label}: {kept:,}{(' ' + why) if why else ''}")
    else:
        dropped = total - kept
        share = f" ({dropped / total:.1%})" if total else ""
        print(f"{pad}{label}: {kept:,} of {total:,} kept, {dropped:,} dropped{share}"
              f"{(' - ' + why) if why else ''}")


plt.rcParams.update({
    "figure.dpi": 110,
    "savefig.dpi": 240,
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "legend.frameon": False,
    "legend.fontsize": 8,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

print("Project:      ", PROJECT_DIR)
print("Input:        ", DPT_OUTPUT_PATH)
print("Outputs:      ", OUT_DIR)
print("Stage cache:  ", STAGE_CACHE_DIR, f"({'on' if STAGE_CACHE_ENABLED else 'off'})")
print("Direction conventions:")
for _column, _meaning in DIRECTION_CONVENTIONS.items():
    print(f"  {_column}: {_meaning}")
print("Landmark panels:", {k: len(v) for k, v in LANDMARK_PANELS.items()})


# %% [markdown]
# ## 0.2 - The PT cohort and the shared coordinate
#
# `cross_species_pt_dpt.h5ad` is 03's PT-only object: its `total_scanpy_dpt` column is the PT-specific
# pseudospace coordinate, not the global nephron one, and it is the coordinate every curve below is
# fitted against. The two species are compared only over the stretch of that coordinate **both**
# species occupy (`common_support_pct`), because outside a species' own support a fitted curve is an
# extrapolation rather than an observation.
#
# The cohort is what it is: two control mouse specimens and two human slices from a single donor.
# The cell below records it, and nothing downstream multiplies those units into a p-value.
#

# %%
# Purpose: load 03's PT object, define the cohort and the shared coordinate every fit below uses.

if not DPT_OUTPUT_PATH.exists():
    raise FileNotFoundError(f"{DPT_OUTPUT_PATH} is missing - run 03_human_vs_healthy_mouse.ipynb first.")

adata = sc.read_h5ad(DPT_OUTPUT_PATH)

if "lognorm" not in adata.layers:
    raise KeyError(f"{DPT_OUTPUT_PATH.name} has no layers['lognorm'].")

required_obs = {"comparison_species", "sample", "total_scanpy_dpt"}
missing_obs = sorted(required_obs - set(adata.obs.columns))

if missing_obs:
    raise KeyError(f"Missing required obs columns in {DPT_OUTPUT_PATH.name}: {missing_obs}")

mask = np.asarray(
    adata.obs["comparison_species"].astype(str).isin(["mouse", "human"])
    & np.isfinite(adata.obs["total_scanpy_dpt"].to_numpy(dtype=float))
).ravel()

adata_pt = adata[mask].copy()
adata_pt.X = adata_pt.layers["lognorm"].copy()

species = adata_pt.obs["comparison_species"].astype(str).to_numpy()
samples = adata_pt.obs["sample"].astype(str).to_numpy()
s = adata_pt.obs["total_scanpy_dpt"].to_numpy(dtype=float)

sample_species = (
    adata_pt.obs[["sample", "comparison_species"]]
    .drop_duplicates()
    .sort_values(["comparison_species", "sample"])
    .reset_index(drop=True)
)

mouse_samples = sample_species.loc[sample_species["comparison_species"].eq("mouse"), "sample"].tolist()
human_samples = sample_species.loc[sample_species["comparison_species"].eq("human"), "sample"].tolist()

if not mouse_samples or not human_samples:
    raise ValueError(f"Both species must be present: mouse={mouse_samples}, human={human_samples}")

region_column = "region" if "region" in adata_pt.obs.columns else None

aggregation = {
    "n_pt_structures": ("total_scanpy_dpt", "size"),
    "dpt_min": ("total_scanpy_dpt", "min"),
    "dpt_max": ("total_scanpy_dpt", "max"),
}

if region_column:
    aggregation["region"] = (region_column, "first")

cohort_summary = (
    adata_pt.obs.groupby(["comparison_species", "sample"], observed=True)
    .agg(**aggregation)
    .reset_index()
)

cohort_summary = cohort_summary[
    ["comparison_species", "sample"] + (["region"] if region_column else [])
    + ["n_pt_structures", "dpt_min", "dpt_max"]
]

cohort_summary["inference_unit"] = np.where(
    cohort_summary["comparison_species"].eq("human"),
    "slice (one donor - sensitivity view only)",
    "specimen",
)

display(cohort_summary.round(4))

_report("PT structures", adata_pt.n_obs)
print(f"  mouse specimens: {mouse_samples}   (independent units)")
print(f"  human slices:    {human_samples}   (one donor: sensitivity views, not replicates)")

# --------------------------------------------------------------------------------------
# Shared support: the stretch of pseudospace both species actually occupy.
# --------------------------------------------------------------------------------------

p_lo, p_hi = CONFIG["common_support_pct"]

lo = max(np.percentile(s[species == specimen_species], p_lo)
         for specimen_species in ("mouse", "human"))
hi = min(np.percentile(s[species == specimen_species], p_hi)
         for specimen_species in ("mouse", "human"))

if not lo < hi:
    raise ValueError(f"No shared PT support: {lo:.4f} >= {hi:.4f}")

grid = np.linspace(lo, hi, CONFIG["grid_points"])

# Unit-scaled position along the shared support. Horizontal displacement is reported in these units
# (fraction of PT), because the raw DPT scale is not comparable between tissues.
grid_unit = (grid - grid[0]) / (grid[-1] - grid[0])

knots = gam_internal_knots(s, basis_df=3 + CONFIG["n_internal_knots"])

_save_table(cohort_summary, "pt_cohort_summary.csv")

print(f"  shared PT support: [{lo:.4f}, {hi:.4f}] on the PT DPT coordinate")
print(f"  grid: {grid.size} points; internal knots: {len(knots)}")


# %% [markdown]
# ## 0.3 - Tested orthologs, per-specimen detection and discovery flags
#
# The analysis universe is the one 03/05 use: one-to-one orthologs detected in at least 2% of PT
# structures and measured in **both** inputs. A gene the ortholog map created as an unfed column is a
# structural zero, not a measured zero, so `measured_in_both_inputs` is required and reported here.
#
# Two flags decide what may appear in an unbiased discovery list:
#
# * `axis_basis_gene` - the gene helped build or orient the PT coordinate (`LANDMARK_PANELS`). Keeping
#   them would make the coordinate's own basis look like the strongest finding in it. They stay in the
#   tables and are used in the registration and validation plots.
# * `technical_gene` - mitochondrial and cytoplasmic/mitochondrial ribosomal protein genes.
#
# Neither flag removes a gene from the object: everything is reported with its flags so the exclusion
# is auditable rather than silent.
#

# %%
# Purpose: the tested ortholog universe, per-specimen detection / abundance, and the two flags.

if "measured_in_both_inputs" in adata_pt.var.columns:
    measured_in_both = adata_pt.var["measured_in_both_inputs"].to_numpy(dtype=bool)
    measured_source = "var['measured_in_both_inputs'] written by 03"
else:
    measured_in_both = np.ones(adata_pt.n_vars, dtype=bool)
    measured_source = "column absent - every retained gene assumed measured in both inputs"
    print("WARNING: measured_in_both_inputs absent; an unfed column would be read as a measured zero.")

Y_all = as_csr(adata_pt.layers["lognorm"])
detected_all = np.asarray((Y_all > 0).sum(axis=0)).ravel()
min_detected = int(np.ceil(CONFIG["min_detected_fraction_all"] * adata_pt.n_obs))

tested = (detected_all >= min_detected) & measured_in_both
tested_index = np.flatnonzero(tested)

gene_names = adata_pt.var_names.to_numpy().astype(str)[tested]
Y_genes = Y_all[:, tested_index].tocsr().astype(np.float64)

_report(
    "tested orthologs", len(gene_names), adata_pt.n_vars,
    why=f"detection >= {CONFIG['min_detected_fraction_all']:.0%} of {adata_pt.n_obs:,} PT structures "
        f"and measured in both inputs ({measured_source})",
)
_report(
    "excluded for low detection", int((detected_all < min_detected).sum()),
    why=f"detected in fewer than {CONFIG['min_detected_fraction_all']:.0%} of PT structures",
)
_report(
    "excluded for measurement", int((~measured_in_both).sum()),
    why="structural zeros created by the accepted ortholog map",
)

gene_lookup_tested = {str(gene).upper(): index for index, gene in enumerate(gene_names)}

# --------------------------------------------------------------------------------------
# Per-specimen detection, mean abundance and bulk level, all on the tested universe.
# --------------------------------------------------------------------------------------

specimen_order = sorted(set(samples.tolist()))
detection_by_specimen = np.vstack([
    np.asarray((Y_genes[samples == specimen] > 0).mean(axis=0)).ravel() for specimen in specimen_order
])
abundance_by_specimen = np.vstack([
    np.asarray(Y_genes[samples == specimen].mean(axis=0)).ravel() for specimen in specimen_order
])

per_specimen_detection = pd.DataFrame(detection_by_specimen, index=specimen_order, columns=gene_names)
per_specimen_abundance = pd.DataFrame(abundance_by_specimen, index=specimen_order, columns=gene_names)

detection_mouse = per_specimen_detection.loc[mouse_samples].mean(axis=0)
detection_human = per_specimen_detection.loc[human_samples].mean(axis=0)
abundance_mouse = per_specimen_abundance.loc[mouse_samples].mean(axis=0)
abundance_human = per_specimen_abundance.loc[human_samples].mean(axis=0)

# --------------------------------------------------------------------------------------
# Flags.
# --------------------------------------------------------------------------------------

axis_basis_gene = np.array([gene in set(AXIS_BASIS_GENES) for gene in gene_names])
technical_gene = np.array([
    any(re.match(pattern, str(gene).upper()) for pattern in TECHNICAL_GENE_PATTERNS)
    for gene in gene_names
])

detected_in_both_species = (
    (detection_mouse.to_numpy() >= CONFIG["min_detection_each_species"])
    & (detection_human.to_numpy() >= CONFIG["min_detection_each_species"])
)

_report("axis-basis genes in the universe", int(axis_basis_gene.sum()), len(gene_names),
        why=f"of {len(AXIS_BASIS_GENES)} landmark-panel genes")
_report("technical genes in the universe", int(technical_gene.sum()), len(gene_names),
        why="mitochondrial / ribosomal protein")
_report("detected in both species", int((detected_in_both_species & ~axis_basis_gene & ~technical_gene).sum()),
        len(gene_names), why=f"detection >= {CONFIG['min_detection_each_species']:.0%} per species, "
                             "axis and technical genes excluded")

print(f"  axis-basis genes present as landmarks: "
      f"{[gene for gene in AXIS_BASIS_GENES if gene.upper() in gene_lookup_tested]}")


# %% [markdown]
# ## 0.4 - Curve helpers
#
# Every metric below is a row-wise operation on a `(genes x grid)` matrix of fitted curves, sharing
# one contract:
#
# * **NaN means "no support here", not zero.** A specimen's curve is NaN outside its own pseudospace
#   support, so a metric is always taken over a finite mask (`_finite_row_mask`) and reported as NaN
#   when too few shared grid points remain.
# * **Nothing extrapolates.** The phase search samples the human curve only inside its own support;
#   grid points a displacement pushes outside it are dropped from that comparison, and out-of-support
#   values never enter a correlation.
# * **Amplitude-free and level-free by construction** where the metric claims to be: `_row_zscore`
#   removes the row mean and the row SD before any correlation or RMS difference.
#
# Positional summaries come in two flavours on purpose: a **weighted centroid** (stable - it uses the
# whole curve) and a **peak** (legible - it uses one grid point). The centroid is the primary reading;
# the peak is reported beside it, never instead of it.
#

# %%
# Purpose: row-wise curve helpers shared by every metric in this notebook.

import warnings


def _finite_row_mask(*curves):
    """Grid points where every supplied row is finite (i.e. supported by all of them)."""
    mask = np.ones(np.asarray(curves[0], dtype=float).shape, dtype=bool)
    for curve in curves:
        mask &= np.isfinite(np.asarray(curve, dtype=float))
    return mask


def _nanmean_safe(values, axis=1, keepdims=False):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanmean(values, axis=axis, keepdims=keepdims)


def _nansum_safe(values, axis=1):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nansum(values, axis=axis)


def _row_pearson(left, right, mask=None, min_points=8):
    """Row-wise Pearson correlation over the shared finite points, else NaN."""
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    if mask is None:
        mask = _finite_row_mask(left, right)

    a = np.where(mask, left, np.nan)
    b = np.where(mask, right, np.nan)
    da = a - _nanmean_safe(a, keepdims=True)
    db = b - _nanmean_safe(b, keepdims=True)

    numerator = _nansum_safe(da * db)
    denominator = np.sqrt(_nansum_safe(da ** 2) * _nansum_safe(db ** 2))
    enough = mask.sum(axis=1) >= min_points

    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(enough & (denominator > 0), numerator / denominator, np.nan)


def _row_spearman(left, right, mask=None, min_points=8):
    """Row-wise Spearman correlation. Rows with any unsupported point are reported as NaN."""
    from scipy.stats import rankdata

    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    if mask is None:
        mask = _finite_row_mask(left, right)

    complete = mask.all(axis=1) & (mask.sum(axis=1) >= min_points)
    out = np.full(left.shape[0], np.nan)
    if complete.any():
        out[complete] = _row_pearson(
            rankdata(left[complete], axis=1), rankdata(right[complete], axis=1)
        )
    return out


def _row_zscore(curve, mask=None, min_points=8):
    """Standardise each row over its supported points: mean 0, SD 1, NaN where the row is flat."""
    curve = np.asarray(curve, dtype=float)
    if mask is None:
        mask = np.isfinite(curve)

    values = np.where(mask, curve, np.nan)
    mean = _nanmean_safe(values, keepdims=True)
    spread = np.sqrt(_nanmean_safe((values - mean) ** 2, keepdims=True))

    with np.errstate(invalid="ignore", divide="ignore"):
        z = np.where(spread > 0, (values - mean) / spread, np.nan)
    z[mask.sum(axis=1) < min_points] = np.nan
    return z


def _row_amplitude(curve, mask=None):
    """Peak-to-trough of each row over its supported points."""
    curve = np.asarray(curve, dtype=float)
    if mask is None:
        mask = np.isfinite(curve)

    values = np.where(mask, curve, np.nan)
    return np.nanmax(values, axis=1) - np.nanmin(values, axis=1)


def _row_weighted_centroid(curve, x, mask=None):
    """Positional centroid weighted by the row above its own minimum (stable, uses every point)."""
    curve = np.asarray(curve, dtype=float)
    if mask is None:
        mask = np.isfinite(curve)

    values = np.where(mask, curve, np.nan)
    weights = np.clip(values - np.nanmin(values, axis=1, keepdims=True), 0.0, None)

    numerator = _nansum_safe(weights * np.asarray(x, dtype=float))
    denominator = _nansum_safe(weights)

    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(denominator > 0, numerator / denominator, np.nan)


def _row_peak_position(curve, x, mask=None):
    """Grid position of the row maximum (legible, single-point; reported beside the centroid)."""
    curve = np.asarray(curve, dtype=float)
    if mask is None:
        mask = np.isfinite(curve)

    index = np.argmax(np.where(mask, curve, -np.inf), axis=1)
    return np.where(mask.any(axis=1), np.asarray(x, dtype=float)[index], np.nan)


def _row_early_to_late(curve, x, mask=None, quantile=0.25):
    """Late-territory mean minus early-territory mean: positive = rises from early to late PT."""
    curve = np.asarray(curve, dtype=float)
    x = np.asarray(x, dtype=float)
    if mask is None:
        mask = np.isfinite(curve)

    early = x <= np.quantile(x, quantile)
    late = x >= np.quantile(x, 1.0 - quantile)

    early_mean = _nanmean_safe(np.where(mask & early, curve, np.nan))
    late_mean = _nanmean_safe(np.where(mask & late, curve, np.nan))
    return late_mean - early_mean


def _row_halfmax_window(curve, x, mask=None):
    """Onset, offset and width of the stretch above the row's half-maximum (fraction of PT)."""
    curve = np.asarray(curve, dtype=float)
    x = np.asarray(x, dtype=float)
    if mask is None:
        mask = np.isfinite(curve)

    values = np.where(mask, curve, np.nan)
    low = np.nanmin(values, axis=1, keepdims=True)
    high = np.nanmax(values, axis=1, keepdims=True)
    above = mask & (values >= low + 0.5 * (high - low))

    onset = np.where(above.any(axis=1), x[np.argmax(above, axis=1)], np.nan)
    offset = np.where(above.any(axis=1), x[x.size - 1 - np.argmax(above[:, ::-1], axis=1)], np.nan)
    return onset, offset, offset - onset


def _interp_rows(values, source_x, target_x):
    """Resample each row onto `target_x` using only its supported points.

    Interpolation runs over the row's finite points alone, so a target position outside the row's own
    finite range is NaN instead of an invented value, and an unsupported stretch contributes nothing
    to a later comparison. An interior gap between supported points is bridged linearly.
    """
    values = np.asarray(values, dtype=float)
    source_x = np.asarray(source_x, dtype=float)
    target_x = np.asarray(target_x, dtype=float)

    out = np.full((values.shape[0], target_x.size), np.nan, dtype=float)
    for row in range(values.shape[0]):
        supported = np.flatnonzero(np.isfinite(values[row]))
        if supported.size < 2:
            continue
        out[row] = np.interp(
            target_x, source_x[supported], values[row][supported], left=np.nan, right=np.nan
        )
    return out


def _row_shift_search(mouse_z, human_z, x_unit, max_shift, step, min_points=8):
    """Constrained horizontal registration of each human z-curve onto its mouse z-curve.

    Positive shift = the human positional program lies **later** along PT than the mouse program, so
    the human curve is sampled at ``x + shift`` and compared with the mouse curve at ``x``. The
    search is capped at ``max_shift`` because beyond that "shifted" stops being a description of an
    aligned gradient. Returns (best shift, best correlation, standardised RMS after the best shift);
    the delta=0 correlation is the curve correlation the caller computes separately.
    """
    mouse_z = np.asarray(mouse_z, dtype=float)
    human_z = np.asarray(human_z, dtype=float)
    x_unit = np.asarray(x_unit, dtype=float)

    shifts = np.arange(-max_shift, max_shift + step / 2.0, step)
    n_genes = mouse_z.shape[0]

    best_shift = np.full(n_genes, np.nan)
    best_corr = np.full(n_genes, np.nan)

    for shift in shifts:
        moved = _interp_rows(human_z, x_unit, x_unit + shift)
        corr = _row_pearson(mouse_z, moved, min_points=min_points)
        improves = np.isfinite(corr) & (~np.isfinite(best_corr) | (corr > best_corr))
        best_corr = np.where(improves, corr, best_corr)
        best_shift = np.where(improves, shift, best_shift)

    residual_rms = np.full(n_genes, np.nan)
    for shift in shifts:
        rows = np.flatnonzero(np.isclose(best_shift, shift))
        if rows.size == 0:
            continue
        moved = _interp_rows(human_z[rows], x_unit, x_unit + shift)
        mask = _finite_row_mask(mouse_z[rows], moved)
        difference = np.where(mask, mouse_z[rows] - moved, np.nan)
        residual_rms[rows] = np.sqrt(_nanmean_safe(difference ** 2))

    return best_shift, best_corr, residual_rms



# %% [markdown]
# ## 1 - Balanced species curves
#
# **Why this metric.** Every cross-species statement in this notebook is a statement about the pair
# `(balanced mouse curve, balanced human curve)` for one gene, so how those two curves are built is the
# first thing that has to be defensible.
#
# **What is done.** The pooled nested fit from 03 is reused (it is cached: same stage name, grid,
# knots, lambda grid, `s`, `c` and Y fingerprint, and `cached_run_level_shape` keys on the fitter's
# source, so this is a cache hit rather than a recomputation). Each specimen is then fitted
# separately with the shared-smooth model at the **pooled per-gene smoothing parameter**, and the
# species curve is the **equal-weight mean of its specimens' curves** - `specimen_balanced_curves`.
#
# **Why not weight by tubules.** A pooled fit weights every tubule equally, so a specimen that
# contributes three times as many structures shapes the curve three times as much. With two mice and
# one donor that is not a nuisance parameter, it is the whole comparison.
#
# **What counts as evidence here.** Nothing yet: this section produces the curves, the per-specimen
# agreement that later sections gate on, and records for each gene how many grid points both species
# actually support. A gene whose balanced curves overlap on too few points cannot support any
# positional reading, however interesting its mean difference is.
#
# **Limits stated up front.** A specimen's curve is NaN outside its own pseudospace support, so the
# balanced curves are read over the shared finite grid only (`n_shared_grid_points`). The per-gene
# smoothing parameter comes from the pooled fit and is reused per specimen, so between-specimen spread
# reflects biology rather than different amounts of smoothing.
#

# %%
# Purpose: pooled fit (shared cache), per-specimen fits (06's cache stage), balanced species curves.

c = (species == "human").astype(float)

gene_fit = cached_run_level_shape(
    Y_genes,
    s,
    c,
    knots,
    grid,
    CONFIG["lambda_grid"],
    stage="section4_gene_fit",          # 03/05's stage name: keyed on the fitter, so this is a hit
    root=STAGE_CACHE_DIR,
    y_fingerprint=digest(Y_genes),
    enabled=STAGE_CACHE_ENABLED,
)

pooled_mouse = gene_fit["curve_healthy"]     # c == 0
pooled_human = gene_fit["curve_aki"]         # c == 1

print(f"pooled fit: mouse {pooled_mouse.shape}, human {pooled_human.shape}, "
      f"{len(np.unique(gene_fit['lam_idx']))} distinct smoothing parameters")


def _fit_specimen_curves():
    """Shared-smooth fit per specimen at the pooled per-gene lambda (04/03 convention)."""
    fitted = {}
    for specimen in specimen_order:
        selected = samples == specimen
        curves, _ = fit_single_condition_curves(
            Y_genes[selected],
            s[selected],
            knots,
            grid,
            CONFIG["lambda_grid"],
            gene_fit["lam_idx"],
            support_pct=CONFIG["common_support_pct"],
        )
        fitted[f"curve__{specimen}"] = curves
    return fitted


specimen_payload = cached_payload(
    "spatial_rewiring_specimen_curves",
    _fit_specimen_curves,
    root=STAGE_CACHE_DIR,
    params={
        "grid": grid,
        "knots": knots,
        "lambda_grid": CONFIG["lambda_grid"],
        "support_pct": np.asarray(CONFIG["common_support_pct"], dtype=float),
        "specimens": np.asarray(specimen_order),
        "lam_idx": gene_fit["lam_idx"],
    },
    inputs={"s": s, "samples": samples, "y": digest(Y_genes)},
    code=fit_single_condition_curves,
    enabled=STAGE_CACHE_ENABLED,
)

specimen_curves = {name.split("__", 1)[1]: specimen_payload[name] for name in specimen_payload}

balanced_mouse = specimen_balanced_curves({name: specimen_curves[name] for name in mouse_samples})
balanced_human = specimen_balanced_curves({name: specimen_curves[name] for name in human_samples})

np.savez_compressed(
    OUT_DIR / "balanced_gene_curves.npz",
    genes=gene_names, grid=grid, grid_unit=grid_unit,
    mouse=balanced_mouse, human=balanced_human,
)
np.savez_compressed(
    OUT_DIR / "specimen_gene_curves.npz",
    genes=gene_names, grid=grid, grid_unit=grid_unit,
    **{f"curve__{name}": specimen_curves[name] for name in specimen_order},
)

print(f"saved {(OUT_DIR / 'balanced_gene_curves.npz').relative_to(PROJECT_DIR)}")
print(f"saved {(OUT_DIR / 'specimen_gene_curves.npz').relative_to(PROJECT_DIR)}")

# --------------------------------------------------------------------------------------
# How much of the grid each specimen and each balanced curve actually supports.
# --------------------------------------------------------------------------------------

support_report = pd.DataFrame({
    "specimen": specimen_order,
    "species": [str(sample_species.set_index("sample").loc[name, "comparison_species"])
                for name in specimen_order],
    "median_supported_grid_points": [int(np.median(np.isfinite(specimen_curves[name]).sum(axis=1)))
                                     for name in specimen_order],
    "median_amplitude": [float(np.nanmedian(_row_amplitude(specimen_curves[name])))
                         for name in specimen_order],
})
display(support_report)

mouse_support = np.isfinite(balanced_mouse)
human_support = np.isfinite(balanced_human)
shared_support = mouse_support & human_support

_report("genes with a finite balanced mouse curve", int(mouse_support.any(axis=1).sum()), len(gene_names))
_report("genes with a finite balanced human curve", int(human_support.any(axis=1).sum()), len(gene_names))
_report("genes over the minimum shared grid points", int((shared_support.sum(axis=1) >= CONFIG["min_pattern_grid_points"]).sum()),
        len(gene_names), why=f">= {CONFIG['min_pattern_grid_points']} points supported by both species")


# %%
# Purpose: within-species reproducibility, per-specimen amplitude, and the cross-check against 05.

from itertools import combinations


def _within_species_reproducibility(specimens):
    """Mean pairwise correlation of the specimen curves of one species (NaN with <2 specimens)."""
    pairs = list(combinations(sorted(specimens), 2))
    if not pairs:
        return np.full(len(gene_names), np.nan)
    correlations = np.vstack([
        _row_pearson(specimen_curves[first], specimen_curves[second]) for first, second in pairs
    ])
    return _nanmean_safe(correlations, axis=0)


mouse_reproducibility = _within_species_reproducibility(mouse_samples)
human_reproducibility = _within_species_reproducibility(human_samples)

amplitude_by_specimen = pd.DataFrame(
    {name: _row_amplitude(specimen_curves[name]) for name in specimen_order},
    index=gene_names,
)

common_balanced = mouse_support & human_support

curve_table = pd.DataFrame({
    "gene": gene_names,
    "n_shared_grid_points": shared_support.sum(axis=1),
    "level_effect_balanced_human_minus_mouse":
        _nanmean_safe(np.where(common_balanced, balanced_human - balanced_mouse, np.nan)),
    "mouse_amplitude_balanced": _row_amplitude(balanced_mouse, common_balanced),
    "human_amplitude_balanced": _row_amplitude(balanced_human, common_balanced),
    "mouse_reproducibility": mouse_reproducibility,
    "human_reproducibility": human_reproducibility,
    "detection_mouse": np.asarray(detection_mouse, dtype=float),
    "detection_human": np.asarray(detection_human, dtype=float),
    "abundance_mouse": np.asarray(abundance_mouse, dtype=float),
    "abundance_human": np.asarray(abundance_human, dtype=float),
    "detected_in_both_species": detected_in_both_species,
    "axis_basis_gene": axis_basis_gene,
    "technical_gene": technical_gene,
})

curve_table["amplitude_log2_ratio_balanced"] = np.where(
    (curve_table["mouse_amplitude_balanced"] > 0) & (curve_table["human_amplitude_balanced"] > 0),
    np.log2(curve_table["human_amplitude_balanced"].to_numpy())
    - np.log2(curve_table["mouse_amplitude_balanced"].to_numpy()),
    np.nan,
)

for name in specimen_order:
    curve_table[f"amplitude__{name}"] = amplitude_by_specimen[name].to_numpy()
    curve_table[f"balanced_corr__{name}"] = np.where(
        str(sample_species.set_index("sample").loc[name, "comparison_species"]) == "mouse",
        _row_pearson(specimen_curves[name], balanced_mouse),
        _row_pearson(specimen_curves[name], balanced_human),
    )

curve_table["median_specimen_amplitude_mouse"] = amplitude_by_specimen[mouse_samples].median(axis=1).to_numpy()
curve_table["median_specimen_amplitude_human"] = amplitude_by_specimen[human_samples].median(axis=1).to_numpy()

_report("median within-mouse reproducibility", float(np.nanmedian(mouse_reproducibility)).__round__(3),
        why="correlation of the two mouse specimen curves")
_report("median within-human reproducibility", float(np.nanmedian(human_reproducibility)).__round__(3),
        why="correlation of the two human slice curves (one donor)")

# --------------------------------------------------------------------------------------
# Cross-check: 05's saved balanced curves, recomputed here from the same object.
#
# This is a check, not an input: 05 wrote these from its own run of the same fitter. Agreement
# means this notebook's curve layer is reproducible from the saved object; disagreement would have
# to be explained before any metric below is believed.
# --------------------------------------------------------------------------------------

cross_check_rows = []

if PREVIOUS_CURVES_PATH.exists():
    with np.load(PREVIOUS_CURVES_PATH, allow_pickle=False) as stored:
        previous_genes = stored["genes"].astype(str)
        previous_mouse = stored["mouse"]
        previous_human = stored["human"]
        previous_grid = stored["grid"]

    shared_genes, previous_index, current_index = np.intersect1d(
        previous_genes, gene_names, return_indices=True
    )

    same_grid = previous_grid.shape == grid.shape and np.allclose(previous_grid, grid)

    for label, previous, current in (
        ("balanced mouse", previous_mouse, balanced_mouse),
        ("balanced human", previous_human, balanced_human),
    ):
        both_finite = np.isfinite(previous[previous_index]) & np.isfinite(current[current_index])
        difference = np.abs(previous[previous_index] - current[current_index])[both_finite]
        cross_check_rows.append({
            "object": label,
            "n_shared_genes": shared_genes.size,
            "grid_identical": bool(same_grid),
            "median_abs_difference": float(np.median(difference)) if difference.size else np.nan,
            "p99_abs_difference": float(np.percentile(difference, 99)) if difference.size else np.nan,
        })

    cross_check = pd.DataFrame(cross_check_rows)
    display(cross_check.round(6))
else:
    cross_check = pd.DataFrame()
    print(f"{PREVIOUS_CURVES_PATH} not found - cross-check against 05 skipped "
          "(run 05_human_vs_healthy_mouse.ipynb to produce it)")

if PREVIOUS_ATLAS_PATH.exists():
    previous_atlas = pd.read_csv(PREVIOUS_ATLAS_PATH)

    # atlas column -> this notebook's balanced-curve column. Names differ because the primary
    # columns of the atlas written in section 13 are the *registered* metrics, not these.
    atlas_pairs = {
        "mouse_amplitude": "mouse_amplitude_balanced",
        "human_amplitude": "human_amplitude_balanced",
        "level_effect_human_minus_mouse": "level_effect_balanced_human_minus_mouse",
    }

    merged = curve_table.merge(previous_atlas[["gene"] + list(atlas_pairs)], on="gene")
    atlas_comparison = pd.DataFrame([
        {
            "column": atlas_column,
            "our_median": float(np.nanmedian(merged[our_column])),
            "atlas05_median": float(np.nanmedian(merged[atlas_column])),
            "spearman": float(spearmanr(merged[our_column], merged[atlas_column],
                                        nan_policy="omit").correlation),
        }
        for atlas_column, our_column in atlas_pairs.items()
    ])
    display(atlas_comparison.round(4))
else:
    print(f"{PREVIOUS_ATLAS_PATH} not found - atlas comparison skipped")


# %% [markdown]
# ## 2 - Cross-species registration of PT position
#
# **Why this step exists.** The mouse and human curves are both functions of *this object's* DPT
# coordinate. That coordinate is oriented by conserved PT markers in each species separately, so equal
# numerical position in the two species need not mean the same anatomical position. Reading a
# horizontal displacement off the unregistered axis would then be reading a coordinate artefact as
# biology. So a monotone human -> mouse mapping is established from landmark programs first, and the
# main analysis runs on the registered coordinate. The unregistered axis is kept as a sensitivity
# analysis, and a gene is only called *shifted* if the effect survives both.
#
# **How the mapping is built.** Canonical PT landmark programs (`LANDMARK_PANELS`: S1, S2, S3, taken
# from 03's own marker definitions) are scored per species as the mean of their members' within-species
# z-scored balanced curves. A landmark is used only if at least two of its members are detected in
# `>= 5%` of both species' structures - a marker that is barely expressed in one species cannot anchor
# anything. Each used landmark contributes its weighted centroid in the two species as one anchor;
# monotone piecewise-linear interpolation through the anchors, plus the endpoints `(0, 0)` and `(1, 1)`,
# is the mapping. Anchor positions are sorted along the mouse axis and forced monotone, so the mapping
# is monotone by construction.
#
# **What counts as evidence.** For the mapping itself: a small number of landmarks whose human
# positions are ordered consistently with their mouse positions. For the rest of the notebook: every
# "shift" claim has to survive the *unregistered* comparison too, which is exactly what section 5
# tests. Registration is deliberately coarse: three programs, not the whole transcriptome, so that it
# cannot smooth away the rewiring we want to measure.
#
# **Limits.** With S1/S2/S3 the mapping is anchored by at most three points and is only as good as
# those programs; residual misalignment inside a segment is not corrected. Landmark genes are flagged
# `axis_basis_gene` and are excluded from unbiased discovery below, because a program used to align
# the two species cannot also be presented as a finding about them.
#

# %%
# Purpose: landmark programs, their positions in each species, and the usable-landmark audit.

detection_mouse_array = np.asarray(detection_mouse, dtype=float)
detection_human_array = np.asarray(detection_human, dtype=float)


def _program_curve(curves, member_index, x_unit):
    """Landmark program curve: mean of the members' within-species z-scored curves."""
    z_scored = _row_zscore(curves[member_index])
    return _nanmean_safe(z_scored, axis=0)


def _build_landmarks(mouse_curves, human_curves, verbose=True):
    """Score every landmark panel in both species and report which panels can anchor the mapping.

    A panel is usable when at least two members are detected in >= `min_detection_each_species` of
    both species, and its own program curve reaches `amplitude_patterned` in both species.
    """
    program_rows = []
    gene_rows = []
    programs = {}

    for panel, requested in LANDMARK_PANELS.items():
        present_index = []
        members = []

        for gene in requested:
            index = gene_lookup_tested.get(gene.upper())
            if index is None:
                gene_rows.append({"panel": panel, "gene": gene, "present": False,
                                  "detection_mouse": np.nan, "detection_human": np.nan,
                                  "usable": False})
                continue
            detected_both = (
                detection_mouse_array[index] >= CONFIG["min_detection_each_species"]
                and detection_human_array[index] >= CONFIG["min_detection_each_species"]
            )
            gene_rows.append({
                "panel": panel, "gene": gene, "present": True,
                "detection_mouse": float(detection_mouse_array[index]),
                "detection_human": float(detection_human_array[index]),
                "usable": bool(detected_both),
            })
            if detected_both:
                present_index.append(index)
                members.append(gene)

        if len(present_index) < 2:
            program_rows.append({
                "panel": panel, "species": "both", "n_members_used": len(present_index),
                "members": ";".join(members), "mouse_centroid": np.nan, "human_centroid": np.nan,
                "mouse_peak": np.nan, "human_peak": np.nan, "mouse_amplitude": np.nan,
                "human_amplitude": np.nan, "used": False,
            })
            continue

        member_index = np.asarray(present_index, dtype=int)
        mouse_program = _program_curve(mouse_curves, member_index, grid_unit)
        human_program = _program_curve(human_curves, member_index, grid_unit)
        programs[panel] = {"mouse": mouse_program, "human": human_program, "members": members}

        mouse_amplitude = float(_row_amplitude(mouse_program[None, :]))
        human_amplitude = float(_row_amplitude(human_program[None, :]))
        used = bool(
            mouse_amplitude >= CONFIG["amplitude_patterned"]
            and human_amplitude >= CONFIG["amplitude_patterned"]
        )

        program_rows.append({
            "panel": panel, "species": "both", "n_members_used": len(present_index),
            "members": ";".join(members),
            "mouse_centroid": float(_row_weighted_centroid(mouse_program[None, :], grid_unit)),
            "human_centroid": float(_row_weighted_centroid(human_program[None, :], grid_unit)),
            "mouse_peak": float(_row_peak_position(mouse_program[None, :], grid_unit)),
            "human_peak": float(_row_peak_position(human_program[None, :], grid_unit)),
            "mouse_amplitude": mouse_amplitude, "human_amplitude": human_amplitude,
            "used": used,
        })

    landmark_table = pd.DataFrame(program_rows)
    landmark_genes = pd.DataFrame(gene_rows)

    if verbose:
        _report("usable landmark panels", int(landmark_table["used"].sum()), len(landmark_table),
                why=">= 2 members detected in both species and an amplitude in both")
        display(landmark_table.round(4))
        print("  landmark gene audit:")
        display(landmark_genes)

    return programs, landmark_table, landmark_genes


landmark_programs, landmark_table, landmark_genes = _build_landmarks(balanced_mouse, balanced_human)

_save_table(landmark_table, "registration_landmarks.csv")
print(f"saved {(OUT_DIR / 'registration_landmarks.csv').relative_to(PROJECT_DIR)}")


# %%
# Purpose: the monotone landmark mapping, the registered human curves, and the registration figure.


def _register(mouse_curves, human_curves, verbose=False):
    """Landmark-anchored monotone human -> mouse registration of the pseudospace coordinate.

    Returns the mapping evaluated on the unit grid (`mapped`, human position -> mouse position), the
    human curves resampled onto the mouse axis, and the support mask that resampling preserves.
    With fewer than two usable landmarks the mapping falls back to the identity and says so.
    """
    programs, table, genes = _build_landmarks(mouse_curves, human_curves, verbose=False)
    used = table[table["used"]].sort_values("mouse_centroid")

    if len(used) < 2:
        mapped = grid_unit.copy()
        anchors = pd.DataFrame(columns=["panel", "human_centroid", "mouse_centroid"])
    else:
        human_anchors = np.concatenate([[0.0], used["human_centroid"].to_numpy(), [1.0]])
        mouse_anchors = np.concatenate([[0.0], used["mouse_centroid"].to_numpy(), [1.0]])
        # Monotone by construction: the human anchors are forced non-decreasing along the mouse axis.
        human_anchors = np.maximum.accumulate(human_anchors)
        mapped = np.interp(grid_unit, human_anchors, mouse_anchors)
        anchors = used[["panel", "human_centroid", "mouse_centroid"]].reset_index(drop=True)

    human_registered = _interp_rows(human_curves, mapped, grid_unit)
    support_indicator = _interp_rows(
        np.isfinite(human_curves).astype(float), mapped, grid_unit
    )
    human_registered_mask = support_indicator >= 0.999

    if verbose:
        print(f"  landmarks used: {len(used)} ({', '.join(anchors['panel']) if len(anchors) else 'none'})")
        print(f"  max |mapped - identity|: {np.max(np.abs(mapped - grid_unit)):.4f} of PT")

    return {
        "mapped": mapped,
        "human_registered": human_registered,
        "human_registered_mask": human_registered_mask,
        "anchors": anchors,
        "n_landmarks": int(len(used)),
        "programs": programs,
        "landmark_table": table,
        "landmark_genes": genes,
    }


primary_registration = _register(balanced_mouse, balanced_human, verbose=True)

human_registered = primary_registration["human_registered"]
human_registered_mask = primary_registration["human_registered_mask"]
registration_shift = primary_registration["mapped"] - grid_unit

registration_summary = pd.DataFrame({
    "quantity": [
        "landmarks_used",
        "max_mapping_deviation_from_identity",
        "median_mapping_deviation_from_identity",
        "grid_points_supported_before",
        "grid_points_supported_after",
    ],
    "value": [
        primary_registration["n_landmarks"],
        float(np.max(np.abs(registration_shift))),
        float(np.median(np.abs(registration_shift))),
        int(np.isfinite(balanced_human).sum(axis=1).mean()),
        int(human_registered_mask.sum(axis=1).mean()),
    ],
})
display(registration_summary)

# --------------------------------------------------------------------------------------
# Figure: canonical programs before registration, landmark positions, the mapping, and the
# canonical gradient genes after registration.
# --------------------------------------------------------------------------------------

validation_genes = ["Slc5a2", "Slc22a6", "Slc22a7"]

figure, axes = plt.subplots(2, 2, figsize=(12, 7.6))
(ax_before, ax_positions), (ax_mapping, ax_after) = axes

for panel, curves in primary_registration["programs"].items():
    ax_before.plot(grid_unit, curves["mouse"], color=SPECIES_COLORS["mouse"], lw=1.8,
                   label=f"{panel} (mouse)" if panel == next(iter(primary_registration["programs"])) else None)
    ax_before.plot(grid_unit, curves["human"], color=SPECIES_COLORS["human"], lw=1.8, ls="--",
                   label=f"{panel} (human)" if panel == next(iter(primary_registration["programs"])) else None)
ax_before.set_title("Landmark programs before registration\n(solid = mouse, dashed = human)", loc="left")
ax_before.set_xlabel("PT position (unit): early -> late")
ax_before.set_ylabel("Mean member z-score")
ax_before.legend()

ax_positions.plot([0, 1], [0, 1], color="0.7", lw=1, ls=":", label="identity")
ax_positions.scatter(landmark_table["mouse_centroid"], landmark_table["human_centroid"],
                     s=60, color="#444444", zorder=3, label="landmark centroid")
for row in landmark_table.itertuples():
    ax_positions.annotate(row.panel, (row.mouse_centroid, row.human_centroid),
                          xytext=(4, 4), textcoords="offset points", fontsize=8)
ax_positions.set_title("Landmark positions: consistent ordering = usable mapping", loc="left")
ax_positions.set_xlabel("Mouse landmark centroid")
ax_positions.set_ylabel("Human landmark centroid")
ax_positions.legend()

ax_mapping.plot(grid_unit, grid_unit, color="0.7", lw=1, ls=":", label="identity")
ax_mapping.plot(grid_unit, primary_registration["mapped"], color="#444444", lw=2.0,
                label="mapping m(human) -> mouse")
if len(primary_registration["anchors"]):
    ax_mapping.scatter(primary_registration["anchors"]["human_centroid"],
                       primary_registration["anchors"]["mouse_centroid"], s=45, color="#D55E00", zorder=3)
ax_mapping.set_title("Monotone human -> mouse mapping", loc="left")
ax_mapping.set_xlabel("Human PT position")
ax_mapping.set_ylabel("Registered mouse PT position")
ax_mapping.legend()

for gene in validation_genes:
    index = gene_lookup_tested.get(gene.upper())
    if index is None:
        continue
    ax_after.plot(grid_unit, balanced_mouse[index], color=SPECIES_COLORS["mouse"], lw=1.8,
                  label="mouse" if gene == validation_genes[0] else None)
    ax_after.plot(grid_unit, human_registered[index], color=SPECIES_COLORS["human"], lw=1.8, ls="--",
                  label="human (registered)" if gene == validation_genes[0] else None)
    ax_after.annotate(gene, (grid_unit[int(np.nanargmax(balanced_mouse[index]))],
                             float(np.nanmax(balanced_mouse[index]))),
                      xytext=(2, 2), textcoords="offset points", fontsize=8)
ax_after.set_title("Canonical PT gradients on the registered axis", loc="left")
ax_after.set_xlabel("PT position (unit, registered)")
ax_after.set_ylabel("Fitted lognorm")
ax_after.legend()

figure.suptitle("Cross-species registration of PT position from S1/S2/S3 landmark programs\n"
                "descriptive alignment for this cohort (2 mouse specimens, 1 human donor)",
                fontsize=12)
figure.tight_layout()
_save_figure(figure, "pt_cross_species_registration.png")
plt.show()


# %% [markdown]
# ## 3 - Gene-level curve metrics
#
# Six quantities per gene, each answering a different question, in the order the decomposition is
# written (level + amplitude + phase + residual):
#
# | | metric | question | invariant to |
# | --- | --- | --- | --- |
# | A | `level_effect_human_minus_mouse` | is the gene's overall abundance different? | nothing (this is the conventional DE-like component) |
# | B | `mouse_amplitude`, `human_amplitude`, `amplitude_log2_ratio_human_over_mouse` | how strong is each species' zonation? | level |
# | C | `shape_corr` (standardised), `pattern_rms_z` | do the two curves have the same spatial pattern? | level **and** amplitude |
# | D | `mouse_position_centroid`, `human_position_centroid`, peak, early-to-late gradient, half-max width | where does expression sit along PT in each species? | level, amplitude |
# | E | `best_shift_human_minus_mouse`, `registered_shape_corr`, `shift_improvement`, `residual_rms_after_shift` | is a modified *position* enough to explain the difference? | level, amplitude |
# | F | `mouse_early_to_late`, `human_early_to_late` | do the two species run in opposite directions? | level, amplitude |
#
# **How to read the phase metric.** `best_shift` is the horizontal displacement, searched only up to
# `max_registration_shift` (20% of PT) and only in the direction that *improves* agreement, that best
# aligns the human curve to the mouse curve after both are standardised. Positive = the human program
# sits later along PT. Its trustworthiness is not the shift value alone: a high `registered_shape_corr`
# *and* a real `shift_improvement` over the unshifted correlation are both required before section 4
# will call a gene shifted, and section 5 additionally asks whether the same shift survives on the
# unregistered axis.
#
# **Primary axis and sensitivity axis.** Every metric is computed twice: on the registered coordinate
# (the section-2 mapping applied to the human curves) and on the unregistered shared coordinate. The
# registered values are the primary ones; the `*_unregistered` columns exist to test them.
#
# **What would count as evidence.** Only the combination: reproducibility within each species, enough
# shared support, adequate amplitude in both species, and — for phase claims — the shift surviving
# both coordinate choices. A gene that is flat or noisy in one species cannot be called displaced,
# however large its level difference.
#

# %%
# Purpose: level, amplitude, shape, position, phase and inversion metrics, registered and unregistered.

mouse_mask = np.isfinite(balanced_mouse)
human_mask_unregistered = np.isfinite(balanced_human)
human_mask_registered = human_registered_mask


def _curve_metrics(mouse, mouse_mask_values, human, human_mask_values, x_unit, suffix=""):
    """Every curve metric for one (mouse, human) balanced-curve pair, on one coordinate choice."""
    common = mouse_mask_values & human_mask_values
    n_shared = common.sum(axis=1)

    level = _nanmean_safe(np.where(common, human - mouse, np.nan))
    mouse_amplitude = _row_amplitude(mouse, common)
    human_amplitude = _row_amplitude(human, common)

    z_mouse = _row_zscore(mouse, common)
    z_human = _row_zscore(human, common)

    with np.errstate(invalid="ignore", divide="ignore"):
        amplitude_difference = human_amplitude - mouse_amplitude
        amplitude_log2_ratio = np.where(
            (mouse_amplitude > 0) & (human_amplitude > 0),
            np.log2(human_amplitude) - np.log2(mouse_amplitude),
            np.nan,
        )

    # C - level- and amplitude-free: both curves are standardised over the shared points first.
    shape_corr = _row_pearson(z_mouse, z_human, mask=common)
    # Rank-based counterpart. It needs every shared grid point to be supported in both species
    # (ranks are meaningless across an unsupported gap), so it is NaN for edge-truncated genes and
    # is reported as a diagnostic beside the Pearson correlation, never instead of it.
    shape_spearman = _row_spearman(z_mouse, z_human, mask=common)
    pattern_rms_z = np.sqrt(_nanmean_safe(np.where(common, z_human - z_mouse, np.nan) ** 2))

    # D - where expression sits, in each species.
    mouse_centroid = _row_weighted_centroid(mouse, x_unit, common)
    human_centroid = _row_weighted_centroid(human, x_unit, common)
    mouse_peak = _row_peak_position(mouse, x_unit, common)
    human_peak = _row_peak_position(human, x_unit, common)
    mouse_early_to_late = _row_early_to_late(mouse, x_unit, common)
    human_early_to_late = _row_early_to_late(human, x_unit, common)
    mouse_onset, mouse_offset, mouse_width = _row_halfmax_window(mouse, x_unit, common)
    human_onset, human_offset, human_width = _row_halfmax_window(human, x_unit, common)

    # E - the constrained phase search.
    best_shift, best_corr, residual_rms = _row_shift_search(
        z_mouse, z_human, x_unit,
        CONFIG["max_registration_shift"], CONFIG["shift_step"],
        min_points=CONFIG["min_pattern_grid_points"],
    )

    return pd.DataFrame({
        f"n_shared_grid_points{suffix}": n_shared,
        f"enough_shared_support{suffix}": n_shared >= CONFIG["min_pattern_grid_points"],
        f"level_effect_human_minus_mouse{suffix}": level,
        f"mouse_amplitude{suffix}": mouse_amplitude,
        f"human_amplitude{suffix}": human_amplitude,
        f"amplitude_difference_human_minus_mouse{suffix}": amplitude_difference,
        f"amplitude_log2_ratio_human_over_mouse{suffix}": amplitude_log2_ratio,
        f"shape_corr{suffix}": shape_corr,
        f"shape_spearman{suffix}": shape_spearman,
        f"pattern_rms_z{suffix}": pattern_rms_z,
        f"best_shift_human_minus_mouse{suffix}": best_shift,
        f"registered_shape_corr{suffix}": best_corr,
        f"shift_improvement{suffix}": best_corr - shape_corr,
        f"residual_rms_after_shift{suffix}": residual_rms,
        f"mouse_position_centroid{suffix}": mouse_centroid,
        f"human_position_centroid{suffix}": human_centroid,
        f"centroid_shift_human_minus_mouse{suffix}": human_centroid - mouse_centroid,
        f"mouse_peak_position{suffix}": mouse_peak,
        f"human_peak_position{suffix}": human_peak,
        f"peak_shift_human_minus_mouse{suffix}": human_peak - mouse_peak,
        f"mouse_early_to_late{suffix}": mouse_early_to_late,
        f"human_early_to_late{suffix}": human_early_to_late,
        f"gradient_change_human_minus_mouse{suffix}": human_early_to_late - mouse_early_to_late,
        f"mouse_halfmax_onset{suffix}": mouse_onset,
        f"mouse_halfmax_offset{suffix}": mouse_offset,
        f"mouse_halfmax_width{suffix}": mouse_width,
        f"human_halfmax_onset{suffix}": human_onset,
        f"human_halfmax_offset{suffix}": human_offset,
        f"human_halfmax_width{suffix}": human_width,
    })


# The unregistered axis carries the metrics a phase claim has to survive on, plus the amplitudes the
# registration is allowed to change; everything else stays in the raw metric frames.
unregistered_columns = [
    "n_shared_grid_points_unregistered",
    "level_effect_human_minus_mouse_unregistered",
    "mouse_amplitude_unregistered",
    "human_amplitude_unregistered",
    "amplitude_log2_ratio_human_over_mouse_unregistered",
    "shape_corr_unregistered",
    "shape_spearman_unregistered",
    "best_shift_human_minus_mouse_unregistered",
    "registered_shape_corr_unregistered",
    "shift_improvement_unregistered",
    "centroid_shift_human_minus_mouse_unregistered",
    "gradient_change_human_minus_mouse_unregistered",
]


def _assemble_metrics(mouse_curves, human_curves, register=True, registration=None):
    """The full per-gene metric table for one curve pair, on the registered or the raw coordinate.

    Shared by the primary analysis and by every sensitivity variant in section 5, so a variant can
    only differ through the curves, the coordinate choice, or the thresholds - never through a
    different code path.
    """
    mouse_mask_values = np.isfinite(mouse_curves)

    if register:
        used_registration = registration if registration is not None else _register(mouse_curves, human_curves)
        human_used = used_registration["human_registered"]
        human_mask_values = used_registration["human_registered_mask"]
    else:
        used_registration = None
        human_used = human_curves
        human_mask_values = np.isfinite(human_curves)

    registered_metrics = _curve_metrics(
        mouse_curves, mouse_mask_values, human_used, human_mask_values, grid_unit
    )
    unregistered_metrics = _curve_metrics(
        mouse_curves, mouse_mask_values, human_curves, np.isfinite(human_curves), grid_unit,
        suffix="_unregistered",
    )

    assembled = (
        pd.DataFrame({"gene": gene_names})
        .join(registered_metrics)
        .join(unregistered_metrics[unregistered_columns])
        .merge(curve_table, on="gene", how="left")
    )
    return assembled, used_registration


gene_metrics, primary_registration_used = _assemble_metrics(
    balanced_mouse, balanced_human, register=True, registration=primary_registration
)

_report("genes with comparable curves on the registered axis",
        int(gene_metrics["enough_shared_support"].sum()), len(gene_names),
        why=f"both species support >= {CONFIG['min_pattern_grid_points']} shared grid points")
_report("genes with a rank-based shape correlation",
        int(np.isfinite(gene_metrics["shape_spearman"]).sum()), len(gene_names),
        why="Spearman needs every shared grid point supported by both species; NaN elsewhere")
_report("genes where registration changed the shared support",
        int((gene_metrics["n_shared_grid_points"]
             != gene_metrics["n_shared_grid_points_unregistered"]).sum()), len(gene_names))

display(gene_metrics[[
    "gene", "level_effect_human_minus_mouse", "mouse_amplitude", "human_amplitude",
    "amplitude_log2_ratio_human_over_mouse", "shape_corr", "best_shift_human_minus_mouse",
    "registered_shape_corr", "centroid_shift_human_minus_mouse",
]].head(10).round(4))


# %% [markdown]
# ## 4 - Spatial phenotypes
#
# **The classes are a reading aid, not a test.** They compress the six metric families into the
# categories the biological question needs — conserved zonation, weaker/stronger zonation in human,
# mouse-zonated/human-flat, human-zonated/mouse-flat, shifted earlier/later in human, gradient
# inversion, complex shape rewiring, weak/uncertain — and every one of them is a conjunction of
# thresholds, evaluated in a fixed priority order. The component metrics stay in the atlas, so a
# reader who disagrees with a threshold can re-derive a different partition without refitting
# anything.
#
# **Priority order** (first match wins, highest priority listed last in the code so it overwrites):
# excluded (not measured in both species, technical, or axis-basis) -> insufficient shared support ->
# mouse-zonated/human-flat -> human-zonated/mouse-flat -> gradient inversion -> shifted
# earlier/later -> conserved zonation -> weaker/stronger zonation in human -> complex shape rewiring ->
# weak/uncertain.
#
# **How the thresholds were chosen.** The next cell reports the empirical distributions the thresholds
# act on (amplitude, within-species reproducibility, shape correlation, displacement) and how many
# genes each one selects. `CONFIG` holds the values; the distributions are printed so the choice is
# auditable rather than inherited. They start from 03/05's operating values and were then checked
# against these distributions.
#
# **Guardrails encoded here.** A gene is never called spatially rewired when one species is
# essentially flat or noisy: every phase/shape/inversion class requires reproducibility within both
# species and adequate amplitude in both. Flatness for the species-specific classes must be visible in
# **each** specimen/slice, not only in the species mean — averaging two slices can flatten a curve that
# neither slice has flat. And a large level effect never blocks a spatial label: level is not used in
# the classification at all.
#

# %%
# Purpose: the empirical distributions the section-4 thresholds act on, and how many genes each selects.

DIAGNOSTIC_QUANTILES = [0.5, 0.75, 0.9, 0.95, 0.99]

diagnostic_columns = [
    "mouse_amplitude", "human_amplitude", "mouse_reproducibility", "human_reproducibility",
    "shape_corr", "registered_shape_corr", "amplitude_log2_ratio_human_over_mouse",
    "best_shift_human_minus_mouse", "level_effect_human_minus_mouse",
]

threshold_distributions = pd.DataFrame({
    column: gene_metrics[column].abs().quantile(DIAGNOSTIC_QUANTILES)
    for column in diagnostic_columns
}).T

threshold_distributions.columns = [f"q{int(q * 100)}" for q in DIAGNOSTIC_QUANTILES]
threshold_distributions["median"] = gene_metrics[diagnostic_columns].abs().median()
threshold_distributions["n_finite"] = gene_metrics[diagnostic_columns].notna().sum()

print("Absolute-value quantiles of the metric distributions (thresholds are read against these):")
display(threshold_distributions.round(4))

enough_support = gene_metrics["enough_shared_support"]
discovery_eligible = (
    enough_support
    & gene_metrics["detected_in_both_species"]
    & ~gene_metrics["technical_gene"]
    & ~gene_metrics["axis_basis_gene"]
)

threshold_selection = pd.DataFrame([
    {"criterion": "patterned in mouse (amplitude >= amplitude_patterned)",
     "n_genes": int((gene_metrics["mouse_amplitude"] >= CONFIG["amplitude_patterned"]).sum())},
    {"criterion": "patterned in human (amplitude >= amplitude_patterned)",
     "n_genes": int((gene_metrics["human_amplitude"] >= CONFIG["amplitude_patterned"]).sum())},
    {"criterion": "flat in mouse (amplitude <= amplitude_flat)",
     "n_genes": int((gene_metrics["mouse_amplitude"] <= CONFIG["amplitude_flat"]).sum())},
    {"criterion": "flat in human (amplitude <= amplitude_flat)",
     "n_genes": int((gene_metrics["human_amplitude"] <= CONFIG["amplitude_flat"]).sum())},
    {"criterion": "reproducible within mouse (r >= within_species_corr)",
     "n_genes": int((gene_metrics["mouse_reproducibility"] >= CONFIG["within_species_corr"]).sum())},
    {"criterion": "reproducible within human (r >= within_species_corr)",
     "n_genes": int((gene_metrics["human_reproducibility"] >= CONFIG["within_species_corr"]).sum())},
    {"criterion": "shape agreement (shape_corr >= conserved_corr)",
     "n_genes": int((gene_metrics["shape_corr"] >= CONFIG["conserved_corr"]).sum())},
    {"criterion": "major amplitude change (|log2 ratio| >= amplitude_change_major)",
     "n_genes": int((gene_metrics["amplitude_log2_ratio_human_over_mouse"].abs()
                     >= CONFIG["amplitude_change_major"]).sum())},
    {"criterion": "displacement found (|best shift| >= minimum_shift)",
     "n_genes": int((gene_metrics["best_shift_human_minus_mouse"].abs()
                     >= CONFIG["minimum_shift"]).sum())},
    {"criterion": "displacement improves agreement (improvement >= shift_improvement)",
     "n_genes": int((gene_metrics["shift_improvement"] >= CONFIG["shift_improvement"]).sum())},
    {"criterion": "discovery-eligible (support + both species + no axis/technical flag)",
     "n_genes": int(discovery_eligible.sum())},
])

display(threshold_selection)
_report("discovery-eligible genes", int(discovery_eligible.sum()), len(gene_metrics),
        why="insufficient shared support, not detected in both species, technical or axis-basis")


# %%
# Purpose: the spatial phenotype classes, the flags behind them, and the per-class summary.

PHENOTYPE_CLASSES = [
    "conserved zonation",
    "weaker zonation in human",
    "stronger zonation in human",
    "mouse-zonated / human-flat",
    "human-zonated / mouse-flat",
    "shifted earlier in human",
    "shifted later in human",
    "gradient inversion",
    "complex shape rewiring",
    "weak / uncertain zonation",
]

# Fine class -> broad family, used for the sensitivity scoring in section 5.
BROAD_PHENOTYPE = {
    "conserved zonation": "conserved",
    "weaker zonation in human": "amplitude-change",
    "stronger zonation in human": "amplitude-change",
    "mouse-zonated / human-flat": "species-specific zonation",
    "human-zonated / mouse-flat": "species-specific zonation",
    "shifted earlier in human": "phase",
    "shifted later in human": "phase",
    "gradient inversion": "inversion",
    "complex shape rewiring": "complex rewiring",
    "weak / uncertain zonation": "no spatial signal",
    "excluded": "excluded",
}


def _phenotype_flags(metrics, config, include_axis_genes=False):
    """Every boolean the classification uses, so the classes can be rebuilt differently."""
    mouse_amplitude = metrics["mouse_amplitude"]
    human_amplitude = metrics["human_amplitude"]
    amplitude_log2 = metrics["amplitude_log2_ratio_human_over_mouse"]
    shift = metrics["best_shift_human_minus_mouse"]

    mouse_slice_columns = [f"amplitude__{name}" for name in mouse_samples]
    human_slice_columns = [f"amplitude__{name}" for name in human_samples]

    # "Flat" must be visible in every specimen/slice, not only in the averaged species curve.
    no_slice_zonated_mouse = metrics[mouse_slice_columns].max(axis=1) < config["amplitude_patterned"]
    no_slice_zonated_human = metrics[human_slice_columns].max(axis=1) < config["amplitude_patterned"]

    flags = {}
    flags["excluded"] = (
        ~metrics["detected_in_both_species"]
        | metrics["technical_gene"]
        | (metrics["axis_basis_gene"] & ~include_axis_genes)
    )
    flags["enough_support"] = metrics["enough_shared_support"]

    flags["strong_mouse"] = mouse_amplitude >= config["amplitude_patterned"]
    flags["strong_human"] = human_amplitude >= config["amplitude_patterned"]

    # NaN reproducibility is "not assessable" (a species with a single specimen) and is permissive.
    flags["reproducible_mouse"] = (
        metrics["mouse_reproducibility"].isna()
        | metrics["mouse_reproducibility"].ge(config["within_species_corr"])
    )
    flags["reproducible_human"] = (
        metrics["human_reproducibility"].isna()
        | metrics["human_reproducibility"].ge(config["within_species_corr"])
    )

    flags["flat_mouse_every_slice"] = (
        (mouse_amplitude <= config["amplitude_flat"]) & no_slice_zonated_mouse
    )
    flags["flat_human_every_slice"] = (
        (human_amplitude <= config["amplitude_flat"]) & no_slice_zonated_human
    )

    flags["patterned_both"] = flags["strong_mouse"] & flags["strong_human"] & flags["enough_support"]
    flags["positional_usable"] = (
        flags["patterned_both"] & flags["reproducible_mouse"] & flags["reproducible_human"]
    )

    flags["inversion"] = (
        flags["positional_usable"]
        & (np.sign(metrics["mouse_early_to_late"]) != np.sign(metrics["human_early_to_late"]))
        & metrics["mouse_early_to_late"].abs().ge(config["gradient_min_abs"])
        & metrics["human_early_to_late"].abs().ge(config["gradient_min_abs"])
        & metrics["shape_corr"].le(config["inversion_max_corr"])
    )

    # A displacement counts only when it survives both coordinate choices with the same sign and
    # improves agreement on both.
    flags["shift_usable"] = (
        flags["positional_usable"]
        & metrics["registered_shape_corr"].ge(config["registered_corr"])
        & shift.abs().ge(config["minimum_shift"])
        & metrics["shift_improvement"].ge(config["shift_improvement"])
        & metrics["registered_shape_corr_unregistered"].ge(config["registered_corr"])
        & metrics["shift_improvement_unregistered"].ge(config["shift_improvement"])
        & (np.sign(shift) == np.sign(metrics["best_shift_human_minus_mouse_unregistered"]))
        & metrics["best_shift_human_minus_mouse_unregistered"].abs().ge(config["minimum_shift"])
    )
    flags["shifted_earlier"] = flags["shift_usable"] & shift.lt(0)
    flags["shifted_later"] = flags["shift_usable"] & shift.gt(0)

    flags["conserved"] = (
        flags["positional_usable"]
        & metrics["shape_corr"].ge(config["conserved_corr"])
        & shift.abs().lt(config["minimum_shift"])
        & amplitude_log2.abs().le(config["amplitude_change_small"])
    )
    # The amplitude classes deliberately do NOT require both curves to clear the "patterned"
    # threshold: a human curve that has dropped below it is the clearest case of weakening. They
    # require reproducibility in both species, a reproducible stronger-shaped curve, and shape
    # agreement, so a flat or noisy curve cannot be given an amplitude label (a flat curve has no
    # usable shape correlation at all, and the species-specific classes take priority anyway).
    flags["amplitude_class_usable"] = (
        flags["enough_support"]
        & flags["reproducible_mouse"]
        & flags["reproducible_human"]
        & (np.maximum(mouse_amplitude, human_amplitude) >= config["amplitude_patterned"])
    )
    flags["weaker_in_human"] = (
        flags["amplitude_class_usable"]
        & metrics["shape_corr"].ge(config["conserved_corr"])
        & amplitude_log2.le(-config["amplitude_change_major"])
    )
    flags["stronger_in_human"] = (
        flags["amplitude_class_usable"]
        & metrics["shape_corr"].ge(config["conserved_corr"])
        & amplitude_log2.ge(config["amplitude_change_major"])
    )

    flags["mouse_zonated_human_flat"] = (
        flags["strong_mouse"] & flags["reproducible_mouse"] & flags["flat_human_every_slice"]
    )
    flags["human_zonated_mouse_flat"] = (
        flags["strong_human"] & flags["reproducible_human"] & flags["flat_mouse_every_slice"]
    )

    flags["complex_rewiring"] = (
        flags["positional_usable"]
        & metrics["shape_corr"].lt(config["conserved_corr"])
        & ~flags["shift_usable"]
        & ~flags["inversion"]
        & amplitude_log2.abs().lt(config["amplitude_change_major"])
    )

    return pd.DataFrame(flags, index=metrics.index)


def _assign_phenotypes(flags):
    """Phenotype per gene. Later entries overwrite earlier ones: the list runs in priority order."""
    labels = np.full(len(flags), "weak / uncertain zonation", dtype=object)

    assignments = [
        (flags["complex_rewiring"], "complex shape rewiring"),
        (flags["stronger_in_human"], "stronger zonation in human"),
        (flags["weaker_in_human"], "weaker zonation in human"),
        (flags["conserved"], "conserved zonation"),
        (flags["shifted_earlier"], "shifted earlier in human"),
        (flags["shifted_later"], "shifted later in human"),
        (flags["inversion"], "gradient inversion"),
        (flags["human_zonated_mouse_flat"], "human-zonated / mouse-flat"),
        (flags["mouse_zonated_human_flat"], "mouse-zonated / human-flat"),
        (~flags["enough_support"], "weak / uncertain zonation"),
        (flags["excluded"], "excluded"),
    ]

    for mask, label in assignments:
        labels[np.asarray(mask, dtype=bool)] = label
    return labels


phenotype_flags = _phenotype_flags(gene_metrics, CONFIG)
gene_metrics = gene_metrics.join(phenotype_flags)
gene_metrics["spatial_phenotype"] = _assign_phenotypes(phenotype_flags)
gene_metrics["broad_phenotype"] = gene_metrics["spatial_phenotype"].map(BROAD_PHENOTYPE)

# A divergence score, not a probability: each term is in units of its own threshold, so 1 unit of the
# score is "one threshold's worth" of amplitude change, displacement or shape divergence.
gene_metrics["spatial_discovery_score"] = (
    (gene_metrics["amplitude_log2_ratio_human_over_mouse"].abs() / CONFIG["amplitude_change_major"]).fillna(0.0)
    + (gene_metrics["best_shift_human_minus_mouse"].abs() / CONFIG["max_registration_shift"]).fillna(0.0)
    + (1.0 - gene_metrics["shape_corr"].clip(-1.0, 1.0).fillna(1.0))
)

phenotype_summary = pd.DataFrame([
    {
        "spatial_phenotype": label,
        "n_genes": int((gene_metrics["spatial_phenotype"] == label).sum()),
        "median_abundance": float(np.nanmedian(np.mean(
            [gene_metrics.loc[gene_metrics["spatial_phenotype"] == label, "abundance_mouse"],
             gene_metrics.loc[gene_metrics["spatial_phenotype"] == label, "abundance_human"]], axis=0))),
        "median_detection": float(np.nanmedian(np.mean(
            [gene_metrics.loc[gene_metrics["spatial_phenotype"] == label, "detection_mouse"],
             gene_metrics.loc[gene_metrics["spatial_phenotype"] == label, "detection_human"]], axis=0))),
        "median_within_species_reproducibility": float(np.nanmedian(np.mean(
            [gene_metrics.loc[gene_metrics["spatial_phenotype"] == label, "mouse_reproducibility"],
             gene_metrics.loc[gene_metrics["spatial_phenotype"] == label, "human_reproducibility"]], axis=0))),
        "median_level_effect": float(np.nanmedian(
            gene_metrics.loc[gene_metrics["spatial_phenotype"] == label, "level_effect_human_minus_mouse"])),
        "median_mouse_amplitude": float(np.nanmedian(
            gene_metrics.loc[gene_metrics["spatial_phenotype"] == label, "mouse_amplitude"])),
        "median_human_amplitude": float(np.nanmedian(
            gene_metrics.loc[gene_metrics["spatial_phenotype"] == label, "human_amplitude"])),
        "median_shape_corr": float(np.nanmedian(
            gene_metrics.loc[gene_metrics["spatial_phenotype"] == label, "shape_corr"])),
        "median_registered_shape_corr": float(np.nanmedian(
            gene_metrics.loc[gene_metrics["spatial_phenotype"] == label, "registered_shape_corr"])),
        "median_abs_shift": float(np.nanmedian(np.abs(
            gene_metrics.loc[gene_metrics["spatial_phenotype"] == label, "best_shift_human_minus_mouse"]))),
        "median_discovery_score": float(np.nanmedian(
            gene_metrics.loc[gene_metrics["spatial_phenotype"] == label, "spatial_discovery_score"])),
    }
    for label in PHENOTYPE_CLASSES + ["excluded"]
])

display(phenotype_summary.round(4))
_save_table(phenotype_summary, "gene_spatial_phenotype_summary.csv")

print("Spatial phenotype counts (all tested orthologs):")
print(gene_metrics["spatial_phenotype"].value_counts().to_string())


# %%
# Purpose: synthetic self-check of the phenotype logic - the classes must come out as intended.

def _phenotype_selfcheck():
    """Build synthetic metric rows with known expected classes and check the classifier.

    The classification is the one piece of this notebook that is easy to get subtly wrong (a
    threshold in the wrong place, a priority order inverted, an amplitude class that silently
    requires the very thing it is meant to detect), and it is not covered by the curve-level
    self-check at the end. Every row here is a constructed metric profile, not data.
    """
    reference = dict(
        mouse_amplitude=0.4, human_amplitude=0.4, amplitude_log2_ratio_human_over_mouse=0.0,
        best_shift_human_minus_mouse=0.0, best_shift_human_minus_mouse_unregistered=0.0,
        registered_shape_corr=0.95, registered_shape_corr_unregistered=0.95,
        shift_improvement=0.01, shift_improvement_unregistered=0.01,
        shape_corr=0.9, mouse_reproducibility=0.8, human_reproducibility=0.8,
        mouse_early_to_late=0.3, human_early_to_late=0.3,
        detected_in_both_species=True, technical_gene=False, axis_basis_gene=False,
        enough_shared_support=True,
        **{f"amplitude__{sample}": 0.4 for sample in mouse_samples + human_samples},
    )

    def gene(**overrides):
        row = dict(reference)
        row.update(overrides)
        return row

    synthetic = {
        "expected conserved": (gene(), "conserved zonation"),
        "expected weaker in human": (
            gene(human_amplitude=0.1, amplitude_log2_ratio_human_over_mouse=-2.0),
            "weaker zonation in human"),
        "expected stronger in human": (
            gene(human_amplitude=0.8, amplitude_log2_ratio_human_over_mouse=1.0),
            "stronger zonation in human"),
        "expected mouse-zonated / human-flat": (
            gene(human_amplitude=0.03, amplitude_log2_ratio_human_over_mouse=-3.7, shape_corr=0.1,
                 registered_shape_corr=0.2,
                 **{f"amplitude__{sample}": 0.03 for sample in human_samples}),
            "mouse-zonated / human-flat"),
        "expected human-zonated / mouse-flat": (
            gene(mouse_amplitude=0.03, amplitude_log2_ratio_human_over_mouse=3.7, shape_corr=0.1,
                 registered_shape_corr=0.2,
                 **{f"amplitude__{sample}": 0.03 for sample in mouse_samples}),
            "human-zonated / mouse-flat"),
        "expected shifted later": (
            gene(shape_corr=0.4, best_shift_human_minus_mouse=0.15,
                 best_shift_human_minus_mouse_unregistered=0.14, registered_shape_corr=0.93,
                 registered_shape_corr_unregistered=0.92, shift_improvement=0.4,
                 shift_improvement_unregistered=0.35),
            "shifted later in human"),
        "expected shifted earlier": (
            gene(shape_corr=0.4, best_shift_human_minus_mouse=-0.15,
                 best_shift_human_minus_mouse_unregistered=-0.14, registered_shape_corr=0.93,
                 registered_shape_corr_unregistered=0.92, shift_improvement=0.4,
                 shift_improvement_unregistered=0.35),
            "shifted earlier in human"),
        "shift on one coordinate only stays complex": (
            gene(shape_corr=0.4, best_shift_human_minus_mouse=0.15,
                 best_shift_human_minus_mouse_unregistered=0.0, registered_shape_corr=0.93,
                 registered_shape_corr_unregistered=0.95, shift_improvement=0.4,
                 shift_improvement_unregistered=0.01),
            "complex shape rewiring"),
        "expected gradient inversion": (
            gene(shape_corr=-0.8, mouse_early_to_late=0.4, human_early_to_late=-0.4),
            "gradient inversion"),
        "expected complex rewiring": (
            gene(shape_corr=0.1, registered_shape_corr=0.3), "complex shape rewiring"),
        "flat in both species": (
            gene(mouse_amplitude=0.02, human_amplitude=0.02, shape_corr=np.nan,
                 **{f"amplitude__{sample}": 0.02 for sample in mouse_samples + human_samples}),
            "weak / uncertain zonation"),
        "flat in one human slice only is not human-flat": (
            gene(human_amplitude=0.03, amplitude_log2_ratio_human_over_mouse=-3.7, shape_corr=0.1,
                 registered_shape_corr=0.2, human_reproducibility=np.nan,
                 **{f"amplitude__{human_samples[0]}": 0.5,
                    f"amplitude__{human_samples[-1]}": 0.02}),
            "weak / uncertain zonation"),
        "axis-basis gene": (gene(axis_basis_gene=True), "excluded"),
        "technical gene": (gene(technical_gene=True), "excluded"),
        "not measured in both inputs": (gene(detected_in_both_species=False), "excluded"),
        "too little shared support": (gene(enough_shared_support=False), "weak / uncertain zonation"),
    }

    frame = pd.DataFrame({name: row for name, (row, _) in synthetic.items()}).T
    frame = frame.reset_index().rename(columns={"index": "gene"})
    for column in frame.columns:
        if column == "gene":
            continue
        frame[column] = frame[column].astype(float)

    for flag in ("detected_in_both_species", "technical_gene", "axis_basis_gene",
                 "enough_shared_support"):
        frame[flag] = frame[flag].astype(bool)

    flags = _phenotype_flags(frame, CONFIG)
    assigned = pd.Series(_assign_phenotypes(flags), index=frame["gene"])

    outcome = pd.DataFrame({
        "synthetic_case": list(synthetic),
        "expected_class": [expected for _, expected in synthetic.values()],
        "assigned_class": [str(assigned[name]) for name in synthetic],
    })
    outcome["passed"] = outcome["expected_class"] == outcome["assigned_class"]
    return outcome


phenotype_selfcheck = _phenotype_selfcheck()
display(phenotype_selfcheck)

if not phenotype_selfcheck["passed"].all():
    failed = phenotype_selfcheck.loc[~phenotype_selfcheck["passed"], "synthetic_case"].tolist()
    raise AssertionError(f"synthetic phenotype self-check failed for: {failed}")

print(f"All {len(phenotype_selfcheck)} synthetic phenotype self-checks passed.")


# %% [markdown]
# ## 5 - Sensitivity of every phenotype to the analysis choices
#
# **Why.** A phenotype class is a conjunction of thresholds applied to curves that were themselves
# built under choices (specimen balancing, coordinate registration). Reporting one partition and
# calling it the result would hide that. So the whole classification is re-derived under each choice
# and the per-gene agreement is kept.
#
# **The variants.** (1) pooled tubule-weighted curves instead of equal-weight specimen curves;
# (2) the unregistered coordinate, i.e. no landmark mapping - for this variant the cross-axis phase
# requirement collapses onto that single axis by construction, so a gene that keeps its shift class here
# kept it without any registration help; (3) leave out each mouse specimen; (4) leave out each human
# slice; (5) amplitude thresholds scaled by 0.8 and 1.25; (6) the displacement threshold scaled by 0.75
# and 1.25; (7) detection threshold at 2% and 10%; (8) axis-basis genes allowed into discovery.
#
# **How the variants are run.** Every variant re-derives landmarks and the mapping from *its own*
# curves - a leave-one-out variant does not inherit the primary registration. All variants share one
# code path (`_assemble_metrics` + `_phenotype_flags` + `_assign_phenotypes`), so a variant can only
# differ through the curves, the coordinate choice or the thresholds.
#
# **What the numbers mean, and what they do not.** `robustness_score` is the fraction of variants that
# keep a gene's **broad** phenotype (conserved / amplitude-change / species-specific zonation / phase /
# inversion / complex rewiring / no spatial signal / excluded). It is a stability score, not a
# significance measure. A gene that leaves the universe in a stricter-detection variant is counted as
# not retaining its phenotype; that is deliberate, because "no longer measured" is a real change in
# what can be claimed.
#
# **The one-donor caveat, stated plainly.** Leaving out a human slice is *slice sensitivity*. With a
# single donor it cannot be biological replication, and no variant here produces a p-value.
#

# %%
# Purpose: re-derive every phenotype under each analysis choice, and score per-gene stability.


def _sensitivity_variants():
    """Every variant the phenotypes are re-derived under: curves, coordinate choice and config."""
    variants = {
        "primary (registered, specimen-balanced)": {
            "mouse": balanced_mouse, "human": balanced_human, "register": True, "config": CONFIG,
            "note": "the primary analysis",
        },
        "pooled tubule-weighted curves": {
            "mouse": pooled_mouse, "human": pooled_human, "register": True, "config": CONFIG,
            "note": "pooled fit instead of equal-weight specimen curves",
        },
        "unregistered coordinate": {
            "mouse": balanced_mouse, "human": balanced_human, "register": False, "config": CONFIG,
            "note": "no landmark registration",
        },
    }

    for name in mouse_samples:
        remaining = [specimen for specimen in mouse_samples if specimen != name]
        variants[f"leave out mouse specimen {name}"] = {
            "mouse": specimen_balanced_curves(
                {specimen: specimen_curves[specimen] for specimen in remaining}
            ),
            "human": balanced_human, "register": True, "config": CONFIG,
            "note": "specimen sensitivity",
        }

    for name in human_samples:
        remaining = [specimen for specimen in human_samples if specimen != name]
        variants[f"leave out human slice {name}"] = {
            "mouse": balanced_mouse,
            "human": specimen_balanced_curves(
                {specimen: specimen_curves[specimen] for specimen in remaining}
            ),
            "register": True, "config": CONFIG,
            "note": "slice sensitivity, not replication (one donor)",
        }

    threshold_variants = {
        "amplitude thresholds -20%": {
            "amplitude_patterned": CONFIG["amplitude_patterned"] * 0.8,
            "amplitude_flat": CONFIG["amplitude_flat"] * 0.8,
        },
        "amplitude thresholds +25%": {
            "amplitude_patterned": CONFIG["amplitude_patterned"] * 1.25,
            "amplitude_flat": CONFIG["amplitude_flat"] * 1.25,
        },
        "shift threshold -25%": {"minimum_shift": CONFIG["minimum_shift"] * 0.75},
        "shift threshold +25%": {"minimum_shift": CONFIG["minimum_shift"] * 1.25},
        "detection threshold 2%": {"min_detection_each_species": 0.02},
        "detection threshold 10%": {"min_detection_each_species": 0.10},
        "axis-basis genes included": {},
    }

    for label, updates in threshold_variants.items():
        variant_config = dict(CONFIG)
        variant_config.update(updates)
        is_axis_variant = label == "axis-basis genes included"
        variants[label] = {
            "mouse": balanced_mouse, "human": balanced_human, "register": True,
            "config": variant_config, "include_axis": is_axis_variant,
            "note": ("axis-basis genes treated as ordinary discoverable genes - the reverse of the "
                     "primary setting, so it shows what the exclusion is hiding"
                     if is_axis_variant else "threshold / universe sensitivity"),
        }

    return variants


sensitivity_variants = _sensitivity_variants()
sensitivity_summary_rows = []
variant_phenotypes = {}

for variant_name, variant_spec in sensitivity_variants.items():
    variant_metrics, variant_registration = _assemble_metrics(
        variant_spec["mouse"], variant_spec["human"], register=variant_spec["register"]
    )

    # The detection and support thresholds act on the universe, not on the curves.
    variant_metrics["detected_in_both_species"] = (
        (variant_metrics["detection_mouse"] >= variant_spec["config"]["min_detection_each_species"])
        & (variant_metrics["detection_human"] >= variant_spec["config"]["min_detection_each_species"])
    )
    variant_metrics["enough_shared_support"] = (
        variant_metrics["n_shared_grid_points"] >= variant_spec["config"]["min_pattern_grid_points"]
    )

    variant_flags = _phenotype_flags(
        variant_metrics, variant_spec["config"],
        include_axis_genes=variant_spec.get("include_axis", False),
    )
    variant_labels = pd.Series(_assign_phenotypes(variant_flags), index=gene_names)
    variant_phenotypes[variant_name] = variant_labels

    primary_labels = gene_metrics["spatial_phenotype"].to_numpy()
    primary_broad = gene_metrics["broad_phenotype"].to_numpy()
    variant_broad = variant_labels.map(BROAD_PHENOTYPE).to_numpy()

    primary_score = gene_metrics["spatial_discovery_score"].to_numpy()
    variant_score = variant_metrics["spatial_discovery_score"].to_numpy()
    order_primary = np.argsort(-np.nan_to_num(primary_score))
    order_variant = np.argsort(-np.nan_to_num(variant_score))

    sensitivity_summary_rows.append({
        "variant": variant_name,
        "note": variant_spec["note"],
        "n_landmarks_used": (np.nan if variant_registration is None
                             else variant_registration["n_landmarks"]),
        "same_fine_phenotype_share": float((variant_labels.to_numpy() == primary_labels).mean()),
        "same_broad_phenotype_share": float((variant_broad == primary_broad).mean()),
        "n_discovery_eligible": int((~variant_flags["excluded"] & variant_flags["enough_support"]).sum()),
        "n_conserved": int((variant_labels == "conserved zonation").sum()),
        "n_amplitude_change": int(variant_labels.isin(
            ["weaker zonation in human", "stronger zonation in human"]).sum()),
        "n_species_specific": int(variant_labels.isin(
            ["mouse-zonated / human-flat", "human-zonated / mouse-flat"]).sum()),
        "n_phase": int(variant_labels.isin(
            ["shifted earlier in human", "shifted later in human"]).sum()),
        "n_inversion": int((variant_labels == "gradient inversion").sum()),
        "n_complex": int((variant_labels == "complex shape rewiring").sum()),
        "spearman_discovery_score_vs_primary": float(
            spearmanr(primary_score, variant_score, nan_policy="omit").correlation),
        "top50_discovery_overlap": len(
            set(gene_names[order_primary[:50]]) & set(gene_names[order_variant[:50]])
        ) / 50.0,
    })

sensitivity_summary = pd.DataFrame(sensitivity_summary_rows)
display(sensitivity_summary.round(4))
_save_table(sensitivity_summary, "sensitivity_summary.csv")


# %%
# Purpose: per-gene stability of the phenotype across the sensitivity variants.

# --------------------------------------------------------------------------------------
# Per-gene stability across the variants.
# --------------------------------------------------------------------------------------

robustness = pd.DataFrame({
    "gene": gene_metrics["gene"],
    "spatial_phenotype": gene_metrics["spatial_phenotype"],
    "broad_phenotype": gene_metrics["broad_phenotype"],
})

for variant_name, variant_labels in variant_phenotypes.items():
    robustness[f"phenotype__{variant_name}"] = variant_labels.to_numpy()
    robustness[f"broad__{variant_name}"] = variant_labels.map(BROAD_PHENOTYPE).to_numpy()

broad_columns = [column for column in robustness if column.startswith("broad__")]
robustness["n_variants_tested"] = robustness[broad_columns].notna().sum(axis=1)
robustness["n_variants_same_broad_phenotype"] = (
    robustness[broad_columns].to_numpy() == robustness["broad_phenotype"].to_numpy()[:, None]
).sum(axis=1)
robustness["robustness_score"] = (
    robustness["n_variants_same_broad_phenotype"]
    / robustness["n_variants_tested"].replace(0, np.nan)
)

gene_metrics = gene_metrics.merge(
    robustness[["gene", "n_variants_tested", "n_variants_same_broad_phenotype", "robustness_score"]],
    on="gene", how="left",
)

_save_table(robustness, "gene_spatial_phenotype_robustness.csv")

print("Median robustness_score per primary phenotype (1.0 = unchanged broad class in every variant):")
display(
    robustness.groupby("spatial_phenotype")["robustness_score"]
    .agg(genes="count", median_robustness="median", min_robustness="min")
    .sort_values("median_robustness", ascending=False)
    .round(3)
)

discovery_genes = ~gene_metrics["axis_basis_gene"] & ~gene_metrics["technical_gene"]
_report("genes unchanged in their broad class under every variant",
        int((gene_metrics["robustness_score"] == 1).sum()), len(gene_names))
_report("discovery-eligible genes robust under >= 80% of variants",
        int((gene_metrics["robustness_score"].ge(0.8) & discovery_genes).sum()), len(gene_names),
        why="used to prioritise the genes shown in sections 6 and 9")


# %% [markdown]
# ## 6 - Primary figures
#
# Five figures, each answering one part of the question:
#
# * **1A** — are homologous PT territories being compared at all? Canonical S1/S2/S3 programs and their
#   member genes on the registered coordinate.
# * **1B** — where does zonation strength sit relative to the identity line, and which genes have strong
#   zonation in both species but poor shape agreement (the rewiring corner)?
# * **1C** — phase versus residual mismatch: is a gene explained by displacement, or is it genuinely
#   reshaped after the best possible displacement?
# * **1D** — representative genes per class, with **both** specimens/slices drawn faintly behind the
#   balanced curves, so within-species consistency is visible rather than asserted.
# * **1E** — shape only: each gene z-scored *within each species*, so the heatmap cannot be driven by
#   an abundance difference between human and mouse.
#
# Genes are chosen by `spatial_discovery_score` within their class, preferring genes that kept their
# broad class in at least 80% of the section-5 variants. Figures are readable without the code: every
# panel states its class, sample sizes and direction convention.
#

# %%
# Purpose: figure 1A - PT alignment on the registered coordinate.

PHENOTYPE_COLORS = {
    "conserved zonation": "#1B7837",
    "weaker zonation in human": "#8C6BB1",
    "stronger zonation in human": "#762A83",
    "mouse-zonated / human-flat": "#0072B2",
    "human-zonated / mouse-flat": "#D55E00",
    "shifted earlier in human": "#E69F00",
    "shifted later in human": "#9E4F00",
    "gradient inversion": "#C2185B",
    "complex shape rewiring": "#4D4D4D",
    "weak / uncertain zonation": "#BDBDBD",
    "excluded": "#E8E8E8",
}


def _top_genes_for_class(label, n, prefer_robust=True):
    """Top `n` genes of a class by divergence score, preferring genes stable across variants."""
    block = gene_metrics[
        gene_metrics["spatial_phenotype"].eq(label)
        & ~gene_metrics["axis_basis_gene"]
        & ~gene_metrics["technical_gene"]
        & gene_metrics["detected_in_both_species"]
    ]
    if prefer_robust and len(block):
        stable = block[block["robustness_score"].ge(0.8)]
        if len(stable) >= n:
            block = stable
    return block.sort_values("spatial_discovery_score", ascending=False).head(n)


discovery_eligible = (
    ~gene_metrics["axis_basis_gene"]
    & ~gene_metrics["technical_gene"]
    & gene_metrics["detected_in_both_species"]
    & gene_metrics["enough_shared_support"]
)

figure, axes = plt.subplots(1, 3, figsize=(13, 4.3), sharey=True)

for axis, (panel, requested) in zip(axes, LANDMARK_PANELS.items()):
    members = [gene for gene in requested if gene.upper() in gene_lookup_tested]
    for gene in members:
        index = gene_lookup_tested[gene.upper()]
        axis.plot(grid_unit, balanced_mouse[index], color=SPECIES_COLORS["mouse"], lw=1.0, alpha=0.55)
        axis.plot(grid_unit, human_registered[index], color=SPECIES_COLORS["human"], lw=1.0,
                  alpha=0.55, ls="--")
        axis.annotate(gene, (grid_unit[int(np.nanargmax(balanced_mouse[index]))],
                             float(np.nanmax(balanced_mouse[index]))),
                      xytext=(2, 2), textcoords="offset points", fontsize=7)

    if panel in primary_registration["programs"]:
        programs = primary_registration["programs"][panel]
        axis.plot(grid_unit, programs["mouse"], color=SPECIES_COLORS["mouse"], lw=2.4)
        axis.plot(grid_unit, programs["human"], color=SPECIES_COLORS["human"], lw=2.4, ls="--")

    axis.set_title(f"{panel} program: {', '.join(members)}", loc="left")
    axis.set_xlabel("Registered PT position (unit): early -> late")

axes[0].set_ylabel("Fitted lognorm (genes) / mean z-score (program)")
axes[0].legend(
    handles=[
        Line2D([], [], color=SPECIES_COLORS["mouse"], lw=2.4, label="mouse (balanced)"),
        Line2D([], [], color=SPECIES_COLORS["human"], lw=2.4, ls="--", label="human (registered)"),
    ], loc="upper right",
)
figure.suptitle("Figure 1A - Canonical S1/S2/S3 territory alignment on the registered coordinate\n"
                "landmark genes thin, program curves bold; two human slices, one donor",
                fontsize=11)
figure.tight_layout()
_save_figure(figure, "fig1a_pt_alignment.png")
plt.show()


# %%
# Purpose: figure 1B - the zonation conservation landscape.

landscape = gene_metrics[discovery_eligible].dropna(
    subset=["mouse_amplitude", "human_amplitude", "shape_corr"]
)

figure, axes = plt.subplots(1, 2, figsize=(13.5, 6.0), width_ratios=[1.15, 1.0])
ax_main, ax_margin = axes

points = ax_main.scatter(
    landscape["mouse_amplitude"], landscape["human_amplitude"],
    c=landscape["shape_corr"], cmap="coolwarm", vmin=-1, vmax=1,
    s=6, alpha=0.6, linewidths=0,
)
limit = float(np.nanpercentile(
    np.concatenate([landscape["mouse_amplitude"], landscape["human_amplitude"]]), 99.5))
limit = max(limit, 0.5)

ax_main.plot([0, limit], [0, limit], color="0.4", lw=1, ls=":", label="equal zonation strength")
ax_main.axvline(CONFIG["amplitude_patterned"], color="0.6", lw=0.8, ls="--")
ax_main.axhline(CONFIG["amplitude_patterned"], color="0.6", lw=0.8, ls="--")
ax_main.set_xlim(-0.02, limit)
ax_main.set_ylim(-0.02, limit)
ax_main.set_xlabel("Mouse zonation amplitude (peak-to-trough lognorm)")
ax_main.set_ylabel("Human zonation amplitude (peak-to-trough lognorm)")

label_candidates = pd.concat([
    _top_genes_for_class("conserved zonation", 4),
    _top_genes_for_class("weaker zonation in human", 3),
    _top_genes_for_class("stronger zonation in human", 3),
    _top_genes_for_class("complex shape rewiring", 3),
    _top_genes_for_class("gradient inversion", 2),
])
for row in label_candidates.itertuples():
    ax_main.annotate(row.gene, (row.mouse_amplitude, row.human_amplitude),
                     xytext=(4, 3), textcoords="offset points", fontsize=7.5)

ax_main.text(0.02, 0.97, "above the line:\nstronger zonation in human", transform=ax_main.transAxes,
             fontsize=8, va="top")
ax_main.text(0.02, 0.03, "below the line:\nweaker zonation in human", transform=ax_main.transAxes,
             fontsize=8, va="bottom")
ax_main.legend(loc="lower right")
ax_main.set_title(f"Figure 1B - Zonation strength, {len(landscape):,} discovery-eligible orthologs\n"
                  "colour = standardised curve correlation (level- and amplitude-free)",
                  loc="left", fontsize=10)

colour_bar = figure.colorbar(points, ax=ax_main, shrink=0.85)
colour_bar.set_label("Standardised mouse-human curve correlation")

# Right panel: the rewiring corner, stated numerically.
strong_both = landscape[
    landscape["mouse_amplitude"].ge(CONFIG["amplitude_patterned"])
    & landscape["human_amplitude"].ge(CONFIG["amplitude_patterned"])
]
quadrants = pd.DataFrame([
    {"group": "strong in both, shape conserved (r >= conserved_corr)",
     "n": int((strong_both["shape_corr"] >= CONFIG["conserved_corr"]).sum())},
    {"group": "strong in both, shape intermediate",
     "n": int(((strong_both["shape_corr"] < CONFIG["conserved_corr"])
               & (strong_both["shape_corr"] > CONFIG["rewired_max_corr"])).sum())},
    {"group": "strong in both, poor shape (r <= rewired_max_corr)",
     "n": int((strong_both["shape_corr"] <= CONFIG["rewired_max_corr"]).sum())},
    {"group": "mouse-dominant (mouse patterned, human flat)",
     "n": int(((landscape["mouse_amplitude"] >= CONFIG["amplitude_patterned"])
               & (landscape["human_amplitude"] <= CONFIG["amplitude_flat"])).sum())},
    {"group": "human-dominant (human patterned, mouse flat)",
     "n": int(((landscape["human_amplitude"] >= CONFIG["amplitude_patterned"])
               & (landscape["mouse_amplitude"] <= CONFIG["amplitude_flat"])).sum())},
])

ax_margin.barh(np.arange(len(quadrants)), quadrants["n"], color="#555555", alpha=0.85)
ax_margin.set_yticks(np.arange(len(quadrants)))
ax_margin.set_yticklabels([textwrap.fill(group, 34) for group in quadrants["group"]], fontsize=8)
ax_margin.set_xlabel("Genes")
ax_margin.set_title("Where the signal sits", loc="left", fontsize=10)
for index, value in enumerate(quadrants["n"]):
    ax_margin.text(value, index, f"  {value:,}", va="center", fontsize=8)

figure.tight_layout()
_save_figure(figure, "zonation_amplitude_human_vs_mouse.png")
plt.show()

_save_table(landscape, "zonation_landscape_eligible_genes.csv")


# %%
# Purpose: figure 1C - positional displacement versus residual shape mismatch.

phase_frame = gene_metrics[discovery_eligible].dropna(
    subset=["shape_corr", "registered_shape_corr", "best_shift_human_minus_mouse",
            "residual_rms_after_shift"]
)

figure, axes = plt.subplots(1, 3, figsize=(16, 5.2), width_ratios=[1.0, 1.0, 1.15])
ax_raw_registered, ax_shift_residual, ax_summary = axes

# Panel 1: does a bounded displacement rescue the agreement?
ax_raw_registered.plot([-1, 1], [-1, 1], color="0.5", lw=1, ls=":", label="no change from shifting")
ax_raw_registered.scatter(
    phase_frame["shape_corr"], phase_frame["registered_shape_corr"],
    s=6, alpha=0.45, linewidths=0, color="#666666",
)
rescued = phase_frame[phase_frame["shift_improvement"] >= CONFIG["shift_improvement"]]
ax_raw_registered.scatter(
    rescued["shape_corr"], rescued["registered_shape_corr"],
    s=9, alpha=0.8, linewidths=0, color=PHENOTYPE_COLORS["shifted later in human"],
)
ax_raw_registered.set_xlabel("Standardised curve correlation (registered axis)")
ax_raw_registered.set_ylabel("Best correlation after displacement\n(|shift| <= 20% of PT)")
ax_raw_registered.set_title("Phase: a bounded displacement explains the difference\n"
                            f"orange = improvement >= {CONFIG['shift_improvement']} "
                            f"(n={len(rescued):,} of {len(phase_frame):,})", loc="left", fontsize=10)
ax_raw_registered.legend(loc="upper left")

# Panel 2: displacement against what is left after it.
for label in ["shifted earlier in human", "shifted later in human", "gradient inversion",
              "complex shape rewiring", "conserved zonation"]:
    block = phase_frame[phase_frame["spatial_phenotype"].eq(label)]
    if not len(block):
        continue
    ax_shift_residual.scatter(
        block["best_shift_human_minus_mouse"], block["residual_rms_after_shift"],
        s=9, alpha=0.7, linewidths=0, color=PHENOTYPE_COLORS[label], label=f"{label} (n={len(block):,})",
    )
ax_shift_residual.axvline(0, color="black", lw=0.8)
ax_shift_residual.set_xlabel("Best displacement of the human program (fraction of PT)\n"
                             "<- earlier in human     later in human ->")
ax_shift_residual.set_ylabel("Standardised RMS difference after the best displacement")
ax_shift_residual.set_title("Phase versus true rewiring\n"
                            "high residual at small displacement = reshaped, not moved",
                            loc="left", fontsize=10)
ax_shift_residual.legend(loc="upper center", fontsize=7.5)

class_counts = (
    gene_metrics[gene_metrics["spatial_phenotype"].isin(PHENOTYPE_CLASSES)]
    .groupby("spatial_phenotype").size().reindex(PHENOTYPE_CLASSES).fillna(0).astype(int)
)
robust_counts = (
    gene_metrics[gene_metrics["spatial_phenotype"].isin(PHENOTYPE_CLASSES)
                 & gene_metrics["robustness_score"].ge(0.8)]
    .groupby("spatial_phenotype").size().reindex(PHENOTYPE_CLASSES).fillna(0).astype(int)
)

y_positions = np.arange(len(PHENOTYPE_CLASSES))
ax_summary.barh(y_positions, class_counts.to_numpy(), color="#BBBBBB", label="all genes in class")
ax_summary.barh(y_positions, robust_counts.to_numpy(), color="#333333",
                label=f"broad class kept in >= 80% of variants")
ax_summary.set_yticks(y_positions)
ax_summary.set_yticklabels([textwrap.fill(label, 26) for label in PHENOTYPE_CLASSES], fontsize=8)
ax_summary.set_xlabel("Genes")
ax_summary.set_title("Class sizes, and how many survive the sensitivity panel", loc="left", fontsize=10)
ax_summary.legend(loc="lower right")
for index, (total, robust) in enumerate(zip(class_counts.to_numpy(), robust_counts.to_numpy())):
    ax_summary.text(total, index, f"  {total:,}", va="center", fontsize=7.5)

figure.suptitle("Figure 1C - Separating positional displacement from residual trajectory rewiring",
                fontsize=11)
figure.tight_layout()
_save_figure(figure, "phase_vs_shape_rewiring.png")
plt.show()


# %%
# Purpose: figure 1D - representative genes per class, with the individual specimens/slices shown.

REPRESENTATIVE_CLASSES = [
    "conserved zonation",
    "weaker zonation in human",
    "stronger zonation in human",
    "mouse-zonated / human-flat",
    "human-zonated / mouse-flat",
    "shifted earlier in human",
    "shifted later in human",
    "gradient inversion",
    "complex shape rewiring",
]
GENES_PER_CLASS = CONFIG["representative_genes_per_class"]

representative_table = pd.concat(
    [_top_genes_for_class(label, GENES_PER_CLASS) for label in REPRESENTATIVE_CLASSES],
    ignore_index=True,
)

representative_class_counts = (
    representative_table.groupby("spatial_phenotype").size()
    .reindex(REPRESENTATIVE_CLASSES).fillna(0).astype(int)
)
print("Representative genes available per class (a class with no eligible gene is left blank):")
display(representative_class_counts.to_frame("n_genes_shown"))

figure, axes = plt.subplots(
    len(REPRESENTATIVE_CLASSES), GENES_PER_CLASS,
    figsize=(3.4 * GENES_PER_CLASS, 1.85 * len(REPRESENTATIVE_CLASSES)),
    sharex=True,
)
axes = np.atleast_2d(axes)

for row_index, label in enumerate(REPRESENTATIVE_CLASSES):
    block = representative_table[representative_table["spatial_phenotype"].eq(label)]

    for column_index in range(GENES_PER_CLASS):
        axis = axes[row_index, column_index]

        if column_index >= len(block):
            axis.set_visible(False)
            continue

        record = block.iloc[column_index]
        index = int(np.flatnonzero(gene_names == record["gene"])[0])

        for specimen in mouse_samples:
            axis.plot(grid_unit, specimen_curves[specimen][index], color=SPECIES_COLORS["mouse"],
                      lw=0.8, alpha=0.35, ls=":")
        for specimen in human_samples:
            axis.plot(grid_unit, specimen_curves[specimen][index], color=SPECIES_COLORS["human"],
                      lw=0.8, alpha=0.35, ls=":")

        axis.plot(grid_unit, balanced_mouse[index], color=SPECIES_COLORS["mouse"], lw=1.9)
        axis.plot(grid_unit, human_registered[index], color=SPECIES_COLORS["human"], lw=1.9, ls="--")

        axis.set_title(f"{record['gene']}\n"
                       f"A_m={record['mouse_amplitude']:.2f} A_h={record['human_amplitude']:.2f} "
                       f"r={record['shape_corr']:.2f}\n"
                       f"shift={record['best_shift_human_minus_mouse']:+.2f} "
                       f"robust={record['robustness_score']:.2f}",
                       fontsize=7.5)

        if column_index == 0:
            axis.set_ylabel(textwrap.fill(label, 18), fontsize=8)
        if row_index == len(REPRESENTATIVE_CLASSES) - 1:
            axis.set_xlabel("PT position (unit)", fontsize=8)

    for axis in axes[row_index]:
        axis.tick_params(labelsize=7)

axes[0, 0].legend(
    handles=[
        Line2D([], [], color=SPECIES_COLORS["mouse"], lw=1.9, label="mouse (balanced)"),
        Line2D([], [], color=SPECIES_COLORS["human"], lw=1.9, ls="--", label="human (registered)"),
        Line2D([], [], color="0.6", lw=0.8, ls=":", label="individual specimen/slice"),
    ], loc="upper right", fontsize=7,
)

figure.suptitle(
    "Figure 1D - Representative genes per spatial phenotype\n"
    "dotted = each mouse specimen / human slice fitted separately; solid/dashed = equal-weight species curves\n"
    "A_m / A_h = zonation amplitude, r = standardised curve correlation, shift = best displacement (+ = later in human)",
    fontsize=10.5,
)
figure.tight_layout(rect=(0, 0, 1, 0.975))
_save_figure(figure, "representative_gene_spatial_phenotypes.png")
plt.show()

_save_table(representative_table, "representative_genes_per_class.csv")


# %%
# Purpose: figure 1E - shape-only heatmap: within-species z-scored curves, mouse and human blocks.

heatmap_classes = [label for label in REPRESENTATIVE_CLASSES if label != "weak / uncertain zonation"]
GENES_PER_CLASS_HEATMAP = CONFIG["heatmap_genes_per_class"]

heatmap_table = pd.concat(
    [_top_genes_for_class(label, GENES_PER_CLASS_HEATMAP) for label in heatmap_classes],
    ignore_index=True,
)

heatmap_index = np.asarray([int(np.flatnonzero(gene_names == gene)[0])
                            for gene in heatmap_table["gene"]], dtype=int)

mouse_shape = _row_zscore(balanced_mouse[heatmap_index])
human_shape = _row_zscore(human_registered[heatmap_index])

figure, axes = plt.subplots(
    1, 2, figsize=(9.5, max(6.0, 0.21 * len(heatmap_index))), sharey=True
)

for axis, matrix, title in zip(
    axes, [mouse_shape, human_shape], ["Healthy mouse", "Human (registered)"]
):
    image = axis.imshow(
        matrix, aspect="auto", interpolation="nearest", cmap="bwr", vmin=-2.5, vmax=2.5,
        extent=[0, 1, len(heatmap_index) - 0.5, -0.5],
    )
    axis.set_xlabel("PT position (registered, unit): early -> late")
    axis.set_title(title, fontsize=10)

    # Class boundaries.
    boundaries = np.cumsum([len(heatmap_table[heatmap_table["spatial_phenotype"].eq(label)])
                            for label in heatmap_classes])[:-1]
    for boundary in boundaries:
        axis.axhline(boundary - 0.5, color="black", lw=0.6, alpha=0.5)

axes[0].set_yticks(np.arange(len(heatmap_index)))
axes[0].set_yticklabels(
    [f"{gene}   [{label}]" for gene, label in
     zip(heatmap_table["gene"], heatmap_table["spatial_phenotype"])],
    fontsize=6.5,
)

colour_bar = figure.colorbar(image, ax=axes, label="Within-species curve z-score", shrink=0.6)

figure.suptitle(
    "Figure 1E - Spatial organisation only\n"
    "each gene z-scored within each species, so abundance differences cannot drive the pattern\n"
    "rows grouped by spatial phenotype",
    fontsize=10.5,
)
figure.tight_layout(rect=(0, 0, 1, 0.94))
_save_figure(figure, "shape_only_gene_heatmap.png")
plt.show()

_save_table(heatmap_table, "shape_only_heatmap_gene_selection.csv")


# %% [markdown]
# ## 7 - Clustering the residual rewiring, and only that
#
# **Why clustering appears here and not earlier.** The classes in section 4 are predefined and
# interpretable; most of the biology the paper will quote is already named by them. What they cannot
# describe is the *shape* of a divergence that is neither a displacement nor an amplitude change — the
# `complex shape rewiring` class. That class is the one place where unsupervised structure can add
# something, so it is the only place clustering is used, exactly as the question demands.
#
# **What is clustered.** The difference curves `D_g(s) = z_human,g(s) - z_mouse,g(s)` on the registered
# axis: each species standardised within itself first, so a cluster cannot be created by abundance or
# by overall zonation strength, only by *where* the two species' patterns disagree.
#
# **How.** Functional PCA on the difference curves, hierarchical clustering (average linkage, 1 -
# correlation distance) on them, `k` chosen in 4–8 by silhouette score. Every cluster is then named
# from its own mean difference curve - which third of PT gains or loses, whether it is biphasic, and
# whether the difference runs against the mouse gradient - never from a textbook pathway name.
# Stability is checked by recomputing the clustering on each leave-one-specimen / leave-one-slice
# variant and reporting the adjusted Rand index against the primary labels. If the class is too small
# to cluster (`cluster_min_genes`), the section says so and is skipped rather than forced.
#

# %%
# Purpose: clustering of the residual complex-rewiring difference curves, with stability and naming.

from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score, silhouette_score

rewiring_genes = gene_metrics.loc[
    gene_metrics["spatial_phenotype"].eq("complex shape rewiring")
    & ~gene_metrics["axis_basis_gene"]
    & ~gene_metrics["technical_gene"],
    "gene",
].to_numpy()
rewiring_index = np.asarray([gene_lookup_tested[str(gene).upper()] for gene in rewiring_genes], dtype=int)

_report("genes in the complex shape rewiring class", len(rewiring_genes),
        why=f"clustering runs only if there are >= {CONFIG['cluster_min_genes']}")

if len(rewiring_genes) < CONFIG["cluster_min_genes"]:
    print("Too few genes to cluster: the section is skipped, not forced.")
    complex_rewiring_clusters = pd.DataFrame()
    complex_cluster_profiles = pd.DataFrame()
else:
    difference_curves = (
        _row_zscore(human_registered[rewiring_index]) - _row_zscore(balanced_mouse[rewiring_index])
    )

    n_components = int(min(10, difference_curves.shape[0] - 1, difference_curves.shape[1] - 1))
    pca = PCA(n_components=n_components, random_state=CONFIG["cluster_seed"])
    component_scores = pca.transform(difference_curves)

    functional_pca = pd.DataFrame({
        "component": np.arange(1, n_components + 1),
        "explained_variance_ratio": pca.explained_variance_ratio_,
        "cumulative": np.cumsum(pca.explained_variance_ratio_),
    })
    display(functional_pca.round(4))

    correlation_distance = 1.0 - np.corrcoef(difference_curves)
    correlation_distance = np.clip(
        np.where(np.isfinite(correlation_distance), correlation_distance, 1.0), 0.0, 2.0
    )
    np.fill_diagonal(correlation_distance, 0.0)
    linkage_matrix = linkage(squareform(correlation_distance, checks=False), method="average")

    choice_rows = []
    for k in range(CONFIG["cluster_k_min"], CONFIG["cluster_k_max"] + 1):
        candidate = fcluster(linkage_matrix, k, criterion="maxclust")
        if np.unique(candidate).size < 2:
            continue
        sizes = np.bincount(candidate)[1:]
        choice_rows.append({
            "k": k,
            "silhouette": float(silhouette_score(difference_curves, candidate, metric="correlation")),
            "smallest_cluster": int(sizes.min()),
            "largest_cluster": int(sizes.max()),
        })

    k_choice = pd.DataFrame(choice_rows)
    display(k_choice.round(4))

    chosen_k = int(k_choice.loc[k_choice["silhouette"].idxmax(), "k"])
    print(f"chosen k = {chosen_k} (best silhouette in [{CONFIG['cluster_k_min']}, {CONFIG['cluster_k_max']}])")

    cluster_ids = fcluster(linkage_matrix, chosen_k, criterion="maxclust")

    mouse_reference = _nanmean_safe(
        _row_zscore(balanced_mouse[rewiring_index]), axis=0
    )
    mouse_pattern_scale = float(np.ptp(mouse_reference))

    def _name_difference_curve(mean_curve, x_unit, reference_curve, reference_scale):
        """Name a cluster from the morphology of its own mean difference curve."""
        amplitude = float(np.ptp(mean_curve))
        thirds = np.array_split(np.arange(mean_curve.size), 3)
        lobe_means = {
            name: float(np.mean(mean_curve[part]))
            for name, part in zip(("early", "mid", "late"), thirds)
        }
        strongest_lobe = max(lobe_means, key=lambda key: abs(lobe_means[key]))
        strongest_value = lobe_means[strongest_lobe]

        with np.errstate(invalid="ignore"):
            reversal = float(np.corrcoef(mean_curve, reference_curve)[0, 1])
        if np.isfinite(reversal) and reversal <= -0.5:
            return "gradient reversal (human pattern opposes the mouse pattern)"

        positive = max(lobe_means.values())
        negative = min(lobe_means.values())
        if amplitude > 0 and min(abs(positive), abs(negative)) > 0.15 * amplitude:
            return "biphasic rewiring (gain in one third, loss in another)"

        if amplitude <= 0 or amplitude <= 0.10 * max(reference_scale, 1e-9):
            return "weak / flat difference"

        direction = "gain" if strongest_value > 0 else "loss"
        return f"{strongest_lobe}-PT human {direction}"

    profile_rows = []
    names_by_cluster = {}
    for cluster in np.unique(cluster_ids):
        member_rows = np.flatnonzero(cluster_ids == cluster)
        mean_curve = _nanmean_safe(difference_curves[member_rows], axis=0)
        cluster_name = _name_difference_curve(
            mean_curve, grid_unit, mouse_reference, mouse_pattern_scale
        )
        names_by_cluster[int(cluster)] = cluster_name
        profile_rows.append({
            "cluster": int(cluster),
            "cluster_name": cluster_name,
            "n_genes": int(member_rows.size),
            "mean_pc1": float(component_scores[member_rows, 0].mean()),
            "mean_amplitude_log2_ratio": float(
                gene_metrics.loc[gene_metrics["gene"].isin(rewiring_genes[member_rows]),
                                 "amplitude_log2_ratio_human_over_mouse"].mean()),
            "mean_shape_corr": float(
                gene_metrics.loc[gene_metrics["gene"].isin(rewiring_genes[member_rows]),
                                 "shape_corr"].mean()),
            "mean_abs_shift": float(np.nanmean(np.abs(
                gene_metrics.loc[gene_metrics["gene"].isin(rewiring_genes[member_rows]),
                                 "best_shift_human_minus_mouse"]))),
        })

    complex_cluster_profiles = pd.DataFrame(profile_rows).sort_values("n_genes", ascending=False)
    display(complex_cluster_profiles.round(3))

    complex_rewiring_clusters = pd.DataFrame({
        "gene": rewiring_genes,
        "cluster": cluster_ids,
        "cluster_name": [names_by_cluster[int(cluster)] for cluster in cluster_ids],
        "difference_amplitude": np.nanmax(difference_curves, axis=1) - np.nanmin(difference_curves, axis=1),
        "difference_mean": _nanmean_safe(difference_curves, axis=1),
        "functional_pc1": component_scores[:, 0],
        "functional_pc2": component_scores[:, 1] if n_components > 1 else np.nan,
        "functional_pc3": component_scores[:, 2] if n_components > 2 else np.nan,
    }).merge(
        gene_metrics[["gene", "mouse_amplitude", "human_amplitude", "amplitude_log2_ratio_human_over_mouse",
                      "shape_corr", "best_shift_human_minus_mouse", "robustness_score",
                      "level_effect_human_minus_mouse"]],
        on="gene", how="left",
    )

    _save_table(complex_rewiring_clusters, "complex_rewiring_clusters.csv")

    # ----------------------------------------------------------------------------------
    # Stability of the clustering under the specimen/slice variants.
    # ----------------------------------------------------------------------------------

    stability_rows = []
    variant_curve_sets = []
    for name in mouse_samples:
        remaining = [specimen for specimen in mouse_samples if specimen != name]
        variant_curve_sets.append((
            f"leave out mouse specimen {name}",
            specimen_balanced_curves({specimen: specimen_curves[specimen] for specimen in remaining}),
            balanced_human,
        ))
    for name in human_samples:
        remaining = [specimen for specimen in human_samples if specimen != name]
        variant_curve_sets.append((
            f"leave out human slice {name}",
            balanced_mouse,
            specimen_balanced_curves({specimen: specimen_curves[specimen] for specimen in remaining}),
        ))

    for variant_name, variant_mouse, variant_human_raw in variant_curve_sets:
        variant_registration = _register(variant_mouse, variant_human_raw)
        variant_difference = (
            _row_zscore(variant_registration["human_registered"][rewiring_index])
            - _row_zscore(variant_mouse[rewiring_index])
        )
        variant_distance = 1.0 - np.corrcoef(variant_difference)
        variant_distance = np.clip(
            np.where(np.isfinite(variant_distance), variant_distance, 1.0), 0.0, 2.0
        )
        np.fill_diagonal(variant_distance, 0.0)
        variant_linkage = linkage(squareform(variant_distance, checks=False), method="average")
        variant_ids = fcluster(variant_linkage, chosen_k, criterion="maxclust")

        stability_rows.append({
            "variant": variant_name,
            "adjusted_rand_index_vs_primary": float(adjusted_rand_score(cluster_ids, variant_ids)),
            "n_clusters_recovered": int(np.unique(variant_ids).size),
        })

    complex_cluster_stability = pd.DataFrame(stability_rows)
    display(complex_cluster_stability.round(3))
    _save_table(complex_cluster_stability, "complex_rewiring_cluster_stability.csv")

    # ----------------------------------------------------------------------------------
    # Figure: cluster mean difference curves and member genes.
    # ----------------------------------------------------------------------------------

    figure, axes = plt.subplots(1, 2, figsize=(13, 4.6), width_ratios=[1.0, 1.25])

    for cluster in np.unique(cluster_ids):
        member_rows = np.flatnonzero(cluster_ids == cluster)
        mean_curve = _nanmean_safe(difference_curves[member_rows], axis=0)
        axes[0].plot(grid_unit, mean_curve, lw=2.0,
                     label=f"{names_by_cluster[int(cluster)]} (n={member_rows.size})")
    axes[0].axhline(0, color="black", lw=0.8)
    axes[0].set_xlabel("Registered PT position (unit): early -> late")
    axes[0].set_ylabel("Mean difference of z-scored curves\nz_human - z_mouse")
    axes[0].set_title("Cluster mean difference curves", loc="left")
    axes[0].legend(fontsize=7)

    ordered_rows = np.argsort(cluster_ids)
    image = axes[1].imshow(
        difference_curves[ordered_rows], aspect="auto", interpolation="nearest", cmap="bwr",
        vmin=-3, vmax=3, extent=[0, 1, ordered_rows.size - 0.5, -0.5],
    )
    boundaries = np.cumsum(np.bincount(cluster_ids)[1:])[:-1]
    for boundary in boundaries:
        axes[1].axhline(boundary - 0.5, color="black", lw=0.7)
    axes[1].set_xlabel("Registered PT position (unit)")
    axes[1].set_ylabel("Genes, grouped by cluster")
    axes[1].set_title("Difference curves of the clustered genes", loc="left")
    figure.colorbar(image, ax=axes[1], label="z_human - z_mouse", shrink=0.85)

    figure.suptitle("Residual (complex shape rewiring) difference curves\n"
                    "z-scored within each species first, so this is spatial disagreement only",
                    fontsize=10.5)
    figure.tight_layout()
    _save_figure(figure, "complex_rewiring_clusters.png")
    plt.show()


# %% [markdown]
# ## 8 - Which functional programs sit in which spatial phenotype
#
# **Why not mean pathway expression.** Averaging a pathway's member genes along PT discards exactly the
# thing under study: two pathways with identical mean trajectories can have members that all sit in the
# same territory or members split across opposing territories. So the test here is whether a pathway's
# member genes are **over-represented in a spatial phenotype** (or in "excluded"), never whether the
# pathway's mean expression differs between species.
#
# **The test.** Hypergeometric over-representation against the discovery-eligible universe: `N` = tested
# orthologs that could be assigned a phenotype, `K` = the pathway's members inside that universe,
# `n` = genes holding the phenotype, `k` = pathway members holding it. Background is the tested
# ortholog universe, and pathways are restricted to `[pathway_min_members, pathway_max_members]`
# members so that tiny curated sets and transcriptome-wide sets cannot dominate. BH correction is
# applied across the **complete** tested pathway x phenotype family (reported), so the FDRs are not
# per-phenotype claims.
#
# **What counts as evidence.** A fold enrichment with a family-wide FDR below `enrichment_fdr`, plus the
# leading genes being named in the table so a reader can see whether the signal is a coherent programme
# or one strong gene. Redundant databases overlap heavily, so the heatmap collapses redundancy groups
# (03's `summarize_pathway_redundancy`, overlap >= 0.6) and shows one representative per group.
#
# **Limits.** This is exploratory: the units of replication are specimens, and enrichment here is a
# property of *this* gene set in *this* cohort, not a population claim. A phenotype with few genes makes
# a weak universe to test against, so `n_universe_with_phenotype` is reported in every row.
#

# %%
# Purpose: pathway x spatial-phenotype over-representation, with family-wide BH correction.

import ast


def _parse_member_list(value):
    """03 stores each pathway's members as a Python or ';'-delimited string."""
    if isinstance(value, (list, tuple)):
        return [str(gene) for gene in value]
    if pd.isna(value):
        return []

    text = str(value)
    try:
        parsed = ast.literal_eval(text)
        if isinstance(parsed, (list, tuple)):
            return [str(gene) for gene in parsed]
    except (ValueError, SyntaxError):
        pass

    return [gene.strip() for gene in text.split(";") if gene.strip()]


if not PATHWAY_MEMBERSHIP_PATH.exists():
    raise FileNotFoundError(
        f"{PATHWAY_MEMBERSHIP_PATH} is missing - run 03_human_vs_healthy_mouse.ipynb first."
    )

pathway_coverage = pd.read_csv(PATHWAY_MEMBERSHIP_PATH)
pathway_coverage["member_list"] = pathway_coverage["genes_present"].map(_parse_member_list)

if "retained" in pathway_coverage.columns:
    membership = pathway_coverage[pathway_coverage["retained"].astype(bool)].copy()
else:
    membership = pathway_coverage.copy()
    print("WARNING: no 'retained' column in the membership file - using every row.")

membership["universe_members"] = membership["member_list"].map(
    lambda genes: sorted({
        gene_lookup_tested[str(gene).upper()] for gene in genes
        if str(gene).upper() in gene_lookup_tested
    })
)
membership["n_universe_members"] = membership["universe_members"].map(len)

_report("pathways retained by 03", len(membership))
_report(
    "pathways inside the tested size range",
    int(membership["n_universe_members"].between(
        CONFIG["pathway_min_members"], CONFIG["pathway_max_members"]).sum()),
    len(membership),
    why=f"{CONFIG['pathway_min_members']} <= members in the tested universe "
        f"<= {CONFIG['pathway_max_members']}",
)

universe_index = np.flatnonzero(discovery_eligible.to_numpy())
universe_size = int(universe_index.size)
universe_mask = np.zeros(len(gene_names), dtype=bool)
universe_mask[universe_index] = True

phenotype_by_row = gene_metrics["spatial_phenotype"].to_numpy()
discovery_score_by_row = np.nan_to_num(gene_metrics["spatial_discovery_score"].to_numpy())

universe_phenotype_counts = pd.Series(phenotype_by_row[universe_index]).value_counts()
print("Phenotype composition of the tested universe:")
display(universe_phenotype_counts.rename("n_genes").to_frame())

# --------------------------------------------------------------------------------------
# The over-representation tests.
# --------------------------------------------------------------------------------------

pathway_tested = membership[
    membership["n_universe_members"].between(
        CONFIG["pathway_min_members"], CONFIG["pathway_max_members"])
].copy()

enrichment_rows = []

for row in pathway_tested.itertuples():
    members = np.asarray(
        [index for index in row.universe_members if universe_mask[index]], dtype=int
    )
    if members.size < CONFIG["pathway_min_members"]:
        continue

    member_phenotypes = phenotype_by_row[members]

    for label in PHENOTYPE_CLASSES:
        n_with_phenotype = int(universe_phenotype_counts.get(label, 0))
        if n_with_phenotype == 0:
            continue

        observed = int((member_phenotypes == label).sum())
        expected = n_with_phenotype * members.size / universe_size
        overlap = members[member_phenotypes == label]
        leading = overlap[np.argsort(-discovery_score_by_row[overlap])][:10]

        enrichment_rows.append({
            "library": row.library,
            "pathway": row.pathway,
            "spatial_phenotype": label,
            "n_pathway_members_in_universe": int(members.size),
            "n_universe_with_phenotype": n_with_phenotype,
            "observed": observed,
            "expected": expected,
            "fold_enrichment": (observed / expected) if expected > 0 else np.nan,
            "p_value": float(hypergeom.sf(observed - 1, universe_size, members.size, n_with_phenotype)),
            "leading_genes": ";".join(gene_names[leading]),
        })

pathway_enrichment = pd.DataFrame(enrichment_rows)
pathway_enrichment["fdr"] = bh_adjust(pathway_enrichment["p_value"].to_numpy())
pathway_enrichment["minus_log10_fdr"] = -np.log10(
    np.maximum(pathway_enrichment["fdr"], 1e-300)
)
# Sign encodes direction of the effect: positive = over-represented, negative = depleted.
pathway_enrichment["signed_strength"] = np.where(
    pathway_enrichment["fold_enrichment"] >= 1.0,
    pathway_enrichment["minus_log10_fdr"],
    -pathway_enrichment["minus_log10_fdr"],
)
pathway_enrichment = (
    pathway_enrichment
    .sort_values(["fdr", "fold_enrichment"], ascending=[True, False])
    .reset_index(drop=True)
)

_report(
    "tested pathway x phenotype pairs", len(pathway_enrichment),
    int(len(pathway_tested) * len(PHENOTYPE_CLASSES)),
    why="BH was applied across this complete family",
)
_report(
    "pairs with FDR < enrichment_fdr",
    int((pathway_enrichment["fdr"] < CONFIG["enrichment_fdr"]).sum()), len(pathway_enrichment),
)

display(pathway_enrichment.head(20).round(4))
_save_table(pathway_enrichment, "pathway_spatial_phenotype_enrichment.csv")


# %%
# Purpose: the pathway x phenotype heatmap, collapsed to non-redundant representative terms.

redundancy_pairs, redundancy_groups = summarize_pathway_redundancy(
    membership.assign(genes_present=membership["member_list"]),
    gene_column="genes_present",
    overlap_threshold=0.6,
)

group_of_pathway = {
    (row.library, row.pathway): int(row.redundancy_group)
    for row in redundancy_groups.itertuples()
}
group_size_of_pathway = {
    (row.library, row.pathway): int(row.group_size)
    for row in redundancy_groups.itertuples()
}

_report("pathway overlap pairs at >= 60% of the smaller set", len(redundancy_pairs))
_report("redundancy groups", int(redundancy_groups["redundancy_group"].nunique()), len(redundancy_groups),
        why="one representative per group is plotted")

significant_pairs = pathway_enrichment[
    pathway_enrichment["fdr"] < CONFIG["enrichment_fdr"]
].sort_values("fdr").reset_index(drop=True)

representatives = []
seen_groups = set()
MAX_HEATMAP_PATHWAYS = 30

for row in significant_pairs.itertuples():
    group = group_of_pathway.get((row.library, row.pathway), None)
    if group is not None and group in seen_groups:
        continue
    if group is not None:
        seen_groups.add(group)
    representatives.append(row)
    if len(representatives) >= MAX_HEATMAP_PATHWAYS:
        break

heatmap_frame = pd.DataFrame([{
    "library": row.library,
    "pathway": row.pathway,
    "spatial_phenotype": row.spatial_phenotype,
    "signed_strength": row.signed_strength,
    "fdr": row.fdr,
    "fold_enrichment": row.fold_enrichment,
    "redundancy_group": group_of_pathway.get((row.library, row.pathway), np.nan),
    "group_size": group_size_of_pathway.get((row.library, row.pathway), np.nan),
    "leading_genes": row.leading_genes,
} for row in representatives])

if heatmap_frame.empty:
    print(f"No pathway x phenotype pair reached FDR < {CONFIG['enrichment_fdr']}: "
          "the heatmap is skipped rather than drawn empty. The full enrichment table is still saved.")
else:
    # A pathway can be significant for more than one phenotype; show its strongest per phenotype.
    heatmap_matrix = (
        heatmap_frame
        .pivot_table(index=["library", "pathway"], columns="spatial_phenotype",
                     values="signed_strength", aggfunc="max")
        .reindex(columns=[label for label in PHENOTYPE_CLASSES
                          if label in set(heatmap_frame["spatial_phenotype"])])
    )
    heatmap_matrix = heatmap_matrix.reindex(
        heatmap_matrix.abs().max(axis=1).sort_values(ascending=False).index
    )

    heatmap_labels = [
        f"{pathway}  [{library}]"
        + ("" if pd.isna(group_of_pathway.get((library, pathway)))
           else f" (group of {int(group_size_of_pathway[(library, pathway)])})")
        for library, pathway in heatmap_matrix.index
    ]
    scale = float(np.nanmax(np.abs(heatmap_matrix.to_numpy())))

    figure, axis = plt.subplots(
        1, 1,
        figsize=(1.15 * heatmap_matrix.shape[1] + 8.0, max(5.0, 0.32 * heatmap_matrix.shape[0])),
    )

    image = axis.imshow(
        heatmap_matrix.to_numpy(), aspect="auto", interpolation="nearest", cmap="PuOr_r",
        vmin=-scale, vmax=scale,
    )

    axis.set_xticks(np.arange(heatmap_matrix.shape[1]))
    axis.set_xticklabels([textwrap.fill(str(label), 18) for label in heatmap_matrix.columns],
                         rotation=45, ha="right", fontsize=8)
    axis.set_yticks(np.arange(heatmap_matrix.shape[0]))
    axis.set_yticklabels([textwrap.fill(str(label), 46) for label in heatmap_labels], fontsize=7.5)

    for row_index in range(heatmap_matrix.shape[0]):
        for column_index in range(heatmap_matrix.shape[1]):
            value = heatmap_matrix.to_numpy()[row_index, column_index]
            if np.isfinite(value):
                axis.text(column_index, row_index, f"{value:.1f}", ha="center", va="center", fontsize=6.5)

    colour_bar = figure.colorbar(image, ax=axis, shrink=0.7)
    colour_bar.set_label("Signed -log10(FDR):  + enriched,  - depleted\n(one term per redundancy group)")

    axis.set_title(
        "Pathway spatial-phenotype enrichment (non-redundant representatives)\n"
        "hypergeometric against the tested-ortholog universe, BH across the complete "
        "pathway x phenotype family",
        loc="left", fontsize=10,
    )
    figure.tight_layout()
    _save_figure(figure, "pathway_spatial_phenotype_heatmap.png")
    plt.show()

    _save_table(heatmap_frame, "pathway_spatial_phenotype_heatmap_terms.csv")


# %% [markdown]
# ## 9 - Inspecting whole pathways along the coordinate
#
# **Why a pathway-level heatmap.** Enrichment tells you that a pathway's genes are over-represented in a
# phenotype; it does not show whether the pathway is internally coherent. Plotting its members as
# within-species z-scored trajectories, mouse and human side by side, answers the questions the
# phenotype framework is meant to settle: does the programme keep its shape, shift, gain or lose
# zonation, or split into opposing sub-programmes?
#
# **Which pathways.** The themes carried by earlier results (peroxisome, fatty acid metabolism, lipid
# metabolism, bile acid metabolism, oxidative phosphorylation/TCA, TNF-alpha/NF-kB, hypoxia) are **looked
# up, not assumed**: for each theme the best-supported matching terms are found in the membership file,
# and the table printed first records the matched terms, their member counts and their strongest
# phenotype enrichment. A theme with no significant enrichment is still shown - that is a result in its
# own right - but its panel says so.
#
# **How the panels are built.** Members are ordered by positional centroid, each gene's curve is z-scored
# within each species, mouse and human blocks are drawn next to each other on the same gene order, and
# every gene row is annotated with its own spatial phenotype. Rows are genes, never samples.
# The candidate themes are shown first, then the strongest non-redundant enrichment findings
# overall in the same layout, so a pathway outside the theme list can still appear.
#

# %%
# Purpose: gene x pseudospace panels for the candidate pathway themes, with per-gene phenotypes.

PATHWAY_THEMES = {
    "Peroxisome": ["peroxisom"],
    "Fatty acid metabolism": ["fatty acid"],
    "Lipid metabolism": ["lipid"],
    "Bile acid metabolism": ["bile acid"],
    "Oxidative phosphorylation / TCA": [
        "oxidative phosphorylation", "citric acid", "tca cycle", "respiratory electron transport",
    ],
    "TNF-alpha / NF-kB": ["tnf", "nf-kb", "nfkb"],
    "Hypoxia": ["hypoxia"],
}

MAX_MEMBERS_PER_THEME = 25

theme_rows = []
theme_selections = {}

for theme, patterns in PATHWAY_THEMES.items():
    matches = membership[
        membership["pathway"].astype(str).str.lower().map(
            lambda name: any(pattern in name for pattern in patterns)
        )
        & membership["n_universe_members"].ge(CONFIG["pathway_min_members"])
    ].sort_values("n_universe_members", ascending=False)

    if not len(matches):
        theme_rows.append({
            "theme": theme, "matched_terms": "none", "chosen_term": "",
            "n_members_in_universe": 0, "best_phenotype": "", "best_fdr": np.nan,
            "best_fold_enrichment": np.nan, "leading_genes": "",
        })
        continue

    chosen = matches.iloc[0]
    term_enrichment = pathway_enrichment[
        pathway_enrichment["pathway"].eq(chosen["pathway"])
        & pathway_enrichment["library"].eq(chosen["library"])
    ].sort_values("fdr")

    best = term_enrichment.iloc[0] if len(term_enrichment) else None

    theme_selections[theme] = {
        "library": chosen["library"],
        "pathway": chosen["pathway"],
        "members": [gene_lookup_tested[str(gene).upper()] for gene in chosen["member_list"]
                    if str(gene).upper() in gene_lookup_tested],
    }

    theme_rows.append({
        "theme": theme,
        "matched_terms": "; ".join(matches["pathway"].astype(str).head(6)),
        "chosen_term": f"{chosen['library']}: {chosen['pathway']}",
        "n_members_in_universe": int(chosen["n_universe_members"]),
        "best_phenotype": "" if best is None else best["spatial_phenotype"],
        "best_fdr": np.nan if best is None else float(best["fdr"]),
        "best_fold_enrichment": np.nan if best is None else float(best["fold_enrichment"]),
        "leading_genes": "" if best is None else best["leading_genes"],
    })

theme_summary = pd.DataFrame(theme_rows)
display(theme_summary.round(4))
_save_table(theme_summary, "pathway_theme_summary.csv")

# --------------------------------------------------------------------------------------
# One row per theme: mouse block and human block, members ordered by positional centroid.
# --------------------------------------------------------------------------------------

plotted_themes = [theme for theme in PATHWAY_THEMES if theme in theme_selections]
height_ratios = []
for theme in plotted_themes:
    members = theme_selections[theme]["members"]
    selected = members[:MAX_MEMBERS_PER_THEME]
    height_ratios.append(max(len(selected), 3))

def _ordered_member_blocks(members, limit):
    """Limit members, order them by human positional centroid, and z-score each species' block.

    Shared by the candidate-theme panels and the strongest-findings panels so the two figures cannot
    drift apart in ordering, scaling or colour limits.
    """
    selected = np.asarray(members[:limit], dtype=int)
    centroids = _row_weighted_centroid(human_registered[selected], grid_unit)
    ordered = selected[np.argsort(np.nan_to_num(centroids, nan=1.0))]
    return ordered, _row_zscore(balanced_mouse[ordered]), _row_zscore(human_registered[ordered])


if not plotted_themes:
    print("No candidate theme matched a pathway in the membership file: the panel figure is skipped.")
else:
    figure, axes = plt.subplots(
        len(plotted_themes), 2,
        figsize=(10.5, 0.30 * sum(height_ratios) + 1.6),
        gridspec_kw={"height_ratios": height_ratios, "hspace": 0.75, "wspace": 0.05},
    )
    axes = np.atleast_2d(axes)

    for row_index, theme in enumerate(plotted_themes):
        selection = theme_selections[theme]
        ordered, mouse_block, human_block = _ordered_member_blocks(
            selection["members"], MAX_MEMBERS_PER_THEME
        )

        for column_index, (block, title) in enumerate([
            (mouse_block, "mouse"), (human_block, "human (registered)"),
        ]):
            axis = axes[row_index, column_index]
            image = axis.imshow(
                block, aspect="auto", interpolation="nearest", cmap="bwr", vmin=-2.5, vmax=2.5,
                extent=[0, 1, len(ordered) - 0.5, -0.5],
            )
            axis.set_yticks(np.arange(len(ordered)))
            if column_index == 0:
                axis.set_yticklabels(
                    [f"{gene_names[index]}  [{phenotype_by_row[index]}]" for index in ordered],
                    fontsize=6.5,
                )
                axis.set_ylabel(textwrap.fill(theme, 26), fontsize=8.5)
            else:
                axis.set_yticklabels([])

            if row_index == 0:
                axis.set_title(title, fontsize=9)
            if row_index == len(plotted_themes) - 1:
                axis.set_xlabel("PT position (unit): early -> late", fontsize=8)
            axis.tick_params(labelsize=7)

        best_fdr = theme_summary.loc[theme_summary["theme"].eq(theme), "best_fdr"].iloc[0]
        best_label = theme_summary.loc[theme_summary["theme"].eq(theme), "best_phenotype"].iloc[0]
        axes[row_index, 0].text(
            0.01, -0.30,
            f"{selection['library']}: {selection['pathway']}   "
            f"(strongest enrichment: {best_label or 'none'}, FDR={best_fdr:.3g})"
            if np.isfinite(best_fdr) else
            f"{selection['library']}: {selection['pathway']}   (no significant phenotype enrichment)",
            transform=axes[row_index, 0].transAxes, fontsize=7.5, va="top",
        )

    figure.suptitle(
        "Candidate pathway themes along the registered coordinate\n"
        "each gene z-scored within species; rows ordered by human positional centroid; "
        "row labels give the gene's spatial phenotype",
        fontsize=10.5,
    )
    figure.tight_layout(rect=(0, 0, 1, 0.965))
    _save_figure(figure, "representative_pathway_spatial_heatmaps.png")
    plt.show()


# %%
# Purpose: the same panels for the strongest pathway findings overall, not only the candidate themes.

STRONGEST_PATHWAYS_TO_PLOT = 6
MAX_MEMBERS_PER_STRONGEST = 20

if heatmap_frame.empty:
    print(f"No pathway x phenotype pair reached FDR < {CONFIG['enrichment_fdr']}: "
          "the strongest-findings panels are skipped.")
    strongest_pathway_findings = pd.DataFrame()
else:
    strongest_rows = (
        heatmap_frame
        .sort_values("fdr")
        .drop_duplicates(subset=["pathway"])
        .head(STRONGEST_PATHWAYS_TO_PLOT)
        .reset_index(drop=True)
    )

    strongest_selections = []
    for row in strongest_rows.itertuples():
        match = membership[
            membership["library"].eq(row.library) & membership["pathway"].eq(row.pathway)
        ]
        if not len(match):
            continue
        members = [
            gene_lookup_tested[str(gene).upper()] for gene in match.iloc[0]["member_list"]
            if str(gene).upper() in gene_lookup_tested
        ]
        if len(members) < 2:
            continue
        strongest_selections.append({
            "library": row.library,
            "pathway": row.pathway,
            "spatial_phenotype": row.spatial_phenotype,
            "fdr": row.fdr,
            "fold_enrichment": row.fold_enrichment,
            "n_members_in_universe": len(members),
            "leading_genes": row.leading_genes,
            "members": members,
        })

    strongest_pathway_findings = pd.DataFrame(
        [{key: value for key, value in entry.items() if key != "members"}
         for entry in strongest_selections]
    )
    display(strongest_pathway_findings.round(4))
    _save_table(strongest_pathway_findings, "strongest_pathway_findings.csv")

    if not strongest_selections:
        print("No plottable pathway in the strongest-findings selection.")
    else:
        height_ratios = [max(min(len(entry["members"]), MAX_MEMBERS_PER_STRONGEST), 3)
                         for entry in strongest_selections]

        figure, axes = plt.subplots(
            len(strongest_selections), 2,
            figsize=(10.5, 0.30 * sum(height_ratios) + 1.6),
            gridspec_kw={"height_ratios": height_ratios, "hspace": 0.85, "wspace": 0.05},
        )
        axes = np.atleast_2d(axes)

        for row_index, entry in enumerate(strongest_selections):
            ordered, mouse_block, human_block = _ordered_member_blocks(
                entry["members"], MAX_MEMBERS_PER_STRONGEST
            )

            for column_index, (block, title) in enumerate([
                (mouse_block, "mouse"), (human_block, "human (registered)"),
            ]):
                axis = axes[row_index, column_index]
                axis.imshow(
                    block, aspect="auto", interpolation="nearest", cmap="bwr", vmin=-2.5, vmax=2.5,
                    extent=[0, 1, len(ordered) - 0.5, -0.5],
                )
                axis.set_yticks(np.arange(len(ordered)))
                if column_index == 0:
                    axis.set_yticklabels(
                        [f"{gene_names[index]}  [{phenotype_by_row[index]}]" for index in ordered],
                        fontsize=6.5,
                    )
                else:
                    axis.set_yticklabels([])
                if row_index == 0:
                    axis.set_title(title, fontsize=9)
                if row_index == len(strongest_selections) - 1:
                    axis.set_xlabel("PT position (unit): early -> late", fontsize=8)
                axis.tick_params(labelsize=7)

            axes[row_index, 0].text(
                0.01, -0.32,
                f"{entry['library']}: {entry['pathway']}   ->   {entry['spatial_phenotype']}   "
                f"(fold={entry['fold_enrichment']:.2f}, FDR={entry['fdr']:.3g}, "
                f"{entry['n_members_in_universe']} members shown {len(ordered)})",
                transform=axes[row_index, 0].transAxes, fontsize=7.5, va="top",
            )

        figure.suptitle(
            "Strongest pathway findings along the registered coordinate\n"
            "non-redundant terms, each gene z-scored within species and row-labelled with its own "
            "spatial phenotype",
            fontsize=10.5,
        )
        figure.tight_layout(rect=(0, 0, 1, 0.955))
        _save_figure(figure, "strongest_pathway_spatial_heatmaps.png")
        plt.show()


# %% [markdown]
# ## 10 - Conventional whole-PT signal versus continuous spatial signal
#
# **The comparison this notebook exists to make.** A whole-PT comparison asks only about the vertical
# offset of the two species' expression; the continuous analysis asks about amplitude, position and
# shape. A gene or pathway can be strong in one and absent in the other, and the two-by-two of those
# outcomes is the methodological result:
#
# | | high spatial divergence | low spatial divergence |
# | --- | --- | --- |
# | **high whole-PT effect** | both abundance and spatial divergence | mostly a level difference |
# | **low whole-PT effect** | **continuous-only biology** | little species difference |
#
# **How the two axes are defined.** Conventional signal = `|level_effect_human_minus_mouse|` from the
# balanced curves (the same quantity a whole-PT DE would report; when 05's pseudobulk whole-PT log fold
# change is present it is shown alongside as a continuity check, not as the primary). Spatial signal =
# `spatial_discovery_score`, which is built only from amplitude change, displacement and shape
# divergence, and therefore cannot be moved by a level difference.
#
# **Descriptive, by construction.** No p-values are used on the conventional axis: with one human donor
# a species-wide test would not be calibrated. The quadrant counts and pathway table are effect-size
# statements that hold for this cohort.
#
# **What counts as evidence.** For a pathway to count as *continuous-only*, its member genes must be
# enriched in a spatial phenotype at family-wide FDR < `enrichment_fdr` while its conventional effect
# sits below the median of all tested pathways. Both conditions are reported in the table.
#

# %%
# Purpose: the gene-level conventional-versus-continuous comparison and its quadrant figure.

conventional_level = gene_metrics["level_effect_human_minus_mouse"].abs()
spatial_divergence = gene_metrics["spatial_discovery_score"]
conventional_percentile = conventional_level.rank(pct=True)
spatial_percentile = spatial_divergence.rank(pct=True)

# Continuity check: 05's pseudobulk whole-PT log fold change, if that notebook has been run. It is a
# sanity link between the conventional axis here and the one 05 reports, not an input.
WHOLE_PT_LOGFOLD_PATH = UPSTREAM_DIR / "gene_rewiring" / "whole_PT_gene_logFC.csv"

if WHOLE_PT_LOGFOLD_PATH.exists():
    whole_pt_logfold = pd.read_csv(WHOLE_PT_LOGFOLD_PATH)
    continuity = gene_metrics[["gene", "level_effect_human_minus_mouse"]].merge(
        whole_pt_logfold[["gene", "log2fc_human_vs_mouse"]], on="gene", how="inner"
    ).dropna()
    print(f"05's whole-PT pseudobulk log2FC versus this notebook's level effect "
          f"({len(continuity):,} shared genes): "
          f"spearman={spearmanr(continuity['level_effect_human_minus_mouse'], continuity['log2fc_human_vs_mouse']).correlation:.3f}")
else:
    print(f"{WHOLE_PT_LOGFOLD_PATH.name} not found - continuity check against 05 skipped")

quadrant_frame = gene_metrics[discovery_eligible].dropna(
    subset=["level_effect_human_minus_mouse"]
).copy()
quadrant_frame["conventional_abs"] = quadrant_frame["level_effect_human_minus_mouse"].abs()

conventional_median = float(np.nanmedian(conventional_level[discovery_eligible]))
spatial_median = float(np.nanmedian(spatial_divergence[discovery_eligible]))

quadrant_labels = {
    "high conventional / high spatial":
        quadrant_frame["conventional_abs"].ge(conventional_median)
        & quadrant_frame["spatial_discovery_score"].ge(spatial_median),
    "high conventional / low spatial":
        quadrant_frame["conventional_abs"].ge(conventional_median)
        & quadrant_frame["spatial_discovery_score"].lt(spatial_median),
    "low conventional / high spatial (continuous-only)":
        quadrant_frame["conventional_abs"].lt(conventional_median)
        & quadrant_frame["spatial_discovery_score"].ge(spatial_median),
    "low conventional / low spatial":
        quadrant_frame["conventional_abs"].lt(conventional_median)
        & quadrant_frame["spatial_discovery_score"].lt(spatial_median),
}

quadrant_frame["conventional_vs_continuous_quadrant"] = np.select(
    list(quadrant_labels.values()), list(quadrant_labels), default="unassigned"
)

quadrant_counts = pd.DataFrame([
    {
        "quadrant": name,
        "n_genes": int(np.asarray(mask).sum()),
        "share": float(np.asarray(mask).mean()),
        "most_common_spatial_phenotype": (
            str(quadrant_frame.loc[np.asarray(mask), "spatial_phenotype"].mode().iloc[0])
            if np.asarray(mask).any() else ""
        ),
    }
    for name, mask in quadrant_labels.items()
])
display(quadrant_counts.round(4))
_save_table(quadrant_frame, "conventional_vs_continuous_gene_quadrants.csv")

figure, axis = plt.subplots(1, 1, figsize=(8.2, 6.2))

for label in PHENOTYPE_CLASSES:
    block = quadrant_frame[quadrant_frame["spatial_phenotype"].eq(label)]
    if not len(block):
        continue
    axis.scatter(
        block["conventional_abs"], block["spatial_discovery_score"],
        s=6, alpha=0.5, linewidths=0, color=PHENOTYPE_COLORS[label], label=f"{label} (n={len(block):,})",
    )

axis.axvline(conventional_median, color="black", lw=0.9, ls="--")
axis.axhline(spatial_median, color="black", lw=0.9, ls="--")
axis.set_yscale("symlog", linthresh=0.2)
axis.set_xlabel("Conventional whole-PT effect: |level effect| (lognorm)")
axis.set_ylabel("Continuous spatial divergence score")
axis.set_title("Figure: conventional whole-PT effect versus continuous spatial signal\n"
               "upper left = continuous-only biology; dashed lines = medians",
               loc="left", fontsize=10)
axis.legend(fontsize=6.5, loc="upper right")

for row in quadrant_frame.nlargest(8, "spatial_discovery_score").itertuples():
    axis.annotate(row.gene, (row.conventional_abs, row.spatial_discovery_score),
                  xytext=(3, 3), textcoords="offset points", fontsize=7)

figure.tight_layout()
_save_figure(figure, "conventional_vs_continuous_genes.png")
plt.show()


# %%
# Purpose: the pathway-level continuous-only table and its figure.

pathway_conventional_effect = {}
pathway_member_stats = {}

for row in pathway_tested.itertuples():
    members = np.asarray([index for index in row.universe_members if universe_mask[index]], dtype=int)
    if members.size == 0:
        continue

    label = f"{row.library}: {row.pathway}"
    pathway_conventional_effect[label] = float(np.nanmedian(conventional_level.to_numpy()[members]))
    pathway_member_stats[label] = {
        "n_members_in_universe": int(members.size),
        "median_member_abs_level_effect": float(np.nanmedian(conventional_level.to_numpy()[members])),
        "median_member_amplitude_log2_ratio": float(np.nanmedian(
            np.abs(gene_metrics["amplitude_log2_ratio_human_over_mouse"].to_numpy()[members]))),
        "median_member_abs_shift": float(np.nanmedian(np.abs(
            gene_metrics["best_shift_human_minus_mouse"].to_numpy()[members]))),
        "median_member_shape_corr": float(np.nanmedian(
            gene_metrics["shape_corr"].to_numpy()[members])),
        "median_member_spatial_score": float(np.nanmedian(
            spatial_divergence.to_numpy()[members])),
    }

conventional_pathway_ranks = pd.Series(pathway_conventional_effect).rank(pct=True)

pathway_spatial_rewiring_summary = pd.DataFrame([
    {
        "library": row.library,
        "pathway": row.pathway,
        "spatial_phenotype": row.spatial_phenotype,
        "fold_enrichment": row.fold_enrichment,
        "fdr": row.fdr,
        "leading_genes": row.leading_genes,
        "conventional_effect": pathway_conventional_effect.get(f"{row.library}: {row.pathway}", np.nan),
        "conventional_percentile": conventional_pathway_ranks.get(f"{row.library}: {row.pathway}", np.nan),
        **pathway_member_stats.get(f"{row.library}: {row.pathway}", {}),
    }
    for row in pathway_enrichment.itertuples()
]).sort_values(["fdr", "fold_enrichment"], ascending=[True, False]).reset_index(drop=True)

_save_table(pathway_spatial_rewiring_summary, "pathway_spatial_rewiring_summary.csv")

significant_pathways = pathway_spatial_rewiring_summary[
    pathway_spatial_rewiring_summary["fdr"] < CONFIG["enrichment_fdr"]
].copy()

best_per_pathway = (
    significant_pathways
    .sort_values("fdr")
    .groupby(["library", "pathway"], as_index=False)
    .first()
)

continuous_only_pathways = best_per_pathway[
    best_per_pathway["conventional_percentile"] < 0.5
].sort_values("fdr").rename(columns={
    "spatial_phenotype": "dominant_spatial_phenotype",
    "fold_enrichment": "enrichment_fold_enrichment",
    "fdr": "enrichment_fdr",
    "median_member_amplitude_log2_ratio": "median_member_abs_amplitude_log2_ratio",
})[[
    "pathway", "library", "conventional_effect", "conventional_percentile",
    "dominant_spatial_phenotype", "enrichment_fdr", "enrichment_fold_enrichment",
    "median_member_abs_amplitude_log2_ratio", "median_member_abs_shift",
    "median_member_shape_corr", "leading_genes",
]]

_report("pathways with a significant spatial enrichment", len(best_per_pathway), len(pathway_tested))
_report("of those, below the median conventional effect (continuous-only)",
        len(continuous_only_pathways), len(best_per_pathway),
        why="a whole-PT level comparison alone would have missed them")

_save_table(continuous_only_pathways, "continuous_only_pathways.csv")

top_spatial = best_per_pathway.nlargest(12, "median_member_spatial_score")

if top_spatial.empty:
    print("No pathway reached the enrichment threshold: the pathway panel is skipped.")
else:
    figure, axis = plt.subplots(1, 1, figsize=(7.5, 5.4))

    positions = np.arange(len(top_spatial))
    axis.barh(
        positions, top_spatial["median_member_spatial_score"],
        color=[PHENOTYPE_COLORS.get(label, "#888888")
               for label in top_spatial["dominant_spatial_phenotype"]],
        alpha=0.9,
    )
    axis.set_yticks(positions)
    axis.set_yticklabels([textwrap.fill(str(row.pathway), 40) for row in top_spatial.itertuples()],
                         fontsize=7.5)
    axis.set_xlabel("Median member continuous spatial divergence score")
    axis.set_title("Pathways carrying the strongest continuous spatial signal\n"
                   "bar colour = dominant spatial phenotype", loc="left", fontsize=10)

    figure.suptitle("Figure: continuous spatial signal that a whole-PT comparison does not rank\n"
                    "descriptive effect sizes for this cohort (2 mouse specimens, 1 human donor)",
                    fontsize=10.5)
    figure.tight_layout()
    _save_figure(figure, "conventional_vs_continuous_pathway_signal.png")
    plt.show()


# %% [markdown]
# ## 11 - Spatially conserved functions
#
# **Why this section exists.** Divergence is the more eye-catching result, but the more transferable
# statement about mouse-to-human translation is often the opposite one: these are the proximal-tubule
# functions whose spatial organisation is *preserved*. This section asks which pathways are enriched for
# `conserved zonation` — patterned in both species, reproducible within both, strong shape agreement,
# little displacement and no major amplitude change — and contrasts them with the pathways that are
# enriched for any divergent phenotype.
#
# **What counts as evidence.** A pathway whose member genes are enriched for conserved zonation at
# family-wide FDR < `enrichment_fdr`, reported with its member count and the median shape correlation of
# its members, so that "conserved" is not a label applied to a large set with a few coherent genes. The
# contrast figure puts conserved enrichment against divergent enrichment on the same axes: a pathway in
# the upper-left quadrant keeps its spatial programme; one in the lower-right has been remodelled.
#
# **Limits.** Conservation here means conservation of the *spatial programme* in this cohort, not of
# expression level (level is deliberately excluded from the conserved class) and not of function.
#

# %%
# Purpose: pathways enriched for conserved zonation, contrasted with pathways enriched for divergence.

DIVERGENT_PHENOTYPES = [
    "weaker zonation in human",
    "stronger zonation in human",
    "mouse-zonated / human-flat",
    "human-zonated / mouse-flat",
    "shifted earlier in human",
    "shifted later in human",
    "gradient inversion",
    "complex shape rewiring",
]

conserved_enrichment = (
    pathway_enrichment[pathway_enrichment["spatial_phenotype"].eq("conserved zonation")]
    .sort_values("fdr").reset_index(drop=True)
)
conserved_enrichment["median_member_shape_corr"] = [
    pathway_member_stats.get(f"{row.library}: {row.pathway}", {}).get("median_member_shape_corr", np.nan)
    for row in conserved_enrichment.itertuples()
]
conserved_enrichment["median_member_abs_amplitude_log2_ratio"] = [
    pathway_member_stats.get(f"{row.library}: {row.pathway}", {}).get(
        "median_member_amplitude_log2_ratio", np.nan)
    for row in conserved_enrichment.itertuples()
]
conserved_enrichment["median_member_abs_shift"] = [
    pathway_member_stats.get(f"{row.library}: {row.pathway}", {}).get("median_member_abs_shift", np.nan)
    for row in conserved_enrichment.itertuples()
]

spatially_conserved_pathways = conserved_enrichment[
    conserved_enrichment["fdr"] < CONFIG["enrichment_fdr"]
][[
    "library", "pathway", "n_pathway_members_in_universe", "observed", "expected",
    "fold_enrichment", "p_value", "fdr", "median_member_shape_corr",
    "median_member_abs_amplitude_log2_ratio", "median_member_abs_shift", "leading_genes",
]]

display(spatially_conserved_pathways.head(20).round(4))
_report("pathways enriched for conserved zonation", len(spatially_conserved_pathways), len(conserved_enrichment),
        why=f"family-wide FDR < {CONFIG['enrichment_fdr']}")
_save_table(spatially_conserved_pathways, "spatially_conserved_pathways.csv")

divergent_enrichment = (
    pathway_enrichment[pathway_enrichment["spatial_phenotype"].isin(DIVERGENT_PHENOTYPES)]
    .sort_values("fdr")
    .groupby(["library", "pathway"], as_index=False)
    .first()[["library", "pathway", "spatial_phenotype", "fold_enrichment", "fdr",
              "n_pathway_members_in_universe", "leading_genes"]]
    .rename(columns={
        "spatial_phenotype": "dominant_divergent_phenotype",
        "fold_enrichment": "divergent_fold_enrichment",
        "fdr": "divergent_fdr",
    })
)

conserved_vs_divergent = (
    spatially_conserved_pathways[["library", "pathway", "fold_enrichment", "fdr"]]
    .rename(columns={"fold_enrichment": "conserved_fold_enrichment", "fdr": "conserved_fdr"})
    .merge(divergent_enrichment, on=["library", "pathway"], how="outer")
)

conserved_vs_divergent["spatially_conserved"] = (
    conserved_vs_divergent["conserved_fdr"] < CONFIG["enrichment_fdr"]
)
conserved_vs_divergent["spatially_divergent"] = (
    conserved_vs_divergent["divergent_fdr"] < CONFIG["enrichment_fdr"]
)
conserved_vs_divergent["classification"] = np.select(
    [
        conserved_vs_divergent["spatially_conserved"] & ~conserved_vs_divergent["spatially_divergent"],
        ~conserved_vs_divergent["spatially_conserved"] & conserved_vs_divergent["spatially_divergent"],
        conserved_vs_divergent["spatially_conserved"] & conserved_vs_divergent["spatially_divergent"],
    ],
    ["spatially conserved", "spatially divergent", "mixed (each significant)"],
    default="neither significant",
)

_save_table(conserved_vs_divergent, "conserved_vs_divergent_pathways.csv")
print(conserved_vs_divergent["classification"].value_counts().to_string())

# --------------------------------------------------------------------------------------
# Figure: what is conserved, what is remodelled.
# --------------------------------------------------------------------------------------

figure, axes = plt.subplots(1, 2, figsize=(14, 5.8), width_ratios=[1.0, 1.15])
ax_conserved, ax_contrast = axes

top_conserved = spatially_conserved_pathways.nlargest(15, "fold_enrichment")
positions = np.arange(len(top_conserved))
ax_conserved.barh(positions, top_conserved["fold_enrichment"], color=PHENOTYPE_COLORS["conserved zonation"],
                  alpha=0.9)
ax_conserved.set_yticks(positions)
ax_conserved.set_yticklabels(
    [textwrap.fill(f"{row.pathway}", 40) for row in top_conserved.itertuples()], fontsize=7.5
)
ax_conserved.set_xlabel("Fold enrichment for conserved zonation")
ax_conserved.set_title(f"Spatially conserved functions (FDR < {CONFIG['enrichment_fdr']})",
                       loc="left", fontsize=10)
for index, row in enumerate(top_conserved.itertuples()):
    ax_conserved.text(row.fold_enrichment, index,
                      f"  n={row.n_pathway_members_in_universe}, r={row.median_member_shape_corr:.2f}",
                      va="center", fontsize=7)

plot_frame = conserved_vs_divergent.dropna(
    subset=["conserved_fold_enrichment", "divergent_fold_enrichment"]
).copy()

for classification, colour in [
    ("spatially conserved", PHENOTYPE_COLORS["conserved zonation"]),
    ("spatially divergent", PHENOTYPE_COLORS["complex shape rewiring"]),
    ("mixed (each significant)", "#8856A7"),
]:
    block = plot_frame[plot_frame["classification"].eq(classification)]
    if not len(block):
        continue
    ax_contrast.scatter(block["conserved_fold_enrichment"], block["divergent_fold_enrichment"],
                        s=22, alpha=0.75, color=colour, label=f"{classification} (n={len(block):,})")

ax_contrast.axhline(1.0, color="black", lw=0.8, ls="--")
ax_contrast.axvline(1.0, color="black", lw=0.8, ls="--")
ax_contrast.set_xlabel("Fold enrichment for conserved zonation")
ax_contrast.set_ylabel("Fold enrichment for any divergent phenotype")
ax_contrast.set_title("Conserved versus remodelled programmes\n"
                      "upper-left = remodelled, lower-right = conserved", loc="left", fontsize=10)
ax_contrast.legend(fontsize=7.5)

figure.suptitle("Figure: contrast between spatially conserved and spatially divergent pathway families\n"
                "pathway-level over-representation of genes in each phenotype class",
                fontsize=11)
figure.tight_layout()
_save_figure(figure, "conserved_vs_divergent_pathway_summary.png")
plt.show()


# %% [markdown]
# ## 12 - Literature validation (deferred)
#
# **Status: deferred, by decision.** The purpose of this section is external validation, kept strictly
# apart from discovery: score a *published* adaptive/failed-repair PT signature continuously along the
# registered coordinate and ask whether its position is conserved, shifted, expanded or absent. That
# requires a signature whose provenance can be stated exactly.
#
# **What is available today.** No published supplementary signature is stored anywhere in this
# repository (checked: `docs/`, `analysis/`, `pseudospace/`, `configs/`, and the saved results trees).
# What exists is the individual marker genes reported in the human/mouse adaptive-repair literature
# (`VCAM1`, `DCDC2`, `HAVCR1`, `SPP1`, `PROM1`, `VIM`, `S100A6`, `ANXA4`) and the project's own
# `celltyping/epithelial_stress_marker_profiles.csv` panel.
#
# **The decision taken.** Rather than present a handful of marker genes as if it were a signature, the
# markers are registered here as an explicit, auditable list and the scoring is deferred until a
# published signature is available or the panel is reviewed. The cell below records the list and reports
# which of the genes are even present in the tested universe - a check that must pass before any
# signature score can mean anything, and one that is useful on its own.
#
# **Why it is kept separate.** Literature markers must not enter the identification of the phenotypes;
# a gene's class here comes only from its own fitted curves and the section-5 sensitivity panel. When
# this section is completed it will read the atlas, not the reverse.
#

# %%
# Purpose: section 12 - the literature anchor list, recorded and audited, with scoring deferred.

LITERATURE_MARKERS = {
    "adaptive repair / failed repair (human and mouse PT)": [
        "VCAM1", "DCDC2", "HAVCR1", "SPP1", "PROM1", "VIM", "S100A6", "ANXA4",
    ],
}

EPITHELIAL_STRESS_PANEL_PATH = UPSTREAM_DIR / "celltyping" / "epithelial_stress_marker_profiles.csv"

literature_rows = []
for signature_name, markers in LITERATURE_MARKERS.items():
    for gene in markers:
        index = gene_lookup_tested.get(gene.upper())
        literature_rows.append({
            "signature": signature_name,
            "gene_as_published": gene,
            "present_in_tested_universe": index is not None,
            "gene_in_object": "" if index is None else gene_names[index],
            "detection_mouse": np.nan if index is None else float(detection_mouse_array[index]),
            "detection_human": np.nan if index is None else float(detection_human_array[index]),
            "mouse_amplitude": np.nan if index is None else float(gene_metrics["mouse_amplitude"].iloc[index]),
            "human_amplitude": np.nan if index is None else float(gene_metrics["human_amplitude"].iloc[index]),
            "spatial_phenotype": "" if index is None else str(phenotype_by_row[index]),
            "robustness_score": np.nan if index is None else float(gene_metrics["robustness_score"].iloc[index]),
        })

literature_audit = pd.DataFrame(literature_rows)
display(literature_audit.round(4))
_save_table(literature_audit, "literature_marker_audit.csv")

if EPITHELIAL_STRESS_PANEL_PATH.exists():
    print(f"project stress-marker panel available for the deferred signature work: "
          f"{EPITHELIAL_STRESS_PANEL_PATH.relative_to(PROJECT_DIR)}")
else:
    print("project stress-marker panel not found")

_report("literature markers present in the tested universe",
        int(literature_audit["present_in_tested_universe"].sum()), len(literature_audit),
        why="a signature cannot be scored on genes the object does not carry")
print("Section 12 scoring is deferred: no published supplementary signature is stored in this "
      "repository. When one is available, it is scored here along the registered coordinate and "
      "compared with the phenotype assignments above - never used to create them.")


# %% [markdown]
# ## 13 - Tables, the run inventory, and the saved analysis notes
#
# **The atlas** (`tables/gene_spatial_rewiring_atlas.csv`) is the per-gene result: level effect, both
# amplitudes and their ratio, standardised curve correlation on both coordinate choices, best
# displacement with its improvement and residual, positional centroids and peaks, early-to-late
# gradients, half-max widths, within-species reproducibility, detection and abundance, the two exclusion
# flags, the phenotype and its robustness score. Everything a reader needs to disagree with a
# classification without refitting anything is here.
#
# **The shortlist** (`tables/top_genes_by_spatial_phenotype.csv`) ranks each class by robustness first
# and divergence score second, so the genes quoted in the figures are the ones that keep their broad
# class across the sensitivity panel.
#
# **Every table and figure this run produced is listed below with its path.** Nothing is deleted
# silently: the counts of genes entering and leaving each step are printed where the step happens, and
# the exclusion reasons are in the flags.
#
# **The analysis notes** (`spatial_rewiring/analysis_notes.md`) are written from the values computed in
# this run — cohort, support, thresholds actually used, class counts, and the caveats that constrain how
# the numbers may be read. They are the file to quote when this analysis is written up.
#

# %%
# Purpose: the analysis atlas, the per-class shortlist, the run inventory and the saved notes.

atlas_columns = [
    "gene",
    "level_effect_human_minus_mouse", "level_effect_human_minus_mouse_unregistered",
    "mouse_amplitude", "human_amplitude", "amplitude_difference_human_minus_mouse",
    "amplitude_log2_ratio_human_over_mouse", "amplitude_log2_ratio_human_over_mouse_unregistered",
    "shape_corr", "shape_corr_unregistered", "shape_spearman", "pattern_rms_z",
    "registered_shape_corr", "registered_shape_corr_unregistered",
    "best_shift_human_minus_mouse", "best_shift_human_minus_mouse_unregistered",
    "shift_improvement", "shift_improvement_unregistered", "residual_rms_after_shift",
    "mouse_position_centroid", "human_position_centroid", "centroid_shift_human_minus_mouse",
    "centroid_shift_human_minus_mouse_unregistered",
    "mouse_peak_position", "human_peak_position", "peak_shift_human_minus_mouse",
    "mouse_early_to_late", "human_early_to_late", "gradient_change_human_minus_mouse",
    "mouse_halfmax_onset", "mouse_halfmax_width", "human_halfmax_onset", "human_halfmax_width",
    "mouse_reproducibility", "human_reproducibility",
    "detection_mouse", "detection_human", "abundance_mouse", "abundance_human",
    "n_shared_grid_points", "enough_shared_support", "detected_in_both_species",
    "axis_basis_gene", "technical_gene",
    "spatial_phenotype", "broad_phenotype", "spatial_discovery_score",
    "n_variants_tested", "n_variants_same_broad_phenotype", "robustness_score",
]

gene_spatial_rewiring_atlas = gene_metrics[atlas_columns].copy()
_save_table(gene_spatial_rewiring_atlas, "gene_spatial_rewiring_atlas.csv")

TOP_N_PER_CLASS = 25

top_genes_by_spatial_phenotype = (
    gene_metrics[
        discovery_eligible & gene_metrics["spatial_phenotype"].isin(PHENOTYPE_CLASSES)
    ]
    .sort_values(["spatial_phenotype", "robustness_score", "spatial_discovery_score"],
                 ascending=[True, False, False])
    .groupby("spatial_phenotype", group_keys=False)
    .head(TOP_N_PER_CLASS)
    .reset_index(drop=True)
)[[
    "spatial_phenotype", "gene", "robustness_score", "spatial_discovery_score",
    "mouse_amplitude", "human_amplitude", "amplitude_log2_ratio_human_over_mouse",
    "shape_corr", "registered_shape_corr", "best_shift_human_minus_mouse",
    "centroid_shift_human_minus_mouse", "level_effect_human_minus_mouse",
    "mouse_reproducibility", "human_reproducibility",
    "detection_mouse", "detection_human",
]]

_save_table(top_genes_by_spatial_phenotype, "top_genes_by_spatial_phenotype.csv")

# --------------------------------------------------------------------------------------
# Run inventory.
# --------------------------------------------------------------------------------------

table_files = sorted(path.name for path in TABLE_DIR.glob("*.csv"))
figure_files = sorted(path.name for path in FIG_DIR.glob("*.png"))
object_files = sorted(path.name for path in OUT_DIR.glob("*.npz"))

print(f"Tables ({len(table_files)}) in {TABLE_DIR.relative_to(PROJECT_DIR)}:")
for name in table_files:
    print(f"  {name}")
print(f"\nFigures ({len(figure_files)}) in {FIG_DIR.relative_to(PROJECT_DIR)}:")
for name in figure_files:
    print(f"  {name}")
print(f"\nCurve payloads ({len(object_files)}) in {OUT_DIR.relative_to(PROJECT_DIR)}:")
for name in object_files:
    print(f"  {name}")

# --------------------------------------------------------------------------------------
# Analysis notes, written from this run's own values.
# --------------------------------------------------------------------------------------

analysis_notes = f"""# Human versus healthy-mouse PT spatial rewiring (notebook 06)

Written by `analysis/notebooks/06_human_mouse_spatial_rewiring.ipynb` at run time from the values it
computed. Read together with the notebook's own section markdown, which states the same caveats.

## Inputs and cohort

- Input object: `{DPT_OUTPUT_PATH.relative_to(PROJECT_DIR)}` (03's PT-specific DPT; `total_scanpy_dpt`).
- Mouse specimens: {', '.join(mouse_samples)}. Human slices: {', '.join(human_samples)} - **one donor**.
- PT structures compared: {adata_pt.n_obs:,}. Shared pseudospace support: [{lo:.4f}, {hi:.4f}] on a
  {grid.size}-point grid; registration and shifts are reported in unit coordinates over that support.
- Tested orthologs: {len(gene_names):,} (>= {CONFIG['min_detected_fraction_all']:.0%} detection overall
  and measured in both inputs); discovery-eligible: {int(discovery_eligible.sum()):,}.

## Method, in one paragraph

Per-specimen curves are fitted with the shared-smooth model at the pooled per-gene smoothing parameter
and averaged with equal specimen weight; the pooled fit is 03/05's cached nested level/shape fit
(`section4_gene_fit`). Cross-species position is registered by a monotone piecewise-linear mapping
anchored on S1/S2/S3 landmark program centroids ({primary_registration['n_landmarks']} landmarks used).
Metrics are level, amplitude, standardised shape correlation, positional centroid/peak/early-to-late
gradient/half-max window, and a bounded displacement search (|shift| <= {CONFIG['max_registration_shift']}
of PT). Every metric is computed on the registered axis and on the unregistered axis.

## Caveats that constrain every number here

- The two human slices are one donor. Nothing here is a population-level species test; no p-value on a
  cross-species contrast is reported or implied. Slice agreement is robustness, not replication.
- Species is inseparable from sampling and from Harmony batch correction (`sample` is the batch key).
- Both human slices are cortex; the `MED1` label does not make the tissue medullary.
- Alignment depends on {primary_registration['n_landmarks']} landmark programs; a gene is only called
  shifted when the shift survives both coordinate choices with the same sign and improves agreement on
  both. Registration error inside a segment is not corrected.
- Phenotype classes are threshold conjunctions. The thresholds are in the notebook's `CONFIG` and the
  distributions they act on are printed in section 4; the component metrics are in the atlas.

## Thresholds actually used

{chr(10).join(f'- `{key}` = {value}' for key, value in CONFIG.items() if not hasattr(value, '__len__'))}

## Phenotype counts (all tested orthologs)

{chr(10).join(f'- {label}: {int((gene_metrics["spatial_phenotype"] == label).sum()):,}' for label in PHENOTYPE_CLASSES + ["excluded"])}

## Robustness

- Genes whose broad class survived every variant: {int((gene_metrics['robustness_score'] == 1).sum()):,}.
- Discovery-eligible genes robust in >= 80% of variants:
  {int((gene_metrics['robustness_score'].ge(0.8) & discovery_genes).sum()):,}.

## Files

- Tables: `tables/` ({len(table_files)} files), figures: `figures/` ({len(figure_files)} files),
  curve payloads: this directory ({len(object_files)} `.npz`).
- Section 12 (literature validation) is deferred: no published supplementary signature is stored in
  this repository; the marker audit is in `tables/literature_marker_audit.csv`.
"""

(OUT_DIR / "analysis_notes.md").write_text(analysis_notes)
print(f"\nsaved {(OUT_DIR / 'analysis_notes.md').relative_to(PROJECT_DIR)}")


# %% [markdown]
# ## 14 - Closing summary
#
# The cell below prints only what this run measured: how many orthologs entered, how many show
# reproducible zonation in each species, how the classes came out, how many of those survive the
# sensitivity panel, which functional programmes are most conserved and most divergent, and which
# pathways carry a spatial signal that a whole-PT comparison would have missed.
#
# Nothing here is a hard-coded conclusion. If the statement the analysis is meant to support -
# "human and mouse proximal tubules share a continuous transcriptional architecture, while specific
# genes and functional programmes diverge through distinct modes of spatial remodelling, including
# altered zonation strength, positional displacement and trajectory rewiring that are obscured by
# whole-PT expression comparisons" - is not supported by the numbers printed here, then it is not
# supported, and the numbers win.
#

# %%
# Purpose: the auto-generated closing summary - measured values only, no hard-coded conclusions.

def _count(mask):
    return int(np.asarray(mask).sum())


summary_lines = []
summary_lines.append("=== Human vs healthy-mouse PT spatial rewiring: measured summary ===")
summary_lines.append("")
summary_lines.append(f"Tested orthologs: {len(gene_names):,} "
                     f"(of {adata_pt.n_vars:,} measured features on {adata_pt.n_obs:,} PT structures)")
summary_lines.append(f"Discovery-eligible (no axis/technical flag, detected in both species, "
                     f"enough shared support): {_count(discovery_eligible):,}")
summary_lines.append(f"Landmarks used for cross-species registration: "
                     f"{primary_registration['n_landmarks']} "
                     f"(max mapping deviation from identity "
                     f"{float(np.max(np.abs(registration_shift))):.3f} of PT)")
summary_lines.append("")
summary_lines.append("Reproducible zonation, by species "
                     f"(amplitude >= {CONFIG['amplitude_patterned']}, within-species r >= "
                     f"{CONFIG['within_species_corr']}):")
summary_lines.append(f"  mouse: {_count((gene_metrics['mouse_amplitude'] >= CONFIG['amplitude_patterned']) & (gene_metrics['mouse_reproducibility'] >= CONFIG['within_species_corr'])):,}")
summary_lines.append(f"  human: {_count((gene_metrics['human_amplitude'] >= CONFIG['amplitude_patterned']) & (gene_metrics['human_reproducibility'] >= CONFIG['within_species_corr'])):,}")
summary_lines.append("")
summary_lines.append("Spatial phenotype counts (all tested orthologs / discovery-eligible):")
for label in PHENOTYPE_CLASSES:
    in_class = gene_metrics["spatial_phenotype"].eq(label)
    summary_lines.append(
        f"  {label}: {_count(in_class):,} / {_count(in_class & discovery_eligible):,}"
    )
summary_lines.append(f"  excluded (axis, technical, not in both species): {_count(gene_metrics['spatial_phenotype'].eq('excluded')):,}")
summary_lines.append("")
summary_lines.append("Sensitivity panel:")
summary_lines.append(f"  variants tested per gene: {int(gene_metrics['n_variants_tested'].max()):,}")
summary_lines.append(f"  genes keeping the same broad phenotype in every variant: "
                     f"{_count(gene_metrics['robustness_score'] == 1):,}")
summary_lines.append(f"  discovery-eligible genes robust in >= 80% of variants: "
                     f"{_count(gene_metrics['robustness_score'].ge(0.8) & discovery_eligible):,}")
summary_lines.append("")
summary_lines.append(f"Significant pathway x phenotype enrichments (FDR < {CONFIG['enrichment_fdr']}): "
                     f"{len(best_per_pathway):,} pathways of {len(pathway_tested):,} tested")
for label in PHENOTYPE_CLASSES:
    block = best_per_pathway[best_per_pathway["spatial_phenotype"].eq(label)]
    if len(block):
        summary_lines.append(f"  {label}: {len(block):,} pathways")
summary_lines.append("")
summary_lines.append("Most conserved programmes (enriched for conserved zonation):")
for row in spatially_conserved_pathways.head(8).itertuples():
    summary_lines.append(f"  {row.pathway} [{row.library}] fold={row.fold_enrichment:.2f} "
                         f"FDR={row.fdr:.3g} members={row.n_pathway_members_in_universe}")
summary_lines.append("")
summary_lines.append("Most spatially divergent programmes (strongest divergent-phenotype enrichment):")
for row in divergent_enrichment.sort_values("divergent_fdr").head(8).itertuples():
    summary_lines.append(f"  {row.pathway} [{row.library}] -> {row.dominant_divergent_phenotype} "
                         f"fold={row.divergent_fold_enrichment:.2f} FDR={row.divergent_fdr:.3g}")
summary_lines.append("")
summary_lines.append("Continuous-only signal (significant spatial enrichment, conventional effect "
                     "below the pathway median):")
summary_lines.append(f"  {len(continuous_only_pathways):,} pathways")
for row in continuous_only_pathways.head(8).itertuples():
    summary_lines.append(f"  {row.pathway} [{row.library}] -> {row.dominant_spatial_phenotype} "
                         f"FDR={row.enrichment_fdr:.3g} conventional percentile "
                         f"{row.conventional_percentile:.2f}")
summary_lines.append("")
summary_lines.append("Caveats: two mouse specimens and two human slices from ONE donor; no "
                     "population-level species inference is reported; species is confounded with "
                     "sampling and batch correction; alignment rests on "
                     f"{primary_registration['n_landmarks']} landmark programmes.")

summary_text = "\n".join(summary_lines)
print(summary_text)
(OUT_DIR / "analysis_summary.txt").write_text(summary_text + "\n")
print(f"\nsaved {(OUT_DIR / 'analysis_summary.txt').relative_to(PROJECT_DIR)}")


# %%
# Purpose: deterministic synthetic self-checks for the inline metric helpers.

from pseudospace.stats_gam import bh_adjust as _bh_reference


def _metric_selfcheck():
    """Check the notebook's own metrics against constructed curves with known answers."""
    x = np.linspace(0.0, 1.0, 101)
    bump = np.exp(-((x - 0.3) ** 2) / (2 * 0.05 ** 2)).reshape(1, -1)
    flat = np.full((1, x.size), 2.0)
    ramp = x.reshape(1, -1)
    inverted = (-bump + 3.0)

    checks = []

    checks.append(("identical curves correlate at 1",
                   np.isclose(_row_pearson(bump, bump)[0], 1.0)))
    checks.append(("mirrored curves correlate at -1",
                   np.isclose(_row_pearson(bump, -bump)[0], -1.0)))
    checks.append(("ties-safe Spearman on a complete row equals Pearson for a monotone pair",
                   abs(_row_spearman(ramp, ramp)[0] - 1.0) < 1e-9))

    z_bump = _row_zscore(bump)[0]
    checks.append(("z-score has mean 0", abs(np.nanmean(z_bump)) < 1e-12))
    checks.append(("z-score has SD 1", abs(np.nanstd(z_bump) - 1.0) < 1e-9))

    checks.append(("flat curve has zero amplitude",
                   np.isclose(_row_amplitude(flat)[0], 0.0)))
    checks.append(("amplitude of a bump is its peak-to-trough",
                   abs(_row_amplitude(bump)[0] - (bump.max() - bump.min())) < 1e-12))

    checks.append(("weighted centroid recovers the bump centre",
                   abs(_row_weighted_centroid(bump, x)[0] - 0.3) < 0.02))
    checks.append(("peak position recovers the bump centre",
                   abs(_row_peak_position(bump, x)[0] - 0.3) <= 0.02))
    checks.append(("early-to-late contrast is positive for a rising ramp",
                   _row_early_to_late(ramp, x)[0] > 0))
    checks.append(("early-to-late contrast is negative for a falling ramp",
                   _row_early_to_late(-ramp, x)[0] < 0))

    onset, offset, width = _row_halfmax_window(bump, x)
    checks.append(("half-maximum window brackets the peak",
                   onset[0] < 0.3 < offset[0] and width[0] > 0))

    partial = bump.copy()
    partial[0, :40] = np.nan
    checks.append(("unsupported points stay NaN, never zero",
                   bool(np.isnan(_row_zscore(partial)[0, :40]).all())))
    checks.append(("unsupported points are excluded from the correlation",
                   np.isclose(_row_pearson(partial, partial)[0], 1.0)))

    resampled = _interp_rows(bump, x, x + 0.5)[0]
    checks.append(("resampling outside the source range is NaN (no extrapolation)",
                   bool(np.isnan(resampled[-1])) and bool(np.isfinite(resampled[0]))))

    # A curve displaced in human must be recovered with the matching sign:
    # H(s) = M(s - d) puts the human program LATER by d (positive shift), H(s) = M(s + d) earlier.
    later = np.interp(x, x + 0.10, bump[0], left=np.nan, right=np.nan).reshape(1, -1)
    earlier = np.interp(x, x - 0.10, bump[0], left=np.nan, right=np.nan).reshape(1, -1)

    later_shift, later_corr, later_residual = _row_shift_search(
        bump, later, x, 0.20, 0.01, min_points=10
    )
    earlier_shift, _, _ = _row_shift_search(bump, earlier, x, 0.20, 0.01, min_points=10)

    checks.append(("a later human program is recovered as a positive displacement",
                   abs(later_shift[0] - 0.10) <= 0.02))
    checks.append(("an earlier human program is recovered as a negative displacement",
                   abs(earlier_shift[0] + 0.10) <= 0.02))
    checks.append(("recovering the displacement improves the correlation",
                   later_corr[0] > _row_pearson(bump, later)[0]))
    checks.append(("residual after the best displacement is ~0 for a pure move",
                   later_residual[0] < 0.05))

    _, inverted_corr, _ = _row_shift_search(bump, inverted, x, 0.20, 0.01, min_points=10)
    checks.append(("an inverted curve is not rescued by displacement",
                   not np.isfinite(inverted_corr[0]) or inverted_corr[0] < 0.8))

    checks.append(("BH adjustment matches the package implementation",
                   np.allclose(bh_adjust(np.array([0.01, 0.02, 0.5])),
                               _bh_reference(np.array([0.01, 0.02, 0.5])))))

    return pd.DataFrame([{"check": name, "passed": bool(passed)} for name, passed in checks])


metric_selfcheck = _metric_selfcheck()
display(metric_selfcheck)

if not metric_selfcheck["passed"].all():
    failed = metric_selfcheck.loc[~metric_selfcheck["passed"], "check"].tolist()
    raise AssertionError(f"synthetic metric self-check failed: {failed}")

print(f"All {len(metric_selfcheck)} synthetic metric self-checks passed.")

