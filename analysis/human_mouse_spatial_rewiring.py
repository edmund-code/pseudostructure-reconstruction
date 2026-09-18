# %% [markdown]
# # Human versus healthy-mouse proximal tubule analysis
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
# specimens' own fitted curves. Absolute abundance and spatial pattern are analysed separately by
# construction: `level` is the mean vertical offset of the two curves, `amplitude` is peak-to-trough of
# each curve, `shape` compares curves after each is standardised within its species, and `displacement`
# is a constrained horizontal shift that has to earn its interpretation against a fixed maximum.
#
# ## Read this before the numbers
#
# * **The two human slices (`HUK1_COR1`, `HUK1_MED1`) come from one donor.** They are not independent
#   human biological replicates and no column in this notebook is a population-level human-versus-mouse
#   test. Every contrast is a descriptive effect size on a specific balanced-curve pair.
# * Agreement between the two human slices, and between the two mouse specimens, is used as a
#   **robustness** measure — "is this curve reproducible within a species" — never as replication-based
#   inference. Section 3.2 leaves one slice out as *slice sensitivity*, not as a replication test.
# * **Species is inseparable from sampling and from Harmony batch correction** in this cohort
#   (`sample` is the batch key). A cross-species difference in a fitted curve is a difference between
#   these two mice and this one donor, processed this way.
# * Both human slices are anatomically cortex; the `MED1` label does not make the tissue medullary.
# * Genes that were used to build or orient the PT coordinate are flagged (`axis_basis_gene`) and are
#   excluded from every unbiased discovery list, but they are retained for the QC and validation plots.
# * Direction conventions, used everywhere below: **positive level effect = higher in human**;
#   **positive shift = the human positional program occurs later along PT than the mouse program**;
#   positive early-to-late gradient = expression increases from early to late PT.
#
# ## The coordinate this notebook analyses
#
# The coordinate is **03's shared PT DPT**, used exactly as 03 oriented it: no cross-species
# re-registration is applied, so a displacement is read as evidence that the two species' spatial
# programmes differ, not as an anatomical measurement between aligned axes. Two consequences are built
# into the analysis rather than left to the reader. First, the displacement criteria are strict — an
# interior optimum, a positional centroid that agrees, and an improved shape agreement — and a
# displacement is reported as a **supporting measurement** of the genes carrying it rather than as a
# phenotype class, so no gene's class rests on a coordinate difference alone. Second, a gene that is
# flat or noisy in one species is never described as spatially rewired, whatever its level difference.
#
# ## Layout
#
# | Section | Question it answers |
# | --- | --- |
# | 0 | Setup, cohort, gene universe, direction conventions |
# | 1 | Balanced human and mouse PT curves, and the curve-layer QC |
# | 2 | How each gene differs along PT: the five quantities, and the conserved-architecture null |
# | 3 | Spatial phenotypes, their robustness across the analysis choices, and the result figures |
# | 4 | Which functional programs sit in which spatial phenotype |
# | 5 | What the continuous analysis adds beyond the whole-PT comparison |
# | 6 | Tables, the run inventory and the closing summary |
#
# Outputs are written to `results/human_vs_healthy_mouse/spatial_rewiring/` (`figures/`, `tables/`,
# the fitted-curve `.npz` payloads and `analysis_notes.md`). The **stage cache is deliberately
# 03/05's** `results/human_vs_healthy_mouse/stage_cache/` rather than a private subdirectory:
# `cached_run_level_shape` keys on the fitter's source plus the fit parameters, so pointing at the
# shared cache turns the pooled gene fit into a reuse instead of a silent recomputation. This
# notebook's own per-specimen fits are cached there as `spatial_rewiring_specimen_curves`.
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
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
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
# section 3 reports the empirical distributions these values are read against.
# ======================================================================================

CONFIG = {
    # --- curve fitting (identical to 03/05 so the shared cache applies) ------------------
    "grid_points": 101,
    "common_support_pct": (1, 99),
    "n_internal_knots": 9,
    "lambda_grid": np.logspace(-3, 3, 13),

    # --- coordinate -----------------------------------------------------------------------
    # One coordinate throughout: 03's shared PT DPT (the PT-specific `total_scanpy_dpt` column of
    # `cross_species_pt_dpt.h5ad`). Every positional statement below is a statement about that axis
    # in this cohort; no cross-species re-registration is applied to it.

    # --- gene universe -------------------------------------------------------------------
    "min_detected_fraction_all": 0.02,     # as in 03/05: >= 2% of all PT structures
    "min_detection_each_species": 0.05,    # a pathway member must be seen in both species
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

    # --- displacement (a supporting measurement, not a class) ------------------------------
    # The search is deliberately wider than any displacement this notebook will interpret: an
    # optimum that lands on the search boundary has not been estimated at all, it only says that
    # agreement keeps improving outwards. Interpreted displacements must sit inside the search
    # with room to spare, and must agree with the independent positional centroid.
    "shift_search_limit": 0.35,            # how far the diagnostic search looks
    "shift_interior_margin": 0.05,         # an optimum must be at least this far from the boundary
    "shift_step": 0.01,
    "minimum_shift": 0.08,                 # a shift must exceed this to be interpreted
    "centroid_shift_tolerance": 0.10,      # the centroid shift must corroborate the curve shift
    "post_shift_corr": 0.80,               # ... and must reach this shape agreement
    "shift_improvement": 0.10,             # ... and improve on the unshifted correlation by this

    # --- gradient inversion ---------------------------------------------------------------
    "gradient_min_abs": 0.10,              # both species need an appreciable early-to-late slope
    "inversion_max_corr": 0.30,

    # --- conventional versus continuous ------------------------------------------------
    # Strong spatial evidence is a categorical state, not a median split of the score: the gated
    # score is exactly 0 for every gene with nothing interpretable, so a median split would call every
    # gene "high spatial" and empty two quadrants. The score then ranks *within* the strong group.
    "spatial_strong_robustness": 0.8,

    # --- pathway ranking tests ------------------------------------------------------------
    # The conventional comparator is a whole-PT preranked GSEA on the level effect - how a whole-PT
    # paper would analyse these pathways - and the conservation view is the same test on a continuous
    # conservation ranking, because the binary conserved class is too small for over-representation.
    "gsea_permutations": 1000,
    "gsea_seed": 42,
    "gsea_min_size": 10,
    "gsea_max_size": 500,

    # --- conserved-architecture null ------------------------------------------------------
    "null_permutations": 5000,          # this test is cheap; 200 was too few for a stable p

    "null_seed": 0,
    "null_min_genes": 50,                  # below this the conservation question is not answerable

    # --- pathway enrichment ---------------------------------------------------------------
    "pathway_min_members": 10,
    "pathway_max_members": 500,
    "enrichment_fdr": 0.05,
    "enrichment_min_members_flagged": 3,

    # A pathway-level label needs member consistency, not just over-representation: enrichment says
    # the phenotype is commoner among the members than in the background; consistency asks whether the
    # members behave that way overall. Without it, "5 of 82 members are conserved" becomes "this
    # pathway is a conserved programme". The behavioural test is phenotype-specific, because a
    # negative median shape correlation is the *expected* signature of inversion or complex rewiring.
    "pathway_consistency_min_fraction": 0.15,             # share of members holding the phenotype
    "pathway_consistency_min_median_shape_corr": 0.0,     # only applied to conserved zonation
    # Only near-duplicate sets are collapsed: a lenient threshold reduces 1,511 pathways to a handful
    # of groups and hides the associations the section exists to report.
    "pathway_redundancy_overlap": 0.85,

    # --- figures --------------------------------------------------------------------------
    "top_labels_per_panel": 12,
    # The complete three-gene atlas is retained as a supplementary figure.  The main figure uses
    # one purpose-selected exemplar per phenotype so it can be read at a single-column page scale.
    "representative_genes_per_class": 3,
    "main_representative_genes_per_class": 1,
    "heatmap_genes_per_class": 6,
}

# How far from the search boundary an optimum must sit to count as an interior estimate.
SHIFT_INTERIOR_LIMIT = CONFIG["shift_search_limit"] - CONFIG["shift_interior_margin"]

PHENOTYPE_COLORS = {
    "conserved zonation": "#1B7837",
    "weaker zonation in human": "#8C6BB1",
    "stronger zonation in human": "#762A83",
    "mouse-zonated / human-flat": "#0072B2",
    "human-zonated / mouse-flat": "#D55E00",
    "gradient inversion": "#C2185B",
    "complex shape rewiring": "#4D4D4D",
    "weak / uncertain zonation": "#BDBDBD",
    "excluded": "#E8E8E8",
}

SPECIES_COLORS = {
    "mouse": "#0072B2",
    "human": "#D55E00",
}

# The displacement metrics are a supporting measurement, not a class: this colour marks them wherever
# they are drawn, so a displacement can never be mistaken for a phenotype.
DISPLACEMENT_COLOR = "#E69F00"

# Direction conventions, restated once so the sign of every column is unambiguous.
DIRECTION_CONVENTIONS = {
    "level_effect_human_minus_mouse": "positive = higher in human",
    "amplitude_log2_ratio_human_over_mouse": "positive = stronger zonation in human",
    "shift_human_minus_mouse": "positive = human program occurs later along PT",
    "gradient_early_to_late": "positive = expression rises from early to late PT",
    "pattern_rms_z": "standardised RMS difference; 0 = identical spatial pattern",
}

# The canonical PT marker panels 03 uses to build and orient the PT coordinate (its
# NEPHRON_AXIS_MARKERS and PT fine panels, plus the S3 additions used for the mouse-only axis). They
# are used here for two things only: the `axis_basis_gene` flag below, so that the basis of the
# coordinate cannot also be presented as a finding about it, and the curve-layer QC panel (figure 1A).
AXIS_BASIS_PANELS = {
    "S1": ["Slc5a2", "Slc5a12", "Gatm", "Lrp2", "Cubn", "Slc34a1"],
    "S2": ["Slc22a6", "Slc13a3", "Cyp2e1"],
    "S3": ["Slc22a7", "Cyp7b1", "Slc7a13", "Slc6a18", "Acsm3"],
}

# Every gene that contributed to building or orienting the PT coordinate.
AXIS_BASIS_GENES = tuple(sorted({gene for panel in AXIS_BASIS_PANELS.values() for gene in panel}))

# Obvious technical features: mitochondrial and cytoplasmic ribosomal protein genes.
TECHNICAL_GENE_PATTERNS = (
    r"^MT-",            # human mitochondrial
    r"^MT(ND|CO|ATP|CYB|RNR)\d",   # mouse mitochondrial (Mt-Nd1 ...)
    r"^RP[SL]\d+[A-Z]?$",          # cytoplasmic ribosomal proteins
    r"^MRP[SL]\d+$",               # mitochondrial ribosomal proteins
)


def _display_path(path):
    """The path relative to the repository root when it lies inside it, otherwise the full path.

    The PSEUDOSPACE_* roots may legitimately point outside the repository - that is how the project
    keeps private data out of Git - so relative_to() cannot be assumed to succeed.
    """
    path = Path(path)
    try:
        return path.relative_to(PROJECT_DIR)
    except ValueError:
        return path


def _save_figure(fig, name, dpi=600):
    """Write review and manuscript-ready raster/vector versions of a publication figure.

    The PNG keeps notebook review convenient; the uncompressed TIFF, PDF and SVG support journal
    submission and manuscript assembly.  The filename passed by each plotting cell stays the
    canonical figure stem, so the four exports cannot silently drift apart.
    """
    path = FIG_DIR / name
    # Render-time geometry audit is opt-in so the notebook does not depend on development-only
    # tooling.  Set PSEUDOSPACE_FIGURE_QA_DIR to a directory containing audit_panel_alignment.py
    # when preparing a manuscript bundle; a failing alignment then blocks export.
    figure_qa_dir = os.environ.get("PSEUDOSPACE_FIGURE_QA_DIR")
    if figure_qa_dir:
        if figure_qa_dir not in sys.path:
            sys.path.insert(0, figure_qa_dir)
        from audit_panel_alignment import require_matplotlib_panel_alignment
        require_matplotlib_panel_alignment(
            fig,
            json_out=path.with_suffix(".alignment.json"),
            overlay_svg=path.with_suffix(".alignment.svg"),
            tolerance_pt=1.5,
            gutter_tolerance_pt=1.5,
            require_panel_labels=False,
            strict=True,
        )
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    tiff_path = path.with_suffix(".tiff")
    pdf_path = path.with_suffix(".pdf")
    svg_path = path.with_suffix(".svg")
    fig.savefig(tiff_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf_path, bbox_inches="tight", facecolor="white")
    fig.savefig(svg_path, bbox_inches="tight", facecolor="white")
    print(f"saved {_display_path(path)} (+ TIFF/PDF/SVG)")
    return path


def _panel_label(axis, label):
    """Place a stable, manuscript-style panel label outside the data rectangle."""
    axis.text(-0.16, 1.07, label, transform=axis.transAxes, fontweight="bold",
              fontsize=11, va="bottom", ha="left", clip_on=False)


def _style_pt_axis(axis, *, xlabel=True, ylabel=None):
    """Use the same positional vocabulary and lightweight grid across curve panels."""
    if xlabel:
        axis.set_xlabel("PT position (shared DPT): early → late")
    if ylabel:
        axis.set_ylabel(ylabel)
    axis.set_xlim(0, 1)
    axis.grid(axis="y", color="0.90", lw=0.6, zorder=0)


def _save_table(frame, name, **kwargs):
    """Write a table into `tables/` and report the path and the shape that was written."""
    path = TABLE_DIR / name
    frame.to_csv(path, index=False, **kwargs)
    print(f"saved {_display_path(path)}  ({len(frame):,} rows x {frame.shape[1]} cols)")
    return path


def _numeric(frame, columns):
    """Cast named columns to float so an empty (or all-NaN) table can still be compared and ranked.

    A table built from an empty result set has object-dtype columns, and pandas refuses nlargest /
    comparisons on object dtype. Casting keeps the empty case working instead of crashing the cell.
    """
    frame = frame.copy()
    for column in columns:
        if column in frame.columns:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
    return frame


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
    "figure.dpi": 120,
    "savefig.dpi": 600,
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "axes.titlepad": 6,
    "axes.labelpad": 3,
    "axes.linewidth": 0.8,
    "legend.frameon": False,
    "legend.fontsize": 7,
    "legend.handlelength": 1.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

print("Project:      ", PROJECT_DIR)
print("Input:        ", DPT_OUTPUT_PATH)
print("Outputs:      ", OUT_DIR)
print("Stage cache:  ", STAGE_CACHE_DIR, f"({'on' if STAGE_CACHE_ENABLED else 'off'})")
print("Coordinate:   ", "03's shared PT DPT (`total_scanpy_dpt`); no re-registration applied")
print("Direction conventions:")
for _column, _meaning in DIRECTION_CONVENTIONS.items():
    print(f"  {_column}: {_meaning}")
print("Axis-basis genes flagged:", len(AXIS_BASIS_GENES))


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
# * `axis_basis_gene` - the gene helped build or orient the PT coordinate (`AXIS_BASIS_GENES`). Keeping
#   them would make the coordinate's own basis look like the strongest finding in it. They stay in the
#   tables and flag every unbiased discovery list.
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
        why=f"of {len(AXIS_BASIS_GENES)} axis-basis genes")
_report("technical genes in the universe", int(technical_gene.sum()), len(gene_names),
        why="mitochondrial / ribosomal protein")
_report("detected in both species", int((detected_in_both_species & ~axis_basis_gene & ~technical_gene).sum()),
        len(gene_names), why=f"detection >= {CONFIG['min_detection_each_species']:.0%} per species, "
                             "axis and technical genes excluded")

print(f"  axis-basis genes present in the universe: "
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
# * **Nothing extrapolates.** The displacement search samples the human curve only inside its own support;
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


def _nanmin_safe(values, axis=1, keepdims=False):
    """NumPy warns on an all-NaN slice; here that slice is a legitimately unsupported row."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanmin(values, axis=axis, keepdims=keepdims)


def _nanmax_safe(values, axis=1, keepdims=False):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanmax(values, axis=axis, keepdims=keepdims)


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
    return _nanmax_safe(values) - _nanmin_safe(values)


def _row_weighted_centroid(curve, x, mask=None):
    """Positional centroid weighted by the row above its own minimum (stable, uses every point)."""
    curve = np.asarray(curve, dtype=float)
    if mask is None:
        mask = np.isfinite(curve)

    values = np.where(mask, curve, np.nan)
    weights = np.clip(values - _nanmin_safe(values, keepdims=True), 0.0, None)

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
    low = _nanmin_safe(values, keepdims=True)
    high = _nanmax_safe(values, keepdims=True)
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


def _row_shift_search(mouse_z, human_z, x_unit, search_limit, step, min_points=8):
    """Constrained horizontal displacement of each human z-curve onto its mouse z-curve.

    Positive shift = the human positional program lies **later** along PT than the mouse program, so
    the human curve is sampled at ``x + shift`` and compared with the mouse curve at ``x``. The
    search is capped at ``search_limit`` because beyond that "shifted" stops being a description of
    an aligned gradient. Returns (best shift, best correlation, standardised RMS after the best
    shift, whether the optimum sits on the search boundary); the delta=0 correlation is the curve
    correlation the caller computes separately. A boundary optimum is not an estimate: it only
    says agreement kept improving outwards, which is why the caller reports it separately.
    """
    mouse_z = np.asarray(mouse_z, dtype=float)
    human_z = np.asarray(human_z, dtype=float)
    x_unit = np.asarray(x_unit, dtype=float)

    shifts = np.arange(-search_limit, search_limit + step / 2.0, step)
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

    at_boundary = np.isfinite(best_shift) & np.isclose(np.abs(best_shift), search_limit)
    return best_shift, best_corr, residual_rms, at_boundary



# %% [markdown]
# ## 1 - Balanced human and mouse PT curves
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
# **The QC panel.** The section closes with figure 1A: the canonical S1 (early), S2 (mid) and S3 (late)
# PT programmes drawn on the shared coordinate. It is the one check the whole notebook depends on - if
# the canonical programmes did not run early to late in both species, no curve metric below would mean
# anything - so the peak positions and early-to-late slopes are printed with the figure.
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

print(f"saved {_display_path(OUT_DIR / 'balanced_gene_curves.npz')}")
print(f"saved {_display_path(OUT_DIR / 'specimen_gene_curves.npz')}")

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
    "n_shared_grid_points_balanced": shared_support.sum(axis=1),
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

    # atlas column -> this notebook's balanced-curve column. The names differ; the comparison is
    # on the same quantity, so the two runs can be read against each other directly.
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


# %%
# Purpose: figure 1A - curve-layer QC: the canonical PT programmes on the shared coordinate.

# What this panel is for. Every metric below is a property of a *fitted curve*, so the curves have to
# be aligned with the biology they claim to represent before any of them is read. The canonical early
# (S1), mid (S2) and late (S3) PT programmes are therefore plotted on the shared PT DPT first. Each
# species' programme curve is the equal-weight mean of its members' within-species z-scored balanced
# curves, so the panel shows position rather than abundance. The expectation is the obvious one - S1
# peaks early, S3 late, in *both* species - and the peak positions are printed beside the figure so the
# claim is auditable rather than asserted.

AXIS_BASIS_PANEL_ORDER = ["S1", "S2", "S3"]


def _programme_curve(curves, members):
    """The programme's curve: mean of its members' within-species z-scored balanced curves."""
    rows = [gene_lookup_tested[gene.upper()] for gene in members if gene.upper() in gene_lookup_tested]
    if not rows:
        return np.full(curves.shape[1], np.nan)
    return _nanmean_safe(_row_zscore(curves[rows]), axis=0)


programme_members = {
    panel: [gene for gene in AXIS_BASIS_PANELS[panel] if gene.upper() in gene_lookup_tested]
    for panel in AXIS_BASIS_PANEL_ORDER
}
programme_curves = {
    (panel, species_label): _programme_curve(curves, programme_members[panel])
    for panel in AXIS_BASIS_PANEL_ORDER
    for species_label, curves in (("mouse", balanced_mouse), ("human", balanced_human))
}

figure, axes = plt.subplots(1, 3, figsize=(7.25, 2.6), sharey=True, layout="constrained")

for axis, panel in zip(axes, AXIS_BASIS_PANEL_ORDER):
    for gene in programme_members[panel]:
        index = gene_lookup_tested[gene.upper()]
        # Retain the member curves, but make the aggregate programme the visual evidence.
        axis.plot(grid_unit, _row_zscore(balanced_mouse[[index]])[0], color=SPECIES_COLORS["mouse"],
                  lw=0.65, alpha=0.20, zorder=1)
        axis.plot(grid_unit, _row_zscore(balanced_human[[index]])[0], color=SPECIES_COLORS["human"],
                  lw=0.65, alpha=0.20, ls="--", zorder=1)

    for species_label, line_style in (("mouse", "-"), ("human", "--")):
        axis.plot(grid_unit, programme_curves[(panel, species_label)],
                  color=SPECIES_COLORS[species_label], lw=2.1, ls=line_style, zorder=3)
    axis.set_title(f"{panel} programme", loc="left", fontweight="bold")
    axis.text(0.02, 0.05, f"{len(programme_members[panel])} member genes", transform=axis.transAxes,
              fontsize=6.5, color="0.35")
    _style_pt_axis(axis, ylabel="Mean member z-score" if axis is axes[0] else None)

_panel_label(axes[0], "a")
axes[0].legend(
    handles=[
        Line2D([], [], color=SPECIES_COLORS["mouse"], lw=2.1, label="mouse programme"),
        Line2D([], [], color=SPECIES_COLORS["human"], lw=2.1, ls="--", label="human programme"),
        Line2D([], [], color="0.45", lw=0.65, alpha=0.35, label="individual member gene"),
    ], loc="upper right", fontsize=6.5,
)
figure.suptitle(
    "Curve-layer QC: the canonical PT programmes sit early to late on the shared coordinate, in both species\n"
    "Faint traces are member genes; each programme curve is their mean within-species z-scored curve",
    fontsize=9.2,
)
_save_figure(figure, "fig1a_pt_marker_programmes.png")
plt.show()

programme_peak_table = pd.DataFrame([
    {"panel": panel, "species": species_label, "n_members": len(programme_members[panel]),
     "peak_position": float(grid_unit[int(np.nanargmax(programme_curves[(panel, species_label)]))]),
     "early_to_late_slope":
         float(_row_early_to_late(programme_curves[(panel, species_label)][None, :], grid_unit)[0])}
    for panel in AXIS_BASIS_PANEL_ORDER
    for species_label in ("mouse", "human")
])
display(programme_peak_table.round(3))
_save_table(programme_peak_table, "pt_marker_programme_peaks.csv")

for species_label in ("mouse", "human"):
    order = (programme_peak_table.loc[programme_peak_table["species"].eq(species_label)]
             .sort_values("peak_position")["panel"].tolist())
    print(f"  {species_label} programmes, ordered early to late by peak position: {order}")



# %% [markdown]
# ## 2 - How each gene differs along PT
#
# ### 2.1 - The five quantities
#
# **Each answers a different question,** Each answers a different question, in the order the decomposition is written
# (level + amplitude + position + shape):
#
# | | metric | question | invariant to |
# | --- | --- | --- | --- |
# | A | `level_effect_human_minus_mouse` | is the gene's overall abundance different? | nothing (this is the conventional DE-like component) |
# | B | `mouse_amplitude`, `human_amplitude`, `amplitude_log2_ratio_human_over_mouse` | how strong is each species' zonation? | level |
# | C | `shape_corr` (standardised), `pattern_rms_z`, `shape_spearman` | do the two curves have the same spatial pattern? | level **and** amplitude |
# | D | `mouse_position_centroid`, `human_position_centroid`, peak, half-max window | where does expression sit along PT in each species? | level, amplitude |
# | E | `mouse_early_to_late`, `human_early_to_late` | do the two species run in opposite directions? | level, amplitude |
#
# **The one supporting measurement.** The displacement of the human programme relative to the mouse one
# (`best_shift_human_minus_mouse` and its family) is computed for every gene but is deliberately *not*
# one of the five: it answers "would a different position explain the difference?", which is a question
# about A-D rather than a separate way of differing. `best_shift` is the horizontal displacement that
# best aligns the human curve to the mouse curve after both are standardised; positive = the human
# program sits later along PT. The search looks out to +/-0.35 of PT, deliberately **wider** than
# anything this notebook will interpret: an optimum pinned to the search boundary has not been estimated
# at all, it only says that agreement keeps improving outwards - those genes are flagged
# (`shift_at_search_boundary`). An *interpreted* displacement must be (a) interior, with |shift| at or
# below the 0.30 interior limit, (b) at least `minimum_shift` with a real `shift_improvement` and
# `post_shift_shape_corr` >= `post_shift_corr`, and (c) corroborated by the independent positional
# centroid within `centroid_shift_tolerance`. Genes meeting all three carry the `shift_usable` flag; the
# count is reported with the classification, section 3 says what a gene whose *only* spatial difference
# is a displacement is called, and figure 1C shows how rare such genes are.
#
# **What counts as evidence.** Only the combination: reproducibility within each species, enough shared
# support, and adequate amplitude in both species. A gene that is flat or noisy in one species cannot be
# called displaced, however large its level difference.
#
# **2.2 asks the aggregate version of the same question.** The five quantities are per gene; the next
# subsection asks the threshold-free prior question - is there a conserved spatial architecture at all? -
# so that the classification in section 3 is read against an answer that does not depend on section 3's
# thresholds.
#

# %%
# Purpose: level, amplitude, shape, position, displacement and inversion metrics on the shared PT DPT.

mouse_mask = np.isfinite(balanced_mouse)
human_mask_shared = np.isfinite(balanced_human)


def _curve_metrics(mouse, mouse_mask_values, human, human_mask_values, x_unit):
    """Every curve metric for one (mouse, human) balanced-curve pair on the shared PT DPT."""
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

    # D - where expression sits, in each species: centroid, peak and half-maximum window.
    mouse_centroid = _row_weighted_centroid(mouse, x_unit, common)
    human_centroid = _row_weighted_centroid(human, x_unit, common)
    mouse_peak = _row_peak_position(mouse, x_unit, common)
    human_peak = _row_peak_position(human, x_unit, common)
    mouse_onset, mouse_offset, mouse_width = _row_halfmax_window(mouse, x_unit, common)
    human_onset, human_offset, human_width = _row_halfmax_window(human, x_unit, common)

    # E - the early-to-late gradient: the second positional summary, and the one whose *sign* the two
    # species can disagree on.
    mouse_early_to_late = _row_early_to_late(mouse, x_unit, common)
    human_early_to_late = _row_early_to_late(human, x_unit, common)

    # The supporting measurement - the constrained displacement search - is not one of the five.
    best_shift, best_corr, residual_rms, at_boundary = _row_shift_search(
        z_mouse, z_human, x_unit,
        CONFIG["shift_search_limit"], CONFIG["shift_step"],
        min_points=CONFIG["min_pattern_grid_points"],
    )

    return pd.DataFrame({
        "n_shared_grid_points": n_shared,
        "enough_shared_support": n_shared >= CONFIG["min_pattern_grid_points"],
        "level_effect_human_minus_mouse": level,
        "mouse_amplitude": mouse_amplitude,
        "human_amplitude": human_amplitude,
        "amplitude_difference_human_minus_mouse": amplitude_difference,
        "amplitude_log2_ratio_human_over_mouse": amplitude_log2_ratio,
        "shape_corr": shape_corr,
        "shape_spearman": shape_spearman,
        "pattern_rms_z": pattern_rms_z,
        "best_shift_human_minus_mouse": best_shift,
        "post_shift_shape_corr": best_corr,
        "shift_improvement": best_corr - shape_corr,
        "residual_rms_after_shift": residual_rms,
        "shift_at_search_boundary": at_boundary,
        "mouse_position_centroid": mouse_centroid,
        "human_position_centroid": human_centroid,
        "centroid_shift_human_minus_mouse": human_centroid - mouse_centroid,
        "mouse_peak_position": mouse_peak,
        "human_peak_position": human_peak,
        "peak_shift_human_minus_mouse": human_peak - mouse_peak,
        "mouse_early_to_late": mouse_early_to_late,
        "human_early_to_late": human_early_to_late,
        "gradient_change_human_minus_mouse": human_early_to_late - mouse_early_to_late,
        "mouse_halfmax_onset": mouse_onset,
        "mouse_halfmax_offset": mouse_offset,
        "mouse_halfmax_width": mouse_width,
        "human_halfmax_onset": human_onset,
        "human_halfmax_offset": human_offset,
        "human_halfmax_width": human_width,
    })


def _spatial_discovery_score(frame, config):
    """Interpretability-gated divergence score, in units of each component's own threshold.

    Each term counts **only where the evidence for it is interpretable**, which is the difference
    between measuring divergence and measuring noise:

    * amplitude change counts when the comparison is interpretable: both species reproducibly
      patterned, or one side patterned and the other flat in every specimen/slice. The ratio of two
      flat curves is a division artefact - but the *loss* of zonation in one species is real evidence,
      and requiring both curves to correlate would score that case zero;
    * displacement counts only when every displacement criterion passes (interior optimum, agreement
      with the mouse curve, corroboration by the independent positional centroid);
    * shape divergence counts only when both species are patterned and reproducible.

    A gene with nothing interpretable therefore scores 0 rather than accumulating a large value out of
    tiny amplitudes, boundary shifts and noisy correlations. The frame must already carry the
    phenotype flags (section 3), so the score is computed after them.
    """
    amplitude_term = (frame["amplitude_log2_ratio_human_over_mouse"].abs()
                      / config["amplitude_change_major"]).where(frame["amplitude_evidence"], 0.0)
    phase_term = (frame["best_shift_human_minus_mouse"].abs()
                  / SHIFT_INTERIOR_LIMIT).where(frame["shift_usable"], 0.0)
    shape_term = (1.0 - frame["shape_corr"].clip(-1.0, 1.0)).where(frame["positional_usable"], 0.0)

    return (amplitude_term.fillna(0.0) + phase_term.fillna(0.0) + shape_term.fillna(0.0))


def _assemble_metrics(mouse_curves, human_curves):
    """The per-gene metric table on the shared PT DPT, joined onto the curve table.

    The primary analysis and every sensitivity variant run through this one function, so a variant
    can only differ through its curves or its thresholds - never through a different code path.
    """
    mouse_mask_values = np.isfinite(mouse_curves)
    human_mask_values = np.isfinite(human_curves)

    metrics = _curve_metrics(
        mouse_curves, mouse_mask_values, human_curves, human_mask_values, grid_unit
    )

    return (
        pd.DataFrame({"gene": gene_names})
        .join(metrics)
        .merge(curve_table, on="gene", how="left")
    )


gene_metrics = _assemble_metrics(balanced_mouse, balanced_human)

_report("genes with comparable curves on the shared PT DPT",
        int(gene_metrics["enough_shared_support"].sum()), len(gene_names),
        why=f"both species support >= {CONFIG['min_pattern_grid_points']} shared grid points")
_report("genes with a rank-based shape correlation",
        int(np.isfinite(gene_metrics["shape_spearman"]).sum()), len(gene_names),
        why="Spearman needs every shared grid point supported by both species; NaN elsewhere")

display(gene_metrics[[
    "gene", "level_effect_human_minus_mouse", "mouse_amplitude", "human_amplitude",
    "amplitude_log2_ratio_human_over_mouse", "shape_corr", "best_shift_human_minus_mouse",
    "post_shift_shape_corr", "centroid_shift_human_minus_mouse",
]].head(10).round(4))


# %% [markdown]
# ### 2.2 - Is there a conserved cross-species spatial architecture?
#
# **Why this section exists.** Everything up to here has been decomposing differences. Before describing
# the exceptions, the analysis owes a threshold-independent answer to a prior question: *is there a
# conserved spatial architecture at all?* The phenotype classes cannot answer it - "172 genes are
# conserved zonation" depends on where the thresholds were put - so this subsection asks the question
# directly, on the genes that are reproducibly patterned in **both** species.
#
# **How.** Two statistics over that set, each with a null obtained by permuting human gene identities
# (i.e. breaking the ortholog pairing while keeping each species' own distribution of curves intact):
# the median matched human-mouse curve correlation after within-species standardisation, and the rank
# correlation between mouse and human positional centroids. The permutation destroys exactly the thing
# being tested - which mouse curve belongs with which human curve - and nothing else, so it is the right
# null for a conservation claim.
#
# **What counts as evidence.** An observed statistic far outside the permutation distribution, reported
# as an empirical p-value `(1 + #{null >= observed}) / (1 + n_permutations)` and as the effect size
# against the null's spread. A conserved architecture at the level of whole pathways would be a weaker
# claim than this one, which deliberately stays at the level of individual curves.
#
# **How to read the result.** The claim this test supports is narrow and should be stated that way:
# ortholog identity preserves PT positional information **above** what random pairing, or pairing within
# the same broad PT territory, would give - but the absolute concordance is modest and heterogeneous
# across genes, not a demonstration that mouse and human PT trajectories are near-identical. Report the
# observed statistic, both nulls and the effect size together; "conserved architecture" without the
# numbers overstates it.
#
# **Limits.** The null presumes gene exchangeability within the patterned set: it does not model the
# possibility that a few strong genes drive the median, so the per-gene distribution is shown in the
# figure rather than summarised only by its median. It is also strictly a statement about this cohort -
# two mouse specimens and one human donor - and about continuous *shape* and *position*, not about
# levels. Positive results here license the phrase "conserved spatial architecture" for these data and
# nothing wider.
#

# %%
# Purpose: the conserved-architecture test - matched curve correlation and centroid concordance vs a null.

architecture_mask = (
    gene_metrics["mouse_amplitude"].ge(CONFIG["amplitude_patterned"])
    & gene_metrics["human_amplitude"].ge(CONFIG["amplitude_patterned"])
    & gene_metrics["mouse_reproducibility"].ge(CONFIG["within_species_corr"])
    & gene_metrics["human_reproducibility"].ge(CONFIG["within_species_corr"])
    & gene_metrics["enough_shared_support"]
    & gene_metrics["detected_in_both_species"]
    & ~gene_metrics["axis_basis_gene"]
    & ~gene_metrics["technical_gene"]
)

architecture_index = np.flatnonzero(architecture_mask.to_numpy())

_report("genes reproducibly patterned in both species", architecture_index.size, len(gene_names),
        why="the set the conservation question is asked about; axis, technical and non-comparable genes excluded")

if architecture_index.size < CONFIG["null_min_genes"]:
    print(f"Fewer than {CONFIG['null_min_genes']} genes qualify: the conservation test is skipped "
          "rather than run on a set too small to answer the question.")
    architecture_null = pd.DataFrame()
else:
    mouse_z_architecture = _row_zscore(balanced_mouse[architecture_index])
    human_z_architecture = _row_zscore(balanced_human[architecture_index])

    def _matched_median_correlation(human_matrix):
        """Median standardised correlation between each mouse curve and the human curve paired with it."""
        mask = _finite_row_mask(mouse_z_architecture, human_matrix)
        correlation = _row_pearson(mouse_z_architecture, human_matrix, mask=mask,
                                   min_points=CONFIG["min_pattern_grid_points"])
        return float(np.nanmedian(correlation))

    observed_curve_correlation = _matched_median_correlation(human_z_architecture)

    mouse_centroids = gene_metrics["mouse_position_centroid"].to_numpy()[architecture_index]
    human_centroids = gene_metrics["human_position_centroid"].to_numpy()[architecture_index]
    observed_centroid_rho = float(spearmanr(mouse_centroids, human_centroids,
                                            nan_policy="omit").correlation)
    observed_centroid_shift = float(np.nanmedian(np.abs(human_centroids - mouse_centroids)))

    # Two nulls. The exchangeable one breaks the ortholog pairing across the whole set - which mouse
    # curve belongs with which human curve - and nothing else. The *stratified* one permutes human
    # identities only within mouse-centroid terciles, so the question becomes "do orthologs match
    # better than genes that already occupy roughly the same PT territory", which is a much harder bar
    # for a positional claim (and the one a reviewer will ask for).
    mouse_tertile_edges = np.nanquantile(mouse_centroids, [1.0 / 3.0, 2.0 / 3.0])
    mouse_tertiles = np.digitize(mouse_centroids, mouse_tertile_edges)
    tertile_rows = [np.flatnonzero(mouse_tertiles == level)
                    for level in np.unique(mouse_tertiles)]

    null_kinds = ["exchangeable", "stratified within mouse-centroid tertiles"]
    null_rng = np.random.default_rng(CONFIG["null_seed"])
    null_curve = {kind: np.empty(CONFIG["null_permutations"]) for kind in null_kinds}
    null_centroid = {kind: np.empty(CONFIG["null_permutations"]) for kind in null_kinds}

    for permutation in range(CONFIG["null_permutations"]):
        orders = {"exchangeable": null_rng.permutation(architecture_index.size)}
        stratified = np.arange(architecture_index.size)
        for rows in tertile_rows:
            stratified[rows] = rows[null_rng.permutation(rows.size)]
        orders[null_kinds[1]] = stratified

        for kind, order in orders.items():
            null_curve[kind][permutation] = _matched_median_correlation(human_z_architecture[order])
            null_centroid[kind][permutation] = float(spearmanr(
                mouse_centroids, human_centroids[order], nan_policy="omit").correlation)

    def _empirical_p(observed, null):
        return float((1 + int((null >= observed).sum())) / (1 + null.size))

    architecture_null = pd.DataFrame([
        {
            "statistic": statistic,
            "null": kind,
            "observed": observed,
            "null_median": float(np.median(null_values)),
            "null_sd": float(np.std(null_values)),
            "observed_minus_null_median": observed - float(np.median(null_values)),
            "empirical_p": _empirical_p(observed, null_values),
        }
        for statistic, observed, nulls in [
            ("median matched human-mouse curve correlation (standardised)",
             observed_curve_correlation, null_curve),
            ("mouse-human positional centroid rank correlation",
             observed_centroid_rho, null_centroid),
        ]
        for kind, null_values in nulls.items()
    ])
    architecture_null.insert(0, "n_genes", architecture_index.size)
    architecture_null["n_permutations"] = CONFIG["null_permutations"]
    architecture_null["median_abs_centroid_shift"] = observed_centroid_shift

    display(architecture_null.round(4))
    _save_table(architecture_null, "conserved_architecture_null.csv")

    # Explicit equal-width geometry keeps the three null panels' plot rectangles and gutters aligned
    # even though their y-axis labels and legends differ in length.
    figure, axes = plt.subplots(1, 3, figsize=(10.2, 3.5))
    figure.subplots_adjust(left=0.075, right=0.985, bottom=0.20, top=0.76, wspace=0.32)

    for kind, colour in zip(null_kinds, ["0.75", "0.45"]):
        axes[0].hist(null_curve[kind], bins=30, color=colour, alpha=0.9,
                     label=f"{kind} (n={CONFIG['null_permutations']})")
    axes[0].axvline(observed_curve_correlation, color=PHENOTYPE_COLORS["conserved zonation"], lw=2.2,
                    label=f"observed {observed_curve_correlation:.3f}")
    axes[0].set_xlabel("Median matched curve correlation")
    axes[0].set_ylabel("Permutations")
    axes[0].set_title("Curve shape: observed vs both nulls\n"
                      f"p = {_empirical_p(observed_curve_correlation, null_curve[null_kinds[0]]):.3g} "
                      "(exchangeable), "
                      f"{_empirical_p(observed_curve_correlation, null_curve[null_kinds[1]]):.3g} "
                      "(stratified)",
                      loc="left", fontsize=9)
    axes[0].legend(fontsize=7)

    for kind, colour in zip(null_kinds, ["0.75", "0.45"]):
        axes[1].hist(null_centroid[kind], bins=30, color=colour, alpha=0.9, label=kind)
    axes[1].axvline(observed_centroid_rho, color=PHENOTYPE_COLORS["human-zonated / mouse-flat"], lw=2.2,
                    label=f"observed {observed_centroid_rho:.3f}")
    axes[1].set_xlabel("Mouse-human centroid rank correlation")
    axes[1].set_ylabel("Permutations")
    axes[1].set_title("Positional concordance: observed vs both nulls\n"
                      f"p = {_empirical_p(observed_centroid_rho, null_centroid[null_kinds[0]]):.3g} "
                      "(exchangeable), "
                      f"{_empirical_p(observed_centroid_rho, null_centroid[null_kinds[1]]):.3g} "
                      "(stratified)",
                      loc="left", fontsize=9)
    axes[1].legend(fontsize=7)

    axes[2].plot([0, 1], [0, 1], color="0.7", lw=1, ls=":", label="identical position")
    axes[2].scatter(mouse_centroids, human_centroids, s=8, alpha=0.5, linewidths=0, color="#444444")
    axes[2].set_xlabel("Mouse positional centroid")
    axes[2].set_ylabel("Human positional centroid")
    axes[2].set_title(f"Matched positions\nmedian |shift| = {observed_centroid_shift:.3f} of PT, "
                      f"rank correlation {observed_centroid_rho:.3f}",
                      loc="left", fontsize=9)
    axes[2].legend(fontsize=7)

    for label, axis in zip(["a", "b", "c"], axes):
        _panel_label(axis, label)
    figure.suptitle(
        "Permutation test of cross-species spatial architecture\n"
        f"{architecture_index.size:,} genes patterned in both species; human gene identities are permuted globally or within mouse-centroid tertiles",
        fontsize=9.2,
    )
    _save_figure(figure, "conserved_architecture_null.png")
    plt.show()


# %% [markdown]
# ## 3 - Spatial phenotypes
#
# **The classes are a reading aid, not a test.** They compress the metric families into the categories
# the biological question needs — conserved zonation, weaker/stronger zonation in human,
# mouse-zonated/human-flat, human-zonated/mouse-flat, gradient inversion, complex shape rewiring,
# weak/uncertain — and every one of them is a conjunction of thresholds, evaluated in a fixed priority
# order. The component metrics stay in the atlas, so a reader who disagrees with a threshold can
# re-derive a different partition without refitting anything.
#
# **Priority order** (first match wins, highest priority listed last in the code so it overwrites):
# excluded (not measured in both species, technical, or axis-basis) -> insufficient shared support ->
# mouse-zonated/human-flat -> human-zonated/mouse-flat -> gradient inversion -> conserved zonation ->
# weaker/stronger zonation in human -> complex shape rewiring -> weak/uncertain.
#
# **There is no "shifted" class, by decision.** A displacement is a supporting measurement here, not a
# phenotype: the classes above describe *how* a spatial pattern differs (lost, gained, inverted,
# reshaped), and a gene whose only difference is where its programme sits is reported through
# `best_shift_human_minus_mouse` and the `shift_usable` flag instead. Such a gene is classed
# `weak / uncertain zonation` — it is patterned in both species, but none of the classes describes it —
# so that bucket holds two different situations and is always read with the flag beside it. The count of
# corroborated displacements is reported with the class counts.
#
# **How the thresholds were chosen.** The next cell reports the empirical distributions the thresholds
# act on (amplitude, within-species reproducibility, shape correlation, displacement) and how many
# genes each one selects. `CONFIG` holds the values; the distributions are printed so the choice is
# auditable rather than inherited. They start from 03/05's operating values and were then checked
# against these distributions.
#
# **Guardrails encoded here.** A gene is never called spatially rewired when one species is essentially
# flat or noisy: every shape and inversion class requires reproducibility within both species and
# adequate amplitude in both. Flatness for the species-specific classes must be visible in **each**
# specimen/slice, not only in the species mean — averaging two slices can flatten a curve that neither
# slice has flat. And a large level effect never blocks a spatial label: level is not used in the
# classification at all.
#
# **3.2 and 3.3.** Subsection 3.2 re-derives every class under the analysis choices (specimen balancing,
# leave-one-out, threshold perturbations) and reports how stable a gene's class is; 3.3 shows the results
# as figures.
#

# %%
# Purpose: the empirical distributions the section-3 thresholds act on, and how many genes each selects.

DIAGNOSTIC_QUANTILES = [0.5, 0.75, 0.9, 0.95, 0.99]

diagnostic_columns = [
    "mouse_amplitude", "human_amplitude", "mouse_reproducibility", "human_reproducibility",
    "shape_corr", "post_shift_shape_corr", "amplitude_log2_ratio_human_over_mouse",
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
    "gradient inversion",
    "complex shape rewiring",
    "weak / uncertain zonation",
]

# Fine class -> broad family, used for the robustness scoring in section 3.2.
BROAD_PHENOTYPE = {
    "conserved zonation": "conserved",
    "weaker zonation in human": "amplitude-change",
    "stronger zonation in human": "amplitude-change",
    "mouse-zonated / human-flat": "species-specific zonation",
    "human-zonated / mouse-flat": "species-specific zonation",
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
    # An amplitude *difference* is evidence when the comparison is interpretable, which is not the
    # same as requiring both species to correlate. Two ways to be interpretable, and they must be
    # symmetric: (i) both species reproducibly patterned, or (ii) one side reproducibly patterned and
    # the other flat in *every* specimen/slice. Without (ii), a gene that loses its zonation entirely
    # would score zero - the flat side's curve correlation is undefined precisely because it is flat,
    # which is the evidence, not a failure to measure.
    both_sides_patterned_and_reproducible = (
        flags["strong_mouse"] & flags["strong_human"]
        & flags["reproducible_mouse"] & flags["reproducible_human"]
    )
    mouse_zonated_human_flat = (
        flags["strong_mouse"] & flags["reproducible_mouse"] & flags["flat_human_every_slice"]
    )
    human_zonated_mouse_flat = (
        flags["strong_human"] & flags["reproducible_human"] & flags["flat_mouse_every_slice"]
    )
    flags["amplitude_evidence"] = flags["enough_support"] & (
        both_sides_patterned_and_reproducible | mouse_zonated_human_flat | human_zonated_mouse_flat
    )

    flags["inversion"] = (
        flags["positional_usable"]
        & (np.sign(metrics["mouse_early_to_late"]) != np.sign(metrics["human_early_to_late"]))
        & metrics["mouse_early_to_late"].abs().ge(config["gradient_min_abs"])
        & metrics["human_early_to_late"].abs().ge(config["gradient_min_abs"])
        & metrics["shape_corr"].le(config["inversion_max_corr"])
    )

    # A displacement counts only when it is an *interior* optimum (not pinned to the search
    # boundary, where the search has only established that agreement keeps improving outwards) and
    # it is corroborated by the independent positional centroid.
    centroid_shift = metrics["centroid_shift_human_minus_mouse"]

    flags["shift_interior"] = shift.abs().le(SHIFT_INTERIOR_LIMIT)
    flags["shift_at_boundary"] = metrics["shift_at_search_boundary"]
    flags["shift_centroid_agrees"] = (
        (np.sign(centroid_shift) == np.sign(shift))
        & (centroid_shift - shift).abs().le(config["centroid_shift_tolerance"])
    )

    # `shift_usable` is a *supporting-metric* flag, not a phenotype class: it marks the genes whose
    # measured displacement clears every criterion above, so the displacement can be reported beside
    # the partition (and contribute to the divergence score) without being part of it. A gene whose
    # only spatial difference is a displacement therefore falls through to `weak / uncertain
    # zonation` - it is patterned in both species, but the classes below do not describe *how*.
    flags["shift_usable"] = (
        flags["positional_usable"]
        & flags["shift_interior"]
        & ~flags["shift_at_boundary"]
        & flags["shift_centroid_agrees"]
        & metrics["post_shift_shape_corr"].ge(config["post_shift_corr"])
        & shift.abs().ge(config["minimum_shift"])
        & metrics["shift_improvement"].ge(config["shift_improvement"])
    )

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
gene_metrics["spatial_discovery_score"] = _spatial_discovery_score(gene_metrics, CONFIG)

_report("genes carrying a corroborated displacement (supporting metric, not a class)",
        int(gene_metrics["shift_usable"].sum()),
        why="interior optimum, centroid-corroborated, and it improves the shape agreement")

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
        "median_post_shift_shape_corr": float(np.nanmedian(
            gene_metrics.loc[gene_metrics["spatial_phenotype"] == label, "post_shift_shape_corr"])),
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
        best_shift_human_minus_mouse=0.0,
        post_shift_shape_corr=0.95,
        shift_improvement=0.01,
        shift_at_search_boundary=False, centroid_shift_human_minus_mouse=0.0,
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
                 post_shift_shape_corr=0.2,
                 **{f"amplitude__{sample}": 0.03 for sample in human_samples}),
            "mouse-zonated / human-flat"),
        "expected human-zonated / mouse-flat": (
            gene(mouse_amplitude=0.03, amplitude_log2_ratio_human_over_mouse=3.7, shape_corr=0.1,
                 post_shift_shape_corr=0.2,
                 **{f"amplitude__{sample}": 0.03 for sample in mouse_samples}),
            "human-zonated / mouse-flat"),
        # A displacement no longer defines a class. These two rows fix the demotion: the flag is
        # set, and the gene is *not* labelled phase - it falls to the bucket for genes the classes do
        # not describe.
        "a corroborated, interior displacement": (
            gene(shape_corr=0.4, best_shift_human_minus_mouse=0.15,
                 centroid_shift_human_minus_mouse=0.15,
                 post_shift_shape_corr=0.93, shift_improvement=0.4),
            "weak / uncertain zonation"),
        "the same displacement pointing the other way": (
            gene(shape_corr=0.4, best_shift_human_minus_mouse=-0.15,
                 centroid_shift_human_minus_mouse=-0.15,
                 post_shift_shape_corr=0.93, shift_improvement=0.4),
            "weak / uncertain zonation"),
        "an optimum pinned to the search boundary is not a displacement": (
            gene(shape_corr=0.4, best_shift_human_minus_mouse=0.35,
                 centroid_shift_human_minus_mouse=0.35, shift_at_search_boundary=True,
                 post_shift_shape_corr=0.93, shift_improvement=0.4),
            "complex shape rewiring"),
        "an interior shift the centroid contradicts is not interpreted": (
            gene(shape_corr=0.4, best_shift_human_minus_mouse=0.16,
                 centroid_shift_human_minus_mouse=-0.16,
                 post_shift_shape_corr=0.93, shift_improvement=0.4),
            "complex shape rewiring"),
        "expected gradient inversion": (
            gene(shape_corr=-0.8, mouse_early_to_late=0.4, human_early_to_late=-0.4),
            "gradient inversion"),
        "expected complex rewiring": (
            gene(shape_corr=0.1, post_shift_shape_corr=0.3), "complex shape rewiring"),
        "flat in both species": (
            gene(mouse_amplitude=0.02, human_amplitude=0.02, shape_corr=np.nan,
                 **{f"amplitude__{sample}": 0.02 for sample in mouse_samples + human_samples}),
            "weak / uncertain zonation"),
        "flat in one human slice only is not human-flat": (
            gene(human_amplitude=0.03, amplitude_log2_ratio_human_over_mouse=-3.7, shape_corr=0.1,
                 post_shift_shape_corr=0.2, human_reproducibility=np.nan,
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
                 "enough_shared_support", "shift_at_search_boundary"):
        frame[flag] = frame[flag].astype(bool)

    flags = _phenotype_flags(frame, CONFIG)
    # The score check joins the flags onto the metrics frame positionally, so `flags` keeps its default
    # index; these two rows look a flag up by the synthetic case's name instead.
    flags_by_gene = flags.set_axis(frame["gene"].to_numpy())
    assigned = pd.Series(_assign_phenotypes(flags), index=frame["gene"])

    outcome = pd.DataFrame({
        "synthetic_case": list(synthetic),
        "expected_class": [expected for _, expected in synthetic.values()],
        "assigned_class": [str(assigned[name]) for name in synthetic],
    })
    outcome["passed"] = outcome["expected_class"] == outcome["assigned_class"]

    # The divergence score must be gated by interpretability: a gene with nothing interpretable scores
    # exactly 0 instead of accumulating tiny-amplitude ratios, boundary shifts and noisy correlations.
    # The score is gated by the flags, so they are joined onto the metrics frame first -
    # exactly the order the notebook itself uses.
    frame_with_flags = frame.join(flags)
    scores = _spatial_discovery_score(frame_with_flags, CONFIG).set_axis(
        frame["gene"].to_numpy()
    )
    uninterpretable = ["flat in both species", "too little shared support"]
    outcome = pd.concat([outcome, pd.DataFrame([
        {"synthetic_case": "score: no interpretable component scores exactly 0",
         "expected_class": "0", "assigned_class": f"{float(scores[uninterpretable].abs().max()):.3g}",
         "passed": bool(np.isclose(scores[uninterpretable].abs().max(), 0.0))},
        {"synthetic_case": "score: a corroborated displacement contributes",
         "expected_class": "> 0",
         "assigned_class": f"{float(scores['a corroborated, interior displacement']):.3g}",
         "passed": bool(scores["a corroborated, interior displacement"] > 0)},
        {"synthetic_case": "flags: a corroborated displacement sets shift_usable",
         "expected_class": "True",
         "assigned_class": str(bool(flags_by_gene.loc["a corroborated, interior displacement", "shift_usable"])),
         "passed": bool(flags_by_gene.loc["a corroborated, interior displacement", "shift_usable"])},
        {"synthetic_case": "flags: a boundary optimum never sets shift_usable",
         "expected_class": "False",
         "assigned_class": str(bool(flags_by_gene.loc[
             "an optimum pinned to the search boundary is not a displacement", "shift_usable"])),
         "passed": not bool(flags_by_gene.loc[
             "an optimum pinned to the search boundary is not a displacement", "shift_usable"])},
        {"synthetic_case": "score: complex rewiring contributes via shape divergence",
         "expected_class": "> 0", "assigned_class": f"{float(scores['expected complex rewiring']):.3g}",
         "passed": bool(scores["expected complex rewiring"] > 0)},
        {"synthetic_case": "score: mouse-zonated / human-flat keeps its amplitude evidence",
         "expected_class": "> 0",
         "assigned_class": f"{float(scores['expected mouse-zonated / human-flat']):.3g}",
         "passed": bool(scores["expected mouse-zonated / human-flat"] > 0)},
        {"synthetic_case": "score: human-zonated / mouse-flat keeps its amplitude evidence",
         "expected_class": "> 0",
         "assigned_class": f"{float(scores['expected human-zonated / mouse-flat']):.3g}",
         "passed": bool(scores["expected human-zonated / mouse-flat"] > 0)},
        {"synthetic_case": "score: the two species-specific directions are symmetric",
         "expected_class": "equal",
         "assigned_class": f"{float(scores['expected mouse-zonated / human-flat'] / scores['expected human-zonated / mouse-flat']):.3f}",
         "passed": bool(np.isclose(scores["expected mouse-zonated / human-flat"],
                                   scores["expected human-zonated / mouse-flat"],
                                   rtol=0.05))},
        {"synthetic_case": "score: a conserved gene scores below a rewired one",
         "expected_class": "< complex",
         "assigned_class": f"{float(scores['expected conserved']):.3g}",
         "passed": bool(scores["expected conserved"] < scores["expected complex rewiring"])},
    ])], ignore_index=True)
    return outcome


phenotype_selfcheck = _phenotype_selfcheck()
display(phenotype_selfcheck)

if not phenotype_selfcheck["passed"].all():
    failed = phenotype_selfcheck.loc[~phenotype_selfcheck["passed"], "synthetic_case"].tolist()
    raise AssertionError(f"synthetic phenotype self-check failed for: {failed}")

print(f"All {len(phenotype_selfcheck)} synthetic phenotype self-checks passed.")


# %% [markdown]
# ### 3.2 - Robustness of the phenotypes across the analysis choices
#
# **Why.** A phenotype class is a conjunction of thresholds applied to curves that were themselves built
# under choices (specimen balancing, the detection threshold). Reporting one partition and calling it the
# result would hide that, so the whole classification is re-derived under each choice and the per-gene
# agreement is kept.
#
# **The variants.** (1) pooled tubule-weighted curves instead of equal-weight specimen curves; (2) leave
# out each mouse specimen (two variants); (3) leave out each human slice (two variants - *slice
# sensitivity*, not replication); (4) amplitude thresholds scaled by 0.8 and 1.25; (5) the displacement
# threshold scaled by 0.75 and 1.25; (6) detection threshold at 2% and 10%; (7) axis-basis genes allowed
# into discovery. Twelve variants in total, beside the primary analysis.
#
# **How they are run.** All variants share one code path (`_assemble_metrics` + `_phenotype_flags` +
# `_assign_phenotypes`), so a variant can only differ through its curves or its thresholds. The
# per-variant class counts are written to `tables/sensitivity_summary.csv`; the table shown here keeps
# only the stability columns.
#
# **What the numbers mean, and what they do not.** `robustness_score` is the fraction of variants that
# keep a gene's **broad** phenotype (conserved / amplitude-change / species-specific zonation / inversion
# / complex rewiring / no spatial signal / excluded). It is a stability score, not a significance
# measure. A gene that leaves the universe in a stricter-detection variant is counted as not retaining
# its phenotype; that is deliberate, because "no longer measured" is a real change in what can be
# claimed. Nothing here walks through every variant: the per-gene score is in the atlas, and the
# fraction of each class that retains its broad class is reported once.
#
# **The one-donor caveat, stated plainly.** Leaving out a human slice cannot be biological replication
# with a single donor, and no variant here produces a p-value.
#

# %%
# Purpose: re-derive every phenotype under each analysis choice, and score per-gene stability.


def _sensitivity_variants():
    """Every variant the phenotypes are re-derived under: the curves and the thresholds."""
    variants = {
        "primary (shared DPT, specimen-balanced)": {
            "mouse": balanced_mouse, "human": balanced_human, "config": CONFIG,
            "note": "the primary analysis",
        },
        "pooled tubule-weighted curves": {
            "mouse": pooled_mouse, "human": pooled_human, "config": CONFIG,
            "note": "pooled fit instead of equal-weight specimen curves",
        },
    }

    for name in mouse_samples:
        remaining = [specimen for specimen in mouse_samples if specimen != name]
        variants[f"leave out mouse specimen {name}"] = {
            "mouse": specimen_balanced_curves(
                {specimen: specimen_curves[specimen] for specimen in remaining}
            ),
            "human": balanced_human, "config": CONFIG,
            "note": "specimen sensitivity",
        }

    for name in human_samples:
        remaining = [specimen for specimen in human_samples if specimen != name]
        variants[f"leave out human slice {name}"] = {
            "mouse": balanced_mouse,
            "human": specimen_balanced_curves(
                {specimen: specimen_curves[specimen] for specimen in remaining}
            ),
            "config": CONFIG,
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
            "mouse": balanced_mouse, "human": balanced_human,
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
    variant_metrics = _assemble_metrics(variant_spec["mouse"], variant_spec["human"])

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
    # The divergence score is gated by the flags, so it is computed after them.
    variant_metrics = variant_metrics.join(variant_flags)
    variant_metrics["spatial_discovery_score"] = _spatial_discovery_score(
        variant_metrics, variant_spec["config"]
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
        "same_fine_phenotype_share": float((variant_labels.to_numpy() == primary_labels).mean()),
        "same_broad_phenotype_share": float((variant_broad == primary_broad).mean()),
        "n_discovery_eligible": int((~variant_flags["excluded"] & variant_flags["enough_support"]).sum()),
        "n_conserved": int((variant_labels == "conserved zonation").sum()),
        "n_amplitude_change": int(variant_labels.isin(
            ["weaker zonation in human", "stronger zonation in human"]).sum()),
        "n_species_specific": int(variant_labels.isin(
            ["mouse-zonated / human-flat", "human-zonated / mouse-flat"]).sum()),
        "n_inversion": int((variant_labels == "gradient inversion").sum()),
        "n_complex": int((variant_labels == "complex shape rewiring").sum()),
        "spearman_discovery_score_vs_primary": float(
            spearmanr(primary_score, variant_score, nan_policy="omit").correlation),
        "top50_discovery_overlap": len(
            set(gene_names[order_primary[:50]]) & set(gene_names[order_variant[:50]])
        ) / 50.0,
    })

sensitivity_summary = pd.DataFrame(sensitivity_summary_rows)

# The per-class counts of every variant are written to `tables/`; only the stability columns are
# displayed, because they are the ones the robustness claim is read from.
STABILITY_COLUMNS = [
    "variant", "note", "same_broad_phenotype_share", "same_fine_phenotype_share",
    "spearman_discovery_score_vs_primary", "top50_discovery_overlap",
]
display(sensitivity_summary[STABILITY_COLUMNS].round(3))
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
        why="used to prioritise the genes shown in sections 3.3 and 4.2")


# %% [markdown]
# ### 3.3 - The result figures
#
# Four figures, each answering one part of the question (the curve-layer QC panel is figure 1A, at the
# end of section 1):
#
# * **1B** — where does zonation strength sit relative to the identity line, and which genes have strong
#   zonation in both species but poor shape agreement (the rewiring corner)?
# * **1C** — displacement versus residual mismatch: is a gene explained by a positional change, or is it
#   genuinely reshaped after the best possible displacement? A displacement is drawn only for genes whose
#   criteria are met, in the supporting-measurement colour.
# * **1D** — representative genes per class, with **both** specimens/slices drawn faintly behind the
#   balanced curves, so within-species consistency is visible rather than asserted.
# * **1E** — shape only: each gene z-scored *within each species*, so the heatmap cannot be driven by an
#   abundance difference between human and mouse.
#
# Genes are chosen by class-specific evidence (the ranking cell below), preferring genes that kept their
# broad class in at least 80% of the 3.2 variants. Figures are readable without the code: every panel
# states its class, sample sizes and direction convention.
#

# %%
# Purpose: the per-class exemplar ranking that the figure cells below use.

# What each phenotype's exemplars must demonstrate, and therefore how they are ranked. A single
# global divergence score picks genes that score high for *any* reason, so each class ranks by the
# evidence it is supposed to illustrate.
CLASS_EVIDENCE = {
    "conserved zonation": {
        "label": "shape agreement first, then the smallest amplitude and position changes",
        "filter": lambda block: block["shape_corr"].ge(CONFIG["conserved_corr"]),
        "rank": ["shape_corr", "pattern_rms_z", "best_shift_human_minus_mouse"],
        "ascending": [False, True, True],
    },
    "weaker zonation in human": {
        "label": "largest amplitude loss with the shape intact and the position unchanged",
        "filter": lambda block: block["shape_corr"].ge(CONFIG["conserved_corr"])
        & block["best_shift_human_minus_mouse"].abs().lt(CONFIG["minimum_shift"]),
        "rank": ["amplitude_log2_ratio_human_over_mouse", "shape_corr"],
        "ascending": [True, False],
    },
    "stronger zonation in human": {
        "label": "largest amplitude gain with the shape intact and the position unchanged",
        "filter": lambda block: block["shape_corr"].ge(CONFIG["conserved_corr"])
        & block["best_shift_human_minus_mouse"].abs().lt(CONFIG["minimum_shift"]),
        "rank": ["amplitude_log2_ratio_human_over_mouse", "shape_corr"],
        "ascending": [False, False],
    },
    "mouse-zonated / human-flat": {
        "label": "largest mouse-to-human amplitude contrast",
        "filter": lambda block: block["shift_at_search_boundary"].eq(False),
        "rank": ["amplitude_log2_ratio_human_over_mouse", "mouse_amplitude"],
        "ascending": [True, False],
    },
    "human-zonated / mouse-flat": {
        "label": "largest human-to-mouse amplitude contrast",
        "filter": lambda block: block["shift_at_search_boundary"].eq(False),
        "rank": ["amplitude_log2_ratio_human_over_mouse", "human_amplitude"],
        "ascending": [False, False],
    },
    "gradient inversion": {
        "label": "strongest opposition between the two early-to-late gradients",
        "filter": lambda block: block["mouse_early_to_late"].abs().ge(CONFIG["gradient_min_abs"])
        & block["human_early_to_late"].abs().ge(CONFIG["gradient_min_abs"]),
        "rank": ["gradient_change_human_minus_mouse", "mouse_early_to_late"],
        "ascending": [False, False],
    },
    "complex shape rewiring": {
        "label": "largest residual mismatch after the best displacement",
        "filter": lambda block: block["shape_corr"].lt(CONFIG["conserved_corr"]),
        "rank": ["residual_rms_after_shift", "shape_corr"],
        "ascending": [False, True],
    },
}


def _top_genes_for_class(label, n, prefer_robust=True):
    """Top `n` genes of a class by that class's own evidence, preferring variant-stable genes.

    Class-specific rather than one global score: a gene can score high for *any* reason - including an
    amplitude difference - and still be a poor illustration of a displacement. Each class ranks by the
    evidence it is meant to demonstrate; if the stricter filter empties the class, that is reported
    rather than hidden, and the class as classified is used.
    """
    block = gene_metrics[
        gene_metrics["spatial_phenotype"].eq(label)
        & ~gene_metrics["axis_basis_gene"]
        & ~gene_metrics["technical_gene"]
        & gene_metrics["detected_in_both_species"]
    ]

    rule = CLASS_EVIDENCE.get(label)
    if rule is None:
        block = block.sort_values("spatial_discovery_score", ascending=False)
    else:
        strict = block[rule["filter"](block)]
        if len(strict):
            block = strict.sort_values(rule["rank"], ascending=rule["ascending"])
        else:
            print(f"  note: no {label} gene meets the stricter exemplar criteria "
                  f"({rule['label']}); falling back to the class as classified")
            block = block.sort_values("spatial_discovery_score", ascending=False)

    if prefer_robust and len(block):
        stable = block[block["robustness_score"].ge(0.8)]
        if len(stable) >= n:
            block = stable

    return block.head(n)


discovery_eligible = (
    ~gene_metrics["axis_basis_gene"]
    & ~gene_metrics["technical_gene"]
    & gene_metrics["detected_in_both_species"]
    & gene_metrics["enough_shared_support"]
)



# %%
# Purpose: figure 1B - the zonation conservation landscape.

landscape = gene_metrics[discovery_eligible].dropna(
    subset=["mouse_amplitude", "human_amplitude", "shape_corr"]
)

strong_both = landscape[
    landscape["mouse_amplitude"].ge(CONFIG["amplitude_patterned"])
    & landscape["human_amplitude"].ge(CONFIG["amplitude_patterned"])
]
quadrants = pd.DataFrame([
    {"group": "strong in both: conserved shape", "n": int((strong_both["shape_corr"] >= CONFIG["conserved_corr"]).sum()),
     "colour": PHENOTYPE_COLORS["conserved zonation"]},
    {"group": "strong in both: intermediate shape", "n": int(((strong_both["shape_corr"] < CONFIG["conserved_corr"])
               & (strong_both["shape_corr"] > CONFIG["rewired_max_corr"])).sum()), "colour": "#8F8F8F"},
    {"group": "strong in both: rewired shape", "n": int((strong_both["shape_corr"] <= CONFIG["rewired_max_corr"]).sum()),
     "colour": PHENOTYPE_COLORS["complex shape rewiring"]},
    {"group": "mouse patterned / human flat", "n": int(((landscape["mouse_amplitude"] >= CONFIG["amplitude_patterned"])
               & (landscape["human_amplitude"] <= CONFIG["amplitude_flat"])).sum()),
     "colour": PHENOTYPE_COLORS["mouse-zonated / human-flat"]},
    {"group": "human patterned / mouse flat", "n": int(((landscape["human_amplitude"] >= CONFIG["amplitude_patterned"])
               & (landscape["mouse_amplitude"] <= CONFIG["amplitude_flat"])).sum()),
     "colour": PHENOTYPE_COLORS["human-zonated / mouse-flat"]},
])

figure = plt.figure(figsize=(8.5, 4.7), layout="constrained")
grid_spec = figure.add_gridspec(2, 2, width_ratios=[1.30, 0.82], height_ratios=[1.0, 0.85])
ax_main = figure.add_subplot(grid_spec[:, 0])
ax_margin = figure.add_subplot(grid_spec[0, 1])
ax_corr = figure.add_subplot(grid_spec[1, 1])

# All observations remain visible as a quiet background.  Colour is reserved for the subset where
# both amplitudes clear the patterned threshold, where a shape correlation is interpretable.
ax_main.scatter(landscape["mouse_amplitude"], landscape["human_amplitude"], s=5, alpha=0.12,
                linewidths=0, color="0.45", rasterized=True, zorder=1)
points = ax_main.scatter(
    strong_both["mouse_amplitude"], strong_both["human_amplitude"], c=strong_both["shape_corr"],
    cmap="RdYlBu_r", vmin=-1, vmax=1, s=8, alpha=0.70, linewidths=0, rasterized=True, zorder=2,
)
limit = max(float(np.nanmax(np.concatenate([
    landscape["mouse_amplitude"].to_numpy(), landscape["human_amplitude"].to_numpy()
])) * 1.03), 0.5)
ax_main.plot([0, limit], [0, limit], color="0.30", lw=0.8, ls=":", zorder=0)
ax_main.axvline(CONFIG["amplitude_patterned"], color="0.55", lw=0.75, ls="--", zorder=0)
ax_main.axhline(CONFIG["amplitude_patterned"], color="0.55", lw=0.75, ls="--", zorder=0)
ax_main.set(xlim=(-0.02, limit), ylim=(-0.02, limit),
            xlabel="Mouse zonation amplitude (peak-to-trough lognorm)",
            ylabel="Human zonation amplitude (peak-to-trough lognorm)")
ax_main.text(0.98, 0.04, "dotted: equal amplitude\ndashed: patterned threshold",
             transform=ax_main.transAxes, fontsize=6.4, ha="right", va="bottom", color="0.35")
ax_main.set_title(f"Zonation landscape ({len(landscape):,} discovery-eligible orthologs)",
                  loc="left", fontweight="bold")
_panel_label(ax_main, "b")

label_candidates = pd.concat([
    _top_genes_for_class("conserved zonation", 2),
    _top_genes_for_class("mouse-zonated / human-flat", 2),
    _top_genes_for_class("human-zonated / mouse-flat", 2),
    _top_genes_for_class("complex shape rewiring", 2),
])
for row in label_candidates.itertuples():
    ax_main.annotate(row.gene, (row.mouse_amplitude, row.human_amplitude), xytext=(3, 3),
                     textcoords="offset points", fontsize=6.2, color="0.20")
colour_bar = figure.colorbar(points, ax=ax_main, pad=0.015, fraction=0.048)
colour_bar.set_label("Shape correlation\n(both species patterned)", fontsize=7)
colour_bar.ax.tick_params(labelsize=6.5)

margin_positions = np.arange(len(quadrants))
ax_margin.barh(margin_positions, quadrants["n"], color=quadrants["colour"], alpha=0.88)
ax_margin.set_yticks(margin_positions)
ax_margin.set_yticklabels([textwrap.fill(group, 27) for group in quadrants["group"]], fontsize=6.5)
ax_margin.invert_yaxis()
ax_margin.set_xlabel("Genes", fontsize=7)
ax_margin.set_title("Signal composition", loc="left", fontweight="bold", fontsize=8.5)
for index, value in enumerate(quadrants["n"]):
    ax_margin.text(value, index, f" {value:,}", va="center", fontsize=6.5)
ax_margin.grid(axis="x", color="0.92", lw=0.6)

ax_corr.hist(strong_both["shape_corr"].dropna(), bins=np.linspace(-1, 1, 25), color="0.55",
             edgecolor="white", lw=0.3)
ax_corr.axvline(CONFIG["rewired_max_corr"], color=PHENOTYPE_COLORS["complex shape rewiring"],
                lw=1.0, ls="--")
ax_corr.axvline(CONFIG["conserved_corr"], color=PHENOTYPE_COLORS["conserved zonation"], lw=1.0, ls="--")
ax_corr.set(xlim=(-1, 1), xlabel="Shape correlation", ylabel="Genes")
ax_corr.set_title("Among genes patterned in both species", loc="left", fontsize=8.5, fontweight="bold")
ax_corr.tick_params(labelsize=6.5)

figure.suptitle(
    "Amplitude and spatial-shape changes are separable in the primary shared PT DPT\n"
    "Colour is shown only where both species have a patterned curve; all other eligible genes remain in grey",
    fontsize=9.5,
)
_save_figure(figure, "zonation_amplitude_human_vs_mouse.png")
plt.show()

_save_table(landscape, "zonation_landscape_eligible_genes.csv")


# %%
# Purpose: figure 1C - displacement versus residual shape mismatch.

displacement_frame = gene_metrics[discovery_eligible].dropna(
    subset=["shape_corr", "post_shift_shape_corr", "best_shift_human_minus_mouse",
            "residual_rms_after_shift"]
)

figure, axes = plt.subplots(1, 3, figsize=(10.5, 3.65), width_ratios=[1.0, 1.0, 1.0])
figure.subplots_adjust(left=0.065, right=0.99, bottom=0.20, top=0.78, wspace=0.34)
ax_raw, ax_shift_residual, ax_summary = axes

# Panel 1: does a bounded displacement *explain* the difference?
#
# Translation almost always improves the correlation somewhat, so an improvement is not evidence of a
# displacement and must not be drawn as one. Orange marks only genes whose displacement meets every
# criterion - a supporting measurement, not a phenotype class - and the merely-improved cloud stays
# grey. The visual message then supports the result rather than contradicting it: clean displacements
# are rare.
ax_raw.plot([-1, 1], [-1, 1], color="0.5", lw=1, ls=":", label="no change from shifting")
ax_raw.scatter(
    displacement_frame["shape_corr"], displacement_frame["post_shift_shape_corr"],
    s=5, alpha=0.18, linewidths=0, color="0.40", rasterized=True,
    label=f"all usable patterns (n={len(displacement_frame):,})",
)
displacement_usable = displacement_frame[displacement_frame["shift_usable"]]
ax_raw.scatter(
    displacement_usable["shape_corr"], displacement_usable["post_shift_shape_corr"],
    s=21, alpha=0.95, linewidths=0.3, edgecolors="white",
    color=DISPLACEMENT_COLOR,
    label=f"interior displacement passes all criteria (n={len(displacement_usable):,})",
)
ax_raw.set_xlabel("Standardised curve correlation (shared PT DPT)")
ax_raw.set_ylabel(
    f"Best correlation after displacement\n(|shift| ≤ {SHIFT_INTERIOR_LIMIT:.2f} of PT)"
)
ax_raw.set_title("A  |  Does bounded displacement rescue the shape?", loc="left",
                           fontweight="bold")
ax_raw.text(0.03, 0.03,
                        f"Orange: {len(displacement_usable)} genes pass every displacement criterion.",
                        transform=ax_raw.transAxes, fontsize=5.9, va="bottom", color="0.30")

# Panel 2: displacement against what is left after it.
for label in ["gradient inversion", "complex shape rewiring", "conserved zonation"]:
    block = displacement_frame[displacement_frame["spatial_phenotype"].eq(label)]
    if not len(block):
        continue
    ax_shift_residual.scatter(
        block["best_shift_human_minus_mouse"], block["residual_rms_after_shift"],
        s=10, alpha=0.75, linewidths=0, color=PHENOTYPE_COLORS[label],
        label=f"{label} (n={len(block):,})", rasterized=True,
    )
ax_shift_residual.axvline(0, color="black", lw=0.8)
ax_shift_residual.set_xlabel("Best displacement of the human program (fraction of PT)\n"
                             "earlier in human  ←     →  later in human")
ax_shift_residual.set_ylabel("Standardised RMS difference after the best displacement")
ax_shift_residual.set_title("B  |  Displacement versus residual rewiring", loc="left", fontweight="bold")
ax_shift_residual.legend(loc="lower left", fontsize=5.3, markerscale=0.8, ncols=2,
                         labelspacing=0.25, columnspacing=0.8)

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
ax_summary.barh(y_positions, class_counts.to_numpy(), color="#D3D3D3", label="all genes in class")
ax_summary.barh(y_positions, robust_counts.to_numpy(), color="#333333",
                label="class stable in ≥80% of sensitivity variants")
ax_summary.set_yticks(y_positions)
ax_summary.set_yticklabels([textwrap.fill(label, 23) for label in PHENOTYPE_CLASSES], fontsize=6.3)
ax_summary.set_xlabel("Genes")
ax_summary.set_xscale("symlog", linthresh=100)
ax_summary.set_title("C  |  Class size and stability\nlight = all; dark = stable in ≥80% of variants",
                     loc="left", fontweight="bold", fontsize=7.7)
ax_summary.invert_yaxis()
for index, (total, robust) in enumerate(zip(class_counts.to_numpy(), robust_counts.to_numpy())):
    ax_summary.text(total, index, f" {total:,}", va="center", fontsize=5.8)

for axis in axes:
    axis.grid(axis="y", color="0.92", lw=0.55, zorder=0)

figure.suptitle(
    "Positional displacement is uncommon; most low-agreement curves remain mismatched after translation\n"
    "Descriptive classifications on the shared PT DPT; no cross-species population-level inference",
    fontsize=9.2,
)
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

main_representative_table = pd.concat(
    [_top_genes_for_class(label, CONFIG["main_representative_genes_per_class"])
     for label in REPRESENTATIVE_CLASSES], ignore_index=True,
)

# A 3 × 3 plate holds one exemplar per class at a legible physical size; the two panels the seven
# classes leave over stay empty. The full per-class selection is still exported below as supplementary
# material rather than being shrunk until it can no longer support the claim it is meant to illustrate.
figure, axes = plt.subplots(3, 3, figsize=(7.25, 6.45), sharex=True, layout="constrained")
for panel_index, (axis, label) in enumerate(zip(axes.flat, REPRESENTATIVE_CLASSES)):
    row_index, column_index = divmod(panel_index, 3)
    block = main_representative_table[
        main_representative_table["spatial_phenotype"].eq(label)
    ]
    if block.empty:
        axis.set_visible(False)
        continue

    record = block.iloc[0]
    index = int(np.flatnonzero(gene_names == record["gene"])[0])
    for specimen in mouse_samples:
        axis.plot(grid_unit, specimen_curves[specimen][index], color=SPECIES_COLORS["mouse"],
                  lw=0.55, alpha=0.26, ls=":", zorder=1)
    for specimen in human_samples:
        axis.plot(grid_unit, specimen_curves[specimen][index], color=SPECIES_COLORS["human"],
                  lw=0.55, alpha=0.26, ls=":", zorder=1)
    axis.plot(grid_unit, balanced_mouse[index], color=SPECIES_COLORS["mouse"], lw=1.65, zorder=3)
    axis.plot(grid_unit, balanced_human[index], color=SPECIES_COLORS["human"], lw=1.65, ls="--", zorder=3)
    axis.set_title(f"{record['gene']}  |  {label}", loc="left", fontsize=7.1, fontweight="bold")
    axis.text(0.02, 0.04,
              f"Aₘ={record['mouse_amplitude']:.2f}; Aₕ={record['human_amplitude']:.2f}; "
              f"r={record['shape_corr']:.2f}; Δpos={record['best_shift_human_minus_mouse']:+.2f}\n"
              f"stability={record['robustness_score']:.2f}",
              transform=axis.transAxes, fontsize=5.6, va="bottom", color="0.25",
              bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 1.0})
    _style_pt_axis(axis, xlabel=row_index == 2,
                   ylabel="Fitted lognorm" if column_index == 0 else None)
    axis.tick_params(labelsize=6.2)

_panel_label(axes[0, 0], "d")
axes[0, 0].legend(
    handles=[
        Line2D([], [], color=SPECIES_COLORS["mouse"], lw=1.65, label="mouse, balanced"),
        Line2D([], [], color=SPECIES_COLORS["human"], lw=1.65, ls="--", label="human, balanced"),
        Line2D([], [], color="0.45", lw=0.7, ls=":", label="individual specimen/slice"),
    ], loc="upper right", fontsize=5.7,
)
figure.suptitle(
    "One high-evidence gene per spatial phenotype\n"
    "Thin dotted curves show the two within-species fits; solid/dashed curves are equal-weight species summaries",
    fontsize=9.2,
)
_save_figure(figure, "representative_gene_spatial_phenotypes.png")
plt.show()

# Retain every selected observation in a separate, full-resolution supplementary atlas.  It is saved
# but deliberately not displayed inline: a multi-panel figure is unsuitable as a main manuscript panel.
full_figure, full_axes = plt.subplots(
    len(REPRESENTATIVE_CLASSES), GENES_PER_CLASS,
    figsize=(8.0, 1.35 * len(REPRESENTATIVE_CLASSES)), sharex=True, layout="constrained",
)
full_axes = np.atleast_2d(full_axes)
for row_index, label in enumerate(REPRESENTATIVE_CLASSES):
    block = representative_table[representative_table["spatial_phenotype"].eq(label)]
    for column_index, axis in enumerate(full_axes[row_index]):
        if column_index >= len(block):
            axis.set_visible(False)
            continue
        record = block.iloc[column_index]
        index = int(np.flatnonzero(gene_names == record["gene"])[0])
        for specimen in mouse_samples:
            axis.plot(grid_unit, specimen_curves[specimen][index], color=SPECIES_COLORS["mouse"],
                      lw=0.45, alpha=0.25, ls=":")
        for specimen in human_samples:
            axis.plot(grid_unit, specimen_curves[specimen][index], color=SPECIES_COLORS["human"],
                      lw=0.45, alpha=0.25, ls=":")
        axis.plot(grid_unit, balanced_mouse[index], color=SPECIES_COLORS["mouse"], lw=1.1)
        axis.plot(grid_unit, balanced_human[index], color=SPECIES_COLORS["human"], lw=1.1, ls="--")
        axis.set_title(record["gene"], fontsize=6.3, fontweight="bold")
        if column_index == 0:
            axis.set_ylabel(textwrap.fill(label, 18), fontsize=6.5)
        if row_index == len(REPRESENTATIVE_CLASSES) - 1:
            axis.set_xlabel("PT position", fontsize=6.5)
        axis.tick_params(labelsize=5.5)
full_figure.suptitle("Supplementary exemplar atlas: three high-evidence genes per phenotype", fontsize=8.5)
_save_figure(full_figure, "supplementary_representative_gene_atlas.png")
plt.close(full_figure)

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
human_shape = _row_zscore(balanced_human[heatmap_index])

figure = plt.figure(figsize=(7.6, max(6.5, 0.175 * len(heatmap_index) + 0.8)), layout="constrained")
grid_spec = figure.add_gridspec(1, 3, width_ratios=[1, 1, 0.045])
axes = [figure.add_subplot(grid_spec[0, 0]), figure.add_subplot(grid_spec[0, 1])]
colour_axis = figure.add_subplot(grid_spec[0, 2])

for axis, matrix, title in zip(axes, [mouse_shape, human_shape], ["Healthy mouse", "Human"]):
    image = axis.imshow(
        matrix, aspect="auto", interpolation="nearest", cmap="RdBu_r", vmin=-2.5, vmax=2.5,
        extent=[0, 1, len(heatmap_index) - 0.5, -0.5],
    )
    _style_pt_axis(axis)
    axis.set_title(title, fontsize=8.5, fontweight="bold")

    # Class boundaries.
    boundaries = np.cumsum([len(heatmap_table[heatmap_table["spatial_phenotype"].eq(label)])
                            for label in heatmap_classes])[:-1]
    for boundary in boundaries:
        axis.axhline(boundary - 0.5, color="black", lw=0.6, alpha=0.5)

axes[0].set_yticks(np.arange(len(heatmap_index)))
axes[0].set_yticklabels(
    [f"{gene}   [{label}]" for gene, label in
     zip(heatmap_table["gene"], heatmap_table["spatial_phenotype"])],
    fontsize=5.8,
)
axes[0].set_ylabel("Genes, grouped by spatial phenotype", fontsize=7)
axes[1].tick_params(labelleft=False)

colour_bar = figure.colorbar(image, cax=colour_axis, label="Within-species curve z-score")
colour_bar.ax.tick_params(labelsize=6.5)
colour_axis.set_label("<colorbar>")

figure.suptitle(
    "Spatial organisation only: the same gene ordering in mouse and human\n"
    "Each row is z-scored within species, so abundance differences cannot drive the pattern",
    fontsize=9.2,
)
_panel_label(axes[0], "e")
_save_figure(figure, "shape_only_gene_heatmap.png")
plt.show()

_save_table(heatmap_table, "shape_only_heatmap_gene_selection.csv")


# %% [markdown]
# ## 4 - Which functional programs sit in which spatial phenotype
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
# **What counts as evidence - and what it licenses.** A fold enrichment with a family-wide FDR below
# `enrichment_fdr` establishes that a pathway's genes are **over-represented among genes classified as**
# that phenotype. That is not the same claim as "this pathway is a conserved/rewired programme", so each
# association also carries a **member-consistency** check: the share of the pathway's tested members
# that hold the phenotype (`member_fraction_in_phenotype`, at least `pathway_consistency_min_fraction`),
# and whether those members behave the way the phenotype claims. That behavioural test is
# **phenotype-specific** (`member_median_shape_corr`, `member_median_amplitude_log2`,
# `member_median_shift`): a negative median shape correlation is the expected signature of an inversion
# or a complex rewiring, whereas a conserved-zonation pathway whose members are not shape-conserved is
# not conserved. A pathway that fails consistency is labelled "enrichment only", which is a real result
# about its member genes and not a claim about the pathway as a whole. Leading genes are named in every row so the
# signal can be checked against a single strong gene. Redundant databases overlap heavily, so the companion heatmap collapses redundancy groups with 03's
# `summarize_pathway_redundancy` - but only near-duplicate sets (overlap >= `pathway_redundancy_overlap`,
# 0.85): a lenient threshold would reduce 1,511 pathways to a few dozen groups and hide the associations
# this section exists to report. The primary figure shows **every** significant association directly, as
# a dot plot, with the member fraction annotated on each dot.
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


def _member_consistency(label, medians, config):
    """Does this pathway's membership of a phenotype describe the pathway's members, for that phenotype?

    Phenotype-specific on purpose. A single "median shape correlation >= 0" rule would mark every
    inversion and complex-rewiring pathway inconsistent, because a negative or low median shape
    correlation is exactly what those phenotypes *are*; conversely a conserved-zonation pathway that
    is not shape-conserved is not conserved. So consistency means "the members behave the way this
    phenotype claims", and what that requires differs per phenotype.
    """
    fraction = medians["fraction"]
    if not np.isfinite(fraction) or fraction < config["pathway_consistency_min_fraction"]:
        return "enrichment only (too few members hold the phenotype)"

    shape = medians["shape_corr"]
    amplitude = medians["amplitude_log2"]

    if label == "conserved zonation":
        consistent = shape >= config["pathway_consistency_min_median_shape_corr"]
    elif label == "weaker zonation in human":
        consistent = np.isfinite(amplitude) and amplitude <= 0
    elif label == "stronger zonation in human":
        consistent = np.isfinite(amplitude) and amplitude >= 0
    elif label == "gradient inversion":
        consistent = np.isfinite(shape) and shape <= 0
    elif label == "complex shape rewiring":
        consistent = np.isfinite(shape) and shape < config["conserved_corr"]
    elif label in {"mouse-zonated / human-flat", "human-zonated / mouse-flat"}:
        consistent = np.isfinite(amplitude) and abs(amplitude) >= config["amplitude_change_major"]
    else:
        return "not assessed"

    return "members consistent" if consistent else "enrichment only (members not consistent)"


def _phenotype_enrichment(phenotype_labels, universe_members_index, pathway_frame, config,
                          member_evidence=None):
    """Upper-tail hypergeometric over-representation of each phenotype among each pathway's members.

    `phenotype_labels` is the per-gene phenotype array over every tested gene; `universe_members_index`
    are the rows that form the background (the discovery-eligible genes); `pathway_frame` must carry
    `universe_members` (row indices) plus `library`/`pathway`. `member_evidence` is an optional frame of
    per-gene evidence (`shape_corr`, `amplitude_log2`, `abs_shift`) used for the coherence columns:
    enrichment is a statement about member *membership*, coherence is a statement about member
    *behaviour*. Without an evidence frame the coherence columns are reported as not assessed. Only over-representation is tested, so a
    pathway that is *short* of a phenotype returns p ~ 1 rather than a small p - a shortfall is not
    evidence of enrichment, and BH in the caller is applied across the complete pathway x phenotype
    family. Kept as a function so the notebook's self-check can run it on a known-answer case.
    """
    universe_mask = np.zeros(phenotype_labels.size, dtype=bool)
    universe_mask[universe_members_index] = True
    universe_size = int(universe_members_index.size)
    universe_counts = pd.Series(phenotype_labels[universe_members_index]).value_counts()

    rows = []
    for pathway_row in pathway_frame.itertuples():
        members = np.asarray(
            [index for index in pathway_row.universe_members if universe_mask[index]], dtype=int
        )
        if members.size < config["pathway_min_members"]:
            continue

        member_phenotypes = phenotype_labels[members]

        for label in PHENOTYPE_CLASSES:
            n_with_phenotype = int(universe_counts.get(label, 0))
            if n_with_phenotype == 0:
                continue

            observed = int((member_phenotypes == label).sum())
            expected = n_with_phenotype * members.size / universe_size
            overlap = members[member_phenotypes == label]
            leading = overlap[np.argsort(-discovery_score_by_row[overlap])][:10]

            # Coherence: is this phenotype's share of the pathway large enough to describe the
            # pathway, and do the members behave that way overall? Enrichment alone is a statement
            # about membership; coherence is a statement about behaviour.
            fraction_in_phenotype = observed / members.size if members.size else np.nan
            if member_evidence is None:
                member_shape_corr = member_amplitude = member_shift = np.nan
                member_abs_amplitude = member_abs_shift = np.nan
                member_consistency = "not assessed"
            else:
                member_shape_corr = float(np.nanmedian(
                    member_evidence["shape_corr"].to_numpy()[members]))
                member_amplitude = float(np.nanmedian(
                    member_evidence["amplitude_log2"].to_numpy()[members]))
                member_abs_amplitude = float(np.nanmedian(
                    np.abs(member_evidence["amplitude_log2"].to_numpy()[members])))
                member_shift = float(np.nanmedian(
                    member_evidence["shift"].to_numpy()[members]))
                member_abs_shift = float(np.nanmedian(
                    member_evidence["abs_shift"].to_numpy()[members]))
                member_consistency = _member_consistency(label, {
                    "fraction": fraction_in_phenotype,
                    "shape_corr": member_shape_corr,
                    "amplitude_log2": member_amplitude,
                    "shift": member_shift,
                }, config)

            rows.append({
                "library": pathway_row.library,
                "pathway": pathway_row.pathway,
                "spatial_phenotype": label,
                "n_pathway_members_in_universe": int(members.size),
                "n_universe_with_phenotype": n_with_phenotype,
                "observed": observed,
                "expected": expected,
                "fold_enrichment": (observed / expected) if expected > 0 else np.nan,
                "p_value": float(hypergeom.sf(observed - 1, universe_size, members.size,
                                              n_with_phenotype)),
                "leading_genes": ";".join(gene_names[leading]),
                "member_fraction_in_phenotype": fraction_in_phenotype,
                "member_median_shape_corr": member_shape_corr,
                "member_median_amplitude_log2": member_amplitude,
                "member_median_abs_amplitude_log2": member_abs_amplitude,
                "member_median_shift": member_shift,
                "member_median_abs_shift": member_abs_shift,
                "member_consistency": member_consistency,
            })

    # Declared columns matter: with no surviving pair the frame must still be addressable by name
    # rather than collapsing to a bare DataFrame with no columns at all.
    return pd.DataFrame(rows, columns=[
        "library", "pathway", "spatial_phenotype", "n_pathway_members_in_universe",
        "n_universe_with_phenotype", "observed", "expected", "fold_enrichment", "p_value",
        "leading_genes", "member_fraction_in_phenotype", "member_median_shape_corr",
        "member_median_amplitude_log2", "member_median_abs_amplitude_log2", "member_median_shift",
        "member_median_abs_shift", "member_consistency",
    ])


# Per-gene evidence for the coherence columns: how the members behave, not only how many there are.
member_evidence = pd.DataFrame({
    "shape_corr": gene_metrics["shape_corr"].to_numpy(),
    "amplitude_log2": gene_metrics["amplitude_log2_ratio_human_over_mouse"].to_numpy(),
    "shift": gene_metrics["best_shift_human_minus_mouse"].to_numpy(),
    "abs_shift": np.abs(gene_metrics["best_shift_human_minus_mouse"].to_numpy()),
})

enrichment_rows = _phenotype_enrichment(
    phenotype_by_row, universe_index, pathway_tested, CONFIG, member_evidence=member_evidence
)

pathway_enrichment = pd.DataFrame(enrichment_rows)
pathway_enrichment = _numeric(pathway_enrichment, [
    "n_pathway_members_in_universe", "n_universe_with_phenotype", "observed", "expected",
    "fold_enrichment", "p_value", "member_fraction_in_phenotype", "member_median_shape_corr",
    "member_median_amplitude_log2", "member_median_abs_amplitude_log2", "member_median_shift",
    "member_median_abs_shift",
])
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
# Purpose: the significant pathway x phenotype associations, and the redundancy-collapsed view.

redundancy_pairs, redundancy_groups = summarize_pathway_redundancy(
    membership.assign(genes_present=membership["member_list"]),
    gene_column="genes_present",
    overlap_threshold=CONFIG["pathway_redundancy_overlap"],
)

group_of_pathway = {
    (row.library, row.pathway): int(row.redundancy_group) for row in redundancy_groups.itertuples()
}
group_size_of_pathway = {
    (row.library, row.pathway): int(row.group_size) for row in redundancy_groups.itertuples()
}

_report("pathway overlap pairs at >= the redundancy threshold", len(redundancy_pairs))
_report("redundancy groups", int(redundancy_groups["redundancy_group"].nunique()),
        len(redundancy_groups),
        why=f"only near-duplicate sets are collapsed (overlap >= {CONFIG['pathway_redundancy_overlap']})")

significant_associations = pathway_enrichment[
    pathway_enrichment["fdr"] < CONFIG["enrichment_fdr"]
].sort_values("fdr").reset_index(drop=True)

significant_associations = significant_associations.assign(
    redundancy_group=[group_of_pathway.get((library, pathway), np.nan)
                      for library, pathway in zip(significant_associations["library"],
                                                  significant_associations["pathway"])],
    redundancy_group_size=[group_size_of_pathway.get((library, pathway), np.nan)
                           for library, pathway in zip(significant_associations["library"],
                                                       significant_associations["pathway"])],
    pathway_label=[
        f"{pathway}  [{library.replace('_2022', '').replace('MSigDB_Hallmark_2020', 'Hallmark').replace('KEGG_2019_Mouse', 'KEGG')}]"
        for library, pathway in zip(significant_associations["library"],
                                    significant_associations["pathway"])
    ],
)

_report("significant pathway x phenotype associations",
        len(significant_associations), len(pathway_enrichment),
        why=f"family-wide FDR < {CONFIG['enrichment_fdr']}; shown directly, not collapsed away")

_save_table(significant_associations, "pathway_significant_associations.csv")

# Bound unconditionally: section 4.2's strongest-findings panels read these whether or not this run
# produced any significant association.
heatmap_frame = pd.DataFrame()

if significant_associations.empty:
    print(f"No pathway x phenotype pair reached FDR < {CONFIG['enrichment_fdr']}: "
          "both figures are skipped rather than drawn empty. The full table is still saved.")
else:
    # ----------------------------------------------------------------------------------
    # Main figure: every significant association, as a dot plot.
    #
    # Cell value = the member fraction, which is the honest way to read an enrichment: how much of
    # the pathway actually holds the phenotype. Dot size carries fold enrichment, colour carries
    # significance, and the row label carries whether the pathway is *coherent* for that phenotype.
    # ----------------------------------------------------------------------------------
    phenotype_columns = [label for label in PHENOTYPE_CLASSES
                         if label in set(significant_associations["spatial_phenotype"])]
    row_order = (significant_associations.groupby("pathway_label")["fdr"].min()
                 .sort_values().index.tolist())

    figure, axis = plt.subplots(
        1, 1, figsize=(7.3, max(3.25, 0.52 * len(row_order) + 1.55)), layout="constrained"
    )

    significance_scale = float(significant_associations["minus_log10_fdr"].max())
    consistency_x = len(phenotype_columns)
    for row_index, pathway_label in enumerate(row_order):
        block = significant_associations[significant_associations["pathway_label"].eq(pathway_label)]
        for entry in block.itertuples():
            x_position = phenotype_columns.index(entry.spatial_phenotype)
            axis.scatter(
                x_position, row_index,
                s=28 + 18 * float(np.clip(entry.fold_enrichment, 0, 10)),
                c=[entry.minus_log10_fdr], cmap="viridis",
                vmin=0, vmax=max(significance_scale, 1e-9),
                edgecolors="#202020", linewidths=0.55, zorder=3,
            )
            axis.annotate(
                f"{entry.member_fraction_in_phenotype:.0%}", (x_position, row_index),
                xytext=(0, -11), textcoords="offset points", ha="center", fontsize=5.9,
            )
        consistent = all(entry.member_consistency == "members consistent" for entry in block.itertuples())
        axis.text(consistency_x, row_index, "yes" if consistent else "no", ha="center", va="center",
                  fontsize=6.8, color=PHENOTYPE_COLORS["conserved zonation"] if consistent else "0.40",
                  fontweight=("bold" if consistent else "normal"))

    axis.set_xticks(np.arange(len(phenotype_columns) + 1))
    axis.set_xticklabels([textwrap.fill(label, 17) for label in phenotype_columns] + ["Members\nconsistent?"],
                         fontsize=7)
    axis.set_yticks(np.arange(len(row_order)))
    axis.set_yticklabels([textwrap.fill(label, 39) for label in row_order], fontsize=6.6)
    axis.set_xlim(-0.55, consistency_x + 0.55)
    axis.set_ylim(-0.6, len(row_order) - 0.4)
    axis.invert_yaxis()
    axis.set_xlabel("Spatial phenotype enriched among pathway members")
    axis.grid(axis="x", color="0.92", lw=0.6, zorder=0)
    axis.set_title(
        "Pathway-by-phenotype enrichment evidence", loc="left", fontweight="bold", fontsize=9,
    )
    colour_bar = figure.colorbar(axis.collections[0], ax=axis, pad=0.015, fraction=0.045)
    colour_bar.set_label("−log₁₀(FDR)", fontsize=6.8)
    colour_bar.ax.tick_params(labelsize=6.1)
    _panel_label(axis, "a")
    figure.suptitle(
        f"Six pathway–phenotype associations pass family-wide FDR < {CONFIG['enrichment_fdr']}\n"
        "Dot area = fold enrichment; colour = −log₁₀(FDR); label = member fraction. “Members consistent” is a stricter behavioural check, not a second significance test.",
        fontsize=8.4,
    )
    _save_figure(figure, "pathway_spatial_phenotype_dotplot.png")
    plt.show()

    # ----------------------------------------------------------------------------------
    # Companion figure: the same associations collapsed to redundancy groups, so a pathway
    # family appears once.
    # ----------------------------------------------------------------------------------
    representatives = []
    seen_groups = set()
    for row in significant_associations.itertuples():
        group = row.redundancy_group
        if np.isfinite(group) and int(group) in seen_groups:
            continue
        if np.isfinite(group):
            seen_groups.add(int(group))
        representatives.append(row)

    heatmap_frame = pd.DataFrame([{
        "library": row.library, "pathway": row.pathway,
        "spatial_phenotype": row.spatial_phenotype,
        "signed_strength": row.signed_strength, "fdr": row.fdr,
        "fold_enrichment": row.fold_enrichment,
        "member_fraction_in_phenotype": row.member_fraction_in_phenotype,
        "member_consistency": row.member_consistency,
        "redundancy_group": row.redundancy_group,
        "group_size": row.redundancy_group_size,
        "leading_genes": row.leading_genes,
    } for row in representatives])

    heatmap_matrix = (
        heatmap_frame
        .pivot_table(index=["library", "pathway"], columns="spatial_phenotype",
                     values="signed_strength", aggfunc="max")
        .reindex(columns=phenotype_columns)
    )
    heatmap_matrix = heatmap_matrix.reindex(
        heatmap_matrix.abs().max(axis=1).sort_values(ascending=False).index
    )
    scale = float(np.nanmax(np.abs(heatmap_matrix.to_numpy())))

    figure, axis = plt.subplots(
        1, 1, figsize=(1.15 * heatmap_matrix.shape[1] + 6.5,
                       max(3.0, 0.48 * heatmap_matrix.shape[0] + 1.2)), layout="constrained"
    )
    image = axis.imshow(heatmap_matrix.to_numpy(), aspect="auto", interpolation="nearest",
                        cmap="PuOr_r", vmin=-scale, vmax=scale)
    axis.set_xticks(np.arange(heatmap_matrix.shape[1]))
    axis.set_xticklabels([textwrap.fill(str(label), 18) for label in heatmap_matrix.columns],
                         rotation=45, rotation_mode="anchor", ha="right", fontsize=8)
    axis.set_yticks(np.arange(heatmap_matrix.shape[0]))
    axis.set_yticklabels([textwrap.fill(f"{pathway}  [{library}]", 46)
                          for library, pathway in heatmap_matrix.index], fontsize=7.5)
    for row_index in range(heatmap_matrix.shape[0]):
        for column_index in range(heatmap_matrix.shape[1]):
            value = heatmap_matrix.to_numpy()[row_index, column_index]
            if np.isfinite(value):
                axis.text(column_index, row_index, f"{value:.1f}", ha="center", va="center",
                          fontsize=6.5)
    figure.colorbar(image, ax=axis, shrink=0.7,
                    label="Signed -log10(FDR): + enriched, - depleted\n(one term per redundancy group)")
    axis.set_title(
        "The same associations, collapsed to redundancy groups\n"
        "near-duplicate sets only; the dot plot above is the primary view",
        loc="left", fontsize=9.5,
    )
    _save_figure(figure, "pathway_spatial_phenotype_heatmap.png")
    plt.show()

    _save_table(heatmap_frame, "pathway_spatial_phenotype_heatmap_terms.csv")


# %%
# Purpose: known-answer self-check of the over-representation test.

def _enrichment_selfcheck():
    """Run the enrichment test on a synthetic phenotype assignment with hand-computed answers.

    The enrichment step decides which functional programmes the notebook reports, and its arithmetic
    (background size, eligible members, upper tail, family-wide correction) is easy to get subtly
    wrong while still producing plausible-looking numbers. This checks it against values computed by
    hand: 60 background genes, a pathway holding 20 of them, 15 of which share a phenotype gives
    expected = 5, fold = 3, p = hypergeom.sf(14, 60, 20, 15).
    """
    labels = np.array(
        ["conserved zonation"] * 15
        + ["gradient inversion"] * 5
        + ["weak / uncertain zonation"] * 30
        + ["excluded"] * int(len(gene_names) - 50)
    )
    background = np.arange(60)
    pathway_frame = pd.DataFrame({
        "library": ["synthetic", "synthetic"],
        "pathway": ["enriched_case", "shortfall_case"],
        "universe_members": [list(range(20)), list(range(20, 40))],
    })

    table = _phenotype_enrichment(labels, background, pathway_frame, CONFIG)

    def _pair(pathway, phenotype):
        block = table[table["pathway"].eq(pathway) & table["spatial_phenotype"].eq(phenotype)]
        return block.iloc[0] if len(block) else None

    enriched = _pair("enriched_case", "conserved zonation")
    shortfall = _pair("shortfall_case", "conserved zonation")

    # Derive which pathway holds every gene of the inverted phenotype from the construction itself,
    # instead of assuming it: the answer follows from the indices above.
    inverted_rows = [row for row in
                     (pathway_frame.iloc[0], pathway_frame.iloc[1])
                     if set(np.flatnonzero(labels == "gradient inversion"))
                     <= set(row["universe_members"])]
    inverted = _pair(inverted_rows[0]["pathway"], "gradient inversion") if inverted_rows else None

    checks = [
        ("enriched pair: background size is the discovery-eligible universe",
         enriched is not None and int(enriched["n_universe_with_phenotype"]) == 15),
        ("enriched pair: eligible members are counted inside the universe",
         enriched is not None and int(enriched["n_pathway_members_in_universe"]) == 20),
        ("enriched pair: observed overlap is counted exactly",
         enriched is not None and int(enriched["observed"]) == 15),
        ("enriched pair: expected is n*K/N = 5",
         enriched is not None and abs(float(enriched["expected"]) - 5.0) < 1e-9),
        ("enriched pair: fold enrichment is 3",
         enriched is not None and abs(float(enriched["fold_enrichment"]) - 3.0) < 1e-9),
        ("enriched pair: p is the hypergeometric upper tail",
         enriched is not None and abs(float(enriched["p_value"])
                                      - float(hypergeom.sf(14, 60, 20, 15))) < 1e-12),
        ("shortfall is not reported as enrichment (p ~ 1)",
         shortfall is None or float(shortfall["p_value"]) > 0.99),
        ("a pathway holding every member of a phenotype returns a small p",
         inverted is not None and float(inverted["p_value"]) < 0.01),
        ("only phenotypes present in the background are tested",
         "excluded" not in set(table["spatial_phenotype"])),
    ]

    return pd.DataFrame([{"check": name, "passed": bool(passed)} for name, passed in checks])


enrichment_selfcheck = _enrichment_selfcheck()
display(enrichment_selfcheck)

if not enrichment_selfcheck["passed"].all():
    failed = enrichment_selfcheck.loc[~enrichment_selfcheck["passed"], "check"].tolist()
    raise AssertionError(f"synthetic enrichment self-check failed for: {failed}")

print(f"All {len(enrichment_selfcheck)} synthetic enrichment self-checks passed.")


# %% [markdown]
# ### 4.2 - Inspecting whole pathways along the coordinate
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

def _ordered_member_blocks(members, limit):
    """Limit members, order them by human positional centroid, and z-score each species' block.

    Shared by the candidate-theme panels and the strongest-findings panels so the two figures cannot
    drift apart in ordering, scaling or colour limits.
    """
    selected = np.asarray(members[:limit], dtype=int)
    centroids = _row_weighted_centroid(balanced_human[selected], grid_unit)
    ordered = selected[np.argsort(np.nan_to_num(centroids, nan=1.0))]
    return ordered, _row_zscore(balanced_mouse[ordered]), _row_zscore(balanced_human[ordered])


# One pathway per page: no raster strip is silently compressed to make a seven-theme dashboard fit.
def _save_pathway_member_pair(*, members, title, subtitle, filename, limit):
    """Save a legible mouse/human trajectory plate for one selected pathway term."""
    ordered, mouse_block, human_block = _ordered_member_blocks(members, limit)
    figure = plt.figure(figsize=(7.45, max(2.85, 0.18 * len(ordered) + 1.8)), layout="constrained")
    grid_spec = figure.add_gridspec(1, 3, width_ratios=[1, 1, 0.045])
    axes = [figure.add_subplot(grid_spec[0, 0]), figure.add_subplot(grid_spec[0, 1])]
    colour_axis = figure.add_subplot(grid_spec[0, 2])
    for axis, block, species in zip(axes, [mouse_block, human_block], ["Healthy mouse", "Human"]):
        image = axis.imshow(block, aspect="auto", interpolation="nearest", cmap="RdBu_r", vmin=-2.5, vmax=2.5,
                            extent=[0, 1, len(ordered) - 0.5, -0.5])
        _style_pt_axis(axis)
        axis.set_title(species, loc="left", fontsize=8.5, fontweight="bold")
        axis.tick_params(labelsize=5.8)
    axes[0].set_yticks(np.arange(len(ordered)))
    axes[0].set_yticklabels([f"{gene_names[index]}  [{phenotype_by_row[index]}]" for index in ordered],
                             fontsize=5.6)
    axes[0].set_ylabel("Genes (ordered by human positional centroid)", fontsize=6.8)
    axes[1].tick_params(labelleft=False)
    colour_bar = figure.colorbar(image, cax=colour_axis, label="Within-species\ncurve z-score")
    colour_bar.ax.tick_params(labelsize=5.8)
    colour_axis.set_label("<colorbar>")
    figure.suptitle(f"{title}\n{subtitle}; {len(ordered)} of {len(members)} tested members displayed",
                    fontsize=8.7)
    _save_figure(figure, filename)
    plt.close(figure)


# A compact index says which candidate themes have a significant association; full trajectory plates
# are written one term per page below.  Candidate-theme inspection is not elevated to enrichment evidence.
theme_display = theme_summary.copy()
theme_display["minus_log10_best_fdr"] = -np.log10(theme_display["best_fdr"])
theme_display = theme_display.sort_values("minus_log10_best_fdr", na_position="last")
figure, axis = plt.subplots(1, 1, figsize=(7.0, 3.4), layout="constrained")
positions = np.arange(len(theme_display))
colours = [PHENOTYPE_COLORS.get(label, "0.65") if np.isfinite(fdr) and fdr < CONFIG["enrichment_fdr"] else "0.78"
           for label, fdr in zip(theme_display["best_phenotype"], theme_display["best_fdr"])]
axis.barh(positions, theme_display["minus_log10_best_fdr"].fillna(0), color=colours)
axis.axvline(-np.log10(CONFIG["enrichment_fdr"]), color="0.35", lw=0.8, ls="--")
axis.set_yticks(positions)
axis.set_yticklabels(theme_display["theme"], fontsize=7)
axis.invert_yaxis()
axis.set_xlabel("−log₁₀(best family-wide FDR)")
axis.set_title("Candidate pathway themes: evidence screen", loc="left", fontweight="bold")
for position, row in enumerate(theme_display.itertuples()):
    label = row.best_phenotype if np.isfinite(row.best_fdr) else "no tested association"
    axis.text(max(row.minus_log10_best_fdr if np.isfinite(row.minus_log10_best_fdr) else 0, 0), position,
              f"  {label}", va="center", fontsize=6.1, color="0.30")
axis.text(0.0, 1.03, "Dashed line: FDR = 0.05. Colours denote the phenotype only for associations that pass it.",
          transform=axis.transAxes, fontsize=6.2, va="bottom", color="0.30")
_save_figure(figure, "representative_pathway_spatial_heatmaps.png")
plt.show()

plotted_themes = [theme for theme in PATHWAY_THEMES if theme in theme_selections]
if not plotted_themes:
    print("No candidate theme matched a pathway in the membership file: detailed plates are skipped.")
else:
    for theme in plotted_themes:
        selection = theme_selections[theme]
        best_fdr = theme_summary.loc[theme_summary["theme"].eq(theme), "best_fdr"].iloc[0]
        best_label = theme_summary.loc[theme_summary["theme"].eq(theme), "best_phenotype"].iloc[0]
        subtitle = (f"{selection['library']}: {selection['pathway']}; strongest enrichment: "
                    f"{best_label or 'none'} (FDR={best_fdr:.3g})" if np.isfinite(best_fdr) else
                    f"{selection['library']}: {selection['pathway']}; no significant phenotype enrichment")
        _save_pathway_member_pair(
            members=selection["members"], title=f"Supplementary candidate pathway: {theme}", subtitle=subtitle,
            filename=f"supplementary_candidate_pathway_{re.sub(r'[^a-z0-9]+', '_', theme.lower()).strip('_')}.png",
            limit=MAX_MEMBERS_PER_THEME,
        )


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
        # The paper-facing overview is a readable evidence ranking; detailed member-level evidence is
        # one supplementary plate per pathway so labels and trajectories retain their physical scale.
        ranked = strongest_pathway_findings.sort_values("fdr").reset_index(drop=True)
        figure, axis = plt.subplots(1, 1, figsize=(7.2, 3.35), layout="constrained")
        positions = np.arange(len(ranked))
        values = -np.log10(ranked["fdr"])
        axis.barh(positions, values,
                  color=[PHENOTYPE_COLORS.get(label, "0.55") for label in ranked["spatial_phenotype"]])
        axis.axvline(-np.log10(CONFIG["enrichment_fdr"]), color="0.35", lw=0.8, ls="--")
        axis.set_yticks(positions)
        axis.set_yticklabels([textwrap.fill(row.pathway, 37) for row in ranked.itertuples()], fontsize=6.7)
        axis.invert_yaxis()
        axis.set_xlabel("−log₁₀(family-wide FDR)")
        axis.set_title("Strongest non-redundant pathway associations", loc="left", fontweight="bold")
        for position, row in enumerate(ranked.itertuples()):
            axis.text(-0.02, position, row.spatial_phenotype, transform=axis.get_yaxis_transform(),
                      ha="right", va="center", fontsize=5.8, color="0.30")
        axis.text(0.0, 1.03, "Colour: associated spatial phenotype. Dashed line: FDR = 0.05.",
                  transform=axis.transAxes, fontsize=6.2, va="bottom", color="0.30")
        _save_figure(figure, "strongest_pathway_spatial_heatmaps.png")
        plt.show()

        for rank, entry in enumerate(strongest_selections, start=1):
            subtitle = (f"{entry['library']}: {entry['pathway']} → {entry['spatial_phenotype']}; "
                        f"fold={entry['fold_enrichment']:.2f}, FDR={entry['fdr']:.3g}")
            _save_pathway_member_pair(
                members=entry["members"], title=f"Supplementary pathway finding {rank}", subtitle=subtitle,
                filename=f"supplementary_top_pathway_{rank:02d}.png", limit=MAX_MEMBERS_PER_STRONGEST,
            )


# %% [markdown]
# ## 5 - What the continuous analysis adds beyond the whole-PT comparison
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
# divergence - never from level - and each term counts **only where its evidence is interpretable**:
# amplitude difference needs a reproducibly patterned species, displacement needs every displacement
# criterion, shape divergence needs both species patterned and reproducible. A gene with nothing
# interpretable scores exactly 0 instead of accumulating division-amplified ratios from flat curves,
# boundary shifts and noisy correlations, so the weak/uncertain cloud sits at the origin rather than
# dominating the continuous-only quadrant.
#
# One distinction to keep in mind when reading the quadrants: a gene can hold a non-zero score while
# still being classified `weak / uncertain zonation` - that happens when one component *is*
# interpretable (say a real amplitude difference) but the slice or specimen agreement needed for a
# species-specific class is not there. The class reflects that disagreement; the score only
# reports which components are measurable.
#
# **How the quadrants are drawn - and why not by a median split.** The conventional axis is continuous
# and well behaved, so its median is a fair boundary. The spatial axis is not: the gated score is exactly
# 0 for every gene with no interpretable component (that is the point of the gate), so its median in this
# cohort is 0 and a "score >= median" rule would call *every* gene high-spatial and leave two quadrants
# empty. Strong spatial evidence is therefore treated as a **categorical** state - a non-weak phenotype
# that kept its broad class in at least `spatial_strong_robustness` of the section-3.2 variants - and the
# continuous score is used only to rank genes inside that group. The upper-left quadrant then answers the
# question the section exists for: reproducible spatial remodelling despite a small whole-PT level
# difference.
#
# **One asymmetry the score no longer has.** Amplitude evidence used to require reproducible
# correlation in *both* species, which is undefined for a species that has lost its zonation - so
# "mouse-zonated / human-flat" genes scored zero while their mirror image scored high. Amplitude
# evidence is now symmetric: both species reproducibly patterned, **or** one side patterned and the
# other flat in every specimen/slice (the same evidence the species-specific classes already require).
#
# **What the conventional comparator is, and what a null result means.** The conventional axis here is
# not the median absolute level effect of a pathway's members - that is a descriptive member statistic,
# not how a whole-PT study tests a pathway. The comparator is a **whole-PT preranked GSEA on the
# level-effect ranking**, so the spatial result is compared against the analysis a conventional paper
# would actually run. `spatially_only_pathways.csv` lists what is spatially significant while not
# reaching significance there - and the name matters, because "spatially only" does **not** mean
# "conventionally weak": these pathways can carry large absolute member level effects. What the
# whole-PT GSEA lacks is a *coherent direction*, since reorganising a pathway along PT moves its
# genes different ways and their whole-PT shifts cancel. If that list comes out empty or nearly so, the honest reading is that conventional
# analysis already reports these programmes as different at the level of whole-PT expression, and the
# continuous framework's contribution is *how* they differ - zonation gain or loss, displacement,
# gradient reversal, trajectory rewiring - rather than *that* they differ.
#
# **Descriptive, by construction.** No p-values are used on the conventional axis: with one human donor
# a species-wide test would not be calibrated. The quadrant counts and pathway table are effect-size
# statements that hold for this cohort.
#
# **What counts as evidence, at pathway level.** A pathway is *spatially only* when its member genes are
# enriched in a spatial phenotype at family-wide FDR < `enrichment_fdr` **and** the conventional
# whole-PT GSEA does not reach that threshold for it. Both results are reported in the table, so the
# selection can be checked rather than trusted.
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

# Strong spatial evidence is categorical, and deliberately so. The gated score is exactly 0 for every
# gene with no interpretable component, which in this cohort is most genes - so its median is 0 and a
# median split would make "high spatial" mean "score >= 0", i.e. every gene, leaving two quadrants
# empty. The categorical rule below asks for a non-weak phenotype that survived the sensitivity panel;
# the continuous score is then used only to rank genes inside that group.
quadrant_frame["strong_spatial_evidence"] = (
    quadrant_frame["spatial_phenotype"].ne("weak / uncertain zonation")
    & quadrant_frame["spatial_phenotype"].ne("excluded")
    & quadrant_frame["robustness_score"].ge(CONFIG["spatial_strong_robustness"])
)

_report("genes with strong, reproducible spatial evidence",
        int(quadrant_frame["strong_spatial_evidence"].sum()), len(quadrant_frame),
        why=f"non-weak phenotype kept in >= {CONFIG['spatial_strong_robustness']:.0%} of the "
            "sensitivity variants")

quadrant_labels = {
    "high conventional / strong spatial":
        quadrant_frame["conventional_abs"].ge(conventional_median)
        & quadrant_frame["strong_spatial_evidence"],
    "high conventional / weak spatial":
        quadrant_frame["conventional_abs"].ge(conventional_median)
        & ~quadrant_frame["strong_spatial_evidence"],
    "low conventional / strong spatial (continuous-only)":
        quadrant_frame["conventional_abs"].lt(conventional_median)
        & quadrant_frame["strong_spatial_evidence"],
    "low conventional / weak spatial":
        quadrant_frame["conventional_abs"].lt(conventional_median)
        & ~quadrant_frame["strong_spatial_evidence"],
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

figure, axes = plt.subplots(1, 2, figsize=(9.4, 3.9), width_ratios=[1.35, 0.78], layout="constrained")
axis, count_axis = axes

weak = quadrant_frame[~quadrant_frame["strong_spatial_evidence"]]
axis.scatter(weak["conventional_abs"], np.log10(1 + weak["spatial_discovery_score"]), s=4, alpha=0.12,
             linewidths=0, color="0.55", rasterized=True,
             label=f"no strong spatial evidence (n={len(weak):,})")

for label in PHENOTYPE_CLASSES:
    block = quadrant_frame[
        quadrant_frame["spatial_phenotype"].eq(label) & quadrant_frame["strong_spatial_evidence"]
    ]
    if not len(block):
        continue
    axis.scatter(
        block["conventional_abs"], np.log10(1 + block["spatial_discovery_score"]),
        s=8, alpha=0.72, linewidths=0, color=PHENOTYPE_COLORS[label], rasterized=True,
        label=f"{label} (n={len(block):,})",
    )

axis.axvline(conventional_median, color="0.25", lw=0.8, ls="--")
axis.set_xlabel("Conventional whole-PT effect: |level effect| (lognorm)")
axis.set_ylabel("log₁₀(1 + continuous spatial divergence score)")
axis.set_title("Continuous spatial signal versus whole-PT level effect", loc="left", fontweight="bold")
axis.text(0.98, 0.04, "Dashed: median whole-PT effect", transform=axis.transAxes,
          ha="right", va="bottom", fontsize=6.0, color="0.30")
axis.legend(fontsize=5.5, loc="upper right", markerscale=1.1)
_panel_label(axis, "a")

continuous_only = quadrant_frame[quadrant_labels["low conventional / strong spatial (continuous-only)"]]
for row in continuous_only.nlargest(5, "spatial_discovery_score").itertuples():
    axis.annotate(row.gene, (row.conventional_abs, np.log1p(row.spatial_discovery_score) / np.log(10)),
                  xytext=(2, 2), textcoords="offset points", fontsize=5.8)

count_matrix = np.array([
    [quadrant_counts.loc[quadrant_counts["quadrant"].eq("low conventional / strong spatial (continuous-only)"), "n_genes"].iloc[0],
     quadrant_counts.loc[quadrant_counts["quadrant"].eq("high conventional / strong spatial"), "n_genes"].iloc[0]],
    [quadrant_counts.loc[quadrant_counts["quadrant"].eq("low conventional / weak spatial"), "n_genes"].iloc[0],
     quadrant_counts.loc[quadrant_counts["quadrant"].eq("high conventional / weak spatial"), "n_genes"].iloc[0]],
], dtype=float)
count_image = count_axis.imshow(count_matrix, cmap="Blues", aspect="equal")
for row_index in range(2):
    for column_index in range(2):
        n_genes = int(count_matrix[row_index, column_index])
        count_axis.text(column_index, row_index, f"{n_genes:,}\n({n_genes / len(quadrant_frame):.1%})",
                        ha="center", va="center", fontsize=7.8,
                        color="white" if count_image.norm(n_genes) > 0.55 else "0.20",
                        fontweight="bold")
count_axis.add_patch(plt.Rectangle((-0.5, -0.5), 1, 1, fill=False,
                                   ec=PHENOTYPE_COLORS["human-zonated / mouse-flat"], lw=2.2))
count_axis.set_xticks([0, 1], ["Low\nwhole-PT effect", "High\nwhole-PT effect"], fontsize=6.8)
count_axis.set_yticks([0, 1], ["Strong, stable\nspatial evidence", "Weak / unstable\nspatial evidence"], fontsize=6.8)
count_axis.set_title("Quadrant counts", loc="left", fontweight="bold")
count_axis.tick_params(length=0)
_panel_label(count_axis, "b")

figure.suptitle(
    "Continuous spatial analysis identifies 289 stable-remodelling genes with below-median whole-PT effects\n"
    "The highlighted cell is the descriptive ‘continuous-only’ quadrant; no species-level p-values are implied",
    fontsize=9.1,
)
_save_figure(figure, "conventional_vs_continuous_genes.png")
plt.show()


# %%
# Purpose: the pathway-level spatially-only table, and the conventional comparator.

pathway_member_stats = {}

for row in pathway_tested.itertuples():
    members = np.asarray([index for index in row.universe_members if universe_mask[index]], dtype=int)
    if members.size == 0:
        continue

    label = f"{row.library}: {row.pathway}"
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

# --------------------------------------------------------------------------------------
# The fair conventional comparator: a whole-PT preranked GSEA on the *level effect*.
#
# The median absolute level effect of a pathway's members is a descriptive member statistic, not how
# a conventional analysis would test a pathway. Ranking every gene by its whole-PT level effect and
# running preranked GSEA gives the conventional answer the spatial result has to beat.
# --------------------------------------------------------------------------------------

import gseapy as gp

pathway_gene_sets = {
    f"{row.library}: {row.pathway}": [gene_names[index] for index in row.universe_members]
    for row in pathway_tested.itertuples()
    if len(row.universe_members) >= CONFIG["pathway_min_members"]
}

conventional_ranking = (
    gene_metrics.loc[discovery_eligible, ["gene", "level_effect_human_minus_mouse"]]
    .dropna()
    .sort_values("level_effect_human_minus_mouse", ascending=False)
    .reset_index(drop=True)
)

print(f"whole-PT preranked GSEA on the level-effect ranking: "
      f"{len(conventional_ranking):,} genes, {len(pathway_gene_sets):,} pathway sets")

conventional_gsea = gp.prerank(
    rnk=conventional_ranking,
    gene_sets=pathway_gene_sets,
    min_size=CONFIG["gsea_min_size"],
    max_size=CONFIG["gsea_max_size"],
    permutation_num=CONFIG["gsea_permutations"],
    seed=CONFIG["gsea_seed"],
    verbose=False,
)

conventional_results = (
    conventional_gsea.res2d.reset_index()
    .rename(columns={"Term": "pathway_key", "NES": "conventional_gsea_nes",
                     "NOM p-val": "conventional_gsea_p", "FDR q-val": "conventional_gsea_fdr",
                     "Lead_genes": "conventional_leading_genes"})
)
conventional_results[["library", "pathway"]] = [
    entry.split(": ", 1) if ": " in entry else (entry, entry)
    for entry in conventional_results["pathway_key"]
]
conventional_results["conventional_gsea_fdr"] = pd.to_numeric(
    conventional_results["conventional_gsea_fdr"], errors="coerce"
)
conventional_results["conventional_gsea_nes"] = pd.to_numeric(
    conventional_results["conventional_gsea_nes"], errors="coerce"
)
conventional_results = conventional_results.sort_values("conventional_gsea_fdr").reset_index(drop=True)

_report("pathways tested by the conventional GSEA", len(conventional_results),
        len(pathway_gene_sets),
        why=f"ranked by the whole-PT level effect; {CONFIG['gsea_permutations']} permutations")
_report("of those, conventionally significant (FDR < enrichment_fdr)",
        int((conventional_results["conventional_gsea_fdr"] < CONFIG["enrichment_fdr"]).sum()),
        len(conventional_results))

_save_table(conventional_results, "conventional_wholePT_preranked_gsea.csv")

conventional_by_pathway = {
    (row.library, row.pathway): (row.conventional_gsea_nes, row.conventional_gsea_fdr)
    for row in conventional_results.itertuples()
}

summary_columns = [
    "library", "pathway", "spatial_phenotype", "fold_enrichment", "fdr", "leading_genes",
    "n_members_in_universe",
    "median_member_abs_level_effect", "median_member_amplitude_log2_ratio",
    "median_member_abs_shift", "median_member_shape_corr", "median_member_spatial_score",
    "member_fraction_in_phenotype", "member_consistency",
    "conventional_gsea_nes", "conventional_gsea_fdr", "minus_log10_fdr",
]

pathway_spatial_rewiring_summary = pd.DataFrame([
    {
        "library": row.library,
        "pathway": row.pathway,
        "spatial_phenotype": row.spatial_phenotype,
        "fold_enrichment": row.fold_enrichment,
        "fdr": row.fdr,
        "leading_genes": row.leading_genes,
        # Coherence travels with the association: "enriched for genes classified as X" is the claim,
        # and these columns say how much of the pathway that is.
        "member_fraction_in_phenotype": row.member_fraction_in_phenotype,
        "member_consistency": row.member_consistency,
        "conventional_gsea_nes": conventional_by_pathway.get((row.library, row.pathway), (np.nan, np.nan))[0],
        "conventional_gsea_fdr": conventional_by_pathway.get((row.library, row.pathway), (np.nan, np.nan))[1],
        "minus_log10_fdr": row.minus_log10_fdr,
        **pathway_member_stats.get(f"{row.library}: {row.pathway}", {}),
    }
    for row in pathway_enrichment.itertuples()
], columns=summary_columns).sort_values(
    ["fdr", "fold_enrichment"], ascending=[True, False]
).reset_index(drop=True)

# An empty result set leaves object-dtype columns, and both the ranking below and every downstream
# comparison need them numeric.
pathway_spatial_rewiring_summary = _numeric(
    pathway_spatial_rewiring_summary,
    [column for column in summary_columns
     if column not in {"library", "pathway", "spatial_phenotype", "leading_genes",
                       "member_consistency"}],
)

_save_table(pathway_spatial_rewiring_summary, "pathway_spatial_rewiring_summary.csv")

significant_pathways = pathway_spatial_rewiring_summary[
    pathway_spatial_rewiring_summary["fdr"] < CONFIG["enrichment_fdr"]
].copy()

if significant_pathways.empty:
    # Keep the columns: everything downstream addresses this frame by name.
    best_per_pathway = significant_pathways.reindex(columns=significant_pathways.columns)
else:
    best_per_pathway = (
        significant_pathways
        .sort_values("fdr")
        .groupby(["library", "pathway"], as_index=False)
        .first()
    )

# Spatially significant but *not* conventionally significant, on the fair comparator. If this set is
# empty, that is the result: the conventional level analysis already ranks these pathways, and the
# continuous framework's contribution is how they differ rather than that they differ.
spatially_only_pathways = best_per_pathway[
    best_per_pathway["conventional_gsea_fdr"].isna()
    | best_per_pathway["conventional_gsea_fdr"].ge(CONFIG["enrichment_fdr"])
].sort_values("fdr").rename(columns={
    "spatial_phenotype": "dominant_spatial_phenotype",
    "fold_enrichment": "enrichment_fold_enrichment",
    "fdr": "enrichment_fdr",
    "median_member_amplitude_log2_ratio": "median_member_abs_amplitude_log2_ratio",
})[[
    # The conventional comparator's own result travels with each row, so a reader can see the NES and
    # FDR that the selection was made on rather than taking it on trust. Note what "spatially only"
    # does *not* mean: these pathways can carry large absolute member level effects. What the whole-PT
    # GSEA lacks is a coherent directional signal, which is exactly what reorganisation along PT
    # produces - the genes move in different directions, so their whole-PT shifts cancel.
    "pathway", "library", "conventional_gsea_nes", "conventional_gsea_fdr",
    "dominant_spatial_phenotype", "enrichment_fdr", "enrichment_fold_enrichment",
    "member_fraction_in_phenotype", "member_consistency",
    "median_member_abs_amplitude_log2_ratio", "median_member_abs_shift",
    "median_member_shape_corr", "leading_genes",
]]

_report("pathways with a significant spatial enrichment", len(best_per_pathway), len(pathway_tested))
_report("of those, not significant in the conventional whole-PT GSEA (spatially only)",
        len(spatially_only_pathways), len(best_per_pathway),
        why="the whole-PT level-effect GSEA finds no coherent direction for them")
if spatially_only_pathways.empty:
    print("  every spatially significant pathway also reaches significance in the conventional "
          "whole-PT GSEA, so the continuous framework's contribution for them is *how* they differ "
          "(zonation change, displacement, shape) rather than *that* they differ.")

_save_table(spatially_only_pathways, "spatially_only_pathways.csv")

comparator = best_per_pathway.dropna(subset=["conventional_gsea_nes"])

if comparator.empty:
    print("No pathway reached the spatial enrichment threshold, or none has a conventional GSEA "
          "result: the comparison figure is skipped.")
else:
    spatially_only = comparator[
        comparator["conventional_gsea_fdr"].isna()
        | comparator["conventional_gsea_fdr"].ge(CONFIG["enrichment_fdr"])
    ]
    both = comparator.drop(spatially_only.index)
    comparator = comparator.sort_values("fdr").reset_index(drop=True)
    spatial_only_mask = (comparator["conventional_gsea_fdr"].isna()
                         | comparator["conventional_gsea_fdr"].ge(CONFIG["enrichment_fdr"]))
    positions = np.arange(len(comparator))

    figure, axes = plt.subplots(1, 2, figsize=(8.5, 3.7), sharey=True,
                                width_ratios=[1.05, 0.92], layout="constrained")
    nes_axis, fdr_axis = axes
    colours = np.where(spatial_only_mask, PHENOTYPE_COLORS["human-zonated / mouse-flat"], "0.48")
    for position, row in enumerate(comparator.itertuples()):
        nes_axis.hlines(position, min(0, row.conventional_gsea_nes), 0, color="0.82", lw=1.0, zorder=1)
    nes_axis.scatter(comparator["conventional_gsea_nes"], positions, s=38, color=colours,
                     edgecolors="white", linewidths=0.5, zorder=3)
    nes_axis.axvline(0, color="0.35", lw=0.8)
    nes_axis.set_xlabel("Whole-PT level-effect GSEA NES\n(negative = lower human level ranking)")
    nes_axis.set_yticks(positions)
    nes_axis.set_yticklabels([textwrap.fill(row.pathway, 34) for row in comparator.itertuples()], fontsize=6.5)
    nes_axis.invert_yaxis()
    nes_axis.set_title("Conventional whole-PT ranking", loc="left", fontweight="bold")
    _panel_label(nes_axis, "a")

    fdr_axis.barh(positions, comparator["minus_log10_fdr"], color=colours, alpha=0.88)
    fdr_axis.axvline(-np.log10(CONFIG["enrichment_fdr"]), color="0.35", lw=0.8, ls="--")
    fdr_axis.set_xlabel("−log₁₀(spatial FDR)")
    fdr_axis.set_title("Spatial phenotype enrichment", loc="left", fontweight="bold")
    fdr_axis.tick_params(labelleft=False)
    fdr_axis.legend(handles=[
        Patch(facecolor=PHENOTYPE_COLORS["human-zonated / mouse-flat"], label=f"spatial only (n={len(spatially_only)})"),
        Patch(facecolor="0.48", label=f"significant in both (n={len(both)})"),
    ], loc="lower right", fontsize=6.0)
    _panel_label(fdr_axis, "b")

    figure.suptitle(
        "Three spatially significant pathway associations are not significant in whole-PT GSEA\n"
        "All six associations rank toward lower human whole-PT level, so the continuous result specifies their spatial pattern rather than a near-zero NES",
        fontsize=8.9,
    )
_save_figure(figure, "conventional_vs_continuous_pathway_signal.png")
plt.show()


# %% [markdown]
# ## 6 - Tables, the run inventory and the saved analysis notes
#
# ### 6.1 - The atlas, the per-class shortlist and the run inventory
#
# **The atlas** (`tables/gene_spatial_rewiring_atlas.csv`) is the per-gene result: level effect, both
# amplitudes and their ratio, the standardised curve correlation, the displacement with its improvement
# and residual (a supporting measurement, not a class), positional centroids and peaks, early-to-late
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
    "level_effect_human_minus_mouse",
    "mouse_amplitude", "human_amplitude", "amplitude_difference_human_minus_mouse",
    "amplitude_log2_ratio_human_over_mouse",
    "shape_corr", "shape_spearman", "pattern_rms_z",
    "post_shift_shape_corr",
    "best_shift_human_minus_mouse",
    "shift_improvement", "residual_rms_after_shift",
    "shift_at_search_boundary", "shift_usable",
    "mouse_position_centroid", "human_position_centroid", "centroid_shift_human_minus_mouse",
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

# Ranked by each class's own evidence, exactly as the figures are, so the table and the panels agree.
top_genes_by_spatial_phenotype = pd.concat(
    [_top_genes_for_class(label, TOP_N_PER_CLASS) for label in PHENOTYPE_CLASSES],
    ignore_index=True,
)
top_genes_by_spatial_phenotype = top_genes_by_spatial_phenotype.assign(
    exemplar_evidence=[
        CLASS_EVIDENCE[label]["label"] if label in CLASS_EVIDENCE else "global divergence score"
        for label in top_genes_by_spatial_phenotype["spatial_phenotype"]
    ]
)[[
    "spatial_phenotype", "gene", "exemplar_evidence", "robustness_score",
    "spatial_discovery_score", "shift_interior", "shift_at_search_boundary",
    "mouse_amplitude", "human_amplitude", "amplitude_log2_ratio_human_over_mouse",
    "shape_corr", "post_shift_shape_corr", "best_shift_human_minus_mouse",
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

print(f"Tables ({len(table_files)}) in {_display_path(TABLE_DIR)}:")
for name in table_files:
    print(f"  {name}")
print(f"\nFigures ({len(figure_files)}) in {_display_path(FIG_DIR)}:")
for name in figure_files:
    print(f"  {name}")
print(f"\nCurve payloads ({len(object_files)}) in {_display_path(OUT_DIR)}:")
for name in object_files:
    print(f"  {name}")

# --------------------------------------------------------------------------------------
# Analysis notes, written from this run's own values.
# --------------------------------------------------------------------------------------

analysis_notes = f"""# Human versus healthy-mouse PT spatial rewiring (notebook 06)

Written by `analysis/notebooks/06_human_mouse_spatial_rewiring.ipynb` at run time from the values it
computed. Read together with the notebook's own section markdown, which states the same caveats.

## Inputs and cohort

- Input object: `{_display_path(DPT_OUTPUT_PATH)}` (03's PT-specific DPT; `total_scanpy_dpt`).
- Coordinate: 03's shared PT DPT (`total_scanpy_dpt`), used as 03 oriented it. No cross-species
  re-registration is applied, so a displacement is reported as a supporting measurement of the genes
  carrying it (`shift_usable`) rather than as a phenotype class.
- Mouse specimens: {', '.join(mouse_samples)}. Human slices: {', '.join(human_samples)} - **one donor**.
- PT structures compared: {adata_pt.n_obs:,}. Shared pseudospace support: [{lo:.4f}, {hi:.4f}] on a
  {grid.size}-point grid; displacements are reported in unit coordinates over that support.
- Tested orthologs: {len(gene_names):,} (>= {CONFIG['min_detected_fraction_all']:.0%} detection overall
  and measured in both inputs); discovery-eligible: {int(discovery_eligible.sum()):,}.

## Method, in one paragraph

Per-specimen curves are fitted with the shared-smooth model at the pooled per-gene smoothing parameter
and averaged with equal specimen weight; the pooled fit is 03/05's cached nested level/shape fit
(`section4_gene_fit`). Metrics are level, amplitude, standardised shape correlation, positional
centroid/peak/early-to-late gradient/half-max window, and a bounded displacement search (interpreted
only to 0.30 of PT, the interior limit; the search itself looks to {CONFIG['shift_search_limit']} of
PT). Every metric is computed on the single coordinate of record, 03's shared PT DPT.

## Caveats that constrain every number here

- The two human slices are one donor. Nothing here is a population-level species test; no p-value on a
  cross-species contrast is reported or implied. Slice agreement is robustness, not replication.
- Species is inseparable from sampling and from Harmony batch correction (`sample` is the batch key).
- Both human slices are cortex; the `MED1` label does not make the tissue medullary.
- Displacement is read on 03's shared PT DPT without cross-species re-registration: it is evidence
  that the two species' programmes differ in position, not an anatomical distance between aligned
  axes. A displacement is reported only when the optimum is interior, the positional centroid agrees
  and the shape agreement improves (`shift_usable`), and it never defines a phenotype class.
- Phenotype classes are threshold conjunctions. The thresholds are in the notebook's `CONFIG` and the
  distributions they act on are printed in section 3; the component metrics are in the atlas.

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
- Stage cache: `{_display_path(STAGE_CACHE_DIR)}` - 03/05's directory, deliberately shared
  so the pooled gene fit is reused rather than recomputed; this notebook's own stages there are
  `spatial_rewiring_specimen_curves` (per-specimen fits) beside the inherited `section4_gene_fit`.
  `PSEUDOSPACE_STAGE_CACHE=0` forces a full rebuild.
"""

(OUT_DIR / "analysis_notes.md").write_text(analysis_notes)
print(f"\nsaved {_display_path(OUT_DIR / 'analysis_notes.md')}")


# %% [markdown]
# ### 6.2 - Closing summary
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
summary_lines.append("Coordinate: 03's shared PT DPT (`total_scanpy_dpt`), used as 03 oriented it")
summary_lines.append("Displacement is a supporting measurement, not a class: "
                     f"{_count(gene_metrics['shift_usable'])} genes meet every displacement criterion")
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
if len(architecture_null):
    for row in architecture_null.itertuples():
        summary_lines.append(
            f"Conserved architecture ({row.n_genes:,} genes patterned in both species): "
            f"{row.statistic} = {row.observed:.3f} vs null {row.null_median:.3f} "
            f"(empirical p = {row.empirical_p:.3g})"
        )
else:
    summary_lines.append("Conserved architecture: not assessed (too few both-species-patterned genes)")
summary_lines.append("")
summary_lines.append(f"Significant pathway x phenotype enrichments (FDR < {CONFIG['enrichment_fdr']}): "
                     f"{len(best_per_pathway):,} pathways of {len(pathway_tested):,} tested")
for label in PHENOTYPE_CLASSES:
    block = best_per_pathway[best_per_pathway["spatial_phenotype"].eq(label)]
    if len(block):
        summary_lines.append(f"  {label}: {len(block):,} pathways")
summary_lines.append("")
summary_lines.append("Spatially significant, conventional-GSEA non-significant pathways:")
summary_lines.append(f"  {len(spatially_only_pathways):,} of {len(best_per_pathway):,} spatially "
                     "significant pathways are not significant in the whole-PT level-effect GSEA")
for row in spatially_only_pathways.head(8).itertuples():
    summary_lines.append(
        f"  {row.pathway} [{row.library}] -> {row.dominant_spatial_phenotype} "
        f"spatial FDR={row.enrichment_fdr:.3g}, conventional GSEA NES={row.conventional_gsea_nes:+.2f} "
        f"(FDR {row.conventional_gsea_fdr:.2g}), members {row.member_fraction_in_phenotype:.0%} in "
        f"phenotype"
    )
summary_lines.append("  (these pathways can carry large absolute member level effects; what the "
                     "whole-PT GSEA lacks is a coherent direction, which reorganisation along PT "
                     "produces)")
summary_lines.append("")
summary_lines.append("Caveats: two mouse specimens and two human slices from ONE donor; no "
                     "population-level species inference is reported; species is confounded with "
                     "sampling and batch correction; the coordinate is 03's shared PT DPT with no "
                     "cross-species re-registration, and displacement is reported as a supporting "
                     "measurement rather than as a phenotype class.")

summary_text = "\n".join(summary_lines)
print(summary_text)
(OUT_DIR / "analysis_summary.txt").write_text(summary_text + "\n")
print(f"\nsaved {_display_path(OUT_DIR / 'analysis_summary.txt')}")


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

    later_shift, later_corr, later_residual, later_at_boundary = _row_shift_search(
        bump, later, x, 0.35, 0.01, min_points=10
    )
    earlier_shift, _, _, _ = _row_shift_search(bump, earlier, x, 0.35, 0.01, min_points=10)

    checks.append(("a later human program is recovered as a positive displacement",
                   abs(later_shift[0] - 0.10) <= 0.02))
    checks.append(("an earlier human program is recovered as a negative displacement",
                   abs(earlier_shift[0] + 0.10) <= 0.02))
    checks.append(("recovering the displacement improves the correlation",
                   later_corr[0] > _row_pearson(bump, later)[0]))
    checks.append(("residual after the best displacement is ~0 for a pure move",
                   later_residual[0] < 0.05))

    _, inverted_corr, _, _ = _row_shift_search(bump, inverted, x, 0.35, 0.01, min_points=10)
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
