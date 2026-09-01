"""Human versus healthy-mouse nephron pseudospace reconstruction.

This is the cross-species counterpart to ``analysis/mouse_only_pseudospace.py``.
It uses the two mouse control samples (Ctrl1A2/Ctrl1A4) and the available human
samples (HUK1_COR1/HUK1_MED1), converts human expression into a documented
one-to-one mouse-ortholog space, reconstructs one shared trajectory, and compares
gene/module curves descriptively.

The human cohort has one donor and cortex/medulla samples are confounded with
species.  Consequently, species curves and level/shape quantities are effect-size
descriptions, not powered species tests.  The generated ``analysis_notes.md``
repeats these caveats next to the results.

Run from the repository root, for example::

    python analysis/human_vs_healthy_mouse.py \
        --data-root /private/project-data \
        --results-root /private/project-results

The notebook ``analysis/notebooks/03_human_vs_healthy_mouse.ipynb`` runs this same
script from an interactive kernel.
"""
# %% [markdown]
# # Human versus healthy-mouse pseudospace
#
# This workflow mirrors the current mouse-only reconstruction while making the
# cross-species conversion explicit. Raw matrices and segmentations remain outside
# Git; use `--data-root` and `--results-root` when running locally.

# %%
from __future__ import annotations

import argparse
import json
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=FutureWarning)


def _project_dir() -> Path:
    starts = [Path(__file__).resolve().parent, Path.cwd().resolve()]
    for start in starts:
        for candidate in (start, *start.parents):
            if (candidate / "pseudospace").is_dir() and (candidate / "data").is_dir():
                return candidate
    raise RuntimeError("Could not locate repository root containing pseudospace/ and data/")


def _roots(project_dir: Path) -> tuple[Path, Path]:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--data-root", type=Path)
    parser.add_argument("--results-root", type=Path)
    args, _ = parser.parse_known_args()
    data = args.data_root or os.environ.get("PSEUDOSPACE_DATA_ROOT") or project_dir / "data"
    results = args.results_root or os.environ.get("PSEUDOSPACE_RESULTS_ROOT") or project_dir / "results"
    return Path(data).expanduser().resolve(), Path(results).expanduser().resolve()


PROJECT_DIR = _project_dir()
if str(PROJECT_DIR) not in sys.path:
    sys.path.insert(0, str(PROJECT_DIR))
DATA_ROOT, RESULTS_ROOT = _roots(PROJECT_DIR)
RESULTS_DIR = RESULTS_ROOT / "human_vs_healthy_mouse"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
for subdir in ("heatmaps", "curves", "diagnostics"):
    (RESULTS_DIR / subdir).mkdir(exist_ok=True)

CACHE_ROOT = Path(os.environ.get("PSEUDOSPACE_CACHE_ROOT", "/tmp/pseudospace_human_mouse_cache"))
for subdir in ("matplotlib", "numba", "xdg"):
    (CACHE_ROOT / subdir).mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(CACHE_ROOT / "matplotlib"))
os.environ.setdefault("NUMBA_CACHE_DIR", str(CACHE_ROOT / "numba"))
os.environ.setdefault("XDG_CACHE_HOME", str(CACHE_ROOT / "xdg"))

TUBULE_BY_GENE_DIR = DATA_ROOT / "tubule_by_gene"
ORTHOLOG_TABLE_PATH = DATA_ROOT / "human_mouse_hcop_fifteen_column.txt.gz"
PATHWAY_LIBRARY_DIR = DATA_ROOT / "mouse_vs_human" / "pathway_gene_sets"
MOUSE_SAMPLES = ("Ctrl1A2", "Ctrl1A4")
HUMAN_SAMPLES = ("HUK1_COR1", "HUK1_MED1")
HUMAN_REGIONS = {"HUK1_COR1": "cortex", "HUK1_MED1": "medulla"}
USE_HARMONY = os.environ.get("PSEUDOSPACE_USE_HARMONY", "1").lower() not in {"0", "false", "no"}
RANDOM_STATE = 0
MIN_GENES_PER_TUBULE = 100
MIN_GENE_TOTAL_COUNTS = 20
HARMONY_THETA = 6
N_PCS = 50
N_NEIGHBORS = 30
N_BINS = 120
GAM_LAMBDA_GRID = np.logspace(-3, 3, 13)

# The panel is intentionally in mouse symbols: human data are converted to this
# space before a marker is ever used for labels, axes, or heatmaps.
SEGMENT_MARKERS = {
    "Podocyte": ["Nphs1", "Nphs2", "Podxl"],
    "PT": ["Slc5a2", "Slc22a6", "Slc13a3", "Slc7a13"],
    "Thin_limb": ["Slc14a2", "Clcnka", "Cldn10"],
    "TAL": ["Slc12a1", "Umod", "Cldn16"],
    "DCT": ["Slc12a3", "Trpm6", "Trpv5"],
    "CD": ["Aqp2", "Aqp4", "Atp6v1b1"],
}
POSITION_MARKERS = {
    "early": ["Nphs1", "Slc5a2", "Slc22a6"],
    "late": ["Aqp2", "Aqp4", "Atp6v1b1"],
}
HEATMAP_MARKERS = {
    "PT-S1": ["Slc5a2", "Slc5a12"],
    "PT-S2": ["Slc22a6", "Slc13a3"],
    "PT-S3": ["Slc7a13", "Slc22a7"],
    "TAL": ["Slc12a1", "Umod"],
    "DCT": ["Slc12a3", "Trpm6"],
    "CD": ["Aqp2", "Aqp4"],
}
SUBSEGMENT_MARKERS = {
    "PT": {
        "PT-S1": ["Slc5a2", "Slc5a12"],
        "PT-S2": ["Slc22a6", "Slc13a3"],
        "PT-S3": ["Slc7a13", "Slc22a7"],
    },
    "TAL": {
        "mTAL": ["Slc12a1", "Umod", "Cldn10"],
        "cTAL": ["Slc12a1", "Umod", "Cldn16"],
    },
    "DCT": {
        "DCT1": ["Slc12a3", "Trpm6"],
        "DCT2": ["Slc12a3", "Trpv5"],
    },
}


def _marker_dict_for_heatmap(markers: dict[str, list[str]]) -> dict[str, list[dict[str, object]]]:
    return {
        group: [{"label": gene, "candidates": [gene]} for gene in genes]
        for group, genes in markers.items()
    }


def _score_and_label(adata, marker_panels: dict[str, list[str]], cluster_col: str = "leiden"):
    """Score available markers and label clusters by their mean panel score."""
    from scipy import sparse

    present = {group: [g for g in genes if g.upper() in {str(v).upper() for v in adata.var_names}]
               for group, genes in marker_panels.items()}
    present = {group: genes for group, genes in present.items() if genes}
    if len(present) < 2:
        raise ValueError("Fewer than two segment marker panels resolved in shared gene space")
    X = adata.layers["lognorm"] if "lognorm" in adata.layers else adata.X
    scores = {}
    lookup = {str(g).upper(): i for i, g in enumerate(adata.var_names)}
    for group, genes in present.items():
        idx = [lookup[g.upper()] for g in genes]
        values = X[:, idx]
        if sparse.issparse(values):
            values = values.toarray()
        scores[group] = np.asarray(values, dtype=float).mean(axis=1)
        adata.obs[f"marker_score_{group}"] = scores[group]
    score_df = pd.DataFrame(scores, index=adata.obs_names)
    cluster_means = score_df.groupby(adata.obs[cluster_col].astype(str)).mean()
    labels = cluster_means.idxmax(axis=1)
    adata.obs["broad_tubule_marker_call"] = adata.obs[cluster_col].astype(str).map(labels).fillna("Unknown").to_numpy()
    adata.obs["celltype_margin"] = score_df.apply(lambda row: np.partition(row.to_numpy(), -2)[-1] - np.partition(row.to_numpy(), -2)[-2], axis=1)
    report = cluster_means.reset_index(names="cluster")
    report["assigned_label"] = report.drop(columns="cluster").idxmax(axis=1)
    return present, report


def _binned_profiles(adata, genes: list[str], species: str, n_bins: int = N_BINS):
    from scipy import sparse

    mask = adata.obs["comparison_species"].astype(str).eq(species).to_numpy()
    sub = adata[mask]
    x = sub.obs["shared_pseudospace"].to_numpy(dtype=float)
    valid = np.isfinite(x)
    x = x[valid]
    sub = sub[valid]
    if len(x) < 10:
        raise ValueError(f"Too few finite pseudotimes for {species} profile")
    lo, hi = float(np.min(x)), float(np.max(x))
    if np.isclose(lo, hi):
        lo, hi = lo - 0.005, hi + 0.005
    grid = np.linspace(lo, hi, n_bins)
    edges = np.linspace(lo, hi, n_bins + 1)
    X = sub[:, genes].layers["lognorm"] if "lognorm" in sub.layers else sub[:, genes].X
    if sparse.issparse(X):
        X = X.toarray()
    values = np.full((len(genes), n_bins), np.nan)
    ids = np.clip(np.digitize(x, edges) - 1, 0, n_bins - 1)
    for b in range(n_bins):
        hit = ids == b
        if hit.any():
            values[:, b] = np.asarray(X[hit]).mean(axis=0)
    for i in range(values.shape[0]):
        good = np.isfinite(values[i])
        if good.sum() >= 2:
            values[i] = np.interp(grid, grid[good], values[i, good])
    return grid, values


# %% [markdown]
# ## 1. Load, QC, and convert to shared ortholog space

# %%
from pseudospace.cross_species import (
    build_one_to_one_ortholog_map,
    combine_cross_species,
    load_cross_species_samples,
    prepare_shared_expression,
    read_ortholog_table,
)

if not ORTHOLOG_TABLE_PATH.exists():
    raise FileNotFoundError(
        f"Missing {ORTHOLOG_TABLE_PATH}. Acquire the private HCOP table and place it under --data-root."
    )
ortholog_table = read_ortholog_table(ORTHOLOG_TABLE_PATH)
ambiguous_human_symbols = ortholog_table.loc[
    ortholog_table.groupby(ortholog_table["human_symbol"].str.upper())["mouse_symbol"]
    .transform("nunique").gt(1), "human_symbol"
].astype(str).unique().tolist()
# Curated overrides may be supplied as a two-column CSV through the environment;
# keeping them external prevents unreviewed symbols from entering Git.
override_path = os.environ.get("PSEUDOSPACE_ORTHOLOG_OVERRIDES")
overrides = {}
if override_path:
    override_df = pd.read_csv(override_path)
    overrides = dict(zip(override_df["human_symbol"], override_df["mouse_symbol"]))
ortholog_map = build_one_to_one_ortholog_map(ortholog_table, overrides=overrides)
ortholog_map.to_csv(RESULTS_DIR / "ortholog_map_used.csv", index=False)

files = sorted(TUBULE_BY_GENE_DIR.glob("*_tubule_by_gene_caleb.h5ad"))
if not files:
    raise FileNotFoundError(f"No tubule matrices found in {TUBULE_BY_GENE_DIR}")
adatas, human_mapping_report = load_cross_species_samples(
    TUBULE_BY_GENE_DIR,
    files,
    ortholog_map,
    mouse_samples=MOUSE_SAMPLES,
    human_samples=HUMAN_SAMPLES,
    human_regions=HUMAN_REGIONS,
    ambiguous_human_symbols=ambiguous_human_symbols,
)
human_mapping_report.to_csv(RESULTS_DIR / "human_gene_mapping_report.csv", index=False)
adata = prepare_shared_expression(
    combine_cross_species(adatas),
    min_genes=MIN_GENES_PER_TUBULE,
    min_gene_total_counts=MIN_GENE_TOTAL_COUNTS,
)
cohort_summary = (
    adata.obs.groupby(["comparison_species", "sample", "region"], observed=True)
    .size().rename("n_tubules").reset_index()
)
cohort_summary.to_csv(RESULTS_DIR / "cohort_summary.csv", index=False)
print(f"Shared matrix: {adata.n_obs:,} tubules x {adata.n_vars:,} mouse-ortholog genes")
print(cohort_summary.to_string(index=False))

# %% [markdown]
# ## 2. Shared embedding, integration, and cell/segment labels
#
# HVGs are selected separately in mouse and human and intersected before PCA. Harmony
# integrates sample labels when R/harmony is available. If it is unavailable, the
# script records a PCA fallback rather than silently pretending integration occurred.

# %%
import scanpy as sc

adata.obs["integration_group"] = adata.obs["comparison_species"].astype(str)
from pseudospace.harmony import run_harmony_rpy2, select_harmony_hvgs_by_condition

adata = select_harmony_hvgs_by_condition(
    adata,
    group_key="integration_group",
    groups=("mouse", "human"),
    mode="intersection",
    min_mean=0.0125,
    max_mean=3,
    min_disp=0.5,
)
hvg = adata.var["highly_variable_for_harmony"].to_numpy()
if int(hvg.sum()) < 50:
    raise ValueError(f"Only {int(hvg.sum())} shared Harmony HVGs resolved")
adata_hvg = adata[:, hvg].copy()
adata_hvg.X = adata_hvg.layers["lognorm"].copy()
sc.pp.scale(adata_hvg, zero_center=True, max_value=10)
sc.tl.pca(adata_hvg, n_comps=min(N_PCS, adata.n_obs - 1, int(hvg.sum()) - 1), random_state=RANDOM_STATE)
adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"].copy()
integration_method = "PCA fallback"
if USE_HARMONY:
    try:
        harmony_result = run_harmony_rpy2(
            adata_hvg, batch_key="sample", n_pcs=min(N_PCS, adata_hvg.obsm["X_pca"].shape[1]),
            theta=HARMONY_THETA, lambda_val=1, max_iter=30, tau=0,
        )
        adata.obsm["X_harmony"] = harmony_result.obsm["X_harmony"].copy()
        integration_method = "Harmony (sample batch)"
    except Exception as exc:
        print(f"WARNING: Harmony unavailable; retaining PCA embedding ({type(exc).__name__}: {exc})")
rep_key = "X_harmony" if "X_harmony" in adata.obsm else "X_pca"
adata.uns["cross_species_integration"] = {
    "method": integration_method,
    "batch_key": "sample" if rep_key == "X_harmony" else None,
    "species_is_confounded_with_donor": True,
}
sc.pp.neighbors(adata, n_neighbors=min(N_NEIGHBORS, adata.n_obs - 1), use_rep=rep_key, random_state=RANDOM_STATE)
sc.tl.leiden(adata, resolution=0.5, key_added="leiden", random_state=RANDOM_STATE)
resolved_panels, cluster_report = _score_and_label(adata, SEGMENT_MARKERS)
cluster_report.to_csv(RESULTS_DIR / "segment_cluster_marker_scores.csv", index=False)
pd.DataFrame([
    {"panel": panel, "resolved_genes": "; ".join(genes), "n_resolved": len(genes)}
    for panel, genes in resolved_panels.items()
]).to_csv(RESULTS_DIR / "segment_marker_resolution.csv", index=False)
adata.write_h5ad(RESULTS_DIR / "cross_species_harmony_pass1.h5ad")
print(f"Embedding: {rep_key}; integration: {integration_method}")
print(adata.obs["broad_tubule_marker_call"].value_counts().to_string())

# %% [markdown]
# ## 2b. Pass-2 integration on the retained tubule continuum
#
# Labels are intentionally not recomputed after this pass, matching the mouse-only workflow.
# Human and mouse are still represented in the shared ortholog space, and sample is the batch
# variable for Harmony. The pass-2 embedding is used only for neighbors/DPT.

# %%
tubule_mask = adata.obs["broad_tubule_marker_call"].astype(str).ne("Unknown").to_numpy()
adata_tubule = adata[tubule_mask].copy()
adata_tubule_hvg = adata_tubule[:, hvg].copy()
adata_tubule_hvg.X = adata_tubule_hvg.layers["lognorm"].copy()
sc.pp.scale(adata_tubule_hvg, zero_center=True, max_value=10)
sc.tl.pca(
    adata_tubule_hvg,
    n_comps=min(N_PCS, adata_tubule.n_obs - 1, int(hvg.sum()) - 1),
    random_state=RANDOM_STATE,
)
adata_tubule.obsm["X_pca"] = adata_tubule_hvg.obsm["X_pca"].copy()
pass2_method = "PCA fallback"
if USE_HARMONY:
    try:
        pass2_result = run_harmony_rpy2(
            adata_tubule_hvg,
            batch_key="sample",
            n_pcs=min(N_PCS, adata_tubule_hvg.obsm["X_pca"].shape[1]),
            theta=HARMONY_THETA,
            lambda_val=1,
            max_iter=30,
            tau=0,
        )
        adata_tubule.obsm["X_harmony"] = pass2_result.obsm["X_harmony"].copy()
        pass2_method = "Harmony (sample batch)"
    except Exception as exc:
        print(f"WARNING: pass-2 Harmony unavailable; retaining PCA ({type(exc).__name__}: {exc})")
pass2_rep = "X_harmony" if "X_harmony" in adata_tubule.obsm else "X_pca"
for key in ("X_pca", "X_harmony"):
    if key in adata_tubule.obsm:
        adata.obsm[key] = np.full((adata.n_obs, adata_tubule.obsm[key].shape[1]), np.nan, dtype=float)
        adata.obsm[key][tubule_mask] = adata_tubule.obsm[key]
adata.uns["cross_species_pass2_integration"] = {"method": pass2_method, "representation": pass2_rep}
adata_tubule.write_h5ad(RESULTS_DIR / "cross_species_harmony_pass2.h5ad")
print(f"Pass-2 embedding: {pass2_rep}; retained tubules: {adata_tubule.n_obs:,}")

# %% [markdown]
# ## 3. Shared diffusion pseudotime and orientation checks

# %%
from pseudospace.markers import compute_total_marker_axis
from pseudospace.trajectory import choose_root_global, orient_and_normalize

axis = compute_total_marker_axis(
    adata,
    POSITION_MARKERS,
    output_col="total_marker_axis",
    layer="lognorm",
    min_markers_per_group=2,
)
rep_key = pass2_rep
sc.pp.neighbors(adata, n_neighbors=min(N_NEIGHBORS, adata.n_obs - 1), use_rep=rep_key, random_state=RANDOM_STATE)
sc.tl.diffmap(adata)
root = choose_root_global(adata.obs["total_marker_axis"].to_numpy(), adata, rep_key, bottom_quantile=0.01)
adata.uns["iroot"] = int(root)
sc.tl.dpt(adata)
adata.obs["total_scanpy_dpt"] = orient_and_normalize(
    adata.obs["dpt_pseudotime"].to_numpy(dtype=float),
    adata.obs["total_marker_axis"].to_numpy(dtype=float),
    min_valid=20,
)
adata.obs["shared_pseudospace"] = adata.obs["total_scanpy_dpt"].to_numpy(dtype=float)
adata.write_h5ad(RESULTS_DIR / "cross_species_dpt.h5ad")
pd.DataFrame({
    "sample": adata.obs["sample"].astype(str),
    "species": adata.obs["comparison_species"].astype(str),
    "region": adata.obs["region"].astype(str),
    "segment": adata.obs["broad_tubule_marker_call"].astype(str),
    "shared_pseudospace": adata.obs["shared_pseudospace"].to_numpy(),
}).to_csv(RESULTS_DIR / "dpt_by_tubule.csv", index=False)

# Recompute family-specific axes where the marker panel has at least two ordered modules.
# Failures are recorded rather than hidden: a small/low-quality family should not invalidate
# the global trajectory or be mistaken for a biological absence of an axis.
from pseudospace.trajectory import recompute_subset_dpt

subset_dpt_rows = []
for family, family_markers in SUBSEGMENT_MARKERS.items():
    try:
        family_report = recompute_subset_dpt(
            adata,
            family,
            _marker_dict_for_heatmap(family_markers),
            list(family_markers),
            output_col=f"{family.lower()}_subset_dpt",
            n_neighbors=min(N_NEIGHBORS, adata.n_obs - 1),
            random_state=RANDOM_STATE,
        )
        subset_dpt_rows.append(family_report)
    except (KeyError, ValueError) as exc:
        subset_dpt_rows.append({"segment": family, "status": f"skipped: {exc}"})
        print(f"WARNING: {family} subset DPT skipped: {exc}")
pd.DataFrame(subset_dpt_rows).to_csv(RESULTS_DIR / "dpt_by_segment.csv", index=False)

# %% [markdown]
# ## 4. Marker heatmaps and species profile similarity

# %%
import matplotlib.pyplot as plt
from pseudospace.heatmaps import plot_marker_heatmap

heatmap_panel = _marker_dict_for_heatmap(HEATMAP_MARKERS)
heatmap_colors = {name: color for name, color in zip(HEATMAP_MARKERS, ["#2166ac", "#67a9cf", "#d1e5f0", "#f4a582", "#d6604d", "#b2182b"])}
for species in ("mouse", "human"):
    try:
        plot_marker_heatmap(
            adata,
            heatmap_panel,
            list(HEATMAP_MARKERS),
            heatmap_colors,
            pseudotime_col="shared_pseudospace",
            cluster_col=None,
            mask=adata.obs["comparison_species"].astype(str).eq(species).to_numpy(),
            title=f"{species.title()} marker gradients on shared pseudospace",
            output_name=f"{species}_marker_heatmap.png",
            strip_col="broad_tubule_marker_call",
            n_bins=N_BINS,
            output_dir=RESULTS_DIR / "heatmaps",
            project_dir=PROJECT_DIR,
        )
    except ValueError as exc:
        print(f"WARNING: {species} heatmap skipped: {exc}")

for family, family_markers in SUBSEGMENT_MARKERS.items():
    pseudotime_col = f"{family.lower()}_subset_dpt"
    if pseudotime_col not in adata.obs:
        continue
    for species in ("mouse", "human"):
        try:
            plot_marker_heatmap(
                adata,
                _marker_dict_for_heatmap(family_markers),
                list(family_markers),
                {name: color for name, color in zip(family_markers, ["#2166ac", "#67a9cf", "#b2182b"])},
                pseudotime_col=pseudotime_col,
                cluster_col=None,
                mask=(adata.obs["comparison_species"].astype(str).eq(species)
                      & adata.obs["broad_tubule_marker_call"].astype(str).eq(family)).to_numpy(),
                title=f"{species.title()} {family} subset marker gradients",
                output_name=f"{species}_{family.lower()}_subset_heatmap.png",
                strip_col="sample",
                n_bins=min(N_BINS, 80),
                output_dir=RESULTS_DIR / "heatmaps",
                project_dir=PROJECT_DIR,
            )
        except ValueError as exc:
            print(f"WARNING: {species} {family} subset heatmap skipped: {exc}")

profile_genes = [gene for genes in HEATMAP_MARKERS.values() for gene in genes if gene in adata.var_names]
profile_genes = list(dict.fromkeys(profile_genes))
if len(profile_genes) >= 3:
    profile_rows = []
    mouse_grid, mouse_values = _binned_profiles(adata, profile_genes, "mouse")
    human_grid, human_values = _binned_profiles(adata, profile_genes, "human")
    common_grid = np.linspace(max(mouse_grid.min(), human_grid.min()), min(mouse_grid.max(), human_grid.max()), N_BINS)
    mouse_common = np.vstack([np.interp(common_grid, mouse_grid, row) for row in mouse_values])
    human_common = np.vstack([np.interp(common_grid, human_grid, row) for row in human_values])
    for gene_i, gene in enumerate(profile_genes):
        for x, mouse_value, human_value in zip(common_grid, mouse_common[gene_i], human_common[gene_i]):
            profile_rows.append({"gene": gene, "pseudospace": x, "mouse_lognorm": mouse_value, "human_lognorm": human_value})
    pd.DataFrame(profile_rows).to_csv(RESULTS_DIR / "curves" / "marker_profile_similarity.csv", index=False)

# %% [markdown]
# ## 5. Gene-level level/shape curves
#
# The nested GAM is retained as a descriptive decomposition. Its cellwise F quantities are
# not presented as confirmatory p-values because one human donor is compared with two mouse
# specimens and species is inseparable from donor/region in this cohort.

# %%
from scipy import sparse
from pseudospace.levelshape import run_level_shape
from pseudospace.stats_gam import gam_internal_knots

pt_mask = adata.obs["broad_tubule_marker_call"].astype(str).eq("PT").to_numpy()
pt_mask &= np.isfinite(adata.obs["shared_pseudospace"].to_numpy(dtype=float))
adata_pt = adata[pt_mask].copy()
if adata_pt.n_obs < 30 or adata_pt.obs["comparison_species"].nunique() < 2:
    raise ValueError("PT cohort is too small for a cross-species trajectory comparison")
s = adata_pt.obs["shared_pseudospace"].to_numpy(dtype=float)
c = adata_pt.obs["comparison_species"].astype(str).eq("human").astype(int).to_numpy()
gene_matrix = adata_pt.layers["lognorm"] if "lognorm" in adata_pt.layers else adata_pt.X
if not sparse.issparse(gene_matrix):
    gene_matrix = sparse.csr_matrix(np.asarray(gene_matrix))
gene_detected = np.asarray((gene_matrix > 0).sum(axis=0)).ravel()
gene_mean = np.asarray(gene_matrix.mean(axis=0)).ravel()
gene_keep = (gene_detected >= max(10, int(np.ceil(0.02 * adata_pt.n_obs)))) & (gene_mean > 0)
gene_names = adata_pt.var_names.to_numpy()[gene_keep]
Y_genes = gene_matrix[:, gene_keep].tocsr().astype(np.float64)
support = []
for species in ("mouse", "human"):
    values = s[c == (species == "human")]
    support.append((np.nanmin(values), np.nanmax(values)))
lo, hi = max(x[0] for x in support), min(x[1] for x in support)
if not lo < hi:
    raise ValueError(f"No overlapping pseudospace support between species: {support}")
grid = np.linspace(lo, hi, 101)
knots = gam_internal_knots(s)
ls = run_level_shape(Y_genes, s, c, knots, grid, GAM_LAMBDA_GRID)
gene_results = pd.DataFrame({
    "gene": gene_names,
    "species_effect_rms": ls["condition_effect_rms"],
    "species_effect_max_abs": ls["condition_effect_max_abs"],
    "descriptive_level_effect_human_minus_mouse": ls["level_effect"],
    "descriptive_shape_rms": ls["shape_rms"],
    "curve_spearman": ls["curve_spearman"],
    "mouse_amplitude": ls["amplitude_healthy"],
    "human_amplitude": ls["amplitude_aki"],
})
gene_results["trajectory_summary"] = np.select(
    [gene_results["curve_spearman"] >= 0.7, gene_results["curve_spearman"] < 0],
    ["similar fitted shape", "opposite fitted directions"],
    default="dissimilar fitted shape",
)
gene_results = gene_results.sort_values("species_effect_rms", ascending=False)
gene_results.to_csv(RESULTS_DIR / "curves" / "gene_trajectory_comparison.csv", index=False)
top_genes = gene_results.head(40)["gene"].tolist()
top_idx = [int(np.flatnonzero(gene_names == gene)[0]) for gene in top_genes]
curve_rows = []
for species, curves in (("mouse", ls["curve_healthy"]), ("human", ls["curve_aki"])):
    for gene_i, gene in zip(top_idx, top_genes):
        for x, value in zip(grid, curves[gene_i]):
            curve_rows.append({"gene": gene, "species": species, "pseudospace": x, "fitted_lognorm": value})
pd.DataFrame(curve_rows).to_csv(RESULTS_DIR / "curves" / "top_gene_fitted_curves.csv", index=False)

from pseudospace.levelshape import fit_single_condition_curves

sample_curve_rows = []
for sample in sorted(adata_pt.obs["sample"].astype(str).unique()):
    sample_mask = adata_pt.obs["sample"].astype(str).eq(sample).to_numpy()
    sample_curves, _ = fit_single_condition_curves(
        Y_genes[sample_mask], s[sample_mask], knots, grid, GAM_LAMBDA_GRID, ls["lam_idx"], support_pct=(0, 100)
    )
    species = str(adata_pt.obs.loc[sample_mask, "comparison_species"].iloc[0])
    for gene_i, gene in zip(top_idx, top_genes):
        for x, value in zip(grid, sample_curves[gene_i]):
            sample_curve_rows.append({"gene": gene, "sample": sample, "species": species, "pseudospace": x, "fitted_lognorm": value})
pd.DataFrame(sample_curve_rows).to_csv(RESULTS_DIR / "curves" / "top_gene_per_sample_curves.csv", index=False)

# %% [markdown]
# ## 6. Pathway/module curves and sample sensitivity

# %%
pathway_rows = []
if PATHWAY_LIBRARY_DIR.exists():
    for library_path in sorted(PATHWAY_LIBRARY_DIR.glob("*.json")):
        try:
            library = json.loads(library_path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            print(f"WARNING: cannot read pathway library {library_path.name}: {exc}")
            continue
        lookup = {str(g).upper(): i for i, g in enumerate(adata_pt.var_names)}
        for pathway, members in library.items():
            if isinstance(members, dict):
                members = members.get("genes", [])
            genes = [lookup[g.upper()] for g in members if str(g).upper() in lookup]
            if len(genes) >= 10:
                pathway_rows.append({"library": library_path.stem, "pathway": pathway, "genes": genes})
pathway_results = []
for row in pathway_rows[:500]:
    values = gene_matrix[:, row["genes"]].mean(axis=1)
    values = np.asarray(values).ravel() if sparse.issparse(values) else np.asarray(values)
    module_ls = run_level_shape(values[:, None], s, c, knots, grid, GAM_LAMBDA_GRID)
    pathway_results.append({
        "library": row["library"], "pathway": row["pathway"], "n_genes": len(row["genes"]),
        "species_effect_rms": float(module_ls["condition_effect_rms"][0]),
        "descriptive_level_effect_human_minus_mouse": float(module_ls["level_effect"][0]),
        "curve_spearman": float(module_ls["curve_spearman"][0]),
    })
pathway_df = pd.DataFrame(pathway_results).sort_values("species_effect_rms", ascending=False) if pathway_results else pd.DataFrame(columns=["library", "pathway", "n_genes", "species_effect_rms", "descriptive_level_effect_human_minus_mouse", "curve_spearman"])
pathway_df.to_csv(RESULTS_DIR / "curves" / "pathway_trajectory_comparison.csv", index=False)

sample_summary = (
    adata_pt.obs.groupby(["sample", "comparison_species", "region"], observed=True)
    .agg(n_pt_tubules=("sample", "size"), pseudospace_min=("shared_pseudospace", "min"), pseudospace_max=("shared_pseudospace", "max"))
    .reset_index()
)
sample_summary.to_csv(RESULTS_DIR / "diagnostics" / "sample_support.csv", index=False)

# %% [markdown]
# ## 6b. Three-axis concordance (secondary fidelity check)
#
# When centroid coordinates and Podocyte calls are present, distance to the nearest Podocyte
# provides a coarse physical cortex-to-medulla proxy. It is not used to fit expression curves;
# it only checks whether the molecular axis agrees with an independent spatial ordering.

# %%
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

concordance_rows = []
if {"x_centroid", "y_centroid"}.issubset(adata.obs.columns):
    xy = adata.obs[["x_centroid", "y_centroid"]].to_numpy(dtype=float)
    podocyte = adata.obs["broad_tubule_marker_call"].astype(str).eq("Podocyte").to_numpy()
    if podocyte.any() and np.isfinite(xy).all(axis=1).any():
        for sample in sorted(adata.obs["sample"].astype(str).unique()):
            sample_mask = adata.obs["sample"].astype(str).eq(sample).to_numpy()
            valid = sample_mask & np.isfinite(xy).all(axis=1)
            anchors = valid & podocyte
            if anchors.sum() < 2 or valid.sum() < 10:
                continue
            distances = np.full(adata.n_obs, np.nan)
            distances[valid] = cKDTree(xy[anchors]).query(xy[valid], k=1)[0]
            # Larger distance is the putative later/deeper direction; this is a diagnostic only.
            for axis_name, values in {
                "marker_axis": adata.obs["total_marker_axis"].to_numpy(dtype=float),
                "shared_pseudospace": adata.obs["shared_pseudospace"].to_numpy(dtype=float),
            }.items():
                finite = valid & np.isfinite(values) & np.isfinite(distances)
                rho = spearmanr(values[finite], distances[finite]).correlation if finite.sum() >= 10 else np.nan
                concordance_rows.append({
                    "sample": sample,
                    "axis": axis_name,
                    "physical_proxy": "distance_to_nearest_podocyte",
                    "n_tubules": int(finite.sum()),
                    "spearman": float(rho),
                })
pd.DataFrame(concordance_rows).to_csv(RESULTS_DIR / "diagnostics" / "three_axis_concordance.csv", index=False)

# %% [markdown]
# ## 7. Summary and publication caveats

# %%
notes = f"""# Human versus healthy-mouse run

Generated by `analysis/human_vs_healthy_mouse.py`.

- Matrix: {adata.n_obs:,} tubules × {adata.n_vars:,} shared mouse-ortholog genes.
- PT comparison: {adata_pt.n_obs:,} tubules; mouse controls = {', '.join(MOUSE_SAMPLES)}; human = {', '.join(HUMAN_SAMPLES)}.
- Integration: {integration_method}; batch key is sample when Harmony is available.
- Human donor count: 1. Mouse specimen count: 2.
- Human cortex/medulla and mouse mixed-kidney anatomy are not separable from species here.
- The nested GAM, curve differences, and pathway rankings are descriptive effect sizes; no
  confirmatory species p-values are reported.
- Re-run sensitivity checks after excluding low-capture `HUK1_MED1` and after restricting
  the mouse cohort to a matched anatomical region if such annotations become available.
- Human genes were mapped with the HCOP table plus any explicitly supplied override CSV;
  see `ortholog_map_used.csv` and `human_gene_mapping_report.csv`.
"""
(RESULTS_DIR / "analysis_notes.md").write_text(notes)
print("\nCross-species workflow complete.")
print(f"Results: {RESULTS_DIR}")
print(f"Top descriptive gene effects:\n{gene_results.head(10)[['gene', 'species_effect_rms', 'curve_spearman']].to_string(index=False)}")
