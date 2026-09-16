"""Cross-species loading and ortholog-space helpers.

The human-versus-healthy-mouse workflow deliberately converts human expression into
mouse-symbol space before integration.  The conversion is explicit and auditable:
the returned mapping report records every source gene as mapped, ambiguous, or
unmapped, and the expression matrix is aggregated with a sparse transformation
matrix rather than silently dropping duplicate columns.
"""
from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from .io_qc import calculate_qc_metrics, sample_name_from_path


def _symbol(value: object) -> str:
    value = str(value).strip()
    return "" if value.lower() in {"", "nan", "none", "-"} else value


def read_ortholog_table(path: Path) -> pd.DataFrame:
    """Read an HCOP-style tab-separated table and validate its symbol columns."""
    table = pd.read_csv(path, sep="\t", compression="infer", dtype=str)
    required = {"human_symbol", "mouse_symbol"}
    missing = sorted(required - set(table.columns))
    if missing:
        raise ValueError(f"Ortholog table {path} is missing columns: {missing}")
    table = table.copy()
    table["human_symbol"] = table["human_symbol"].map(_symbol)
    table["mouse_symbol"] = table["mouse_symbol"].map(_symbol)
    table = table.loc[(table["human_symbol"] != "") & (table["mouse_symbol"] != "")].copy()
    if table.empty:
        raise ValueError(f"Ortholog table {path} contains no usable symbol pairs")
    return table


def build_one_to_one_ortholog_map(
    table: pd.DataFrame,
    overrides: Mapping[str, str] | None = None,
    min_support: int = 3,
) -> pd.DataFrame:
    """Build a deterministic one-to-one human->mouse symbol map.

    HCOP can contain one-to-many and many-to-one assertions.  Those pairs are not
    used automatically: only symbols with exactly one partner on both sides are
    retained.  ``overrides`` is intended for reviewed kidney-specific exceptions
    and takes precedence after removing any conflicting HCOP rows.
    When a ``support`` column is present (standard HCOP tables), consensus support-weighted
    reciprocal best matching is applied:
      1. Low-confidence noise pairs supported by fewer than ``min_support`` databases are pruned.
      2. A candidate score is computed from database support count with a tie-breaking bonus for
         identical gene symbols.
      3. Only mutual reciprocal best pairs (the top mouse partner for the human gene is also
         the top human partner for that mouse gene) are retained.

    When ``support`` is absent (e.g. synthetic test tables), an unweighted bidirectional
    uniqueness filter is applied as fallback.

    ``overrides`` allows reviewed kidney-specific exceptions and takes precedence after
    removing any conflicting HCOP rows.
    """
    pairs = table[["human_symbol", "mouse_symbol"]].drop_duplicates().copy()
    pairs = table.copy()
    pairs["human_symbol"] = pairs["human_symbol"].map(_symbol)
    pairs["mouse_symbol"] = pairs["mouse_symbol"].map(_symbol)
    pairs = pairs.loc[(pairs["human_symbol"] != "") & (pairs["mouse_symbol"] != "")].copy()
    pairs["human_key"] = pairs["human_symbol"].str.upper()
    pairs["mouse_key"] = pairs["mouse_symbol"].str.upper()
    h_counts = pairs.groupby("human_key")["mouse_key"].nunique()
    m_counts = pairs.groupby("mouse_key")["human_key"].nunique()
    pairs["mapping_status"] = np.where(
        (pairs["human_key"].map(h_counts) == 1) & (pairs["mouse_key"].map(m_counts) == 1),
        "hcop_one_to_one",
        "ambiguous_hcop",
    )
    clean = pairs.loc[pairs["mapping_status"].eq("hcop_one_to_one"),
                      ["human_symbol", "mouse_symbol", "mapping_status"]].copy()

    if "support" in pairs.columns:
        pairs["n_support"] = pairs["support"].apply(
            lambda x: len(str(x).split(",")) if pd.notna(x) and str(x) != "-" else 1
        )
        pairs["same_name"] = (pairs["human_key"] == pairs["mouse_key"]).astype(int)
        pairs["score"] = pairs["n_support"] * 10 + pairs["same_name"]

        filtered = pairs[pairs["n_support"] >= min_support].copy()
        best_h = filtered.sort_values(["human_key", "score"], ascending=[True, False]).drop_duplicates(subset=["human_key"], keep="first")
        best_m = filtered.sort_values(["mouse_key", "score"], ascending=[True, False]).drop_duplicates(subset=["mouse_key"], keep="first")

        clean = pd.merge(
            best_h[["human_symbol", "mouse_symbol", "human_key", "mouse_key"]],
            best_m[["human_symbol", "mouse_symbol", "human_key", "mouse_key"]],
            on=["human_symbol", "mouse_symbol", "human_key", "mouse_key"],
        )
        clean["mapping_status"] = "hcop_one_to_one"
    else:
        pairs = pairs[["human_symbol", "mouse_symbol", "human_key", "mouse_key"]].drop_duplicates()
        h_counts = pairs.groupby("human_key")["mouse_key"].nunique()
        m_counts = pairs.groupby("mouse_key")["human_key"].nunique()
        pairs["mapping_status"] = np.where(
            (pairs["human_key"].map(h_counts) == 1) & (pairs["mouse_key"].map(m_counts) == 1),
            "hcop_one_to_one",
            "ambiguous_hcop",
        )
        clean = pairs.loc[pairs["mapping_status"].eq("hcop_one_to_one"),
                          ["human_symbol", "mouse_symbol", "mapping_status"]].copy()

    if overrides:
        override_rows = []
        for human, mouse in overrides.items():
            human, mouse = _symbol(human), _symbol(mouse)
            if human and mouse:
                override_rows.append({
                    "human_symbol": human,
                    "mouse_symbol": mouse,
                    "mapping_status": "curated_override",
                })
        if override_rows:
            override_df = pd.DataFrame(override_rows)
            h_keys = override_df["human_symbol"].str.upper()
            m_keys = override_df["mouse_symbol"].str.upper()
            clean = clean.loc[
                ~clean["human_symbol"].str.upper().isin(h_keys)
                & ~clean["mouse_symbol"].str.upper().isin(m_keys)
            ]
            clean = pd.concat([clean, override_df], ignore_index=True)

    clean = clean[["human_symbol", "mouse_symbol", "mapping_status"]].copy()
    if clean["human_symbol"].str.upper().duplicated().any() or clean["mouse_symbol"].str.upper().duplicated().any():
        raise ValueError("Ortholog overrides create a non-bijective mapping")
    return clean.sort_values("human_symbol", key=lambda s: s.str.upper()).reset_index(drop=True)


def map_human_to_mouse_space(
    adata: ad.AnnData,
    ortholog_map: pd.DataFrame,
    ambiguous_human_symbols: Sequence[str] | None = None,
) -> tuple[ad.AnnData, pd.DataFrame]:
    """Convert an AnnData object to mouse-symbol space and return a gene report."""
    required = {"human_symbol", "mouse_symbol", "mapping_status"}
    missing = sorted(required - set(ortholog_map.columns))
    if missing:
        raise ValueError(f"ortholog_map is missing columns: {missing}")

    human_to_mouse = {
        str(h).upper(): str(m)
        for h, m in zip(ortholog_map["human_symbol"], ortholog_map["mouse_symbol"])
    }
    target_symbols = list(dict.fromkeys(ortholog_map["mouse_symbol"].astype(str)))
    target_index = {symbol.upper(): i for i, symbol in enumerate(target_symbols)}
    source_symbols = pd.Index(adata.var_names.astype(str))

    report = pd.DataFrame({"source_human_symbol": source_symbols})
    report["mouse_symbol"] = report["source_human_symbol"].str.upper().map(human_to_mouse)
    ambiguous = {str(symbol).upper() for symbol in (ambiguous_human_symbols or ())}
    report["mapping_status"] = np.select(
        [report["mouse_symbol"].notna(), report["source_human_symbol"].str.upper().isin(ambiguous)],
        ["mapped", "ambiguous_hcop"],
        default="unmapped",
    )
    report["source_column_index"] = np.arange(adata.n_vars, dtype=int)
    mapped_mask = report["mouse_symbol"].notna().to_numpy()
    source_idx = np.flatnonzero(mapped_mask)
    target_idx = np.array([target_index[str(x).upper()] for x in report.loc[mapped_mask, "mouse_symbol"]], dtype=int)
    transform = sparse.csr_matrix(
        (np.ones(len(source_idx), dtype=float), (source_idx, target_idx)),
        shape=(adata.n_vars, len(target_symbols)),
    )

    def convert(matrix):
        if sparse.issparse(matrix):
            return matrix.tocsr() @ transform
        return np.asarray(matrix) @ transform.toarray()

    out = ad.AnnData(
        X=convert(adata.X),
        obs=adata.obs.copy(),
        var=pd.DataFrame(index=pd.Index(target_symbols, name="mouse_symbol")),
    )
    for layer_name, layer in adata.layers.items():
        if getattr(layer, "shape", None) == adata.shape:
            out.layers[layer_name] = convert(layer)
    out.uns["cross_species_mapping"] = {
        "n_source_human_genes": int(adata.n_vars),
        "n_mapped_source_genes": int(mapped_mask.sum()),
        "n_target_mouse_genes": int(len(target_symbols)),
    }
    out.var["source_human_symbol"] = [
        next((h for h, m in human_to_mouse.items() if m.upper() == target.upper()), "")
        for target in target_symbols
    ]
    # A target column exists for every accepted ortholog, but it is only a MEASUREMENT when some
    # source gene actually fed it. The rest are structural zeros that must not be read as "this gene
    # has no expression here" - the distinction the shared space has to carry forward.
    measured_targets = np.zeros(len(target_symbols), dtype=bool)
    if len(target_idx):
        measured_targets[target_idx] = True
    out.var["measured_in_source_input"] = measured_targets
    out.uns["cross_species_mapping"]["n_targets_measured_in_source_input"] = int(
        measured_targets.sum())
    out.uns["cross_species_mapping"]["n_targets_without_source_feature"] = int(
        (~measured_targets).sum())
    return out, report


def load_cross_species_samples(
    data_dir: Path,
    files: Sequence[str | Path],
    ortholog_map: pd.DataFrame,
    mouse_samples: Sequence[str] = ("Ctrl1A2", "Ctrl1A4"),
    human_samples: Sequence[str] = ("HUK1_COR1", "HUK1_MED1"),
    human_regions: Mapping[str, str] | None = None,
    ambiguous_human_symbols: Sequence[str] | None = None,
) -> tuple[dict[str, ad.AnnData], pd.DataFrame]:
    """Load two mouse controls and human samples into a shared gene space."""
    mouse_samples, human_samples = set(mouse_samples), set(human_samples)
    human_regions = dict(human_regions or {})
    adatas: dict[str, ad.AnnData] = {}
    reports: list[pd.DataFrame] = []
    for file in files:
        path = Path(file)
        name = sample_name_from_path(path)
        if name not in mouse_samples | human_samples:
            continue
        native = sc.read_h5ad(data_dir / path if not path.is_absolute() else path)
        if "counts" in native.layers:
            native.X = native.layers["counts"].copy()
        species = "mouse" if name in mouse_samples else "human"
        native = calculate_qc_metrics(native, species=species)
        native.obs["original_obs_name"] = native.obs_names.astype(str)
        native.obs_names = pd.Index([f"{name}_{obs}" for obs in native.obs_names])
        native.obs["sample"] = name
        native.obs["species"] = species
        native.obs["comparison_species"] = species
        native.obs["condition"] = "Control" if species == "mouse" else "Healthy"
        native.obs["region"] = ("mixed_kidney" if species == "mouse"
                                 else human_regions.get(name, "unspecified"))

        if species == "human":
            converted, report = map_human_to_mouse_space(
                native, ortholog_map, ambiguous_human_symbols=ambiguous_human_symbols
            )
            converted.obs = native.obs.copy()
            report.insert(0, "sample", name)
            reports.append(report)
            adata = converted
        else:
            adata = native
            report = pd.DataFrame({
                "sample": name,
                "source_human_symbol": pd.Series(dtype=str),
                "mouse_symbol": pd.Series(dtype=str),
                "mapping_status": pd.Series(dtype=str),
            })
        if species == "mouse":
            # Mouse objects are already in the target space, so every column is a measurement.
            adata.var["measured_in_source_input"] = True
            adata.uns["cross_species_mapping"] = {
                "n_targets_measured_in_source_input": int(adata.n_vars),
                "n_targets_without_source_feature": 0,
            }
        adatas[name] = adata
        print(f"{name}: {species}, {adata.n_obs:,} tubules x {adata.n_vars:,} shared genes "
              f"({int(adata.var['measured_in_source_input'].sum()):,} measured in the source input)")

    expected = mouse_samples | human_samples
    missing = sorted(expected - set(adatas))
    if missing:
        raise FileNotFoundError(f"Missing cross-species samples: {missing}")
    mapping_report = pd.concat(reports, ignore_index=True) if reports else pd.DataFrame()
    return adatas, mapping_report


def combine_cross_species(adatas: Mapping[str, ad.AnnData], *,
                          require_measured_in_both: bool = True) -> ad.AnnData:
    """Concatenate already-converted objects over the genes measured in every input.

    An ortholog mapping guarantees a unique counterpart, not that both input matrices contain that
    gene. Converting human symbols creates a target column for every accepted pair and fills zeros
    where the human input has no such feature, so a plain symbol intersection keeps those structural
    zeros and turns "not supplied in this measurement" into apparent evidence of no expression.
    Genes measured in only one species are therefore dropped by default, with the dropped symbols
    recorded in ``uns['cross_species_availability']``. ``require_measured_in_both=False`` restores
    the older symbol-only intersection for comparison.
    """
    measured_sets = []
    for name, adata in adatas.items():
        flags = adata.var.get("measured_in_source_input")
        if flags is None:
            measured_sets.append(set(map(str, adata.var_names)))
        else:
            measured_sets.append(set(map(str, adata.var_names[np.asarray(flags, dtype=bool)])))
    shared = set.intersection(*measured_sets) if measured_sets else set()
    # What each object holds but cannot contribute to the joint analysis: a gene missing from the
    # shared set was not measured in some input, whether or not this particular object has it.
    per_species_unmeasured = {
        str(name): sorted(set(map(str, adata.var_names)) - shared)[:500]
        for name, adata in adatas.items()
    }
    all_symbols = set().union(*(set(map(str, adata.var_names)) for adata in adatas.values()))
    dropped_union = sorted(all_symbols - shared)
    if require_measured_in_both:
        adatas = {
            name: adata[:, [gene for gene in adata.var_names if str(gene) in shared]]
            for name, adata in adatas.items()
        }
    availability = {
        "n_shared_genes": len(shared),
        "n_dropped_not_measured_in_every_input": len(dropped_union),
        "dropped_symbols": dropped_union[:500],
        "unmeasured_symbols_by_sample": per_species_unmeasured,
        "rule": ("measured_in_source_input in every input"
                 if require_measured_in_both else
                 "measured_in_source_input flagged per gene; the object keeps every gene"),
        "genes_removed_from_the_object": bool(require_measured_in_both),
    }
    out = sc.concat(dict(adatas), join="inner", label="sample_from_concat", index_unique=None)
    out.obs["sample"] = out.obs["sample"].astype(str)
    out.obs["species"] = out.obs["species"].astype(str)
    if not out.obs_names.is_unique:
        out.obs_names_make_unique()
    if not out.var_names.is_unique:
        out.var_names_make_unique()
    # Per gene, in BOTH modes: True only where every input actually measured the gene. The notebooks
    # gate their analysis gene set on this, so the availability rule can be applied without changing
    # which genes feed the HVG/PCA/Harmony step (a selection step that is sensitive to the gene set,
    # and therefore to the clustering and its hand-reviewed labels).
    out.var["measured_in_both_inputs"] = [str(gene) in shared for gene in out.var_names]
    out.uns["cross_species_availability"] = availability
    return out


def prepare_shared_expression(
    adata: ad.AnnData,
    min_genes: int = 100,
    min_gene_total_counts: int = 20,
) -> ad.AnnData:
    """Filter, normalize, and annotate the shared mouse-symbol matrix."""
    out = adata.copy()
    if "counts" in out.layers:
        out.X = out.layers["counts"].copy()
    sc.pp.filter_cells(out, min_genes=min_genes)
    sc.pp.filter_genes(out, min_counts=min_gene_total_counts)
    out.layers["counts"] = out.X.copy()
    sc.pp.normalize_total(out, target_sum=1e4)
    sc.pp.log1p(out)
    out.layers["lognorm"] = out.X.copy()
    out.var["mt"] = out.var_names.astype(str).str.startswith("mt-")
    out.var["ribo"] = out.var_names.astype(str).str.startswith(("Rpl", "Rps", "Mrpl", "Mrps"))
    out.var["exclude_from_harmony_hvg"] = out.var["mt"] | out.var["ribo"]
    return out


def plot_cross_species_marker_alignment(
    adata: ad.AnnData,
    output_path: Path | str | None = None,
    markers: Sequence[tuple[str, str, str]] | None = None,
    basis: str = "umap",
    species_key: str = "comparison_species",
    segment_key: str = "segment_class",
    dpi: int = 200,
):
    """Plot side-by-side (Human vs Mouse) canonical marker alignment on embedding coordinates.

    Row 0 shows the categorical segment annotation for biological context.
    Subsequent rows display canonical nephron segment markers with synchronized
    color scaling, sorted so high-expressing structures are rendered on top.
    """
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    if markers is None:
        markers = [
            ("Slc5a2", "PT-S1", "Slc5a2 (SGLT2)"),
            ("Slc22a6", "PT-S2", "Slc22a6 (OAT1)"),
            ("Slc22a7", "PT-S3", "Slc22a7 (OAT2)"),
            ("Slc12a1", "mTAL", "Slc12a1 (NKCC2)"),
            ("Slc12a3", "DCT", "Slc12a3 (NCC)"),
            ("Aqp2", "Collecting Duct", "Aqp2 (Aquaporin-2)"),
        ]

    rep_key = f"X_{basis}" if not basis.startswith("X_") else basis
    if rep_key not in adata.obsm:
        raise KeyError(f"Embedding coordinates {rep_key!r} not found in adata.obsm")

    coords = np.asarray(adata.obsm[rep_key])
    expr_mat = adata.layers["lognorm"] if "lognorm" in adata.layers else adata.X
    species = adata.obs[species_key].to_numpy()
    is_human = (species == "human")
    is_mouse = (species == "mouse")

    segments = ["PT-S1", "PT-S2", "PT-S3", "mTAL", "DCT2", "CNT_CD", "IMCD"]
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2"]
    seg_color_map = dict(zip(segments, colors))
    canonical_palette = {
        "PT": "#1f77b4",
        "PT-S1": "#4C9BD3",
        "PT-S2": "#5B8E55",
        "PT-S3": "#D6B48A",
        "ATL": "#26A69A",
        "DTL": "#7E57C2",
        "mTAL": "#F58518",
        "cTAL": "#E8A33D",
        "AL": "#F58518",
        "DCT": "#D81B60",
        "DCT1": "#D81B60",
        "DCT2": "#EC6EA5",
        "CNT": "#76B7B2",
        "CCD": "#4E79A7",
        "OMCD": "#9C755F",
        "IMCD": "#593C8F",
        "CNT_CD": "#8D6E63",
    }
    seg_series = adata.obs[segment_key].astype(str).to_numpy()
    unique_present = [s for s in pd.unique(adata.obs[segment_key]) if str(s) not in ("nan", "None", "")]
    canonical_order = list(canonical_palette.keys())
    segments = [s for s in canonical_order if s in unique_present] + [s for s in unique_present if s not in canonical_palette]

    default_colors = plt.cm.tab20.colors
    seg_color_map = {}
    for i, seg in enumerate(segments):
        if seg in canonical_palette:
            seg_color_map[seg] = canonical_palette[seg]
        else:
            seg_color_map[seg] = default_colors[i % len(default_colors)]

    n_rows = 1 + len(markers)
    fig = plt.figure(figsize=(11, 3.4 * n_rows))
    gs = fig.add_gridspec(
        n_rows, 2, width_ratios=[1, 1], wspace=0.08, hspace=0.18,
        left=0.05, right=0.88, top=0.96, bottom=0.02
    )

    # Row 0: Segment Annotations
    ax_h0 = fig.add_subplot(gs[0, 0])
    ax_m0 = fig.add_subplot(gs[0, 1])

    for ax, mask, sp_name in [
        (ax_h0, is_human, f"Human ({is_human.sum():,} structures)"),
        (ax_m0, is_mouse, f"Mouse ({is_mouse.sum():,} structures)"),
    ]:
        ax.scatter(coords[~mask, 0], coords[~mask, 1], c="#ececec", s=1.0, alpha=0.3, rasterized=True)
        for seg in segments:
            seg_mask = mask & (seg_series == seg)
            if seg_mask.any():
                ax.scatter(
                    coords[seg_mask, 0], coords[seg_mask, 1],
                    c=seg_color_map[seg], s=3.0, alpha=0.85, label=seg, rasterized=True
                )
        ax.set_title(f"{sp_name}\nSegment Annotation", fontsize=11, fontweight="bold", pad=8)
        ax.axis("off")

    legend_elements = [
        Line2D([0], [0], marker="o", color="w", label=s, markerfacecolor=seg_color_map[s], markersize=8)
        for s in segments
    ]
    fig.legend(
        handles=legend_elements, loc="upper left", bbox_to_anchor=(0.89, 0.96),
        title="Segment Class", frameon=False, fontsize=9, title_fontsize=10
    )

    cmap = plt.cm.inferno

    for row_idx, (gene, seg_label, title) in enumerate(markers, start=1):
        if gene not in adata.var_names:
            continue
        g_idx = adata.var_names.get_loc(gene)
        expr = expr_mat[:, g_idx]
        if hasattr(expr, "toarray"):
            expr = expr.toarray().ravel()
        else:
            expr = np.asarray(expr).ravel()

        pos_expr = expr[expr > 0]
        vmax = float(np.percentile(pos_expr, 99)) if len(pos_expr) > 0 else 1.0

        ax_h = fig.add_subplot(gs[row_idx, 0])
        ax_m = fig.add_subplot(gs[row_idx, 1])

        for ax, mask, sp_name in [(ax_h, is_human, "Human"), (ax_m, is_mouse, "Mouse")]:
            ax.scatter(coords[~mask, 0], coords[~mask, 1], c="#f0f0f0", s=1.0, alpha=0.3, rasterized=True)

            sp_coords = coords[mask]
            sp_expr = expr[mask]

            zero_mask = (sp_expr == 0)
            ax.scatter(sp_coords[zero_mask, 0], sp_coords[zero_mask, 1], c="#cfd4d8", s=2.0, alpha=0.4, rasterized=True)

            pos_sub = ~zero_mask
            if pos_sub.any():
                sort_order = np.argsort(sp_expr[pos_sub])
                sc_plot = ax.scatter(
                    sp_coords[pos_sub][sort_order, 0], sp_coords[pos_sub][sort_order, 1],
                    c=sp_expr[pos_sub][sort_order], cmap=cmap, vmin=0, vmax=vmax,
                    s=3.5, alpha=0.9, rasterized=True
                )

            pct_pos = (sp_expr > 0).mean() * 100
            ax.set_title(f"{sp_name}: {title} ({pct_pos:.1f}% > 0)", fontsize=10.5, fontweight="bold", pad=6)
            ax.axis("off")

        cbar_ax = fig.add_axes([
            0.90, gs[row_idx, 1].get_position(fig).y0 + 0.015,
            0.015, gs[row_idx, 1].get_position(fig).height * 0.7
        ])
        cbar = fig.colorbar(sc_plot, cax=cbar_ax)
        cbar.ax.tick_params(labelsize=8)
        cbar.set_label("log-norm counts", fontsize=8)

    if output_path is not None:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
        print(f"Saved side-by-side marker alignment to {output_path}")
    return fig


def plot_pre_filtering_tubule_qc(
    adata: ad.AnnData,
    output_path: Path | str | None = None,
    basis: str = "umap",
    leiden_key: str = "leiden_coarse",
    coarse_key: str = "coarse_class",
    segment_key: str = "segment_class",
    species_key: str = "comparison_species",
    sample_key: str = "sample",
    remove_classes: Sequence[str] = ("Unassigned",),
    stroma_gene: str = "Vim",
    doublet_gene: str = "Cryab",
    dpi: int = 200,
):
    """Plot a 6-panel pre-filtering diagnostic on the Pass-1 manifold.

    Audits which segmented structures are kept vs removed before filtering for Pass 2:
    - Panel A: UMAP of Kept vs Removed structures.
    - Panel B: Leiden clusters (0-8) with kept/removed status labeled.
    - Panel C: Species distribution (Human vs Mouse).
    - Panel D: Non-tubular stroma / glomerular signature (e.g. Vim).
    - Panel E: Mixed doublet / stress bridge signature (e.g. Cryab).
    - Panel F: Retention rate bar chart by sample.
    """
    import matplotlib.pyplot as plt

    rep_key = f"X_{basis}" if not basis.startswith("X_") else basis
    if rep_key not in adata.obsm:
        raise KeyError(f"Embedding coordinates {rep_key!r} not found in adata.obsm")

    coords = np.asarray(adata.obsm[rep_key])
    expr_mat = adata.layers["lognorm"] if "lognorm" in adata.layers else adata.X

    is_kept = ~adata.obs[coarse_key].isin(remove_classes).to_numpy()
    leiden = adata.obs[leiden_key].astype(str).to_numpy()
    species = adata.obs[species_key].astype(str).to_numpy()
    sample = adata.obs[sample_key].astype(str).to_numpy()

    fig = plt.figure(figsize=(19, 12))
    gs = fig.add_gridspec(2, 3, wspace=0.22, hspace=0.25, left=0.04, right=0.96, top=0.93, bottom=0.06)

    # 1. Kept vs Removed
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.scatter(coords[is_kept, 0], coords[is_kept, 1], c="#2b7bba", s=2.5, alpha=0.7,
                label=f"Kept: Nephron Tubules (n={is_kept.sum():,})", rasterized=True)
    ax1.scatter(coords[~is_kept, 0], coords[~is_kept, 1], c="#e6550d", s=3.5, alpha=0.85,
                label=f"Removed: Non-tubule / Mixed (n={(~is_kept).sum():,})", rasterized=True)
    ax1.set_title("A. Pass-1 Filtering Decision: Kept vs Removed", fontsize=12, fontweight="bold", pad=8)
    ax1.legend(loc="lower left", frameon=True, facecolor="white", framealpha=0.92, fontsize=9.5)
    ax1.axis("off")

    # 2. Leiden Clusters with removed clusters highlighted
    ax2 = fig.add_subplot(gs[0, 1])
    unique_clusters = sorted(np.unique(leiden), key=lambda x: int(x) if x.isdigit() else x)
    palette = ["#1f77b4", "#33a02c", "#6baed6", "#9467bd", "#ff7f0e", "#e31a1c", "#ff0055", "#7f3b08", "#e377c2", "#17becf", "#bcbd22"]

    import matplotlib.patheffects as PathEffects

    cluster_labels = {}
    for c in unique_clusters:
        c_mask = (leiden == c)
        seg = adata.obs.loc[c_mask, segment_key].iloc[0] if segment_key in adata.obs else ""
        c_removed = bool(adata.obs.loc[c_mask, coarse_key].isin(remove_classes).mean() > 0.5)
        status = "REMOVED" if c_removed else "kept"
        cluster_labels[c] = f"{c}: {seg} ({status})" if seg else f"{c} ({status})"

    for i, c in enumerate(unique_clusters):
        mask = (leiden == c)
        color = palette[i % len(palette)]
        ax2.scatter(coords[mask, 0], coords[mask, 1], c=color, s=2.5, alpha=0.75,
                    label=cluster_labels.get(c, c), rasterized=True)
        if mask.any():
            x_c, y_c = np.median(coords[mask], axis=0)
            txt = ax2.text(x_c, y_c, str(c), ha="center", va="center", fontsize=9.5, fontweight="bold", color="black")
            txt.set_path_effects([PathEffects.withStroke(linewidth=3.0, foreground="white")])

    ax2.set_title("B. Pass-1 Leiden Clusters (0–8)", fontsize=12, fontweight="bold", pad=8)
    ax2.legend(loc="lower left", frameon=True, facecolor="white", framealpha=0.92, fontsize=8, markerscale=2)
    ax2.axis("off")

    # 3. Species breakdown on Pass 1
    ax3 = fig.add_subplot(gs[0, 2])
    sp_colors = {"human": "#1f77b4", "mouse": "#ff7f0e"}
    for sp in sorted(np.unique(species)):
        mask = (species == sp)
        ax3.scatter(coords[mask, 0], coords[mask, 1], c=sp_colors.get(sp, "#333333"),
                    s=2.5, alpha=0.7, label=f"{sp.capitalize()} (n={mask.sum():,})", rasterized=True)
    ax3.set_title("C. Species Distribution (Pass 1)", fontsize=12, fontweight="bold", pad=8)
    ax3.legend(loc="lower left", frameon=True, facecolor="white", framealpha=0.92, fontsize=10, markerscale=2)
    ax3.axis("off")

    # Helper for gene expression panels
    def _plot_gene(ax, gene_name, title):
        if gene_name not in adata.var_names:
            ax.text(0.5, 0.5, f"Gene '{gene_name}' not found", ha="center", va="center", fontsize=11)
            ax.axis("off")
            return
        g_idx = adata.var_names.get_loc(gene_name)
        expr = expr_mat[:, g_idx]
        if hasattr(expr, "toarray"):
            expr = expr.toarray().ravel()
        else:
            expr = np.asarray(expr).ravel()
        pos = expr[expr > 0]
        vmax = float(np.percentile(pos, 99)) if len(pos) > 0 else 1.0
        sort_order = np.argsort(expr)
        sc = ax.scatter(
            coords[sort_order, 0], coords[sort_order, 1],
            c=expr[sort_order], cmap="magma", s=2.5, alpha=0.85,
            vmin=0, vmax=vmax, rasterized=True
        )
        ax.set_title(title, fontsize=12, fontweight="bold", pad=8)
        cb = plt.colorbar(sc, ax=ax, fraction=0.035, pad=0.03)
        cb.set_label("log-norm expression", fontsize=8.5)
        ax.axis("off")

    # 4. Cluster 7 Marker: Vim (Stroma / Glomerular / Interstitial)
    ax4 = fig.add_subplot(gs[1, 0])
    _plot_gene(ax4, stroma_gene, f"D. Cluster 7 Signature: {stroma_gene} (Stroma / Podocyte / Interstitial)")

    # 5. Cluster 6 Marker: Cryab (Mixed / Stress / Edge Program)
    ax5 = fig.add_subplot(gs[1, 1])
    _plot_gene(ax5, doublet_gene, f"E. Cluster 6 Signature: {doublet_gene} (Mixed Doublet / Stress Bridge)")

    # 6. Breakdown table/bar chart by sample
    ax6 = fig.add_subplot(gs[1, 2])
    ct = pd.crosstab(
        adata.obs[sample_key],
        pd.Series(np.where(is_kept, "Kept", "Removed"), index=adata.obs_names)
    )
    for col in ["Kept", "Removed"]:
        if col not in ct.columns:
            ct[col] = 0
    ct["Total"] = ct["Kept"] + ct["Removed"]
    ct["% Kept"] = (ct["Kept"] / ct["Total"] * 100).round(1)

    samples = ct.index.tolist()
    y_pos = np.arange(len(samples))
    bar_height = 0.55

    ax6.barh(y_pos, ct["Kept"], height=bar_height, color="#2b7bba", label="Kept Tubules", edgecolor="black", linewidth=0.5)
    ax6.barh(y_pos, ct["Removed"], left=ct["Kept"], height=bar_height, color="#e6550d", label="Removed", edgecolor="black", linewidth=0.5)

    for i, (k, r, pct) in enumerate(zip(ct["Kept"], ct["Removed"], ct["% Kept"])):
        ax6.text(k / 2, i, f"{k:,}", va="center", ha="center", color="white", fontweight="bold", fontsize=9)
        if r > 300:
            ax6.text(k + r / 2, i, f"{r:,}", va="center", ha="center", color="white", fontweight="bold", fontsize=8.5)
        ax6.text(k + r + 200, i, f"{pct}% kept", va="center", ha="left", fontweight="bold", fontsize=9, color="#333333")

    ax6.set_yticks(y_pos)
    ax6.set_yticklabels(samples, fontsize=9.5, fontweight="bold")
    ax6.set_xlabel("Number of Segmented Structures", fontsize=10)
    ax6.set_xlim(0, ct["Total"].max() * 1.30)
    ax6.set_title("F. Tubule Retention by Sample", fontsize=12, fontweight="bold", pad=8)
    ax6.legend(loc="lower right", frameon=True, fontsize=9.5)
    ax6.spines["top"].set_visible(False)
    ax6.spines["right"].set_visible(False)

    plt.suptitle("Pre-Filtering Tubule QC & Filtering Audit (Pass-1 Harmony Manifold)", fontsize=15, fontweight="bold", y=0.98)

    if output_path is not None:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
        print(f"Saved pre-filtering tubule QC to {output_path}")
    return fig
