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
) -> pd.DataFrame:
    """Build a deterministic one-to-one human->mouse symbol map.

    HCOP can contain one-to-many and many-to-one assertions.  Those pairs are not
    used automatically: only symbols with exactly one partner on both sides are
    retained.  ``overrides`` is intended for reviewed kidney-specific exceptions
    and takes precedence after removing any conflicting HCOP rows.
    """
    pairs = table[["human_symbol", "mouse_symbol"]].drop_duplicates().copy()
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
        adatas[name] = adata
        print(f"{name}: {species}, {adata.n_obs:,} tubules x {adata.n_vars:,} shared genes")

    expected = mouse_samples | human_samples
    missing = sorted(expected - set(adatas))
    if missing:
        raise FileNotFoundError(f"Missing cross-species samples: {missing}")
    mapping_report = pd.concat(reports, ignore_index=True) if reports else pd.DataFrame()
    return adatas, mapping_report


def combine_cross_species(adatas: Mapping[str, ad.AnnData]) -> ad.AnnData:
    """Concatenate already-converted objects using the shared-gene intersection."""
    out = sc.concat(dict(adatas), join="inner", label="sample_from_concat", index_unique=None)
    out.obs["sample"] = out.obs["sample"].astype(str)
    out.obs["species"] = out.obs["species"].astype(str)
    if not out.obs_names.is_unique:
        out.obs_names_make_unique()
    if not out.var_names.is_unique:
        out.var_names_make_unique()
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
