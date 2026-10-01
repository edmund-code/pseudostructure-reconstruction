"""Load current mouse structure-level counts with validated segmentation annotations.

This loader deliberately stops on stale ``feature_index`` joins.  The feature index was
assigned by enumerating a particular GeoJSON, so a range check alone is not sufficient.
"""
from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from shapely.geometry import shape

from .io_qc import combine_adatas, load_and_process_mouse_samples
from .vocabulary import coarse_for

MOUSE_SAMPLES = ("Ctrl1A2", "Ctrl1A4", "IR2A2", "IR2A4")
SPOT_QUANTILE = 0.05
CENTROID_TOLERANCE_PX = 1.0  # same contract as notebook 01 verify_segmentation_matches_obs


def _counts_layer(adata: ad.AnnData, sample: str):
    """Return validated raw counts, refusing transformed matrices as count input."""
    if not adata.var_names.is_unique:
        raise ValueError(f"{sample}: duplicate gene names would be silently renamed during combine")
    matrix = adata.layers["counts"] if "counts" in adata.layers else adata.X
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError(f"{sample}: raw counts must be finite and nonnegative")
    if not np.equal(values, np.floor(values)).all():
        raise ValueError(
            f"{sample}: count source is not integer-valued; refusing to treat normalized/logged X as raw counts"
        )
    row_sums = np.asarray(matrix.sum(axis=1)).ravel()
    if not np.isfinite(row_sums).all() or (row_sums <= 0).any():
        raise ValueError(f"{sample}: every structure must have a positive finite raw-count sum")
    return matrix


def _read_and_validate_segmentation(data_root: Path, sample: str, obs: pd.DataFrame):
    path = data_root / f"{sample}_kept_tubules_labeled_fine.geojson"
    if not path.is_file():
        raise FileNotFoundError(f"{sample}: expected current fine segmentation at {path}")
    payload = json.loads(path.read_text())
    features = payload.get("features", [])
    required = {"feature_index", "x_centroid", "y_centroid"}
    missing = sorted(required - set(obs.columns))
    if missing:
        raise ValueError(f"{sample}: expression matrix lacks segmentation join fields: {missing}")
    indices = pd.to_numeric(obs.feature_index, errors="coerce").to_numpy()
    if not np.isfinite(indices).all() or not np.equal(indices, np.floor(indices)).all():
        raise ValueError(f"{sample}: feature_index contains missing or non-integer values")
    indices = indices.astype(int)
    if len(np.unique(indices)) != len(indices):
        raise ValueError(f"{sample}: feature_index is not unique")
    if ((indices < 0) | (indices >= len(features))).any():
        raise ValueError(f"{sample}: feature_index outside current GeoJSON feature range")

    centroids = obs[["x_centroid", "y_centroid"]].to_numpy(dtype=float)
    if not np.isfinite(centroids).all():
        raise ValueError(f"{sample}: source x_centroid/y_centroid contains non-finite values")
    labels, max_dev = [], 0.0
    for row, fi in zip(obs.itertuples(), indices):
        feature = features[fi]
        props = feature.get("properties") or {}
        label = props.get("segment_class")
        if not label:
            raise ValueError(f"{sample}: feature {fi} has no segment_class property")
        geom = feature.get("geometry")
        if not geom:
            raise ValueError(f"{sample}: feature {fi} has no geometry")
        polygon = shape(geom)
        if not polygon.is_valid:
            polygon = polygon.buffer(0)
        if polygon.is_empty:
            raise ValueError(f"{sample}: feature {fi} has empty geometry")
        centroid = polygon.centroid
        dev = float(np.hypot(centroid.x - float(row.x_centroid), centroid.y - float(row.y_centroid)))
        if not np.isfinite([centroid.x, centroid.y, dev]).all():
            raise ValueError(f"{sample}: feature {fi} has non-finite geometry/source centroid")
        max_dev = max(max_dev, dev)
        if dev > CENTROID_TOLERANCE_PX:
            raise ValueError(
                f"{sample}: feature_index/centroid mismatch at feature {fi}: {dev:.3f}px "
                f"(tolerance {CENTROID_TOLERANCE_PX}px); refusing stale segmentation join"
            )
        labels.append(str(label))
    return np.asarray(labels, dtype=object), {
        "path": path.name,
        "feature_count": len(features),
        "observations_checked": len(obs),
        "max_centroid_deviation_px": max_dev,
        "centroid_tolerance_px": CENTROID_TOLERANCE_PX,
        "size_bytes": path.stat().st_size,
        "mtime_ns": path.stat().st_mtime_ns,
    }


def load_repeated_mouse_inputs(data_root: str | Path) -> ad.AnnData:
    """Return PT structures from the four current mouse count matrices.

    ``layers['counts']`` contains raw aggregated counts and ``layers['lognorm']`` contains
    Scanpy total-count normalization to 10,000 followed by natural log1p. ``X`` is lognorm.
    All genes shared by the four inputs are retained; no HVG or tubule gene-count filter runs.
    The notebook 02 pooled mouse p5 ``n_spots`` structural floor is applied before PT selection.
    """
    root = Path(data_root).expanduser().resolve()
    matrix_dir = root / "tubule_by_gene"
    files = [f"{sample}_tubule_by_gene_caleb.h5ad" for sample in MOUSE_SAMPLES]
    absent = [str(matrix_dir / f) for f in files if not (matrix_dir / f).is_file()]
    if absent:
        raise FileNotFoundError("Missing mouse tubule count matrices: " + ", ".join(absent))

    sample_adatas = load_and_process_mouse_samples(matrix_dir, files)
    # Attach current labels only after validating every source observation against its exact
    # feature index and centroid.  The whole current input set is required to pass.
    validation = {}
    for sample in MOUSE_SAMPLES:
        a = sample_adatas[sample]
        # load_and_process_mouse_samples materializes layers['counts'] into X when present.
        # Validate before combine_adatas can make duplicate feature names unique.
        _counts_layer(a, sample)
        labels, evidence = _read_and_validate_segmentation(root, sample, a.obs)
        a.obs["feature_index"] = pd.to_numeric(a.obs["feature_index"], errors="raise").astype("int64")
        a.obs["segment_class"] = labels
        a.obs["coarse_class"] = [coarse_for(x) for x in labels]
        a.obs["broad_tubule_marker_call"] = a.obs["coarse_class"].astype(str)
        validation[sample] = evidence

    combined = combine_adatas(sample_adatas)
    if "n_spots" not in combined.obs:
        raise ValueError("Mouse count matrices lack required n_spots structural support")
    spots = pd.to_numeric(combined.obs.n_spots, errors="coerce").to_numpy(float)
    if not np.isfinite(spots).all() or (spots < 0).any():
        raise ValueError("n_spots must contain finite, nonnegative values")
    floor = float(np.quantile(spots, SPOT_QUANTILE))
    before = int(combined.n_obs)
    keep = spots >= floor
    combined = combined[keep].copy()
    after_spot_floor = int(combined.n_obs)

    # Match the analysis-wide denominator: normalize each retained structure using all common
    # measured genes before the PT view is selected. This keeps PT scaling independent of labels.
    counts = combined.X.copy()
    if sparse.issparse(counts):
        counts = counts.tocsr()
    combined.layers["counts"] = counts.copy()
    lognorm = combined.copy()
    sc.pp.normalize_total(lognorm, target_sum=1e4)
    sc.pp.log1p(lognorm)
    combined.layers["lognorm"] = lognorm.X.copy()
    combined.X = combined.layers["lognorm"].copy()

    pt = combined.obs.coarse_class.astype(str).eq("PT").to_numpy()
    combined = combined[pt].copy()
    if combined.n_obs == 0:
        raise ValueError("No PT structures remain after pooled p5 n_spots floor")

    # Avoid calculating or applying per-structure gene QC.
    combined.obs["structure_id"] = (
        combined.obs["sample"].astype(str) + "|" + combined.obs["feature_index"].astype(str)
    ).to_numpy()
    if not combined.obs.structure_id.is_unique:
        raise ValueError("sample|feature_index structure IDs are not unique")
    combined.uns["repeated_mouse_input_manifest"] = {
        "samples": list(MOUSE_SAMPLES),
        "source_matrices": {
            s: {
                "path": f"tubule_by_gene/{s}_tubule_by_gene_caleb.h5ad",
                "size_bytes": (matrix_dir / f"{s}_tubule_by_gene_caleb.h5ad").stat().st_size,
                "mtime_ns": (matrix_dir / f"{s}_tubule_by_gene_caleb.h5ad").stat().st_mtime_ns,
            }
            for s in MOUSE_SAMPLES
        },
        "segmentation_validation": validation,
        "centroid_validation": "all matrix observations checked against GeoJSON feature_index",
        "counts_layer": "raw structure-aggregated X from input matrices",
        "normalization": "after p5 spot floor and before PT selection: Scanpy normalize_total(target_sum=10000) across all common measured genes, then natural log1p",
        "genes": "inner intersection across four samples; no HVG selection",
        "n_spots_quantile": SPOT_QUANTILE,
        "n_spots_floor": floor,
        "structures_before_floor": before,
        "structures_after_floor_all_classes": after_spot_floor,
        "structures_after_pt_selection": int(combined.n_obs),
        "no_tubule_gene_count_filter": True,
    }
    return combined
