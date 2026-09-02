"""No-private-data regression coverage for reusable pseudospace helpers."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from pathlib import Path

ad = pytest.importorskip("anndata")
pytest.importorskip("scanpy")

from pseudospace.heatmaps import binned_interpolated_expression
from pseudospace.harmony import require_supported_harmony_version
from pseudospace.io_qc import annotate_mito_ribo_mouse_symbols, sample_name_from_path
from pseudospace.levelshape import loso_shape_stability, run_level_shape
from pseudospace.markers import assign_cluster_labels, resolve_available_marker_groups
from pseudospace.trajectory import orient_and_normalize


def test_harmony_version_guard_rejects_legacy_r_package():
    require_supported_harmony_version("2.0.5")
    require_supported_harmony_version("[1] 2 0 5")
    with pytest.raises(RuntimeError, match="requires harmony >= 2.0.5"):
        require_supported_harmony_version("1.2.4")


def _adata(n_obs: int = 24):
    x = np.column_stack([
        np.linspace(1, 4, n_obs),
        np.linspace(4, 1, n_obs),
        np.full(n_obs, 2.0),
        np.arange(n_obs) % 3,
    ])
    obj = ad.AnnData(X=x)
    obj.var_names = ["mt-Nd1", "Rpl3", "Slc5a2", "Aqp2"]
    obj.obs["pseudotime"] = np.linspace(0, 1, n_obs)
    obj.layers["lognorm"] = x.copy()
    return obj


def test_qc_marker_annotation_and_filename_contract():
    obj = annotate_mito_ribo_mouse_symbols(_adata())
    assert obj.var.loc["mt-Nd1", "exclude_from_harmony_hvg"]
    assert obj.var.loc["Rpl3", "exclude_from_harmony_hvg"]
    assert not obj.var.loc["Slc5a2", "exclude_from_harmony_hvg"]
    assert sample_name_from_path(Path("Ctrl1A2_tubule_by_gene_caleb.h5ad")) == "Ctrl1A2"


def test_marker_groups_and_manual_labels_are_explicit():
    obj = _adata()
    groups = {
        "PT": [{"label": "PT: Slc5a2", "candidates": ["slc5a2"]}],
        "CD": [{"label": "CD: Aqp2", "candidates": ["AQP2"]}],
        "missing": [{"label": "Missing", "candidates": ["NotAGene"]}],
    }
    resolved, labels, indices, missing = resolve_available_marker_groups(
        obj, groups, ["PT", "CD", "missing"]
    )
    assert resolved == {"PT": ["PT: Slc5a2"], "CD": ["CD: Aqp2"]}
    assert labels == ["PT: Slc5a2", "CD: Aqp2"]
    assert indices == [2, 3]
    assert missing == ["Missing"]

    obj.obs["cluster"] = ["0", "1"] * 12
    assign_cluster_labels(obj, "cluster", {"0": "PT", "1": "CD"}, "segment")
    assert set(obj.obs["segment"].astype(str)) == {"PT", "CD"}


def test_trajectory_orientation_and_heatmap_binning():
    marker_axis = np.linspace(-2, 2, 24)
    oriented = orient_and_normalize(np.linspace(1, 0, 24), marker_axis, min_valid=10)
    assert np.isclose(oriented.min(), 0)
    assert np.isclose(oriented.max(), 1)
    assert np.corrcoef(oriented, marker_axis)[0, 1] > 0.99

    raw, zscore, centers, edges, occupied = binned_interpolated_expression(
        _adata(), [0, 1], "pseudotime", n_bins=8, smooth_sigma=0
    )
    assert raw.shape == zscore.shape == (2, 8)
    assert len(centers) == 8 and len(edges) == 9
    assert occupied.all()


def test_level_shape_and_loso_stability_on_synthetic_replicates():
    n_per_sample = 12
    samples = np.repeat(["Ctrl1", "Ctrl2", "IR1", "IR2"], n_per_sample)
    condition = np.isin(samples, ["IR1", "IR2"]).astype(float)
    position = np.tile(np.linspace(0.02, 0.98, n_per_sample), 4)
    # One condition-dependent shape feature and one stable feature.
    response = np.column_stack([
        np.sin(np.pi * position) + condition * position**2,
        0.2 + np.cos(np.pi * position),
    ])
    knots = np.linspace(0.15, 0.85, 4)
    grid = np.linspace(0.05, 0.95, 25)
    lambdas = np.array([0.01, 0.1, 1.0])

    fitted = run_level_shape(response, position, condition, knots, grid, lambdas)
    assert fitted["shape_rms"][0] > fitted["shape_rms"][1]

    stability = loso_shape_stability(
        response, position, condition, samples, knots, grid, lambdas,
        top_idx=[0], feature_names=["shape_feature", "stable_feature"],
    )
    assert len(stability) == 4
    assert stability["loso_shape_rms"].notna().all()
    assert (stability["feature"] == "shape_feature").all()
