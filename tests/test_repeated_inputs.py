"""Small input-integrity regressions for measured repeated-structure loading."""
import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from pseudospace.repeated_inputs import (
    MOUSE_SAMPLES,
    _counts_layer,
    _read_and_validate_segmentation,
    load_repeated_mouse_inputs,
)
import pseudospace.repeated_inputs as repeated_inputs


def _square_feature():
    return {
        "type": "Feature",
        "properties": {"segment_class": "PT-S1"},
        "geometry": {
            "type": "Polygon",
            "coordinates": [[[0, 0], [2, 0], [2, 2], [0, 2], [0, 0]]],
        },
    }


@pytest.mark.parametrize("samples", [(), ("Ctrl1A2", "Ctrl1A2"), ("unknown",), "Ctrl1A2"])
def test_loader_rejects_invalid_sample_selection_before_input_io(tmp_path, samples):
    with pytest.raises(ValueError, match="nonempty unique subset"):
        load_repeated_mouse_inputs(tmp_path / "does-not-exist", samples=samples)


def test_mouse_sample_contract_is_explicit():
    assert MOUSE_SAMPLES == ("Ctrl1A2", "Ctrl1A4", "IR2A2", "IR2A4")


def test_loader_can_load_healthy_subset_without_aki_inputs(tmp_path, monkeypatch):
    healthy = ("Ctrl1A2", "Ctrl1A4")
    matrix_dir = tmp_path / "tubule_by_gene"
    matrix_dir.mkdir()
    for sample in healthy:
        (matrix_dir / f"{sample}_tubule_by_gene_caleb.h5ad").touch()
    captured = {}

    def feature(index, segment):
        x0 = index * 3
        return {
            "type": "Feature",
            "properties": {"segment_class": segment},
            "geometry": {
                "type": "Polygon",
                "coordinates": [[[x0, 0], [x0 + 2, 0], [x0 + 2, 2], [x0, 2], [x0, 0]]],
            },
        }

    def fake_load(data_dir, files):
        captured["files"] = list(files)
        sample_adatas = {}
        for sample in healthy:
            (tmp_path / f"{sample}_kept_tubules_labeled_fine.geojson").write_text(
                json.dumps({"features": [feature(0, "PT-S1"), feature(1, "PT-S2")]})
            )
            obs = pd.DataFrame({
                "feature_index": [0, 1],
                "x_centroid": [1.0, 4.0],
                "y_centroid": [1.0, 1.0],
                "n_spots": [10, 20],
                "sample": [sample, sample],
                "condition": ["Control", "Control"],
            })
            sample_adatas[sample] = ad.AnnData(
                X=sparse.csr_matrix([[2, 1], [1, 3]]),
                obs=obs,
                var=pd.DataFrame(index=["GeneA", "GeneB"]),
            )
            sample_adatas[sample].obs_names = [f"{sample}_unit_0", f"{sample}_unit_1"]
        return sample_adatas

    monkeypatch.setattr(repeated_inputs, "load_and_process_mouse_samples", fake_load)
    loaded = load_repeated_mouse_inputs(tmp_path, samples=healthy)

    assert captured["files"] == [f"{sample}_tubule_by_gene_caleb.h5ad" for sample in healthy]
    manifest = loaded.uns["repeated_mouse_input_manifest"]
    assert manifest["samples"] == list(healthy)
    assert manifest["n_spots_cohort"] == list(healthy)
    assert set(manifest["source_matrices"]) == set(healthy)
    assert set(manifest["segmentation_validation"]) == set(healthy)
    assert set(loaded.obs["sample"].astype(str)) == set(healthy)
    assert not any(sample.startswith("IR") for sample in manifest["source_matrices"])


def test_segmentation_validation_rejects_nonfinite_source_centroid(tmp_path):
    (tmp_path / "sample_kept_tubules_labeled_fine.geojson").write_text(
        json.dumps({"features": [_square_feature()]})
    )
    obs = pd.DataFrame({"feature_index": [0], "x_centroid": [np.nan], "y_centroid": [1.]})
    with pytest.raises(ValueError, match="non-finite"):
        _read_and_validate_segmentation(tmp_path, "sample", obs)


def test_segmentation_validation_rejects_centroid_mismatch(tmp_path):
    (tmp_path / "sample_kept_tubules_labeled_fine.geojson").write_text(
        json.dumps({"features": [_square_feature()]})
    )
    obs = pd.DataFrame({"feature_index": [0], "x_centroid": [3.], "y_centroid": [1.]})
    with pytest.raises(ValueError, match="centroid mismatch"):
        _read_and_validate_segmentation(tmp_path, "sample", obs)


def test_segmentation_validation_rejects_duplicate_feature_indices(tmp_path):
    (tmp_path / "sample_kept_tubules_labeled_fine.geojson").write_text(
        json.dumps({"features": [_square_feature(), _square_feature()]})
    )
    obs = pd.DataFrame({"feature_index": [0, 0], "x_centroid": [1., 1.], "y_centroid": [1., 1.]})
    with pytest.raises(ValueError, match="not unique"):
        _read_and_validate_segmentation(tmp_path, "sample", obs)


@pytest.mark.parametrize(
    "matrix, var_names, message",
    [
        (np.array([[1., np.nan]]), ["a", "b"], "finite and nonnegative"),
        (np.array([[1., -1.]]), ["a", "b"], "finite and nonnegative"),
        (np.array([[1., 0.5]]), ["a", "b"], "integer-valued"),
        (np.array([[1., 2.]]), ["a", "a"], "duplicate gene names"),
    ],
)
def test_raw_count_validation_rejects_bad_inputs(matrix, var_names, message):
    a = ad.AnnData(X=sparse.csr_matrix(matrix), var=pd.DataFrame(index=var_names))
    with pytest.raises(ValueError, match=message):
        _counts_layer(a, "sample")


def test_raw_count_validation_rejects_empty_structures():
    a = ad.AnnData(X=sparse.csr_matrix([[0, 0], [1, 0]]))
    with pytest.raises(ValueError, match="positive finite"):
        _counts_layer(a, "sample")
