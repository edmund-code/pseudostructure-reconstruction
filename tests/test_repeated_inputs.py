"""Small input-integrity regressions for measured repeated-structure loading."""
import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from pseudospace.repeated_inputs import _counts_layer, _read_and_validate_segmentation


def _square_feature():
    return {
        "type": "Feature",
        "properties": {"segment_class": "PT-S1"},
        "geometry": {
            "type": "Polygon",
            "coordinates": [[[0, 0], [2, 0], [2, 2], [0, 2], [0, 0]]],
        },
    }


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
