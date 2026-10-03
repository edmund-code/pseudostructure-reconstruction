"""Synthetic regressions for the canonical molecular mean atlas."""
import inspect

import numpy as np
import pytest

from pseudospace.canonical_mean_atlas import (
    conditional_gene_predictions,
    fit_global_offset,
    fit_mean_atlas,
    project_mean_atlas,
    segment_means,
)


def _curve(z):
    z = np.asarray(z)
    return np.column_stack([z, np.sin(2 * np.pi * z), np.cos(np.pi * z)])


def test_projection_recovers_positions_and_returns_grid_posterior():
    rng = np.random.default_rng(22)
    z_train = np.linspace(0, 1, 501)
    y_train = _curve(z_train) + rng.normal(0, 0.01, (len(z_train), 3))
    atlas = fit_mean_atlas(y_train, z_train)
    z_test = np.linspace(0.04, 0.96, 31)
    projected = project_mean_atlas(atlas, _curve(z_test))

    assert np.mean(np.abs(projected["z_MAP"] - z_test)) < 0.02
    assert projected["posterior"].shape == (len(z_test), 101)
    assert np.allclose(projected["posterior"].sum(axis=1), 1)
    assert np.all(projected["entropy"] >= 0)
    assert np.all(projected["residual_norm"] >= 0)
    assert np.mean(np.abs(projected["z_mean"] - z_test)) < 0.03
    assert set(inspect.signature(project_mean_atlas).parameters) == {"atlas", "Y"}


def test_global_offset_uses_fixed_atlas_and_returns_one_shrunk_vector():
    z = np.linspace(0, 1, 301)
    train = _curve(z)
    atlas = fit_mean_atlas(train, z)
    mean_before = atlas["mean"].copy()
    offset_true = np.array([0.5, -0.25, 0.1])
    estimated, history = fit_global_offset(atlas, _curve(z[::4]) + offset_true)

    assert estimated.shape == (train.shape[1],)
    assert np.allclose(estimated, offset_true / 5, atol=0.025)
    assert len(history) == 5
    assert np.array_equal(atlas["mean"], mean_before)


def test_conditional_gene_prediction_uses_training_profiles_only():
    z = np.linspace(0, 1, 201)
    values = np.column_stack([z, 3 - 2 * z])
    predicted = conditional_gene_predictions(values, z, [0.5], bandwidth=0.03)
    assert predicted.shape == (1, 2)
    assert np.allclose(predicted[0], [0.5, 2.0], atol=0.01)


def test_tiny_bandwidth_keeps_nearest_kernel_support_finite():
    z = np.array([0.0, 0.5, 1.0])
    values = np.array([[1.0], [2.0], [3.0]])
    predicted = conditional_gene_predictions(values, z, [0.49], bandwidth=1e-300)
    assert np.isfinite(predicted).all()
    assert predicted[0, 0] == pytest.approx(2.0)
    atlas = fit_mean_atlas(values, z, grid=np.array([0.0, 0.5, 1.0]),
                           bandwidth=1e-300, min_effective=1)
    assert np.isfinite(atlas["mean"]).all()


def test_segment_mean_oracle_returns_exact_training_means_and_counts():
    x_train = np.array([[1., 0.], [3., 2.], [10., 4.], [14., 8.], [20., 9.]])
    train_labels = np.array(["S1", "S1", "S2", "S2", "S3"])
    out = segment_means(x_train, train_labels, ["S3", "S1", "S2"])
    assert np.array_equal(out["prediction"], [[20., 9.], [2., 1.], [12., 6.]])
    assert out["counts"] == {"S1": 2, "S2": 2, "S3": 1}
    with pytest.raises(ValueError, match="support each"):
        segment_means(x_train[:4], train_labels[:4], ["S3"])


def test_atlas_refuses_insufficient_effective_support():
    z = np.linspace(0, 1, 20)
    with pytest.raises(ValueError, match="insufficient effective training support"):
        fit_mean_atlas(_curve(z), z, min_effective=12)


def test_atlas_rejects_nonfinite_and_mismatched_inputs():
    z = np.linspace(0, 1, 301)
    x = _curve(z)
    with pytest.raises(ValueError, match="one finite value per profile"):
        fit_mean_atlas(x, z[:-1])
    with pytest.raises(ValueError, match="finite"):
        project_mean_atlas(fit_mean_atlas(x, z), [[np.nan, 0, 0]])
