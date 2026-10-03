import numpy as np
import pytest

from pseudospace.residual_metric_atlas import (
    fit_residual_metric_atlas, project_residual_metric_atlas,
)


def test_correlated_residual_structure_changes_direction_weighting():
    rng = np.random.default_rng(31)
    n = 700
    z = rng.uniform(0, 1, n)
    common = rng.normal(0, 1.1, n)
    contrast = rng.normal(0, .16, n)
    Y = np.column_stack([z, common + contrast, common - contrast])
    atlas = fit_residual_metric_atlas(Y, z, bandwidth=.12, min_effective=12)
    center = atlas['mean'][50]
    high_variance = center + np.array([0., 1 / np.sqrt(2), 1 / np.sqrt(2)])
    low_variance = center + np.array([0., 1 / np.sqrt(2), -1 / np.sqrt(2)])
    projected = project_residual_metric_atlas(atlas, np.vstack([high_variance, low_variance]))
    assert projected['log_density'][0] > projected['log_density'][1]
    assert np.linalg.norm(high_variance - center) == pytest.approx(np.linalg.norm(low_variance - center))
    assert not np.allclose(atlas['precision'], np.eye(3) / np.trace(atlas['covariance']) * 3)


def test_singular_residual_directions_are_floored_and_projection_is_finite():
    z = np.linspace(0, 1, 240)
    Y = np.column_stack([z, 2 * z, np.ones_like(z)])
    atlas = fit_residual_metric_atlas(Y, z, bandwidth=.12, min_effective=12)
    eig = np.linalg.eigvalsh(atlas['covariance'])
    assert np.all(eig > 0)
    result = project_residual_metric_atlas(atlas, Y[[0, 80, 239]])
    assert np.isfinite(result['posterior']).all()
    np.testing.assert_allclose(result['posterior'].sum(axis=1), 1.)
    assert np.isfinite(result['entropy']).all()
    assert np.isfinite(result['residual_norm']).all()
    assert np.isfinite(result['log_density']).all()
    assert np.isfinite(result['log_density_at_MAP']).all()


def test_constant_training_profiles_get_finite_floored_covariance():
    z = np.linspace(0, 1, 240)
    Y = np.ones((len(z), 3)) * np.array([2., -1., .5])
    atlas = fit_residual_metric_atlas(Y, z, bandwidth=.15, min_effective=12)
    assert np.all(np.linalg.eigvalsh(atlas['covariance']) > 0)
    result = project_residual_metric_atlas(atlas, np.vstack([Y[0], Y[-1]]))
    np.testing.assert_allclose(result['posterior'].sum(axis=1), 1.)
    assert np.isfinite(result['posterior']).all()
    assert np.isfinite(result['log_density']).all()
    assert np.isfinite(result['log_density_at_MAP']).all()


def test_projection_is_invariant_to_query_batching():
    rng = np.random.default_rng(9)
    z = rng.uniform(0, 1, 260)
    Y = np.column_stack([z, np.sin(3 * z), np.cos(2 * z)]) + rng.normal(0, .04, (260, 3))
    atlas = fit_residual_metric_atlas(Y, z, bandwidth=.15, min_effective=12)
    query = Y[[4, 35, 130, 251]] + .01
    full = project_residual_metric_atlas(atlas, query)
    left = project_residual_metric_atlas(atlas, query[:2])
    right = project_residual_metric_atlas(atlas, query[2:])
    np.testing.assert_allclose(full['posterior'], np.vstack([left['posterior'], right['posterior']]))
    np.testing.assert_allclose(full['z_mean'], np.r_[left['z_mean'], right['z_mean']])


@pytest.mark.parametrize('bad', [np.array([[np.nan, 0.]]), np.array([[1.]])])
def test_projection_rejects_nonfinite_or_wrong_feature_count(bad):
    z = np.linspace(0, 1, 200)
    Y = np.column_stack([z, z**2])
    atlas = fit_residual_metric_atlas(Y, z, bandwidth=.15, min_effective=12)
    with pytest.raises(ValueError):
        project_residual_metric_atlas(atlas, bad)
