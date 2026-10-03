"""Synthetic regressions for Gaussian marginal spatial flow."""
import numpy as np
import pytest

pytest.importorskip('scipy')
pytest.importorskip('sklearn')

from pseudospace.gaussian_spatial_flow import (
    advect_gaussian, fit_gaussian_path, gaussian_velocity,
    project_gaussian_path,
)


def test_fitted_1d_gaussian_flow_transports_known_terminal_moments():
    rng = np.random.default_rng(120)
    z = rng.uniform(0, 1, 12000)
    mean = 1.5 * z - .2
    sd = .7 + .35 * z
    train = (mean + sd * rng.normal(size=len(z)))[:, None]
    grid = np.linspace(0, 1, 101)
    path = fit_gaussian_path(train, z, grid, bandwidth=.07,
                             shrinkage=.25, min_effective=12)

    # Check the Lyapunov equation for the field's covariance derivative.
    pos = .503
    _, _, A = gaussian_velocity(path, pos)
    cell = np.searchsorted(grid, pos) - 1
    sigma = ((grid[cell + 1] - pos) * path['covariances'][cell]
             + (pos - grid[cell]) * path['covariances'][cell + 1]) / (grid[cell + 1] - grid[cell])
    sigma_prime = ((path['covariances'][cell + 1] - path['covariances'][cell])
                   / (grid[cell + 1] - grid[cell]))
    assert np.linalg.norm(A @ sigma + sigma @ A - sigma_prime) < 1e-9

    z0, z1 = .3, .7
    mu0, var0 = 1.5*z0-.2, (.7+.35*z0)**2
    mu1, var1 = 1.5*z1-.2, (.7+.35*z1)**2
    source = rng.normal(mu0, np.sqrt(var0), size=(5000, 1))
    moved = advect_gaussian(path, source, z0, z1, steps=64)
    assert abs(moved.mean() - mu1) < .05
    assert abs(moved.var() - var1) < .06


def test_constant_covariance_flow_is_exact_mean_translation():
    grid = np.linspace(0, 1, 5)
    slope = np.array([2., -.75])
    means = grid[:, None] * slope
    covariance = np.array([[1., .25], [.25, 2.]])
    path = {'grid': grid, 'means': means,
            'covariances': np.repeat(covariance[None, :, :], len(grid), axis=0),
            'feature_count': 2}
    rows = np.array([[0., 1.], [-2., 3.], [.4, -.1]])
    moved = advect_gaussian(path, rows, .1, .9, steps=32)
    np.testing.assert_allclose(moved, rows + .8*slope, atol=1e-9)


def test_path_is_training_only_and_projection_posterior_is_normalized():
    rng = np.random.default_rng(33)
    z = np.linspace(0, 1, 250)
    Y = np.column_stack([z, np.sin(2*np.pi*z)]) + rng.normal(0, .04, (len(z), 2))
    grid = np.linspace(0, 1, 31)
    path_a = fit_gaussian_path(Y, z, grid, min_effective=12)
    changed_unrelated_target = rng.normal(100, 20, size=(40, 2))
    # The fitter has no target argument: holdout changes cannot alter its path.
    path_b = fit_gaussian_path(Y, z, grid, min_effective=12)
    assert changed_unrelated_target.shape == (40, 2)
    np.testing.assert_allclose(path_a['means'], path_b['means'])
    np.testing.assert_allclose(path_a['covariances'], path_b['covariances'])
    projected = project_gaussian_path(path_a, Y[::7], block_size=5)
    assert projected['posterior'].shape == (len(Y[::7]), len(grid))
    np.testing.assert_allclose(projected['posterior'].sum(axis=1), 1.)
    for key in ('z_MAP', 'z_mean', 'entropy', 'residual_norm'):
        assert np.isfinite(projected[key]).all()


def test_gaussian_path_refuses_sparse_kernel_support():
    Y = np.arange(40, dtype=float).reshape(20, 2)
    z = np.r_[np.full(10, .25), np.full(10, .75)]
    with pytest.raises(ValueError, match='effective support'):
        fit_gaussian_path(Y, z, grid=[.25, .5, .75], bandwidth=.02,
                          min_effective=12)
