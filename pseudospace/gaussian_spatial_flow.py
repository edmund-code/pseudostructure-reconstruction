"""Smooth Gaussian marginal probability flow for spatial molecular states.

The fitted velocity is the Gaussian probability-flow field for the estimated
mean and covariance path. It is a distributional interpolation, not an
identified biological mechanism or a unique particle-level trajectory.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
from scipy.linalg import cholesky, solve_sylvester, solve_triangular

from pseudospace.spatial_probability_flow import energy_distance, window_clouds


def _matrix(X, name, min_rows=2):
    X = np.asarray(X, dtype=float)
    if X.ndim != 2 or X.shape[0] < min_rows or X.shape[1] < 1:
        raise ValueError(f'{name} must be a 2D matrix with at least {min_rows} rows and one feature')
    if not np.isfinite(X).all():
        raise ValueError(f'{name} must contain only finite values')
    return X


def fit_gaussian_path(Y, z, grid=None, bandwidth=.08, shrinkage=.25,
                      variance_floor=1e-3, min_effective=12):
    """Estimate a kernel-weighted Gaussian marginal at each grid position.

    Each position uses all training structures with Gaussian kernel weights.
    Effective support is reported and enforced before covariance estimation.
    Covariances are shrunk toward their diagonal and receive a small global
    variance floor to keep every marginal positive definite.
    """
    Y = _matrix(Y, 'Y')
    z = np.asarray(z, dtype=float).ravel()
    if len(z) != len(Y) or not np.isfinite(z).all():
        raise ValueError('z must be finite and match Y rows')
    if not np.isfinite(bandwidth) or bandwidth <= 0:
        raise ValueError('bandwidth must be positive and finite')
    if not 0 <= shrinkage <= 1 or not np.isfinite(shrinkage):
        raise ValueError('shrinkage must be finite and in [0,1]')
    if not np.isfinite(variance_floor) or variance_floor <= 0:
        raise ValueError('variance_floor must be positive and finite')
    if not np.isfinite(min_effective) or min_effective < 2:
        raise ValueError('min_effective must be at least 2')
    if grid is None:
        grid = np.linspace(0., 1., 101)
    else:
        grid = np.asarray(grid, dtype=float).ravel()
    if len(grid) < 2 or not np.isfinite(grid).all() or np.any(np.diff(grid) <= 0):
        raise ValueError('grid must be finite, strictly increasing, and contain at least two values')

    global_center = Y.mean(axis=0)
    global_cov = np.cov(Y, rowvar=False, ddof=1)
    global_cov = np.atleast_2d(global_cov)
    global_scale = max(float(np.trace(global_cov) / Y.shape[1]), np.finfo(float).eps)
    floor = global_scale * variance_floor
    means, covariances, support = [], [], []
    for position in grid:
        weights = np.exp(-.5 * ((z - position) / bandwidth) ** 2)
        sum_w, sum_w2 = weights.sum(), np.dot(weights, weights)
        n_eff = sum_w**2 / sum_w2 if sum_w2 > 0 else 0.
        support.append(n_eff)
        if n_eff < min_effective:
            raise ValueError(f'Gaussian path at z={position:g} has effective support {n_eff:.2f}; '
                             f'requires {min_effective}')
        weights /= sum_w
        mean = weights @ Y
        centered = Y - mean
        cov = (centered * weights[:, None]).T @ centered
        diagonal = np.diag(np.diag(cov))
        cov = (1 - shrinkage) * cov + shrinkage * diagonal
        cov = (cov + cov.T) * .5 + floor * np.eye(Y.shape[1])
        means.append(mean)
        covariances.append(cov)
    covariances = np.asarray(covariances)
    factors = np.asarray([cholesky(cov, lower=True) for cov in covariances])
    return {'grid': grid, 'means': np.asarray(means),
            'covariances': covariances, 'cholesky': factors,
            'support': np.asarray(support),
            'bandwidth': float(bandwidth), 'shrinkage': float(shrinkage),
            'variance_floor': float(floor), 'feature_count': Y.shape[1]}


def project_gaussian_path(path, Y, block_size=256):
    """Project rows onto the frozen Gaussian atlas using full-covariance logpdfs.

    Returns MAP and posterior-mean grid coordinates, the grid posterior,
    entropy, and Euclidean residual from the posterior-averaged molecular mean.
    A uniform prior is used over the fitted grid; no query labels or DPT enter.
    """
    Y = _matrix(Y, 'Y', min_rows=1)
    if Y.shape[1] != path['feature_count']:
        raise ValueError('Y feature count does not match Gaussian path')
    if not isinstance(block_size, (int, np.integer)) or block_size < 1:
        raise ValueError('block_size must be a positive integer')
    means = path['means']
    factors = path['cholesky']
    grid = path['grid']
    d, k = Y.shape[1], len(grid)
    logp = np.empty((len(Y), k), dtype=float)
    for j in range(k):
        L = factors[j]
        logdet = 2 * np.log(np.diag(L)).sum()
        for start in range(0, len(Y), int(block_size)):
            stop = min(start + int(block_size), len(Y))
            residual = (Y[start:stop] - means[j]).T
            whitened = solve_triangular(L, residual, lower=True, check_finite=False)
            logp[start:stop, j] = -.5 * (d*np.log(2*np.pi) + logdet
                                         + np.sum(whitened**2, axis=0))
    logp -= np.log(k)  # uniform grid prior
    row_max = np.max(logp, axis=1, keepdims=True)
    log_norm = row_max + np.log(np.exp(logp - row_max).sum(axis=1, keepdims=True))
    log_post = logp - log_norm
    probabilities = np.exp(log_post)
    z_map = grid[np.argmax(logp, axis=1)]
    z_mean = probabilities @ grid
    entropy = -np.sum(probabilities * log_post, axis=1)
    expected_mean = probabilities @ means
    residual_norm = np.linalg.norm(Y - expected_mean, axis=1)
    return {'z_MAP': z_map, 'z_mean': z_mean, 'posterior': probabilities,
            'entropy': entropy, 'residual_norm': residual_norm}


def gaussian_velocity(path, z):
    """Interpolate Gaussian moments and solve the symmetric Lyapunov flow."""
    grid = path['grid']
    z = float(z)
    if not np.isfinite(z) or z < grid[0] or z > grid[-1]:
        raise ValueError('z must lie within the fitted Gaussian path grid')
    mean_spline = CubicSpline(grid, path['means'], axis=0)
    mean = np.asarray(mean_spline(z))
    mean_prime = np.asarray(mean_spline(z, 1))
    covs = path['covariances']
    i = min(int(np.searchsorted(grid, z, side='right') - 1), len(grid) - 2)
    i = max(i, 0)
    dz = grid[i + 1] - grid[i]
    frac = (z - grid[i]) / dz
    cov = (1 - frac) * covs[i] + frac * covs[i + 1]
    cov_prime = (covs[i + 1] - covs[i]) / dz
    A = solve_sylvester(cov, cov, cov_prime)
    A = (A + A.T) * .5
    return mean, mean_prime, A


def _velocity(path, X, z):
    mean, mean_prime, A = gaussian_velocity(path, z)
    return mean_prime + (X - mean) @ A.T


def advect_gaussian(path, Y, z0, z1, steps=32):
    """RK4 advection under the fitted affine Gaussian probability flow."""
    Y = _matrix(Y, 'Y', min_rows=1)
    if Y.shape[1] != path['feature_count']:
        raise ValueError('Y feature count does not match Gaussian path')
    if not np.isfinite([z0, z1]).all() or steps < 1:
        raise ValueError('z endpoints must be finite and steps positive')
    if min(z0, z1) < path['grid'][0] or max(z0, z1) > path['grid'][-1]:
        raise ValueError('z endpoints must lie within fitted Gaussian path grid')
    h = (float(z1) - float(z0)) / int(steps)
    X = Y.copy()
    for j in range(int(steps)):
        z = float(z0) + j * h
        k1 = _velocity(path, X, z)
        k2 = _velocity(path, X + h*k1/2, z + h/2)
        k3 = _velocity(path, X + h*k2/2, z + h/2)
        k4 = _velocity(path, X + h*k3, z + h)
        X += h * (k1 + 2*k2 + 2*k3 + k4) / 6
    if not np.isfinite(X).all():
        raise RuntimeError('Gaussian flow integration diverged')
    return X


def _sample_cov(X):
    cov = np.cov(X, rowvar=False, ddof=1)
    return np.atleast_2d(cov)


def gaussian_transport_benchmark(Ytrain, ztrain, Ytest, ztest, centers,
                                 halfwidth, min_count=12, max_points=128,
                                 seed=15, bandwidth=.08, shrinkage=.25,
                                 variance_floor=1e-3, min_effective=12,
                                 steps=32):
    """Evaluate frozen Gaussian spatial flow against three empirical baselines."""
    Ytrain = _matrix(Ytrain, 'Ytrain')
    Ytest = _matrix(Ytest, 'Ytest')
    ztrain = np.asarray(ztrain, float).ravel()
    ztest = np.asarray(ztest, float).ravel()
    centers = np.asarray(centers, float).ravel()
    if len(ztrain) != len(Ytrain) or len(ztest) != len(Ytest):
        raise ValueError('ztrain and ztest must match their molecular matrices')
    if Ytrain.shape[1] != Ytest.shape[1]:
        raise ValueError('train and test feature counts must match')
    if len(centers) < 2 or not np.isfinite(centers).all() or np.any(np.diff(centers) <= 0):
        raise ValueError('centers must be finite and strictly increasing')
    train_clouds = window_clouds(Ytrain, ztrain, centers, halfwidth,
                                 min_count, max_points, seed)
    test_clouds = window_clouds(Ytest, ztest, centers, halfwidth,
                                min_count, max_points, seed + 1)
    grid = np.linspace(float(np.min(ztrain)), float(np.max(ztrain)), 101)
    path = fit_gaussian_path(Ytrain, ztrain, grid, bandwidth, shrinkage,
                             variance_floor, min_effective)
    rows = []
    for i, (z0, z1) in enumerate(zip(centers[:-1], centers[1:])):
        source, target = test_clouds[i], test_clouds[i + 1]
        train_shift = train_clouds[i + 1].mean(axis=0) - train_clouds[i].mean(axis=0)
        mu0 = gaussian_velocity(path, z0)[0]
        mu1 = gaussian_velocity(path, z1)[0]
        predictions = {
            'identity': source,
            'centroid': source + train_shift,
            'smooth_mean': source + (mu1 - mu0),
            'gaussian_flow': advect_gaussian(path, source, z0, z1, steps),
        }
        target_centered = target - target.mean(axis=0)
        for method, pred in predictions.items():
            pred_centered = pred - pred.mean(axis=0)
            rows.append({
                'transition': i, 'z_from': z0, 'z_to': z1, 'method': method,
                'energy_distance': energy_distance(pred, target),
                'mean_error': float(np.linalg.norm(pred.mean(axis=0) - target.mean(axis=0))),
                'covariance_frobenius_error': float(np.linalg.norm(_sample_cov(pred) - _sample_cov(target), ord='fro')),
                'shape_energy_distance': energy_distance(pred_centered, target_centered),
                'n_source': len(source), 'n_target': len(target),
            })
    return pd.DataFrame(rows)
