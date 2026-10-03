"""Mean atlas with one frozen global covariance for residual-aware projection."""
from __future__ import annotations

import numpy as np
from scipy.special import logsumexp, xlogy
from sklearn.covariance import LedoitWolf

from pseudospace.canonical_mean_atlas import (
    _coordinates, _grid, _matrix, fit_mean_atlas,
)


def fit_residual_metric_atlas(Ytrain, ztrain, grid=None, bandwidth=.08,
                              min_effective=12):
    """Fit a kernel mean path and one Ledoit-Wolf covariance of residuals.

    Residuals are measured from the interpolated training mean path at each
    supplied training coordinate. Covariance is global across positions; no
    local covariance, flow, or query-dependent fitting is used.
    """
    y = _matrix(Ytrain, 'Ytrain')
    z = _coordinates(ztrain, len(y), 'ztrain')
    mean_atlas = fit_mean_atlas(y, z, grid=grid, bandwidth=bandwidth,
                                min_effective=min_effective)
    fitted_mean = np.column_stack([
        np.interp(z, mean_atlas['grid'], mean_atlas['mean'][:, j])
        for j in range(y.shape[1])
    ])
    residual = y - fitted_mean
    estimator = LedoitWolf(assume_centered=True).fit(residual)
    covariance = (estimator.covariance_ + estimator.covariance_.T) * .5
    values, vectors = np.linalg.eigh(covariance)
    floor = max(float(np.trace(covariance)) / y.shape[1] * 1e-6, 1e-8)
    values = np.maximum(values, floor)
    covariance = (vectors * values) @ vectors.T
    covariance = (covariance + covariance.T) * .5
    precision = (vectors * (1. / values)) @ vectors.T
    logdet = float(np.log(values).sum())
    result = dict(mean_atlas)
    result.update(covariance=covariance, precision=precision,
                  covariance_shrinkage=float(estimator.shrinkage_),
                  covariance_eigenvalue_floor=floor, logdet=logdet,
                  feature_count=y.shape[1], residual_centered=False)
    return result


def project_residual_metric_atlas(atlas, Y):
    """Project profiles against the frozen mean path with global Mahalanobis distance."""
    y = _matrix(Y, 'Y')
    mean = np.asarray(atlas['mean'], dtype=float)
    grid = _grid(atlas['grid'])
    covariance = np.asarray(atlas['covariance'], dtype=float)
    precision = np.asarray(atlas['precision'], dtype=float)
    if y.shape[1] != mean.shape[1] or atlas.get('feature_count') != y.shape[1]:
        raise ValueError('Y feature count does not match fitted residual metric atlas')
    if mean.shape != (len(grid), y.shape[1]):
        raise ValueError('fitted mean path has invalid dimensions')
    if covariance.shape != (y.shape[1], y.shape[1]) or precision.shape != covariance.shape:
        raise ValueError('fitted residual covariance or precision has invalid dimensions')
    if not np.isfinite(mean).all() or not np.isfinite(covariance).all() or not np.isfinite(precision).all():
        raise ValueError('fitted atlas contains non-finite values')
    logdet = float(atlas['logdet'])
    if not np.isfinite(logdet):
        raise ValueError('fitted covariance log determinant must be finite')

    n, d = y.shape
    n_grid = len(grid)
    log_density = np.empty((n, n_grid), dtype=float)
    for start in range(0, n, 256):
        block = y[start:start + 256]
        for k, center in enumerate(mean):
            residual = block - center
            d2 = np.einsum('ij,jk,ik->i', residual, precision, residual, optimize=True)
            log_density[start:start + len(block), k] = -.5 * (
                np.maximum(d2, 0.) + logdet + d * np.log(2. * np.pi))
    log_posterior = log_density - logsumexp(log_density, axis=1, keepdims=True)
    posterior = np.exp(log_posterior)
    posterior /= posterior.sum(axis=1, keepdims=True)
    map_index = np.argmax(log_density, axis=1)
    entropy = -np.sum(xlogy(posterior, posterior), axis=1)
    residual_norm = np.sqrt(np.maximum(-2. * (log_density[np.arange(n), map_index]
        + .5 * (logdet + d * np.log(2. * np.pi))), 0.))
    return {
        'z_MAP': grid[map_index],
        'z_mean': np.clip(posterior @ grid, 0., 1.),
        'posterior': posterior,
        'entropy': entropy,
        'residual_norm': residual_norm,
        'log_density': logsumexp(log_density, axis=1) - np.log(n_grid),
        'log_density_at_MAP': log_density[np.arange(n), map_index],
    }
