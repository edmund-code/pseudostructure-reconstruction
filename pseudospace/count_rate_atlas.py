"""Frozen conditional-count atlas with equal-specimen spatial rates.

This is a projection pilot: gene counts are independent negative-binomial
observations conditional on all-gene library exposure and a finite spatial
grid. Its fixed dispersion and posterior are approximations, not calibrated
uncertainty or anatomical ground truth.
"""
from __future__ import annotations

import numpy as np
from scipy.special import gammaln, logsumexp, xlogy

from pseudospace.canonical_mean_atlas import (
    _bandwidth, _coordinates, _grid, _kernel_weights,
)
from pseudospace.count_representation import _counts_matrix, _exposure, _labels


def _theta(value):
    if not np.isscalar(value):
        raise ValueError('theta must be positive finite or positive infinity')
    try:
        theta = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError('theta must be positive finite or positive infinity') from exc
    if np.isnan(theta) or theta <= 0:
        raise ValueError('theta must be positive finite or positive infinity')
    return theta


def _positive_finite(value, name):
    if not np.isscalar(value):
        raise ValueError(f'{name} must be positive and finite')
    try:
        value = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f'{name} must be positive and finite') from exc
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f'{name} must be positive and finite')
    return value


def fit_count_rate_atlas(counts, library, z, specimens, grid=None,
                         bandwidth=.08, theta=100., min_effective=12.,
                         rate_floor=1e-12):
    """Fit equal-specimen conditional gene rates over a supplied training axis."""
    counts = _counts_matrix(counts, 'counts')
    library = _exposure(library, counts, 'library')
    z = _coordinates(z, len(counts))
    specimens = _labels(specimens, len(counts))
    grid = _grid(grid)
    bandwidth = _bandwidth(bandwidth)
    theta = _theta(theta)
    min_effective = _positive_finite(min_effective, 'min_effective')
    rate_floor = _positive_finite(rate_floor, 'rate_floor')
    if rate_floor >= 1 / counts.shape[1]:
        raise ValueError('rate_floor must be smaller than 1 / number of genes')

    training_specimens = np.unique(specimens)
    by_specimen = []
    support = np.empty((len(grid), len(training_specimens)), dtype=float)
    for sidx, specimen in enumerate(training_specimens):
        use = specimens == specimen
        weights = _kernel_weights(grid, z[use], bandwidth)
        weights /= weights.sum(axis=1, keepdims=True)
        n_eff = 1. / np.sum(weights**2, axis=1)
        support[:, sidx] = n_eff
        if np.any(n_eff < min_effective):
            j = int(np.flatnonzero(n_eff < min_effective)[0])
            raise ValueError(f'grid position {grid[j]:g}, specimen {specimen!r} has '
                             f'effective support {n_eff[j]:.2f}; requires {min_effective}')
        numerator = weights @ counts[use]
        denominator = weights @ library[use]
        by_specimen.append(numerator / denominator[:, None])
    rates = np.mean(by_specimen, axis=0)
    rates = np.maximum(rates, rate_floor)
    row_sums = rates.sum(axis=1)
    too_large = row_sums > 1.
    if np.any(too_large):
        rates[too_large] /= row_sums[too_large, None]
    return {'grid': grid, 'rates': rates, 'theta': theta,
            'n_effective': support, 'training_specimens': training_specimens,
            'feature_count': counts.shape[1], 'bandwidth': bandwidth,
            'rate_floor': float(rate_floor)}


def _nb_log_density(counts, library, rates, theta):
    """Full profile-by-grid independent NB/Poisson log likelihood."""
    n, genes = counts.shape
    n_grid = len(rates)
    logp = np.empty((n, n_grid), dtype=float)
    count_sum = counts.sum(axis=1)
    if np.isinf(theta):
        constant = count_sum * np.log(library) - gammaln(counts + 1).sum(axis=1)
        log_rates = np.log(rates)
        rate_sums = rates.sum(axis=1)
        for k in range(n_grid):
            logp[:, k] = counts @ log_rates[k] - library * rate_sums[k] + constant
    else:
        log_theta = np.log(theta)
        constant = (gammaln(counts + theta).sum(axis=1)
                    - genes*gammaln(theta) - gammaln(counts + 1).sum(axis=1)
                    + count_sum * (np.log(library) - log_theta))
        log_rates = np.log(rates)
        for k in range(n_grid):
            penalty = np.sum((counts + theta)
                             * np.log1p(library[:, None] * rates[k][None, :] / theta), axis=1)
            logp[:, k] = counts @ log_rates[k] - penalty + constant
    return logp


def project_count_rate_atlas(atlas, counts, library, block_size=128):
    """Project counts/exposures only onto a frozen finite-grid rate atlas."""
    counts = _counts_matrix(counts, 'counts', min_rows=1)
    if counts.shape[1] != atlas['feature_count']:
        raise ValueError('query count feature count does not match atlas')
    library = _exposure(library, counts, 'library')
    if (not isinstance(block_size, (int, np.integer))
            or isinstance(block_size, (bool, np.bool_)) or block_size < 1):
        raise ValueError('block_size must be a positive integer')
    posterior = np.empty((len(counts), len(atlas['grid'])), dtype=float)
    log_density = np.empty(len(counts), dtype=float)
    z_mean = np.empty(len(counts), dtype=float)
    z_map = np.empty(len(counts), dtype=float)
    entropy = np.empty(len(counts), dtype=float)
    log_prior = np.log(len(atlas['grid']))
    for start in range(0, len(counts), int(block_size)):
        stop = min(start + int(block_size), len(counts))
        logp = _nb_log_density(counts[start:stop], library[start:stop],
                               atlas['rates'], atlas['theta'])
        normalizer = logsumexp(logp, axis=1)
        log_post = logp - normalizer[:, None]
        p = np.exp(log_post)
        p /= p.sum(axis=1, keepdims=True)
        posterior[start:stop] = p
        log_density[start:stop] = normalizer - log_prior
        z_mean[start:stop] = np.clip(p @ atlas['grid'], atlas['grid'][0], atlas['grid'][-1])
        z_map[start:stop] = atlas['grid'][np.argmax(logp, axis=1)]
        entropy[start:stop] = -np.sum(xlogy(p, p), axis=1)
    return {'posterior': posterior, 'z_mean': z_mean, 'z_MAP': z_map,
            'entropy': entropy, 'log_density': log_density}
