"""Tubule-to-tubule variability per gene, beyond the mean (notebook 107; docs/results/pt-injury-heterogeneity.md).

* ``moment_dispersion``: per gene, the NB method-of-moments dispersion over a set of units with an exposure (library)
  offset: alpha = sum[(y - mu)^2 - mu] / sum[mu^2], mu_u = N_u * sum(y) / sum(N); floored.
* ``residual_dispersion``: log alpha minus a lowess trend on log mean count across genes (the residual
  over-dispersion of Eling et al. 2018, here by moments).
* ``detection``: share of units with a nonzero count, per gene.

Ponytail: moments, not a likelihood; adequate at the >= 30 units and >= 0.5 mean counts the notebook requires. The
upgrade is an NB2 location-scale fit with a unit-size term (alpha = tau^2 + kappa / m).
"""
from __future__ import annotations

import numpy as np
from scipy import sparse
from statsmodels.nonparametric.smoothers_lowess import lowess

FLOOR = 1e-3


def moment_dispersion(counts, library, *, floor=FLOOR):
    """(alpha, mean count) per gene for units x genes ``counts`` with per-unit ``library``."""
    y = sparse.csr_matrix(counts, dtype=float)
    n = np.asarray(library, float)
    share = np.asarray(y.sum(axis=0)).ravel() / n.sum()                        # rate per library unit
    mu_sq = (n ** 2).sum() * share ** 2                                         # sum_u mu_u^2
    y_sq = np.asarray(y.multiply(y).sum(axis=0)).ravel()
    y_mu = np.asarray((sparse.diags(n) @ y).sum(axis=0)).ravel() * share        # sum_u y_u mu_u
    sum_mu = n.sum() * share
    ss = y_sq - 2 * y_mu + mu_sq                                                # sum (y - mu)^2
    alpha = np.divide(ss - sum_mu, mu_sq, out=np.full(len(share), np.nan), where=mu_sq > 0)
    mean = np.asarray(y.mean(axis=0)).ravel()
    return np.maximum(alpha, floor), mean


def residual_dispersion(alpha, mean, *, frac=.3):
    """(epsilon, trend): log alpha minus a lowess trend of log alpha on log mean, over finite genes."""
    a, m = np.log(np.asarray(alpha, float)), np.log(np.asarray(mean, float))
    ok = np.isfinite(a) & np.isfinite(m)
    trend = np.full(len(a), np.nan)
    if ok.sum() >= 10:
        fit = lowess(a[ok], m[ok], frac=frac, return_sorted=False)
        trend[ok] = fit
    return a - trend, trend


def detection(counts):
    """Share of units with a nonzero count, per gene."""
    y = sparse.csr_matrix(counts)
    return np.asarray((y > 0).mean(axis=0)).ravel()
