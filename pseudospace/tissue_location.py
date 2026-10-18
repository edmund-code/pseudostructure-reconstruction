"""Nephron position against tissue location for unit-level injury (notebook 96; docs/results/pt-injury-tissue-location.md).

* ``program_score``: per unit, the mean over a program's genes of (expression - reference mean) / max(reference SD,
  floor), with the mean and SD over reference units.
* ``position_expectation``: the reference units' binned mean score, interpolated at each unit's position; with
  ``cross_fit`` a reference unit's expectation comes from the other reference specimens only.
* ``glomerular_depth``: mean distance to the k nearest glomerulus centroids, and its percentile within the section.
* ``equal_count_bins``, ``spline_basis``: position strata for nulls, and a natural cubic spline design in position.
* ``depth_effect``: OLS of the effect on spline(x) + depth (Frisch-Waugh), with a null permuting depth within position
  bins; returns beta (scaled to the 0-1 depth range), partial R2 and the null betas.
* ``knn_weights``, ``morans_i``, ``moran_null``: row-standardised k-nearest-neighbour weights, Moran's I and its null
  permuting values within position bins.
* ``mixture_bic``: BIC of one against two Gaussian components.
"""
from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree

from .regulator_activity import stratified_permutations


def program_score(expression, reference, *, floor=.25):
    """Units x genes expression -> per unit mean of genes standardised on the ``reference`` units (SD floored)."""
    X = np.asarray(expression.toarray() if sparse.issparse(expression) else expression, float)
    ref = np.asarray(reference, bool)
    mu = X[ref].mean(axis=0)
    sd = np.maximum(X[ref].std(axis=0, ddof=1), floor)
    return ((X - mu) / sd).mean(axis=1)


def _binned_curve(x, y, edges):
    idx = np.clip(np.searchsorted(edges, x, side='right') - 1, 0, len(edges) - 2)
    centres, means = [], []
    for b in range(len(edges) - 1):
        m = idx == b
        if m.any():
            centres.append(.5 * (edges[b] + edges[b + 1]))
            means.append(float(np.mean(y[m])))
    return np.array(centres), np.array(means)


def position_expectation(x, score, specimen, reference_specimens, *, width=.05, cross_fit=True, builders=None):
    """Expected score at each unit's position from reference units' binned means (linear between bin centres).

    Non-reference units use all reference specimens; with ``cross_fit`` a reference unit uses the other reference
    specimens only. ``builders`` (bool per unit) limits which reference units build the curves (default: all).
    """
    x, score, specimen = np.asarray(x, float), np.asarray(score, float), np.asarray(specimen)
    build = np.ones(len(x), bool) if builders is None else np.asarray(builders, bool)
    refs = list(reference_specimens)
    edges = np.arange(0, 1 + width, width)
    edges[-1] = max(edges[-1], x.max() + 1e-9)
    out = np.full(len(x), np.nan)
    ref_mask = np.isin(specimen, refs)
    c, m = _binned_curve(x[ref_mask & build], score[ref_mask & build], edges)
    out[~ref_mask] = np.interp(x[~ref_mask], c, m)
    for s in refs:
        own = specimen == s
        others = (np.isin(specimen, [r for r in refs if r != s]) if cross_fit else ref_mask) & build
        c, m = _binned_curve(x[others], score[others], edges)
        out[own] = np.interp(x[own], c, m)
    return out


def glomerular_depth(units_xy, glomeruli_xy, *, k=3):
    """(mean distance to the k nearest glomeruli, its percentile among the units given)."""
    d, _ = cKDTree(np.asarray(glomeruli_xy, float)).query(np.asarray(units_xy, float), k=k)
    d = d.mean(axis=1) if d.ndim == 2 else d
    pct = (np.argsort(np.argsort(d)) + .5) / len(d)
    return d, pct


def equal_count_bins(x, n=20):
    """Bin index 0..n-1 of each value, by equal-count quantiles."""
    ranks = np.argsort(np.argsort(np.asarray(x, float)))
    return np.minimum((ranks * n) // len(ranks), n - 1)


def spline_basis(x, df=5):
    """Intercept plus a natural cubic regression spline of x with ``df`` degrees of freedom (patsy's ``cr``)."""
    from patsy import dmatrix

    return np.asarray(dmatrix(f'cr(x, df={df})', {'x': np.asarray(x, float)}, return_type='matrix'))


def _residualise(X, Y):
    coef, *_ = np.linalg.lstsq(X, Y, rcond=None)
    return Y - X @ coef


def depth_effect(effect, x, depth, strata, *, n_perm=1000, seed=96, df=5):
    """beta of depth (0-1 scale) on the effect at fixed spline(x), its partial R2 and the null betas."""
    e, x, d = (np.asarray(v, float) for v in (effect, x, depth))
    X = spline_basis(x, df)
    e_r = _residualise(X, e)
    d_r = _residualise(X, d)
    beta = float(d_r @ e_r / (d_r @ d_r))
    partial = float((d_r @ e_r) ** 2 / ((d_r @ d_r) * (e_r @ e_r)))
    perms = stratified_permutations(strata, n_perm, seed)
    D = _residualise(X, d[perms].T)
    null = (D * e_r[:, None]).sum(axis=0) / (D ** 2).sum(axis=0)
    fitted_residual = e_r - beta * d_r
    return {'beta': beta, 'partial R2': partial, 'null': null,
            'p': float((1 + (np.abs(null) >= abs(beta)).sum()) / (1 + n_perm)), 'residual': fitted_residual}


def knn_weights(xy, k=6):
    """Row-standardised binary weights to the k nearest other units."""
    xy = np.asarray(xy, float)
    _, idx = cKDTree(xy).query(xy, k=k + 1)
    rows = np.repeat(np.arange(len(xy)), k)
    return sparse.csr_matrix((np.full(rows.size, 1. / k), (rows, idx[:, 1:].ravel())), shape=(len(xy), len(xy)))


def morans_i(values, W):
    """Moran's I with row-standardised W (S0 = n): z'Wz / z'z."""
    z = np.asarray(values, float)
    z = z - z.mean()
    return float(z @ (W @ z) / (z @ z))


def moran_null(values, W, strata, *, n_perm=1000, seed=96):
    """Moran's I over permutations of the values within strata."""
    z = np.asarray(values, float)
    perms = stratified_permutations(strata, n_perm, seed)
    Z = z[perms].T
    Z = Z - Z.mean(axis=0)
    return np.asarray((Z * (W @ Z)).sum(axis=0) / (Z ** 2).sum(axis=0)).ravel()


def mixture_bic(values, *, seed=96):
    """BIC of one and two Gaussian components, the two-component weights and means (sorted by mean)."""
    from sklearn.mixture import GaussianMixture

    v = np.asarray(values, float).reshape(-1, 1)
    g1 = GaussianMixture(1, random_state=seed).fit(v)
    g2 = GaussianMixture(2, random_state=seed, n_init=5).fit(v)
    order = np.argsort(g2.means_.ravel())
    return {'bic 1': float(g1.bic(v)), 'bic 2': float(g2.bic(v)), 'weights': g2.weights_[order].tolist(),
            'means': g2.means_.ravel()[order].tolist()}
