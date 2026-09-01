"""GAM (P-spline) design/penalty utilities and generic statistics helpers.

Extracted verbatim from ``6_mouse_only_pseudospace.ipynb`` (Section 4.0), which
copied them from ``5_mouse_healthy_vs_aki.ipynb``. ``single_smooth_slice``,
``sparse_gam_sse`` and ``dense_gam_sse`` are unused by notebook 6 but retained
here because they are live utilities in notebooks 4/5.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.stats import spearmanr
from patsy import dmatrix
from statsmodels.stats.multitest import multipletests


def as_csr(matrix):
    if sp.issparse(matrix):
        return matrix.tocsr().astype(np.float64)
    return sp.csr_matrix(np.asarray(matrix, dtype=np.float64))


def zscore_rows(values):
    values = np.asarray(values, dtype=float)
    mean = np.nanmean(values, axis=1, keepdims=True)
    std = np.nanstd(values, axis=1, keepdims=True)
    std[~np.isfinite(std) | (std == 0)] = 1.0
    return (values - mean) / std


def make_gam_design(x, knots):
    smooth = np.asarray(
        dmatrix(
            'bs(x, knots=knots, degree=3, include_intercept=False, '
            'lower_bound=0, upper_bound=1) - 1',
            {'x': np.asarray(x, dtype=float), 'knots': np.asarray(knots, dtype=float)},
        )
    )
    return np.column_stack([np.ones(len(x)), smooth])


def gam_internal_knots(x, basis_df=6):
    degree = 3
    n_internal = max(int(basis_df) - degree, 0)
    if n_internal == 0:
        return np.array([], dtype=float)
    values = np.clip(np.asarray(x, dtype=float), 0, 1)
    probs = np.linspace(0, 1, n_internal + 2)[1:-1]
    knots = np.unique(np.clip(np.quantile(values, probs), 1e-6, 1 - 1e-6))
    if knots.size != n_internal:
        knots = np.linspace(0, 1, n_internal + 2)[1:-1]
    return knots


def second_difference_penalty(n_coefficients):
    if n_coefficients <= 2:
        return np.eye(n_coefficients)
    diff = np.diff(np.eye(n_coefficients), n=2, axis=0)
    return diff.T @ diff


def make_gam_penalty(n_columns, smooth_slices, n_observations, alpha=1.0):
    penalty = np.zeros((n_columns, n_columns), dtype=float)
    for smooth_slice in smooth_slices:
        start, stop, step = smooth_slice.indices(n_columns)
        if step != 1:
            raise ValueError('GAM smooth slices must be contiguous')
        n_smooth = stop - start
        if n_smooth <= 0:
            continue
        scale = float(alpha) * max(int(n_observations), 1) / max(n_smooth, 1)
        penalty[start:stop, start:stop] += scale * second_difference_penalty(n_smooth)
    return penalty


def single_smooth_slice(design):
    return [slice(1, design.shape[1])]


def sparse_gam_sse(matrix, design, smooth_slices=None):
    values = as_csr(matrix)
    x = np.asarray(design, dtype=float)
    smooth_slices = [] if smooth_slices is None else smooth_slices
    xtx = x.T @ x
    penalty = make_gam_penalty(x.shape[1], smooth_slices, x.shape[0])
    system_inv = np.linalg.pinv(xtx + penalty)
    xty = np.asarray(x.T @ values)
    beta = np.asarray(system_inv @ xty)
    yty = np.asarray(values.multiply(values).sum(axis=0)).ravel()
    fitted_cross = np.sum(beta * xty, axis=0)
    model_cross = np.sum(beta * (xtx @ beta), axis=0)
    sse = yty - 2 * fitted_cross + model_cross
    edf = float(np.trace(system_inv @ xtx))
    return np.maximum(sse, 0.0), beta, edf


def dense_gam_sse(values, design, smooth_slices=None):
    values = np.asarray(values, dtype=np.float64)
    x = np.asarray(design, dtype=np.float64)
    smooth_slices = [] if smooth_slices is None else smooth_slices
    xtx = x.T @ x
    penalty = make_gam_penalty(x.shape[1], smooth_slices, x.shape[0])
    system_inv = np.linalg.pinv(xtx + penalty)
    xty = x.T @ values
    beta = system_inv @ xty
    yty = np.sum(values * values, axis=0)
    fitted_cross = np.sum(beta * xty, axis=0)
    model_cross = np.sum(beta * (xtx @ beta), axis=0)
    sse = yty - 2 * fitted_cross + model_cross
    edf = float(np.trace(system_inv @ xtx))
    return np.maximum(sse, 0.0), beta, edf


def bh_adjust(pvalues):
    pvalues = np.asarray(pvalues, dtype=float)
    adjusted = np.ones_like(pvalues)
    valid = np.isfinite(pvalues)
    adjusted[valid] = multipletests(pvalues[valid], method='fdr_bh')[1]
    return adjusted


def safe_spearman(left, right):
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    valid = np.isfinite(left) & np.isfinite(right)
    if valid.sum() < 3 or np.ptp(left[valid]) == 0 or np.ptp(right[valid]) == 0:
        return np.nan
    return spearmanr(left[valid], right[valid]).correlation


def resolve_present(requested, universe):
    lookup = {str(gene).upper(): gene for gene in universe}
    return [lookup[str(gene).upper()] for gene in requested if str(gene).upper() in lookup]
