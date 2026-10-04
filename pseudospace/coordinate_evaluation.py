"""Coordinate validation helpers for the cross-species PT axis (notebook 35).

Every metric is gauge-free where possible: Spearman correlations do not change under a monotone
reparameterization of the coordinate. Perturbations act on raw counts, so a physical section keeps
its identity while its measurement changes. Nothing here estimates a coordinate.
"""
from __future__ import annotations

from hashlib import sha256

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata


def sha256_folds(genes, n_folds=3):
    """Deterministic gene folds that do not depend on gene order or the Python hash seed."""
    if n_folds < 2:
        raise ValueError('Need at least two folds.')
    return np.array([int(sha256(str(g).encode()).hexdigest(), 16) % n_folds for g in genes], dtype=int)


def log_normalize(counts, library, target_sum=1e4):
    """log1p(count / library * target_sum), sparse in and out."""
    counts = sparse.csr_matrix(counts, dtype=float)
    library = np.asarray(library, float).ravel()
    if len(library) != counts.shape[0] or not np.isfinite(library).all() or (library <= 0).any():
        raise ValueError('Need one positive finite library size per row.')
    out = counts.multiply((target_sum / library)[:, None]).tocsr()
    out.data = np.log1p(out.data)
    return out


def group_rank_correlations(expression, z, groups, *, min_rows=10):
    """Spearman correlation of every gene with ``z`` inside each group (genes x groups).

    Ties take average ranks. A gene that is constant inside a group, or a group with fewer than
    ``min_rows`` rows or a constant ``z``, gives NaN for that entry.
    """
    x = sparse.csr_matrix(expression, dtype=float)
    z, groups = np.asarray(z, float).ravel(), np.asarray(groups).astype(str)
    if x.shape[0] != len(z) or len(groups) != len(z) or not np.isfinite(z).all():
        raise ValueError('Expression rows, coordinate and groups must align; coordinate must be finite.')
    labels = sorted(set(groups))
    out = np.full((x.shape[1], len(labels)), np.nan)
    for j, label in enumerate(labels):
        rows = np.flatnonzero(groups == label)
        if len(rows) < min_rows or np.ptp(z[rows]) == 0:
            continue
        zr = rankdata(z[rows])
        zr = (zr - zr.mean()) / np.linalg.norm(zr - zr.mean())
        # Ponytail: dense block per group (<= a few thousand rows); chunk genes if memory binds.
        values = rankdata(x[rows].toarray(), axis=0)
        values -= values.mean(axis=0)
        norm = np.linalg.norm(values, axis=0)
        with np.errstate(invalid='ignore', divide='ignore'):
            out[:, j] = np.where(norm > 0, (values.T @ zr) / norm, np.nan)
    return pd.DataFrame(out, columns=labels)


def concordance(values, reference):
    """Spearman correlation across genes, ignoring genes missing in either vector."""
    values, reference = np.asarray(values, float), np.asarray(reference, float)
    keep = np.isfinite(values) & np.isfinite(reference)
    if keep.sum() < 3:
        return np.nan, int(keep.sum())
    a, b = rankdata(values[keep]), rankdata(reference[keep])
    if np.ptp(a) == 0 or np.ptp(b) == 0:
        return np.nan, int(keep.sum())
    return float(np.corrcoef(a, b)[0, 1]), int(keep.sum())


def absorption_ratio(stat_original, stat_frozen, stat_refit, *, floor=.02):
    """How much of an injected change survives once the coordinate is refitted.

    ``stat_frozen`` uses perturbed data with the unperturbed coordinate, ``stat_refit`` the same
    data with the refitted coordinate. A ratio of 1 means the coordinate did not absorb the
    effect; 0 means it absorbed all of it. Genes whose frozen change is below ``floor`` are
    excluded from the per-gene ratios but kept in the slope through the origin.
    """
    original, frozen, refit = (np.asarray(v, float) for v in (stat_original, stat_frozen, stat_refit))
    if not (original.shape == frozen.shape == refit.shape):
        raise ValueError('Statistics must align.')
    d_frozen, d_refit = frozen - original, refit - original
    keep = np.isfinite(d_frozen) & np.isfinite(d_refit)
    big = keep & (np.abs(d_frozen) >= floor)
    ratios = d_refit[big] / d_frozen[big]
    denominator = np.sum(d_frozen[keep] ** 2)
    return {'median_ratio': float(np.median(ratios)) if big.any() else np.nan,
            'slope': float(np.sum(d_refit[keep] * d_frozen[keep]) / denominator) if denominator > 0 else np.nan,
            'n_genes': int(keep.sum()), 'n_ratio_genes': int(big.sum()),
            'median_abs_frozen_change': float(np.median(np.abs(d_frozen[keep]))) if keep.any() else np.nan}


def _check_counts(counts):
    counts = sparse.csr_matrix(counts)
    data = counts.data
    if (not np.isfinite(data).all() or (data < 0).any() or not np.allclose(data, np.round(data))):
        raise ValueError('Expected nonnegative integer counts.')
    return counts


def binomial_thin_to_depth(counts, library, target, rng):
    """Thin each row with p = min(1, target / library); rows at or below target stay unchanged."""
    counts = _check_counts(counts)
    library = np.asarray(library, float).ravel()
    if len(library) != counts.shape[0] or not np.isfinite(target) or target <= 0 or (library <= 0).any():
        raise ValueError('Need a positive target and one positive library size per row.')
    p = np.minimum(1.0, target / library)
    out = counts.astype(np.int64).tocsr(copy=True)
    row_of_entry = np.repeat(np.arange(out.shape[0]), np.diff(out.indptr))
    out.data = rng.binomial(out.data, p[row_of_entry])
    out.eliminate_zeros()
    return out


def poisson_count_split(counts, epsilon, rng):
    """Split counts into A ~ Binomial(x, epsilon) and B = x - A (count splitting / thinning).

    Under a Poisson model A and B are independent; overdispersed counts make that approximate.
    """
    counts = _check_counts(counts)
    if not 0 < epsilon < 1:
        raise ValueError('epsilon must lie strictly between 0 and 1.')
    a = counts.astype(np.int64).tocsr(copy=True)
    a.data = rng.binomial(a.data, epsilon)
    b = counts.astype(np.int64).tocsr() - a
    a.eliminate_zeros()
    b.eliminate_zeros()
    return a, b


def scale_counts(counts, rows, cols, log_factor, rng):
    """Multiply the expected counts of selected entries by exp(log_factor[row, col]).

    Factors above 1 add Poisson(x * (factor - 1)); factors below 1 binomially thin. Each entry
    keeps expectation x * factor given the observed x. ``log_factor`` has shape
    (len(rows), len(cols)). Unselected entries are returned unchanged.
    """
    counts = _check_counts(counts).astype(np.int64)
    rows, cols = np.asarray(rows, int), np.asarray(cols, int)
    log_factor = np.asarray(log_factor, float)
    if log_factor.shape != (len(rows), len(cols)) or not np.isfinite(log_factor).all():
        raise ValueError('log_factor must be finite with shape (len(rows), len(cols)).')
    if len(np.unique(rows)) != len(rows) or len(np.unique(cols)) != len(cols):
        raise ValueError('Rows and columns must be unique.')
    block = counts[rows][:, cols].toarray()
    factor = np.exp(log_factor)
    up = factor > 1
    new = block.copy()
    new[up] = block[up] + rng.poisson(block[up] * (factor[up] - 1))
    new[~up] = rng.binomial(block[~up], factor[~up])
    delta = new - block
    r, c = np.nonzero(delta)
    full = sparse.csr_matrix((delta[r, c], (rows[r], cols[c])), shape=counts.shape)
    out = (counts + full).tocsr()
    out.eliminate_zeros()
    if (out.data < 0).any():
        raise RuntimeError('Negative counts after scaling.')
    return out


def fold_noise_variance(coordinates, specimen, segment, species):
    """Positional noise per species from disagreement between gene-fold coordinates.

    For every fold pair and specimen x segment cell, Var(z_a - z_b) / 2 estimates the variance of
    one fold coordinate's panel noise if the two errors are independent. Cells are averaged with
    equal weight within each species. This captures gene-panel noise only, so it is a lower bound
    on total positional error.
    """
    names = list(coordinates)
    if len(names) < 2:
        raise ValueError('Need at least two fold coordinates.')
    specimen, segment, species = (np.asarray(v).astype(str) for v in (specimen, segment, species))
    cells = pd.DataFrame({'specimen': specimen, 'segment': segment, 'species': species})
    rows = []
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            diff = np.asarray(coordinates[a], float) - np.asarray(coordinates[b], float)
            frame = cells.assign(diff=diff)
            var = frame.groupby(['species', 'specimen', 'segment'], observed=True)['diff'].var(ddof=1) / 2
            rows.append(var.rename(f'{a}-{b}'))
    table = pd.concat(rows, axis=1)
    per_cell = table.mean(axis=1)
    return per_cell.groupby(level='species').mean().to_dict(), table


def reflect_unit(z):
    """Reflect values into [0, 1] (repeatedly, so any finite value lands inside)."""
    z = np.mod(np.asarray(z, float), 2.0)
    return np.where(z > 1, 2.0 - z, z)


def interpolate_curves(curves, grid, points):
    """Evaluate gene x grid curves at given coordinate points by linear interpolation."""
    curves, grid = np.asarray(curves, float), np.asarray(grid, float)
    return np.vstack([np.interp(points, grid, row) for row in curves])
