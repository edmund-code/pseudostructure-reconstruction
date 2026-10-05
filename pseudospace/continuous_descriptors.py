"""What a continuous PT coordinate adds beyond S1/S2/S3 segment means (notebooks 56-60).

Segment pseudobulk is a projection of a gene's positional profile onto segment-wise constants. Two
things lie outside that projection: structure inside a segment (gradients, interior extrema, where a
switch happens) and the precise location of a switch relative to a segment boundary. These helpers
measure such features per specimen so that they can be replicated between specimens, between
coordinates and in external data. They never select genes by a species contrast.

Two artefacts can create apparent within-segment structure: segment labels that disagree with the
coordinate (label noise) and positional noise in the coordinate (smears a step into a ramp).
``label_agreeing_thresholds`` and ``blockwise_permutation`` build the semi-synthetic nulls for both.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from .stats_gam import make_gam_design

SEGMENTS = ('PT-S1', 'PT-S2', 'PT-S3')


def label_agreeing_thresholds(position, segment, order=SEGMENTS, max_candidates=600):
    """Two cut points t12 < t23 on one specimen's coordinate that best reproduce its ordered labels.

    Maximises the number of structures with label order[0] below t12, order[1] in [t12, t23) and
    order[2] at or above t23. Candidates are midpoints between sorted positions, thinned to at most
    ``max_candidates`` quantiles (ceiling: a cut can be off by one candidate spacing). Returns
    (t12, t23, agreement fraction).
    """
    s = np.asarray(position, float)
    lab = np.asarray(segment).astype(str)
    if s.shape != lab.shape or s.size < 3 or not np.isfinite(s).all():
        raise ValueError('Need aligned finite positions and labels.')
    code = np.full(len(lab), -1)
    for k, name in enumerate(order):
        code[lab == name] = k
    if (code < 0).any() or len(np.unique(code)) < 3:
        raise ValueError('Every structure needs one of the three ordered labels, and all three must occur.')
    sort = np.argsort(s, kind='stable')
    s_sorted, c_sorted = s[sort], code[sort]
    mids = (s_sorted[1:] + s_sorted[:-1]) / 2
    mids = mids[np.diff(s_sorted) > 0]
    cuts = np.unique(np.r_[s_sorted[0] - 1e-9, mids, s_sorted[-1] + 1e-9])
    if len(cuts) > max_candidates:
        cuts = np.unique(np.quantile(cuts, np.linspace(0, 1, max_candidates)))
    # below[k][i] = structures of class k with position < cuts[i]
    below = np.stack([np.searchsorted(s_sorted[c_sorted == k], cuts, side='left') for k in range(3)])
    n_last = int((c_sorted == 2).sum())
    # agreement(i, j) = S1 below cut i + S2 in [cut i, cut j) + S3 at or above cut j, for i <= j
    first = below[0][:, None]
    middle = below[1][None, :] - below[1][:, None]
    last = (n_last - below[2])[None, :]
    score = first + middle + last
    score[np.tril_indices(len(cuts), -1)] = -1
    i, j = np.unravel_index(np.argmax(score), score.shape)
    return float(cuts[i]), float(cuts[j]), float(score[i, j] / len(s))


def coordinate_segments(position, segment, specimen, order=SEGMENTS):
    """Per-specimen coordinate-threshold segments: labels redrawn by the coordinate alone.

    Returns (labels named like ``order``, table of specimen thresholds and label agreement).
    """
    s = np.asarray(position, float)
    lab, spec = np.asarray(segment).astype(str), np.asarray(specimen).astype(str)
    out = np.empty(len(s), dtype=object)
    rows = []
    for name in np.unique(spec):
        m = spec == name
        t12, t23, agree = label_agreeing_thresholds(s[m], lab[m], order)
        out[m] = np.where(s[m] < t12, order[0], np.where(s[m] < t23, order[1], order[2]))
        rows.append({'specimen': name, 't12': t12, 't23': t23, 'label agreement': agree})
    return out.astype(str), pd.DataFrame(rows)


def segment_units(position, t12, t23, lo=0., hi=1.):
    """Map positions to segment units: [lo, t12) -> [0, 1), [t12, t23) -> [1, 2), [t23, hi] -> [2, 3]."""
    s = np.asarray(position, float)
    if not lo < t12 < t23 < hi:
        raise ValueError('Need lo < t12 < t23 < hi.')
    return np.where(s < t12, (s - lo) / (t12 - lo),
                    np.where(s < t23, 1 + (s - t12) / (t23 - t12), 2 + (s - t23) / (hi - t23)))


def logistic_transitions(position, segment, specimen, order=SEGMENTS):
    """Per-specimen S1->S2 and S2->S3 transitions by logistic regression on position (notebook 37)."""
    from sklearn.linear_model import LogisticRegression

    s = np.asarray(position, float)
    lab, spec = np.asarray(segment).astype(str), np.asarray(specimen).astype(str)
    rows = []
    for name in np.unique(spec):
        for first, second in zip(order[:-1], order[1:]):
            keep = (spec == name) & np.isin(lab, [first, second])
            model = LogisticRegression(C=1e4).fit(s[keep][:, None], (lab[keep] == second).astype(int))
            rows.append({'specimen': name, 'transition': f'{first[-2:]}→{second[-2:]}',
                         'position': float(-model.intercept_[0] / model.coef_[0, 0])})
    return pd.DataFrame(rows)


def blockwise_permutation(blocks, rng):
    """Row index that permutes rows only within each block (same permutation for every gene).

    Applied as ``Y[index]`` with positions and labels kept, it removes every dependence of
    expression on position and on labels inside a block while keeping block means, the joint
    distribution across genes and therefore inter-gene correlation.
    """
    blocks = np.asarray(blocks)
    index = np.arange(len(blocks))
    for name in np.unique(blocks):
        rows = np.flatnonzero(blocks == name)
        index[rows] = rows[rng.permutation(len(rows))]
    return index


def specimen_curves(expression, position, specimen, grid, knots):
    """Unweighted regression-spline fit per specimen, evaluated on a shared grid.

    Returns {specimen: array (genes, len(grid))}. ``expression`` is structures x genes (dense or
    sparse). The basis is notebook 12's cubic B-spline with shared knots.
    """
    from scipy import sparse

    s = np.asarray(position, float)
    spec = np.asarray(specimen).astype(str)
    grid_design = make_gam_design(np.asarray(grid, float), knots)
    out = {}
    for name in np.unique(spec):
        m = spec == name
        design = make_gam_design(s[m], knots)
        block = expression[m]
        block = block.toarray() if sparse.issparse(block) else np.asarray(block, float)
        coef, *_ = np.linalg.lstsq(design, block, rcond=None)
        out[name] = (grid_design @ coef).T
    return out


def _crossing(u, y, level):
    """First position where an increasing-on-average curve reaches ``level`` (linear interpolation)."""
    above = np.flatnonzero(y >= level)
    if not len(above):
        return np.nan
    j = above[0]
    if j == 0:
        return float(u[0])
    y0, y1 = y[j - 1], y[j]
    frac = 0. if y1 == y0 else (level - y0) / (y1 - y0)
    return float(u[j - 1] + frac * (u[j] - u[j - 1]))


def monotone_onsets(curves, units, *, min_amplitude=.5, min_abs_rho=.9):
    """Onset (half-amplitude crossing) and width (10-90%) for monotone-like curves.

    ``curves`` is genes x grid, ``units`` the grid in segment units (increasing). A curve is
    monotone-like when its amplitude (max - min) reaches ``min_amplitude`` and its Spearman
    correlation with the grid reaches ``min_abs_rho`` in absolute value. Falling curves are flipped,
    so the onset is where the change is half done in either direction.
    """
    y = np.asarray(curves, float)
    u = np.asarray(units, float)
    if y.ndim != 2 or y.shape[1] != len(u) or np.any(np.diff(u) <= 0):
        raise ValueError('Need genes x grid curves and increasing grid units.')
    grid_rank = rankdata(u)
    amp = y.max(axis=1) - y.min(axis=1)
    ry = rankdata(y, axis=1)
    rc = ry - ry.mean(axis=1, keepdims=True)
    gc = grid_rank - grid_rank.mean()
    rho = (rc @ gc) / np.maximum(np.linalg.norm(rc, axis=1) * np.linalg.norm(gc), 1e-12)
    monotone = (amp >= min_amplitude) & (np.abs(rho) >= min_abs_rho)
    onset, width = np.full(len(y), np.nan), np.full(len(y), np.nan)
    for i in np.flatnonzero(monotone):
        z = y[i] * np.sign(rho[i])
        z = (z - z.min()) / (z.max() - z.min())
        onset[i] = _crossing(u, z, .5)
        width[i] = _crossing(u, z, .9) - _crossing(u, z, .1)
    return pd.DataFrame({'amplitude': amp, 'rho': rho, 'monotone': monotone,
                         'direction': np.sign(rho), 'onset': onset, 'width': width})


def segment_means(expression, segment, order=SEGMENTS):
    """Mean expression per segment (genes x 3) for one specimen's rows."""
    from scipy import sparse

    lab = np.asarray(segment).astype(str)
    cols = []
    for name in order:
        block = expression[lab == name]
        cols.append(np.asarray(block.mean(axis=0)).ravel() if sparse.issparse(block) else np.asarray(block).mean(axis=0))
    return np.column_stack(cols)


def switch_fraction(means, min_range=.1):
    """Discrete onset summary r = (m2 - m1) / (m3 - m1) from three segment means.

    For a monotone switch, larger r means an earlier switch, in either direction of change. It
    confounds switch location with switch width. Genes with |m3 - m1| < ``min_range`` get NaN.
    """
    m = np.asarray(means, float)
    if m.ndim != 2 or m.shape[1] != 3:
        raise ValueError('Need genes x 3 segment means.')
    span = m[:, 2] - m[:, 0]
    r = (m[:, 1] - m[:, 0]) / np.where(np.abs(span) >= min_range, span, np.nan)
    return r


def partial_spearman(x, y, covariates):
    """Rank-based partial correlation of x and y given covariates (Pearson on rank residuals)."""
    x, y = rankdata(np.asarray(x, float)), rankdata(np.asarray(y, float))
    z = np.asarray(covariates, float)
    z = z[:, None] if z.ndim == 1 else z
    if not (len(x) == len(y) == len(z)) or len(x) < z.shape[1] + 3:
        raise ValueError('Need aligned x, y and covariates with enough rows.')
    design = np.column_stack([np.ones(len(x)), np.column_stack([rankdata(c) for c in z.T])])
    rx = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
    ry = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    return float(rx @ ry / np.sqrt((rx @ rx) * (ry @ ry)))


def bootstrap_partial_spearman(x, y, covariates, *, n_boot=2000, seed=0, level=.95):
    """Estimate with a percentile bootstrap over rows (genes treated as independent)."""
    x, y, z = np.asarray(x, float), np.asarray(y, float), np.asarray(covariates, float)
    estimate = partial_spearman(x, y, z)
    rng = np.random.default_rng(seed)
    draws = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, len(x), len(x))
        draws[b] = partial_spearman(x[idx], y[idx], z[idx])
    tail = (1 - level) / 2
    return estimate, float(np.nanquantile(draws, tail)), float(np.nanquantile(draws, 1 - tail))


def group_difference_curves(expression, position, group, specimen, knots, grid, *, nuisance_shape=None):
    """Group-1 minus group-0 smooth difference Δ(s) on ``grid`` from notebook 31's full design.

    Uses ``pathway_calibration.contrast_designs`` (specimen contrasts, equal specimen weights, and an
    optional nuisance shape such as species in a relabeled comparison) and returns genes x grid.
    """
    from scipy import sparse

    from .pathway_calibration import contrast_designs

    designs, weights = contrast_designs(position, group, specimen, knots, nuisance_shape=nuisance_shape)
    x = designs['full']
    p_b = make_gam_design(np.asarray(grid[:1], float), knots).shape[1]
    root = np.sqrt(weights)
    y = expression.toarray() if sparse.issparse(expression) else np.asarray(expression, float)
    beta, *_ = np.linalg.lstsq(root[:, None] * x, root[:, None] * y, rcond=None)
    basis = make_gam_design(np.asarray(grid, float), knots)[:, 1:]
    # Column layout of the full design: [intercept, basis (p_b - 1), group, group x basis, ...].
    level, shape = beta[p_b], beta[p_b + 1:2 * p_b]
    return (level[None, :] + basis @ shape).T


def within_segment_slopes(expression, position, specimen, segment, order=SEGMENTS, min_rows=30):
    """OLS slope t-statistic of each gene on position inside each specimen x label segment.

    Returns {specimen: DataFrame genes x segments}. A segment with fewer than ``min_rows``
    structures gets NaN. Structures are treated as independent (descriptive t).
    """
    from scipy import sparse

    s = np.asarray(position, float)
    spec, lab = np.asarray(specimen).astype(str), np.asarray(segment).astype(str)
    out = {}
    for name in np.unique(spec):
        cols = {}
        for seg in order:
            rows = (spec == name) & (lab == seg)
            if rows.sum() < min_rows:
                cols[seg] = np.full(expression.shape[1], np.nan)
                continue
            y = expression[rows]
            y = y.toarray() if sparse.issparse(y) else np.asarray(y, float)
            x = s[rows] - s[rows].mean()
            sxx = x @ x
            yc = y - y.mean(axis=0)
            slope = (x @ yc) / sxx
            resid = np.maximum((yc ** 2).sum(axis=0) - slope ** 2 * sxx, 0) / (rows.sum() - 2)
            cols[seg] = slope / np.sqrt(np.maximum(resid / sxx, 1e-30))
        out[name] = pd.DataFrame(cols)
    return out


def interior_extrema(curves, units, *, min_amplitude=.5, min_prominence=.2, inner=.2):
    """Segment holding a curve's interior maximum or minimum, or -1.

    The global extremum (maximum, else minimum) must lie in a segment's inner part (at least
    ``inner`` of the segment length from both boundaries) and exceed both of that segment's boundary
    values by ``min_prominence``; the curve's amplitude must reach ``min_amplitude``. Returns
    (segment index or -1, sign +1 for a peak / -1 for a trough).
    """
    y = np.asarray(curves, float)
    u = np.asarray(units, float)
    seg_of = np.full(len(y), -1)
    sign_of = np.zeros(len(y), int)
    amp = y.max(axis=1) - y.min(axis=1)
    for i in np.flatnonzero(amp >= min_amplitude):
        for sign in (1, -1):
            z = sign * y[i]
            j = int(np.argmax(z))
            k = int(np.floor(u[j]))
            frac = u[j] - k
            if k > 2 or frac < inner or frac > 1 - inner:
                continue
            ends = np.interp([k, k + 1], u, z)
            if z[j] - ends.max() >= min_prominence:
                seg_of[i], sign_of[i] = k, sign
                break
    return seg_of, sign_of
