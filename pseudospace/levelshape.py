"""Section 4 nested level/shape GAM decomposition and honest 2v2 sample-permutation.

Extracted verbatim from ``6_mouse_only_pseudospace.ipynb`` (Section 4). Depends on
the spline design/penalty utilities in :mod:`pseudospace.stats_gam`.
"""
from __future__ import annotations

import itertools
import warnings

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .stats_gam import make_gam_design, make_gam_penalty, safe_spearman


def build_ls_designs(s, c, knots):
    base = make_gam_design(s, knots)                       # (n, p_b) = [intercept, basis]
    p_b = base.shape[1]
    k = p_b - 1
    basis = base[:, 1:]
    c = np.asarray(c, dtype=float).reshape(-1)
    x0 = base
    x1 = np.column_stack([base, c])
    x2 = np.column_stack([base, c, c[:, None] * basis])
    sl0 = [slice(1, p_b)]
    sl1 = [slice(1, p_b)]
    sl2 = [slice(1, p_b), slice(p_b + 1, p_b + 1 + k)]
    return x0, x1, x2, sl0, sl1, sl2, p_b, k


def ls_grid_designs(grid, knots):
    base_g = make_gam_design(grid, knots)
    return base_g, base_g[:, 1:]


def fit_sse_grid(xtx, xty, yty, n_obs, smooth_slices, lambda_grid, keep_beta=False):
    p = xtx.shape[0]
    g = xty.shape[1]
    n_lam = len(lambda_grid)
    sse = np.empty((n_lam, g))
    edf = np.empty(n_lam)
    betas = [] if keep_beta else None
    for li, lam in enumerate(lambda_grid):
        penalty = make_gam_penalty(p, smooth_slices, n_obs, alpha=float(lam))
        ainv = np.linalg.pinv(xtx + penalty)
        beta = ainv @ xty
        fitted_cross = np.sum(beta * xty, axis=0)
        model_cross = np.sum(beta * (xtx @ beta), axis=0)
        sse[li] = np.maximum(yty - 2 * fitted_cross + model_cross, 0.0)
        edf[li] = float(np.trace(ainv @ xtx))
        if keep_beta:
            betas.append(beta)
    return sse, edf, betas


def m2_curves(beta2, base_g, basis_g, p_b):
    # beta2: (p2, g); columns [intercept, basis, c, c*basis]
    curve_h = base_g @ beta2[:p_b]                          # (G, g) condition c=0
    cond_coef = beta2[p_b]                                  # (g,)
    interaction = basis_g @ beta2[p_b + 1:]                 # (G, g)
    curve_a = curve_h + cond_coef[None, :] + interaction    # condition c=1
    return curve_h.T, curve_a.T


def shape_metrics(curve_h, curve_a):
    diff = curve_a - curve_h                                # (g, G)
    mean_gap = diff.mean(axis=1)
    shape_rms = np.sqrt(np.mean((diff - mean_gap[:, None]) ** 2, axis=1))
    cond_rms = np.sqrt(np.mean(diff ** 2, axis=1))
    cond_max = np.max(np.abs(diff), axis=1)
    return mean_gap, shape_rms, cond_rms, cond_max


# Reading thresholds for `summarize_curve_effects`. They label the numeric columns; they are not
# tests, and each one is documented so a reader can disagree with the cut explicitly.
EFFECT_RMS_FLOOR = 0.05          # lognorm units: below this there is nothing to describe
PATTERN_RMS_Z_TOL = 0.25         # z-units: the curves keep the same shape below this
AMPLITUDE_LOG2_TOL = 0.5         # a ~1.4x amplitude change
LEVEL_FRACTION_TOL = 0.5         # the difference is mostly a vertical offset
# Below this gradient amplitude a z-standardised comparison magnifies noise: standardising a nearly
# flat curve turns tiny fluctuations into an apparently large pattern difference. Such features get
# `pattern_status = 'insufficient amplitude'` instead of a pattern claim.
PATTERN_AMPLITUDE_FLOOR = 0.15


def _classify_difference(frame: pd.DataFrame) -> np.ndarray:
    """Dominant kind of difference: level, amplitude, redistribution, or nothing detectable.

    Reported separately from :func:`_pattern_status`, because a large level offset can coexist with an
    unsupported pattern comparison, and a single label that lets the pattern term win hides that.
    """
    detectable = frame['condition_effect_rms'].to_numpy() >= EFFECT_RMS_FLOOR
    level_dominant = frame['level_fraction'].to_numpy() >= LEVEL_FRACTION_TOL
    amplitude_moves = np.abs(frame['amplitude_log2_ratio'].to_numpy()) >= AMPLITUDE_LOG2_TOL
    pattern_supported = frame['pattern_status'].to_numpy() == 'supported'

    labels = np.full(len(frame), 'undetermined', dtype=object)
    labels[level_dominant] = 'level'
    labels[~level_dominant & amplitude_moves & ~pattern_supported] = 'amplitude'
    labels[~level_dominant & pattern_supported] = 'redistribution'
    labels[~detectable] = 'none'
    return labels


def _pattern_status(frame: pd.DataFrame) -> np.ndarray:
    """Whether a pattern comparison is meaningful for this feature.

    ``supported`` requires both curves to carry a gradient and the standardised difference to exceed
    the tolerance. Anything flatter is ``insufficient amplitude``, where standardising magnifies noise;
    a measurable gradient that does not move is ``unchanged``.
    """
    measurable = ((frame['amplitude_reference'].to_numpy() >= PATTERN_AMPLITUDE_FLOOR)
                  & (frame['amplitude_comparison'].to_numpy() >= PATTERN_AMPLITUDE_FLOOR))
    moved = frame['pattern_rms_z'].to_numpy() >= PATTERN_RMS_Z_TOL
    return np.where(~measurable, 'insufficient amplitude', np.where(moved, 'supported', 'unchanged'))


def summarize_curve_effects(curve_reference, curve_comparison, feature_names=None):
    """Split a fitted reference/comparison curve pair into level, amplitude and spatial pattern.

    ``run_level_shape``'s ``shape_rms`` is the RMS of the mean-centred difference, so a gene whose
    peak stays put while its gradient flattens is reported as a shape effect. Three questions are
    worth separating:

    * level - a constant vertical offset: ``level_effect`` is the mean difference over the grid.
    * amplitude - how strong the gradient is: peak-to-trough of each curve and their ratio.
    * pattern - where expression sits along pseudospace. ``pattern_rms_z`` compares the two curves
      after each is standardised over the grid, and ``curve_spearman`` ranks them; both are
      invariant to level and amplitude, so they only move when peak position, width or
      monotonicity changes.

    ``level_fraction`` / ``shape_fraction`` use the exact orthogonal split
    ``RMS(diff)**2 = mean(diff)**2 + RMS(diff - mean(diff))**2``: the share of the total fitted
    separation that is vertical offset rather than spatial redistribution.

    NaN cells (a group evaluated outside its own support) are ignored per feature. Returns one row
    per feature; ``difference_type`` is a reading aid and the numeric columns are the result.
    """
    lhs = np.asarray(curve_reference, dtype=float)
    rhs = np.asarray(curve_comparison, dtype=float)
    if lhs.shape != rhs.shape:
        raise ValueError(f'curve shapes differ: {lhs.shape} vs {rhs.shape}')
    if lhs.ndim != 2:
        raise ValueError('curves must be 2-D (feature x grid)')

    def _masked(values):
        return np.where(np.isfinite(values), values, np.nan)

    # A curve with no support on the grid (or zero spread) legitimately produces all-NaN
    # aggregations; those rows are reported as NaN rather than warned about.
    warnings.simplefilter('ignore', RuntimeWarning)
    diff = _masked(rhs - lhs)
    n_finite = np.sum(np.isfinite(diff), axis=1)
    level_effect = np.nanmean(diff, axis=1)
    centred = diff - level_effect[:, None]
    shape_rms = np.sqrt(np.nanmean(centred ** 2, axis=1))
    total_rms = np.sqrt(np.nanmean(diff ** 2, axis=1))
    total_sq = total_rms ** 2
    with np.errstate(invalid='ignore', divide='ignore'):
        level_fraction = np.where(total_sq > 0, level_effect ** 2 / total_sq, 0.0)
        shape_fraction = np.where(total_sq > 0, shape_rms ** 2 / total_sq, 0.0)

    def _amplitude(curve):
        finite = _masked(curve)
        return np.nanmax(finite, axis=1) - np.nanmin(finite, axis=1)

    amplitude_reference = _amplitude(lhs)
    amplitude_comparison = _amplitude(rhs)
    with np.errstate(invalid='ignore', divide='ignore'):
        amplitude_ratio = np.where(
            amplitude_reference > 0, amplitude_comparison / amplitude_reference, np.nan
        )
        amplitude_log2_ratio = np.where(
            (amplitude_reference > 0) & (amplitude_comparison > 0),
            np.log2(amplitude_comparison) - np.log2(amplitude_reference),
            np.nan,
        )

    def _standardized(curve):
        finite = _masked(curve)
        mean = np.nanmean(finite, axis=1, keepdims=True)
        spread = np.nanstd(finite, axis=1, keepdims=True)
        return np.where(spread > 0, (finite - mean) / spread, np.nan)

    z_diff = _standardized(rhs) - _standardized(lhs)
    pattern_rms_z = np.sqrt(np.nanmean(_masked(z_diff) ** 2, axis=1))

    frame = pd.DataFrame({
        'feature': (list(feature_names) if feature_names is not None
                    else [f'feature_{i}' for i in range(lhs.shape[0])]),
        'n_grid_points': n_finite,
        'level_effect': level_effect,
        'shape_rms': shape_rms,
        'condition_effect_rms': total_rms,
        'level_fraction': level_fraction,
        'shape_fraction': shape_fraction,
        'amplitude_reference': amplitude_reference,
        'amplitude_comparison': amplitude_comparison,
        'amplitude_ratio': amplitude_ratio,
        'amplitude_log2_ratio': amplitude_log2_ratio,
        'pattern_rms_z': pattern_rms_z,
    })
    # A pattern verdict is only offered when both curves carry a gradient, and the dominant kind of
    # difference is reported separately, so a large vertical offset cannot hide inside a pattern label.
    frame['pattern_status'] = _pattern_status(frame)
    frame['difference_type'] = _classify_difference(frame)
    frame['dominant_difference'] = frame['difference_type']
    return frame


def run_level_shape_summary(Y, s, c, knots, grid, lambda_grid, feature_names=None):
    """``run_level_shape`` plus the level/amplitude/pattern table for its fitted curves."""
    fit = run_level_shape(Y, s, c, knots, grid, lambda_grid)
    summary = summarize_curve_effects(fit['curve_healthy'], fit['curve_aki'],
                                      feature_names=feature_names)
    summary['curve_spearman'] = fit['curve_spearman']
    return fit, summary


def _xty_of(design, Y, sparse):
    return np.asarray(design.T @ Y) if sparse else (np.asarray(design, dtype=float).T @ np.asarray(Y, dtype=float))


def run_level_shape(Y, s, c, knots, grid, lambda_grid):
    """Nested M0/M1/M2 GAM with per-response GCV lambda (selected on M2, reused for M0/M1).
    Y: (n, g) sparse CSR (genes) or dense ndarray (pathway module scores)."""
    n = Y.shape[0]
    g = Y.shape[1]
    x0, x1, x2, sl0, sl1, sl2, p_b, k = build_ls_designs(s, c, knots)
    base_g, basis_g = ls_grid_designs(grid, knots)
    sparse = sp.issparse(Y)
    yty = (np.asarray(Y.multiply(Y).sum(axis=0)).ravel() if sparse
           else np.sum(np.asarray(Y, dtype=float) ** 2, axis=0))
    xty0, xty1, xty2 = _xty_of(x0, Y, sparse), _xty_of(x1, Y, sparse), _xty_of(x2, Y, sparse)
    sse0, edf0, _ = fit_sse_grid(x0.T @ x0, xty0, yty, n, sl0, lambda_grid)
    sse1, edf1, _ = fit_sse_grid(x1.T @ x1, xty1, yty, n, sl1, lambda_grid)
    sse2, edf2, betas2 = fit_sse_grid(x2.T @ x2, xty2, yty, n, sl2, lambda_grid, keep_beta=True)

    # GCV on the full M2 model -> one lambda per response, reused across M0/M1/M2.
    gcv2 = n * sse2 / np.maximum(n - edf2[:, None], 1e-6) ** 2
    lam_idx = np.argmin(gcv2, axis=0)
    ar = np.arange(g)
    s0, s1, s2 = sse0[lam_idx, ar], sse1[lam_idx, ar], sse2[lam_idx, ar]
    e0, e1, e2 = edf0[lam_idx], edf1[lam_idx], edf2[lam_idx]

    level_F = (np.maximum(s0 - s1, 0) / np.maximum(e1 - e0, 1e-6)) / np.maximum(s1 / np.maximum(n - e1, 1.0), 1e-12)
    shape_F = (np.maximum(s1 - s2, 0) / np.maximum(e2 - e1, 1e-6)) / np.maximum(s2 / np.maximum(n - e2, 1.0), 1e-12)

    beta2_sel = np.empty((x2.shape[1], g))
    for li in range(len(lambda_grid)):
        m = lam_idx == li
        if m.any():
            beta2_sel[:, m] = betas2[li][:, m]
    curve_h, curve_a = m2_curves(beta2_sel, base_g, basis_g, p_b)
    mean_gap, shape_rms, cond_rms, cond_max = shape_metrics(curve_h, curve_a)
    curve_sp = np.array([safe_spearman(curve_h[i], curve_a[i]) for i in range(g)])
    return dict(
        lam_idx=lam_idx, sl2=sl2, p_b=p_b, k=k,
        level_F=level_F, shape_F=shape_F,
        level_effect=mean_gap, m2_condition_coef=beta2_sel[p_b],
        shape_rms=shape_rms, condition_effect_rms=cond_rms, condition_effect_max_abs=cond_max,
        curve_spearman=curve_sp,
        amplitude_healthy=np.ptp(curve_h, axis=1), amplitude_aki=np.ptp(curve_a, axis=1),
        curve_healthy=curve_h, curve_aki=curve_a,
    )


def fit_single_condition_curves(Y, s, knots, grid, lambda_grid, lam_idx, support_pct=(1, 99)):
    """Fit intercept + f(s) (the M0 form) to ONE homogeneous group of cells, e.g. a single
    specimen, and evaluate on ``grid``.

    ``run_level_shape`` cannot be used per specimen: with one specimen the condition indicator is
    constant, so the M1/M2 designs are rank-deficient. This fits the shared-smooth model only.

    ``lam_idx`` is the per-feature index into ``lambda_grid`` chosen by GCV on the POOLED M2 fit;
    reusing it keeps every specimen curve on the same smoothing scale as the pooled condition
    curves it is drawn against, so visible spread is between-specimen variation and not a
    difference in how hard each curve was smoothed.

    Returns ``(curves, in_support)``: ``curves`` is (g, len(grid)) with NaN outside this group's
    own [p_lo, p_hi] pseudospace support -- a specimen is not extrapolated into a stretch of
    pseudospace where it has no cells. ``in_support`` is the boolean grid mask.
    """
    s = np.asarray(s, dtype=float).reshape(-1)
    grid = np.asarray(grid, dtype=float)
    sparse = sp.issparse(Y)
    g = Y.shape[1]

    base = make_gam_design(s, knots)
    p_b = base.shape[1]
    smooth_slices = [slice(1, p_b)]
    base_g = make_gam_design(grid, knots)

    xtx = base.T @ base
    xty = _xty_of(base, Y, sparse)

    curves = np.full((g, grid.size), np.nan)
    lam_idx = np.asarray(lam_idx, dtype=int).reshape(-1)
    for li in np.unique(lam_idx):                 # solve once per distinct lambda, not per feature
        cols = np.flatnonzero(lam_idx == li)
        penalty = make_gam_penalty(p_b, smooth_slices, len(s), alpha=float(lambda_grid[li]))
        beta = np.linalg.pinv(xtx + penalty) @ xty[:, cols]
        curves[cols] = (base_g @ beta).T

    p_lo, p_hi = support_pct
    in_support = (grid >= np.percentile(s, p_lo)) & (grid <= np.percentile(s, p_hi))
    curves[:, ~in_support] = np.nan
    return curves, in_support


def loso_shape_stability(Y, s, c, samples, knots, grid, lambda_grid, top_idx, feature_names=None):
    """Leave-one-specimen-out robustness of the top ``shape_rms`` hits.

    With 2 vs 2 specimens the honest question is not significance but whether a
    top shape hit survives dropping either specimen. For each specimen this refits
    the full nested level/shape model on the remaining tubules (still >=1 per
    condition) and records the recomputed ``shape_rms`` and its global rank for
    every feature in ``top_idx``.

    Y: (n, g) sparse CSR (genes) or dense (module scores) -- same object passed to
    :func:`run_level_shape`. Returns a long DataFrame
    (dropped_sample x feature) for a tornado/slope plot.
    """
    samples = np.asarray(samples).astype(str)
    uniq = sorted(np.unique(samples).tolist())
    top_idx = np.asarray(top_idx, dtype=int)

    base = run_level_shape(Y, s, c, knots, grid, lambda_grid)
    base_rms = base['shape_rms']
    base_rank = (-base_rms).argsort().argsort()          # 0 = largest shape_rms

    def _name(i):
        return feature_names[i] if feature_names is not None else int(i)

    rows = []
    for dropped in uniq:
        keep = samples != dropped
        n_cond = np.unique(np.asarray(c)[keep]).size
        if n_cond < 2:
            # dropping this specimen collapses a condition -- cannot fit level/shape
            for ti in top_idx:
                rows.append({'dropped_sample': dropped, 'feature': _name(ti),
                             'baseline_shape_rms': float(base_rms[ti]),
                             'loso_shape_rms': np.nan,
                             'baseline_rank': int(base_rank[ti]), 'loso_rank': -1})
            continue
        Ysub = Y[keep] if sp.issparse(Y) else np.asarray(Y)[keep]
        ls = run_level_shape(Ysub, np.asarray(s)[keep], np.asarray(c)[keep],
                             knots, grid, lambda_grid)
        rms = ls['shape_rms']
        rank = (-rms).argsort().argsort()
        for ti in top_idx:
            rows.append({'dropped_sample': dropped, 'feature': _name(ti),
                         'baseline_shape_rms': float(base_rms[ti]),
                         'loso_shape_rms': float(rms[ti]),
                         'baseline_rank': int(base_rank[ti]),
                         'loso_rank': int(rank[ti])})
    return pd.DataFrame(rows)


def sample_perm_pvalues(Y, s, samples, knots, grid, lambda_grid, lam_idx, sl2, p_b, control_samples):
    """Honest 2v2 significance: permute the healthy/aki label across the 4 SAMPLES.

    All C(4,2)=6 relabelings form the exhaustive (exact) reference set, and the observed labeling
    is one of them, so p is the inclusive tail ``#{splits with statistic >= observed} / 6``.
    Ties are counted as >=, which is what an exact test requires: if all six splits produced the
    same statistic, every split is "at least as extreme" and p = 1. The 6 splits form 3
    mirror-image (healthy/aki swap) pairs and shape_rms is symmetric under that swap -> a maximal
    feature lands on 2/6, not 1/6; 2/6 is the floor this design can reach. A magnitude-relative
    tolerance keeps the comparison deterministic against float noise on those near-tied twins.
    Each response's GCV-selected lambda* is held FIXED across all 6 relabelings (computational
    choice; mildly conservative). Pseudospace is NEVER shuffled."""
    n = Y.shape[0]
    g = Y.shape[1]
    ar = np.arange(g)
    base_g, basis_g = ls_grid_designs(grid, knots)
    base = make_gam_design(s, knots)
    basis = base[:, 1:]
    sparse = sp.issparse(Y)
    samples = np.asarray(samples).astype(str)
    uniq = sorted(np.unique(samples).tolist())
    splits = list(itertools.combinations(uniq, 2))            # 6 ways to pick 2 healthy samples
    perm_shape = np.empty((len(splits), g))
    perm_level = np.empty((len(splits), g))
    for si, healthy_set in enumerate(splits):
        cp = np.where(np.isin(samples, np.asarray(healthy_set)), 0.0, 1.0)
        x2p = np.column_stack([base, cp, cp[:, None] * basis])
        xtx2p = x2p.T @ x2p
        xty2p = _xty_of(x2p, Y, sparse)
        sh = np.empty((len(lambda_grid), g))
        lv = np.empty((len(lambda_grid), g))
        for li, lam in enumerate(lambda_grid):
            penalty = make_gam_penalty(x2p.shape[1], sl2, n, alpha=float(lam))
            beta = np.linalg.pinv(xtx2p + penalty) @ xty2p
            ch, ca = m2_curves(beta, base_g, basis_g, p_b)
            mg, srms, _, _ = shape_metrics(ch, ca)
            sh[li] = srms
            lv[li] = mg
        perm_shape[si] = sh[lam_idx, ar]
        perm_level[si] = lv[lam_idx, ar]
    true_idx = next(i for i, hs in enumerate(splits) if set(hs) == set(control_samples))
    obs_shape = perm_shape[true_idx]
    obs_level = np.abs(perm_level[true_idx])
    n_splits = len(splits)
    # Inclusive-tail exhaustive p-value: p = #{relabelings with statistic >= observed} / 6, the
    # observed labeling included among the six. Counting ties as "at least as extreme" is the
    # exact-test definition; the previous strict-greater count returned 1/6 for a statistic that
    # all six relabelings reproduced, where the correct value is 1. Because the six splits form
    # three mirror-image (healthy/aki swap) pairs and shape_rms is symmetric under the swap, a
    # maximal feature stops at 2/6 rather than 1/6: that is the floor of this design. A
    # magnitude-relative tolerance keeps the comparison deterministic against float noise on the
    # near-tied twins. Each response's GCV-selected lambda* is held FIXED across all 6 relabelings
    # (computational; mildly conservative). Pseudospace is NEVER shuffled.
    tol_shape = 1e-9 * np.maximum(np.abs(obs_shape), 1.0)
    tol_level = 1e-9 * np.maximum(obs_level, 1.0)
    shape_p = np.sum(perm_shape >= (obs_shape - tol_shape)[None, :], axis=0) / n_splits
    level_p = np.sum(np.abs(perm_level) >= (obs_level - tol_level)[None, :], axis=0) / n_splits
    return shape_p, level_p, splits, true_idx
