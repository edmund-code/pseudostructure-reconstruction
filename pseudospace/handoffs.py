"""Transporter handoffs along an axis (notebook 101; docs/results/pt-transporter-handoffs.md).

* ``scaled_profile``: a curve min-max scaled over its finite grid points (NaN elsewhere).
* ``crossing``: the first position where early - late changes sign from positive to non-positive, linearly
  interpolated between grid points; NaN when it never does.
* ``half_points``: the early profile's last fall below 0.5 and the late profile's first rise to 0.5 (interpolated).
* ``pair_count``: within-family early x late pairs, sum over families of (early members) x (late members).
* ``family_null``: that count over permutations of the family labels among the eligible genes.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def scaled_profile(curve):
    """Min-max scale over finite values; constant or empty curves give NaN."""
    c = np.asarray(curve, float)
    ok = np.isfinite(c)
    out = np.full(c.shape, np.nan)
    if ok.sum() < 2:
        return out
    lo, hi = c[ok].min(), c[ok].max()
    if hi <= lo:
        return out
    out[ok] = (c[ok] - lo) / (hi - lo)
    return out


def _interp_root(x0, x1, y0, y1):
    return x0 if y1 == y0 else x0 + (x1 - x0) * (-y0) / (y1 - y0)


def crossing(early, late, grid):
    """First x where early - late goes from > 0 to <= 0 (interpolated); NaN if never."""
    d = np.asarray(early, float) - np.asarray(late, float)
    g = np.asarray(grid, float)
    for i in range(len(d) - 1):
        if np.isfinite(d[i]) and np.isfinite(d[i + 1]) and d[i] > 0 >= d[i + 1]:
            return float(_interp_root(g[i], g[i + 1], d[i], d[i + 1]))
    return np.nan


def half_points(early, late, grid, *, level=.5):
    """(x where early last falls below ``level``, x where late first reaches it); NaN when absent."""
    e, l, g = (np.asarray(v, float) for v in (early, late, grid))
    fall = np.nan
    for i in range(len(e) - 1):
        if np.isfinite(e[i]) and np.isfinite(e[i + 1]) and e[i] >= level > e[i + 1]:
            fall = float(_interp_root(g[i], g[i + 1], e[i] - level, e[i + 1] - level))
    rise = np.nan
    for i in range(len(l) - 1):
        if np.isfinite(l[i]) and np.isfinite(l[i + 1]) and l[i] < level <= l[i + 1]:
            rise = float(_interp_root(g[i], g[i + 1], l[i] - level, l[i + 1] - level))
            break
    return fall, rise


def pair_count(family, early, late):
    """Sum over families of (early members) x (late members); arguments aligned per gene."""
    f = pd.Series(np.asarray(family))
    e, l = np.asarray(early, bool), np.asarray(late, bool)
    return int((f[e].value_counts() * f[l].value_counts()).fillna(0).sum())


def family_null(family, early, late, *, n=1000, seed=101):
    """``pair_count`` over permutations of the family labels among the genes given."""
    rng = np.random.default_rng(seed)
    f = np.asarray(family)
    return np.array([pair_count(rng.permutation(f), early, late) for _ in range(n)])
