"""Where pathways that a segment-level analysis puts in one PT segment sit inside it (notebook 84).

Protocol: ``docs/results/pt-within-segment-pathways.md``.

- ``microdissection_segment_t``: D1, the moderated t of one segment against the mean of the others in replicated
  microdissection samples.
- ``greedy_representatives``: D4, one representative per group of pathways whose cores overlap.
- ``thirds_profile``: C2, a curve's first, middle and last thirds inside one segment, relative to each other.
- ``replicated_spread`` and ``position_test``: C3, the between-pathway spread of those profiles counted where two
  specimens agree, against a stratified gene permutation null.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .program_localization import moderated_t

SEGMENTS = ('S1', 'S2', 'S3')
FLAT = .25              # log2; a curve whose range inside the segment is below this is flat there (76's O5 threshold)


def microdissection_segment_t(logexpr, groups, segment, *, segments=SEGMENTS, covariate=None):
    """Moderated t of ``segment`` against the mean of the other segments (D1).

    ``logexpr`` is genes x samples (log2), ``groups`` each sample's segment. The residual variance is the one-way
    model's over ``segments``; the effect is divided by sqrt(1/n_k + sum(1/n_o) / len(others)**2) before
    ``moderated_t`` (limma squeezeVar, trend on ``covariate``). Returns moderated_t's table plus ``log2_effect``.
    """
    x = pd.DataFrame(logexpr, dtype=float)
    groups = np.asarray(groups).astype(str)
    if segment not in segments or set(groups) - set(segments):
        raise ValueError('Every sample needs a segment among ``segments``, and ``segment`` must be one of them.')
    means, n, resid = {}, {}, np.zeros(len(x))
    for s in segments:
        block = x.loc[:, groups == s]
        if block.shape[1] < 1:
            raise ValueError(f'No sample of {s}.')
        means[s], n[s] = block.mean(axis=1), block.shape[1]
        resid += ((block.sub(means[s], axis=0)) ** 2).sum(axis=1).to_numpy()
    df = len(groups) - len(segments)
    if df < 1:
        raise ValueError('Need more samples than segments.')
    others = [s for s in segments if s != segment]
    effect = means[segment] - sum(means[s] for s in others) / len(others)
    scale = np.sqrt(1 / n[segment] + sum(1 / n[s] for s in others) / len(others) ** 2)
    out = moderated_t(effect / scale, resid / df, df, covariate=covariate)
    out['log2_effect'] = effect.to_numpy()
    return out


def greedy_representatives(cores, order, *, overlap=.6):
    """Walk ``order``; keep a pathway unless its core overlaps a kept one's by >= ``overlap`` (|A ∩ B| / min).

    Returns {representative: [absorbed pathways]} in the order kept (D4)."""
    reps = {}
    for p in order:
        a = set(cores[p])
        for r in reps:
            b = set(cores[r])
            if a and b and len(a & b) / min(len(a), len(b)) >= overlap:
                reps[r].append(p)
                break
        else:
            reps[p] = []
    return reps


def thirds_profile(curve, units, k, *, n_parts=3, flat=FLAT):
    """The curve inside segment ``k`` (k <= units < k + 1): part means relative to their mean, range, peak (C2).

    Returns profile (n_parts values summing to 0; NaN when a part has no grid point), range, peak position, peak
    part (0-based argmax of the profile) and flat (range < ``flat``)."""
    y, u = np.asarray(curve, float), np.asarray(units, float)
    inside = (u >= k) & (u < k + 1) & np.isfinite(y)
    empty = {'profile': np.full(n_parts, np.nan), 'range': np.nan, 'peak position': np.nan, 'peak part': -1,
             'flat': True}
    if inside.sum() < 2:
        return empty
    uu, yy = u[inside], y[inside]
    edges = k + np.arange(n_parts + 1) / n_parts
    parts = np.array([yy[(uu >= edges[i]) & (uu < edges[i + 1])].mean() if ((uu >= edges[i]) & (uu < edges[i + 1])).any()
                      else np.nan for i in range(n_parts)])
    if np.isnan(parts).all():
        return empty
    profile = parts - np.nanmean(parts)
    r = float(yy.max() - yy.min())
    return {'profile': profile, 'range': r, 'peak position': float(uu[np.argmax(yy)]),
            'peak part': int(np.nanargmax(profile)), 'flat': bool(r < flat)}


def replicated_spread(profiles_a, profiles_b):
    """Σ over rows of (a_p − mean a) · (b_p − mean b): the between-row spread counted only where A and B agree (C3).

    Rows with any NaN in either array are dropped first."""
    a, b = np.asarray(profiles_a, float), np.asarray(profiles_b, float)
    ok = np.isfinite(a).all(axis=1) & np.isfinite(b).all(axis=1)
    a, b = a[ok], b[ok]
    if len(a) < 2:
        return np.nan
    return float(((a - a.mean(axis=0)) * (b - b.mean(axis=0))).sum())


def stratified_permutation(strata, rng):
    """One permutation of ``range(len(strata))`` that moves indices only within their stratum."""
    strata = np.asarray(strata)
    perm = np.arange(len(strata))
    for s in np.unique(strata):
        idx = np.flatnonzero(strata == s)
        perm[idx] = rng.permutation(idx)
    return perm


def position_test(curves_a, curves_b, units, cores, k, strata, *, n_perm=2000, seed=84, n_parts=3):
    """T and its permutation p for pathway cores inside segment ``k`` (C3).

    ``curves_a``/``curves_b``: pool genes x grid (log2) for the two specimens, rows in pool order; ``cores``: one
    index array into the pool per pathway; ``strata``: one label per pool gene. A pathway curve is the nan-median over
    its core rows. The same permutation is used for both specimens. Returns T, p, the null sample and the profiles."""
    a, b = np.asarray(curves_a, float), np.asarray(curves_b, float)
    cores = [np.asarray(c, int) for c in cores]

    def profiles(perm):
        pa = np.array([thirds_profile(np.nanmedian(a[perm[c]], axis=0), units, k, n_parts=n_parts)['profile'] for c in cores])
        pb = np.array([thirds_profile(np.nanmedian(b[perm[c]], axis=0), units, k, n_parts=n_parts)['profile'] for c in cores])
        return pa, pb

    identity = np.arange(len(a))
    pa, pb = profiles(identity)
    t = replicated_spread(pa, pb)
    rng = np.random.default_rng(seed)
    null = np.empty(n_perm)
    for i in range(n_perm):
        null[i] = replicated_spread(*profiles(stratified_permutation(strata, rng)))
    p = (1 + np.sum(null >= t)) / (1 + n_perm)
    return {'T': t, 'p': float(p), 'null': null, 'profiles': (pa, pb)}
