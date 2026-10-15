"""Glucose handling along the PT, and the size of a drug experiment that could test it (notebook 92).

Protocol: ``docs/results/pt-glucose-handling.md``.

- ``bootstrap_deviation``: deviation profiles of target genes on units resampled within one specimen.
- ``switch_fit``: notebook 90's switch model alone (intercept plus a logistic), for many profiles at once.
- ``balanced_relabelings``: how many other two-group splits a design with n specimens per group has, and the smallest
  permutation p it can reach.
- ``within_group_sd``: the between-specimen SD of a per-specimen statistic from replicates of each condition, pooled
  within expression strata.
- ``minimum_detectable_effect``: the smallest group difference a two-sample t-test detects at n per group.
- ``luminal_glucose``, ``calibrate_capacities``: a steady-state glucose balance along the PT with transporter
  densities that follow measured profiles (an illustration with stated assumptions, not an inference).
"""
from __future__ import annotations

from math import comb

import numpy as np
from scipy import sparse, stats
from scipy.optimize import brentq

from .positional_programs import MIN_UNITS, bin_profiles, reference_profile


def bootstrap_deviation(counts, library, bins, n_bins, targets, usable, *, n_boot, seed, min_units=MIN_UNITS):
    """Resample one specimen's units with replacement ``n_boot`` times. Each draw rebins the units, centres every gene
    on its mean over valid bins, and subtracts the reference profile (the median ``usable`` gene per bin), as the
    unresampled profiles do. Returns n_boot x len(targets) x n_bins (NaN where a bin has < ``min_units`` units).

    Units of one tubule can neighbour each other, so this understates the uncertainty; the specimens are the
    replicates."""
    X = sparse.csr_matrix(counts)
    lib = np.asarray(library, float)
    b = np.asarray(bins)
    keep = np.flatnonzero(b >= 0)
    rng = np.random.default_rng(seed)
    out = np.full((n_boot, len(targets), n_bins), np.nan)
    for k in range(n_boot):
        idx = keep[rng.integers(0, len(keep), len(keep))]
        prof = bin_profiles(X[idx], lib[idx], b[idx], n_bins, np.zeros(len(idx), int), [0], min_units=min_units)[0]
        out[k] = prof[targets] - reference_profile(prof, usable)
    return out


def switch_fit(profiles, weights, centres, *, m_grid=None, w_grid=None):
    """90's switch model S (``continuous_axis.fit_ladder``, one specimen) on every row of ``profiles``: weighted least
    squares of the profile on an intercept and a logistic with midpoint m and 10–90% width w, over 90's grids. Bins
    that are NaN in any row are dropped (in a bootstrap draw the thin bins are the same for every gene). Returns
    (midpoint, width, amplitude) arrays."""
    Y = np.atleast_2d(np.asarray(profiles, float))
    u = np.asarray(centres, float)
    w = np.asarray(weights, float)
    ok = np.isfinite(Y).all(axis=0) & (w > 0)
    Y, u, w = Y[:, ok], u[ok], w[ok] / w[ok].mean()
    m_grid = np.linspace(.15, 2.35, 45) if m_grid is None else np.asarray(m_grid, float)
    w_grid = np.geomspace(.1, 4., 25) if w_grid is None else np.asarray(w_grid, float)
    mm, ww = (g.ravel() for g in np.meshgrid(m_grid, w_grid, indexing='ij'))
    X = 1. / (1. + np.exp(-np.clip((u[None, :] - mm[:, None]) * np.log(81.) / ww[:, None], -60, 60)))
    keep = np.ptp(X, axis=1) >= 1e-6
    X, mm, ww = X[keep], mm[keep], ww[keep]
    Xc = X - (X * w).sum(axis=1, keepdims=True) / w.sum()
    Yc = Y - (Y * w).sum(axis=1, keepdims=True) / w.sum()
    num = (Xc * w) @ Yc.T                                   # combos x rows
    den = (Xc ** 2 * w).sum(axis=1)[:, None]
    rss = (Yc ** 2 * w).sum(axis=1)[None, :] - num ** 2 / den
    best = np.argmin(rss, axis=0)
    cols = np.arange(Y.shape[0])
    return mm[best], ww[best], num[best, cols] / den[best, 0]


def balanced_relabelings(n):
    """(other splits, smallest p): splits of 2n specimens into two groups of n other than the true one, and
    1 / all splits, the smallest permutation p the design can reach."""
    if n < 1:
        raise ValueError('Need at least one specimen per group.')
    total = comb(2 * n, n) // 2
    return total - 1, 1. / total


def within_group_sd(values, groups, strata):
    """``values``: genes x specimens of a per-specimen statistic; ``groups``: lists of column indices, one per
    condition (replicates); ``strata``: an expression stratum per gene. Each gene's pooled within-group variance has
    df = Σ (n_g − 1). Its median within a stratum, divided by the median of χ²_df / df, estimates the stratum's
    variance (the median resists genes with a real specimen difference). Returns (per-gene variance, df, {stratum: SD})."""
    V = np.asarray(values, float)
    ss = np.zeros(len(V))
    df = 0
    for g in groups:
        x = V[:, list(g)]
        ss += ((x - x.mean(axis=1, keepdims=True)) ** 2).sum(axis=1)
        df += len(g) - 1
    if df < 1:
        raise ValueError('Need a condition with at least two specimens.')
    var = ss / df
    scale = stats.chi2.median(df) / df
    strata = np.asarray(strata)
    sd = {s: float(np.sqrt(np.nanmedian(var[strata == s]) / scale)) for s in np.unique(strata)}
    return var, df, sd


def minimum_detectable_effect(sd, n, *, alpha=.05, power=.8):
    """(t_{1−α/2} + t_{power}) · sd · √(2/n) with df = 2n − 2: the true difference a two-sided two-sample t-test
    detects with the given power at n specimens per group."""
    if n < 2:
        raise ValueError('A t-test needs at least two specimens per group.')
    df = 2 * n - 2
    return float((stats.t.ppf(1 - alpha / 2, df) + stats.t.ppf(power, df)) * sd * np.sqrt(2. / n))


def luminal_glucose(s2, s1, a2, a1, *, c0, k2=4.9, k1=1.8, flow_end=1 / 3, substeps=40):
    """Glucose along the PT at steady state, one row per bin of equal length.

    The luminal glucose flux F (a fraction of the filtered load) falls as dF/dx = −(a2·s2·c/(k2 + c) + a1·s1·c/(k1 + c))
    / c0, with c = c0·F/Q the luminal concentration (mM) and Q the tubular flow, falling linearly from 1 to
    ``flow_end`` along the PT. ``s2`` and ``s1`` are relative SGLT2 and SGLT1 densities per bin (constant within a
    bin); a2 and a1 scale them to capacities (filtered loads per unit length at saturation, times c0). k2 and k1 are
    the human K0.5 values for glucose (Hummel 2011, 37 °C). Returns (F at the bin edges, uptake per bin by SGLT2,
    uptake per bin by SGLT1), all as fractions of the filtered load.

    Ceiling: mRNA stands for transporter density, every bin is the same length, flow falls linearly, and paracellular
    flux, GLUTs and the Na+ gradient are ignored. Upgrade path: a published nephron model (Layton's) with these
    profiles as its axial densities."""
    s2, s1 = np.asarray(s2, float), np.asarray(s1, float)
    n = len(s2)
    if len(s1) != n or not (np.all(s2 >= 0) and np.all(s1 >= 0)):
        raise ValueError('s2 and s1 must be non-negative and the same length.')
    h = 1. / substeps
    F = np.empty(n + 1)
    F[0] = 1.
    up2, up1 = np.zeros(n), np.zeros(n)

    def rates(f, x, j):
        q = 1. - (1. - flow_end) * x / n
        c = c0 * max(f, 0.) / q
        return a2 * s2[j] * c / (k2 + c) / c0, a1 * s1[j] * c / (k1 + c) / c0

    f = 1.
    for j in range(n):
        for k in range(substeps):          # RK4 in position x (bins)
            x = j + k * h
            r = [sum(rates(f, x, j))]
            r.append(sum(rates(f - h / 2 * r[0], x + h / 2, j)))
            r.append(sum(rates(f - h / 2 * r[1], x + h / 2, j)))
            r.append(sum(rates(f - h * r[2], x + h, j)))
            step = h * (r[0] + 2 * r[1] + 2 * r[2] + r[3]) / 6
            r2, r1 = rates(f, x + h / 2, j)
            share = r2 / (r2 + r1) if r2 + r1 > 0 else 0.
            up2[j] += step * share
            up1[j] += step * (1 - share)
            f = max(f - step, 0.)
        F[j + 1] = f
    return F, up2, up1


def calibrate_capacities(s2, s1, *, block_reabsorbed, threshold, threshold_excreted=.01, c_block=8., **kw):
    """(a2, a1): a1 so that with SGLT2 absent the PT reabsorbs ``block_reabsorbed`` of the load at plasma glucose
    ``c_block`` mM; then a2 so that ``threshold_excreted`` of the load escapes at plasma glucose ``threshold`` mM."""
    def reabsorbed(a2, a1, c0):
        return 1. - luminal_glucose(s2, s1, a2, a1, c0=c0, **kw)[0][-1]

    a1 = brentq(lambda a: reabsorbed(0., a, c_block) - block_reabsorbed, 1e-6, 1e4)
    a2 = brentq(lambda a: (1. - reabsorbed(a, a1, threshold)) - threshold_excreted, 1e-6, 1e5)
    return a2, a1
