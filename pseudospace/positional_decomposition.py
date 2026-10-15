"""A discovery screen along the PT: level, location and gradient survival (notebook 91).

Protocol: ``docs/results/pt-positional-decomposition.md``.

- ``level_location``, ``row_corr``: a condition effect Δ(x) split into its mean over the axis (level) and the rest
  (location), with the location's size.
- ``survival_parts``, ``cross_fitted_survival``, ``survival_class``: how much of a reference gradient survives in
  another condition, with the template and the denominator taken from different reference specimens.
- ``set_mean_z``, ``set_profile_z``, ``set_ratio_z``: set statistics against member-matched random sets.
- ``gradient_cosine``: agreement of two gradient vectors over a set (species).
- ``interior_dip``: a replicated interior minimum of a profile along the axis.
- ``pseudobulk_log_cpm``: per-specimen pseudobulk log2 CPM.
"""
from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.stats import norm

SURVIVAL_CUTS = (-.25, .25, .6, 1.5)
SURVIVAL_CLASSES = ('reversed', 'erased', 'attenuated', 'kept', 'amplified')


def level_location(delta):
    """Rows of ``delta`` (genes x bins; NaN outside the valid bins): level L (mean over valid bins), location
    P = Δ − L, and its size M (root mean square of P)."""
    D = np.asarray(delta, float)
    with np.errstate(invalid='ignore'):
        L = np.nanmean(D, axis=1)
        P = D - L[:, None]
        M = np.sqrt(np.nanmean(P ** 2, axis=1))
    return L, P, M


def row_corr(a, b):
    """Pearson r of each row of ``a`` with the same row of ``b`` over the bins finite in both (NaN when < 3 bins or
    a row is constant)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    n = ok.sum(axis=1)
    a0, b0 = np.where(ok, a, 0.), np.where(ok, b, 0.)
    with np.errstate(invalid='ignore', divide='ignore'):
        ma, mb = a0.sum(axis=1) / n, b0.sum(axis=1) / n
        da, db = np.where(ok, a - ma[:, None], 0.), np.where(ok, b - mb[:, None], 0.)
        r = (da * db).sum(axis=1) / np.sqrt((da ** 2).sum(axis=1) * (db ** 2).sum(axis=1))
    return np.where(n >= 3, r, np.nan)


def survival_parts(case, template, held):
    """Numerator Σ case·template and denominator Σ held·template per row, over bins finite in all three."""
    c, t, h = (np.asarray(x, float) for x in (case, template, held))
    ok = np.isfinite(c) & np.isfinite(t) & np.isfinite(h)
    num = np.where(ok, c * t, 0.).sum(axis=1)
    den = np.where(ok, h * t, 0.).sum(axis=1)
    return num, den


def cross_fitted_survival(cases, refs):
    """Gradient survival of the ``cases`` (list of deviation arrays) against two reference specimens ``refs``.

    For each case i and template t (the other reference h held out), A_{i,t} = Σ d_i·d_t / Σ d_h·d_t. A gradient
    taken from the same specimens as the denominator would share their noise and make the cases look flattened.
    Returns the per-gene median over the (i, t) combinations, and the numerator and denominator arrays
    (genes x combinations) for pooled and set ratios."""
    if len(refs) != 2:
        raise ValueError('Cross-fitted survival needs exactly two reference specimens.')
    nums, dens = [], []
    for t, h in ((0, 1), (1, 0)):
        for c in cases:
            n_, d_ = survival_parts(c, refs[t], refs[h])
            nums.append(n_)
            dens.append(d_)
    num, den = np.column_stack(nums), np.column_stack(dens)
    with np.errstate(invalid='ignore', divide='ignore'):
        per = np.where(den > 0, num / den, np.nan)
    med = np.full(len(per), np.nan)
    some = np.isfinite(per).any(axis=1)
    med[some] = np.nanmedian(per[some], axis=1)
    return med, num, den


def pooled_ratio(num, den, rows=None):
    """The median over combinations of Σ num / Σ den over ``rows`` (all rows when None)."""
    rows = slice(None) if rows is None else rows
    return float(np.median(np.asarray(num)[rows].sum(axis=0) / np.asarray(den)[rows].sum(axis=0)))


def survival_class(a, cuts=SURVIVAL_CUTS, classes=SURVIVAL_CLASSES):
    """Survival classes: reversed < −0.25 ≤ erased ≤ 0.25 < attenuated ≤ 0.6 < kept ≤ 1.5 < amplified."""
    a = np.asarray(a, float)
    lo, e, att, keep = cuts
    with np.errstate(invalid='ignore'):
        conds = [a < lo, a <= e, a <= att, a <= keep, a > keep]
    return np.select(conds, list(classes), default='').astype(object)


def _z_p(obs, null, tail):
    mu, sd = null.mean(), null.std(ddof=1)
    z = (obs - mu) / sd if sd > 0 else np.nan
    if tail == 'upper':
        p = norm.sf(z)
    elif tail == 'lower':
        p = norm.cdf(z)
    else:
        p = 2 * norm.sf(abs(z))
    return z, p


def set_mean_z(values, members, perms, *, tail='two'):
    """Mean of ``values`` over ``members`` (row indices) against ``perms`` (n_perm x genes index maps from
    ``strata_permutations``). Returns (observed, z, p)."""
    v = np.asarray(values, float)
    obs = float(v[members].mean())
    null = v[perms[:, members]].mean(axis=1)
    z, p = _z_p(obs, null, tail)
    return obs, z, p


def set_profile_z(profiles, members, perms):
    """Size (RMS over bins) of the members' mean profile against random sets; upper tail. ``profiles`` is
    genes x bins with NaN for invalid bins (the same bins for every gene). Returns (mean profile, size, z, p)."""
    P = np.asarray(profiles, float)
    ok = np.isfinite(P).all(axis=0)
    P = P[:, ok]
    prof = P[members].mean(axis=0)
    obs = float(np.sqrt(np.mean(prof ** 2)))
    null = np.sqrt((P[perms[:, members]].mean(axis=1) ** 2).mean(axis=1))
    z, p = _z_p(obs, null, 'upper')
    full = np.full(ok.shape, np.nan)
    full[ok] = prof
    return full, obs, z, p


def set_ratio_z(num, den, members, perms, *, tail='two'):
    """Pooled ratio (median over combinations of Σ num / Σ den) over ``members`` against random sets."""
    num, den = np.asarray(num, float), np.asarray(den, float)
    obs = pooled_ratio(num, den, members)
    null = np.median(num[perms[:, members]].sum(axis=1) / den[perms[:, members]].sum(axis=1), axis=1)
    z, p = _z_p(obs, null, tail)
    return obs, z, p


def gradient_cosine(a, b):
    """Σ a·b / √(Σ a² Σ b²) over finite pairs: +1 the same gradient direction, −1 reversed."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    den = np.sqrt((a[ok] ** 2).sum() * (b[ok] ** 2).sum())
    return float((a[ok] * b[ok]).sum() / den) if den > 0 else np.nan


def set_cosine_z(a, b, members, perms, *, tail='two'):
    """``gradient_cosine`` over ``members`` against random sets (a and b finite for every row)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    obs = gradient_cosine(a[members], b[members])
    pa, pb = a[perms[:, members]], b[perms[:, members]]
    null = (pa * pb).sum(axis=1) / np.sqrt((pa ** 2).sum(axis=1) * (pb ** 2).sum(axis=1))
    z, p = _z_p(obs, null, tail)
    return obs, z, p


def set_log_rms_ratio_z(a, b, members, perms, *, tail='two'):
    """log2(RMS b / RMS a) over ``members`` against random sets."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    obs = float(.5 * np.log2((b[members] ** 2).mean() / (a[members] ** 2).mean()))
    null = .5 * np.log2((b[perms[:, members]] ** 2).mean(axis=1) / (a[perms[:, members]] ** 2).mean(axis=1))
    z, p = _z_p(obs, null, tail)
    return obs, z, p


def interior_dip(profile, centres, *, lo=.35, hi=1.95, frac=.1):
    """Is the minimum of ``profile`` over centres in [lo, hi] an interior dip? Both side maxima must exceed it by
    ≥ ``frac`` × the profile's range. Returns (dip, position of the minimum, smaller side excess / range)."""
    y, c = np.asarray(profile, float), np.asarray(centres, float)
    ok = np.isfinite(y)
    inner = np.flatnonzero(ok & (c >= lo - 1e-9) & (c <= hi + 1e-9))
    if len(inner) == 0:
        return False, np.nan, np.nan
    j = inner[np.argmin(y[inner])]
    left, right = y[ok & (c < c[j])], y[ok & (c > c[j])]
    rng = np.ptp(y[ok])
    if len(left) == 0 or len(right) == 0 or rng <= 0:
        return False, float(c[j]), np.nan
    excess = min(left.max(), right.max()) - y[j]
    return bool(excess >= frac * rng), float(c[j]), float(excess / rng)


def pseudobulk_log_cpm(counts, library, specimen, names, *, pseudo=1.):
    """Specimens x genes: log2((Σ counts + pseudo) / Σ library × 1e6) for each specimen in ``names``."""
    X = sparse.csr_matrix(counts)
    specimen = np.asarray(specimen)
    out = []
    for n in names:
        rows = specimen == n
        lib = float(np.asarray(library, float)[rows].sum())
        with np.errstate(divide='ignore'):
            out.append(np.log2((np.asarray(X[rows].sum(axis=0)).ravel() + pseudo) / lib * 1e6))
    return np.vstack(out)
