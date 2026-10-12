"""Published human vs mouse findings re-read along the PT (notebooks 81 and 82).

Protocol: ``docs/results/pt-published-findings.md``. The helpers here are the rules that protocol fixes:
Cheval 2012 SAGE counts and S3 − S1 contrasts, agreement by gradient strength, the class rules for a program's
contrast in independent data (one library, B2c; or donors, B3), the gene class of a CellSpectra block member (A2), a
Mantel–Haenszel test (E1), the continuous descriptors of a program curve with the species-level agreement (G), and
the pathway label (L). Notebook 79's statistics come from ``program_localization``.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import binomtest, norm, rankdata, spearmanr, wilcoxon

from .continuous_descriptors import interior_extrema, monotone_onsets

EFFECT = .25            # log2; 76's O5 threshold, used throughout
SEGMENTS = ('S1', 'S2', 'S3')


# --- Cheval 2012 SAGE (B1, B2) ----------------------------------------------------------------------------------

def sage_counts(values):
    """Tag counts and library sizes from SAGE values per 10,000 tags.

    One tag is the column's smallest non-zero value, so a column's library is 10,000 / that value. Returns
    (counts frame, library series). Missing values stay missing.
    """
    values = pd.DataFrame(values, dtype=float)
    unit = values.where(values > 0).min()
    if unit.isna().any():
        raise ValueError(f'columns without a non-zero value: {list(unit.index[unit.isna()])}')
    return (values / unit).round(), 1e4 / unit


def sage_contrast(c_hi, c_lo, n_hi, n_lo, pseudo=.5):
    """log2 rate of ``hi`` over ``lo`` (pseudocount ``pseudo``) and its Poisson SE."""
    c_hi, c_lo = np.asarray(c_hi, float), np.asarray(c_lo, float)
    d = np.log2((c_hi + pseudo) / n_hi) - np.log2((c_lo + pseudo) / n_lo)
    se = np.sqrt(1 / (c_hi + pseudo) + 1 / (c_lo + pseudo)) / np.log(2)
    return d, se


def zonated(d, se, *, min_effect=.5, min_z=2.):
    d, se = np.asarray(d, float), np.asarray(se, float)
    return (np.abs(d) >= min_effect) & (np.abs(d) >= min_z * se)


def cochran_armitage(successes, totals, scores=None):
    """One-sided (increasing) Cochran–Armitage trend test of proportions. Returns (z, p)."""
    x, n = np.asarray(successes, float), np.asarray(totals, float)
    s = np.arange(len(n), dtype=float) if scores is None else np.asarray(scores, float)
    keep = n > 0
    x, n, s = x[keep], n[keep], s[keep]
    if len(n) < 2:
        return np.nan, np.nan
    p = x.sum() / n.sum()
    t = (s * (x - n * p)).sum()
    var = p * (1 - p) * ((n * s * s).sum() - (n * s).sum() ** 2 / n.sum())
    if var <= 0:
        return np.nan, np.nan
    z = t / np.sqrt(var)
    return float(z), float(norm.sf(z))


def strength_agreement(d_select, d_other, selected, bins=(.5, 1., 2., np.inf)):
    """Share of selected genes whose other-species contrast has the same sign, by |d_select| bin (B1a).

    Returns (per-bin frame with n, agree, share and binomial p against 0.5 one-sided greater; trend z and p).
    """
    a, b = np.asarray(d_select, float), np.asarray(d_other, float)
    sel = np.asarray(selected, bool)
    rows = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        m = sel & (np.abs(a) >= lo) & (np.abs(a) < hi)
        n = int(m.sum())
        k = int((np.sign(a[m]) == np.sign(b[m])).sum())
        rows.append({'bin': f'[{lo:g}, {hi:g})', 'n': n, 'agree': k, 'share': k / n if n else np.nan,
                     'p (> 0.5)': binomtest(k, n, .5, alternative='greater').pvalue if n else np.nan})
    frame = pd.DataFrame(rows)
    return frame, cochran_armitage(frame.agree, frame.n)


def slope_through_origin(x, y, *, n_boot=2000, seed=0, level=.95):
    """Least-squares slope of y on x through the origin, with a percentile bootstrap over genes."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    est = float((x * y).sum() / (x * x).sum())
    rng = np.random.default_rng(seed)
    draws = np.empty(n_boot)
    for b in range(n_boot):
        i = rng.integers(0, len(x), len(x))
        draws[b] = (x[i] * y[i]).sum() / (x[i] * x[i]).sum()
    tail = (1 - level) / 2
    return est, float(np.quantile(draws, tail)), float(np.quantile(draws, 1 - tail))


def _rank_residual(v, covariates):
    r = rankdata(v)
    if covariates is None:
        return r - r.mean()
    z = np.column_stack([np.ones(len(r))] + [rankdata(c) for c in np.asarray(covariates, float).T])
    return r - z @ np.linalg.lstsq(z, r, rcond=None)[0]


def permutation_spearman(x, y, *, covariates=None, n_perm=10000, seed=0):
    """Spearman (or rank-partial, given covariates genes x k) of x and y; one-sided permutation p (greater).

    Both rank residuals are computed once; the permutation shuffles y's residual against x's (residual
    permutation), so the null keeps each variable's dependence on the covariates.
    """
    rx, ry = _rank_residual(np.asarray(x, float), covariates), _rank_residual(np.asarray(y, float), covariates)
    norm_ = np.sqrt((rx @ rx) * (ry @ ry))
    rho = float(rx @ ry / norm_)
    rng = np.random.default_rng(seed)
    null = np.array([rx @ ry[rng.permutation(len(ry))] for _ in range(n_perm)]) / norm_
    return rho, float((1 + (null >= rho).sum()) / (1 + n_perm))


# --- Class rules in independent data (B2c, B3) ------------------------------------------------------------------

def library_class(contrasts, *, effect=EFFECT, share_hi=.7, share_lo=.3, alpha=.05):
    """B2c: class of a program's own-segment member contrasts in one library per species (Cheval)."""
    c = np.asarray(contrasts, float)
    c = c[np.isfinite(c)]
    if len(c) < 5:
        return {'class': 'not testable', 'n': len(c), 'median': np.nan, 'share': np.nan, 'p': np.nan}
    m, pos = float(np.median(c)), int((c > 0).sum())
    nonzero = int((c != 0).sum())
    share = pos / nonzero if nonzero else np.nan
    p_up = binomtest(pos, nonzero, .5, alternative='greater').pvalue if nonzero else 1.
    p_dn = binomtest(pos, nonzero, .5, alternative='less').pvalue if nonzero else 1.
    if m >= effect and share >= share_hi and p_up <= alpha:
        cls = 'zonated'
    elif m <= -effect and share <= share_lo and p_dn <= alpha:
        cls = 'reversed'
    elif abs(m) < effect:
        cls = 'absent'
    else:
        cls = 'attenuated'
    return {'class': cls, 'n': len(c), 'median': m, 'share': share, 'p': float(min(p_up, p_dn))}


def donor_class(values, *, effect=EFFECT, fraction=.8, alpha=.05):
    """B3: class of a program's per-donor median contrasts."""
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    if len(v) < 3:
        return {'class': 'not testable', 'donors': len(v), 'median': np.nan, 'share': np.nan, 'p': np.nan}
    m = float(np.median(v))
    up, dn = float((v >= effect).mean()), float((v <= -effect).mean())
    p_up = float(wilcoxon(v, alternative='greater').pvalue) if np.any(v != 0) else 1.
    p_dn = float(wilcoxon(v, alternative='less').pvalue) if np.any(v != 0) else 1.
    if m >= effect and up >= fraction and p_up <= alpha:
        cls, share, p = 'zonated', up, p_up
    elif m <= -effect and dn >= fraction and p_dn <= alpha:
        cls, share, p = 'reversed', dn, p_dn
    elif abs(m) < effect:
        cls, share, p = 'absent', up, p_up
    else:
        cls, share, p = 'attenuated', up, p_up
    return {'class': cls, 'donors': len(v), 'median': m, 'share': share, 'p': p}


def segment_contrasts(rates, segments=SEGMENTS):
    """Each segment minus the mean of the others, from a frame whose columns are segments (genes x segments)."""
    rates = pd.DataFrame(rates)
    return pd.DataFrame({s: rates[s] - rates[[t for t in segments if t != s]].mean(axis=1) for s in segments})


def donor_segment_rates(sums, totals, meta, genes, *, min_cells=50, segments=SEGMENTS, cpm_floor=1.):
    """Per donor, log2 CPM (+ floor) per segment from pseudobulk sums (rows) with meta columns donor, segment, cells.

    Only donors with at least ``min_cells`` in every listed segment are kept. Returns {donor: genes x segments}.
    """
    meta = pd.DataFrame(meta).reset_index(drop=True)
    frame = pd.DataFrame(np.asarray(sums, float), columns=list(genes))
    totals = np.asarray(totals, float)
    cells = meta.pivot_table(index='donor', columns='segment', values='cells', aggfunc='sum', fill_value=0)
    donors = cells.index[(cells.reindex(columns=list(segments), fill_value=0) >= min_cells).all(axis=1)]
    out = {}
    for donor in donors:
        cols = {}
        for seg in segments:
            rows = meta.index[meta.donor.eq(donor) & meta.segment.eq(seg)]
            cols[seg] = np.log2(frame.loc[rows].sum() / totals[rows].sum() * 1e6 + cpm_floor)
        out[donor] = pd.DataFrame(cols)
    return out


# --- A2: gene class of a published block member -----------------------------------------------------------------

def gene_class(q_condition, q_pairs, listed, *, effect=EFFECT):
    """Class of one gene (A2). ``q_condition`` is length 3 (S1..S3, signed by the claim), ``q_pairs`` pairs x 3,
    ``listed`` the listed segment's index. Rules in the protocol's order."""
    qc, qp = np.asarray(q_condition, float), np.asarray(q_pairs, float)
    if not np.isfinite(qc).all() or not np.isfinite(qp).all():
        return 'not measured'
    all_pairs = lambda j, op: op(qc[j]) and all(op(v) for v in qp[:, j])    # noqa: E731
    above = lambda v: v >= effect                                            # noqa: E731
    below = lambda v: v <= -effect                                           # noqa: E731
    if all(all_pairs(j, above) for j in range(3)):
        return 'present in all segments'
    others = [j for j in range(3) if j != listed]
    if all_pairs(listed, above) and any(qc[j] < effect for j in others):
        return 'restricted to the listed segment'
    if qc[listed] < effect and any(all_pairs(j, above) for j in others):
        return 'elsewhere only'
    if any(all_pairs(j, below) for j in range(3)):
        return 'opposite somewhere'
    return 'not reproduced'


# --- E1: Mantel–Haenszel ----------------------------------------------------------------------------------------

def mantel_haenszel(exposed, outcome, strata):
    """Common odds ratio (Mantel–Haenszel) of outcome given exposure across strata, with the CMH test (no continuity
    correction); p is one-sided for OR > 1. Returns dict(or_mh, z, p, n, a)."""
    e, o, s = np.asarray(exposed, bool), np.asarray(outcome, bool), np.asarray(strata)
    num = den = t = var = 0.
    a_tot = 0
    for level in np.unique(s):
        m = s == level
        a, b = np.sum(e[m] & o[m]), np.sum(e[m] & ~o[m])
        c, d = np.sum(~e[m] & o[m]), np.sum(~e[m] & ~o[m])
        n = a + b + c + d
        if n < 2:
            continue
        num += a * d / n
        den += b * c / n
        t += a - (a + b) * (a + c) / n
        var += (a + b) * (c + d) * (a + c) * (b + d) / (n * n * (n - 1))
        a_tot += a
    z = t / np.sqrt(var) if var > 0 else np.nan
    return {'or_mh': num / den if den > 0 else np.inf, 'z': float(z), 'p': float(norm.sf(z)) if np.isfinite(z) else np.nan,
            'n': int(len(e)), 'a': int(a_tot)}


# --- G: continuous descriptors ----------------------------------------------------------------------------------

def reference_curve(curves, usable):
    """The typical gene's shape: the median over usable genes of their mean-centred curves (natural log units).

    Subtracting it is the within-species analogue of notebook 76's ``centre_rows`` (the median gene at each position):
    it removes the shape every gene shares in that specimen, which random gene sets otherwise reproduce.
    """
    y = np.asarray(curves, float)[np.asarray(usable, bool)]
    return np.median(y - y.mean(axis=1, keepdims=True), axis=0)


def program_curve(curves, usable, *, reference=None, log=np.log(2), min_members=5):
    """Median over usable members of their mean-centred curves minus ``reference`` (if given), in log2.

    ``curves`` is members x grid in natural log; ``reference`` a grid-length curve in the same units."""
    y = np.asarray(curves, float)[np.asarray(usable, bool)]
    if len(y) < min_members:
        return np.full(np.asarray(curves).shape[1], np.nan), len(y)
    y = y - y.mean(axis=1, keepdims=True)
    if reference is not None:
        y = y - np.asarray(reference, float)
    return np.median(y, axis=0) / log, len(y)


def segment_profile(rates, members, *, segments=SEGMENTS, min_members=1):
    """The cluster-resolution analogue of ``program_curve`` for one specimen or donor (log2 in and out).

    ``rates`` is genes x segments. Each gene's segment levels minus its mean over the segments, centred by the median
    gene per segment; the median over the ``members`` present. Genes missing a segment are dropped first."""
    r = pd.DataFrame(rates)[list(segments)].dropna()
    c = r.sub(r.mean(axis=1), axis=0)
    c = c - c.median()
    present = [g for g in dict.fromkeys(members) if g in c.index]
    if len(present) < min_members:
        return pd.Series(np.nan, index=list(segments))
    return c.loc[present].median()


def within_segment_change(curve, units, *, segments=SEGMENTS, n_end=3, trim=5):
    """Per segment, the curve's change from its first to its last ``n_end`` grid points inside that segment, and the
    position (segment units) and slope of the steepest point, ignoring ``trim`` grid points at each end (spline
    edges). What a segment mean cannot show."""
    y, u = np.asarray(curve, float), np.asarray(units, float)
    out = {}
    for k, seg in enumerate(segments):
        inside = np.flatnonzero((u >= k) & (u < k + 1) & np.isfinite(y))
        out[seg] = y[inside[-n_end:]].mean() - y[inside[:n_end]].mean() if len(inside) >= 2 * n_end else np.nan
    slope = np.gradient(y, u)[trim:len(y) - trim]
    i = int(np.nanargmax(np.abs(slope)))
    out['steepest position'], out['steepest slope'] = float(u[trim + i]), float(slope[i])
    return out


def shape_descriptors(curve, units, *, flat=EFFECT, prominence=EFFECT / 2, inner=.2, min_rho=.9, sharp=.5,
                      gradient=1.):
    """Shape class, position and width of one curve on segment units (G)."""
    y = np.asarray(curve, float)
    u = np.asarray(units, float)
    out = {'shape': 'not measured', 'amplitude': np.nan, 'direction': 0, 'position': np.nan, 'width': np.nan,
           'switch': ''}
    if not np.isfinite(y).all():
        return out
    amp = float(y.max() - y.min())
    out['amplitude'] = amp
    if amp < flat:
        out['shape'] = 'flat'
        return out
    seg, sign = interior_extrema(y[None], u, min_amplitude=flat, min_prominence=prominence, inner=inner)
    if seg[0] >= 0:
        j = int(np.argmax(sign[0] * y))
        out.update(shape='interior peak' if sign[0] > 0 else 'interior trough', direction=int(sign[0]),
                   position=float(u[j]))
        return out
    mono = monotone_onsets(y[None], u, min_amplitude=flat, min_abs_rho=min_rho).iloc[0]
    if bool(mono.monotone):
        width = float(mono.width)
        out.update(shape='rising' if mono.direction > 0 else 'falling', direction=int(mono.direction),
                   position=float(mono.onset), width=width,
                   switch='sharp switch' if width <= sharp else 'gradient' if width >= gradient else 'graded switch')
        return out
    out['shape'] = 'mixed'
    return out


def species_statement(first, second, *, tolerance=EFFECT, weak=(1., 2.), weak_label='human S2: weak cell'):
    """Species-level description from two specimens' descriptors (dicts from ``shape_descriptors``)."""
    shapes = {first['shape'], second['shape']}
    out = {'shape': '', 'position': np.nan, 'width': np.nan, 'switch': '', 'agree': False, 'flag': ''}
    if len(shapes) != 1 or 'not measured' in shapes:
        out['shape'] = 'specimens disagree' if 'not measured' not in shapes else 'not measured'
        return out
    shape = first['shape']
    out['shape'] = shape
    if shape in ('flat', 'mixed'):
        out['agree'] = True
        return out
    if abs(first['position'] - second['position']) > tolerance:
        out['shape'] = f'{shape}, positions disagree'
        return out
    pos = (first['position'] + second['position']) / 2
    out.update(agree=True, position=pos)
    if shape in ('rising', 'falling'):
        width = (first['width'] + second['width']) / 2
        out['width'] = width
        out['switch'] = first['switch'] if first['switch'] == second['switch'] else 'graded switch'
    if weak is not None and weak[0] <= pos < weak[1]:
        out['flag'] = weak_label
    return out


def segment_position(units):
    """'S2 0.62' for 1.62 segment units (S1 starts at 0)."""
    if not np.isfinite(units):
        return ''
    k = int(np.clip(np.floor(units), 0, 2))
    return f'{SEGMENTS[k]} {units - k:.2f}'


def cross_species(mouse, human, *, gap=.37):
    """Compare two species statements (G): same monotone class required; |Δposition| ≥ gap to call a shift."""
    if not (mouse['agree'] and human['agree']):
        return 'not comparable'
    if mouse['shape'] != human['shape']:
        return f'different shape (mouse {mouse["shape"]}, human {human["shape"]})'
    if mouse['shape'] in ('flat', 'mixed'):
        return f'both {mouse["shape"]}'
    delta = human['position'] - mouse['position']
    if abs(delta) < gap:
        return 'same place within registration error'
    return f'{"later" if delta > 0 else "earlier"} in human by {abs(delta):.2f} segment units'


# --- L: pathway label -------------------------------------------------------------------------------------------

def pathway_label(genes, libraries, *, min_genes=3):
    """The library term holding most of ``genes`` (>= min_genes). ``libraries`` is an ordered {library: {term:
    members}}; ties go to the earlier library, then the term name. Returns (library, term, genes in it) or None."""
    genes = set(genes)
    best = None
    for rank, (lib, sets) in enumerate(libraries.items()):
        for term in sorted(sets):
            hit = sorted(genes & set(sets[term]))
            key = (-len(hit), rank, term)
            if len(hit) >= min_genes and (best is None or key < best[0]):
                best = (key, lib, term, hit)
    return None if best is None else best[1:]


def spearman(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    r = spearmanr(x[ok], y[ok])
    return float(r.statistic), float(r.pvalue), int(ok.sum())
