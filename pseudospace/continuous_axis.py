"""Six uses of the continuous PT axis from the literature review (notebook 90).

Protocol: ``docs/results/pt-continuous-axis-experiments.md``.

- ``label_fractions``, ``fit_ladder``, ``ladder_class``, ``ladder_curve``, ``half_point``: A, gene shape classes with
  segment steps (constant within each label) as a competing model.
- ``centre_of_mass``, ``aligned_amplitude``: B, set positions and sub-trends a median curve flattens.
- ``nb_bin_loglik``, ``log_posterior``, ``plugin_mi``, ``corrected_mi``, ``within_label_mi``: C, decoding position.
- ``end_contrast``, ``reliability``, ``disattenuated_r``, ``graded_categories``: D, level-free species comparison.
- ``local_change``, ``discrete_profiles``: E, change density against a discrete transcriptome.
- ``matched_fold_changes``: F, segment-level against position-matched fold changes.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.special import logsumexp

from .transition_sharpness import binned_log_rates

EFFECT = .25
BIC_MARGIN = 6.
LN81 = np.log(81.)
MODELS = ('F', 'L', 'S', 'B', 'D', 'DS')
CONTINUOUS = ('L', 'S', 'B')


def label_fractions(bins, labels, n_bins, label_order):
    """Bins x labels: the share of the units in each bin that carry each label (NaN rows for empty bins)."""
    bins, labels = np.asarray(bins), np.asarray(labels)
    out = np.zeros((n_bins, len(label_order)))
    keep = bins >= 0
    for j, lab in enumerate(label_order):
        out[:, j] = np.bincount(bins[keep & (labels == lab)], minlength=n_bins)[:n_bins]
    total = out.sum(axis=1, keepdims=True)
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(total > 0, out / total, np.nan)


def _sig(u, m, w):
    return 1. / (1. + np.exp(-np.clip((u - m) * LN81 / w, -60, 60)))


def _gauss(u, c, fwhm):
    return np.exp(-4 * np.log(2) * (u - c) ** 2 / fwhm ** 2)


def _stack(profiles, weights, label_frac, centres):
    """Rows of the joint fit: finite bins of every mouse. Returns y (genes x rows), sqrt weights, mouse index, u and
    the label fractions of each row."""
    ys, ws, ms, us, fs = [], [], [], [], []
    for k, (P, w, f) in enumerate(zip(profiles, weights, label_frac)):
        P = np.asarray(P, float)
        ok = np.isfinite(P).all(axis=0) & np.isfinite(np.asarray(f, float)).all(axis=1) & (np.asarray(w, float) > 0)
        ys.append(P[:, ok])
        ws.append(np.asarray(w, float)[ok])
        ms.append(np.full(ok.sum(), k))
        us.append(np.asarray(centres, float)[ok])
        fs.append(np.asarray(f, float)[ok])
    w = np.concatenate(ws)
    return np.hstack(ys), np.sqrt(w / w.mean()), np.concatenate(ms), np.concatenate(us), np.vstack(fs)


def _intercepts(mouse):
    return np.column_stack([(mouse == k).astype(float) for k in np.unique(mouse)])


def _discrete_columns(mouse, u, frac, slopes):
    """Label columns beyond the first (the intercepts absorb it) and, with ``slopes``, one within-label slope each."""
    cols = [frac[:, j] for j in range(1, frac.shape[1])]
    if slopes:
        for j in range(frac.shape[1]):
            share = frac[:, j]
            centre = np.sum(share * u) / max(share.sum(), 1e-12)
            cols.append(share * (u - centre))
    return np.column_stack(cols)


def _solve(X, sw, Y):
    """Weighted least squares of every row of Y on X: (WRSS per row, coefficients rows x k)."""
    Xw = X * sw[:, None]
    Yw = Y * sw[None, :]
    q, r = np.linalg.qr(Xw)
    proj = Yw @ q
    rss = np.maximum((Yw ** 2).sum(axis=1) - (proj ** 2).sum(axis=1), 0.)
    coef = np.linalg.solve(r, proj.T).T
    return rss, coef


def fit_ladder(profiles, weights, centres, label_frac, *, m_grid=None, w_grid=None, c_grid=None, fwhm_grid=None):
    """Fit the six shape models of protocol A to every gene.

    ``profiles``, ``weights`` and ``label_frac`` are lists with one entry per mouse (genes x bins, bins, bins x labels).
    A single-mouse fit is a list of length one. Returns a DataFrame with one row per gene: BIC of each model, the
    nonlinear parameters and amplitudes, and the range of each model's fitted shape (the mean over mice)."""
    centres = np.asarray(centres, float)
    m_grid = np.linspace(.15, 2.35, 45) if m_grid is None else np.asarray(m_grid, float)
    w_grid = np.geomspace(.1, 4., 25) if w_grid is None else np.asarray(w_grid, float)
    c_grid = np.linspace(.15, 2.35, 45) if c_grid is None else np.asarray(c_grid, float)
    fwhm_grid = np.geomspace(.2, 1.2, 25) if fwhm_grid is None else np.asarray(fwhm_grid, float)
    Y, sw, mouse, u, frac = _stack(profiles, weights, label_frac, centres)
    n_genes, n = Y.shape
    n_mice = len(np.unique(mouse))
    I = _intercepts(mouse)
    out = {}

    def shape_range(col_values):
        """Mean over mice of the range of a fitted shape (rows x genes)."""
        return np.mean([np.ptp(col_values[mouse == k], axis=0) for k in range(n_mice)], axis=0)

    rss, _ = _solve(I, sw, Y)
    out['F'] = (rss, n_mice)
    out['range F'] = np.zeros(n_genes)
    rss, coef = _solve(np.column_stack([I, u]), sw, Y)
    out['L'] = (rss, n_mice + 1)
    out['L slope'] = coef[:, -1]
    out['range L'] = np.abs(coef[:, -1]) * np.mean([np.ptp(u[mouse == k]) for k in range(n_mice)])

    for name, grid_a, grid_b, fn in (('S', m_grid, w_grid, _sig), ('B', c_grid, fwhm_grid, _gauss)):
        best = np.full(n_genes, np.inf)
        par = np.full((n_genes, 3), np.nan)
        for a in grid_a:
            for b in grid_b:
                col = fn(u, a, b)
                if np.ptp(col) < 1e-6:
                    continue
                rss, coef = _solve(np.column_stack([I, col]), sw, Y)
                better = rss < best
                best[better] = rss[better]
                par[better] = np.column_stack([np.full(better.sum(), a), np.full(better.sum(), b), coef[better, -1]])
        out[name] = (best, n_mice + 3)
        out[f'{name} a'], out[f'{name} b'], out[f'{name} amp'] = par.T
        rng_shape = np.array([np.mean([np.ptp(fn(u[mouse == k], pa, pb)) for k in range(n_mice)]) if np.isfinite(pa)
                              else np.nan for pa, pb in zip(par[:, 0], par[:, 1])])
        out[f'range {name}'] = np.abs(par[:, 2]) * rng_shape

    for name, slopes in (('D', False), ('DS', True)):
        cols = _discrete_columns(mouse, u, frac, slopes)
        rss, coef = _solve(np.column_stack([I, cols]), sw, Y)
        out[name] = (rss, n_mice + cols.shape[1])
        fitted = cols @ coef[:, n_mice:].T
        out[f'range {name}'] = shape_range(fitted)
        if name == 'D':
            out['D S2-S1'], out['D S3-S1'] = coef[:, n_mice], coef[:, n_mice + 1]
        else:
            out['DS slopes'] = coef[:, n_mice + frac.shape[1] - 1:]
            out['DS jumps'] = coef[:, n_mice:n_mice + frac.shape[1] - 1]

    table = pd.DataFrame(index=np.arange(n_genes))
    for name in MODELS:
        rss, k = out[name]
        table[f'bic {name}'] = n * np.log(np.maximum(rss, 1e-300) / n) + k * np.log(n)
        table[f'range {name}'] = out[f'range {name}']
    table['n'] = n
    table['L slope'] = out['L slope']
    for name in ('S', 'B'):
        table[f'{name} position'], table[f'{name} width'], table[f'{name} amp'] = (out[f'{name} a'], out[f'{name} b'],
                                                                                 out[f'{name} amp'])
    table['D S2-S1'], table['D S3-S1'] = out['D S2-S1'], out['D S3-S1']
    table['DS jump'] = np.abs(out['DS jumps'][:, 0]) + np.abs(out['DS jumps'][:, 1] - out['DS jumps'][:, 0])
    table['DS within slope'] = np.abs(out['DS slopes']).mean(axis=1)
    return table


def ladder_class(table, *, flat=EFFECT, margin=BIC_MARGIN, edge_window=.15, sharp=.5, gradual=1.):
    """Class of every gene from :func:`fit_ladder`'s table: the lowest-BIC model, 'flat' when that is F or its fitted
    range is below ``flat``. Adds ΔBIC continuous-vs-discrete (BIC(D) − best continuous BIC), its reading, and the
    location and width of switches and bumps."""
    t = table.copy()
    bic = t[[f'bic {m}' for m in MODELS]].to_numpy()
    best = np.array(MODELS)[np.argmin(bic, axis=1)]
    rng_best = np.array([t.at[i, f'range {b}'] for i, b in zip(t.index, best)])
    t['best'] = best
    t['range best'] = rng_best
    flat_mask = (best == 'F') | ~(rng_best >= flat)
    names = {'L': 'gradient', 'S': 'switch', 'D': 'discrete', 'DS': 'discrete + within'}
    cls = []
    for b, amp in zip(best, t['B amp']):
        cls.append('bump' if b == 'B' and amp > 0 else 'trough' if b == 'B' else names.get(b, 'flat'))
    t['class'] = np.where(flat_mask, 'flat', cls)
    cont = t[[f'bic {m}' for m in CONTINUOUS]].min(axis=1)
    t['dBIC continuous'] = t['bic D'] - cont
    t['continuous model'] = np.array(CONTINUOUS)[np.argmin(t[[f'bic {m}' for m in CONTINUOUS]].to_numpy(), axis=1)]
    reading = np.where(t['dBIC continuous'] > margin, 'continuous',
                       np.where(t['dBIC continuous'] < -margin, 'discrete', 'undecided'))
    t['continuous vs discrete'] = np.where(t['class'].eq('flat'), 'flat', reading)
    loc = []
    for c, ps, pb in zip(t['class'], t['S position'], t['B position']):
        p = ps if c == 'switch' else pb if c in ('bump', 'trough') else np.nan
        if not np.isfinite(p):
            loc.append('')
            continue
        near = [k for k in (1, 2) if abs(p - k) <= edge_window]
        loc.append(f'at S{near[0]}/S{near[0] + 1}' if near else f'inside S{int(np.clip(np.floor(p), 0, 2)) + 1}')
    t['location'] = loc
    t['switch width'] = np.where(t['class'].eq('switch'),
                                 pd.cut(t['S width'], [0, sharp, gradual - 1e-12, np.inf],
                                        labels=['sharp', 'graded', 'gradual']).astype(str), '')
    return t


def ladder_curve(model, row, centres, label_frac=None):
    """The fitted shape (no intercept) of one gene's ``model`` on ``centres``; ``row`` is its ladder-table row.
    D needs ``label_frac`` (bins x labels) of the mouse whose curve is drawn."""
    u = np.asarray(centres, float)
    if model == 'F':
        return np.zeros_like(u)
    if model == 'L':
        return row['L slope'] * (u - u.mean())
    if model == 'S':
        return row['S amp'] * _sig(u, row['S position'], row['S width'])
    if model == 'B':
        return row['B amp'] * _gauss(u, row['B position'], row['B width'])
    if model == 'D':
        f = np.asarray(label_frac, float)
        return f[:, 1] * row['D S2-S1'] + f[:, 2] * row['D S3-S1']
    raise ValueError(f'No curve for model {model}.')


def half_point(centres, curve):
    """Position where a curve first crosses halfway between its values at the two ends (NaN for a flat curve); for a
    curve with an interior extremum larger than both ends, the extremum's position."""
    u, y = np.asarray(centres, float), np.asarray(curve, float)
    ok = np.isfinite(y)
    u, y = u[ok], y[ok]
    if len(y) < 3 or np.ptp(y) <= 0:
        return np.nan
    lo_end, hi_end = y[0], y[-1]
    j = int(np.argmax(np.abs(y - (lo_end + hi_end) / 2)))
    if abs(y[j] - (lo_end + hi_end) / 2) > abs(hi_end - lo_end) / 2 + 1e-12:
        return float(u[j])
    target = (lo_end + hi_end) / 2
    s = np.sign(y - target)
    cross = np.flatnonzero(s[:-1] * s[1:] <= 0)
    if not len(cross):
        return np.nan
    i = cross[0]
    if y[i + 1] == y[i]:
        return float(u[i])
    return float(u[i] + (target - y[i]) * (u[i + 1] - u[i]) / (y[i + 1] - y[i]))


def centre_of_mass(deviation, centres):
    """Per row: Σ u·2^d / Σ 2^d over finite bins (positions in ``centres``)."""
    D = np.atleast_2d(np.asarray(deviation, float))
    e = np.where(np.isfinite(D), np.exp2(np.nan_to_num(D)), 0.)
    u = np.asarray(centres, float)[None, :]
    return (e * u).sum(axis=1) / e.sum(axis=1)


def pc1_signs(profiles):
    """Signs of the members' loadings on the first principal component of their standardised profiles
    (members x bins), oriented so that the aligned median rises along the axis."""
    P = np.asarray(profiles, float)
    P = P[:, np.isfinite(P).all(axis=0)]
    sd = P.std(axis=1, keepdims=True)
    Z = (P - P.mean(axis=1, keepdims=True)) / np.where(sd > 0, sd, 1.)
    u, s, vt = np.linalg.svd(Z, full_matrices=False)
    signs = np.sign(u[:, 0])
    signs[signs == 0] = 1.
    trend = np.median(signs[:, None] * P, axis=0)
    if np.corrcoef(trend, np.arange(len(trend)))[0, 1] < 0:
        signs = -signs
    return signs


def aligned_amplitude(train, test):
    """Cross-validated aligned amplitude: signs from ``train`` (members x bins), the range of the median sign-aligned
    profile in ``test``. Returns (amplitude, signs)."""
    signs = pc1_signs(train)
    aligned = np.nanmedian(signs[:, None] * np.asarray(test, float), axis=0)
    return float(np.nanmax(aligned) - np.nanmin(aligned)), signs


def nb_bin_loglik(counts, library, rates, alpha, *, min_alpha=1e-4):
    """Units x bins NB log-likelihood (bin-dependent terms only) of each unit's ``counts`` (units x genes) under mean
    library x rate_gb, with per-gene dispersion ``alpha``."""
    C = counts.toarray() if sparse.issparse(counts) else np.asarray(counts, float)
    L = np.asarray(library, float)
    R = np.maximum(np.asarray(rates, float), 1e-12)
    size = 1. / np.maximum(np.asarray(alpha, float), min_alpha)
    out = C @ np.log(R)
    cr = C + size[None, :]
    for b in range(R.shape[1]):
        out[:, b] -= (cr * np.log(size[None, :] + L[:, None] * R[None, :, b])).sum(axis=1)
    return out


def log_posterior(loglik, log_prior):
    """Normalised log posterior over bins (units x bins); ``log_prior`` is bins or units x bins."""
    lp = np.asarray(loglik, float) + np.asarray(log_prior, float)
    return lp - logsumexp(lp, axis=1, keepdims=True)


def plugin_mi(x, y, nx, ny):
    """Plug-in mutual information (bits) of two integer label vectors."""
    x, y = np.asarray(x, int), np.asarray(y, int)
    joint = np.zeros((nx, ny))
    np.add.at(joint, (x, y), 1.)
    joint /= joint.sum()
    px, py = joint.sum(axis=1, keepdims=True), joint.sum(axis=0, keepdims=True)
    nz = joint > 0
    return float((joint[nz] * np.log2(joint[nz] / (px @ py)[nz])).sum())


def corrected_mi(x, y, nx, ny, rng, n_shuffle=20):
    """Plug-in MI minus the mean MI after shuffling ``y`` (the plug-in bias)."""
    raw = plugin_mi(x, y, nx, ny)
    y = np.asarray(y, int)
    bias = np.mean([plugin_mi(x, rng.permutation(y), nx, ny) for _ in range(n_shuffle)])
    return raw - bias


def within_label_mi(x, y, labels, nx, ny, rng, n_shuffle=20):
    """Σ_l p(l)·MI(x; y | label l), each term shuffle-corrected."""
    labels = np.asarray(labels)
    total = 0.
    for lab in np.unique(labels):
        r = labels == lab
        if r.sum() < 10:
            continue
        total += r.mean() * corrected_mi(np.asarray(x)[r], np.asarray(y)[r], nx, ny, rng, n_shuffle)
    return total


def end_contrast(deviation, n_end=8):
    """Per row: mean over the last ``n_end`` finite bins minus the first ``n_end``."""
    D = np.atleast_2d(np.asarray(deviation, float))
    ok = np.flatnonzero(np.isfinite(D).all(axis=0))
    return D[:, ok[-n_end:]].mean(axis=1) - D[:, ok[:n_end]].mean(axis=1)


def reliability(r):
    """Spearman–Brown reliability of the mean of two replicates whose correlation is ``r``."""
    return 2 * r / (1 + r)


def disattenuated_r(r_between, r_a, r_b):
    """Correlation of two two-replicate means corrected for each side's replicate noise."""
    return r_between / np.sqrt(reliability(r_a) * reliability(r_b))


def graded_categories(e1a, e1b, e2a, e2b, *, graded=.5, flat=.2):
    """Per gene: 'conserved', 'reversed', 'first only', 'second only' or '' from the end contrasts of two groups
    of two specimens each (protocol D)."""
    e1a, e1b, e2a, e2b = (np.asarray(v, float) for v in (e1a, e1b, e2a, e2b))
    g1 = (np.abs(e1a) >= graded) & (np.abs(e1b) >= graded) & (np.sign(e1a) == np.sign(e1b))
    g2 = (np.abs(e2a) >= graded) & (np.abs(e2b) >= graded) & (np.sign(e2a) == np.sign(e2b))
    f1 = (np.abs(e1a) < flat) & (np.abs(e1b) < flat)
    f2 = (np.abs(e2a) < flat) & (np.abs(e2b) < flat)
    out = np.full(len(e1a), '', dtype=object)
    out[g1 & g2 & (np.sign(e1a) == np.sign(e2a))] = 'conserved'
    out[g1 & g2 & (np.sign(e1a) != np.sign(e2a))] = 'reversed'
    out[g1 & f2] = 'first only'
    out[g2 & f1] = 'second only'
    return out


def local_change(deviation, half=2):
    """Genes x positions: mean of the ``half`` bins after each interior bin edge minus the ``half`` bins before it.
    Position j (j = half .. n_bins - half) sits between bins j - 1 and j."""
    D = np.asarray(deviation, float)
    n = D.shape[1]
    return np.stack([D[:, j:j + half].mean(axis=1) - D[:, j - half:j].mean(axis=1) for j in range(half, n - half + 1)],
                    axis=1)


def discrete_profiles(counts, library, labels, bins, n_bins, specimen, names, label_order, *, min_units=20, pseudo=.5):
    """Per specimen: the genes x bins log2 rates a constant-within-label gene would give (each bin Σ_l f_bl·R_l, R_l
    the specimen's label pseudobulk rate), centred like ``positional_programs.bin_profiles``."""
    counts = sparse.csr_matrix(counts)
    labels, bins, specimen = np.asarray(labels), np.asarray(bins), np.asarray(specimen)
    library = np.asarray(library, float)
    out = {}
    for name in names:
        rows = (specimen == name) & (bins >= 0)
        frac = label_fractions(bins[rows], labels[rows], n_bins, label_order)
        R = np.column_stack([(np.asarray(counts[rows & (labels == lab)].sum(axis=0)).ravel() + pseudo)
                             / max(library[rows & (labels == lab)].sum(), 1e-12) for lab in label_order])
        y = np.log2(R @ np.nan_to_num(frac).T)
        valid = np.bincount(bins[rows], minlength=n_bins)[:n_bins] >= min_units
        y[:, ~valid] = np.nan
        out[name] = y - np.nanmean(y, axis=1, keepdims=True)
    return out


def matched_fold_changes(counts, library, bins, labels, condition, n_bins, label, *, case='case', ref='ref',
                         pseudo=.5):
    """Per gene, for the units carrying ``label``: the pooled log2 FC case − ref (what a segment pseudobulk reports)
    and the position-matched log2 FC (both sides' binned rates weighted by the case units' library share per bin,
    over bins holding units of both). Returns (pooled, matched, case mean position − ref mean position in bins)."""
    counts = sparse.csr_matrix(counts)
    bins, labels, condition = np.asarray(bins), np.asarray(labels), np.asarray(condition)
    library = np.asarray(library, float)
    sel = (labels == label) & (bins >= 0)

    def side(group):
        r = sel & (condition == group)
        with np.errstate(divide='ignore'):
            c = binned_log_rates(counts[r], library[r], bins[r], n_bins, pseudo=pseudo)
        lib = np.bincount(bins[r], weights=library[r], minlength=n_bins)[:n_bins]
        pooled = np.log2((np.asarray(counts[r].sum(axis=0)).ravel() + pseudo) / library[r].sum())
        return c, lib, pooled, float(np.average(bins[r], weights=library[r]))

    c_case, lib_case, pooled_case, pos_case = side(case)
    c_ref, lib_ref, pooled_ref, pos_ref = side(ref)
    both = (lib_case > 0) & (lib_ref > 0)
    w = lib_case[both] / lib_case[both].sum()
    matched = np.log2(np.exp2(c_case[:, both]) @ w) - np.log2(np.exp2(c_ref[:, both]) @ w)
    return pooled_case - pooled_ref, matched, pos_case - pos_ref
