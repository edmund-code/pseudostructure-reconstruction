"""Dilution-free comparison of axial zonation amplitude between species (notebook 39).

Our human and mouse segment gradients (for example S2 - S1) are both regressed on an
**independent** reference gradient measured in other animals or donors. Each slope is the
covariance with the reference divided by the reference variance. Reference noise therefore
attenuates both slopes by the same factor, and their ratio cov(human, ref) / cov(mouse, ref) is
free of regression dilution. Regressing human on our own mouse gradient is not: our mouse noise
sits in the predictor, so that slope shrinks toward zero and human looks flat.

A mouse reference biases the ratio downward when zonation differs between species, because genes
zonated only in mouse inflate the mouse covariance. A human reference biases it upward for the same
reason. The two kinds of reference therefore bracket the global amplitude ratio.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse


def group_sums(counts, labels):
    """Sum rows of a (structures x genes) count matrix per label. Returns (labels, sums)."""
    labels = np.asarray(labels, dtype=str)
    names, codes = np.unique(labels, return_inverse=True)
    indicator = sparse.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))), shape=(len(names), len(codes)))
    matrix = counts if sparse.issparse(counts) else np.asarray(counts, float)
    sums = indicator @ matrix
    return names, np.asarray(sums.toarray() if sparse.issparse(sums) else sums, float)


def log_ratio(counts_late, library_late, counts_early, library_early, *, prior=.5):
    """log2 of the CPM ratio late / early from summed counts (a pseudobulk log fold change)."""
    late = np.log2((np.asarray(counts_late, float) + prior) / (float(library_late) + 1) * 1e6)
    early = np.log2((np.asarray(counts_early, float) + prior) / (float(library_early) + 1) * 1e6)
    return late - early


def ols_slope(x, y):
    """Least-squares slope of y on x with an intercept."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    xc = x - x.mean()
    return float(xc @ (y - y.mean()) / (xc @ xc))


def deming_slope(x, y, delta):
    """Deming regression slope of y on x; delta = error variance of y / error variance of x.

    delta = 1 is orthogonal (total least squares) regression; delta -> inf recovers OLS.
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    sxx, syy = np.var(x, ddof=1), np.var(y, ddof=1)
    sxy = np.cov(x, y, ddof=1)[0, 1]
    if sxy == 0:
        return 0.
    d = float(delta)
    return float((syy - d * sxx + np.sqrt((syy - d * sxx) ** 2 + 4 * d * sxy ** 2)) / (2 * sxy))


def amplitude_ratio(human, mouse, reference, *, method='ols', delta_human=None, delta_mouse=None):
    """Slopes of our human and mouse gradients on a reference gradient, and their ratio.

    ``method``: 'ols' (dilution-free ratio), 'deming' (needs delta_human and delta_mouse, the
    error-variance ratios of each response to the reference) or 'orthogonal' (delta = 1).
    """
    if method == 'ols':
        b_h, b_m = ols_slope(reference, human), ols_slope(reference, mouse)
    elif method == 'deming':
        if delta_human is None or delta_mouse is None:
            raise ValueError('Deming regression needs both error-variance ratios.')
        b_h, b_m = deming_slope(reference, human, delta_human), deming_slope(reference, mouse, delta_mouse)
    elif method == 'orthogonal':
        b_h, b_m = deming_slope(reference, human, 1.), deming_slope(reference, mouse, 1.)
    else:
        raise ValueError(f'Unknown method: {method}')
    return {'slope_human': b_h, 'slope_mouse': b_m, 'ratio': b_h / b_m if b_m != 0 else np.nan}


def bootstrap_ratio(human, mouse, reference_draw, *, n_boot=2000, seed=0, method='ols', deltas=None):
    """Percentile interval of the amplitude ratio, resampling genes and (optionally) reference donors.

    ``human`` and ``mouse`` are gene-indexed Series. ``reference_draw(rng)`` returns a
    gene-indexed Series of the reference gradient over the selected genes. Return the observed
    data with ``rng=None``, and a donor-resampled version otherwise. Genes are resampled with
    replacement within each draw. Returns (2.5%, 50%, 97.5%) and the draws.
    """
    rng = np.random.default_rng(seed)
    ratios = []
    for _ in range(n_boot):
        reference = reference_draw(rng).dropna()
        genes = reference.index.intersection(human.index).intersection(mouse.index)
        pick = genes[rng.integers(0, len(genes), len(genes))]
        kwargs = {} if deltas is None else {'delta_human': deltas[0], 'delta_mouse': deltas[1]}
        ratios.append(amplitude_ratio(human.loc[pick].to_numpy(), mouse.loc[pick].to_numpy(),
                                      reference.loc[pick].to_numpy(), method=method, **kwargs)['ratio'])
    ratios = np.asarray(ratios, float)
    return np.nanpercentile(ratios, [2.5, 50, 97.5]), ratios


def flattening_class(human_pairs, mouse_pairs, ratio, *, residual=.5, signed=.3, specific=.5, absent=.25):
    """Classify one gene contrast against global flattening (human = ratio x mouse).

    ``human_pairs`` / ``mouse_pairs``: per-specimen gradients (log2). Beyond flattening needs
    |mean residual| >= ``residual`` and the same residual sign in every human-mouse specimen
    pairing. Subclasses: reversal (opposite signs, both |gradient| >= ``signed``),
    human-specific gradient (|human| >= ``specific`` while |mouse| < ``absent``),
    human-stronger (same sign, |human| > |mouse|), other.
    """
    h, m = np.asarray(human_pairs, float), np.asarray(mouse_pairs, float)
    if not (np.isfinite(h).all() and np.isfinite(m).all()):
        return 'not testable', np.nan
    residuals = np.array([a - ratio * b for a in h for b in m])
    mean_h, mean_m = h.mean(), m.mean()
    mean_residual = mean_h - ratio * mean_m
    consistent = np.all(np.sign(residuals) == np.sign(mean_residual))
    if abs(mean_residual) < residual or not consistent:
        return 'consistent with flattening', mean_residual
    if np.sign(mean_h) != np.sign(mean_m) and abs(mean_h) >= signed and abs(mean_m) >= signed:
        return 'reversal', mean_residual
    if abs(mean_h) >= specific and abs(mean_m) < absent:
        return 'human-specific gradient', mean_residual
    if np.sign(mean_h) == np.sign(mean_m) and abs(mean_h) > abs(mean_m):
        return 'human-stronger', mean_residual
    return 'other beyond flattening', mean_residual
