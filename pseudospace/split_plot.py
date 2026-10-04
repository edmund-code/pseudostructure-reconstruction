"""Specimen-level split-plot tests on positional pseudobulk (notebook 37).

Species is a property of whole specimens; position varies within each specimen. In this
split-plot layout, a species offset can only be judged against between-specimen variation. With
2 + 2 specimens that leaves 2 df. A species-by-position difference is judged against the
specimen-by-position deviations instead, which leave 2(k - 1) df for k position bins
(Altman & Krzywinski 2015, Nat Methods 12:377).

limma is not installed in this environment, so this module carries a small Python version of
three limma pieces: voom-style precision weights (Law et al. 2014), empirical-Bayes variance
moderation (fitFDist without a trend; Smyth 2004), and moderated F-tests. It omits limma's
robust and trend options.
"""
from __future__ import annotations

import numpy as np
from scipy import stats
from scipy.special import digamma, polygamma


def log_cpm(counts, library_size):
    """voom's log2-CPM: log2((count + 0.5) / (library + 1) * 1e6). Rows are samples."""
    counts = np.asarray(counts, float)
    library_size = np.asarray(library_size, float)
    if counts.ndim != 2 or library_size.shape != (counts.shape[0],) or (library_size <= 0).any():
        raise ValueError('Need a samples x genes count matrix and positive library sizes.')
    return np.log2((counts + .5) / (library_size[:, None] + 1) * 1e6)


def _ols(y, design):
    beta, *_ = np.linalg.lstsq(design, y, rcond=None)
    fitted = design @ beta
    rank = np.linalg.matrix_rank(design)
    resid_df = design.shape[0] - rank
    if resid_df < 1:
        raise ValueError('The design leaves no residual degrees of freedom.')
    sigma = np.sqrt(((y - fitted) ** 2).sum(axis=0) / resid_df)
    return fitted, sigma


def voom_weights(logcpm, library_size, design, *, span=.5):
    """Observation weights from the mean-variance trend of an unweighted fit (voom).

    sqrt(residual SD) is smoothed against average log2 count; each observation's weight is the
    inverse fourth power of the trend at its fitted log2 count.
    """
    from statsmodels.nonparametric.smoothers_lowess import lowess

    y = np.asarray(logcpm, float)
    design = np.asarray(design, float)
    fitted, sigma = _ols(y, design)
    offset = np.log2(np.asarray(library_size, float) + 1) - np.log2(1e6)
    sx = y.mean(axis=0) + offset.mean()
    sy = np.sqrt(sigma)
    trend = lowess(sy, sx, frac=span, return_sorted=True)
    xs, index = np.unique(trend[:, 0], return_index=True)
    ys = trend[index, 1]
    predicted = np.interp(fitted + offset[:, None], xs, ys)
    return 1 / np.maximum(predicted, 1e-8) ** 4


def _trigamma_inverse(value):
    """Solve trigamma(x) = value for x > 0 (Newton iteration as in limma)."""
    value = float(value)
    if value > 1e7:
        return 1 / np.sqrt(value)
    if value < 1e-6:
        return 1 / value
    x = .5 + 1 / value
    for _ in range(50):
        tri = polygamma(1, x)
        step = tri * (1 - tri / value) / polygamma(2, x)
        x += step
        if -step / x < 1e-8:
            break
    return x


def squeeze_var(s2, df):
    """limma's fitFDist + squeezeVar without a trend. Returns (posterior s2, d0, s0^2)."""
    s2 = np.asarray(s2, float)
    if (s2 <= 0).any() or not np.isfinite(s2).all() or df <= 0:
        raise ValueError('Variances must be positive and finite with positive df.')
    z = np.log(s2)
    e = z - digamma(df / 2) + np.log(df / 2)
    emean = e.mean()
    evar = ((e - emean) ** 2).sum() / (len(e) - 1) - polygamma(1, df / 2)
    if evar > 0:
        d0 = 2 * _trigamma_inverse(evar)
        s02 = np.exp(emean + digamma(d0 / 2) - np.log(d0 / 2))
        posterior = (d0 * s02 + df * s2) / (d0 + df)
    else:
        d0, s02 = np.inf, np.exp(emean)
        posterior = np.full_like(s2, s02)
    return posterior, d0, s02


def moderated_f(y, weights, design, test_columns):
    """Weighted least squares per gene with a moderated F for a block of coefficients.

    ``y`` and ``weights`` are samples x genes; ``test_columns`` index the tested block of
    ``design``. F = beta_I' V_II^-1 beta_I / (q * s2_post) on (q, d0 + d) df, with V the
    unscaled covariance (equivalent to limma's F from moderated t). Returns a dict with F, p,
    the df, raw and posterior variances, and sqrt(weight)-scaled residuals for
    correlation estimates.
    """
    y, w, x = (np.asarray(v, float) for v in (y, weights, design))
    n, p = x.shape
    test = np.asarray(test_columns, int)
    if y.shape != w.shape or y.shape[0] != n or (w <= 0).any() or np.linalg.matrix_rank(x) < p:
        raise ValueError('Need aligned positive weights and a full-rank design.')
    resid_df = n - p
    if resid_df < 1:
        raise ValueError('No residual degrees of freedom.')
    xtwx = np.einsum('np,ng,nq->gpq', x, w, x)
    xtwy = np.einsum('np,ng,ng->gp', x, w, y)
    cov = np.linalg.inv(xtwx)
    beta = np.einsum('gpq,gq->gp', cov, xtwy)
    resid = y - x @ beta.T
    s2 = (w * resid ** 2).sum(axis=0) / resid_df
    s2 = np.maximum(s2, 1e-12)
    posterior, d0, s02 = squeeze_var(s2, resid_df)
    block = cov[:, test][:, :, test]
    b = beta[:, test]
    quad = np.einsum('gp,gp->g', b, np.linalg.solve(block, b[..., None])[..., 0])
    q = len(test)
    total_df = min(d0 + resid_df, 1e6)
    f_stat = quad / (q * posterior)
    return {'F': f_stat, 'p': stats.f.sf(f_stat, q, total_df), 'df1': q, 'df_resid': resid_df,
            'd0': d0, 's02': s02, 's2': s2, 's2_post': posterior, 'beta': beta,
            'weighted_residuals': np.sqrt(w) * resid}
