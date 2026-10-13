"""Audit of the NB gene model (notebook 87): its tests under other settings, and NB diagnostics at unit level.

Protocol: ``docs/results/pt-gene-model-audit.md``.

* ``positional_lrs`` recomputes 71's T_spatial, T_total and LR_spec with explicit knots, weights, offset and
  dispersion, by the same two ``nested_nb_lr`` calls as ``standard_gene_model`` and ``specimen_shape_lr``, so 71's
  settings reproduce 71's statistics and one setting can be changed at a time;
* ``zero_check``, ``pearson_dispersion``, ``quantile_residuals``: is a fitted NB adequate for a gene's counts;
* ``grouped_alpha``: ML dispersions per group of units at fixed fitted means, with the LR against one dispersion;
* ``dispersion_trend`` / ``shrunken_alpha``: the trended and DESeq2-style shrunken dispersions;
* ``poscounts_size_factors``: DESeq2's size factors for sparse counts;
* ``cv_folds`` / ``heldout_deviance``: unit-level cross-validation of the spline's degrees of freedom.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.special import polygamma
from scipy.stats import chi2, kstest, nbinom, norm

from .nb_gam import _loglik, _ml_alpha, nb_dispersion, nb_fit, nested_nb_lr
from .pathway_calibration import contrast_designs
from .pathway_pipelines import specimen_shape_designs


def positional_lrs(counts, library, position, group, specimen, knots, *, dispersion=None, nuisance_shape=None,
                   unit_weights=False, offset=None, block_size=200, n_jobs=1):
    """71's gene statistics of one partition with explicit settings.

    T_level, T_spatial and T_total as ``standard_gene_model`` (designs base, level, full) and LR_spec as
    ``specimen_shape_lr`` (full, spec), at one fixed dispersion. ``dispersion`` None estimates it under this
    partition's M_full (71: the species split; pass it to the relabelings). ``unit_weights`` drops 71's specimen
    balancing; ``offset`` replaces log(library). Returns a dict of per-gene arrays plus df1, df2 and the dispersion.
    """
    s = np.asarray(position, float)
    lib = np.asarray(library, float)
    designs, w = contrast_designs(s, group, specimen, knots, nuisance_shape=nuisance_shape)
    spec, _ = specimen_shape_designs(s, group, specimen, knots, nuisance_shape=nuisance_shape)
    if unit_weights:
        w = np.ones(len(s))
    off = np.log(lib) if offset is None else np.asarray(offset, float)
    run = {'block_size': block_size, 'n_jobs': n_jobs}
    if dispersion is None:
        dispersion = nb_dispersion(counts, designs['full'], w, off, **run)
    model = nested_nb_lr(counts, {k: designs[k] for k in ('base', 'level', 'full')}, w, off, dispersion,
                         {'T_level': ('base', 'level'), 'T_spatial': ('level', 'full'), 'T_total': ('base', 'full')},
                         **run)
    shape = nested_nb_lr(counts, {'full': spec['full'], 'spec': spec['spec']}, w, off, dispersion,
                         {'LR_spec': ('full', 'spec')}, **run)
    return {'T_level': model['T_level'], 'T_spatial': model['T_spatial'], 'T_total': model['T_total'],
            'LR_spec': shape['LR_spec'], 'converged': model['converged'] & shape['converged'],
            'separated': model['separated'] | shape['separated'], 'df1': model['df']['T_spatial'][0],
            'df_total': model['df']['T_total'][0], 'df2': shape['df']['LR_spec'][0], 'dispersion': dispersion}


def _dense(y):
    return y.toarray() if sparse.issparse(y) else np.asarray(y, float)


def zero_check(y, mu, alpha):
    """Observed against expected zeros per gene (columns). Returns observed, NB and Poisson expected fractions, the
    excess-zero z under the NB (observed - expected zeros over sqrt(sum p0 (1 - p0))) and its one-sided p."""
    y, mu = _dense(y), np.asarray(mu, float)
    a = np.asarray(alpha, float)[None, :]
    p0 = np.exp(-np.log1p(a * mu) / a)
    observed = (y == 0).sum(0)
    expected = p0.sum(0)
    z = (observed - expected) / np.sqrt(np.maximum((p0 * (1 - p0)).sum(0), 1e-12))
    n = y.shape[0]
    return pd.DataFrame({'observed': observed / n, 'nb': expected / n, 'poisson': np.exp(-mu).mean(0), 'z': z,
                         'p': norm.sf(z)})


def pearson_dispersion(y, mu, alpha, n_params):
    """Sum of squared Pearson residuals over residual df, under the NB and under a Poisson variance."""
    y, mu = _dense(y), np.asarray(mu, float)
    a = np.asarray(alpha, float)[None, :]
    df = y.shape[0] - n_params
    return pd.DataFrame({'nb': ((y - mu) ** 2 / (mu + a * mu ** 2)).sum(0) / df,
                         'poisson': ((y - mu) ** 2 / mu).sum(0) / df})


def quantile_residuals(y, mu, alpha, rng):
    """Randomized quantile residuals (Dunn & Smyth 1996) of a NB(mu, alpha) fit: units x genes, N(0, 1) when the
    model holds."""
    y, mu = _dense(y), np.asarray(mu, float)
    size = 1. / np.broadcast_to(np.asarray(alpha, float)[None, :], mu.shape)
    prob = size / (size + mu)
    lo = np.where(y > 0, nbinom.cdf(y - 1, size, prob), 0.)
    hi = nbinom.cdf(y, size, prob)
    u = lo + rng.uniform(size=y.shape) * (hi - lo)
    return norm.ppf(np.clip(u, 1e-12, 1 - 1e-12))


def residual_ks(residuals):
    """Per gene (column) Kolmogorov-Smirnov D and p of the residuals against N(0, 1)."""
    out = [kstest(r, 'norm') for r in np.asarray(residuals, float).T]
    return pd.DataFrame({'D': [o.statistic for o in out], 'p': [o.pvalue for o in out]})


def grouped_alpha(y, mu, alpha0, groups, *, weights=None):
    """ML dispersion per group of units at fixed means, and the LR of the grouped against one dispersion.

    ``groups`` labels each unit (row). Returns (DataFrame genes x groups of alpha, DataFrame with LR, df, p, the common
    alpha). The common alpha is re-estimated on the same means, so the LR compares like with like.
    """
    y, mu = _dense(y), np.asarray(mu, float)
    w = np.ones(y.shape[0]) if weights is None else np.asarray(weights, float)
    labels = np.asarray(groups)
    names = list(dict.fromkeys(labels))
    common, _ = _ml_alpha(y, mu, w, np.asarray(alpha0, float))
    ll0 = _loglik(y, mu, common, w)
    per, ll1 = {}, np.zeros(y.shape[1])
    for name in names:
        rows = labels == name
        per[name], _ = _ml_alpha(y[rows], mu[rows], w[rows], common)
        ll1 += _loglik(y[rows], mu[rows], per[name], w[rows])
    lr = np.maximum(2 * (ll1 - ll0), 0.)
    df = len(names) - 1
    return pd.DataFrame(per), pd.DataFrame({'LR': lr, 'df': df, 'p': chi2.sf(lr, df), 'alpha_common': common})


def dispersion_trend(alpha, mean):
    """Least-squares fit of log alpha on the DESeq2 parametric form alpha = a + b / mean (a, b >= 0, grid over a)."""
    alpha, mean = np.asarray(alpha, float), np.asarray(mean, float)
    ok = np.isfinite(alpha) & (alpha > 0) & (mean > 0)
    la, inv = np.log(alpha[ok]), 1. / mean[ok]
    best = None
    for a in np.concatenate([[0.], np.geomspace(1e-4, np.quantile(alpha[ok], .9), 200)]):
        bs = np.geomspace(1e-6, 100, 400)
        loss = ((la[None, :] - np.log(a + bs[:, None] * inv[None, :])) ** 2).sum(1)
        k = int(np.argmin(loss))
        if best is None or loss[k] < best[0]:
            best = (loss[k], a, bs[k])
    _, a, b = best
    return a, b, a + b / np.maximum(mean, 1e-12)


def shrunken_alpha(alpha, trend, n_units, n_params):
    """DESeq2-style MAP on log alpha with a normal prior around the trend.

    Sampling variance of log alpha: trigamma((n - p) / 2); prior variance: the variance of log(alpha / trend) minus it,
    floored at 0.25 (DESeq2's floor). Returns the shrunken alpha and the weight on the gene's own estimate.
    """
    la, lt = np.log(np.asarray(alpha, float)), np.log(np.asarray(trend, float))
    sampling = float(polygamma(1, (n_units - n_params) / 2))
    prior = max(float(np.var(la - lt)) - sampling, .25)
    weight = prior / (prior + sampling)
    return np.exp(weight * la + (1 - weight) * lt), weight


def poscounts_size_factors(counts):
    """DESeq2's ``poscounts`` size factors: per gene the geometric mean over its positive counts (scaled by all units),
    per unit the median over its positive genes of log(count / geometric mean), exponentiated and normalised to a
    geometric mean of 1."""
    y = sparse.csc_matrix(counts, dtype=float)
    n = y.shape[0]
    logs = y.copy()
    logs.data = np.log(logs.data)
    lgm = np.asarray(logs.sum(axis=0)).ravel() / n                    # log of the poscounts geometric mean
    ratios = sparse.csr_matrix(logs)
    ratios.data = ratios.data - lgm[ratios.indices]
    sf = np.array([np.median(ratios.data[ratios.indptr[i]:ratios.indptr[i + 1]])
                   if ratios.indptr[i + 1] > ratios.indptr[i] else np.nan for i in range(n)])
    if np.isnan(sf).any():
        raise ValueError('A unit has no positive count.')
    sf = np.exp(sf)
    return sf / np.exp(np.mean(np.log(sf)))


def cv_folds(specimen, n_folds, rng):
    """Fold label per unit, drawn within each specimen (equal sizes up to one)."""
    spec = np.asarray(specimen).astype(str)
    out = np.empty(len(spec), int)
    for name in np.unique(spec):
        rows = np.flatnonzero(spec == name)
        out[rows] = rng.permutation(np.arange(len(rows)) % n_folds)
    return out


def heldout_deviance(y, design, weights, offset, folds, *, block_size=200, n_jobs=1):
    """Summed held-out NB deviance per gene: per fold, alpha and beta on the training units (weighted), the unweighted
    deviance on the held-out units. Returns (deviance per gene, all folds converged per gene)."""
    y = sparse.csc_matrix(y, dtype=float)
    x = np.asarray(design, float)
    w, off = np.asarray(weights, float), np.asarray(offset, float)
    total, ok = np.zeros(y.shape[1]), np.ones(y.shape[1], bool)
    for k in np.unique(folds):
        train, test = folds != k, folds == k
        disp = nb_dispersion(y[train], x[train], w[train], off[train], block_size=block_size, n_jobs=n_jobs)
        fit = nb_fit(y[train], x[train], w[train], off[train], disp, block_size=block_size, n_jobs=n_jobs)
        mu = np.clip(np.exp(x[test] @ fit['beta'].T + off[test][:, None]), 1e-12, 1e8)
        yt = y[test].toarray()
        a = disp['alpha'][None, :]
        term = (np.where(yt > 0, yt * np.log(np.where(yt > 0, yt, 1.) / mu), 0.)
                - (yt + 1. / a) * (np.log1p(a * yt) - np.log1p(a * mu)))
        total += 2. * term.sum(0)
        ok &= fit['converged']
    return total, ok
