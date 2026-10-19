"""Joining per-gene results from many notebooks into one atlas (notebook 112; docs/results/pt-gene-atlas.md).

* ``join_lines``: a gene x line boolean frame, with a separate testable mask, from (call, universe) pairs.
* ``convergence_counts``: how many genes are called by at least k testable lines, k = 1..n.
* ``stratified_line_null``: the same counts after permuting each line's calls inside expression strata, which keeps
  every line's size and its expression advantage and destroys gene-specific agreement.
* ``line_similarity``: pairwise Jaccard and phi between lines, on the genes testable on both.
* ``yield_profile``: per calling method, the clean / supported / anchored / unique shares and the two yields.
* ``evidence_model``: an L2 logistic fit of "this gene is in a named set" on the line calls (notebook 113), giving
  each line a weight on one scale and a scorer for held-out genes.

Ponytail: lines are binary and unweighted, so a strong and a marginal call count the same; the upgrade is a
weighted count using each line's standardized statistic.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def join_lines(index, lines):
    """(calls, testable): two gene x line boolean frames.

    ``lines`` maps a line name to (called, universe): iterables of gene names. A gene is testable on a line when it
    is in that line's universe; a call outside the universe is kept and marks the gene testable, since a method
    cannot call a gene it did not look at.
    """
    index = pd.Index(index, name='gene')
    calls = pd.DataFrame(False, index=index, columns=list(lines))
    testable = pd.DataFrame(False, index=index, columns=list(lines))
    for name, (called, universe) in lines.items():
        called = index.isin(pd.Index(called))
        testable[name] = index.isin(pd.Index(universe)) | called
        calls[name] = called & testable[name].to_numpy()
    return calls, testable


def convergence_counts(calls, testable=None, *, kmax=None):
    """Series k -> number of genes called by at least k lines (genes with no testable line are dropped)."""
    c = calls.to_numpy(bool)
    if testable is not None:
        keep = testable.to_numpy(bool).any(axis=1)
        c = c[keep]
    n = c.sum(axis=1)
    kmax = calls.shape[1] if kmax is None else kmax
    return pd.Series({k: int((n >= k).sum()) for k in range(1, kmax + 1)}, name='genes')


def stratified_line_null(calls, testable, strata, *, draws=2000, seed=0, kmax=None):
    """DataFrame draw x k of ``convergence_counts`` after permuting each line inside ``strata``.

    Within one line and one stratum the called genes are redrawn, without replacement, from the testable genes of
    that stratum, keeping the number of calls. Genes testable on no line are dropped, as in the observed count.
    """
    rng = np.random.default_rng(seed)
    strata = pd.Series(strata).reindex(calls.index).to_numpy()
    c, t = calls.to_numpy(bool), testable.to_numpy(bool)
    keep = t.any(axis=1)
    kmax = calls.shape[1] if kmax is None else kmax
    groups = [np.flatnonzero((strata == s) & keep) for s in pd.unique(strata[keep])]
    out = np.zeros((draws, kmax), int)
    for d in range(draws):
        drawn = np.zeros_like(c)
        for j in range(c.shape[1]):
            for g in groups:
                pool = g[t[g, j]]
                k = int(c[g, j].sum())
                if k and len(pool):
                    drawn[rng.choice(pool, min(k, len(pool)), replace=False), j] = True
        n = drawn[keep].sum(axis=1)
        out[d] = [(n >= k).sum() for k in range(1, kmax + 1)]
    return pd.DataFrame(out, columns=range(1, kmax + 1))


def line_similarity(calls, testable):
    """(jaccard, phi): line x line frames computed on the genes testable on both lines."""
    names = list(calls.columns)
    c, t = calls.to_numpy(bool), testable.to_numpy(bool)
    jac = pd.DataFrame(np.nan, index=names, columns=names)
    phi = pd.DataFrame(np.nan, index=names, columns=names)
    for i in range(len(names)):
        for j in range(len(names)):
            both = t[:, i] & t[:, j]
            a, b = c[both, i], c[both, j]
            union = (a | b).sum()
            jac.iloc[i, j] = (a & b).sum() / union if union else np.nan
            n11, n10, n01, n00 = (a & b).sum(), (a & ~b).sum(), (~a & b).sum(), (~a & ~b).sum()
            den = np.sqrt(float(n11 + n10) * (n01 + n00) * (n11 + n01) * (n10 + n00))
            phi.iloc[i, j] = (n11 * n00 - n10 * n01) / den if den > 0 else np.nan
    return jac, phi


def yield_profile(calls, testable, methods, *, clean, supported, anchored, min_lines=2):
    """One row per method in ``methods`` with the shares and counts of the frozen yield criteria.

    ``clean``, ``supported`` and ``anchored`` are boolean Series over the atlas index.
    """
    clean, supported, anchored = (pd.Series(s).reindex(calls.index).fillna(False).astype(bool)
                                  for s in (clean, supported, anchored))
    n_lines = calls.sum(axis=1)
    rows = []
    for m in methods:
        called = calls[m]
        others = calls[[o for o in methods if o != m]].any(axis=1)
        cl = called & clean
        su = cl & supported
        an = su & anchored
        rows.append({
            'method': m,
            'called': int(called.sum()),
            'clean': int(cl.sum()),
            'clean share': cl.sum() / called.sum() if called.sum() else np.nan,
            'supported': int(su.sum()),
            'supported share': su.sum() / cl.sum() if cl.sum() else np.nan,
            'anchored share': an.sum() / su.sum() if su.sum() else np.nan,
            'unique share': (called & ~others).sum() / called.sum() if called.sum() else np.nan,
            'yield': int(su.sum()),
            'novel yield': int((su & ~anchored & (n_lines >= min_lines)).sum()),
        })
    return pd.DataFrame(rows)


def evidence_model(features, outcome, *, C=1.0, seed=0):
    """(coefficients, score): an L2 logistic fit of ``outcome`` on the numeric frame ``features``.

    ``score`` takes a frame with the same columns and returns the linear predictor, so a model fitted on one set of
    genes can rank another. Columns that are constant in the training genes get a zero coefficient.
    """
    from sklearn.linear_model import LogisticRegression

    cols = list(features.columns)
    X = features[cols].to_numpy(float)
    y = np.asarray(outcome, int)
    keep = X.std(axis=0) > 0
    fit = LogisticRegression(C=C, solver='liblinear', max_iter=2000,      # liblinear is L2 by default
                             random_state=seed).fit(X[:, keep], y)
    coef = pd.Series(0., index=cols)
    coef[np.array(cols)[keep]] = fit.coef_[0]

    def score(frame):
        return frame[cols].to_numpy(float) @ coef.to_numpy() + fit.intercept_[0]

    return coef, score
