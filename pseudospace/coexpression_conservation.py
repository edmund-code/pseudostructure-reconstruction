"""Co-expression conservation (notebook 97; docs/results/pt-coexpression-conservation.md).

Crow, Suresh, Lee & Gillis 2022 (NAR, PMID 35451481): for gene g, take its top k co-expression partners in one network
and ask how well g's co-expression ranks in another network recover them (AUROC); average both directions.

* ``spearman_network``: genes x genes Spearman correlation over samples (constant genes give 0 off the diagonal).
* ``rank_standardise`` and ``aggregate``: Crow's aggregation, every off-diagonal edge ranked over the network and divided
  by the maximum rank, then averaged over networks.
* ``partner_auroc``: per gene, the AUROC of A's top-k partners within B's row (self excluded; average ranks).
* ``conservation``: the mean of ``partner_auroc(A, B)`` and ``partner_auroc(B, A)``.
* ``binned_residuals``: expression minus the mean of its position bin (equal-count bins), per column.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import rankdata


def spearman_network(X):
    """Samples x genes -> genes x genes Spearman correlation; a constant gene correlates 0 with every other gene."""
    R = rankdata(np.asarray(X, float), axis=0)
    R -= R.mean(axis=0)
    norm = np.sqrt((R ** 2).sum(axis=0))
    R = np.divide(R, norm, out=np.zeros_like(R), where=norm > 0)
    C = R.T @ R
    np.fill_diagonal(C, 1.)
    return C


def rank_standardise(net):
    """Off-diagonal edges ranked over the whole (symmetric) network and divided by the maximum rank; diagonal 1."""
    A = np.asarray(net, float)
    iu = np.triu_indices_from(A, k=1)
    r = rankdata(A[iu])
    out = np.ones_like(A)
    out[iu] = r / r.max()
    out[(iu[1], iu[0])] = out[iu]
    return out


def aggregate(nets):
    """Mean of rank-standardised networks."""
    return np.mean([rank_standardise(n) for n in nets], axis=0)


def partner_auroc(A, B, *, k=10):
    """Per gene, the AUROC with which B's row ranks A's top-k partners above the other genes (self excluded)."""
    A, B = np.asarray(A, float).copy(), np.asarray(B, float).copy()
    n = A.shape[0]
    if A.shape != B.shape or n < k + 2:
        raise ValueError('Need two square networks over the same genes, larger than k + 1.')
    np.fill_diagonal(A, -np.inf)
    np.fill_diagonal(B, -np.inf)
    partners = np.argpartition(-A, k, axis=1)[:, :k]
    ranks = rankdata(B, axis=1) - 1.              # the diagonal holds rank 1; the others rank 1..n-1 among n-1 genes
    got = np.take_along_axis(ranks, partners, axis=1).sum(axis=1)
    return (got - k * (k + 1) / 2) / (k * (n - 1 - k))


def conservation(A, B, *, k=10):
    """Crow's co-expression conservation: the mean AUROC of both directions."""
    return .5 * (partner_auroc(A, B, k=k) + partner_auroc(B, A, k=k))


def binned_residuals(X, x, *, n_bins=20):
    """Each column minus the mean of the samples in its position bin (equal-count bins of x)."""
    X = np.asarray(X, float)
    order = np.argsort(np.argsort(np.asarray(x, float)))
    bins = np.minimum(order * n_bins // len(order), n_bins - 1)
    out = X.copy()
    for b in range(n_bins):
        m = bins == b
        out[m] -= X[m].mean(axis=0)
    return out
