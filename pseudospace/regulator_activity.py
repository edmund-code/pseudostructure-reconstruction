"""Regulator and signalling activity from target genes (notebook 95; protocol in docs/results/pt-regulator-activity.md).

* ``collectri_network``: OmniPath's CollecTRI rows -> signed TF -> target edges in mouse symbols. NF-kB dimers collapse to
  ``NFKB`` and AP-1 dimers to ``AP1`` (decoupleR's convention); the sign is -1 for an inhibition-only edge and +1
  otherwise (CollecTRI's default is activation); human symbols map through a one-to-one ortholog table.
* ``progeny_network``: OmniPath's PROGENy annotation rows -> the ``top`` genes per pathway by p value, with their weights.
* ``weight_matrix``: genes x sources sparse weights over a gene universe, sources with >= ``min_targets`` targets.
* ``ulm``: decoupleR's univariate linear model, vectorised: for every source and statistic vector, the t value of the
  slope of y ~ 1 + w over all genes in the universe (non-targets have w = 0).
* ``ulm_null``: the same t over ``n`` permutations of each statistic vector within gene strata (expression deciles),
  as z against the permutation mean and SD. Genes are treated as independent, so z is a ranking device; specificity
  comes from relabeled contrasts and replication.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse

NFKB = {'NFKB1', 'NFKB2', 'RELA', 'RELB', 'REL'}
AP1 = {'JUN', 'JUNB', 'JUND', 'FOS', 'FOSB', 'FOSL1', 'FOSL2'}


def _collapse(source):
    parts = set(str(source).replace('COMPLEX:', '').split('_'))
    if len(parts) > 1 and parts <= NFKB:
        return 'NFKB'
    if len(parts) > 1 and parts <= AP1:
        return 'AP1'
    return source if len(parts) == 1 else None


def collectri_network(rows, human_to_mouse):
    """Signed edges (source, target, weight) from CollecTRI rows; targets mapped to mouse, sources kept human-named.

    ``rows``: DataFrame with source_genesymbol, target_genesymbol, is_stimulation, is_inhibition.
    ``human_to_mouse``: mapping of one-to-one orthologs. Complexes other than NF-kB and AP-1 dimers are dropped; an
    edge listed twice (e.g. two NF-kB dimers on one target) keeps one row, with -1 only if every copy is -1.
    """
    df = pd.DataFrame({'source': rows.source_genesymbol.map(_collapse), 'target': rows.target_genesymbol.map(human_to_mouse),
                       'weight': np.where(rows.is_inhibition.astype(bool) & ~rows.is_stimulation.astype(bool), -1., 1.)})
    df = df.dropna(subset=['source', 'target'])
    return df.groupby(['source', 'target'], as_index=False).weight.max()


def progeny_network(rows, human_to_mouse, *, top=500):
    """PROGENy weights (pathway, target, weight): the ``top`` genes per pathway by p value, mapped to mouse."""
    wide = rows.pivot_table(index=['record_id', 'genesymbol'], columns='label', values='value', aggfunc='first').reset_index()
    wide = wide.dropna(subset=['pathway', 'weight', 'p_value'])
    wide['weight'] = wide.weight.astype(float)
    wide['p_value'] = wide.p_value.astype(float)
    wide = wide.sort_values(['pathway', 'p_value', 'genesymbol']).groupby('pathway').head(top)
    out = pd.DataFrame({'source': wide.pathway, 'target': wide.genesymbol.map(human_to_mouse), 'weight': wide.weight})
    out = out.dropna(subset=['target'])
    return out.groupby(['source', 'target'], as_index=False).weight.mean()


def weight_matrix(network, genes, *, min_targets=10):
    """Sparse genes x sources weights over ``genes``; sources with >= ``min_targets`` targets in ``genes``."""
    index = pd.Index(genes)
    net = network[network.target.isin(index)]
    sizes = net.groupby('source').size()
    sources = sorted(sizes.index[sizes >= min_targets])
    net = net[net.source.isin(sources)]
    col = pd.Index(sources)
    W = sparse.csr_matrix((net.weight.to_numpy(float), (index.get_indexer(net.target), col.get_indexer(net.source))),
                          shape=(len(index), len(col)))
    return W, sources


def ulm(W, Y):
    """t of the slope of y ~ 1 + w for every source (columns of W) and vector (columns of Y): sources x vectors."""
    W = sparse.csc_matrix(W, dtype=float)
    Y = np.asarray(Y, float)
    if Y.ndim == 1:
        Y = Y[:, None]
    if not np.isfinite(Y).all():
        raise ValueError('Statistic vectors must be finite.')
    n = Y.shape[0]
    sw = np.asarray(W.sum(axis=0)).ravel()
    sww = np.asarray(W.multiply(W).sum(axis=0)).ravel()
    sy, syy = Y.sum(axis=0), (Y ** 2).sum(axis=0)
    cov = np.asarray(W.T @ Y) - np.outer(sw, sy) / n
    var_w = sww - sw ** 2 / n
    var_y = syy - sy ** 2 / n
    r = cov / np.sqrt(np.outer(var_w, var_y))
    r = np.clip(r, -1 + 1e-12, 1 - 1e-12)
    return r * np.sqrt((n - 2) / (1 - r ** 2))


def stratified_permutations(strata, n, seed):
    """``n`` permutations (rows) of range(len(strata)) that move genes only within their stratum."""
    strata = np.asarray(strata)
    rng = np.random.default_rng(seed)
    out = np.tile(np.arange(len(strata)), (n, 1))
    for s in np.unique(strata):
        members = np.flatnonzero(strata == s)
        keys = rng.random((n, len(members)))
        out[:, members] = members[np.argsort(keys, axis=1)]
    return out


def ulm_null(W, y, permutations):
    """Observed ulm t and its z against the t over permuted copies of ``y`` (rows of ``permutations``)."""
    y = np.asarray(y, float)
    observed = ulm(W, y)[:, 0]
    null = ulm(W, y[np.asarray(permutations)].T)
    mu, sd = null.mean(axis=1), null.std(axis=1, ddof=1)
    return observed, (observed - mu) / np.where(sd > 0, sd, np.nan)
