"""Overlap-preserving nulls for collections of gene sets (notebook 48).

A list of called pathways shares genes, so per-pathway scores are not independent and a
Mann–Whitney test across pathways treats one shared program as many observations. Here one random
relabeling of gene identities, drawn within matching strata, is applied to every set at once. Each
permuted set keeps its size and strata composition, and every pairwise overlap between sets is kept
exactly. The pathway score is notebook 37's joint replication z: the covariate-matched rank-AUC z
(exact null moments, which a within-strata relabeling leaves unchanged) divided by the square root of
a CAMERA variance inflation from residual correlations of the set's current members.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata

from .specimen_null import matched_auc_null, membership


def stratified_permutation(strata, rng):
    """A random gene relabeling that only exchanges genes within the same stratum."""
    strata = np.asarray(strata)
    perm = np.arange(len(strata))
    for label in np.unique(strata):
        idx = np.flatnonzero(strata == label)
        perm[idx] = rng.permutation(idx)
    return perm


class JointReplicationZ:
    """Joint z (matched rank-AUC z ÷ √VIF) for many gene sets, under any within-strata relabeling.

    ``statistic``: Series over genes; ``strata``: matching strata aligned to it; ``gene_sets``: dict
    of member lists; ``residuals``: replicates x genes residual matrix whose within-set correlations
    give the VIF (one block, as notebook 37's donor-level residuals).
    """

    def __init__(self, statistic, strata, gene_sets, residuals):
        statistic = pd.Series(statistic, dtype=float)
        residuals = np.asarray(residuals, float)
        if residuals.shape[1] != len(statistic) or residuals.shape[0] < 4:
            raise ValueError('Residuals must be replicates x genes, with at least four replicates.')
        self.ids = list(gene_sets)
        moments = matched_auc_null(statistic.to_frame('a'), gene_sets, strata).set_index('pathway_id').loc[self.ids]
        self.null_mean = moments.null_auc_mean.to_numpy()
        self.null_sd = np.maximum(moments.null_auc_sd.to_numpy(), 1e-12)
        self.m = membership(gene_sets, statistic.index)
        self.size = np.asarray(self.m.sum(axis=1)).ravel()
        self.total = len(statistic)
        self.ranks = rankdata(statistic.to_numpy())
        ties = np.unique(statistic.to_numpy(), return_counts=True)[1].astype(float)
        self.tie_factor = 1 - np.sum(ties ** 3 - ties) / (self.total ** 3 - self.total)
        centred = residuals - residuals.mean(axis=0)
        lengths = np.linalg.norm(centred, axis=0)
        self.variable = lengths > 1e-12 * np.sqrt(len(centred))
        self.columns = np.where(self.variable, centred / np.where(self.variable, lengths, 1.), 0.)

    def z(self, perm=None):
        """Joint z of every set after mapping each member g to perm[g] (identity when None)."""
        ranks = self.ranks if perm is None else self.ranks[perm]
        columns = self.columns if perm is None else self.columns[:, perm]
        variable = self.variable if perm is None else self.variable[perm]
        n, rest = self.size, self.total - self.size
        auc = (self.m @ ranks - n * (n + 1) / 2) / (n * rest)
        z_matched = (auc - self.null_mean) / self.null_sd
        summed = np.asarray(self.m @ columns.T)
        k = np.asarray(self.m @ variable.astype(float)).ravel()
        with np.errstate(invalid='ignore', divide='ignore'):
            rho = np.clip(((summed ** 2).sum(axis=1) - k) / (k * (k - 1)), -1., 1.)
        rho = np.where(k >= 2, np.maximum(rho, 0.), np.nan)
        variance = ((np.arcsin(1.) * n * rest + np.arcsin(.5) * n * rest * (rest - 1)
                     + np.arcsin(rho / 2) * n * (n - 1) * rest * (rest - 1)
                     + np.arcsin((1 + rho) / 2) * n * (n - 1) * rest) / (2 * np.pi) * self.tie_factor)
        independent = n * rest * (self.total + 1) / 12 * self.tie_factor
        return z_matched / np.sqrt(variance / independent)


def collection_null(scorer, collections, strata, n_perm, seed):
    """Observed and permuted mean joint z of each collection (dict name -> list of set ids).

    Every permutation relabels genes within ``strata`` and rescores all sets at once, so overlaps
    within and between collections are kept. Returns (observed Series, permuted DataFrame).
    """
    index = {p: i for i, p in enumerate(scorer.ids)}
    rows = {name: np.array([index[p] for p in ids]) for name, ids in collections.items()}
    observed_z = scorer.z()
    observed = pd.Series({name: float(np.mean(observed_z[r])) for name, r in rows.items()})
    rng = np.random.default_rng(seed)
    draws = np.empty((n_perm, len(rows)))
    for b in range(n_perm):
        z = scorer.z(stratified_permutation(strata, rng))
        draws[b] = [np.mean(z[r]) for r in rows.values()]
    return observed, pd.DataFrame(draws, columns=list(rows))


def upper_p(observed, null):
    """One-sided permutation p with the observed value counted: (1 + #{null >= obs}) / (B + 1)."""
    null = np.asarray(null, float)
    return (1 + np.count_nonzero(null >= observed - 1e-12)) / (len(null) + 1)


def label_permutation_p(first, second, n_perm=9999, seed=0, exact_limit=50000):
    """One-sided p that mean(first) - mean(second) is this large when unit labels are exchangeable.

    Exact enumeration when the number of splits is at most ``exact_limit``; otherwise random splits.
    """
    from itertools import combinations
    from math import comb

    values = np.concatenate([np.asarray(first, float), np.asarray(second, float)])
    k, total = len(first), len(values)
    observed = values[:k].mean() - values[k:].mean()
    if comb(total, k) <= exact_limit:
        diffs = []
        for chosen in combinations(range(total), k):
            mask = np.zeros(total, bool)
            mask[list(chosen)] = True
            diffs.append(values[mask].mean() - values[~mask].mean())
        diffs = np.array(diffs)
        return float(np.mean(diffs >= observed - 1e-12)), observed
    rng = np.random.default_rng(seed)
    diffs = np.empty(n_perm)
    for b in range(n_perm):
        mask = np.zeros(total, bool)
        mask[rng.choice(total, k, replace=False)] = True
        diffs[b] = values[mask].mean() - values[~mask].mean()
    return upper_p(observed, diffs), observed
