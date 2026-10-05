"""Many-draw specimen nulls for PT pathway screens at segment resolution (notebook 45).

Notebook 37 judged pathway strategies by two balanced relabelings of four specimens. Public
snRNA atlases with many donors give many independent 2 + 2 designs instead. Their units (nuclei)
carry segment labels but no coordinate, so the screens are run at segment resolution:

- ``T_level``: a group offset (average level), notebook 37's strategy 4;
- ``T_steps``: group x segment beyond the offset (position, steps), notebook 37's strategy 6.

Every design here is constant within specimen x segment cells, so the weighted partial F of
notebook 37 follows exactly from cell sizes, cell means and within-cell sums of squares. One draw
then costs a regression on a few dozen rows instead of tens of thousands of units.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata

from .split_plot import log_cpm, moderated_f, voom_weights

SEGMENT_COMPARISONS = {'T_level': ('base', 'level'), 'T_steps': ('level', 'steps')}


def cell_summaries(expression, donor, segment, counts=None, library=None):
    """Sizes, means and within-cell sums of squares per donor x segment cell.

    ``expression`` (units x genes, sparse or dense) is the normalized matrix the models use;
    ``donor``/``segment`` label every unit. With ``counts`` and ``library`` the raw-count sums per
    cell (for pseudobulk) are kept as well. Returns a dict of arrays indexed [donor, segment, gene].
    """
    donor, segment = np.asarray(donor, str), np.asarray(segment, str)
    if expression.shape[0] != len(donor) or donor.shape != segment.shape:
        raise ValueError('Units, donors and segments must align.')
    donors, segments = np.unique(donor), np.unique(segment)
    d_code, s_code = np.searchsorted(donors, donor), np.searchsorted(segments, segment)
    cell = d_code * len(segments) + s_code
    n_cells = len(donors) * len(segments)
    indicator = sparse.csr_matrix((np.ones(len(cell)), (cell, np.arange(len(cell)))),
                                  shape=(n_cells, len(cell)))
    x = sparse.csr_matrix(expression, dtype=float)
    n = np.asarray(indicator.sum(axis=1)).ravel()
    sums = np.asarray((indicator @ x).todense())
    squares = np.asarray((indicator @ x.multiply(x)).todense())
    with np.errstate(invalid='ignore', divide='ignore'):
        mean = sums / n[:, None]
    within = np.maximum(squares - n[:, None] * np.nan_to_num(mean) ** 2, 0.)
    shape = (len(donors), len(segments), x.shape[1])
    out = {'donors': donors, 'segments': segments, 'n': n.reshape(shape[:2]),
           'mean': mean.reshape(shape), 'within': within.reshape(shape)}
    if counts is not None:
        raw = sparse.csr_matrix(counts, dtype=float)
        if raw.shape[0] != len(donor) or library is None or len(library) != len(donor):
            raise ValueError('Counts and library sizes must align with the units.')
        out['counts'] = np.asarray((indicator @ raw).todense()).reshape(shape[0], shape[1], raw.shape[1])
        out['library'] = (indicator @ np.asarray(library, float)).reshape(shape[:2])
    return out


def _orthobasis(design, root):
    u, singular, _ = np.linalg.svd(root[:, None] * design, full_matrices=False)
    tol = np.finfo(float).eps * max(design.shape) * (singular[0] if len(singular) else 0.)
    return u[:, singular > tol]


def segment_designs(donor_of_row, segment_of_row, group, nuisance=None):
    """Notebook 37's nested designs with segment intercepts as the shared shape.

    Rows are specimen x segment cells. ``group`` and ``nuisance`` map each specimen to 0/1.
    Specimen intercepts enter as within-group sum-to-zero contrasts; the nuisance factor's segment
    shape beyond its offset enters every model (relabelings), as in ``contrast_designs``.
    """
    d, s = np.asarray(donor_of_row, str), np.asarray(segment_of_row, str)
    labels = np.unique(s)
    shape = np.column_stack([(s == k).astype(float) for k in labels])
    g = np.array([group[x] for x in d], float)
    columns = [shape]
    for level in (0., 1.):
        names = sorted({x for x in d if group[x] == level})
        if not names:
            raise ValueError('Both groups need at least one specimen.')
        columns += [((d == name).astype(float) - (d == names[-1]).astype(float))[:, None] for name in names[:-1]]
    if nuisance is not None:
        q = np.array([nuisance[x] for x in d], float)
        columns.append(q[:, None] * shape[:, 1:])
    base = np.column_stack(columns)
    level = np.column_stack([base, g])
    return {'base': base, 'level': level, 'steps': np.column_stack([level, g[:, None] * shape])}


def summary_partial_f(summary, donors, group, nuisance=None, comparisons=None):
    """Exact unit-level weighted partial F for segment designs, from cell summaries.

    ``donors``: the specimens in this design; ``group``/``nuisance``: dicts specimen -> 0/1. Unit
    weights are notebook 37's (each group half the weight, each specimen an equal share). Returns
    {statistic: F per gene} plus 'df'.
    """
    pairs = comparisons or SEGMENT_COMPARISONS
    index = pd.Index(summary['donors'])
    rows_d, rows_s, means, within, n_rows = [], [], [], [], []
    for name in donors:
        i = index.get_loc(name)
        for j, seg in enumerate(summary['segments']):
            if summary['n'][i, j] > 0:
                rows_d.append(name), rows_s.append(seg)
                means.append(summary['mean'][i, j]), within.append(summary['within'][i, j])
                n_rows.append(summary['n'][i, j])
    rows_d, n_rows = np.array(rows_d), np.array(n_rows, float)
    total = n_rows.sum()
    unit_weight = np.empty(len(rows_d))
    for level in (0., 1.):
        names = [x for x in donors if group[x] == level]
        for name in names:
            mask = rows_d == name
            unit_weight[mask] = total / (2 * len(names) * n_rows[mask].sum())
    designs = segment_designs(rows_d, rows_s, group, nuisance)
    root = np.sqrt(unit_weight * n_rows)
    y = root[:, None] * np.vstack(means)
    common = (unit_weight[:, None] * np.vstack(within)).sum(axis=0)
    sse, rank = {}, {}
    for name in {x for pair in pairs.values() for x in pair}:
        q = _orthobasis(designs[name], root)
        resid = y - q @ (q.T @ y)
        sse[name], rank[name] = common + np.einsum('ij,ij->j', resid, resid), q.shape[1]
    out, df = {}, {}
    for name, (reduced, full) in pairs.items():
        df[name] = (rank[full] - rank[reduced], total - rank[full])
        if df[name][0] <= 0 or df[name][1] <= 2:
            raise ValueError(f'{name}: no added columns or no residual degrees of freedom.')
        out[name] = (np.maximum(sse[reduced] - sse[full], 0.) / df[name][0]
                     / np.maximum(sse[full] / df[name][1], 1e-12))
    out['df'] = df
    return out


def membership(gene_sets, genes):
    """Sparse pathways x genes indicator in the order of ``genes``."""
    index = pd.Index(genes)
    rows, cols = [], []
    for i, members in enumerate(gene_sets.values()):
        idx = index.get_indexer(sorted(set(members)))
        if (idx < 0).any():
            raise ValueError('Every member must be in the gene universe.')
        rows.extend([i] * len(idx)), cols.extend(idx)
    return sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(gene_sets), len(index)))


def matched_auc_null(statistics, gene_sets, strata):
    """Rank-AUC with the exact mean and SD of notebook 12's covariate-matched random-set null.

    ``matched_pathway_tests`` draws equally sized random sets within each stratum (without
    replacement, from the whole universe) and uses their AUC mean and SD in a normal approximation
    (``normal_matched_q``). Those moments are exact here: a stratum contributing c of its N genes
    adds c * mean(ranks) to the expected rank sum and c * var(ranks) * (N - c) / (N - 1) to its
    variance. ``statistics``: DataFrame genes x statistics. Returns long rows with auc, effect,
    null_auc_mean, null_auc_sd and z_matched.
    """
    values = statistics.to_numpy(float)
    strata = np.asarray(strata)
    if not np.isfinite(values).all() or len(strata) != len(values):
        raise ValueError('Statistics must be finite and aligned with strata.')
    m = membership(gene_sets, statistics.index)
    labels = np.unique(strata)
    onehot = sparse.csr_matrix((np.ones(len(strata)), (np.arange(len(strata)), np.searchsorted(labels, strata))),
                               shape=(len(strata), len(labels)))
    count = np.asarray((m @ onehot).todense())                       # pathways x strata
    size = np.asarray(onehot.sum(axis=0)).ravel()
    n = count.sum(axis=1)
    total = len(values)
    scale = n * (total - n)
    frames = []
    for j, name in enumerate(statistics.columns):
        ranks = rankdata(values[:, j])
        mean_k = np.asarray(onehot.T @ ranks).ravel() / size
        var_k = np.asarray(onehot.T @ ranks ** 2).ravel() / size - mean_k ** 2
        finite = np.where(size > 1, (size[None, :] - count) / np.maximum(size[None, :] - 1, 1), 0.)
        null_var = (count * var_k[None, :] * finite).sum(axis=1)
        observed = np.asarray(m @ ranks).ravel()
        offset = n * (n + 1) / 2
        auc = (observed - offset) / scale
        null_mean = (count @ mean_k - offset) / scale
        null_sd = np.sqrt(np.maximum(null_var, 0.)) / scale
        frames.append(pd.DataFrame({'pathway_id': list(gene_sets), 'statistic': name, 'n_genes': n.astype(int),
                                    'auc': auc, 'effect': auc - .5, 'null_auc_mean': null_mean,
                                    'null_auc_sd': null_sd,
                                    'z_matched': (auc - null_mean) / np.maximum(null_sd, 1e-12)}))
    return pd.concat(frames, ignore_index=True)


def pseudobulk_t(counts, library, design, coefficient):
    """voom-style moderated t for one coefficient of a sample-level design (limma port).

    ``counts``: samples x genes raw sums; ``library``: samples. Returns (t, p) per gene.
    """
    y = log_cpm(counts, library)
    weights = voom_weights(y, library, design)
    fit = moderated_f(y, weights, design, [coefficient])
    t = np.sign(fit['beta'][:, coefficient]) * np.sqrt(fit['F'])
    return t, fit['p']


def donor_deviations(summary, donors, groups=None):
    """Specimen offsets and specimen x segment deviations of cell means (two-way decomposition).

    With ``groups`` (specimen -> label, for example species), each group's mean segment profile is
    removed first, so only deviations between specimens of the same group remain. Returns
    (offsets: donors x genes, interactions: (donors * segments) x genes), centred over specimens,
    with equal segment weights. Inputs to the pathway-coherence check.
    """
    index = pd.Index(summary['donors'])
    mean = summary['mean'][index.get_indexer(donors)].copy()
    if not np.isfinite(mean).all():
        raise ValueError('Every specimen needs every segment.')
    if groups is not None:
        labels = np.array([groups[d] for d in donors])
        for label in np.unique(labels):
            mean[labels == label] -= mean[labels == label].mean(axis=0, keepdims=True)
    grand = mean.mean(axis=(0, 1), keepdims=True)
    donor_effect = mean.mean(axis=1, keepdims=True) - grand
    segment_effect = mean.mean(axis=0, keepdims=True) - grand
    interaction = mean - grand - donor_effect - segment_effect
    return donor_effect[:, 0, :], interaction.reshape(-1, mean.shape[2])
