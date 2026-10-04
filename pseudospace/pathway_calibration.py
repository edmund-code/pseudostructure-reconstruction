"""Negative-control calibration of PT pathway screens by balanced specimen relabeling.

With two specimens per species, a species contrast and a specimen contrast use
the same structure-level machinery. Relabeling the four specimens into two
mouse+human groups keeps the weights, gene universe and pathway test but removes
the species difference (the random-partition control of Lamian, Hou et al. 2023).
A screen that still reports many pathways is measuring specimen heterogeneity or
structure-level pseudoreplication, not species biology.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import f as f_dist, norm
from statsmodels.stats.multitest import multipletests

from .levelshape import build_ls_designs

COMPARISONS = {
    'T_level': ('base', 'level'),                      # constant species offset
    'T_spatial': ('level', 'full'),                    # smooth species-by-position change
    'T_total': ('base', 'full'),                       # either
    'T_discrete_total': ('base', 'discrete'),          # species-by-segment steps
    'T_discrete_spatial': ('level', 'discrete_shape'),  # steps beyond a constant offset
    'T_position_given_segments': ('discrete', 'union'),  # smooth position beyond steps
}


def balanced_partitions(specimen, species):
    """The species split of four specimens plus its two balanced relabelings.

    Returns {label: (group-1 specimens)}; 'species' puts the second species in
    group 1, and each 'swap:...' groups one specimen of each species.
    """
    specimen, species = np.asarray(specimen, str), np.asarray(species, str)
    owner = pd.Series(species).groupby(specimen).unique()
    if len(owner) != 4 or owner.map(len).ne(1).any():
        raise ValueError('Expected four specimens, each from exactly one species.')
    owner = owner.map(lambda v: v[0])
    kinds = sorted(owner.unique())
    if len(kinds) != 2 or owner.value_counts().ne(2).any():
        raise ValueError('Expected two specimens from each of two species.')
    first, second = (sorted(owner.index[owner.eq(k)]) for k in kinds)
    partitions = {'species': tuple(second)}
    for partner in second:
        partitions['swap:' + first[0] + '+' + partner] = (first[0], partner)
    return partitions


def contrast_designs(position, group, specimen, knots, *, segments=None, nuisance_shape=None):
    """Nested weighted designs for one two-group comparison along a shared coordinate.

    Specimen intercepts enter as within-group sum-to-zero contrasts, as in
    notebook 12, so the group level stays identifiable. ``nuisance_shape`` is an
    optional second 0/1 partition (species, when ``group`` is a relabeling) whose
    smooth position-dependent difference enters every model. ``segments`` adds
    group-by-segment step models. Each group gets half the weight and each of its
    specimens an equal share.
    """
    s, g = np.asarray(position, float), np.asarray(group, float)
    samples = np.asarray(specimen, str)
    if not len(s) == len(g) == len(samples) or set(np.unique(g)) != {0., 1.}:
        raise ValueError('Need aligned positions, a 0/1 group and specimens.')
    n = len(s)
    weights, nuisance = np.zeros(n), []
    for level in (0., 1.):
        names = np.unique(samples[g == level])
        for name in names:
            mask = samples == name
            if not np.all(g[mask] == level):
                raise ValueError('A specimen cannot belong to both groups.')
            weights[mask] = n / (2 * len(names) * mask.sum())
        nuisance.extend((samples == name).astype(float) - (samples == names[-1]).astype(float)
                        for name in names[:-1])
    x0, x1, x2 = build_ls_designs(s, g, knots)[:3]
    if nuisance_shape is not None:
        q = np.asarray(nuisance_shape, float)
        if q.shape != s.shape:
            raise ValueError('nuisance_shape must align with positions.')
        nuisance.extend((q[:, None] * x0[:, 1:]).T)
    extra = np.column_stack(nuisance) if nuisance else np.empty((n, 0))
    designs = {name: np.column_stack([x, extra]) for name, x in
               (('base', x0), ('level', x1), ('full', x2))}
    if segments is not None:
        labels = np.asarray(segments, str)
        if labels.shape != s.shape:
            raise ValueError('segments must align with positions.')
        steps = np.column_stack([g * (labels == name) for name in np.unique(labels)])
        designs['discrete'] = np.column_stack([designs['base'], steps])
        designs['discrete_shape'] = np.column_stack([designs['level'], steps])
        designs['union'] = np.column_stack([designs['discrete'], designs['full']])
    return designs, weights


def _orthobasis(design, root_weights):
    u, singular, _ = np.linalg.svd(root_weights[:, None] * design, full_matrices=False)
    tol = np.finfo(float).eps * max(design.shape) * (singular[0] if len(singular) else 0.)
    return u[:, singular > tol]


def nested_partial_f(expression, designs, weights, comparisons=None, *, block_size=512):
    """Weighted partial-F scores for nested designs (redundant columns allowed).

    Each comparison is (reduced, full) design names; comparisons whose designs are
    absent are skipped. Returns {name: scores} plus 'df' = {name: (added, residual)}.
    """
    pairs = {k: v for k, v in (comparisons or COMPARISONS).items()
             if v[0] in designs and v[1] in designs}
    w = np.asarray(weights, float)
    n = expression.shape[0]
    if not pairs or w.shape != (n,) or not np.isfinite(w).all() or (w <= 0).any():
        raise ValueError('Need comparisons and positive finite weights aligned to rows.')
    root = np.sqrt(w)
    used = sorted({name for pair in pairs.values() for name in pair})
    bases = {name: _orthobasis(np.asarray(designs[name], float), root) for name in used}
    df = {}
    for name, (reduced, full) in pairs.items():
        qr, qf = bases[reduced], bases[full]
        if np.linalg.norm(qr - qf @ (qf.T @ qr)) > 1e-8 * max(1., np.sqrt(qr.shape[1])):
            raise ValueError(f'{name}: {reduced} is not nested in {full}.')
        df[name] = (qf.shape[1] - qr.shape[1], n - qf.shape[1])
        if df[name][0] <= 0 or df[name][1] <= 2:
            raise ValueError(f'{name}: no added columns or no residual degrees of freedom.')
    sse = {name: np.empty(expression.shape[1]) for name in used}
    for start in range(0, expression.shape[1], block_size):
        block = expression[:, start:start + block_size]
        block = block.toarray() if sparse.issparse(block) else np.asarray(block, float)
        if not np.isfinite(block).all():
            raise ValueError('Expression must be finite.')
        yw = root[:, None] * block
        for name in used:
            resid = yw - bases[name] @ (bases[name].T @ yw)
            sse[name][start:start + block.shape[1]] = np.einsum('ij,ij->j', resid, resid)
    scores = {name: np.maximum(sse[r] - sse[f], 0.) / df[name][0]
              / np.maximum(sse[f] / df[name][1], 1e-12) for name, (r, f) in pairs.items()}
    scores['df'] = df
    return scores


def normal_matched_q(tests):
    """Add a normal-approximation p/q from the matched null mean and SD.

    With B draws an empirical p cannot fall below 1/(B+1), so BH over m pathways
    calls nothing unless about m/(alpha*(B+1)) of them sit at that floor: a few
    strong pathways read as zero. The z-score keeps the matched null, not its floor.
    """
    out = tests.copy()
    out['z_matched'] = (out.auc - out.null_auc_mean) / out.null_auc_sd.clip(lower=1e-12)
    out['p_normal'] = norm.sf(out.z_matched)
    keys = [c for c in ('partition', 'statistic') if c in out]
    out['q_normal'] = out.groupby(keys).p_normal.transform(
        lambda p: multipletests(p, method='fdr_bh')[1])
    return out


def pathway_scores(expression, gene_sets, genes, weights):
    """Per-structure mean of weighted-standardized member expression (a module score)."""
    index = pd.Index(genes)
    x = expression.toarray() if sparse.issparse(expression) else np.asarray(expression, float)
    w = np.asarray(weights, float) / np.sum(weights)
    mean = w @ x
    z = (x - mean) / np.sqrt(np.maximum(w @ (x - mean) ** 2, 1e-12))
    members = {p: index.get_indexer(v) for p, v in gene_sets.items()}
    if any((idx < 0).any() for idx in members.values()):
        raise ValueError('Every pathway member must be in the gene universe.')
    return pd.DataFrame({p: z[:, idx].mean(axis=1) for p, idx in members.items()})


def parametric_calls(scores, df, alpha=.05):
    """Structure-level F-test p/q for pathway scores: treats structures as replicates."""
    frames = []
    for name, values in scores.items():
        p = f_dist.sf(np.asarray(values, float), *df[name])
        frames.append(pd.DataFrame({'statistic': name, 'pathway_id': list(values.index),
                                    'F': values.to_numpy(), 'p': p,
                                    'q': multipletests(p, method='fdr_bh')[1]}))
    out = pd.concat(frames, ignore_index=True)
    out['hit'] = out.q.le(alpha)
    return out


def negative_control_summary(hits, *, true_label='species'):
    """Discoveries per strategy for the true and relabeled partitions.

    ``hits``: one row per (strategy, partition, pathway_id) with boolean 'hit'.
    The negative-control discovery ratio (NCDR) is the mean relabeled count over
    the true count: near 0 is specific; near or above 1 means the screen reports
    as much under relabeling as for the real contrast.
    """
    counts = hits.groupby(['strategy', 'partition']).hit.sum().unstack('partition').fillna(0)
    nulls = [c for c in counts if c != true_label]
    if true_label not in counts or not nulls:
        raise ValueError('Need the true partition and at least one relabeling.')
    summary = pd.DataFrame({'true_hits': counts[true_label], 'null_mean': counts[nulls].mean(axis=1),
                            'null_max': counts[nulls].max(axis=1)})
    summary['ncdr'] = summary.null_mean / summary.true_hits.where(summary.true_hits > 0)
    flagged = (hits[hits.partition.isin(nulls) & hits.hit]
               .groupby('strategy').pathway_id.apply(set))
    true_sets = hits[hits.partition.eq(true_label) & hits.hit].groupby('strategy').pathway_id.apply(set)
    summary['true_hits_also_in_null'] = [len(true_sets.get(s, set()) & flagged.get(s, set()))
                                         for s in summary.index]
    summary['true_only'] = summary.true_hits - summary.true_hits_also_in_null
    return summary.join(counts[nulls].add_prefix('hits_'))
