"""Gene-first PT pathway remodeling for notebook 12; no coordinate estimation.

Statistics are conditional on the observed structures. Neither structural HC3
errors nor competitive gene-set nulls supply independent human donors.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import mannwhitneyu, rankdata
from statsmodels.stats.multitest import multipletests

from .levelshape import build_ls_designs
from .stats_gam import gam_internal_knots, make_gam_design


def fit_nested_trajectories(expression, position, human, specimen, grid, *, basis_df=6):
    """Fixed-basis Gaussian regression-spline GAMs with balanced specimen weights.

    The same within-species sum-to-zero specimen intercept contrasts enter all
    three models. This avoids confounding an unconstrained specimen dummy with
    species. No penalty is used: fixed low basis dimension regularizes the fit
    and preserves exact nesting. Return partial-F *ranking* statistics and HC3
    structural uncertainty for the full-model human-minus-mouse contrast.
    """
    y = sparse.csr_matrix(expression, dtype=float)
    s, c = np.asarray(position, float), np.asarray(human, float)
    samples, grid = np.asarray(specimen, str), np.asarray(grid, float)
    if y.shape[0] != len(s) or s.shape != c.shape or s.shape != samples.shape:
        raise ValueError('Expression rows, coordinate, species and specimens must align.')
    if (y.shape[1] == 0 or not np.isfinite(y.data).all() or
            not np.isfinite(s).all() or not np.isfinite(grid).all() or
            s.min() < 0 or s.max() > 1 or grid.size < 2 or
            grid.min() < s.min() or grid.max() > s.max() or np.any(np.diff(grid) <= 0)):
        raise ValueError('Need finite expression and increasing grid within observed [0, 1] support.')
    if set(np.unique(c)) != {0., 1.} or basis_df < 3:
        raise ValueError('Need both species (mouse=0, human=1) and basis_df >= 3.')
    n = len(s)
    weights = np.zeros(n)
    nuisance = []
    for group in (0., 1.):
        names = np.unique(samples[c == group])
        for name in names:
            mask = samples == name
            if not np.all(c[mask] == group):
                raise ValueError('A specimen cannot belong to both species.')
            weights[mask] = n / (2 * len(names) * mask.sum())
        # Contrasts predict zero at the equal-specimen species mean.
        nuisance.extend((samples == name).astype(float) - (samples == names[-1]).astype(float)
                        for name in names[:-1])
    knots = gam_internal_knots(s, basis_df=basis_df)
    x0, x1, x2, _, _, _, p_base, _ = build_ls_designs(s, c, knots)
    designs = [np.column_stack([x, *nuisance]) for x in (x0, x1, x2)]
    solvers = []
    for x in designs:
        if np.linalg.matrix_rank(x) != x.shape[1] or n <= x.shape[1] + 2:
            raise ValueError('Insufficient coordinate support for the requested spline/specimen design.')
        inv = np.linalg.inv(x.T @ (weights[:, None] * x))
        solvers.append(inv @ (x.T * weights))
    base = make_gam_design(grid, knots)
    mouse_design = np.column_stack([base, np.zeros((len(grid), p_base + len(nuisance)))])
    human_design = mouse_design.copy()
    human_design[:, p_base] = 1
    human_design[:, p_base + 1:2 * p_base] = base[:, 1:]
    contrast = human_design - mouse_design
    projection = contrast @ solvers[2]
    leverage = np.einsum('ij,ji->i', designs[2], solvers[2])
    if np.any(leverage >= .99):
        raise ValueError('Near-unit leverage: reduce spline dimension or improve support.')
    g = y.shape[1]
    result = {key: np.empty(g) for key in ('T_level', 'T_spatial', 'T_total')}
    result.update({key: np.empty((g, len(grid))) for key in ('mouse', 'human', 'delta', 'se', 'z')})
    # Ponytail: 256-gene blocks bound dense residual memory; tune only after profiling.
    for start in range(0, g, 256):
        block = slice(start, min(start + 256, g))
        values = y[:, block].toarray()
        betas = [solver @ values for solver in solvers]
        residuals = [values - x @ beta for x, beta in zip(designs, betas)]
        sse = [np.sum(weights[:, None] * residual ** 2, axis=0) for residual in residuals]
        for key, reduced, full in (('T_level', 0, 1), ('T_spatial', 1, 2), ('T_total', 0, 2)):
            added = designs[full].shape[1] - designs[reduced].shape[1]
            noise = np.maximum(sse[full] / (n - designs[full].shape[1]), 1e-12)
            result[key][block] = np.maximum(sse[reduced] - sse[full], 0) / added / noise
        delta = contrast @ betas[2]
        # HC3 sandwich: balance weights are not assumed to be inverse variances.
        se = np.sqrt(np.maximum(projection ** 2 @ (residuals[2] / (1 - leverage[:, None])) ** 2, 1e-12))
        result['mouse'][block] = (mouse_design @ betas[2]).T
        result['human'][block] = (human_design @ betas[2]).T
        result['delta'][block], result['se'][block], result['z'][block] = delta.T, se.T, (delta / se).T
    result['knots'], result['grid'] = knots, grid
    result['design_df'] = np.array([x.shape[1] for x in designs])
    return result


def gene_covariates(expression, position, human, specimen, genes, *, n_bins=10):
    """Equal-specimen mean/detection; fraction of occupied bins with any detection."""
    x = sparse.csr_matrix(expression)
    s, c, samples = np.asarray(position), np.asarray(human), np.asarray(specimen)
    names = np.unique(samples)
    means = np.array([np.asarray(x[samples == name].mean(axis=0)).ravel() for name in names])
    detections = np.array([np.asarray((x[samples == name] > 0).mean(axis=0)).ravel() for name in names])
    edges = np.linspace(s.min(), s.max(), n_bins + 1)
    bins = np.clip(np.searchsorted(edges, s, side='right') - 1, 0, n_bins - 1)
    coverage = []
    for name in names:
        present = [np.asarray((x[(samples == name) & (bins == b)] > 0).sum(axis=0)).ravel() > 0
                   for b in np.unique(bins[samples == name])]
        coverage.append(np.mean(present, axis=0))
    frame = pd.DataFrame({'gene': genes, 'mean_expression': means.mean(axis=0),
                          'detection': detections.mean(axis=0), 'coverage': np.mean(coverage, axis=0)})
    for label, group in (('mouse', 0), ('human', 1)):
        frame[f'detection_{label}'] = np.asarray((x[c == group] > 0).mean(axis=0)).ravel()
    return frame


def matching_strata(covariates, *, bins=3):
    """Coarse quantile matching, with ties kept together; no outcome used."""
    codes = []
    for name in ('mean_expression', 'detection', 'coverage'):
        values = covariates[name].to_numpy(float)
        if not np.isfinite(values).all():
            raise ValueError(f'Non-finite matching covariate: {name}')
        edges = np.unique(np.quantile(values, np.linspace(0, 1, bins + 1)[1:-1]))
        codes.append(np.searchsorted(edges, values, side='right'))
    return np.unique(np.column_stack(codes), axis=0, return_inverse=True)[1]


def rank_auc(scores, members):
    """Probability that a member exceeds a nonmember, counting a tie as one half."""
    values = np.asarray(scores, float)
    idx = np.asarray(members, int)
    n, total = len(idx), len(values)
    if not 0 < n < total or len(np.unique(idx)) != n or (idx < 0).any() or (idx >= total).any():
        raise ValueError('Need unique in-range members and a nonempty complement.')
    if not np.isfinite(values).all():
        raise ValueError('Ranking scores must be finite.')
    ranks = rankdata(values, axis=0, method='average')
    return (ranks[idx].sum(axis=0) - n * (n + 1) / 2) / (n * (total - n))


def matched_pathway_tests(statistics, gene_sets, strata, *, n_null=9999, seed=12):
    """One-sided rank-AUC tests against equally sized covariate-matched gene sets.

    Draw without replacement within every stratum, from the entire eligible
    universe (including original members). Reuse stratum/count draws across
    pathways and statistics. BH spans all supplied pathway/statistic pairs.
    Does NOT preserve within-pathway gene correlation or calibrate donor inference.
    """
    if n_null < 1 or not statistics.index.is_unique or not gene_sets:
        raise ValueError('Need permutations, unique genes and at least one gene set.')
    values, strata = statistics.to_numpy(float), np.asarray(strata)
    if not np.isfinite(values).all() or len(strata) != len(values):
        raise ValueError('Scores and matching strata must be finite and aligned.')
    ranks = rankdata(values, axis=0)
    rng = np.random.default_rng(seed)
    pools = {label: np.flatnonzero(strata == label) for label in np.unique(strata)}
    rank_sums, rows = {}, []
    for pathway, members in gene_sets.items():
        idx = statistics.index.get_indexer(sorted(set(members)))
        if (idx < 0).any() or not 0 < len(idx) < len(values):
            raise ValueError(f'{pathway}: invalid membership for this universe.')
        n, total = len(idx), len(values)
        null_sum = np.zeros((n_null, values.shape[1]))
        fixed_members = 0
        for label, count in zip(*np.unique(strata[idx], return_counts=True)):
            key = (label, count)
            if key not in rank_sums:
                pool = pools[label]
                if count == len(pool):
                    sums = np.broadcast_to(ranks[pool].sum(axis=0), null_sum.shape).copy()
                else:
                    draws = np.array([rng.choice(pool, count, replace=False) for _ in range(n_null)])
                    sums = ranks[draws].sum(axis=1)
                rank_sums[key] = sums
            null_sum += rank_sums[key]
            fixed_members += count if count == len(pools[label]) else 0
        observed = (ranks[idx].sum(axis=0) - n * (n + 1) / 2) / (n * (total - n))
        null = (null_sum - n * (n + 1) / 2) / (n * (total - n))
        outside = np.ones(total, bool)
        outside[idx] = False
        for j, statistic in enumerate(statistics.columns):
            hits = np.count_nonzero(null[:, j] >= observed[j] - 1e-12)
            p = (hits + 1) / (n_null + 1)
            rows.append({'pathway_id': pathway, 'statistic': statistic, 'n_genes': n,
                         'auc': observed[j], 'effect': observed[j] - .5,
                         'p_nominal': mannwhitneyu(values[idx, j], values[outside, j], alternative='greater').pvalue,
                         'p_empirical': p, 'mc_se': np.sqrt(p * (1 - p) / (n_null + 1)),
                         'null_auc_mean': null[:, j].mean(), 'null_auc_sd': null[:, j].std(),
                         'fixed_member_fraction': fixed_members / n, 'n_null': n_null})
    frame = pd.DataFrame(rows)
    frame['q_empirical'] = multipletests(frame.p_empirical, method='fdr_bh')[1]
    return frame


def local_pathway_curves(z, genes, gene_sets, grid):
    """Local divergence and signed relative-rank direction, using every eligible gene."""
    index = pd.Index(genes)
    z = np.asarray(z, float)
    if z.shape != (len(genes), len(grid)) or not np.isfinite(z).all():
        raise ValueError('Local Z must be finite and aligned to genes and grid.')
    magnitude_ranks, signed_ranks = rankdata(np.abs(z), axis=0), rankdata(z, axis=0)
    rows = []
    for pathway, members in gene_sets.items():
        idx = index.get_indexer(members)
        n, total = len(idx), len(index)
        if (idx < 0).any() or len(np.unique(idx)) != n or not 0 < n < total:
            raise ValueError(f'{pathway}: invalid local pathway membership.')
        offset, denominator = n * (n + 1) / 2, n * (total - n)
        divergence = (magnitude_ranks[idx].sum(axis=0) - offset) / denominator - .5
        direction = 2 * ((signed_ranks[idx].sum(axis=0) - offset) / denominator - .5)
        # An absolute signed summary disambiguates relative rank from absolute direction.
        median_z = np.median(z[idx], axis=0)
        rows.append(pd.DataFrame({'pathway_id': pathway, 'position': grid,
                                  'divergence': divergence, 'direction': direction,
                                  'median_member_z': median_z,
                                  'human_high_fraction': np.mean(z[idx] > 0, axis=0)}))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(
        columns=['pathway_id', 'position', 'divergence', 'direction', 'median_member_z', 'human_high_fraction'])


def spatial_descriptions(curves):
    """Post-discovery descriptions; half-height widths and sign deadband are not tests."""
    rows = []
    for pathway, frame in curves.groupby('pathway_id', sort=False):
        frame = frame.sort_values('position')
        s, d, direction = (frame[key].to_numpy() for key in ('position', 'divergence', 'direction'))
        peak = int(np.argmax(d))
        affected = (d > 0) & (d >= max(d[peak], 0) / 2)
        signs = np.sign(direction[np.abs(direction) >= .1])
        positive = np.maximum(d, 0)
        thirds = np.minimum(((s - s.min()) / (s.max() - s.min()) * 3).astype(int), 2)
        mass = np.array([positive[thirds == b].sum() for b in range(3)])
        shares = mass / mass.sum() if mass.sum() else np.full(3, np.nan)
        rows.append({'pathway_id': pathway, 'peak_position': s[peak], 'peak_divergence': d[peak],
                     'direction_at_peak': direction[peak],
                     'median_member_z_at_peak': frame.median_member_z.iloc[peak],
                     'affected_grid_fraction': affected.mean(),
                     'affected_width': affected.mean() * (s.max() - s.min()),
                     'direction_sign_changes': int(np.count_nonzero(np.diff(signs))) if len(signs) else 0,
                     'early_share': shares[0], 'mid_share': shares[1], 'late_share': shares[2]})
    return pd.DataFrame(rows).reindex(columns=['pathway_id', 'peak_position', 'peak_divergence',
        'direction_at_peak', 'median_member_z_at_peak', 'affected_grid_fraction', 'affected_width',
        'direction_sign_changes', 'early_share', 'mid_share', 'late_share'])


def collapse_programs(candidates, gene_sets, leading_edges, curves, *, distance_cut=.45):
    """Complete-linkage descriptive grouping after discovery; retain every term.

    Equal-weight similarities: gene Jaccard, leading-edge Jaccard, divergence
    correlation, direction correlation. Missing leading edges score zero.
    """
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform

    ids = list(candidates)
    if not ids:
        return pd.DataFrame(columns=['pathway_id', 'program'])
    distances = np.zeros((len(ids), len(ids)))
    traces = {p: curves[curves.pathway_id.eq(p)].sort_values('position') for p in ids}
    def jaccard(a, b):
        a, b = set(a), set(b)
        return len(a & b) / len(a | b) if a | b else 0.
    for i, left in enumerate(ids):
        for j in range(i):
            right = ids[j]
            similarities = [jaccard(gene_sets[left], gene_sets[right]),
                            jaccard(leading_edges.get(left, []), leading_edges.get(right, []))]
            for column in ('divergence', 'direction'):
                a, b = traces[left][column].to_numpy(), traces[right][column].to_numpy()
                corr = np.corrcoef(a, b)[0, 1] if np.std(a) > 1e-10 and np.std(b) > 1e-10 else 0.
                similarities.append(max(float(corr), 0.))
            distances[i, j] = distances[j, i] = 1 - np.mean(similarities)
    groups = (fcluster(linkage(squareform(distances), method='complete'), distance_cut, criterion='distance')
              if len(ids) > 1 else np.ones(1, int))
    return pd.DataFrame({'pathway_id': ids, 'program': groups})
