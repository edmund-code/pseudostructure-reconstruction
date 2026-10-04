"""Coordinate arms and their comparison summaries (notebook 36).

The arms reuse notebook 13's machinery: species-intersection HVGs, PCA, Harmony, a nonbranching
scFates curve and DPT from the same root. They differ in which genes and structures enter.
scFates is numerically sensitive (ulp-level input changes move some structures), so callers should
report seed stability next to any scFates result. Positions are a gauge on [0, 1], not lengths.
"""
from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.stats import spearmanr

from .stats_gam import gam_internal_knots, make_gam_design


def scanpy_lognorm(counts):
    """Notebook 13's normalization: sorted CSR, normalize_total(1e4) by division, then log1p."""
    block = sparse.csr_matrix(counts, dtype=np.float64)
    block.sort_indices()
    block.eliminate_zeros()
    obj = ad.AnnData(X=block)
    sc.pp.normalize_total(obj, target_sum=1e4)
    sc.pp.log1p(obj)
    return sparse.csr_matrix(obj.X)


def species_intersection_hvgs(lognorm, obs, genes, exclude):
    """Notebook 13's Harmony gene rule: HVG (seurat, batch = sample) in both species."""
    from .harmony import select_harmony_hvgs_by_condition
    var = pd.DataFrame({'exclude_from_harmony_hvg': np.asarray(exclude, bool)}, index=pd.Index(genes))
    obj = ad.AnnData(X=lognorm, obs=obs.copy(), var=var)
    obj.layers['lognorm'] = obj.X
    obj = select_harmony_hvgs_by_condition(obj, group_key='comparison_species', groups=('mouse', 'human'),
                                           mode='intersection', min_mean=.0125, max_mean=3, min_disp=.5)
    return obj.var.highly_variable_for_harmony.to_numpy(bool)


def harmony_trajectory(lognorm, obs, columns, *, pca_mask=None, n_comps=50, theta=6, seed=0, dims=5):
    """PCA on the selected columns, Harmony over ``sample``, first ``dims`` corrected components.

    ``pca_mask`` (boolean over ``columns``) reproduces scanpy's ``mask_var='highly_variable'``
    default when the object carries such a column; ``None`` uses every selected column.
    """
    import rpy2.robjects as ro
    from .harmony import run_harmony_rpy2
    columns = np.asarray(columns)
    obj = ad.AnnData(X=sparse.csr_matrix(lognorm[:, columns]), obs=obs.copy())
    if pca_mask is not None:
        obj.var['highly_variable'] = np.asarray(pca_mask, bool)
        sc.tl.pca(obj, n_comps=n_comps, random_state=seed)
    else:
        sc.tl.pca(obj, n_comps=n_comps, random_state=seed, mask_var=None)
    ro.r(f'set.seed({seed})')
    obj = run_harmony_rpy2(obj, batch_key='sample', n_pcs=n_comps, theta=theta, lambda_val=1, max_iter=30, tau=0)
    return np.asarray(obj.obsm['X_harmony'][:, :dims], dtype=float)


def marker_tip_score(values, species, early_columns, late_columns):
    """Notebook 13's early-tip score: per-species z-scored early minus late marker means."""
    values = np.asarray(values, float)
    species = np.asarray(species).astype(str)
    z = np.empty_like(values)
    for sp in np.unique(species):
        rows = species == sp
        z[rows] = (values[rows] - values[rows].mean(axis=0)) / np.maximum(values[rows].std(axis=0), 1e-6)
    return z[:, early_columns].mean(axis=1) - z[:, late_columns].mean(axis=1)


def fit_scfates_curve(X, tip_score, *, nodes=30, seed=0, n_map=0, agree_groups=None):
    """Nonbranching scFates curve with notebook 13's root rule, scaled to [0, 1].

    Without ``agree_groups`` the root is scFates' automatic tip choice on ``tip_score``. With it,
    each tip is scored by the mean ``tip_score`` of its top-1% members within every group, and all
    groups must pick the same tip. The DPT root cell is the near-root structure (top 1% membership)
    with the highest score.
    """
    import scFates as scf
    X = np.asarray(X, float)
    score = np.asarray(tip_score, float)
    obj = ad.AnnData(X=np.zeros((len(X), 1)))
    obj.obsm['X_trajectory'] = X
    obj.obs['early_tip_score'] = score
    scf.tl.curve(obj, Nodes=nodes, use_rep='X_trajectory', ndims_rep=X.shape[1], seed=seed)
    graph = obj.uns['graph']
    degree = np.asarray(graph['B']).astype(bool).sum(axis=0)
    if not (len(graph['tips']) == 2 and len(graph['forks']) == 0 and np.all(degree <= 2)):
        raise ValueError('scFates graph is not one nonbranching path.')
    membership = np.asarray(obj.obsm['X_R'])
    if agree_groups is None:
        scf.tl.root(obj, 'early_tip_score', tips_only=True)
    else:
        groups = np.asarray(agree_groups).astype(str)
        choices = {}
        for g in np.unique(groups):
            rows = np.flatnonzero(groups == g)
            means = []
            for tip in graph['tips']:
                m = membership[rows, tip]
                top = rows[m >= np.quantile(m, .99)]
                means.append(score[top].mean())
            choices[g] = int(np.asarray(graph['tips'])[int(np.argmax(means))])
        if len(set(choices.values())) != 1:
            raise ValueError(f'Groups disagree on the early tip: {choices}')
        scf.tl.root(obj, next(iter(choices.values())))
    root_node = int(obj.uns['graph']['root'])
    # scFates' tips_only rule zeroes interior nodes, so all-negative tip scores pick an interior node.
    if root_node not in set(np.asarray(graph['tips']).tolist()):
        raise ValueError('Root is not a tip; the tip score must be positive at the early tip.')
    m = membership[:, root_node]
    near = np.flatnonzero(m >= np.quantile(m, .99))
    root_cell = int(near[np.argmax(score[near])])
    scf.tl.pseudotime(obj, n_jobs=2, n_map=1, seed=seed)
    t = obj.obs.t.to_numpy(float)
    result = {'z': (t - t.min()) / np.ptp(t), 'root_cell': root_cell, 'root_node': root_node}
    if n_map:
        scf.tl.pseudotime(obj, n_jobs=8, n_map=n_map, seed=seed)
        maps = obj.uns['pseudotime_list']
        result['maps'] = np.vstack([(maps[str(i)].loc[obj.obs_names, 't'].to_numpy(float) - t.min()) / np.ptp(t)
                                    for i in range(n_map)])
    return result


def fit_dpt(X, root_cell, *, n_neighbors=30, seed=0):
    """Scanpy DPT on the same representation and root, scaled to [0, 1] (NaN if disconnected)."""
    obj = ad.AnnData(X=np.zeros((len(X), 1)))
    obj.obsm['X_trajectory'] = np.asarray(X, float)
    sc.pp.neighbors(obj, n_neighbors=n_neighbors, use_rep='X_trajectory', random_state=seed,
                    key_added='trajectory_neighbors')
    sc.tl.diffmap(obj, neighbors_key='trajectory_neighbors', random_state=seed)
    obj.uns['iroot'] = int(root_cell)
    sc.tl.dpt(obj, neighbors_key='trajectory_neighbors')
    d = obj.obs.dpt_pseudotime.to_numpy(float)
    finite = np.isfinite(d)
    out = np.full(len(d), np.nan)
    out[finite] = (d[finite] - d[finite].min()) / np.ptp(d[finite])
    return out


def loso_positions(X, tip_score, specimen, held_out, *, nodes=30, seed=0, k=15, agree_groups=None):
    """Fit the curve without one specimen and place it by frozen kNN interpolation (P3b)."""
    from .trajectory_benchmark import project_neighbor_axis
    specimen = np.asarray(specimen).astype(str)
    train = specimen != str(held_out)
    if train.all() or not train.any():
        raise ValueError('held_out must name one present specimen.')
    groups = None if agree_groups is None else np.asarray(agree_groups)[train]
    fit = fit_scfates_curve(np.asarray(X)[train], np.asarray(tip_score)[train], nodes=nodes, seed=seed,
                            agree_groups=groups)
    return project_neighbor_axis(np.asarray(X)[train], fit['z'], np.asarray(X)[~train], k=k)


def within_segment_agreement(a, b, specimen, segment):
    """Spearman of two coordinates inside every specimen x segment cell."""
    frame = pd.DataFrame({'a': a, 'b': b, 'specimen': np.asarray(specimen).astype(str),
                          'segment': np.asarray(segment).astype(str)})
    rows = []
    for (name, seg), cell in frame.groupby(['specimen', 'segment']):
        rho = spearmanr(cell.a, cell.b, nan_policy='omit').statistic if len(cell) > 2 else np.nan
        rows.append({'specimen': name, 'segment': seg, 'rho': float(rho)})
    return pd.DataFrame(rows)


def registration_gap(z, species, specimen, segment, *, first='human', second='mouse'):
    """Per-segment difference of species medians (specimen medians averaged within species)."""
    frame = pd.DataFrame({'z': z, 'species': np.asarray(species).astype(str),
                          'specimen': np.asarray(specimen).astype(str), 'segment': np.asarray(segment).astype(str)})
    medians = frame.groupby(['species', 'specimen', 'segment']).z.median().groupby(['species', 'segment']).mean()
    gap = medians.loc[first] - medians.loc[second]
    return gap, float(gap.abs().max())


def heldout_prediction(train_y, train_z, train_segment, test_y, test_z, test_segment, *, basis_df=6):
    """Per-gene test MSE for a coordinate spline and for training segment means.

    The spline (notebook 12's cubic B-spline basis) is fitted on the training specimen's positions;
    test positions are clipped to the training range. Returns MSEs and the training variance used
    for scaling.
    """
    ytr = np.asarray(train_y.toarray() if sparse.issparse(train_y) else train_y, float)
    yte = np.asarray(test_y.toarray() if sparse.issparse(test_y) else test_y, float)
    ztr, zte = np.asarray(train_z, float), np.asarray(test_z, float)
    keep_tr, keep_te = np.isfinite(ztr), np.isfinite(zte)
    ytr, ztr, seg_tr = ytr[keep_tr], ztr[keep_tr], np.asarray(train_segment).astype(str)[keep_tr]
    yte, zte, seg_te = yte[keep_te], zte[keep_te], np.asarray(test_segment).astype(str)[keep_te]
    knots = gam_internal_knots(ztr, basis_df=basis_df)
    design_tr = make_gam_design(ztr, knots)
    design_te = make_gam_design(np.clip(zte, ztr.min(), ztr.max()), knots)
    beta, *_ = np.linalg.lstsq(design_tr, ytr, rcond=None)
    coordinate_mse = np.mean((yte - design_te @ beta) ** 2, axis=0)
    means = {s: ytr[seg_tr == s].mean(axis=0) for s in np.unique(seg_tr)}
    oracle = np.vstack([means[s] for s in seg_te])
    segment_mse = np.mean((yte - oracle) ** 2, axis=0)
    return {'coordinate_mse': coordinate_mse, 'segment_mse': segment_mse, 'train_var': ytr.var(axis=0)}
