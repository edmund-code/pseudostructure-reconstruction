"""Training-only trajectory baselines for repeated-section atlas benchmarks."""
from __future__ import annotations

import numpy as np
from importlib.metadata import PackageNotFoundError, version
from sklearn.neighbors import NearestNeighbors
from scipy.sparse.csgraph import connected_components

from pseudospace.repeated_atlas_baselines import _anatomy, _matrix


def _version(package):
    try:
        return version(package)
    except PackageNotFoundError:
        return 'unknown'


def _seed(seed):
    if not isinstance(seed, (int, np.integer)) or isinstance(seed, (bool, np.bool_)):
        raise ValueError('seed must be an integer')
    return int(seed)


def _orient(z, anatomy):
    flipped = bool(np.median(z[anatomy == 2]) < np.median(z[anatomy == 0]))
    if flipped:
        z = 1. - z
    lo, hi = float(np.min(z)), float(np.max(z))
    if not np.isfinite([lo, hi]).all() or hi - lo <= np.finfo(float).eps:
        raise ValueError('trajectory produced a non-finite or constant coordinate')
    return np.clip((z - lo) / (hi - lo), 0., 1.), flipped


def fit_dpt_axis(Y, anatomy, *, n_neighbors=30, seed=15):
    """Fit Scanpy DPT on all training rows, rooted at the S1-mean-nearest row."""
    Y = _matrix(Y, 'Y')
    anatomy = _anatomy(anatomy, len(Y), require_all=True)
    seed = _seed(seed)
    if not isinstance(n_neighbors, (int, np.integer)) or isinstance(n_neighbors, (bool, np.bool_)) or n_neighbors < 2:
        raise ValueError('n_neighbors must be an integer >= 2')
    if len(Y) < 4:
        raise ValueError('DPT requires at least four training observations')
    try:
        import anndata as ad
        import scanpy as sc
    except ImportError as exc:  # pragma: no cover - analysis environment dependent
        raise ImportError('fit_dpt_axis requires scanpy and anndata') from exc

    obj = ad.AnnData(X=Y.copy())
    obj.obsm['X_benchmark'] = Y.copy()
    k = min(int(n_neighbors), len(Y) - 1)
    sc.pp.neighbors(obj, n_neighbors=k, use_rep='X_benchmark', random_state=seed)
    n_comps = min(15, len(Y) - 1)
    sc.tl.diffmap(obj, n_comps=n_comps)
    obj.uns['iroot'] = int(np.flatnonzero(anatomy == 0)[np.argmin(
        np.linalg.norm(Y[anatomy == 0] - Y[anatomy == 0].mean(axis=0), axis=1))])
    sc.tl.dpt(obj, n_dcs=min(10, n_comps))
    raw = np.asarray(obj.obs['dpt_pseudotime'], dtype=float)
    n_nonfinite = int((~np.isfinite(raw)).sum())
    if n_nonfinite:
        raise ValueError(f'DPT returned {n_nonfinite} non-finite training coordinates')
    z, flipped = _orient(raw, anatomy)
    return z, {'method': 'scanpy_dpt', 'root_index': int(obj.uns['iroot']),
               'nonfinite_count': n_nonfinite, 'n_neighbors': k,
               'n_components': n_comps, 'n_dcs': min(10, n_comps),
               'seed': seed, 'orientation_flipped': flipped,
               'scanpy_version': _version('scanpy'), 'anndata_version': _version('anndata')}


def _scfates_graph_guard(adata):
    """Refuse branched inferred graphs when ScFates exposes its graph topology."""
    graph_info = adata.uns.get('graph', {})
    adjacency = graph_info.get('B') if isinstance(graph_info, dict) else None
    if adjacency is None:
        raise ValueError('ScFates did not expose an inferred graph for topology validation')
    adjacency = np.asarray(adjacency, dtype=float)
    if adjacency.ndim != 2 or adjacency.shape[0] != adjacency.shape[1]:
        raise ValueError('ScFates graph adjacency is not a square matrix')
    degrees = np.count_nonzero(adjacency, axis=1)
    tips = int(np.sum(degrees == 1))
    forks = int(np.sum(degrees > 2))
    n_components = connected_components(adjacency, directed=False, return_labels=False)
    if tips != 2 or forks != 0 or np.any((degrees < 1) | (degrees > 2)) or n_components != 1:
        raise ValueError(f'ScFates graph must be unbranched with two tips; found {tips} tips and {forks} forks')
    return tips, forks


def fit_scfates_axis(Y, anatomy, *, nodes=30, seed=15):
    """Fit an unbranched ScFates curve and return training pseudotime."""
    Y = _matrix(Y, 'Y')
    anatomy = _anatomy(anatomy, len(Y), require_all=True)
    seed = _seed(seed)
    if not isinstance(nodes, (int, np.integer)) or isinstance(nodes, (bool, np.bool_)) or nodes < 3:
        raise ValueError('nodes must be an integer >= 3')
    try:
        import anndata as ad
        import scFates as scf
    except ImportError as exc:  # pragma: no cover - analysis environment dependent
        raise ImportError('fit_scfates_axis requires scFates and anndata') from exc

    obj = ad.AnnData(X=Y.copy())
    obj.obsm['X_benchmark'] = Y.copy()
    scf.tl.curve(obj, use_rep='X_benchmark', Nodes=int(nodes), seed=seed,
                 epg_lambda=0.01, epg_mu=0.1, basis=None, device='cpu')
    tips, forks = _scfates_graph_guard(obj)
    s1_mean = Y[anatomy == 0].mean(axis=0)
    tip_coords = np.asarray(obj.uns['graph']['tips'], dtype=int)
    curve = np.asarray(obj.uns['graph']['F'], dtype=float)
    root_tip = int(tip_coords[np.argmin(np.linalg.norm(curve[:, tip_coords] - s1_mean[:, None], axis=0))])
    scf.tl.root(obj, int(root_tip))
    scf.tl.pseudotime(obj, n_jobs=2, n_map=1, seed=seed)
    if 't' in obj.obs:
        raw = np.asarray(obj.obs['t'], dtype=float)
    elif 'pseudotime' in obj.obs:
        raw = np.asarray(obj.obs['pseudotime'], dtype=float)
    else:
        raise ValueError('ScFates did not return a recognized pseudotime column')
    n_nonfinite = int((~np.isfinite(raw)).sum())
    if n_nonfinite:
        raise ValueError(f'ScFates returned {n_nonfinite} non-finite training coordinates')
    z, flipped = _orient(raw, anatomy)
    return z, {'method': 'scfates_curve', 'root_tip': root_tip,
               'nonfinite_count': n_nonfinite, 'nodes': int(nodes),
               'seed': seed, 'tips': tips, 'forks': forks,
               'orientation_flipped': flipped,
               'scfates_version': _version('scFates'), 'anndata_version': _version('anndata'),
               'root_rule': 'curve tip nearest training S1 molecular mean'}


def project_neighbor_axis(train_Y, train_z, query_Y, *, k=15):
    """Project query rows by frozen Euclidean kNN inverse-distance weights."""
    train_Y = _matrix(train_Y, 'train_Y')
    query_Y = np.asarray(query_Y, dtype=float)
    if query_Y.ndim != 2 or query_Y.shape[1] != train_Y.shape[1] or not np.isfinite(query_Y).all():
        raise ValueError('query_Y must be finite, 2D, and match training features')
    train_z = np.asarray(train_z, dtype=float).ravel()
    if len(train_z) != len(train_Y) or not np.isfinite(train_z).all() or np.any((train_z < 0) | (train_z > 1)):
        raise ValueError('train_z must be finite in [0,1] and match training rows')
    if not isinstance(k, (int, np.integer)) or isinstance(k, (bool, np.bool_)) or k < 1:
        raise ValueError('k must be a positive integer')
    if not len(query_Y):
        return np.empty(0, dtype=float)
    k = min(int(k), len(train_Y))
    distances, indices = NearestNeighbors(n_neighbors=k).fit(train_Y).kneighbors(query_Y)
    result = np.empty(len(query_Y), dtype=float)
    for i, (d, idx) in enumerate(zip(distances, indices)):
        zero = d <= 1e-12
        if np.any(zero):
            result[i] = train_z[idx[zero]].mean()
        else:
            w = 1. / d
            result[i] = w @ train_z[idx] / w.sum()
    return np.clip(result, 0., 1.)
