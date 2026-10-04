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


def _trajectory_anatomy(anatomy, n, anchor_mask=None):
    """Validate anchor labels and mark non-anchor rows with the -1 sentinel."""
    if anchor_mask is None:
        labels = _anatomy(anatomy, n, require_all=True)
        return labels, np.ones(n, dtype=bool)

    mask = np.asarray(anchor_mask)
    if mask.ndim != 1 or len(mask) != n or mask.dtype.kind != 'b':
        raise ValueError('anchor_mask must be a one-dimensional boolean vector matching rows')
    labels = np.asarray(anatomy)
    if labels.ndim != 1 or len(labels) != n:
        raise ValueError('anatomy must be a one-dimensional vector matching rows')
    anchors = _anatomy(labels[mask], int(mask.sum()), require_all=True)
    validated = np.full(n, -1, dtype=int)
    validated[mask] = anchors
    return validated, mask


def fit_dpt_axis(Y, anatomy, *, anchor_mask=None, n_neighbors=30, seed=15):
    """Fit DPT on all rows, calibrating its root and orientation from anatomy anchors.

    With ``anchor_mask``, unlabeled rows still contribute to the DPT graph and
    receive coordinates, while only masked rows calibrate segment rooting and
    orientation. Without it, every row must retain the strict 0/1/2 labels.
    """
    Y = _matrix(Y, 'Y')
    anatomy, anchor_mask = _trajectory_anatomy(anatomy, len(Y), anchor_mask)
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
    s1_rows = np.flatnonzero(anchor_mask & (anatomy == 0))
    obj.uns['iroot'] = int(s1_rows[np.argmin(
        np.linalg.norm(Y[s1_rows] - Y[s1_rows].mean(axis=0), axis=1))])
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
               'anchor_count': int(anchor_mask.sum()),
               'unlabeled_count': int((~anchor_mask).sum()),
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


def _anatomy_curve_initialization(Y, anatomy, anchor_mask):
    """Build an ordered three-node path from anchor-only feature centroids."""
    positions = np.vstack([
        Y[anchor_mask & (anatomy == label)].mean(axis=0)
        for label in (0, 1, 2)
    ])
    if any(np.array_equal(positions[i], positions[j])
           for i in range(3) for j in range(i + 1, 3)):
        raise ValueError('anatomy initialization requires three distinct anchor centroids')
    return {
        'InitNodePositions': positions,
        'InitEdges': np.array([[0, 1], [1, 2]], dtype=int),
        'InitNodes': 3,
    }


def fit_scfates_axis(Y, anatomy, *, anchor_mask=None, nodes=30, seed=15,
                     initialization='default'):
    """Fit an unbranched curve with source-only anatomical anchors.

    Non-anchor rows remain curve observations and receive coordinates, but do not
    contribute anatomy labels to calibration. By default all rows are strict
    0/1/2 anatomy anchors. ``initialization='anatomy'`` initializes the curve
    with an ordered path through the three anchor-only feature centroids; it
    does not remove any observed fit rows.
    """
    if not isinstance(initialization, str) or initialization not in ('default', 'anatomy'):
        raise ValueError("initialization must be 'default' or 'anatomy'")
    Y = _matrix(Y, 'Y')
    anatomy, anchor_mask = _trajectory_anatomy(anatomy, len(Y), anchor_mask)
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
    if initialization == 'default':
        # Preserve the original backend call and its defaults exactly.
        scf.tl.curve(obj, use_rep='X_benchmark', Nodes=int(nodes), seed=seed,
                     epg_lambda=0.01, epg_mu=0.1, basis=None, device='cpu')
    else:
        init = _anatomy_curve_initialization(Y, anatomy, anchor_mask)
        scf.tl.curve(obj, use_rep='X_benchmark', Nodes=int(nodes), seed=seed,
                     epg_lambda=0.01, epg_mu=0.1, basis=None, device='cpu', **init)
    tips, forks = _scfates_graph_guard(obj)
    s1_mean = Y[anchor_mask & (anatomy == 0)].mean(axis=0)
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
               'initialization': initialization,
               'anchor_count': int(anchor_mask.sum()),
               'unlabeled_count': int((~anchor_mask).sum()),
               'scfates_version': _version('scFates'), 'anndata_version': _version('anndata'),
               'root_rule': 'curve tip nearest training S1 molecular mean'}


def fit_unoriented_scfates_axis(Y, *, nodes=30, seed=15):
    """Fit an unlabeled unbranched curve with an arbitrary native orientation.

    This returns coordinates normalized over this fit's own observed rows.
    That scale describes the fitted query curve only; it does not assert that
    the observations cover the full canonical interval.
    """
    Y = _matrix(Y, 'Y')
    seed = _seed(seed)
    if seed < 0:
        raise ValueError('seed must be non-negative')
    if (not isinstance(nodes, (int, np.integer)) or isinstance(nodes, (bool, np.bool_))
            or nodes < 3):
        raise ValueError('nodes must be an integer >= 3')
    try:
        import anndata as ad
        import scFates as scf
    except ImportError as exc:  # pragma: no cover - analysis environment dependent
        raise ImportError('fit_unoriented_scfates_axis requires scFates and anndata') from exc

    obj = ad.AnnData(X=Y.copy())
    obj.obsm['X_benchmark'] = Y.copy()
    scf.tl.curve(obj, use_rep='X_benchmark', Nodes=int(nodes), seed=seed,
                 epg_lambda=0.01, epg_mu=0.1, basis=None, device='cpu')
    tips, forks = _scfates_graph_guard(obj)
    tip_coords = np.asarray(obj.uns['graph']['tips'], dtype=int)
    if tip_coords.ndim != 1 or len(tip_coords) != 2:
        raise ValueError('ScFates graph must expose exactly two tip indices')
    root_tip = int(np.min(tip_coords))
    scf.tl.root(obj, root_tip)
    scf.tl.pseudotime(obj, n_jobs=2, n_map=1, seed=seed)
    if 't' in obj.obs:
        raw = np.asarray(obj.obs['t'], dtype=float)
    elif 'pseudotime' in obj.obs:
        raw = np.asarray(obj.obs['pseudotime'], dtype=float)
    else:
        raise ValueError('ScFates did not return a recognized pseudotime column')
    if raw.ndim != 1 or len(raw) != len(Y) or not np.isfinite(raw).all():
        raise ValueError('ScFates returned invalid unlabeled pseudotime values')
    lo, hi = float(raw.min()), float(raw.max())
    if hi - lo <= np.finfo(float).eps:
        raise ValueError('trajectory produced a constant coordinate')
    z = np.clip((raw - lo) / (hi - lo), 0., 1.)
    return z, {
        'method': 'scfates_unoriented_curve', 'root_tip': root_tip,
        'root_rule': 'lowest-index inferred graph tip; arbitrary orientation, no anatomy',
        'nonfinite_count': 0, 'nodes': int(nodes), 'seed': seed,
        'tips': tips, 'forks': forks,
        'scfates_version': _version('scFates'), 'anndata_version': _version('anndata'),
    }


def project_neighbor_axis(train_Y, train_z, query_Y, *, k=15,
                          return_distribution=False):
    """Project by frozen Euclidean kNN inverse-distance weights.

    The optional distribution returns the existing assignment supports/weights
    alongside their coordinate mean. These weights are not calibrated spatial
    probabilities; retaining them permits fair nonlinear response averaging.
    """
    if not isinstance(return_distribution, (bool, np.bool_)):
        raise ValueError('return_distribution must be boolean')
    train_Y = _matrix(train_Y, 'train_Y')
    query_Y = np.asarray(query_Y, dtype=float)
    if query_Y.ndim != 2 or query_Y.shape[1] != train_Y.shape[1] or not np.isfinite(query_Y).all():
        raise ValueError('query_Y must be finite, 2D, and match training features')
    train_z = np.asarray(train_z, dtype=float).ravel()
    if len(train_z) != len(train_Y) or not np.isfinite(train_z).all() or np.any((train_z < 0) | (train_z > 1)):
        raise ValueError('train_z must be finite in [0,1] and match training rows')
    if not isinstance(k, (int, np.integer)) or isinstance(k, (bool, np.bool_)) or k < 1:
        raise ValueError('k must be a positive integer')
    k = min(int(k), len(train_Y))
    if not len(query_Y):
        empty = np.empty(0, dtype=float)
        return ({'z_mean': empty, 'support': np.empty((0, k)),
                 'weights': np.empty((0, k))} if return_distribution else empty)
    distances, indices = NearestNeighbors(n_neighbors=k).fit(train_Y).kneighbors(query_Y)
    result = np.empty(len(query_Y), dtype=float)
    weights = np.zeros(distances.shape) if return_distribution else None
    for i, (d, idx) in enumerate(zip(distances, indices)):
        zero = d <= 1e-12
        if np.any(zero):
            result[i] = train_z[idx[zero]].mean()
            if return_distribution:
                weights[i, zero] = 1. / zero.sum()
        else:
            w = 1. / d
            result[i] = w @ train_z[idx] / w.sum()
            if return_distribution:
                weights[i] = w / w.sum()
    result = np.clip(result, 0., 1.)
    return ({'z_mean': result, 'support': train_z[indices], 'weights': weights}
            if return_distribution else result)
