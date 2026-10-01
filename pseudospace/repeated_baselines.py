"""Simple baselines and diagnostics for repeated-structure pseudospace work.

These helpers are deliberately conventional.  In particular, ``panel_dpt`` is
an inductive DPT baseline on a fixed gene panel, not the project's canonical
Harmony DPT workflow.
"""
from __future__ import annotations

import numpy as np
import warnings
from scipy.interpolate import UnivariateSpline
from scipy.sparse import csr_matrix, diags
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import MatrixRankWarning, spsolve
from scipy.stats import pearsonr, spearmanr
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler


def _matrix(X, name):
    X = np.asarray(X, dtype=float)
    if X.ndim != 2 or X.shape[0] < 2 or X.shape[1] < 1:
        raise ValueError(f'{name} must be a 2D matrix with at least two rows and one column')
    if not np.isfinite(X).all():
        raise ValueError(f'{name} must contain only finite values')
    return X


def _anatomy_order(labels):
    a = np.asarray(labels).astype(str)
    order = {'S1': 0.0, 'S2': 0.5, 'S3': 1.0}
    unknown = sorted(set(a) - set(order))
    if unknown:
        raise ValueError(f'anatomy labels must be S1/S2/S3; found {unknown}')
    return np.array([order[x] for x in a])


def _orient_scale(train, test, anatomy_train):
    """Orient and scale from training values only; held-out values are clamped."""
    train = np.asarray(train, float).ravel()
    test = np.asarray(test, float).ravel()
    anat = _anatomy_order(anatomy_train)
    if len(anat) != len(train):
        raise ValueError('anatomy_train length must match training observations')
    corr = spearmanr(train, anat).correlation
    if np.isfinite(corr) and corr < 0:
        train, test = -train, -test
    lo, hi = float(np.min(train)), float(np.max(train))
    if hi - lo <= np.finfo(float).eps:
        raise ValueError('training coordinate has no variation')
    return np.clip((train - lo) / (hi - lo), 0, 1), np.clip((test - lo) / (hi - lo), 0, 1)


def baseline_coordinates(Xtrain, Xtest, anatomy_train, method='pca1', seed=0):
    """Return train/test coordinates from a training-only PCA1 or panel-DPT fit.

    The panel-DPT graph and diffusion pseudotime are fitted on training samples.
    Held-out samples are projected by inverse-distance-weighted interpolation
    from their nearest training neighbors in the training PCA space.
    """
    Xtrain = _matrix(Xtrain, 'Xtrain')
    Xtest = np.asarray(Xtest, dtype=float)
    if Xtest.ndim != 2 or Xtest.shape[1] != Xtrain.shape[1] or not np.isfinite(Xtest).all():
        raise ValueError('Xtest must be finite, 2D, and have the same features as Xtrain')
    if len(anatomy_train) != len(Xtrain):
        raise ValueError('anatomy_train length must match Xtrain')
    if method not in {'pca1', 'panel_dpt'}:
        raise ValueError("method must be 'pca1' or 'panel_dpt'")
    scaler = StandardScaler().fit(Xtrain)
    A, B = scaler.transform(Xtrain), scaler.transform(Xtest)
    requested_components = 1 if method == 'pca1' else 20
    n_components = min(requested_components, A.shape[0] - 1, A.shape[1])
    pca = PCA(n_components=max(1, n_components), random_state=seed).fit(A)
    train_embedding, test_embedding = pca.transform(A), pca.transform(B)
    if method == 'pca1':
        return _orient_scale(train_embedding[:, 0], test_embedding[:, 0], anatomy_train)

    # Scanpy is optional for the core package; import only for this baseline.
    try:
        import anndata as ad
        import scanpy as sc
    except ImportError as exc:  # pragma: no cover - depends on analysis environment
        raise ImportError("method='panel_dpt' requires scanpy and anndata") from exc
    if len(Xtrain) < 4:
        raise ValueError('panel_dpt requires at least four training observations')
    obj = ad.AnnData(X=A.astype(np.float32))
    obj.obsm['X_panel_pca'] = train_embedding.astype(np.float32)
    n_neighbors = min(15, len(Xtrain) - 1)
    sc.pp.neighbors(obj, n_neighbors=n_neighbors, use_rep='X_panel_pca', random_state=seed)
    sc.tl.diffmap(obj)
    anatomy_rank = _anatomy_order(anatomy_train)
    s1 = np.flatnonzero(anatomy_rank == 0)
    if not len(s1):
        raise ValueError('panel_dpt requires at least one S1 training observation')
    root_centroid = train_embedding[s1].mean(axis=0)
    root = int(s1[np.argmin(np.linalg.norm(train_embedding[s1] - root_centroid, axis=1))])
    obj.uns['iroot'] = root
    sc.tl.dpt(obj)
    ztrain = np.asarray(obj.obs['dpt_pseudotime'], dtype=float)
    if not np.isfinite(ztrain).all():
        raise ValueError('panel_dpt produced non-finite training coordinates')
    if len(Xtest):
        k = min(5, len(Xtrain))
        nn = NearestNeighbors(n_neighbors=k).fit(train_embedding)
        distances, indices = nn.kneighbors(test_embedding)
        exact = distances[:, 0] <= 1e-12
        weights = 1 / np.maximum(distances, 1e-12)
        ztest = np.sum(weights * ztrain[indices], axis=1) / weights.sum(axis=1)
        ztest[exact] = ztrain[indices[exact, 0]]
    else:
        ztest = np.empty(0, dtype=float)
    return _orient_scale(ztrain, ztest, anatomy_train)


def anatomy_initialize(X, anatomy, seed=0):
    """Initialize rank positions by PCA1 within each S1/S2/S3 interval.

    The within-interval direction is intentionally arbitrary and uses no marker
    genes.  The thirds encode only the supplied coarse anatomy.
    """
    X = _matrix(X, 'X')
    anatomy = np.asarray(anatomy).astype(str)
    if len(anatomy) != len(X):
        raise ValueError('anatomy length must match X observations')
    anat_rank = _anatomy_order(anatomy)
    z = np.full(len(X), np.nan)
    for label, lower, upper in [('S1', 0, 1/3), ('S2', 1/3, 2/3), ('S3', 2/3, 1)]:
        idx = np.flatnonzero(anatomy == label)
        if not len(idx):
            continue
        if len(idx) == 1:
            ranks = np.array([0.5])
        else:
            scaled = StandardScaler().fit_transform(X[idx])
            pc = PCA(n_components=1, random_state=seed).fit_transform(scaled)[:, 0]
            ranks = (np.argsort(np.argsort(pc, kind='mergesort'), kind='mergesort') + .5) / len(idx)
        z[idx] = lower + ranks * (upper - lower)
    return z


def curve_agreement(X, z, specimens, gene_names, n_grid=50):
    """Summarize pairwise agreement of specimen-specific smoothed gene curves.

    Fits specimen-specific smoothing splines with a fixed mild smoothing
    penalty (rather than the workflow atlas fitter). Returns a long DataFrame with Pearson curve agreement over each pair's
    overlapping coordinate support, peak displacement, transition displacement,
    and monotonic direction agreement. A specimen needs four distinct finite
    positions for a spline fit.
    """
    import pandas as pd

    X = _matrix(X, 'X')
    z = np.asarray(z, dtype=float).ravel()
    specimens = np.asarray(specimens).astype(str)
    genes = np.asarray(gene_names).astype(str)
    if len(z) != len(X) or len(specimens) != len(X) or len(genes) != X.shape[1]:
        raise ValueError('z/specimens/gene_names dimensions must match X')
    if n_grid < 5:
        raise ValueError('n_grid must be at least 5')
    groups = np.unique(specimens)
    if len(groups) < 2:
        raise ValueError('at least two specimens are required')
    fits = {}
    for gene_idx, gene in enumerate(genes):
        for specimen in groups:
            mask = (specimens == specimen) & np.isfinite(z)
            xx, yy = z[mask], X[mask, gene_idx]
            if len(xx) < 4 or np.ptp(xx) <= 0:
                continue
            order = np.argsort(xx)
            xx, yy = xx[order], yy[order]
            # Average repeated positions before fitting; use mild smoothing.
            ux, inverse = np.unique(xx, return_inverse=True)
            if len(ux) < 4:
                continue
            uy = np.array([yy[inverse == j].mean() for j in range(len(ux))])
            spline = UnivariateSpline(ux, uy, k=min(3, len(ux)-1), s=len(ux) * np.var(uy) * .05)
            grid = np.linspace(ux.min(), ux.max(), n_grid)
            curve = spline(grid)
            derivative = spline.derivative()(grid)
            # Peak is maximum absolute deviation from specimen-specific mean.
            peak = float(grid[np.argmax(np.abs(curve - np.mean(curve)))])
            transition = float(grid[np.argmax(np.abs(derivative))])
            sign = float(np.sign(np.corrcoef(grid, curve)[0, 1])) if np.std(curve) else 0.0
            fits[(gene, specimen)] = (grid, curve, peak, transition, sign)
    rows = []
    for gene in genes:
        for i, a in enumerate(groups):
            for b in groups[i+1:]:
                if (gene, a) not in fits or (gene, b) not in fits:
                    continue
                ga, ya, pa, ta, sa = fits[(gene, a)]
                gb, yb, pb, tb, sb = fits[(gene, b)]
                lo, hi = max(ga[0], gb[0]), min(ga[-1], gb[-1])
                if hi <= lo:
                    continue
                common = np.linspace(lo, hi, n_grid)
                va = np.interp(common, ga, ya)
                vb = np.interp(common, gb, yb)
                corr = pearsonr(va, vb).statistic if np.std(va) and np.std(vb) else np.nan
                peak_a = common[np.argmax(np.abs(va - np.mean(va)))]
                peak_b = common[np.argmax(np.abs(vb - np.mean(vb)))]
                trans_a = common[np.argmax(np.abs(np.gradient(va, common)))]
                trans_b = common[np.argmax(np.abs(np.gradient(vb, common)))]
                sign_a = np.sign(np.corrcoef(common, va)[0, 1]) if np.std(va) else 0
                sign_b = np.sign(np.corrcoef(common, vb)[0, 1]) if np.std(vb) else 0
                rows.append({'gene': gene, 'specimen_a': a, 'specimen_b': b,
                    'curve_pearson': float(corr), 'peak_difference': float(abs(peak_a-peak_b)),
                    'transition_difference': float(abs(trans_a-trans_b)),
                    'monotonic_sign_agreement': bool(sign_a == sign_b) if sign_a and sign_b else np.nan,
                    'overlap_min': lo, 'overlap_max': hi})
    return pd.DataFrame(rows)


def topology_coordinate(X, anatomy, k=15):
    """Harmonic graph coordinate with S1/S3 boundaries; disconnected gaps stay NaN.

    Returns ``(coordinate, diagnostics)``. S1-S3 edges are prohibited; components
    lacking either boundary are not bridged and retain NaN coordinates.
    """
    X = _matrix(X, 'X')
    anatomy = np.asarray(anatomy).astype(str)
    if len(anatomy) != len(X):
        raise ValueError('anatomy length must match X observations')
    labels = _anatomy_order(anatomy)
    if int(k) < 1:
        raise ValueError('k must be positive')
    emb = PCA(n_components=min(10, X.shape[0]-1, X.shape[1]), random_state=0).fit_transform(StandardScaler().fit_transform(X))
    nn = NearestNeighbors(n_neighbors=min(int(k)+1, len(X))).fit(emb)
    distances, neighbors = nn.kneighbors(emb)
    positive_distances = distances[:, 1:][distances[:, 1:] > 0]
    distance_scale = float(np.median(positive_distances)) if positive_distances.size else 1.0
    distance_scale = max(distance_scale, 1e-12)
    rows, cols, vals = [], [], []
    removed_edges = 0
    for i in range(len(X)):
        for d, j in zip(distances[i, 1:], neighbors[i, 1:]):
            # Do not create shortcuts directly between S1 and S3.
            if {labels[i], labels[j]} == {0.0, 1.0}:
                removed_edges += 1
                continue
            w = float(np.exp(-d*d / distance_scale**2))
            rows.extend([i, int(j)]); cols.extend([int(j), i]); vals.extend([w, w])
    W = csr_matrix((vals, (rows, cols)), shape=(len(X), len(X)))
    ncomp, comp = connected_components(W, directed=False)
    z = np.full(len(X), np.nan)
    component_info = []
    n_harmonic_range_failures = 0
    for c in range(ncomp):
        idx = np.flatnonzero(comp == c)
        lo = idx[labels[idx] == 0]
        hi = idx[labels[idx] == 1]
        valid = bool(len(lo) and len(hi))
        component_info.append({'component': int(c), 'n': int(len(idx)), 'has_s1': bool(len(lo)), 'has_s3': bool(len(hi)), 'coordinate_defined': valid})
        if not valid:
            continue
        boundary = np.unique(np.concatenate([lo, hi]))
        # Resolve any overlap (not expected under disjoint labels).
        interior = np.setdiff1d(idx, boundary, assume_unique=False)
        z[lo] = 0.0; z[hi] = 1.0
        if not len(interior):
            continue
        Wc = W[idx][:, idx].tocsr()
        L = diags(np.asarray(Wc.sum(axis=1)).ravel()) - Wc
        # Indices are local positions because idx is sorted by flatnonzero.
        local = np.full(len(X), -1, dtype=int)
        local[idx] = np.arange(len(idx))
        b_loc = local[boundary]
        i_loc = local[interior]
        b_values = np.where(labels[boundary] == 0, 0.0, 1.0)
        try:
            L_ii = L[i_loc][:, i_loc].tocsr()
            rhs = -L[i_loc][:, b_loc] @ b_values
            with warnings.catch_warnings():
                warnings.simplefilter('error', MatrixRankWarning)
                solution = np.asarray(spsolve(L_ii, rhs), dtype=float).ravel()
            if not np.isfinite(solution).all() or np.any(solution < -1e-6) or np.any(solution > 1 + 1e-6):
                z[interior] = np.nan
                n_harmonic_range_failures += 1
            else:
                # The harmonic maximum principle bounds values by the endpoint
                # values; clip only floating-point solver error at the boundary.
                z[interior] = np.clip(solution, 0.0, 1.0)
        except (np.linalg.LinAlgError, MatrixRankWarning, RuntimeError, ValueError):
            z[interior] = np.nan
    return z, {'n_components': int(ncomp), 'components': component_info,
               'n_unmapped': int(np.isnan(z).sum()), 'n_s1_s3_edges_removed': int(removed_edges),
               'n_harmonic_range_failures': n_harmonic_range_failures}
