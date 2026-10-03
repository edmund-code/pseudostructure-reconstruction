"""Small probability-flow proof-of-concept tools for spatial atlases.

This is a finite-step OT distillation pilot.  It does not enforce the
continuous continuity equation and its alternating latent fit is not an EM
optimizer; it is intended for bounded synthetic and exploratory experiments.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.spatial import cKDTree
from scipy.special import expit, logsumexp
from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsRegressor
from sklearn.preprocessing import StandardScaler


def scaffold_centroid_rows(query_xy, query_samples, reference_xy, reference_samples,
                           tolerance=1.):
    """Match query centroids to unique reference rows within the same specimen.

    Returns reference row indices in query order, using ``-1`` for a query
    with no reference centroid within ``tolerance``. Coordinates are matched
    geometrically; observation names are deliberately not consulted.
    """
    query_xy = np.asarray(query_xy, dtype=float)
    reference_xy = np.asarray(reference_xy, dtype=float)
    if query_xy.ndim != 2 or query_xy.shape[1:] != (2,) or not np.isfinite(query_xy).all():
        raise ValueError('query_xy must be a finite numeric N x 2 array')
    if reference_xy.ndim != 2 or reference_xy.shape[1:] != (2,) or not np.isfinite(reference_xy).all():
        raise ValueError('reference_xy must be a finite numeric N x 2 array')
    query_samples = np.asarray(query_samples, dtype=object).ravel()
    reference_samples = np.asarray(reference_samples, dtype=object).ravel()
    if len(query_samples) != len(query_xy):
        raise ValueError('query_samples length must match query_xy')
    if len(reference_samples) != len(reference_xy):
        raise ValueError('reference_samples length must match reference_xy')
    if any(pd.isna(x) or not str(x).strip() for x in query_samples):
        raise ValueError('query_samples must contain nonempty specimen labels')
    if any(pd.isna(x) or not str(x).strip() for x in reference_samples):
        raise ValueError('reference_samples must contain nonempty specimen labels')
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError('tolerance must be positive and finite')

    matched = np.full(len(query_xy), -1, dtype=int)
    # Build one small tree per specimen so coordinates cannot match across samples.
    for sample in dict.fromkeys(query_samples.tolist()):
        qrows = np.flatnonzero(query_samples == sample)
        rrows = np.flatnonzero(reference_samples == sample)
        if not len(rrows):
            continue
        coords = reference_xy[rrows]
        tree = cKDTree(coords)
        for qrow in qrows:
            candidates = tree.query_ball_point(query_xy[qrow], tolerance)
            if len(candidates) > 1:
                raise ValueError(f'query row {qrow} has ambiguous reference centroids within tolerance')
            if candidates:
                matched[qrow] = rrows[candidates[0]]
        assigned = matched[qrows]
        assigned = assigned[assigned >= 0]
        if len(np.unique(assigned)) != len(assigned):
            raise ValueError(f'multiple query rows reuse a reference centroid in specimen {sample!r}')
    return matched


def _matrix(X, name, min_rows=1):
    X = np.asarray(X, dtype=float)
    if X.ndim != 2 or X.shape[0] < min_rows or X.shape[1] < 1:
        raise ValueError(f'{name} must be a 2D matrix with at least {min_rows} rows and one feature')
    if not np.isfinite(X).all():
        raise ValueError(f'{name} must contain only finite values')
    return X


def representation(Xtrain, Xtest, n_components=15, seed=15):
    """Fit a training-only scaler/PCA representation and transform held-out rows."""
    Xtrain = _matrix(Xtrain, 'Xtrain', 2)
    Xtest = np.asarray(Xtest, dtype=float)
    if Xtest.ndim != 2 or Xtest.shape[1] != Xtrain.shape[1] or not np.isfinite(Xtest).all():
        raise ValueError('Xtest must be finite and have the same number of features as Xtrain')
    if not isinstance(n_components, (int, np.integer)) or n_components < 1:
        raise ValueError('n_components must be a positive integer')
    scaler = StandardScaler().fit(Xtrain)
    A = scaler.transform(Xtrain)
    k = min(int(n_components), len(Xtrain) - 1, Xtrain.shape[1])
    if k < 1:
        raise ValueError('at least two training rows are required for PCA')
    pca = PCA(n_components=k, random_state=seed).fit(A)
    B = (pca.transform(scaler.transform(Xtest)) if len(Xtest)
         else np.empty((0, k)))
    return pca.transform(A), B, {'scaler': scaler, 'pca': pca}


def _centers(centers):
    c = np.asarray(centers, dtype=float).ravel()
    if len(c) < 2 or not np.isfinite(c).all() or np.any(np.diff(c) <= 0):
        raise ValueError('centers must be a finite, strictly increasing vector with at least two entries')
    return c


def window_clouds(Y, z, centers, halfwidth, min_count=12, max_points=128, seed=15,
                  specimens=None):
    """Build bounded empirical clouds, optionally balanced across specimens."""
    Y = _matrix(Y, 'Y')
    z = np.asarray(z, float).ravel()
    c = _centers(centers)
    if len(z) != len(Y) or not np.isfinite(z).all():
        raise ValueError('z must be finite and match Y rows')
    if specimens is not None:
        specimens = np.asarray(specimens, dtype=object).ravel()
        if len(specimens) != len(Y) or any(pd.isna(x) or not str(x).strip() for x in specimens):
            raise ValueError('specimens must contain one nonempty, nonmissing label per Y row')
        specimens = np.asarray([str(x) for x in specimens], dtype=str)
    if not np.isfinite(halfwidth) or halfwidth <= 0 or min_count < 2 or max_points < min_count:
        raise ValueError('halfwidth must be positive, min_count >= 2, and max_points >= min_count')
    rng = np.random.default_rng(seed)
    clouds = []
    for window_i, center in enumerate(c):
        idx = np.flatnonzero(np.abs(z - center) <= halfwidth)
        if specimens is None:
            if len(idx) < min_count:
                raise ValueError(f'window at {center:g} has {len(idx)} rows; requires {min_count}')
            if len(idx) > max_points:
                idx = np.sort(rng.choice(idx, size=max_points, replace=False))
        else:
            groups = [idx[specimens[idx] == s] for s in np.unique(specimens)]
            counts = [len(group) for group in groups]
            if any(n < min_count for n in counts):
                raise ValueError(f'window at {center:g} has per-specimen counts {counts}; each requires {min_count}')
            per_specimen = min(min(counts), max_points // len(groups))
            if per_specimen < min_count:
                raise ValueError(f'window at {center:g} cannot retain {min_count} rows per specimen under max_points={max_points}')
            selected = []
            for group in groups:
                selected.append(np.sort(rng.choice(group, size=per_specimen, replace=False)))
            idx = np.concatenate(selected)
        clouds.append(Y[idx].copy())
    return clouds


def _sinkhorn_barycenters(source, target, epsilon, max_iter=10000, tol=1e-4):
    if epsilon <= 0 or not np.isfinite(epsilon):
        raise ValueError('epsilon must be positive and finite')
    cost = np.sum((source[:, None, :] - target[None, :, :]) ** 2, axis=2)
    # Entropic balanced plan in log space with uniform marginals.
    # epsilon is relative to the median nonzero squared distance, so one
    # setting remains usable across standardized and PCA-scaled representations.
    positive = cost[cost > 0]
    scale = float(np.median(positive)) if len(positive) else 1.
    logK = -cost / (epsilon * scale)
    loga = np.full(len(source), -np.log(len(source)))
    logb = np.full(len(target), -np.log(len(target)))
    u, v = np.zeros_like(loga), np.zeros_like(logb)
    converged = False
    row_resid = col_resid = np.inf
    for iteration in range(max_iter):
        u = loga - logsumexp(logK + v[None, :], axis=1)
        v = logb - logsumexp(logK + u[:, None], axis=0)
        # Dual changes depend on arbitrary gauge shifts. Assess the actual
        # coupling marginals periodically instead, relative to each marginal.
        if iteration % 10 == 9 or iteration == max_iter - 1:
            plan = np.exp(logK + u[:, None] + v[None, :])
            rows_now = plan.sum(axis=1)
            cols_now = plan.sum(axis=0)
            row_resid = np.max(np.abs(rows_now - 1 / len(source))) / (1 / len(source))
            col_resid = np.max(np.abs(cols_now - 1 / len(target))) / (1 / len(target))
            if max(row_resid, col_resid) <= tol:
                converged = True
                break
    plan = np.exp(logK + u[:, None] + v[None, :])
    rows = plan.sum(axis=1, keepdims=True)
    if not converged or not np.isfinite(plan).all() or np.any(rows <= 0):
        cost_stats = np.percentile(cost, [0, 50, 90, 99, 100])
        raise RuntimeError(
            'Sinkhorn coupling failed to converge: '
            f'relative row residual={row_resid:.3g}, column residual={col_resid:.3g}, '
            f'iterations={iteration + 1}, epsilon_scale={epsilon * scale:.3g}, '
            f'cost_quantiles={np.array2string(cost_stats, precision=4)}'
        )
    return plan @ target / rows


def _flow_features(x, z):
    x = np.asarray(x, float)
    z = np.asarray(z, float).reshape(-1, 1)
    return np.column_stack([np.ones(len(x)), z[:, 0], z[:, 0] ** 2, x, z * x])


def fit_flow(clouds, centers, ridge=1., epsilon=.5):
    """Distill adjacent entropic-OT barycentric displacement into affine v(x,z)."""
    c = _centers(centers)
    if len(clouds) != len(c):
        raise ValueError('cloud count must equal center count')
    clouds = [_matrix(x, f'cloud[{i}]', 2) for i, x in enumerate(clouds)]
    d = clouds[0].shape[1]
    if any(x.shape[1] != d for x in clouds):
        raise ValueError('all clouds must have the same feature count')
    if ridge < 0 or not np.isfinite(ridge):
        raise ValueError('ridge must be finite and nonnegative')
    feats, targets = [], []
    for i, (a, b) in enumerate(zip(clouds[:-1], clouds[1:])):
        dz = c[i + 1] - c[i]
        bary = _sinkhorn_barycenters(a, b, epsilon)
        weight = 1 / np.sqrt(len(a))
        feats.append(weight * _flow_features(a, np.full(len(a), (c[i] + c[i + 1]) / 2)))
        targets.append(weight * (bary - a) / dz)
    F, T = np.vstack(feats), np.vstack(targets)
    penalty = np.eye(F.shape[1]) * ridge
    penalty[0, 0] = 0.
    coef = (np.linalg.solve(F.T @ F + penalty, F.T @ T) if ridge > 0
            else np.linalg.lstsq(F, T, rcond=None)[0])
    return {'coef': coef, 'n_features': d, 'ridge': float(ridge), 'epsilon': float(epsilon)}


def _velocity(model, x, z):
    return _flow_features(x, np.full(len(x), z)) @ model['coef']


def advect(model, Y, z0, z1, steps=8):
    """Integrate the fitted field with RK4 along scalar z endpoints."""
    Y = _matrix(Y, 'Y')
    if Y.shape[1] != model['n_features'] or not np.isfinite([z0, z1]).all() or steps < 1:
        raise ValueError('invalid dimension, endpoints, or steps')
    h = (float(z1) - float(z0)) / int(steps)
    x = Y.copy()
    for j in range(int(steps)):
        z = float(z0) + j * h
        k1 = _velocity(model, x, z)
        k2 = _velocity(model, x + h * k1 / 2, z + h / 2)
        k3 = _velocity(model, x + h * k2 / 2, z + h / 2)
        k4 = _velocity(model, x + h * k3, z + h)
        x += h * (k1 + 2*k2 + 2*k3 + k4) / 6
        if not np.isfinite(x).all():
            raise RuntimeError('Flow integration diverged')
    return x


def energy_distance(X, Y):
    """Multivariate energy V-statistic, clamped against roundoff below zero."""
    X, Y = _matrix(X, 'X'), _matrix(Y, 'Y')
    if X.shape[1] != Y.shape[1]:
        raise ValueError('X and Y must have the same feature count')
    from scipy.spatial.distance import cdist
    value = 2 * cdist(X, Y).mean() - cdist(X, X).mean() - cdist(Y, Y).mean()
    return float(max(0., value))


def benchmark(Ytrain, ztrain, Ytest, ztest, centers, halfwidth,
              min_count=12, max_points=128, seed=15, ridge=1., epsilon=.5, steps=8):
    """Cross-specimen one-step prediction benchmark using shared bounded clouds."""
    train = window_clouds(Ytrain, ztrain, centers, halfwidth, min_count, max_points, seed)
    test = window_clouds(Ytest, ztest, centers, halfwidth, min_count, max_points, seed + 1)
    c = _centers(centers)
    model = fit_flow(train, c, ridge, epsilon)
    # Local mean increments are smoothed over the training sequence.
    means = np.array([x.mean(axis=0) for x in train])
    ztrain_arr = np.asarray(ztrain, float).reshape(-1, 1)
    reg = KNeighborsRegressor(n_neighbors=min(64, len(ztrain_arr)), weights='distance').fit(ztrain_arr, np.asarray(Ytrain, float))
    rows = []
    for i in range(len(c) - 1):
        source, target = test[i], test[i + 1]
        predictions = {
            'identity': source,
            'centroid': source + (means[i + 1] - means[i]),
            'dpt_local': source + (reg.predict([[c[i + 1]]]).mean(axis=0) - reg.predict([[c[i]]]).mean(axis=0)),
            'flow': advect(model, source, c[i], c[i + 1], steps),
        }
        for method, prediction in predictions.items():
            rows.append({'transition': i, 'z_from': c[i], 'z_to': c[i + 1],
                         'method': method, 'energy_distance': energy_distance(prediction, target),
                         'n_source': len(source), 'n_target': len(target)})
    return pd.DataFrame(rows)


def _normal_logpdf(Y, mean, var):
    return -.5 * np.sum(np.log(2*np.pi*var) + (Y - mean) ** 2 / var, axis=1)


def fit_density_atlas(Y, z, labels, centers, halfwidth, flow_weight=0., anatomy_strength=.25,
                      endpoint_only=False, min_count=12, max_points=128, seed=15,
                      ridge=1., epsilon=.5, specimens=None):
    """Fit local diagonal Gaussian atlas and optional flow-smoothed moments."""
    Y = _matrix(Y, 'Y', 2)
    z = np.asarray(z, float).ravel()
    labels = np.asarray(labels).astype(str).ravel()
    c = _centers(centers)
    if len(z) != len(Y) or len(labels) != len(Y) or not np.isfinite(z).all():
        raise ValueError('z and labels must match Y rows; z must be finite')
    if not 0 <= flow_weight <= 1 or anatomy_strength < 0:
        raise ValueError('flow_weight must be in [0,1] and anatomy_strength nonnegative')
    clouds = window_clouds(Y, z, c, halfwidth, min_count, max_points, seed, specimens)
    means = np.array([x.mean(axis=0) for x in clouds])
    variances = np.array([x.var(axis=0, ddof=1) for x in clouds])
    floor = max(float(np.median(variances[variances > 0])) * 1e-3 if np.any(variances > 0) else 1e-6, 1e-6)
    variances = np.maximum(variances, floor)
    model = None
    if flow_weight > 0:
        model = fit_flow(clouds, c, ridge, epsilon)
        flowed_means, flowed_vars = [means[0]], [variances[0]]
        for i in range(len(c)-1):
            transported = advect(model, clouds[i], c[i], c[i+1])
            flowed_means.append(transported.mean(axis=0))
            flowed_vars.append(np.maximum(transported.var(axis=0, ddof=1), floor))
        raw_means, raw_vars = means.copy(), variances.copy()
        means = (1-flow_weight)*raw_means + flow_weight*np.asarray(flowed_means)
        variances = ((1-flow_weight) * (raw_vars + (raw_means - means)**2)
                     + flow_weight * (np.asarray(flowed_vars) + (np.asarray(flowed_means) - means)**2))
    # Ordered cutpoints in [0,1], estimated from available labels and positions.
    use = np.isin(labels, ['S1', 'S2', 'S3'])
    if endpoint_only:
        use &= np.isin(labels, ['S1', 'S3'])
    order = {'S1': 0, 'S2': 1, 'S3': 2}
    ranks = np.array([order.get(x, -1) for x in labels])
    cutpoints = np.array([np.nan, np.nan], float)
    tau = .1
    if anatomy_strength > 0:
        if use.sum() < 4 or not np.any(ranks[use] == 0) or not np.any(ranks[use] == 2):
            raise ValueError('anatomy supervision requires at least four labeled rows and both S1 and S3')
        # Parameterize c2=c1+positive gap; mild quadratic regularization prevents extremes.
        def unpack(p):
            a = expit(p[0]); b = a + (1-a)*expit(p[1]); return a, b
        def loss(p):
            a, b = unpack(p)
            zz, rr = z[use], ranks[use]
            p1, p2 = expit((zz-a)/tau), expit((zz-b)/tau)
            probs = np.column_stack([1-p1, p1-p2, p2])
            nll = -np.log(np.maximum(probs[np.arange(len(rr)), rr], 1e-10)).sum()
            return nll + 1e-3*np.sum(p*p)
        res = minimize(loss, np.array([-0.7, 0.]), method='L-BFGS-B')
        if not res.success or not np.isfinite(res.fun):
            raise RuntimeError('ordered anatomy cutpoint fit failed')
        cutpoints[:] = unpack(res.x)
    return {'centers': c, 'means': means, 'variances': variances, 'clouds': clouds,
            'flow_model': model, 'cutpoints': cutpoints, 'tau': tau,
            'anatomy_strength': float(anatomy_strength), 'endpoint_only': bool(endpoint_only),
            'feature_count': Y.shape[1]}


def posterior(atlas, Y, labels=None):
    """Return posterior mean coordinate, grid probabilities, and entropy."""
    Y = _matrix(Y, 'Y')
    if Y.shape[1] != atlas['feature_count']:
        raise ValueError('Y feature count does not match atlas')
    c, means, vars_ = atlas['centers'], atlas['means'], atlas['variances']
    logp = np.column_stack([_normal_logpdf(Y, means[j], vars_[j]) for j in range(len(c))])
    logp -= np.log(len(c))  # uniform grid prior
    if labels is not None and atlas['anatomy_strength'] > 0:
        labels = np.asarray(labels).astype(str).ravel()
        if len(labels) != len(Y):
            raise ValueError('labels must match Y rows')
        a, b = atlas['cutpoints']; tau = atlas['tau']
        for i, lab in enumerate(labels):
            if lab not in {'S1', 'S2', 'S3'} or (atlas['endpoint_only'] and lab == 'S2'):
                continue
            p1, p2 = expit((c-a)/tau), expit((c-b)/tau)
            probs = {'S1': 1-p1, 'S2': p1-p2, 'S3': p2}[lab]
            logp[i] += atlas['anatomy_strength'] * np.log(np.maximum(probs, 1e-10))
    logp -= logsumexp(logp, axis=1, keepdims=True)
    probs = np.exp(logp)
    mean_z = probs @ c
    entropy = -np.sum(probs * logp, axis=1)
    return mean_z, probs, entropy


def latent_fit(Y, z0, labels, centers, halfwidth, flow_weight=0., iterations=4,
               endpoint_only=False, anatomy_strength=.25, min_count=12,
               max_points=128, seed=15, ridge=1., epsilon=.5, specimens=None):
    """Alternating atlas/posterior pilot initialized from supplied positions."""
    Y = _matrix(Y, 'Y', 2)
    z = np.asarray(z0, float).ravel().copy()
    if len(z) != len(Y) or not np.isfinite(z).all() or iterations < 1:
        raise ValueError('z0 must be finite and match Y; iterations must be positive')
    atlas = None
    for _ in range(int(iterations)):
        atlas = fit_density_atlas(Y, z, labels, centers, halfwidth, flow_weight,
                                  anatomy_strength, endpoint_only, min_count,
                                  max_points, seed, ridge, epsilon, specimens)
        z, _, _ = posterior(atlas, Y, labels)
    atlas = fit_density_atlas(Y, z, labels, centers, halfwidth, flow_weight,
                              anatomy_strength, endpoint_only, min_count,
                              max_points, seed, ridge, epsilon, specimens)
    z, probs, entropy = posterior(atlas, Y, labels)
    return {'z': z, 'posterior': probs, 'entropy': entropy, 'atlas': atlas}
