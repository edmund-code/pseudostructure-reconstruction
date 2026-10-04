"""Joint repeated-specimen canonical atlas with latent positions.

Profiles are modeled as noisy cross-sections around a shared smooth mean,
with small constant specimen offsets and a latent finite-grid coordinate.
Nephron identities are marginalized: fitted posterior assignments do not
track a physical tubule. Projection uses expression alone and is exploratory,
with a fixed gauge, finite-grid uniform prior, and uncalibrated posterior.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit, logsumexp
from sklearn.utils.extmath import randomized_svd

from pseudospace.repeated_structure import _basis


def _matrix(X, name='X', min_rows=1):
    X = np.asarray(X, dtype=float)
    if X.ndim != 2 or X.shape[0] < min_rows or X.shape[1] < 1:
        raise ValueError(f'{name} must be a 2D matrix with at least {min_rows} rows and one feature')
    if not np.isfinite(X).all():
        raise ValueError(f'{name} must contain only finite values')
    return X


def _specimen_weights(specimens, balance=True):
    labels, inverse, counts = np.unique(specimens, return_inverse=True, return_counts=True)
    if balance:
        weights = 1. / (len(labels) * counts[inverse])
    else:
        weights = np.full(len(specimens), 1. / len(specimens))
    return labels, inverse, weights


def fit_balanced_representation(X, specimens, n_components=15, seed=15):
    """Fit equal-specimen weighted standardization and randomized-SVD PCs."""
    X = _matrix(X, 'X', min_rows=2)
    specimens = np.asarray(specimens, dtype=object).ravel()
    if len(specimens) != len(X) or any(x is None or (isinstance(x, (float, np.floating)) and np.isnan(x))
                                         or not str(x).strip() for x in specimens):
        raise ValueError('specimens must have one nonempty label per row')
    if (not isinstance(n_components, (int, np.integer)) or isinstance(n_components, (bool, np.bool_))
            or n_components < 1):
        raise ValueError('n_components must be a positive integer')
    specimens = np.asarray([str(x) for x in specimens])
    _, _, weights = _specimen_weights(specimens)
    mean = weights @ X
    scale = np.sqrt(weights @ ((X - mean) ** 2))
    scale[scale < 1e-8] = 1.
    standardized = (X - mean) / scale
    weighted = standardized * np.sqrt(weights * len(X))[:, None]
    k = min(int(n_components), len(X) - 1, X.shape[1])
    if k < 1:
        raise ValueError('n_components must be positive')
    _, _, components = randomized_svd(weighted, n_components=k, random_state=seed)
    transform = {'mean': mean, 'scale': scale, 'components': components}
    return standardized @ components.T, transform


def transform_balanced_representation(X, transform):
    """Apply a frozen balanced representation to new rows."""
    X = _matrix(X, 'X')
    mean = np.asarray(transform['mean'], float)
    scale = np.asarray(transform['scale'], float)
    components = np.asarray(transform['components'], float)
    if mean.ndim != 1 or scale.shape != mean.shape or X.shape[1] != len(mean):
        raise ValueError('X feature count does not match representation')
    if components.ndim != 2 or components.shape[1] != len(mean) or not np.isfinite(components).all():
        raise ValueError('representation components have invalid dimensions or values')
    if not np.isfinite(mean).all() or not np.isfinite(scale).all() or np.any(scale <= 0):
        raise ValueError('representation mean/scale must be finite with positive scales')
    return ((X - mean) / scale) @ components.T


def _ordinal_probs(grid, cutpoints, width):
    c1, c2 = cutpoints
    p1 = expit((grid - c1) / width)
    p2 = expit((grid - c2) / width)
    return np.column_stack([1 - p1, p1 - p2, p2])


def _cutpoint_fit(posterior, labels, weights, grid, width, initial):
    observed = labels >= 0
    if not observed.any():
        return np.asarray(initial, float)

    def unpack(p):
        c1 = .03 + .89 * expit(p[0])
        c2 = c1 + .05 + (.97 - c1 - .05) * expit(p[1])
        return np.array([c1, c2])

    def objective(p):
        cuts = unpack(p)
        probs = _ordinal_probs(grid, cuts, width)
        log_label = np.log(np.maximum(probs[:, labels[observed]].T, 1e-12))
        expected = np.sum(posterior[observed] * log_label, axis=1)
        return -np.dot(weights[observed], expected)

    def logit(q):
        q = np.clip(q, 1e-6, 1-1e-6)
        return np.log(q / (1-q))
    c1, c2 = initial
    p0 = [logit((c1-.03)/.89), logit((c2-c1-.05)/(.97-c1-.05))]
    result = minimize(objective, p0, method='L-BFGS-B')
    if not result.success or not np.isfinite(result.fun) or not np.isfinite(result.x).all():
        raise RuntimeError(f'ordered anatomy cutpoint optimization failed: {result.message}')
    return unpack(result.x)


def _posterior(Y, mean, offsets, specimens_idx, sigma2, grid, labels,
               anatomy_strength, cutpoints, anatomy_width):
    residual = Y[:, None, :] - mean[None, :, :] - offsets[specimens_idx, None, :]
    logp = -np.sum(residual**2, axis=2) / (2 * sigma2)
    if labels is not None and anatomy_strength > 0:
        probs = np.maximum(_ordinal_probs(grid, cutpoints, anatomy_width), 1e-12)
        valid = labels >= 0
        logp[valid] += anatomy_strength * np.log(probs[:, labels[valid]].T)
    logp -= np.log(len(grid))
    logp -= logsumexp(logp, axis=1, keepdims=True)
    posterior = np.exp(np.maximum(logp, np.log(np.finfo(float).tiny)))
    posterior /= posterior.sum(axis=1, keepdims=True)
    return posterior, residual


def fit_joint_atlas(Y, specimens, anatomy, initial_z, grid=None, roughness=.001,
                    offset_penalty=4., anatomy_strength=1., anatomy_width=.1,
                    max_iter=60, tol=1e-5, fit_offsets=True,
                    balance_specimens=True, allow_unlabeled_anatomy=False,
                    prior_scale=1.0):
    """Fit a shared atlas and latent finite-grid positions.

    ``allow_unlabeled_anatomy=True`` permits ``-1`` labels for unlabeled
    specimen point clouds. Those rows receive no ordinal anatomy information;
    this does not replace them with pseudo-S2 labels.

    ``prior_scale`` divides all feature values by one common scalar during
    fitting, then restores returned atlas values to input units. The quadratic
    smoothness and offset penalties therefore act on beta / prior_scale and
    offsets / prior_scale. For fixed Y, changing this scale changes prior strength. Jointly scaling
    Y and prior_scale preserves the fit. The scalar does not reweight features
    or alter upstream PCA/graph inputs. Objective histories use reference
    units; all returned molecular parameters and variances use input units.
    """
    Y = _matrix(Y, 'Y', min_rows=4)
    if not np.isscalar(prior_scale) or isinstance(prior_scale, (bool, np.bool_)):
        raise ValueError('prior_scale must be a positive finite scalar')
    try:
        prior_scale = float(prior_scale)
    except (TypeError, ValueError) as exc:
        raise ValueError('prior_scale must be a positive finite scalar') from exc
    if not np.isfinite(prior_scale) or prior_scale <= 0:
        raise ValueError('prior_scale must be a positive finite scalar')
    with np.errstate(over='ignore', under='ignore', invalid='ignore'):
        scale2 = float(prior_scale * prior_scale)
    if not np.isfinite(scale2) or scale2 <= 0:
        raise ValueError('prior_scale squared must be positive and finite')
    with np.errstate(over='ignore', divide='ignore', invalid='ignore'):
        Y = Y / prior_scale
    if not np.isfinite(Y).all():
        raise ValueError('prior_scale normalization overflowed')
    specimens = np.asarray(specimens, dtype=object).ravel()
    initial_z = np.asarray(initial_z, float).ravel()
    anatomy = np.asarray(anatomy)
    if len(specimens) != len(Y) or any(x is None or (isinstance(x, (float, np.floating)) and np.isnan(x))
                                       or not str(x).strip() for x in specimens):
        raise ValueError('specimens must have one nonempty label per row')
    specimens = np.asarray([str(x) for x in specimens])
    if len(initial_z) != len(Y) or not np.isfinite(initial_z).all() or np.any((initial_z < 0) | (initial_z > 1)):
        raise ValueError('initial_z must be finite, match rows, and lie in [0,1]')
    if anatomy.ndim != 1 or len(anatomy) != len(Y):
        raise ValueError('anatomy must have one label per row')
    if not isinstance(allow_unlabeled_anatomy, (bool, np.bool_)):
        raise ValueError('allow_unlabeled_anatomy must be boolean')
    try:
        numeric_anatomy = anatomy.astype(float)
    except (TypeError, ValueError) as exc:
        raise ValueError('anatomy labels must be integers in {0,1,2}') from exc
    allowed_anatomy = [0, 1, 2, -1] if allow_unlabeled_anatomy else [0, 1, 2]
    if (not np.isfinite(numeric_anatomy).all()
            or not np.equal(numeric_anatomy, np.rint(numeric_anatomy)).all()
            or not np.isin(numeric_anatomy, allowed_anatomy).all()):
        if allow_unlabeled_anatomy:
            raise ValueError('anatomy labels must be integers in {-1,0,1,2}')
        raise ValueError('anatomy labels must be integers in {0,1,2}')
    labels = numeric_anatomy.astype(int)
    if grid is None:
        grid = np.linspace(0., 1., 101)
    else:
        grid = np.asarray(grid, float).ravel()
    if len(grid) < 8 or not np.isfinite(grid).all() or np.any(np.diff(grid) <= 0) or grid[0] != 0 or grid[-1] != 1:
        raise ValueError('grid must be strictly increasing, contain at least 8 points, and span [0,1]')
    if any(not np.isscalar(value) for value in
           (roughness, offset_penalty, anatomy_strength, anatomy_width, tol)):
        raise ValueError('roughness, penalties, anatomy_width, and tol must be scalar values')
    scalar_params = np.asarray([roughness, offset_penalty, anatomy_strength,
                                anatomy_width, tol], dtype=float)
    if not np.isfinite(scalar_params).all() or roughness < 0 or offset_penalty < 0 or anatomy_strength < 0 or anatomy_width <= 0:
        raise ValueError('penalties/strength must be finite and nonnegative; anatomy_width must be finite and positive')
    if (not isinstance(max_iter, (int, np.integer)) or isinstance(max_iter, (bool, np.bool_))
            or max_iter < 1 or tol <= 0):
        raise ValueError('max_iter must be a positive integer and tol must be finite and positive')

    specimens_unique, specimen_idx, weights = _specimen_weights(specimens, balance_specimens)
    n, d = Y.shape
    B = _basis(grid, grid)
    q = B.shape[1]
    B_init = _basis(initial_z, grid)
    d2 = np.diff(np.eye(q), n=2, axis=0)
    roughness_matrix = d2.T @ d2
    penalty = roughness * roughness_matrix
    weighted_y_var = float(np.sum(weights[:, None] * (Y - weights @ Y)**2) / d)
    beta = np.linalg.solve(B_init.T @ (weights[:, None] * B_init)
                           + roughness * max(weighted_y_var, 1e-8) * roughness_matrix
                           + 1e-8*np.eye(q),
                           B_init.T @ (weights[:, None] * Y))
    offsets = np.zeros((len(specimens_unique), d))
    mean = B @ beta
    residual = Y - B_init @ beta
    sigma2 = max(float(np.sum(weights[:, None] * residual**2) / d), 1e-6)
    cutpoints = np.array([.33, .67])
    history = []
    converged = False
    posterior, _ = _posterior(Y, mean, offsets, specimen_idx, sigma2, grid,
                              labels, anatomy_strength, cutpoints, anatomy_width)

    for iteration in range(int(max_iter)):
        if anatomy_strength > 0:
            cutpoints = _cutpoint_fit(posterior, labels, weights, grid,
                                      anatomy_width, cutpoints)
        # Two coordinate sweeps reduce the coupling between mean and offsets.
        for _ in range(2):
            weighted_grid = posterior.T @ weights
            rhs = B.T @ (posterior.T @ (weights[:, None] * (Y - offsets[specimen_idx])))
            lhs = B.T @ (weighted_grid[:, None] * B) + sigma2*penalty + 1e-8*np.eye(q)
            beta = np.linalg.solve(lhs, rhs)
            mean = B @ beta
            if fit_offsets:
                for s in range(len(specimens_unique)):
                    use = specimen_idx == s
                    group_mass = weights[use].sum()
                    expected_mean = posterior[use] @ mean
                    residual_mean = np.sum(weights[use, None] * (Y[use] - expected_mean), axis=0) / group_mass
                    shrink = 1 + offset_penalty * sigma2 / (len(specimens_unique) * group_mass)
                    offsets[s] = residual_mean / shrink
                # Fixed gauge: equal-specimen average offset is exactly zero.
                gauge = offsets.mean(axis=0)
                offsets -= gauge
                beta += gauge[None, :]  # cubic basis partitions unity
                mean = B @ beta

        # The current E-step responsibilities remain fixed through the M-step,
        # including the variance update.
        residual = Y[:, None, :] - mean[None, :, :] - offsets[specimen_idx, None, :]
        expected_sse = np.sum(weights[:, None, None] * posterior[:, :, None] * residual**2)
        sigma2 = max(float(expected_sse / d), 1e-6)
        posterior, _ = _posterior(Y, mean, offsets, specimen_idx, sigma2,
                                  grid, labels, anatomy_strength,
                                  cutpoints, anatomy_width)
        # Observed weighted objective includes weak labels and explicit penalties.
        residual = Y[:, None, :] - mean[None, :, :] - offsets[specimen_idx, None, :]
        logp = (-np.sum(residual**2, axis=2) / (2*sigma2)
                - .5*d*np.log(2*np.pi*sigma2))
        if anatomy_strength > 0:
            lp = np.log(np.maximum(_ordinal_probs(grid, cutpoints, anatomy_width), 1e-12))
            observed = labels >= 0
            logp[observed] += anatomy_strength * lp[:, labels[observed]].T
        logp -= np.log(len(grid))
        objective = -float(np.dot(weights, logsumexp(logp, axis=1)))
        objective += .5*roughness*np.sum((d2 @ beta)**2)
        objective += .5*offset_penalty*np.mean(np.sum(offsets**2, axis=1))
        objective_change = (np.nan if not history else objective - history[-1]['objective'])
        history.append({'iteration': iteration + 1, 'objective': objective,
                        'sigma2': sigma2 * scale2,
                        'objective_change': objective_change})
        if len(history) > 1 and abs(objective_change) <= tol * max(1., abs(history[-2]['objective'])):
            converged = True
            break

    posterior, residual = _posterior(Y, mean, offsets, specimen_idx, sigma2,
                                      grid, labels, anatomy_strength,
                                      cutpoints, anatomy_width)
    z_mean = posterior @ grid
    z_map = grid[np.argmax(posterior, axis=1)]
    offsets_by_specimen = {sample: offsets[i].copy()
                           for i, sample in enumerate(specimens_unique)}
    mean = mean * prior_scale
    beta = beta * prior_scale
    offsets_by_specimen = {sample: value * prior_scale
                           for sample, value in offsets_by_specimen.items()}
    sigma2 = sigma2 * scale2
    if (not np.isfinite(mean).all() or not np.isfinite(beta).all()
            or not np.isfinite(sigma2) or sigma2 <= 0
            or any(not np.isfinite(value).all() for value in offsets_by_specimen.values())
            or any(not np.isfinite(item['sigma2']) or item['sigma2'] <= 0 for item in history)):
        raise ValueError('prior_scale restoration produced nonfinite parameters or nonpositive variance')
    return {'grid': grid, 'mean': mean, 'beta': beta, 'offsets': offsets_by_specimen,
            'specimens': specimens_unique, 'sigma2': sigma2,
            'cutpoints': cutpoints, 'posterior': posterior,
            'z_mean': z_mean, 'z_MAP': z_map, 'history': history,
            'converged': converged,
            'weights': weights, 'feature_count': d,
            'anatomy_strength': anatomy_strength, 'anatomy_width': anatomy_width,
            'prior_scale': prior_scale}


def project_joint_atlas(atlas, Y, offset=None, block_size=512):
    """Project expression alone onto a frozen isotropic Gaussian mean atlas."""
    Y = _matrix(Y, 'Y', min_rows=1)
    if Y.shape[1] != atlas['feature_count']:
        raise ValueError('Y feature count does not match atlas')
    if not isinstance(block_size, (int, np.integer)) or block_size < 1:
        raise ValueError('block_size must be a positive integer')
    if offset is None:
        offset = np.zeros(Y.shape[1])
    offset = np.asarray(offset, float).ravel()
    if len(offset) != Y.shape[1] or not np.isfinite(offset).all():
        raise ValueError('offset must be a finite vector matching the feature count')
    grid, mean, sigma2 = atlas['grid'], atlas['mean'], atlas['sigma2']
    logp = np.empty((len(Y), len(grid)))
    for start in range(0, len(Y), int(block_size)):
        stop = min(start + int(block_size), len(Y))
        residual = Y[start:stop, None, :] - mean[None, :, :] - offset
        logp[start:stop] = -np.sum(residual**2, axis=2) / (2*sigma2)
    logp -= np.log(len(grid))
    normalizer = logsumexp(logp, axis=1, keepdims=True)
    log_post = logp - normalizer
    posterior = np.exp(np.maximum(log_post, np.log(np.finfo(float).tiny)))
    posterior /= posterior.sum(axis=1, keepdims=True)
    z_map = grid[np.argmax(logp, axis=1)]
    z_mean = posterior @ grid
    entropy = -np.sum(posterior * log_post, axis=1)
    expected_mean = posterior @ mean + offset
    residual_norm = np.linalg.norm(Y - expected_mean, axis=1)
    return {'posterior': posterior, 'z_mean': z_mean, 'z_MAP': z_map,
            'entropy': entropy, 'residual_norm': residual_norm}
