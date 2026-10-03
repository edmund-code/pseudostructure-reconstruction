"""Privileged within-segment diagnostics for coordinate validation."""
from __future__ import annotations

import numpy as np


def _labels(values, n, name):
    labels = np.asarray(values, dtype=object)
    if labels.ndim != 1 or len(labels) != n:
        raise ValueError(f'{name} must be a one-dimensional vector matching rows')
    if any(x is None or (isinstance(x, (float, np.floating)) and np.isnan(x))
           or not str(x).strip() for x in labels):
        raise ValueError(f'{name} must contain nonmissing labels')
    return np.asarray([str(x) for x in labels], dtype=str)


def _matrix(X, name, min_rows=1):
    X = np.asarray(X, dtype=float)
    if X.ndim != 2 or X.shape[0] < min_rows or X.shape[1] < 1:
        raise ValueError(f'{name} must be a 2D matrix with at least {min_rows} rows and one column')
    if not np.isfinite(X).all():
        raise ValueError(f'{name} must contain only finite values')
    return X


def _anatomy(anatomy, n):
    raw = np.asarray(anatomy)
    if raw.ndim != 1 or len(raw) != n:
        raise ValueError('anatomy must be a one-dimensional vector matching rows')
    try:
        numeric = raw.astype(float)
    except (TypeError, ValueError) as exc:
        raise ValueError('anatomy labels must be integers in {0,1,2}') from exc
    if not np.isfinite(numeric).all() or not np.isin(numeric, [0, 1, 2]).all():
        raise ValueError('anatomy labels must be integers in {0,1,2}')
    return numeric.astype(int)


def permute_positions_within_segments(z, specimens, anatomy, seed):
    """Shuffle coordinates independently within specimen × segment strata."""
    z = np.asarray(z, dtype=float)
    if z.ndim != 1 or not np.isfinite(z).all() or np.any((z < 0) | (z > 1)):
        raise ValueError('z must be a finite one-dimensional vector in [0,1]')
    specimens = _labels(specimens, len(z), 'specimens')
    anatomy = _anatomy(anatomy, len(z))
    if not isinstance(seed, (int, np.integer)) or isinstance(seed, (bool, np.bool_)):
        raise ValueError('seed must be an integer')
    result = z.copy()
    rng = np.random.default_rng(seed)
    groups = {}
    for i, key in enumerate(zip(specimens, anatomy)):
        groups.setdefault(key, []).append(i)
    for rows in groups.values():
        rows = np.asarray(rows, dtype=int)
        result[rows] = rng.permutation(z[rows])
    return result


def fit_segment_covariate_oracle(expression, covariates, anatomy, specimens,
                                 ridge=.01, degree=1):
    """Fit specimen-balanced ridge predictions separately within segments.

    Covariates are centered and scaled within each segment using equal total
    mass per original specimen. The gene-response intercept is the weighted
    segment mean and is not penalized. This is a privileged diagnostic oracle:
    prediction requires the query's supplied segment label. Degree two adds
    squares and pairwise interactions, with expanded-basis scaling frozen
    from the original-specimen-weighted training observations.
    """
    expression = _matrix(expression, 'expression', min_rows=2)
    covariates = _matrix(covariates, 'covariates', min_rows=2)
    if len(covariates) != len(expression):
        raise ValueError('covariates rows must match expression')
    anatomy = _anatomy(anatomy, len(expression))
    specimens = _labels(specimens, len(expression), 'specimens')
    if not np.isscalar(ridge):
        raise ValueError('ridge must be finite and nonnegative')
    try:
        ridge = float(ridge)
    except (TypeError, ValueError) as exc:
        raise ValueError('ridge must be finite and nonnegative') from exc
    if not np.isfinite(ridge) or ridge < 0:
        raise ValueError('ridge must be finite and nonnegative')
    if (not isinstance(degree, (int, np.integer)) or isinstance(degree, (bool, np.bool_))
            or degree not in (1, 2)):
        raise ValueError('degree must be integer 1 or 2')

    segments = {}
    for segment in (0, 1, 2):
        use = anatomy == segment
        if not use.any():
            raise ValueError(f'training data lack segment {segment}')
        segment_specimens = np.unique(specimens[use])
        local_weights = np.zeros(use.sum(), dtype=float)
        local_specimens = specimens[use]
        for specimen in segment_specimens:
            local = local_specimens == specimen
            local_weights[local] = 1. / (len(segment_specimens) * local.sum())
        X, R = covariates[use], expression[use]
        center = local_weights @ X
        scale = np.sqrt(local_weights @ ((X - center) ** 2))
        scale[scale < 1e-8] = 1.
        standardized = (X - center) / scale
        intercept = local_weights @ R
        attrs = {'covariate_mean': center, 'covariate_scale': scale,
                 'intercept': intercept, 'training_specimens': segment_specimens}
        if degree == 1:
            design = standardized
        else:
            poly_rows, poly_cols = np.triu_indices(X.shape[1])
            polynomial = np.column_stack(
                [standardized[:, i] * standardized[:, j]
                 for i, j in zip(poly_rows, poly_cols)])
            basis = np.column_stack([standardized, polynomial])
            basis_mean = local_weights @ basis
            basis_scale = np.sqrt(local_weights @ ((basis - basis_mean) ** 2))
            basis_scale[basis_scale < 1e-8] = 1.
            design = (basis - basis_mean) / basis_scale
            attrs.update(polynomial_rows=poly_rows, polynomial_cols=poly_cols,
                         basis_mean=basis_mean, basis_scale=basis_scale)
        lhs = design.T @ (local_weights[:, None] * design)
        rhs = design.T @ (local_weights[:, None] * (R - intercept))
        beta = np.linalg.solve(lhs + ridge*np.eye(design.shape[1])
                               + 1e-10*np.eye(design.shape[1]), rhs)
        attrs['beta'] = beta
        segments[segment] = attrs
    return {'segments': segments, 'covariate_count': covariates.shape[1],
            'expression_count': expression.shape[1], 'ridge': ridge,
            'degree': int(degree)}


def predict_segment_covariate_oracle(model, covariates, anatomy):
    """Predict expression from covariates and explicitly supplied segment IDs."""
    covariates = _matrix(covariates, 'covariates')
    if covariates.shape[1] != model['covariate_count']:
        raise ValueError('covariate feature count does not match fitted oracle')
    anatomy = _anatomy(anatomy, len(covariates))
    prediction = np.empty((len(covariates), model['expression_count']), dtype=float)
    degree = model.get('degree', 1)
    for segment in (0, 1, 2):
        use = anatomy == segment
        if not use.any():
            continue
        params = model['segments'][segment]
        scaled = (covariates[use] - params['covariate_mean']) / params['covariate_scale']
        if degree == 1:
            design = scaled
        elif degree == 2:
            polynomial = np.column_stack(
                [scaled[:, i] * scaled[:, j]
                 for i, j in zip(params['polynomial_rows'], params['polynomial_cols'])])
            basis = np.column_stack([scaled, polynomial])
            design = (basis - params['basis_mean']) / params['basis_scale']
        else:
            raise ValueError('fitted oracle degree must be 1 or 2')
        prediction[use] = params['intercept'] + design @ params['beta']
    return prediction
