"""Low-complexity reference baselines for repeated-specimen atlases."""
from __future__ import annotations

import warnings

import numpy as np
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression


def _matrix(X, name, min_rows=2):
    X = np.asarray(X, dtype=float)
    if X.ndim != 2 or X.shape[0] < min_rows or X.shape[1] < 1:
        raise ValueError(f'{name} must be a 2D matrix with at least {min_rows} rows and one feature')
    if not np.isfinite(X).all():
        raise ValueError(f'{name} must contain only finite values')
    return X


def _anatomy(anatomy, n, require_all=False):
    labels = np.asarray(anatomy)
    if labels.ndim != 1 or len(labels) != n:
        raise ValueError('anatomy must be a one-dimensional vector matching rows')
    try:
        numeric = labels.astype(float)
    except (TypeError, ValueError) as exc:
        raise ValueError('anatomy labels must be integers in {0,1,2}') from exc
    if not np.isfinite(numeric).all() or not np.isin(numeric, [0, 1, 2]).all():
        raise ValueError('anatomy labels must be integers in {0,1,2}')
    labels = numeric.astype(int)
    if require_all and not np.array_equal(np.unique(labels), [0, 1, 2]):
        raise ValueError('training anatomy must contain all labels 0, 1, and 2')
    if not np.any(labels == 0) or not np.any(labels == 2):
        raise ValueError('training anatomy must contain endpoint labels 0 and 2')
    return labels


def fit_pc1_reference(Y, anatomy):
    """Orient and robustly scale the existing first balanced-PC coordinate."""
    Y = _matrix(Y, 'Y')
    anatomy = _anatomy(anatomy, len(Y))
    scores = Y[:, 0].copy()
    if np.ptp(scores) <= np.finfo(float).eps:
        raise ValueError('training first coordinate has no variation')
    sign = 1.
    if np.median(scores[anatomy == 2]) < np.median(scores[anatomy == 0]):
        sign = -1.
        scores *= sign
    low, high = np.quantile(scores, [.01, .99])
    if not np.isfinite([low, high]).all() or high <= low:
        raise ValueError('training PC1 quantile range is constant')
    return {'sign': sign, 'low': float(low), 'high': float(high),
            'feature_count': Y.shape[1]}


def project_pc1_reference(reference, Y):
    """Project rows to the frozen training PC1 scale, clipped to [0,1]."""
    Y = _matrix(Y, 'Y', min_rows=1)
    if Y.shape[1] != reference['feature_count']:
        raise ValueError('Y feature count does not match PC1 reference')
    score = Y[:, 0] * reference['sign']
    return np.clip((score - reference['low']) / (reference['high'] - reference['low']), 0., 1.)


def fit_segment_mixture(Y, expression, anatomy, specimens, C=1., seed=15):
    """Fit a specimen-balanced soft segment classifier and canonical means.

    Specimen IDs affect training weights and within-segment mean averaging;
    they are never classifier features. Query projection requires only Y.
    """
    Y = _matrix(Y, 'Y')
    expression = _matrix(expression, 'expression')
    if len(expression) != len(Y):
        raise ValueError('expression rows must match Y')
    labels = _anatomy(anatomy, len(Y), require_all=True)
    specimens = np.asarray(specimens, dtype=object).ravel()
    if len(specimens) != len(Y) or any(
            x is None or (isinstance(x, (float, np.floating)) and np.isnan(x))
            or not str(x).strip() for x in specimens):
        raise ValueError('specimens must contain one nonempty label per row')
    specimens = np.asarray([str(x) for x in specimens])
    if not np.isscalar(C):
        raise ValueError('C must be positive and finite')
    try:
        C = float(C)
    except (TypeError, ValueError) as exc:
        raise ValueError('C must be positive and finite') from exc
    if not np.isfinite(C) or C <= 0:
        raise ValueError('C must be positive and finite')
    if (not isinstance(seed, (int, np.integer)) or isinstance(seed, (bool, np.bool_))):
        raise ValueError('seed must be an integer')

    unique, inverse, counts = np.unique(specimens, return_inverse=True, return_counts=True)
    weights = 1. / (len(unique) * counts[inverse])
    mean = weights @ Y
    scale = np.sqrt(weights @ ((Y - mean) ** 2))
    scale[scale < 1e-8] = 1.
    standardized = (Y - mean) / scale
    classifier = LogisticRegression(C=C, solver='lbfgs', max_iter=1000,
                                    random_state=int(seed))
    with warnings.catch_warnings():
        warnings.simplefilter('error', ConvergenceWarning)
        try:
            classifier.fit(standardized, labels, sample_weight=weights * len(Y))
        except ConvergenceWarning as exc:
            raise RuntimeError('segment mixture classifier did not converge') from exc

    gene_means = np.empty((3, expression.shape[1]), dtype=float)
    for segment in range(3):
        sample_means = []
        for specimen in unique:
            use = (specimens == specimen) & (labels == segment)
            if np.any(use):
                sample_means.append(expression[use].mean(axis=0))
        if not sample_means:
            raise ValueError(f'training anatomy lacks segment {segment}')
        gene_means[segment] = np.mean(sample_means, axis=0)
    return {'classifier': classifier, 'mean': mean, 'scale': scale,
            'gene_means': gene_means, 'feature_count': Y.shape[1],
            'expression_count': expression.shape[1], 'class_order': np.arange(3),
            'training_specimens': unique}


def predict_segment_mixture(model, Y):
    """Predict canonical molecular mean using expression-only soft weights."""
    Y = _matrix(Y, 'Y', min_rows=1)
    if Y.shape[1] != model['feature_count']:
        raise ValueError('Y feature count does not match segment mixture')
    classifier = model['classifier']
    raw = classifier.predict_proba((Y - model['mean']) / model['scale'])
    probabilities = np.zeros((len(Y), 3), dtype=float)
    probabilities[:, classifier.classes_.astype(int)] = raw
    prediction = probabilities @ model['gene_means']
    return {'prediction': prediction, 'probabilities': probabilities}
