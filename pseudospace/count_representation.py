"""Frozen Pearson-residual representation for count-based molecular inputs.

This adapts an observation-depth residualization to aggregate count profiles.
It is a representation experiment, not a depth-corrected spatial coordinate.
Genes with zero training probability carry no information and are assigned
zero residuals in every transformed query, including nonzero query counts.
"""
from __future__ import annotations

import numpy as np
from sklearn.utils.extmath import randomized_svd


def _counts_matrix(counts, name, min_rows=2):
    X = np.asarray(counts, dtype=float)
    if X.ndim != 2 or X.shape[0] < min_rows or X.shape[1] < 1:
        raise ValueError(f'{name} must be a 2D count matrix with at least {min_rows} rows and one gene')
    if not np.isfinite(X).all() or np.any(X < 0):
        raise ValueError(f'{name} must contain finite nonnegative counts')
    if np.any(np.abs(X - np.rint(X)) > 1e-6):
        raise ValueError(f'{name} must contain integer-valued counts')
    return X


def _labels(specimens, n):
    labels = np.asarray(specimens, dtype=object).ravel()
    if len(labels) != n or any(
            x is None or (isinstance(x, (float, np.floating)) and np.isnan(x))
            or not str(x).strip() for x in labels):
        raise ValueError('specimens must contain one nonempty label per count row')
    return np.asarray([str(x) for x in labels], dtype=str)


def _exposure(library, counts, name):
    library = np.asarray(library, dtype=float).ravel()
    if len(library) != len(counts) or not np.isfinite(library).all() or np.any(library <= 0):
        raise ValueError(f'{name} must have one positive finite exposure per row')
    if np.any(library + 1e-8 < counts.sum(axis=1)):
        raise ValueError(f'{name} must be at least the measured modeling-gene count sum')
    return library


def _probability(counts, library, specimens):
    labels = np.unique(specimens)
    per_sample = []
    for label in labels:
        use = specimens == label
        per_sample.append(counts[use].sum(axis=0) / library[use].sum())
    return np.mean(per_sample, axis=0)


def _pearson(counts, library, probability, theta, clip):
    mu = library[:, None] * probability[None, :]
    if np.isinf(theta):
        variance = mu
    else:
        variance = mu + mu**2 / theta
    residual = np.zeros_like(counts, dtype=float)
    informative = probability > 0
    residual[:, informative] = ((counts[:, informative] - mu[:, informative])
                                / np.sqrt(variance[:, informative]))
    return np.clip(residual, -clip, clip)


def fit_count_representation(counts, library, specimens, n_components=15,
                             theta=100., seed=15):
    """Fit specimen-balanced gene probabilities, residual center, and PCs.

    ``library`` is total raw-count exposure over all measured genes, while
    ``counts`` may contain only modeling genes. Equal specimen weight is used
    both for estimating gene probabilities and centering residuals. PCA uses
    randomized SVD without a second per-gene standardization.
    """
    counts = _counts_matrix(counts, 'counts')
    specimens = _labels(specimens, len(counts))
    library = _exposure(library, counts, 'library')
    if (not isinstance(n_components, (int, np.integer))
            or isinstance(n_components, (bool, np.bool_)) or n_components < 1):
        raise ValueError('n_components must be a positive integer')
    if not np.isscalar(theta):
        raise ValueError('theta must be positive finite or positive infinity')
    try:
        theta = float(theta)
    except (TypeError, ValueError) as exc:
        raise ValueError('theta must be positive finite or positive infinity') from exc
    if np.isnan(theta) or theta <= 0:
        raise ValueError('theta must be positive finite or positive infinity')
    k = min(int(n_components), len(counts) - 1, counts.shape[1])
    if k < 1:
        raise ValueError('at least two training rows and one gene are required')

    probability = _probability(counts, library, specimens)
    clip = float(np.sqrt(len(counts)))
    residual = _pearson(counts, library, probability, theta, clip)
    labels, inverse, group_counts = np.unique(specimens, return_inverse=True, return_counts=True)
    weights = 1. / (len(labels) * group_counts[inverse])
    center = weights @ residual
    centered = residual - center
    weighted = centered * np.sqrt(weights * len(counts))[:, None]
    _, _, components = randomized_svd(weighted, n_components=k, random_state=seed)
    transform = {'probability': probability, 'theta': theta, 'clip': clip,
                 'center': center, 'components': components,
                 'n_features': counts.shape[1], 'n_training': len(counts),
                 'training_specimens': labels}
    return centered @ components.T, transform


def transform_count_representation(counts, library, transform):
    """Transform query counts using frozen training probabilities and PCs."""
    counts = _counts_matrix(counts, 'counts', min_rows=1)
    if counts.shape[1] != transform['n_features']:
        raise ValueError('query count feature count does not match fitted representation')
    library = _exposure(library, counts, 'library')
    residual = _pearson(counts, library, transform['probability'],
                        transform['theta'], transform['clip'])
    if transform['center'].shape != (counts.shape[1],):
        raise ValueError('fitted residual center has invalid dimensions')
    return (residual - transform['center']) @ transform['components'].T
