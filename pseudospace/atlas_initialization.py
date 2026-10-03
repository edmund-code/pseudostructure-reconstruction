"""Training-only coordinate initializers for repeated-specimen atlases."""
from __future__ import annotations

import numpy as np

from pseudospace.count_representation import _labels
from pseudospace.repeated_atlas_baselines import (
    _anatomy, _matrix, fit_pc1_reference, fit_segment_mixture,
    predict_segment_mixture, project_pc1_reference,
)


def initialize_atlas_positions(Y, anatomy, specimens, *, method='pc1', seed=15):
    """Initialize atlas z from PC1, an endpoint contrast, or soft segments.

    These are training coordinate initializers only. ``soft_segment`` returns
    an expected ordinal label from expression, rather than known fine
    positions; its probabilities can vary continuously within segments.
    """
    Y = _matrix(Y, 'Y')
    anatomy = _anatomy(anatomy, len(Y))
    specimens = _labels(specimens, len(Y))
    if not isinstance(seed, (int, np.integer)) or isinstance(seed, (bool, np.bool_)):
        raise ValueError('seed must be an integer')
    if not isinstance(method, str) or method not in {'pc1', 'endpoint_contrast', 'soft_segment'}:
        raise ValueError("method must be 'pc1', 'endpoint_contrast', or 'soft_segment'")
    if not np.any(anatomy == 0) or not np.any(anatomy == 2):
        raise ValueError('anatomy must contain endpoint labels 0 and 2')

    if method == 'pc1':
        reference = fit_pc1_reference(Y, anatomy)
        z = project_pc1_reference(reference, Y)
        return z, {'method': method, 'low': reference['low'],
                   'high': reference['high'],
                   'contract': 'existing first PC, oriented S1-to-S3 and scaled by training 1st/99th percentiles'}

    if method == 'endpoint_contrast':
        per_label = {}
        for label in (0, 2):
            specimen_means = []
            for specimen in np.unique(specimens):
                use = (specimens == specimen) & (anatomy == label)
                if np.any(use):
                    specimen_means.append(Y[use].mean(axis=0))
            per_label[label] = np.mean(specimen_means, axis=0)
        direction = per_label[2] - per_label[0]
        norm = float(np.linalg.norm(direction))
        if not np.isfinite(norm) or norm <= np.finfo(float).eps:
            raise ValueError('endpoint contrast is degenerate')
        direction /= norm
        score = Y @ direction
        low, high = np.quantile(score, [.01, .99])
        if not np.isfinite([low, high]).all() or high <= low:
            raise ValueError('endpoint contrast score has no robust training range')
        z = np.clip((score - low) / (high - low), 0., 1.)
        return z, {'method': method, 'direction': direction,
                   'low': float(low), 'high': float(high),
                   'contract': 'equal-specimen S3-minus-S1 mean direction, scaled by training 1st/99th percentiles'}

    mixture = fit_segment_mixture(Y, anatomy[:, None]/2., anatomy, specimens, C=1., seed=seed)
    output = predict_segment_mixture(mixture, Y)
    z = output['probabilities'] @ np.arange(3, dtype=float) / 2.
    return z, {'method': method,
               'contract': 'expected ordinal segment code from training-only soft classifier; no true fine positions'}
