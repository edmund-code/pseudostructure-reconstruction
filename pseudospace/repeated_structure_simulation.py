"""Synthetic unordered repeated-structure count data for method checks."""
from __future__ import annotations

import numpy as np


def simulate_repeated_sections(*, seed, n_per_specimen=500,
                               specimens=('A', 'B'), interval=(0., 1.),
                               spatial_strength=3., nuisance_strength=3.,
                               nuisance_mode='independent', realization_sd=.1,
                               exposure_multiplier=1.):
    """Generate balanced specimens with an unordered spatial/nuisance count program.

    Each row is an independent structure realization; no nephron identities or
    trajectories are generated. ``nuisance_mode`` can make nuisance q
    independent of, aligned with, or reversed relative to the full-axis signal.
    """
    if (not isinstance(n_per_specimen, (int, np.integer))
            or isinstance(n_per_specimen, (bool, np.bool_)) or n_per_specimen < 20):
        raise ValueError('n_per_specimen must be an integer of at least 20')
    raw_specimens = np.asarray(specimens, dtype=object)
    if raw_specimens.ndim != 1 or len(raw_specimens) < 1 or any(
            not isinstance(x, (str, np.str_)) or not str(x).strip() for x in raw_specimens):
        raise ValueError('specimens must contain one or more nonempty strings')
    specimen_labels = np.asarray([str(x) for x in raw_specimens], dtype=str)
    if len(np.unique(specimen_labels)) != len(specimen_labels):
        raise ValueError('specimens must be unique')
    interval = np.asarray(interval, dtype=float)
    if interval.shape != (2,) or not np.isfinite(interval).all() or not 0 <= interval[0] < interval[1] <= 1:
        raise ValueError('interval must satisfy 0 <= lower < upper <= 1')
    if not np.isscalar(spatial_strength) or not np.isfinite(spatial_strength) or spatial_strength <= 0:
        raise ValueError('spatial_strength must be positive and finite')
    if not np.isscalar(nuisance_strength) or not np.isfinite(nuisance_strength) or nuisance_strength < 0:
        raise ValueError('nuisance_strength must be nonnegative and finite')
    if not np.isscalar(realization_sd) or not np.isfinite(realization_sd) or realization_sd < 0:
        raise ValueError('realization_sd must be nonnegative and finite')
    if not np.isscalar(exposure_multiplier) or not np.isfinite(exposure_multiplier) or exposure_multiplier <= 0:
        raise ValueError('exposure_multiplier must be positive and finite')
    if nuisance_mode not in {'independent', 'coupled', 'reversed'}:
        raise ValueError("nuisance_mode must be 'independent', 'coupled', or 'reversed'")
    if not isinstance(seed, (int, np.integer)) or isinstance(seed, (bool, np.bool_)):
        raise ValueError('seed must be an integer')

    rng = np.random.default_rng(seed)
    n_specimens = len(specimen_labels)
    n = n_per_specimen * n_specimens
    z = rng.uniform(interval[0], interval[1], size=n)
    specimens = np.repeat(specimen_labels, n_per_specimen)
    if nuisance_mode == 'independent':
        q = rng.normal(size=n)
    else:
        q = np.sqrt(3.) * (2*z - 1)
        if nuisance_mode == 'reversed':
            q = -q

    spatial_direction = np.tile(np.array([1., -1., .6, -.6]), 10)
    spatial_direction /= np.linalg.norm(spatial_direction)
    nuisance_direction = np.tile(np.array([1., 1., -1., -1.]), 10)
    nuisance_direction /= np.linalg.norm(nuisance_direction)
    logits = (spatial_strength * (2*z - 1)[:, None] * spatial_direction[None, :]
              + nuisance_strength * q[:, None] * nuisance_direction[None, :]
              + rng.normal(0., realization_sd, size=(n, 40)))
    logits -= logits.max(axis=1, keepdims=True)
    probabilities = np.exp(logits)
    probabilities /= probabilities.sum(axis=1, keepdims=True)
    library = np.maximum(1, np.rint(
        rng.lognormal(np.log(8000.), .35, size=n) * exposure_multiplier
    )).astype(np.int64)
    counts = np.vstack([rng.multinomial(int(total), probability)
                        for total, probability in zip(library, probabilities)])
    anatomy = np.where(z < 1/3, 0, np.where(z > 2/3, 2, 1)).astype(int)
    gene_names = np.asarray([f'Gene{i:02d}' for i in range(40)], dtype=str)
    structure_ids = np.asarray([
        f'{specimen}_{i:05d}'
        for specimen in specimen_labels
        for i in range(n_per_specimen)
    ], dtype=str)
    order = rng.permutation(n)
    return {'counts': counts[order], 'library': library[order],
            'specimens': specimens[order], 'anatomy': anatomy[order],
            'gene_names': gene_names, 'structure_ids': structure_ids[order],
            'true_z': z[order], 'nuisance_q': q[order]}
