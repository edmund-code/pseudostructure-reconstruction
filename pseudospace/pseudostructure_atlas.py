"""Prototype wrapper for a frozen repeated-specimen pseudostructure atlas.

This composes the existing count residual, joint Gaussian mean, and
conditional count-rate estimators. Their coordinates are exploratory: neither
Gaussian/NB posterior precision is calibrated nor anatomy established as truth.
"""
from __future__ import annotations

import numpy as np
from scipy import sparse

from pseudospace.count_rate_atlas import fit_count_rate_atlas, project_count_rate_atlas
from pseudospace.count_representation import (
    _counts_matrix, _exposure, _labels, fit_count_representation,
    transform_count_representation,
)
from pseudospace.joint_repeated_atlas import (
    fit_joint_atlas, project_joint_atlas,
)
from pseudospace.repeated_atlas_baselines import _anatomy


def _gene_names(values, expected_length=None):
    names = np.asarray(values, dtype=object)
    if names.ndim != 1 or (expected_length is not None and len(names) != expected_length):
        raise ValueError('gene_names must be a one-dimensional vector matching count columns')
    if any(not isinstance(x, (str, np.str_)) or not str(x).strip() for x in names):
        raise ValueError('gene_names must contain nonempty strings')
    names = np.asarray([str(x) for x in names], dtype=str)
    folded = [x.casefold() for x in names]
    if len(set(folded)) != len(folded):
        raise ValueError('gene_names must be unique after casefolding')
    return names


def _structure_ids(values, n):
    ids = np.asarray(values, dtype=object)
    if ids.ndim != 1 or len(ids) != n:
        raise ValueError('structure_ids must be a one-dimensional vector matching query rows')
    if any(not isinstance(x, (str, np.str_)) or not str(x).strip() for x in ids):
        raise ValueError('structure_ids must contain nonempty strings')
    ids = np.asarray([str(x) for x in ids], dtype=str)
    if len(set(ids.tolist())) != n:
        raise ValueError('structure_ids must be unique')
    return ids


def _gaussian_projection(atlas, Y, specimens=None):
    output = {key: np.empty(len(Y), dtype=float)
              for key in ('z_mean', 'z_MAP', 'entropy', 'residual_norm')}
    output['posterior'] = np.empty((len(Y), len(atlas['grid'])), dtype=float)
    if specimens is None:
        groups = [(np.arange(len(Y)), None)]
    else:
        groups = [(np.flatnonzero(specimens == specimen), specimen)
                  for specimen in np.unique(specimens)]
    for rows, specimen in groups:
        offset = atlas['offsets'].get(specimen, np.zeros(Y.shape[1]))
        projected = project_joint_atlas(atlas, Y[rows], offset=offset)
        for key in output:
            output[key][rows] = projected[key]
    return output


def fit_pseudostructure_atlas(counts, library, specimens, anatomy,
                              gene_names, model_gene_mask, *,
                              n_components=15, seed=15, grid=None,
                              roughness=.001, offset_penalty=4.,
                              anatomy_strength=1., anatomy_width=.1,
                              max_iter=60, tol=1e-5,
                              representation_theta=100., count_theta=100.,
                              rate_bandwidth=.08, min_effective=12.):
    """Fit training-only count-Gaussian and count-rate references.

    ``counts`` contains the supplied measured features. ``library`` is raw
    all-gene exposure, including features excluded from modeling. Sparse inputs
    are densified; this prototype targets filtered structure aggregates. ``model_gene_mask`` is an explicit caller-controlled feature
    contract; marker exclusions and availability decisions remain upstream.
    """
    counts = _counts_matrix(counts.toarray() if sparse.issparse(counts) else counts, 'counts')
    n, g = counts.shape
    library = _exposure(library, counts, 'library')
    specimens = _labels(specimens, n)
    anatomy = _anatomy(anatomy, n)
    names = _gene_names(gene_names, g)
    raw_mask = np.asarray(model_gene_mask)
    if raw_mask.ndim != 1 or len(raw_mask) != g or raw_mask.dtype.kind != 'b':
        raise ValueError('model_gene_mask must be a boolean vector matching genes')
    mask = raw_mask.astype(bool, copy=True)
    if mask.sum() < 2:
        raise ValueError('model_gene_mask must include at least two genes')
    if np.sum(anatomy == 0) < 2 or np.sum(anatomy == 2) < 2:
        raise ValueError('training requires at least two S1 and two S3 labels for orientation')
    if (not isinstance(n_components, (int, np.integer))
            or isinstance(n_components, (bool, np.bool_)) or n_components < 1):
        raise ValueError('n_components must be a positive integer')
    if (not isinstance(seed, (int, np.integer)) or isinstance(seed, (bool, np.bool_))):
        raise ValueError('seed must be an integer')

    fit_grid = None if grid is None else np.array(grid, dtype=float, copy=True)
    modeling_counts = counts[:, mask]
    Y, representation = fit_count_representation(
        modeling_counts, library, specimens, n_components=n_components,
        theta=representation_theta, seed=seed)
    initial = Y[:, 0].copy()
    if np.median(initial[anatomy == 2]) < np.median(initial[anatomy == 0]):
        initial *= -1
    low, high = np.quantile(initial, [.01, .99])
    if not np.isfinite([low, high]).all() or high <= low:
        raise ValueError('training PC1 initialization has no robust range')
    initial = np.clip((initial - low) / (high - low), 0., 1.)
    gaussian = fit_joint_atlas(
        Y, specimens, anatomy, initial_z=initial, grid=fit_grid,
        roughness=roughness, offset_penalty=offset_penalty,
        anatomy_strength=anatomy_strength, anatomy_width=anatomy_width,
        max_iter=max_iter, tol=tol)
    if not gaussian['converged']:
        raise RuntimeError('joint Gaussian atlas did not converge')
    training_projection = _gaussian_projection(gaussian, Y, specimens)
    training_z = training_projection['z_mean']
    rate = fit_count_rate_atlas(
        modeling_counts, library, training_z, specimens, grid=fit_grid,
        bandwidth=rate_bandwidth, theta=count_theta,
        min_effective=min_effective)
    params = {'n_components': int(n_components), 'seed': int(seed),
              'grid': gaussian['grid'].copy(), 'roughness': float(roughness),
              'offset_penalty': float(offset_penalty),
              'anatomy_strength': float(anatomy_strength),
              'anatomy_width': float(anatomy_width), 'max_iter': int(max_iter),
              'tol': float(tol), 'representation_theta': float(representation_theta),
              'count_theta': float(count_theta), 'rate_bandwidth': float(rate_bandwidth),
              'min_effective': float(min_effective),
              'initialization_quantiles': (float(low), float(high))}
    return {'gene_names': names, 'model_gene_mask': mask,
            'representation': representation, 'gaussian': gaussian,
            'rate': rate, 'training_z': training_z,
            'training_specimens': np.unique(specimens), 'params': params}


def project_pseudostructure_atlas(reference, counts, library, gene_names,
                                  structure_ids, *, specimens=None):
    """Project counts/exposure through frozen models; no query labels or refit."""
    counts = _counts_matrix(counts.toarray() if sparse.issparse(counts) else counts, 'counts', min_rows=1)
    names = _gene_names(gene_names, counts.shape[1])
    if not np.array_equal(names, reference['gene_names']):
        raise ValueError('query gene names and order must exactly match reference')
    library = _exposure(library, counts, 'library')
    ids = _structure_ids(structure_ids, len(counts))
    if specimens is None:
        query_specimens = None
    else:
        query_specimens = _labels(specimens, len(counts))
    mask = reference['model_gene_mask']
    model_counts = counts[:, mask]
    Y = transform_count_representation(model_counts, library, reference['representation'])
    gaussian = _gaussian_projection(reference['gaussian'], Y, query_specimens)
    count = project_count_rate_atlas(reference['rate'], model_counts, library)
    disagreement = np.abs(gaussian['z_mean'] - count['z_mean'])
    return {'structure_ids': ids, 'gaussian': gaussian, 'count': count,
            'coordinate_disagreement': disagreement}


def summarize_panel_sensitivity(projections):
    """Report sensitivity in a common numerical gauge, not calibrated intervals.

    Panel references must use the same axis convention and profile order.
    Conditional posteriors are retained separately rather than averaged.
    """
    if not isinstance(projections, (list, tuple)) or len(projections) < 2:
        raise ValueError('at least two aligned panel projections are required')
    ids = np.asarray(projections[0]['structure_ids'], dtype=str)
    n = len(ids)
    ids = _structure_ids(ids, n)
    if n == 0:
        raise ValueError('panel projections must include query profiles')
    for projection in projections:
        current = np.asarray(projection['structure_ids'], dtype=str)
        if current.ndim != 1 or not np.array_equal(current, ids):
            raise ValueError('panel projections must have identical structure ID order')
    summary = {'structure_ids': ids, 'n_panels': len(projections)}
    for model_name in ('gaussian', 'count'):
        vectors = [np.asarray(p[model_name]['z_mean'], float) for p in projections]
        entropy_vectors = [np.asarray(p[model_name]['entropy'], float) for p in projections]
        if any(v.ndim != 1 or len(v) != n for v in vectors + entropy_vectors):
            raise ValueError(f'{model_name} projection vectors must match structure IDs')
        coordinates = np.vstack(vectors)
        entropies = np.vstack(entropy_vectors)
        if coordinates.shape != (len(projections), n) or entropies.shape != coordinates.shape:
            raise ValueError(f'{model_name} projection vectors must match structure IDs')
        if (not np.isfinite(coordinates).all() or np.any((coordinates < 0) | (coordinates > 1))
                or not np.isfinite(entropies).all() or np.any(entropies < -1e-12)):
            raise ValueError(f'{model_name} positions/entropy must be finite and valid')
        summary[model_name] = {
            'z_median': np.median(coordinates, axis=0),
            'z_min': np.min(coordinates, axis=0),
            'z_max': np.max(coordinates, axis=0),
            'z_range': np.ptp(coordinates, axis=0),
            'entropy_mean': np.mean(entropies, axis=0),
            'entropy_min': np.min(entropies, axis=0),
        }
    return summary
