"""Small summaries of pathway trajectories across spatial resolutions."""
from __future__ import annotations

from collections.abc import Mapping

import numpy as np
import pandas as pd

from .levelshape import build_ls_designs


def aggregate_saved_trajectories(fit, position, human, specimen, gene_sets):
    """Average saved gene trajectories and compute pathway-mean HC3 bands.

    Gene residuals are averaged within each pathway before HC3 variance is
    formed, preserving cross-gene residual covariance. This summarizes the
    already-fitted curves; it does not refit expression models. ``gene_sets``
    must contain the caller's common count-supported pathway memberships.
    """
    if not isinstance(gene_sets, Mapping):
        raise ValueError('gene_sets must map pathway IDs to gene members.')
    genes = np.asarray(fit['genes']).astype(str)
    grid = np.asarray(fit['grid'], dtype=float)
    knots = np.asarray(fit['knots'], dtype=float)
    position = np.asarray(position, dtype=float)
    human = np.asarray(human, dtype=float)
    specimen = np.asarray(specimen)
    if pd.isna(specimen).any():
        raise ValueError('Specimen labels must be present for every structure.')
    specimen = specimen.astype(str)
    mouse_fit = np.asarray(fit['mouse'], dtype=float)
    human_fit = np.asarray(fit['human'], dtype=float)
    residuals = np.asarray(fit['residuals'], dtype=float)
    if residuals.ndim != 2:
        raise ValueError('Saved full-model residuals must be a structures-by-genes matrix.')
    n, g = residuals.shape
    if (genes.ndim != 1 or len(genes) != g or len(np.unique(genes)) != g
            or position.shape != (n,) or human.shape != (n,) or specimen.shape != (n,)
            or mouse_fit.shape != (g, len(grid)) or human_fit.shape != mouse_fit.shape
            or not np.isfinite(position).all() or not np.isfinite(human).all()
            or not np.isfinite(grid).all() or not np.isfinite(knots).all()
            or not np.isfinite(mouse_fit).all() or not np.isfinite(human_fit).all()
            or not np.isfinite(residuals).all()):
        raise ValueError('Saved fit arrays, genes, and structure metadata must align and be finite.')
    if set(np.unique(human)) != {0., 1.} or len(np.unique(specimen)) < 2:
        raise ValueError('Need both species (mouse=0, human=1) and multiple specimens.')
    if position.min() >= position.max() or grid.size < 2 or np.any(np.diff(grid) <= 0):
        raise ValueError('Need a nonzero position range and increasing fitted grid.')
    if grid.min() < position.min() - 1e-12 or grid.max() > position.max() + 1e-12:
        raise ValueError('Fitted grid must lie within observed common support.')

    weights = np.zeros(n)
    nuisance = []
    for group in (0., 1.):
        names = np.unique(specimen[human == group])
        for name in names:
            mask = specimen == name
            if not np.all(human[mask] == group):
                raise ValueError('A specimen cannot belong to both species.')
            weights[mask] = n / (2 * len(names) * mask.sum())
        nuisance.extend((specimen == name).astype(float) - (specimen == names[-1]).astype(float)
                        for name in names[:-1])

    x2 = build_ls_designs(position, human, knots)[2]
    design = np.column_stack([x2, *nuisance])
    if np.linalg.matrix_rank(design) != design.shape[1] or n <= design.shape[1] + 2:
        raise ValueError('Insufficient support for the saved full-model design.')
    inverse = np.linalg.inv(design.T @ (weights[:, None] * design))
    solver = inverse @ (design.T * weights)
    leverage = np.einsum('ij,ji->i', design, solver)
    if np.any(leverage >= .99):
        raise ValueError('Near-unit leverage in saved full-model design.')

    grid_mouse = build_ls_designs(grid, np.zeros(len(grid)), knots)[2]
    grid_human = build_ls_designs(grid, np.ones(len(grid)), knots)[2]
    grid_mouse = np.column_stack([grid_mouse, np.zeros((len(grid), len(nuisance)))])
    grid_human = np.column_stack([grid_human, np.zeros((len(grid), len(nuisance)))])
    # Guard the expected Mfull ordering [base, species, species:basis, specimen].
    if grid_mouse.shape[1] != design.shape[1]:
        raise ValueError('Saved full-model parameterization does not match the current spline design.')
    delta_projection = (grid_human - grid_mouse) @ solver
    mouse_projection = grid_mouse @ solver
    human_projection = grid_human @ solver
    correction = 1 / (1 - leverage)

    gene_index = pd.Index(genes)
    rows = []
    for pathway_id, members in gene_sets.items():
        members = list(members)
        if not members or len(set(members)) != len(members):
            raise ValueError(f'{pathway_id}: pathway members must be nonempty and unique.')
        indices = gene_index.get_indexer(members)
        if np.any(indices < 0):
            raise ValueError(f'{pathway_id}: unknown genes in pathway membership.')
        residual_mean = residuals[:, indices].mean(axis=1)
        hc3_score_sq = (residual_mean * correction) ** 2
        se_mouse = np.sqrt(np.maximum(mouse_projection ** 2 @ hc3_score_sq, 1e-12))
        se_human = np.sqrt(np.maximum(human_projection ** 2 @ hc3_score_sq, 1e-12))
        se_delta = np.sqrt(np.maximum(delta_projection ** 2 @ hc3_score_sq, 1e-12))
        rows.append(pd.DataFrame({
            'pathway_id': pathway_id,
            'original_position': grid,
            'position': (grid - grid.min()) / (grid.max() - grid.min()),
            'mouse': mouse_fit[indices].mean(axis=0),
            'human': human_fit[indices].mean(axis=0),
            'delta': human_fit[indices].mean(axis=0) - mouse_fit[indices].mean(axis=0),
            'se_mouse': se_mouse,
            'se_human': se_human,
            'se_delta': se_delta,
            'n_genes': len(indices),
        }))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=[
        'pathway_id', 'original_position', 'position', 'mouse', 'human', 'delta',
        'se_mouse', 'se_human', 'se_delta', 'n_genes'])
