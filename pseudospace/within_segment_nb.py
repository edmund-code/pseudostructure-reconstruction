"""Within-segment structure along the PT under the NB count model (notebook 74).

Notebooks 57-60 asked, with Gaussian models on log1p, what a continuous coordinate adds beyond S1/S2/S3 steps.
Notebook 74 asks it again with the NB model of notebooks 69/71 (raw counts, offset log of the mt-free library, a fixed
per-gene dispersion), so every statistic here is an NB likelihood ratio at a fixed alpha.

Every new gene test follows the positional F-test of notebook 71: an NB LR numerator judged against a specimen-level
LR denominator that **mirrors** it (if the numerator adds group x X beyond group x B, the denominator adds
within-group specimen x X beyond specimen x B), then ``pathway_pipelines.moderated_f``.

- ``beyond_step_designs``: a smooth group difference beyond group x segment steps (57's tests A and B). Shared steps
  enter every model (``shared_steps=True``, notebook 74's protocol), so that each group carries its own full step
  pattern in the reduced model; ``shared_steps=False`` is notebook 57's design, where only the case group had steps
  and the reference group's steps were left to the shared spline.
- ``step_ladder_designs``: group x steps beyond a group offset with the shared spline as nuisance (37's "steps on the
  spline" with labels, and its six coordinate bins).
- ``within_species_designs``: one species' shared gradient beyond its specimens' intercepts (W1), or beyond shared
  steps (W2, W3), against specimen-specific shape.
- ``nb_segment_slopes``: the signed root LR of a linear position term inside each specimen x segment (replaces the
  OLS slope t of notebook 58).
- ``nb_specimen_curves`` and ``pseudobulk_log_rates``: per-specimen NB curves and segment rates for 58's descriptors.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .nb_gam import nb_fit, nested_nb_lr
from .pathway_calibration import contrast_designs
from .stats_gam import make_gam_design

SEGMENTS = ('PT-S1', 'PT-S2', 'PT-S3')
NUMERATOR, DENOMINATOR = 'numerator', 'denominator'


def _indicators(labels):
    """0/1 columns of every level but the first (sorted), so they add steps beyond an intercept."""
    labels = np.asarray(labels).astype(str)
    levels = np.unique(labels)
    if len(levels) < 2:
        raise ValueError('Need at least two segments.')
    return np.column_stack([(labels == k).astype(float) for k in levels[1:]])


def _specimen_contrasts(specimen, group):
    """Within-group sum-to-zero specimen contrasts, coded as ``contrast_designs`` codes its specimen intercepts."""
    samples, g = np.asarray(specimen).astype(str), np.asarray(group, float)
    columns = []
    for level in (0., 1.):
        names = np.unique(samples[g == level])
        columns.extend((samples == name).astype(float) - (samples == names[-1]).astype(float) for name in names[:-1])
    if not columns:
        raise ValueError('Need at least two specimens in a group.')
    return np.column_stack(columns)


def _times(contrasts, block):
    return np.column_stack([c[:, None] * block for c in contrasts.T])


def beyond_step_designs(position, group, specimen, knots, segments, *, nuisance=None, shared_steps=True):
    """Designs of "group x spline beyond group x steps" and its mirrored specimen denominator.

    numerator: reduced (shared spline, specimen contrasts, group x steps [, shared steps, nuisance x steps]) ->
    full (+ group x spline); denominator: full + specimen-within-group x steps -> + specimen-within-group x spline.
    ``nuisance`` (0/1, the species when ``group`` is a relabeling) adds its spline shape to every model, as in
    ``contrast_designs``, and with ``shared_steps`` also its steps. Returns (designs, comparisons, weights).
    """
    s = np.asarray(position, float)
    designs, weights = contrast_designs(s, group, specimen, knots, segments=segments, nuisance_shape=nuisance)
    steps = _indicators(segments)
    extra = []
    if shared_steps:
        extra.append(steps)
        if nuisance is not None:
            extra.append(np.asarray(nuisance, float)[:, None] * steps)
    reduced = np.column_stack([designs['discrete'], *extra])
    full = np.column_stack([designs['union'], *extra])
    contrasts = _specimen_contrasts(specimen, group)
    den_reduced = np.column_stack([full, _times(contrasts, steps)])
    den_full = np.column_stack([den_reduced, _times(contrasts, make_gam_design(s, knots)[:, 1:])])
    out = {'reduced': reduced, 'full': full, 'den_reduced': den_reduced, 'den_full': den_full}
    return out, {NUMERATOR: ('reduced', 'full'), DENOMINATOR: ('den_reduced', 'den_full')}, weights


def step_ladder_designs(position, group, specimen, knots, segments, *, nuisance=None):
    """Designs of "group x steps beyond a group offset", with the shared spline and shared steps in every model.

    numerator: reduced (shared spline, shared steps, specimen contrasts, group [, nuisance spline and steps]) ->
    full (+ group x steps); denominator: full -> + specimen-within-group x steps. With the package labels this is 37's
    "steps on the spline"; with six equal-width position bins it is 37's six-bin row. Returns (designs, comparisons,
    weights).
    """
    s = np.asarray(position, float)
    designs, weights = contrast_designs(s, group, specimen, knots, segments=segments, nuisance_shape=nuisance)
    steps = _indicators(segments)
    extra = [steps] + ([np.asarray(nuisance, float)[:, None] * steps] if nuisance is not None else [])
    reduced = np.column_stack([designs['level'], *extra])
    full = np.column_stack([designs['discrete_shape'], *extra])
    den_full = np.column_stack([full, _times(_specimen_contrasts(specimen, group), steps)])
    return ({'reduced': reduced, 'full': full, 'den_full': den_full},
            {NUMERATOR: ('reduced', 'full'), DENOMINATOR: ('full', 'den_full')}, weights)


def within_species_designs(position, specimen, knots, segments=None):
    """One species' shared smooth gradient (beyond shared steps when ``segments`` is given) against specimen shape.

    Rows are one species' units. reduced: intercept, specimen sum-to-zero contrasts [, shared steps]; full: + spline;
    denominator: full [+ specimen x steps] -> + specimen x spline. Each specimen gets an equal share of the weight.
    Returns (designs, comparisons, weights).
    """
    s = np.asarray(position, float)
    samples = np.asarray(specimen).astype(str)
    names = np.unique(samples)
    if len(names) < 2:
        raise ValueError('Need at least two specimens.')
    n = len(s)
    weights = np.zeros(n)
    for name in names:
        mask = samples == name
        weights[mask] = n / (len(names) * mask.sum())
    contrasts = np.column_stack([(samples == name).astype(float) - (samples == names[-1]).astype(float)
                                 for name in names[:-1]])
    basis = make_gam_design(s, knots)[:, 1:]
    reduced = np.column_stack([np.ones(n), contrasts])
    if segments is not None:
        steps = _indicators(segments)
        reduced = np.column_stack([reduced, steps])
    full = np.column_stack([reduced, basis])
    shape = _times(contrasts, basis)
    if segments is None:
        designs = {'reduced': reduced, 'full': full, 'den_full': np.column_stack([full, shape])}
        return designs, {NUMERATOR: ('reduced', 'full'), DENOMINATOR: ('full', 'den_full')}, weights
    den_reduced = np.column_stack([full, _times(contrasts, steps)])
    designs = {'reduced': reduced, 'full': full, 'den_reduced': den_reduced,
               'den_full': np.column_stack([den_reduced, shape])}
    return designs, {NUMERATOR: ('reduced', 'full'), DENOMINATOR: ('den_reduced', 'den_full')}, weights


def nested_lr(counts, library, designs, comparisons, weights, dispersion, *, block_size=200, n_jobs=1, fill='median'):
    """``nb_gam.nested_nb_lr`` of a (designs, comparisons) pair: numerator and denominator LRs with their df.

    Returns dict numerator, denominator, df1, df2, converged, separated (nb_gam's rules: non-converged genes are filled
    by ``fill``, separated genes are flagged).
    """
    lib = np.asarray(library, float)
    if not np.isfinite(lib).all() or (lib <= 0).any():
        raise ValueError('Library sizes must be positive and finite.')
    lr = nested_nb_lr(counts, designs, weights, np.log(lib), dispersion, comparisons, block_size=block_size,
                      n_jobs=n_jobs, fill=fill)
    return {NUMERATOR: lr[NUMERATOR], DENOMINATOR: lr[DENOMINATOR], 'df1': int(lr['df'][NUMERATOR][0]),
            'df2': int(lr['df'][DENOMINATOR][0]), 'converged': lr['converged'], 'separated': lr['separated']}


def nb_segment_slopes(counts, library, position, specimen, segment, dispersion, *, order=SEGMENTS, min_rows=30,
                      block_size=200, n_jobs=1):
    """Signed root LR of a linear position term inside each specimen x segment (unit weights, fixed dispersion).

    sign(slope) x sqrt(D(intercept) - D(intercept + position)); NaN when the cell has fewer than ``min_rows`` units or
    the gene did not converge or is separated there. Returns {specimen: DataFrame genes x ``order``}.
    """
    from scipy import sparse

    s = np.asarray(position, float)
    spec, lab = np.asarray(specimen).astype(str), np.asarray(segment).astype(str)
    lib = np.asarray(library, float)
    counts = sparse.csr_matrix(counts) if sparse.issparse(counts) else np.asarray(counts)
    out = {}
    for name in np.unique(spec):
        cols = {}
        for seg in order:
            rows = np.flatnonzero((spec == name) & (lab == seg))
            if len(rows) < min_rows:
                cols[seg] = np.full(counts.shape[1], np.nan)
                continue
            y, off, w = counts[rows], np.log(lib[rows]), np.ones(len(rows))
            x = s[rows] - s[rows].mean()
            run = {'block_size': block_size, 'n_jobs': n_jobs}
            f0 = nb_fit(y, np.ones((len(rows), 1)), w, off, dispersion, **run)
            f1 = nb_fit(y, np.column_stack([np.ones(len(rows)), x]), w, off, dispersion, **run)
            lr = np.maximum(f0['deviance'] - f1['deviance'], 0.)
            ok = f0['converged'] & f1['converged'] & ~f1['separated'] & np.isfinite(lr)
            cols[seg] = np.where(ok, np.sign(f1['beta'][:, 1]) * np.sqrt(lr), np.nan)
        out[name] = pd.DataFrame(cols)
    return out


def nb_specimen_curves(counts, library, position, specimen, grid, knots, dispersion, *, quantiles=(.01, .99),
                       block_size=200, n_jobs=1):
    """Per-specimen NB curves (natural-log rate per library unit) on ``grid``, as ``standard_gene_model`` draws them.

    Unit weights, the shared knots and the fixed dispersion; NaN outside the specimen's [1st, 99th] percentile range and
    for genes that did not converge or are separated in that specimen. Returns {specimen: genes x grid}.
    """
    from scipy import sparse

    s = np.asarray(position, float)
    spec = np.asarray(specimen).astype(str)
    lib = np.asarray(library, float)
    grid = np.asarray(grid, float)
    base = make_gam_design(grid, knots)
    counts = sparse.csr_matrix(counts) if sparse.issparse(counts) else np.asarray(counts)
    out = {}
    for name in np.unique(spec):
        rows = np.flatnonzero(spec == name)
        fit = nb_fit(counts[rows], make_gam_design(s[rows], knots), np.ones(len(rows)), np.log(lib[rows]), dispersion,
                     block_size=block_size, n_jobs=n_jobs)
        curve = fit['beta'] @ base.T
        lo, hi = np.quantile(s[rows], quantiles)
        curve[:, (grid < lo) | (grid > hi)] = np.nan
        curve[~fit['converged'] | fit['separated']] = np.nan
        out[name] = curve
    return out


def pseudobulk_log_rates(counts, library, specimen, segment, *, order=SEGMENTS, floor=.5):
    """Specimen x segment pseudobulk natural-log rates ln(max(sum counts, floor) / sum library): {specimen: genes x 3}.

    NaN for a segment without units in that specimen.
    """
    from scipy import sparse

    spec, lab = np.asarray(specimen).astype(str), np.asarray(segment).astype(str)
    lib = np.asarray(library, float)
    counts = sparse.csr_matrix(counts)
    out = {}
    for name in np.unique(spec):
        cols = []
        for seg in order:
            rows = np.flatnonzero((spec == name) & (lab == seg))
            if not len(rows):
                cols.append(np.full(counts.shape[1], np.nan))
                continue
            total = np.asarray(counts[rows].sum(axis=0)).ravel()
            cols.append(np.log(np.maximum(total, floor) / lib[rows].sum()))
        out[name] = np.column_stack(cols)
    return out


def permuted_rows(blocks, rng):
    """Row index that permutes units within each block; apply it to counts and library together (positions and
    labels stay with their units). The same index for every gene keeps segment means and inter-gene correlation."""
    blocks = np.asarray(blocks)
    index = np.arange(len(blocks))
    for name in np.unique(blocks):
        rows = np.flatnonzero(blocks == name)
        index[rows] = rows[rng.permutation(len(rows))]
    return index
