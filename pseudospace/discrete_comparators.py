"""Discrete comparators of the positional F-test: the segment F-test, family-wide ORA and the swap rule (notebook 75).

The positional F-test (notebooks 67, 69, 71) asks whether the group difference changes along the PT beyond a constant
offset, with a 6-df spline in position. Its **segment twin** asks the same question with the S1/S2/S3 labels in place
of the spline, so the only change is continuous against discrete:

- rows are units, the model is an NB GLM on raw counts with offset log(library) and ``contrast_designs``' prior
  weights n / (2 J n_j) (each group half the weight, each of its specimens an equal share);
- ``segment_step_designs``: segment intercepts are the shared shape, with within-group sum-to-zero specimen
  contrasts and, for a relabeling, the species x segment shape as nuisance (``specimen_null.segment_designs``);
  ``level`` adds the group, ``steps`` adds group x segment, and ``spec`` adds within-group specimen contrasts x
  segment, so every specimen gets its own steps;
- the numerator T_step = LR(level -> steps), df1 = 2 for three segments; the denominator LR_spec = LR(steps -> spec),
  df2 = 4 for a condition split of 2 + 2 specimens and 2 for a balanced relabeling (from the design ranks);
- ``segment_f_test`` moderates LR_spec / df2 by limma's squeezeVar and forms F = (T_step / df1) / s2_post, its p and
  z = Phi^-1(1 - p) (``pathway_pipelines.moderated_f``), exactly as the positional F-test does with T_spatial.

The dispersion is fixed (the condition split's alpha of the spline model, the NB note's rule for the step screen), so
no dispersion is estimated here. ``segment_step_lr`` also returns the group difference in each segment: the case -
reference natural-log rate ratio at equal specimen weights, with its HC3 SE.

``family_ora`` is notebook 56's ORA: one-sided hypergeometric p of every set in every gene list, BH over the whole
family of set x list tests, a call when any test passes, and the list-based calibration p (the best p times the
number of lists, notebook 61's convention). ``swap_rule`` is the pass rule every row of notebook 75's method table
uses: calls > 0, NCDR < 0.25, each relabeling below 5% of the tested sets, and a pooled share of relabeled
calibration p <= 0.05 of at most 0.075.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import false_discovery_control, hypergeom

from .nb_gam import nb_fit, nested_nb_lr
from .specimen_null import segment_designs

SEGMENT_LRS = {'T_step': ('level', 'steps'), 'LR_spec': ('steps', 'spec')}
NCDR_MAX, NULL_RATE, CALIBRATION_MAX, ALPHA = .25, .05, .075, .05


def _per_specimen(values, samples, name):
    """{specimen: value} from per-unit values that must be constant within each specimen."""
    frame = pd.Series(np.asarray(values, float), index=samples)
    per = frame.groupby(level=0).agg(['min', 'max'])
    if not np.allclose(per['min'], per['max']):
        raise ValueError(f'{name} must be constant within each specimen.')
    return per['min'].to_dict()


def segment_step_designs(segment, group, specimen, *, nuisance=None):
    """The base, level, steps and spec designs of the segment F-test and their prior weights.

    segment, group (0/1), specimen and the optional nuisance (0/1, the species when ``group`` is a relabeling) are
    per unit; group and nuisance must be constant within specimen. Returns ({name: X}, weights, segment labels).
    ``steps`` is written without a separate group column (base + group x each segment), so its last columns are the
    group difference in each segment; ``level`` is base + group.
    """
    seg, samples = np.asarray(segment).astype(str), np.asarray(specimen).astype(str)
    g = np.asarray(group, float)
    n = len(seg)
    if g.shape != (n,) or samples.shape != (n,) or set(np.unique(g)) != {0., 1.}:
        raise ValueError('Need aligned segments, a 0/1 group and specimens.')
    group_of = _per_specimen(g, samples, 'group')
    nuisance_of = None if nuisance is None else _per_specimen(nuisance, samples, 'nuisance')
    labels = np.unique(seg)
    if len(labels) < 2:
        raise ValueError('Need at least two segments.')
    parts = segment_designs(samples, seg, group_of, nuisance_of)
    shape = np.column_stack([(seg == k).astype(float) for k in labels])
    base = parts['base']
    steps = np.column_stack([base, g[:, None] * shape])
    contrasts = []
    weights = np.zeros(n)
    for level in (0., 1.):
        names = sorted(name for name, v in group_of.items() if v == level)
        for name in names:
            mask = samples == name
            weights[mask] = n / (2 * len(names) * mask.sum())
        contrasts.extend((samples == name).astype(float) - (samples == names[-1]).astype(float) for name in names[:-1])
    spec = np.column_stack([steps, *[c[:, None] * shape[:, 1:] for c in contrasts]])
    return {'base': base, 'level': parts['level'], 'steps': steps, 'spec': spec}, weights, labels


def segment_step_lr(counts, library, segment, group, specimen, dispersion, *, nuisance=None, block_size=200, n_jobs=1,
                    fill='median'):
    """NB LRs of the segment F-test at a fixed dispersion, and the group difference per segment.

    counts: units x genes raw counts (dense or sparse); library: per-unit library (offset log library); dispersion:
    ``nb_dispersion``'s dict or one alpha per gene. Returns a dict: T_step and LR_spec (``nb_gam.nested_nb_lr``,
    non-converged genes filled by ``fill``), df1, df2, residual_df, converged (every used fit), separated, segments,
    delta and se_hc3 (genes x segments: case - reference natural-log rate ratio at equal specimen weights).
    """
    lib = np.asarray(library, float)
    if lib.shape != (np.asarray(segment).shape[0],) or not np.isfinite(lib).all() or (lib <= 0).any():
        raise ValueError('Library sizes must be positive, finite and aligned to units.')
    designs, weights, labels = segment_step_designs(segment, group, specimen, nuisance=nuisance)
    offset = np.log(lib)
    run = {'block_size': block_size, 'n_jobs': n_jobs}
    lr = nested_nb_lr(counts, {k: designs[k] for k in ('level', 'steps', 'spec')}, weights, offset, dispersion,
                      SEGMENT_LRS, fill=fill, **run)
    (df1, _), (df2, residual) = lr['df']['T_step'], lr['df']['LR_spec']
    k = len(labels)
    pick = np.zeros((k, designs['steps'].shape[1]))
    pick[:, -k:] = np.eye(k)
    fit = nb_fit(counts, designs['steps'], weights, offset, dispersion, contrasts={'delta': pick}, **run)
    return {'T_step': lr['T_step'], 'LR_spec': lr['LR_spec'], 'df1': int(df1), 'df2': int(df2),
            'residual_df': int(residual), 'converged': lr['converged'] & fit['converged'],
            'separated': lr['separated'] | fit['separated'], 'segments': [str(s) for s in labels],
            'delta': fit['contrasts']['delta']['value'], 'se_hc3': fit['contrasts']['delta']['se_hc3']}


def segment_f_test(counts, library, segment, group, specimen, dispersion, *, genes=None, nuisance=None, covariate=None,
                   robust=False, block_size=200, n_jobs=1):
    """The segment F-test of one partition: ``segment_step_lr`` then ``pathway_pipelines.moderated_f``.

    covariate: the squeezeVar trend covariate per gene (notebook 67: log mean raw count per unit). Returns one row per
    gene: moderated_f's s2, s2_prior, d0, s2_post, F, df_total, log_p, p, z; q (BH over all genes); T_step, LR_spec,
    df2, converged, separated; delta_<segment> and se_<segment>.
    """
    from .pathway_pipelines import moderated_f

    out = segment_step_lr(counts, library, segment, group, specimen, dispersion, nuisance=nuisance,
                          block_size=block_size, n_jobs=n_jobs)
    index = pd.RangeIndex(len(out['T_step'])) if genes is None else pd.Index(genes, name='gene')
    if len(index) != len(out['T_step']):
        raise ValueError('One gene name per count column.')
    frame = moderated_f(pd.Series(out['T_step'], index=index), out['df1'], out['LR_spec'], out['df2'],
                        covariate=covariate, robust=robust)
    frame['q'] = false_discovery_control(frame.p.to_numpy())
    frame['T_step'] = out['T_step']
    frame['LR_spec'] = out['LR_spec']
    frame['df2'] = out['df2']
    frame['converged'] = out['converged']
    frame['separated'] = out['separated']
    for j, s in enumerate(out['segments']):
        frame[f'delta_{s}'] = out['delta'][:, j]
        frame[f'se_{s}'] = out['se_hc3'][:, j]
    return frame


def membership(gene_sets, universe):
    """Sets x genes 0/1 CSR matrix over ``universe`` (members outside it are ignored), in ``gene_sets`` order."""
    index = pd.Index(list(universe))
    rows, cols = [], []
    for i, members in enumerate(gene_sets.values()):
        idx = index.get_indexer(list(members))
        idx = idx[idx >= 0]
        rows.extend([i] * len(idx))
        cols.extend(idx)
    return sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(gene_sets), len(index)))


def family_ora(selections, member_matrix, *, alpha=ALPHA):
    """Notebook 56's ORA over a family of gene lists.

    selections: list of boolean vectors, one per gene list, each over its background; member_matrix: sets x background
    genes (``membership``), either one matrix shared by every list or a list of matrices, one per list (lists drawn
    from different backgrounds, such as DE genes of different segments; every matrix has the same sets as rows). Each
    set gets a one-sided hypergeometric p in every list (1 for an empty list); BH runs over the whole family of set x
    list p; a set is called when any of its tests has q <= ``alpha``. Returns a DataFrame per set: p_<i> and q_<i> per
    list, best_p, best_list, call, and calibration_p = min(1, number of lists x best_p).
    """
    selections = list(selections)
    if not selections:
        raise ValueError('Need at least one gene list.')
    matrices = member_matrix if isinstance(member_matrix, (list, tuple)) else [member_matrix] * len(selections)
    if len(matrices) != len(selections):
        raise ValueError('One member matrix per gene list, or one for all.')
    matrices = [sparse.csr_matrix(m, dtype=float) for m in matrices]
    if len({m.shape[0] for m in matrices}) != 1:
        raise ValueError('Every member matrix must hold the same sets.')
    columns = []
    for selected, member in zip(selections, matrices):
        selected = np.asarray(selected, bool)
        n_all = member.shape[1]
        if selected.shape != (n_all,):
            raise ValueError('Every selection must be a boolean vector over its background genes.')
        n_sel = int(selected.sum())
        if n_sel == 0:
            columns.append(np.ones(member.shape[0]))
            continue
        k = np.asarray(member @ selected.astype(float)).ravel()
        size = np.asarray(member.sum(axis=1)).ravel()
        columns.append(hypergeom.sf(k - 1, n_all, size, n_sel))
    p = np.column_stack(columns)
    q = false_discovery_control(p.ravel()).reshape(p.shape)
    out = pd.DataFrame({**{f'p_{i}': p[:, i] for i in range(p.shape[1])}, **{f'q_{i}': q[:, i] for i in range(p.shape[1])}})
    out['best_p'] = p.min(axis=1)
    out['best_list'] = p.argmin(axis=1)
    out['call'] = (q <= alpha).any(axis=1)
    out['calibration_p'] = np.minimum(1., p.shape[1] * out.best_p.to_numpy())
    return out


def swap_rule(condition_calls, relabeled_calls, tested, calibration_p, *, ncdr_max=NCDR_MAX, null_rate=NULL_RATE,
              calibration_max=CALIBRATION_MAX, alpha=ALPHA):
    """Notebook 75's pass rule for one row on one dataset.

    condition_calls: the condition split's call count (or several, for a row with several condition runs: their mean
    is used); relabeled_calls: the two relabelings' counts; tested: tested sets; calibration_p: the relabelings'
    per-set calibration p, pooled. Pass = calls > 0, NCDR < ``ncdr_max``, every relabeling < ``null_rate`` x tested
    and a calibration share (pooled share of p <= ``alpha``) <= ``calibration_max``.
    """
    calls = float(np.mean(np.atleast_1d(np.asarray(condition_calls, float))))
    relabeled = np.atleast_1d(np.asarray(relabeled_calls, float))
    p = np.asarray(calibration_p, float).ravel()
    if tested <= 0 or not len(relabeled) or not len(p) or np.isnan(p).any():
        raise ValueError('Need tested sets, relabeled calls and finite calibration p.')
    ncdr = relabeled.mean() / calls if calls > 0 else np.nan
    relabeling = bool(calls > 0 and ncdr < ncdr_max and (relabeled < null_rate * tested).all())
    share = float(np.mean(p <= alpha))
    calibration = share <= calibration_max
    return {'calls': calls, 'NCDR': ncdr, 'relabeling rule': 'pass' if relabeling else 'fail',
            'calibration share': share, 'calibration': 'pass' if calibration else 'fail',
            'pass': 'pass' if relabeling and calibration else 'fail'}
