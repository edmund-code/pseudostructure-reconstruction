"""Exact 2v2 specimen-permutation p-values: inclusive tail, ties counted, mirror-pair floor.

The permutation set is the six ways to split four specimens into two pairs, so the observed
labeling is one of them and the smallest attainable p is 2/6 (a feature and its condition-swapped
mirror tie). These tests pin that behaviour down without private data.
"""
from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("patsy")
pytest.importorskip("statsmodels")

from pseudospace.levelshape import run_level_shape, sample_perm_pvalues
from pseudospace.stats_gam import gam_internal_knots

SAMPLES = ("C1", "C2", "A1", "A2")
CONTROL_SAMPLES = ("C1", "C2")
LAMBDA_GRID = np.logspace(-3, 3, 5)
GRID = np.linspace(0.0, 1.0, 21)
N_PER_SAMPLE = 25


def _panel(curve_by_sample, *, noise=0.0, seed=0):
    """Long-form (Y, s, c, samples) for four specimens and two genes."""
    rng = np.random.default_rng(seed)
    s, c, samples, rows = [], [], [], []
    for name in SAMPLES:
        x = np.linspace(0.0, 1.0, N_PER_SAMPLE)
        gene_1 = curve_by_sample[name](x)
        gene_2 = np.full_like(x, 2.0)                 # flat: no pseudospace structure
        values = np.column_stack([gene_1, gene_2])
        if noise:
            values = values + rng.normal(0.0, noise, size=values.shape)
        rows.append(values)
        s.append(x)
        c.append(np.full(N_PER_SAMPLE, 1.0 if name in ("A1", "A2") else 0.0))
        samples.extend([name] * N_PER_SAMPLE)
    return (np.vstack(rows), np.concatenate(s), np.concatenate(c), np.array(samples))


def _perm_pvalues(Y, s, c, samples, debug=False):
    knots = gam_internal_knots(s, basis_df=6)
    fit = run_level_shape(Y, s, c, knots, GRID, LAMBDA_GRID)
    out = sample_perm_pvalues(
        Y, s, samples, knots, GRID, LAMBDA_GRID, fit["lam_idx"], fit["sl2"], fit["p_b"],
        CONTROL_SAMPLES,
    )
    if debug:
        print("shape_p", out[0], "level_p", out[1], "true_idx", out[3])
    return out


def test_identical_specimen_curves_give_p_one_for_every_gene():
    """All six relabelings reproduce the observed statistic -> exact p is 1, not 1/6."""
    flat = lambda x: 1.0 + 2.0 * x                      # noqa: E731 - every specimen identical
    Y, s, c, samples = _panel({name: flat for name in SAMPLES})
    shape_p, level_p, splits, _ = _perm_pvalues(Y, s, c, samples)
    assert len(splits) == 6
    assert np.allclose(shape_p, 1.0), shape_p
    assert np.allclose(level_p, 1.0), level_p


def test_crossing_conditions_hit_the_two_of_six_mirror_floor():
    """A maximal shape effect ties with its condition-swapped mirror -> p = 2/6, never 1/6."""
    curves = {
        "C1": lambda x: 1.0 + 3.0 * x,
        "C2": lambda x: 1.0 + 3.0 * x,
        "A1": lambda x: 4.0 - 3.0 * x,                  # crossing: level gap ~ 0, shape large
        "A2": lambda x: 4.0 - 3.0 * x,
    }
    Y, s, c, samples = _panel(curves)
    shape_p, level_p, splits, true_idx = _perm_pvalues(Y, s, c, samples)
    assert set(splits[true_idx]) == set(CONTROL_SAMPLES)
    assert shape_p[0] == pytest.approx(2 / 6), shape_p
    # The observed labeling has no mean gap, and only its condition-swapped mirror is equally
    # gap-free, so an exact two-sided level test also lands on the 2/6 floor here.
    assert level_p[0] == pytest.approx(2 / 6), level_p
    assert np.all(shape_p >= 1 / 6) and np.all(shape_p <= 1.0)
    assert np.all(level_p >= 1 / 6) and np.all(level_p <= 1.0)


def test_pure_level_shift_is_seen_by_the_level_statistic_not_the_shape_statistic():
    """A constant offset between conditions must not be reported as a shape change."""
    curves = {
        "C1": lambda x: 1.0 + 3.0 * x, "C2": lambda x: 1.0 + 3.0 * x,
        "A1": lambda x: 3.0 + 3.0 * x, "A2": lambda x: 3.0 + 3.0 * x,
    }
    Y, s, c, samples = _panel(curves)
    shape_p, level_p, _, _ = _perm_pvalues(Y, s, c, samples)
    assert level_p[0] == pytest.approx(2 / 6), level_p
    assert shape_p[0] == pytest.approx(1.0), shape_p


def test_permutation_pvalue_is_a_valid_tail_and_monotone_in_effect_size():
    """Weaker effects cannot be 'more significant' than stronger ones sharing the same design."""
    strong = {
        "C1": lambda x: 1.0 + 3.0 * x, "C2": lambda x: 1.0 + 3.0 * x,
        "A1": lambda x: 4.0 - 3.0 * x, "A2": lambda x: 4.0 - 3.0 * x,
    }
    weak = {
        "C1": lambda x: 1.0 + 3.0 * x, "C2": lambda x: 1.0 + 3.0 * x,
        "A1": lambda x: 2.0 + 1.0 * x, "A2": lambda x: 2.0 + 1.0 * x,
    }
    strong_p = _perm_pvalues(*_panel(strong))[0][0]
    weak_p = _perm_pvalues(*_panel(weak))[0][0]
    assert strong_p <= weak_p, (strong_p, weak_p)
    assert np.isfinite(strong_p) and np.isfinite(weak_p)


def test_ties_are_counted_as_at_least_as_extreme():
    """Noisy-but-identical condition curves tie across relabelings; p stays inclusive."""
    flat = lambda x: 1.0 + 2.0 * x                      # noqa: E731
    Y, s, c, samples = _panel({name: flat for name in SAMPLES}, noise=0.05, seed=3)
    shape_p, _, _, _ = _perm_pvalues(Y, s, c, samples)
    assert np.all(shape_p >= 1 / 6), shape_p
    # Near-identical condition curves: most relabelings must be counted, so p cannot be tiny.
    assert shape_p[0] >= 3 / 6, shape_p
