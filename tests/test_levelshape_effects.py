"""Level / amplitude / spatial-pattern decomposition of fitted curve pairs.

`shape_rms` alone reports an amplitude loss as a shape change; these tests pin the separation of
the three questions without private data.
"""
from __future__ import annotations

import numpy as np
import pytest

from pseudospace.levelshape import summarize_curve_effects

GRID = np.linspace(0.0, 1.0, 41)


def _curves(*rows):
    return np.vstack(rows)


def test_pure_level_shift_is_all_level():
    reference = _curves(1.0 + 2.0 * GRID)
    comparison = _curves(3.0 + 2.0 * GRID)
    out = summarize_curve_effects(reference, comparison, feature_names=["gene"])
    row = out.iloc[0]
    assert row['level_effect'] == pytest.approx(2.0)
    assert row['level_fraction'] == pytest.approx(1.0)
    assert row['shape_fraction'] == pytest.approx(0.0, abs=1e-12)
    assert row['shape_rms'] == pytest.approx(0.0, abs=1e-12)
    assert row['pattern_rms_z'] == pytest.approx(0.0, abs=1e-12)
    assert row['difference_type'] == 'level shift'


def test_amplitude_loss_is_not_reported_as_a_pattern_change():
    """Same peak position, half the gradient: amplitude moves, pattern does not."""
    reference = _curves(1.0 + 2.0 * GRID)
    comparison = _curves(1.5 + 1.0 * GRID)          # same mean, same peak position, half slope
    out = summarize_curve_effects(reference, comparison).iloc[0]
    assert out['amplitude_ratio'] == pytest.approx(0.5)
    assert out['amplitude_log2_ratio'] == pytest.approx(-1.0)
    assert out['pattern_rms_z'] == pytest.approx(0.0, abs=1e-12)
    assert out['shape_rms'] > 0.1                    # the legacy metric does move
    assert out['difference_type'] == 'amplitude change'


def test_peak_relocation_is_a_pattern_change():
    reference = _curves(np.exp(-((GRID - 0.25) ** 2) / 0.01))
    comparison = _curves(np.exp(-((GRID - 0.75) ** 2) / 0.01))
    out = summarize_curve_effects(reference, comparison).iloc[0]
    assert out['pattern_rms_z'] > 1.0
    assert out['difference_type'] == 'pattern shift'


def test_identical_curves_have_no_detectable_difference():
    reference = _curves(1.0 + 2.0 * GRID)
    out = summarize_curve_effects(reference, reference.copy()).iloc[0]
    assert out['condition_effect_rms'] == pytest.approx(0.0, abs=1e-12)
    assert out['level_effect'] == pytest.approx(0.0, abs=1e-12)
    assert out['difference_type'] == 'no detectable difference'


def test_level_and_shape_fractions_are_orthogonal_and_sum_to_one():
    rng = np.random.default_rng(0)
    reference = rng.normal(size=(6, GRID.size))
    comparison = reference + rng.normal(size=(6, GRID.size))
    out = summarize_curve_effects(reference, comparison)
    assert np.allclose(out['level_fraction'] + out['shape_fraction'], 1.0)
    assert np.allclose(
        out['condition_effect_rms'] ** 2,
        out['level_effect'] ** 2 + out['shape_rms'] ** 2,
    )


def test_nan_outside_support_is_ignored_per_feature():
    reference = _curves(1.0 + 2.0 * GRID, 1.0 + 2.0 * GRID)
    comparison = _curves(2.0 + 2.0 * GRID, 2.0 + 2.0 * GRID)
    comparison[0, :10] = np.nan                      # first feature has no support on 10 cells
    out = summarize_curve_effects(reference, comparison)
    assert out.loc[0, 'n_grid_points'] == GRID.size - 10
    assert out.loc[0, 'level_effect'] == pytest.approx(1.0)
    assert out.loc[1, 'n_grid_points'] == GRID.size


def test_zero_variance_curve_does_not_produce_a_pattern_change():
    reference = _curves(np.full(GRID.size, 2.0))
    comparison = _curves(np.full(GRID.size, 5.0))
    out = summarize_curve_effects(reference, comparison).iloc[0]
    assert np.isnan(out['pattern_rms_z'])
    assert out['difference_type'] == 'level shift'


def test_effect_columns_that_collide_with_a_fit_table_are_identifiable():
    """`shape_rms`/`level_effect` exist in both the fit table and this summary.

    The notebooks therefore merge only the new columns; if this set ever gains a name that a fit
    table already carries, the merge would silently rename rather than overwrite.
    """
    reference = _curves(1.0 + 2.0 * GRID)
    comparison = _curves(3.0 + 2.0 * GRID)
    summary = summarize_curve_effects(reference, comparison, feature_names=['gene'])
    fit_table_columns = {'gene', 'shape_rms', 'level_effect', 'condition_effect_rms',
                        'curve_spearman', 'amplitude_healthy', 'amplitude_aki'}
    overlap = fit_table_columns & set(summary.columns)
    assert overlap == {'shape_rms', 'level_effect', 'condition_effect_rms'}
    merged_new_columns = ['level_fraction', 'shape_fraction', 'pattern_rms_z', 'amplitude_reference',
                         'amplitude_comparison', 'amplitude_ratio', 'amplitude_log2_ratio',
                         'difference_type']
    assert not (set(merged_new_columns) & overlap), 'the notebook merge list must avoid these'
    for column in merged_new_columns:
        assert column in summary.columns
