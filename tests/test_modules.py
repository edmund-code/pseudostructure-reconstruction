"""Curve-module discovery: shape-based grouping, descriptors, stability and enrichment."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from pseudospace.modules import (
    curve_descriptors,
    difference_curves,
    discover_curve_modules,
    enrich_modules,
    module_stability,
)

GRID = np.linspace(0.0, 1.0, 51)


def _early_declining(n=4, scale=1.0):
    return [scale * (2.5 - 2.0 * GRID) for _ in range(n)]


def _late_rising(n=4, scale=1.0):
    return [scale * (0.5 + 2.0 * GRID) for _ in range(n)]


def test_difference_curves_is_signed_and_centerable():
    reference = np.vstack([1.0 + GRID])
    comparison = np.vstack([3.0 + GRID])
    delta = difference_curves(reference, comparison)
    assert np.allclose(delta, 2.0)
    centered = difference_curves(reference, comparison, center=True)
    assert np.allclose(centered, 0.0)


def test_curve_descriptors_label_known_shapes():
    curves = np.vstack([
        2.5 - 2.0 * GRID,                                     # declining
        0.5 + 2.0 * GRID,                                     # rising
        np.exp(-((GRID - 0.5) ** 2) / 0.01),                  # mid peak
        np.full(GRID.size, 2.0),                              # flat
        np.exp(-((GRID - 0.95) ** 2) / 0.005),                # late peak
    ])
    described = curve_descriptors(curves, GRID)
    assert described.loc[0, 'direction'] == 'early-declining'
    assert described.loc[1, 'direction'] == 'late-rising'
    assert described.loc[2, 'direction'] == 'mid-peaked'
    assert described.loc[3, 'direction'] == 'flat'
    assert described.loc[4, 'peak_position'] == pytest.approx(0.95, abs=0.03)
    assert described.loc[2, 'amplitude'] == pytest.approx(1.0, abs=0.02)
    assert described.loc[0, 'monotonicity'] < -0.9
    assert described.loc[1, 'monotonicity'] > 0.9


def test_modules_group_by_shape_not_by_amplitude_or_level():
    curves = np.vstack(
        _early_declining(3, scale=1.0)      # same shape, different scale ...
        + _early_declining(2, scale=2.5)    # ... and a level offset
        + _late_rising(4)
        + [np.full(GRID.size, 1.0)]         # flat: unassigned
    )
    curves[5] += 4.0                        # level offset on one declining curve
    labels, table = discover_curve_modules(curves, GRID, max_distance=0.4, min_features=3)
    assert set(labels.iloc[:5]) == {'M1'}
    assert set(labels.iloc[5:9]) == {'M2'}
    assert labels.iloc[9] == 'unassigned'
    assert set(table['module']) == {'M1', 'M2'}
    assert table.set_index('module').loc['M1', 'n_features'] == 5


def test_small_and_too_flat_curves_stay_unassigned():
    curves = np.vstack(_early_declining(3) + _late_rising(2) + [2.0 + 0.001 * GRID])
    labels, _ = discover_curve_modules(curves, GRID, max_distance=0.5, min_features=3,
                                       min_amplitude=0.01)
    assert set(labels.iloc[:3]) == {'M1'}
    assert set(labels.iloc[3:5]) == {'unassigned'}, 'a 2-member cluster is not a module'
    assert labels.iloc[5] == 'unassigned'


def test_n_modules_cut_is_honoured():
    curves = np.vstack(_early_declining(3) + _late_rising(3) + [np.exp(-((GRID - 0.5) ** 2) / 0.01)] * 3)
    labels, table = discover_curve_modules(curves, GRID, n_modules=3, min_features=2)
    assert len(set(labels)) == 3
    assert len(table) == 3


def test_module_stability_reports_agreement_and_disorder():
    base = pd.Series(['M1'] * 5 + ['M2'] * 5, index=[f'g{i}' for i in range(10)])
    same = base.copy()
    renamed = pd.Series(['M2'] * 5 + ['M1'] * 5, index=base.index)   # same partition, other names
    repartitioned = pd.Series(['M1', 'M2'] * 5, index=base.index)    # different membership
    moved = base.copy()
    moved.iloc[0] = 'M2'
    stability = module_stability({'baseline': base, 'identical': same, 'renamed': renamed,
                                  'repartitioned': repartitioned, 'one_feature_moved': moved})
    table = stability.set_index('run')
    assert table.loc['identical', 'adjusted_rand_index_vs_reference'] == pytest.approx(1.0)
    # module names are arbitrary, so an exact relabelling must still count as full agreement
    assert table.loc['renamed', 'adjusted_rand_index_vs_reference'] == pytest.approx(1.0)
    assert table.loc['repartitioned', 'adjusted_rand_index_vs_reference'] < 0.0
    assert table.loc['one_feature_moved', 'fraction_same_module'] == pytest.approx(0.9)


def test_enrichment_uses_the_discovery_background_and_corrects_across_pairs():
    background = [f'bg{i}' for i in range(100)]
    modules = {'M1': background[:10], 'M2': background[10:20]}
    gene_sets = {'hits_M1': background[:8] + ['outside'], 'hits_none': ['outside']}
    result = enrich_modules(modules, gene_sets, background=background)
    best = result[(result['module'] == 'M1') & (result['gene_set'] == 'hits_M1')].iloc[0]
    assert best['n_module_genes'] == 10
    assert best['n_set_genes_in_background'] == 8      # 'outside' is not in the background
    assert best['n_overlap'] == 8
    assert best['p_value'] < 1e-6
    assert best['p_value_adjusted'] >= best['p_value']
    none = result[(result['module'] == 'M1') & (result['gene_set'] == 'hits_none')].iloc[0]
    assert np.isnan(none['p_value']), 'a set with no background genes cannot be tested'
    assert len(result) == 4
