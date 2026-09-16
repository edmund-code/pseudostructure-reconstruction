"""Signed gene rankings and competitive, correlation-aware set enrichment."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from pseudospace.enrichment import camera_like_enrichment, score_signed_signatures, signed_gene_rankings

GRID = np.linspace(0.0, 1.0, 51)


def test_signed_rankings_separate_level_amplitude_and_redistribution():
    reference = np.vstack([
        1.0 + 2.0 * GRID,                        # 0: level up, same shape
        1.0 + 2.0 * GRID,                        # 1: amplitude halved
        np.exp(-((GRID - 0.25) ** 2) / 0.01),    # 2: peak moves towards late
        np.exp(-((GRID - 0.25) ** 2) / 0.01),    # 3: supressed late
    ])
    comparison = np.vstack([
        3.0 + 2.0 * GRID,
        1.5 + 1.0 * GRID,
        np.exp(-((GRID - 0.75) ** 2) / 0.01),
        0.3 * np.exp(-((GRID - 0.25) ** 2) / 0.01),
    ])
    out = signed_gene_rankings(reference, comparison, GRID,
                               gene_names=['level', 'amplitude', 'redistribution', 'late'])
    by_gene = out.set_index('gene')
    assert by_gene.loc['level', 'level_effect'] == pytest.approx(2.0)
    # a pure vertical offset has no rank correlation with pseudospace: report it as undefined
    assert np.isnan(by_gene.loc['level', 'redistribution'])
    assert by_gene.loc['amplitude', 'level_effect'] == pytest.approx(0.0, abs=1e-6)
    assert by_gene.loc['amplitude', 'amplitude_log2_ratio'] == pytest.approx(-1.0)
    assert by_gene.loc['redistribution', 'redistribution'] > 0.5          # towards late pseudospace
    assert by_gene.loc['redistribution', 'early_delta'] < by_gene.loc['redistribution', 'late_delta']
    # suppressing an early peak shows up as a negative *early* contrast, not a late one
    assert by_gene.loc['late', 'early_delta'] < by_gene.loc['late', 'late_delta']
    assert by_gene.loc['late', 'early_delta'] < -0.3


def test_competitive_test_contrasts_with_the_background_not_zero():
    """Every gene shifted up by the same amount must NOT be significant: the null is the background."""
    rng = np.random.default_rng(0)
    names = [f'g{i}' for i in range(200)]
    # every gene carries the same offset: the raw mean is large, the deviation from the background is not
    shifted = {name: 5.0 + float(value) for name, value in zip(names, rng.normal(size=200))}
    sets = {'half': names[:100]}
    result = camera_like_enrichment(shifted, sets, min_set_size=5)
    row = result.set_index('gene_set').loc['half']
    assert row['set_mean_statistic'] > 4.5                       # raw mean is dominated by the offset
    assert abs(row['set_mean_statistic_centred']) < 0.5
    assert row['p_value_competitive_normal'] > 0.05, 'no deviation from the background mean'

    elevated = dict(shifted)
    elevated.update({name: shifted[name] + 4.0 for name in names[:100]})
    result = camera_like_enrichment(elevated, sets, min_set_size=5)
    row = result.set_index('gene_set').loc['half']
    # the background mean includes the elevated members, so the centred effect is ~2, not exactly 2
    assert row['set_mean_statistic_centred'] == pytest.approx(2.0, abs=0.3)
    assert row['p_value_competitive_normal'] < 1e-3


def test_permutation_validation_agrees_with_the_normal_approximation():
    """The analytic competitive p must match a brute-force sampling-without-replacement null."""
    rng = np.random.default_rng(7)
    names = [f'g{i}' for i in range(400)]
    statistics = {name: float(value) for name, value in zip(names, rng.normal(size=400))}
    member_names = names[:12]
    observed = np.mean([statistics[name] for name in member_names])
    all_values = np.array([statistics[name] for name in names])
    background_mean = all_values.mean()
    draws = 20_000
    empirical = np.array([
        all_values[rng.choice(len(all_values), size=12, replace=False)].mean() - background_mean
        for _ in range(draws)
    ])
    empirical_p = float(np.mean(np.abs(empirical) >= abs(observed - background_mean)))
    result = camera_like_enrichment(statistics, {'set': member_names}, min_set_size=5)
    analytic_p = float(result['p_value_competitive_normal'].iloc[0])
    assert 0.0 < empirical_p <= 1.0 and 0.0 < analytic_p <= 1.0
    assert 0.4 < analytic_p / max(empirical_p, 1e-4) < 2.5, (analytic_p, empirical_p)


def test_inflation_uses_residual_correlation_when_a_design_is_given():
    rng = np.random.default_rng(0)
    n_structures, n_genes = 400, 30
    trend = np.linspace(-1, 1, n_structures)[:, None]
    block = trend @ rng.normal(size=(1, 15)) + 0.2 * rng.normal(size=(n_structures, 15))
    expression = np.column_stack([block, rng.normal(size=(n_structures, n_genes - 15))])
    names = [f'g{i}' for i in range(n_genes)]
    statistics = {name: (1.0 if index < 15 else 0.0) for index, name in enumerate(names)}
    design = np.column_stack([np.ones(n_structures), trend.ravel()])
    raw = camera_like_enrichment(statistics, {'correlated': names[:15]}, expression=expression,
                                 gene_names=names, min_set_size=5)
    residual = camera_like_enrichment(statistics, {'correlated': names[:15]}, expression=expression,
                                      gene_names=names, design=design, min_set_size=5)
    raw_inflation = float(raw['variance_inflation'].iloc[0])
    residual_inflation = float(residual['variance_inflation'].iloc[0])
    # the shared trend is what makes the genes look correlated; removing it must lower the inflation
    assert raw_inflation > residual_inflation > 1.0
    assert residual['correlation_basis'].iloc[0] == 'residual'
    assert raw['correlation_basis'].iloc[0] == 'mean-centred'


def test_small_or_out_of_background_sets_are_not_tested():
    statistics = {f'g{i}': float(i) for i in range(30)}
    sets = {'too_small': ['g0', 'g1'], 'outside_background': ['g26', 'g27', 'g28'],
            'partly_inside': ['g20', 'g21', 'g22', 'g23', 'g24']}
    result = camera_like_enrichment(statistics, sets, background=[f'g{i}' for i in range(25)],
                                    min_set_size=5)
    by_set = result.set_index('gene_set')
    assert np.isnan(by_set.loc['too_small', 'p_value_competitive_normal'])
    assert np.isnan(by_set.loc['outside_background', 'p_value_competitive_normal'])
    assert by_set.loc['partly_inside', 'n_genes_tested'] == 5
    assert 'p_value_competitive_normal_adjusted' in result.columns


def test_signature_scores_use_weights_and_report_coverage():
    expression = np.array([
        [0.0, 10.0, 1.0], [1.0, 10.0, 1.0], [2.0, 10.0, 1.0], [3.0, 10.0, 1.0],
    ])
    signatures = {
        'up': {'g1': 1.0, 'g2': 1.0},
        'down': {'g1': -1.0, 'g2': 1.0},
        'absent': {'gX': 1.0},
    }
    scores, coverage = score_signed_signatures(expression, signatures,
                                               gene_names=['g1', 'g2', 'g3'], min_genes=2)
    assert list(scores.columns) == ['up', 'down']
    assert scores['up'].iloc[-1] > scores['up'].iloc[0]        # g1 rises, so 'up' rises
    assert scores['down'].iloc[-1] < scores['down'].iloc[0]
    assert coverage.set_index('signature').loc['absent', 'usable'] == False  # noqa: E712
    assert coverage.set_index('signature').loc['up', 'n_members_present'] == 2
