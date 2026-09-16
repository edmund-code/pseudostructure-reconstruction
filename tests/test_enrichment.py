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


def test_camera_like_reports_a_competitive_p_and_a_correlation_inflation():
    rng = np.random.default_rng(0)
    n_structures, n_genes = 200, 60
    latent = rng.normal(size=(n_structures, 1))
    block = latent + 0.1 * rng.normal(size=(n_structures, 20))        # 20 correlated genes
    noise = rng.normal(size=(n_structures, n_genes - 20))
    expression = np.column_stack([block, noise])
    names = [f'g{i}' for i in range(n_genes)]
    signal = rng.normal(1.0, 0.1, size=20)          # the correlated block carries the signal
    statistics = {name: float(value) for name, value in
                  zip(names, np.concatenate([signal, np.zeros(n_genes - 20)]))}
    result = camera_like_enrichment(statistics, {'correlated_up': names[:20], 'noise': names[20:]},
                                    expression=expression, gene_names=names, n_permutations=200,
                                    seed=1, min_set_size=5)
    by_set = result.set_index('gene_set')
    assert by_set.loc['correlated_up', 'set_mean_statistic'] > 0.9
    assert by_set.loc['correlated_up', 'p_value_permutation'] < 0.05
    # the correlated block must inflate the variance of the set mean
    assert by_set.loc['correlated_up', 'correlation_inflation'] > 1.0
    assert np.isfinite(by_set.loc['correlated_up', 'p_value_correlation_aware'])
    assert by_set.loc['noise', 'set_mean_statistic'] == pytest.approx(0.0)
    assert by_set.loc['noise', 'p_value_permutation'] > 0.5
    assert 'p_value_permutation_adjusted' in result.columns


def test_correlation_estimate_is_capped_without_losing_the_inflation():
    """The capped estimate is a noisier version of the same number, not a different regime."""
    rng = np.random.default_rng(4)
    n_correlated = 30
    latent = rng.normal(size=(400, 1))
    expression = np.column_stack([latent + 0.1 * rng.normal(size=(400, n_correlated)),
                                  rng.normal(size=(400, 20))])
    names = [f'g{i}' for i in range(expression.shape[1])]
    statistics = {name: (1.0 if index < n_correlated else 0.0)
                  for index, name in enumerate(names)}
    gene_sets = {'correlated': names[:n_correlated]}
    capped = camera_like_enrichment(statistics, gene_sets, expression=expression, gene_names=names,
                                    n_permutations=50, min_set_size=5,
                                    correlation_gene_cap=15, correlation_structure_cap=200)
    uncapped = camera_like_enrichment(statistics, gene_sets, expression=expression, gene_names=names,
                                      n_permutations=50, min_set_size=5)
    capped_inflation = float(capped['correlation_inflation'].iloc[0])
    uncapped_inflation = float(uncapped['correlation_inflation'].iloc[0])
    assert capped_inflation > 1.5 and uncapped_inflation > 1.5
    assert 0.5 < capped_inflation / uncapped_inflation < 2.0


def test_camera_like_respects_background_and_minimum_size():
    statistics = {f'g{i}': float(i) for i in range(30)}
    sets = {
        'too_small': ['g0', 'g1'],
        'outside_background': ['g26', 'g27', 'g28', 'g29'],
        'partly_inside': ['g20', 'g21', 'g22', 'g23', 'g24', 'g25', 'g26', 'g27'],
    }
    result = camera_like_enrichment(statistics, sets, background=[f'g{i}' for i in range(25)],
                                    n_permutations=100, min_set_size=5)
    by_set = result.set_index('gene_set')
    assert np.isnan(by_set.loc['too_small', 'p_value_permutation'])
    assert np.isnan(by_set.loc['outside_background', 'p_value_permutation'])
    # only the members that sit in the background window are tested
    assert by_set.loc['partly_inside', 'n_genes_tested'] == 5


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
