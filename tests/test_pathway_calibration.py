import numpy as np
import pandas as pd
import pytest

from pseudospace.pathway_calibration import (balanced_partitions, contrast_designs, nested_partial_f,
                                             negative_control_summary, normal_matched_q, parametric_calls,
                                             pathway_scores)
from pseudospace.stats_gam import gam_internal_knots


def _cohort(seed=3, per=150):
    rng = np.random.default_rng(seed)
    specimen = np.repeat(['m1', 'm2', 'h1', 'h2'], per)
    human = np.isin(specimen, ['h1', 'h2']).astype(float)
    position = rng.uniform(0, 1, len(specimen))
    return rng, specimen, human, position


def test_balanced_partitions_lists_species_split_and_two_mixed_swaps():
    specimen = np.array(['m1', 'm2', 'h1', 'h2'])
    parts = balanced_partitions(specimen, ['mouse', 'mouse', 'human', 'human'])
    assert parts['species'] == ('m1', 'm2')  # 'mouse' sorts after 'human': group 1
    swaps = [v for k, v in parts.items() if k.startswith('swap:')]
    assert len(swaps) == 2 and all(len({s[0] for s in pair}) == 2 for pair in swaps)
    with pytest.raises(ValueError):
        balanced_partitions(specimen[:3], ['mouse', 'mouse', 'human'])


def test_relabeled_shape_test_ignores_species_shape_and_offsets():
    rng, specimen, human, position = _cohort()
    offsets = {'m1': 0., 'm2': .8, 'h1': 1.5, 'h2': 2.4}       # specimen offsets only
    y_species = np.sin(3 * position) + human * np.cos(4 * position)  # true species shape
    y = y_species + np.vectorize(offsets.get)(specimen)
    y = np.column_stack([y + rng.normal(0, .3, len(y)) for _ in range(40)])
    knots = gam_internal_knots(position, basis_df=6)
    swap = np.isin(specimen, ['m1', 'h1']).astype(float)
    designs, w = contrast_designs(position, swap, specimen, knots, nuisance_shape=human)
    scores = nested_partial_f(y, designs, w, {'T_spatial': ('level', 'full')})
    # Under the relabeling there is no group shape: F stays near its null mean of ~1.
    assert np.median(scores['T_spatial']) < 2
    designs, w = contrast_designs(position, human, specimen, knots)
    true = nested_partial_f(y, designs, w, {'T_spatial': ('level', 'full')})
    assert np.median(true['T_spatial']) > 20


def test_nested_partial_f_matches_direct_weighted_least_squares_and_rejects_nonnesting():
    rng, specimen, human, position = _cohort(per=40)
    knots = gam_internal_knots(position, basis_df=4)
    labels = np.where(position < .4, 'S1', np.where(position < .7, 'S2', 'S3'))
    designs, w = contrast_designs(position, human, specimen, knots, segments=labels)
    y = rng.normal(size=(len(position), 3)) + human[:, None] * position[:, None]
    got = nested_partial_f(y, designs, w)

    def sse(x, col):
        r = np.sqrt(w)
        beta = np.linalg.lstsq(r[:, None] * x, r * y[:, col], rcond=None)[0]
        res = r * (y[:, col] - x @ beta)
        return res @ res, np.linalg.matrix_rank(r[:, None] * x)

    for name, (reduced, full) in [('T_spatial', ('level', 'full')),
                                  ('T_position_given_segments', ('discrete', 'union'))]:
        for col in range(3):
            s0, r0 = sse(designs[reduced], col)
            s1, r1 = sse(designs[full], col)
            expect = (s0 - s1) / (r1 - r0) / (s1 / (len(y) - r1))
            assert got[name][col] == pytest.approx(expect, rel=1e-8)
    with pytest.raises(ValueError):
        nested_partial_f(y, designs, w, {'bad': ('full', 'discrete')})


def test_normal_matched_q_and_summary_count_relabeled_discoveries():
    tests = pd.DataFrame({'partition': ['species'] * 3 + ['swap:a'] * 3, 'statistic': 'T',
                          'auc': [.9, .8, .5, .5, .52, .9], 'null_auc_mean': .5, 'null_auc_sd': .02})
    tests = normal_matched_q(tests)
    assert tests.loc[0, 'q_normal'] < 1e-6 and tests.loc[2, 'q_normal'] > .3
    hits = pd.DataFrame({'strategy': 'T', 'partition': tests.partition,
                         'pathway_id': ['p1', 'p2', 'p3', 'p1', 'p2', 'p3'],
                         'hit': tests.q_normal.le(.05)})
    summary = negative_control_summary(hits).loc['T']
    assert summary.true_hits == 2 and summary.null_mean == 1 and summary.ncdr == .5
    assert summary.true_hits_also_in_null == 0 and summary.true_only == 2


def test_pathway_scores_and_parametric_calls_shapes():
    rng = np.random.default_rng(0)
    x = rng.normal(size=(50, 4))
    scores = pathway_scores(x, {'a': ['g0', 'g1'], 'b': ['g2', 'g3']}, ['g0', 'g1', 'g2', 'g3'],
                            np.ones(50))
    assert scores.shape == (50, 2) and abs(scores.a.mean()) < 1e-12
    calls = parametric_calls({'T': pd.Series([0.1, 50.], index=['a', 'b'])}, {'T': (2, 100)})
    assert calls.hit.tolist() == [False, True]
    with pytest.raises(ValueError):
        pathway_scores(x, {'a': ['missing']}, ['g0', 'g1', 'g2', 'g3'], np.ones(50))
