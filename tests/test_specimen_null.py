import numpy as np
import pandas as pd
import pytest

from pseudospace.pathway_calibration import nested_partial_f
from pseudospace.pathway_remodeling import matched_pathway_tests
from pseudospace.specimen_null import (SEGMENT_COMPARISONS, cell_summaries, matched_auc_null,
                                       pseudobulk_t, segment_designs, summary_partial_f)


def _units(seed=0):
    rng = np.random.default_rng(seed)
    donors = np.repeat(['a', 'b', 'c', 'd'], [150, 90, 120, 70])
    segments = rng.choice(['S1', 'S2', 'S3'], len(donors))
    y = rng.gamma(1., 1., (len(donors), 25)) * (rng.random((len(donors), 25)) > .4)
    y[donors == 'b', :5] += 1.
    y[(donors == 'c') & (segments == 'S3'), 5:10] += 1.
    return y, donors, segments


@pytest.mark.parametrize('nuisance', [None, {'a': 0., 'b': 1., 'c': 0., 'd': 1.}])
def test_summary_partial_f_matches_unit_level_fit(nuisance):
    y, donors, segments = _units()
    group = {'a': 0., 'b': 0., 'c': 1., 'd': 1.}
    designs = segment_designs(donors, segments, group, nuisance)
    weights = np.empty(len(donors))
    for level in (0., 1.):
        names = [k for k, v in group.items() if v == level]
        for name in names:
            weights[donors == name] = len(donors) / (2 * len(names) * (donors == name).sum())
    expected = nested_partial_f(y, designs, weights, SEGMENT_COMPARISONS)
    got = summary_partial_f(cell_summaries(y, donors, segments), list('abcd'), group, nuisance)
    for name in SEGMENT_COMPARISONS:
        np.testing.assert_allclose(got[name], expected[name], rtol=1e-8, atol=1e-10)
        assert got['df'][name] == expected['df'][name]


def test_exact_matched_null_agrees_with_random_sets():
    rng = np.random.default_rng(3)
    genes = [f'g{i}' for i in range(300)]
    stats = pd.DataFrame({'T': rng.gamma(2., 1., 300)}, index=genes)
    strata = rng.integers(0, 4, 300)
    sets = {'p1': genes[:20], 'p2': genes[50:90], 'p3': list(np.array(genes)[strata == 2][:15])}
    exact = matched_auc_null(stats, sets, strata).set_index('pathway_id')
    sampled = matched_pathway_tests(stats, sets, strata, n_null=20000, seed=1).set_index('pathway_id')
    np.testing.assert_allclose(exact.auc, sampled.auc, atol=1e-12)
    np.testing.assert_allclose(exact.null_auc_mean, sampled.null_auc_mean, atol=2e-3)
    np.testing.assert_allclose(exact.null_auc_sd, sampled.null_auc_sd, rtol=.03)


def test_pseudobulk_t_sign_follows_group():
    rng = np.random.default_rng(5)
    counts = rng.poisson(50, (4, 200)).astype(float)
    counts[2:, :10] *= 4
    design = np.column_stack([np.ones(4), [0, 0, 1, 1]])
    t, p = pseudobulk_t(counts, counts.sum(axis=1), design, 1)
    assert (t[:10] > 0).all() and np.median(p[:10]) < np.median(p[10:])


def test_within_group_deviations_remove_group_profiles():
    from pseudospace.specimen_null import donor_deviations
    y, donors, segments = _units(seed=2)
    summary = cell_summaries(y, donors, segments)
    shifted = dict(summary, mean=summary['mean'].copy())
    shifted['mean'][2:] += np.arange(3)[None, :, None] * 5.  # a group-wide segment profile for c, d
    groups = {'a': 0, 'b': 0, 'c': 1, 'd': 1}
    for got, expected in zip(donor_deviations(shifted, list('abcd'), groups),
                             donor_deviations(summary, list('abcd'), groups)):
        np.testing.assert_allclose(got, expected, atol=1e-10)
