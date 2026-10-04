import numpy as np
import pandas as pd
import pytest

from pseudospace.pathway_decoys import (averaging_loss, bootstrap_decoy_calls, decoy_calls, decoy_fdp,
                                        decoy_inflation, joint_matched_test, overlap_programs,
                                        shape_class, welch_t)


def _tests(z, sd=.01, vifs=None):
    ids = [f'p{i}' for i in range(len(z))]
    frame = pd.DataFrame({'pathway_id': ids, 'partition': 'species', 'statistic': 'T_spatial',
                          'null_auc_mean': .5, 'null_auc_sd': sd})
    frame['auc'] = .5 + np.asarray(z) * sd
    return frame


def test_joint_test_reduces_to_matched_z_without_inflation():
    tests = _tests([4., 1., 0.])
    out = joint_matched_test(tests, {p: 1. for p in tests.pathway_id})
    np.testing.assert_allclose(out.z_joint, [4., 1., 0.])
    assert out.q_joint.iloc[0] < .001


def test_joint_test_deflates_and_reserves_missing_inflation():
    tests = _tests([4., 4.])
    out = joint_matched_test(tests, {'p0': 4.})
    assert out.z_joint.iloc[0] == pytest.approx(2.)
    assert out.p_joint.iloc[1] == 1. and np.isnan(out.q_joint.iloc[1])
    with pytest.raises(ValueError):
        joint_matched_test(tests, {'p0': 0., 'p1': 1.})


def test_decoy_fdp_counts_and_monotone_q():
    target = pd.Series([5., 4., 3., 2., 1.], index=list('abcde'))
    decoys = pd.DataFrame({'d1': [0., 0., 3.5, 0., 0.], 'd2': [0., 4.5, 0., 0., 0.]}, index=target.index)
    curve = decoy_fdp(target, decoys)
    row = curve.set_index('threshold').loc[3.]
    assert row.target_calls == 3 and row.decoy_calls == 1.
    assert (np.diff(curve.q.to_numpy()) >= -1e-12).all()       # q never falls as the threshold drops
    assert list(decoy_calls(target, decoys, alpha=.1)) == ['a']
    assert list(decoy_calls(target, decoys, alpha=.25)) == list('abcde')   # FDP(1) = 1/5
    assert len(decoy_calls(target, decoys * 10, alpha=.1)) == 0


def test_decoy_calls_separate_spiked_targets_from_noise():
    rng = np.random.default_rng(0)
    target = pd.Series(np.r_[rng.normal(6, 1, 40), rng.normal(0, 1, 960)])
    decoys = pd.DataFrame({'a': rng.normal(0, 1, 1000), 'b': rng.normal(0, 1, 1000)})
    calls = decoy_calls(target, decoys, alpha=.1)
    assert 35 <= len(calls) <= 60 and len(set(range(40)) & set(calls)) >= 38
    low, mid, high = bootstrap_decoy_calls(target, decoys, n_boot=50)
    assert low <= len(calls) <= high


def test_decoy_inflation_recovers_equicorrelation():
    rng = np.random.default_rng(1)
    sizes = pd.Series(rng.integers(10, 300, 3000))
    rho = .01
    z = pd.DataFrame({d: rng.normal(0, np.sqrt(1 + (sizes - 1) * rho)) for d in ('x', 'y')})
    a, rho_hat, vif = decoy_inflation(z, sizes)
    assert a == pytest.approx(1, abs=.15) and rho_hat == pytest.approx(rho, abs=.003)
    assert (vif >= 1).all()


def test_overlap_programs_groups_shared_drivers():
    sets = {'A': list('abcdefg'), 'B': list('abcdxyz'), 'C': list('pqrstuv'), 'D': list('pqrsmno')}
    out = overlap_programs(sets, list(sets), contributing=set('abcdpqrs'), cut=.5)
    groups = out.groupby('program').pathway_id.apply(set).tolist()
    assert {'A', 'B'} in groups and {'C', 'D'} in groups


def test_averaging_loss_and_shape_classes():
    grid = np.linspace(0, 1, 61)
    loss = averaging_loss(np.vstack([np.ones(61), np.sin(2 * np.pi * grid)]))
    assert loss[0] == pytest.approx(0) and loss[1] == pytest.approx(1, abs=1e-6)
    assert shape_class(np.sin(2 * np.pi * grid)) == 'reversing'
    assert shape_class(np.exp(-((grid - .5) / .05) ** 2)) == 'localized'
    assert shape_class(grid) == 'graded'
    assert shape_class(np.zeros(5)) == 'flat'


def test_welch_t_matches_scipy():
    from scipy.stats import ttest_ind
    rng = np.random.default_rng(2)
    a, b = rng.normal(1, 1, (6, 4)), rng.normal(0, 2, (12, 4))
    t, _ = welch_t(a, b)
    np.testing.assert_allclose(t, ttest_ind(a, b, equal_var=False).statistic)
    with pytest.raises(ValueError):
        welch_t(a[:1], b)
