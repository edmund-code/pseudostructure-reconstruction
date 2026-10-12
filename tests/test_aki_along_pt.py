"""aki_along_pt on synthetic data: segment units, centred program profiles, summaries, calibration, placement."""
import numpy as np
import pandas as pd
import pytest

from pseudospace.aki_along_pt import (balanced_cut, centre_of_mass, centre_rows, empirical_p, equal_weight_quantiles,
                                      knn_place, profile_summary, program_profile, segment_cuts, segment_units)


def test_balanced_cut_separates_two_groups():
    x = np.r_[np.linspace(0, .4, 50), np.linspace(.6, 1, 10)]
    assert balanced_cut(x, x > .5) == pytest.approx(.5)
    with pytest.raises(ValueError):
        balanced_cut(x, np.zeros_like(x, bool))


def test_segment_cuts_and_units():
    rng = np.random.default_rng(0)
    x = rng.uniform(0, 1, 3000)
    seg = np.where(x < .3, 'PT-S1', np.where(x < .7, 'PT-S2', 'PT-S3'))
    t12, t23 = segment_cuts(x, seg)
    assert t12 == pytest.approx(.3, abs=.01) and t23 == pytest.approx(.7, abs=.01)
    s = segment_units(np.array([0, t12, (t12 + t23) / 2, t23, 1, 1.15]), t12, t23)
    np.testing.assert_allclose(s[:5], [0, 1, 1.5, 2, 3])
    assert s[5] == pytest.approx(3 + .15 / (1 - t23))          # linear extension beyond hi
    with pytest.raises(ValueError):
        segment_units(x, .7, .3)


def test_centring_removes_a_shared_shift_and_program_median_is_log2():
    grid = 11
    shift = np.linspace(0, 1, grid)                            # a composition tilt shared by every gene
    delta = np.tile(shift, (9, 1))
    delta[:3] += np.log(4)                                     # three genes rise fourfold everywhere
    centred, m = centre_rows(delta, np.ones(9, bool))
    np.testing.assert_allclose(m, shift)
    profile, n = program_profile(centred, [f'g{i}' for i in range(9)], ['g0', 'g1', 'g2', 'absent'])
    assert n == 3
    np.testing.assert_allclose(profile, 2.)                    # log2(4)
    _, n_usable = program_profile(centred, [f'g{i}' for i in range(9)], ['g0', 'g1'], usable=np.r_[False, np.ones(8, bool)])
    assert n_usable == 1


def test_profile_summary_level_drift_centre_and_segments():
    s = np.linspace(0, 3, 101)
    ramp = s / 3                                               # rises toward S3
    out = profile_summary(ramp, s)
    assert out['level'] == pytest.approx(.5)
    assert out['drift'] == pytest.approx(ramp[-21:].mean() - ramp[:21].mean())
    assert out['peak'] == pytest.approx(3)
    assert out['centre_of_mass'] == pytest.approx((s * ramp).sum() / ramp.sum())
    assert out['within_S3'] == pytest.approx(1 / 3, abs=.02)
    bump = np.exp(-(s - 2.5) ** 2 / .05)                       # activation confined to late S3
    assert centre_of_mass(bump, s) == pytest.approx(2.5, abs=.01)
    assert np.isnan(centre_of_mass(-bump, s))


def test_empirical_p_counts_the_null_at_or_below():
    null = np.array([.01, .02, .5, .9])
    np.testing.assert_allclose(empirical_p([.001, .02, 1.], null), [1 / 5, 3 / 5, 5 / 5])


def test_equal_weight_quantiles_ignore_group_size():
    values = np.r_[np.zeros(1000), np.ones(10)]                # a large group at 0 and a small one at 1
    groups = np.r_[np.repeat('a', 1000), np.repeat('b', 10)]
    q = equal_weight_quantiles(values, groups, [.25, .75])
    assert q[0] == pytest.approx(0) and q[1] == pytest.approx(1)


def test_knn_place_reports_median_and_label_shares():
    ref = np.r_[np.zeros((20, 2)), np.full((20, 2), 10.)]
    values = np.r_[np.full(20, .1), np.full(20, .9)]
    labels = np.r_[np.repeat('PT-S1', 20), np.repeat('cTAL', 20)]
    out = knn_place(np.array([[0.1, 0.], [9.9, 10.]]), ref, values, labels, k=5)
    assert list(out.majority) == ['PT-S1', 'cTAL']
    np.testing.assert_allclose(out.median_value, [.1, .9])
    np.testing.assert_allclose(out['share_PT-S1'], [1, 0])
    with pytest.raises(ValueError):
        knn_place(np.zeros((1, 3)), ref, values, labels)
