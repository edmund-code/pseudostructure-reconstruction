import numpy as np
import pandas as pd
import pytest

from pseudospace.continuous_descriptors import (blockwise_permutation, bootstrap_partial_spearman,
                                                coordinate_segments, label_agreeing_thresholds,
                                                monotone_onsets, partial_spearman, segment_means,
                                                segment_units, specimen_curves, switch_fraction)
from pseudospace.stats_gam import gam_internal_knots


def _labels(s, t12=.3, t23=.65):
    return np.where(s < t12, 'PT-S1', np.where(s < t23, 'PT-S2', 'PT-S3'))


def test_thresholds_recover_cuts_despite_label_noise():
    rng = np.random.default_rng(0)
    s = rng.random(3000)
    lab = _labels(s)
    flip = rng.random(len(s)) < .08
    lab[flip] = rng.choice(['PT-S1', 'PT-S2', 'PT-S3'], flip.sum())
    t12, t23, agree = label_agreeing_thresholds(s, lab)
    assert abs(t12 - .3) < .03 and abs(t23 - .65) < .03
    assert .9 < agree < .97


def test_coordinate_segments_are_per_specimen():
    rng = np.random.default_rng(1)
    s = rng.random(2000)
    spec = np.repeat(['a', 'b'], 1000)
    lab = np.where(spec == 'a', _labels(s, .2, .5), _labels(s, .4, .8))
    seg, table = coordinate_segments(s, lab, spec)
    assert (seg == lab).mean() > .99
    assert table.set_index('specimen').loc['b', 't12'] == pytest.approx(.4, abs=.02)


def test_segment_units_and_blockwise_permutation():
    u = segment_units(np.array([0., .15, .3, .65, 1.]), .3, .65)
    assert np.allclose(u, [0, .5, 1, 2, 3])
    blocks = np.repeat(['x', 'y', 'z'], [5, 7, 3])
    index = blockwise_permutation(blocks, np.random.default_rng(2))
    assert sorted(index) == list(range(15))
    assert (blocks[index] == blocks).all()


def test_monotone_onsets_recover_logistic_midpoint_and_width():
    u = np.linspace(0, 3, 121)
    curves = np.vstack([1 / (1 + np.exp(-(u - 1.4) / .1)),             # rising, sharp, onset 1.4
                        -2 / (1 + np.exp(-(u - 1.8) / .2)),            # falling, onset 1.8
                        .1 * np.sin(u),                                # small amplitude
                        np.exp(-(u - 1.5) ** 2 / .1)])                # bump: not monotone
    table = monotone_onsets(curves, u)
    assert list(table.monotone) == [True, True, False, False]
    assert table.onset[0] == pytest.approx(1.4, abs=.02)
    assert table.onset[1] == pytest.approx(1.8, abs=.02)
    assert table.width[0] == pytest.approx(2 * np.log(9) * .1, abs=.05)
    assert table.direction[1] == -1


def test_switch_fraction_orders_early_and_late_switches():
    early, late = [0., .9, 1.], [0., .1, 1.]
    falling = [1., .2, 0.]
    r = switch_fraction(np.array([early, late, falling, [1., 1., 1.02]]))
    assert r[0] > r[1] and r[2] == pytest.approx(.8)
    assert np.isnan(r[3])


def test_partial_spearman_removes_shared_covariate():
    rng = np.random.default_rng(3)
    z = rng.normal(size=500)
    x = z + .3 * rng.normal(size=500)
    y = z + .3 * rng.normal(size=500)
    assert abs(partial_spearman(x, y, z)) < .15
    w = rng.normal(size=500)
    est, lo, hi = bootstrap_partial_spearman(x + w, y + w, z, n_boot=200, seed=0)
    assert est > .3 and lo > 0 and lo <= est <= hi


def test_specimen_curves_and_segment_means():
    rng = np.random.default_rng(4)
    s = rng.random(4000)
    spec = np.repeat(['a', 'b'], 2000)
    truth = np.column_stack([np.sin(3 * s), s ** 2])
    y = truth + .05 * rng.normal(size=truth.shape)
    grid = np.linspace(.05, .95, 31)
    curves = specimen_curves(y, s, spec, grid, gam_internal_knots(s, 6))
    assert set(curves) == {'a', 'b'}
    assert np.abs(curves['a'][0] - np.sin(3 * grid)).max() < .03
    means = segment_means(y, _labels(s))
    assert means.shape == (2, 3) and means[1, 0] < means[1, 1] < means[1, 2]


def test_read_coordinate_and_common_support(tmp_path):
    from pseudospace.coordinate_inputs import common_support_mask, read_coordinate

    path = tmp_path / 'coord.csv'
    pd.DataFrame({'structure_id': ['a', 'b', 'c'], 'position': [.1, .5, np.nan],
                  'other': [1., 2., 3.]}).to_csv(path, index=False)
    assert read_coordinate(path).to_dict() == {'a': .1, 'b': .5}
    assert read_coordinate(path, 'other')['c'] == 3.
    with pytest.raises(ValueError):
        read_coordinate(path, 'missing')
    s = np.r_[np.linspace(0, 1, 101), np.linspace(.2, .8, 101)]
    spec = np.repeat(['x', 'y'], 101)
    keep = common_support_mask(s, spec)
    assert s[keep].min() >= .2 and s[keep].max() <= .8 and keep[101:].mean() > .97


def test_group_difference_curves_recover_shape():
    from pseudospace.continuous_descriptors import group_difference_curves

    rng = np.random.default_rng(5)
    n = 4000
    s = rng.random(n)
    spec = np.repeat(['m1', 'm2', 'h1', 'h2'], n // 4)
    g = np.isin(spec, ['h1', 'h2']).astype(float)
    truth_delta = .5 + np.sin(4 * s)
    y = np.column_stack([np.cos(3 * s) + g * truth_delta + .05 * rng.normal(size=n),
                         s + .05 * rng.normal(size=n)])
    grid = np.linspace(.05, .95, 19)
    knots = gam_internal_knots(s, 6)
    delta = group_difference_curves(y, s, g, spec, knots, grid)
    assert delta.shape == (2, 19)
    assert np.abs(delta[0] - (.5 + np.sin(4 * grid))).max() < .1
    assert np.abs(delta[1]).max() < .05


def test_within_segment_slopes_and_interior_extrema():
    from pseudospace.continuous_descriptors import interior_extrema, within_segment_slopes

    rng = np.random.default_rng(6)
    s = rng.random(3000)
    lab = _labels(s)
    y = np.column_stack([s + .05 * rng.normal(size=3000),                       # rising within every segment
                         (lab == 'PT-S2') * 1. + .05 * rng.normal(size=3000)])   # pure step: flat inside segments
    t = within_segment_slopes(y, s, np.repeat('a', 3000), lab)['a']
    assert (t.iloc[0] > 10).all() and (t.iloc[1].abs() < 4).all()
    u = np.linspace(0, 3, 121)
    curves = np.vstack([np.exp(-(u - 1.5) ** 2 / .02),          # peak in the middle of S2
                        np.exp(-(u - 1.02) ** 2 / .02),         # peak at the S1/S2 boundary: not interior
                        -np.exp(-(u - 2.5) ** 2 / .02)])        # trough inside S3
    seg, sign = interior_extrema(curves, u)
    assert list(seg) == [1, -1, 2] and list(sign) == [1, 0, -1]
