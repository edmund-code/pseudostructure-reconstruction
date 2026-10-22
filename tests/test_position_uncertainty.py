import numpy as np
import pandas as pd

from pseudospace.position_uncertainty import call_share, half_sample, interval, jackknife_se, spatial_blocks, unit_summary


def _sections(rng, scales=(1.0, 40.0)):
    xy, sec = [], []
    for i, s in enumerate(scales):            # same layout at very different pixel scales
        xy.append(rng.uniform(0, 100, (2000, 2)) * s)
        sec += [f's{i}'] * 2000
    return np.vstack(xy), np.array(sec)


def test_blocks_hit_target_size_at_any_pixel_scale():
    rng = np.random.default_rng(0)
    xy, sec = _sections(rng)
    b = spatial_blocks(xy, sec, per_block=25)
    for s in ('s0', 's1'):
        sizes = pd.Series(b[sec == s]).value_counts()
        assert 20 <= sizes.mean() <= 27
        assert all(name.startswith(f'{s}|') for name in sizes.index)


def test_half_sample_keeps_half_the_blocks_of_every_section_and_whole_blocks():
    rng = np.random.default_rng(1)
    xy, sec = _sections(rng)
    b = spatial_blocks(xy, sec)
    keep = half_sample(b, sec, np.random.default_rng(2))
    for s in ('s0', 's1'):
        n_blocks = len(np.unique(b[sec == s]))
        assert len(np.unique(b[(sec == s) & keep])) == n_blocks // 2
    kept_blocks = set(b[keep])
    assert all(keep[b == name].all() for name in kept_blocks)   # a block is kept whole or not at all


def test_delete_half_jackknife_se_matches_the_se_of_a_mean():
    rng = np.random.default_rng(3)
    x = rng.normal(0, 1, 4000)
    halves = [x[rng.choice(4000, 2000, replace=False)].mean() for _ in range(400)]
    assert abs(jackknife_se(halves) - 1 / np.sqrt(4000)) < .2 / np.sqrt(4000)
    iv = interval(x.mean(), halves)
    assert iv['lo'] < 0 < iv['hi'] and iv['n_draws'] == 400


def test_unit_summary_uses_only_the_draws_that_contain_the_unit():
    pos = pd.DataFrame([[.1, np.nan, .3], [np.nan, np.nan, np.nan], [.5, .5, np.nan]], index=['a', 'b', 'c'])
    s = unit_summary(pos)
    assert s.loc['a', 'n_draws'] == 2 and np.isclose(s.loc['a', 'mean'], .2) and np.isclose(s.loc['a', 'se'], .1)
    assert np.isnan(s.loc['b', 'mean']) and np.isnan(s.loc['b', 'se'])
    assert s.loc['c', 'se'] == 0


def test_call_share_ignores_failed_draws():
    assert call_share([1.0, -1.0, np.nan, 2.0, None], lambda v: v > 0) == 2 / 3
