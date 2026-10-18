import numpy as np
import pytest

from pseudospace import handoffs as ho


def test_profiles_crossing_and_half_points():
    g = np.linspace(0, 1, 101)
    early = ho.scaled_profile(np.log(5) - 4 * g)                 # falls linearly
    late = ho.scaled_profile(2 * g)                               # rises linearly
    assert early[0] == 1 and early[-1] == 0 and np.isnan(ho.scaled_profile([np.nan, 1.])).all()
    assert ho.crossing(early, late, g) == pytest.approx(.5, abs=1e-9)
    fall, rise = ho.half_points(early, late, g)
    assert fall == pytest.approx(.5, abs=1e-9) and rise == pytest.approx(.5, abs=1e-9)
    late_shift = ho.scaled_profile(1 / (1 + np.exp(-(g - .8) / .02)))
    early_step = ho.scaled_profile(1 / (1 + np.exp((g - .3) / .02)))
    fall, rise = ho.half_points(early_step, late_shift, g)
    assert rise - fall == pytest.approx(.5, abs=.02)            # a gap: positive width
    assert np.isnan(ho.crossing(late, early, g))                 # never early > late then below


def test_pair_count_and_null():
    family = np.array(['A'] * 4 + ['B'] * 4 + ['C'] * 2)
    early = np.array([1, 1, 0, 0, 1, 0, 0, 0, 0, 0], bool)
    late = np.array([0, 0, 1, 1, 0, 1, 0, 0, 1, 0], bool)
    assert ho.pair_count(family, early, late) == 2 * 2 + 1 * 1
    null = ho.family_null(family, early, late, n=500, seed=1)
    assert null.min() >= 0 and null.mean() < 5
