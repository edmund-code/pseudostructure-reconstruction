import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from pseudospace.zonation_amplitude import (amplitude_ratio, bootstrap_ratio, deming_slope, flattening_class,
                                            group_sums, log_ratio, ols_slope)


def _simulate(n=4000, ratio=.4, seed=0, ref_noise=1., our_noise=.5):
    rng = np.random.default_rng(seed)
    truth = rng.normal(0, 1, n)                                   # true mouse gradient
    mouse = truth + rng.normal(0, our_noise, n)
    human = ratio * truth + rng.normal(0, our_noise, n)
    reference = truth + rng.normal(0, ref_noise, n)                # independent noisy reference
    return truth, human, mouse, reference


def test_reference_ratio_is_dilution_free_while_naive_regression_is_not():
    _, human, mouse, reference = _simulate(ref_noise=1.5)
    out = amplitude_ratio(human, mouse, reference)
    assert out['ratio'] == pytest.approx(.4, abs=.05)
    assert out['slope_mouse'] < .5                                 # each slope is diluted ...
    naive = ols_slope(mouse, human)                                # ... and regressing on our mouse is biased
    assert naive < .35


def test_deming_recovers_slope_with_known_error_ratio():
    rng = np.random.default_rng(1)
    x_true = rng.normal(0, 1, 20000)
    x = x_true + rng.normal(0, .8, x_true.size)
    y = 2 * x_true + rng.normal(0, .4, x_true.size)
    assert deming_slope(x, y, (.4 / .8) ** 2) == pytest.approx(2, abs=.05)
    assert ols_slope(x, y) < 1.6
    assert deming_slope(x, y, 1e9) == pytest.approx(ols_slope(x, y), rel=1e-3)


def test_orthogonal_and_deming_methods_need_inputs():
    _, human, mouse, reference = _simulate()
    with pytest.raises(ValueError):
        amplitude_ratio(human, mouse, reference, method='deming')
    assert np.isfinite(amplitude_ratio(human, mouse, reference, method='orthogonal')['ratio'])


def test_bootstrap_interval_covers_truth_and_resamples_donors():
    truth, human, mouse, _ = _simulate(n=1500, seed=2)
    genes = pd.Index([f'g{i}' for i in range(len(truth))])
    h, m = pd.Series(human, index=genes), pd.Series(mouse, index=genes)
    donors = truth[None, :] + np.random.default_rng(3).normal(0, 2, (6, len(truth)))

    def draw(rng):
        rows = np.arange(6) if rng is None else rng.integers(0, 6, 6)
        return pd.Series(donors[rows].mean(axis=0), index=genes)

    (low, mid, high), draws = bootstrap_ratio(h, m, draw, n_boot=300, seed=4)
    assert low < .4 < high and len(draws) == 300


def test_group_sums_and_log_ratio():
    counts = sparse.csr_matrix(np.array([[1, 0], [3, 2], [0, 4]], float))
    names, sums = group_sums(counts, ['b', 'a', 'b'])
    assert list(names) == ['a', 'b'] and sums.tolist() == [[3, 2], [1, 4]]
    assert log_ratio(np.array([99.5]), 1e6 - 1, np.array([49.5]), 1e6 - 1)[0] == pytest.approx(1)


def test_flattening_classes():
    assert flattening_class([.2, .25], [.5, .6], .4)[0] == 'consistent with flattening'
    assert flattening_class([1., 1.2], [-.8, -1.], .4)[0] == 'reversal'
    assert flattening_class([1., 1.1], [.1, .05], .4)[0] == 'human-specific gradient'
    assert flattening_class([2., 2.2], [1., 1.1], .4)[0] == 'human-stronger'
    assert flattening_class([1., -1.], [.1, .1], .4)[0] == 'consistent with flattening'   # pairings disagree
    assert flattening_class([np.nan, 1.], [1., 1.], .4)[0] == 'not testable'
