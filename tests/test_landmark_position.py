"""Synthetic checks for the landmark count-position model (no private data)."""
import numpy as np
import pytest
from scipy.stats import spearmanr

from pseudospace.landmark_position import (_balanced_aggregates, deviance_gain, fit_lcp, grid_points,
                                           landmark_ratio_score, piecewise_shape, predictive_loglik,
                                           project_lcp, residual_spline_gain, scaled_gain, spline_basis)


def simulate(n_per=(500, 300), genes=40, seed=0, dispersion=20.0):
    rng = np.random.default_rng(seed)
    z = np.concatenate([rng.uniform(0, 1, n) for n in n_per])
    specimen = np.concatenate([[f's{i}'] * n for i, n in enumerate(n_per)])
    library = np.exp(rng.normal(np.log(4000), .5, len(z)))
    slopes = rng.choice([-1, 1], genes) * rng.uniform(1.5, 3.0, genes)
    level = np.log(rng.uniform(2e-4, 2e-3, genes))
    log_rate = level[None, :] + slopes[None, :] * (z[:, None] - .5)
    mu = library[:, None] * np.exp(log_rate)
    counts = rng.negative_binomial(dispersion, dispersion / (dispersion + mu))
    grid = grid_points(50)
    shape = piecewise_shape([1 / 6, 1 / 2, 5 / 6], slopes[:, None] * (np.array([1 / 6, 1 / 2, 5 / 6]) - .5), grid)
    return z, specimen, library, counts, shape, slopes


def test_helpers():
    grid = grid_points(11)
    shape = piecewise_shape([.2, .8], [[0., 1.]], grid)
    assert np.allclose(shape.mean(), 0) and shape[0, 0] == shape[0, 1] and shape[0, -1] == shape[0, -2]
    assert spline_basis(grid_points(50), 6).shape == (50, 6)
    assert landmark_ratio_score([0, 9.5], [0, 0])[1] == pytest.approx(np.log(20))
    with pytest.raises(ValueError):
        landmark_ratio_score([-1], [0])
    with pytest.raises(ValueError):
        piecewise_shape([.5, .2], [[0, 1]], grid)


def test_recovers_order_orientation_and_ignores_library():
    z, specimen, library, counts, shape, _ = simulate()
    fit = fit_lcp(counts, library, specimen, shape, temperature=1.0)
    rho = spearmanr(fit['mean'], z).statistic
    assert rho > .9  # oriented by the external shape, not flipped
    error = fit['mean'] - z
    assert abs(spearmanr(error, library).statistic) < .15
    assert np.all(fit['q10'] <= fit['q90']) and np.all(fit['entropy'] >= 0)
    assert fit['objective'].size == fit['iterations']


def test_projection_reproduces_fit_and_temperature_widens():
    z, specimen, library, counts, shape, _ = simulate(seed=1)
    fit = fit_lcp(counts, library, specimen, shape, temperature=2.0, max_iter=60)
    proj = project_lcp(fit, counts, library)
    assert np.allclose(proj['mean'], fit['mean'], atol=5e-3)
    hot = project_lcp(fit, counts, library, temperature=8.0)
    assert np.mean(hot['entropy']) > np.mean(proj['entropy'])


def test_balanced_aggregates_weight_specimens_equally():
    rng = np.random.default_rng(3)
    Y = np.vstack([np.full((90, 1), 10.0), np.full((10, 1), 30.0)])
    library = np.full(100, 1000.0)
    specimen = np.r_[['a'] * 90, ['b'] * 10]
    posterior = np.ones((100, 1))
    S, E = _balanced_aggregates(Y, library, posterior, specimen, .5)
    assert S[0, 0] / E[0] == pytest.approx((0.01 + 0.03) / 2)  # not the pooled 0.012
    assert rng is not None


def test_predictive_loglik_prefers_correct_model():
    z, specimen, library, counts, shape, _ = simulate(seed=4)
    rng = np.random.default_rng(4)
    a = rng.binomial(counts, .5)
    b = counts - a
    fit = fit_lcp(a, library / 2, specimen, shape, max_iter=40)
    good = predictive_loglik(fit, b, library / 2)
    shuffled = dict(fit, posterior=fit['posterior'][rng.permutation(len(z))])
    assert np.isfinite(good) and good > predictive_loglik(shuffled, b, library / 2)


def test_p4c_and_deviance_gain_detect_positional_signal_only():
    z, specimen, library, counts, shape, _ = simulate(n_per=(600, 600), seed=5)
    tr, te = specimen == 's0', specimen == 's1'
    seg = np.where(z < 1 / 3, 'S1', np.where(z < 2 / 3, 'S2', 'S3'))
    y = np.log1p(counts / library[:, None] * 1e4)
    base = {s: y[tr & (seg == s)].mean(axis=0) for s in ('S1', 'S2', 'S3')}
    btr, bte = np.vstack([base[s] for s in seg[tr]]), np.vstack([base[s] for s in seg[te]])
    noise = np.random.default_rng(99).uniform(size=len(z))
    out = residual_spline_gain(y[tr], z[tr], seg[tr], btr, y[te], z[te], seg[te], bte)
    gain = scaled_gain(out['base_mse'], out['coordinate_mse'], out['train_var'])
    null = residual_spline_gain(y[tr], noise[tr], seg[tr], btr, y[te], noise[te], seg[te], bte)
    gain_null = scaled_gain(null['base_mse'], null['coordinate_mse'], null['train_var'])
    assert gain > .02 and gain_null < .005  # a useless coordinate only overfits
    dev = deviance_gain(counts[tr], library[tr], z[tr], seg[tr], counts[te], library[te], z[te], seg[te])
    dev_null = deviance_gain(counts[tr], library[tr], noise[tr], seg[tr], counts[te], library[te], noise[te], seg[te])
    g = 1 - dev['coordinate_deviance'].sum() / dev['segment_deviance'].sum()
    g0 = 1 - dev_null['coordinate_deviance'].sum() / dev_null['segment_deviance'].sum()
    assert g > .02 and g0 < .005


def test_input_validation():
    z, specimen, library, counts, shape, _ = simulate(n_per=(50, 50), genes=5)
    with pytest.raises(ValueError):
        fit_lcp(counts, library, specimen, shape[:, :10])
    with pytest.raises(ValueError):
        fit_lcp(counts, library, specimen, shape, temperature=.5)
    with pytest.raises(ValueError):
        fit_lcp(-counts, library, specimen, shape)
