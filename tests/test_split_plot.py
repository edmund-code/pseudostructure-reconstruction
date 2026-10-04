import numpy as np
import pytest
from scipy import stats

from pseudospace.split_plot import log_cpm, moderated_f, squeeze_var, voom_weights


def _layout(n_bins=8):
    specimen = np.repeat(np.arange(4), n_bins)
    species = (specimen >= 2).astype(float)
    position = np.tile(np.linspace(0, 1, n_bins), 4)
    poly = np.column_stack([position, position ** 2, position ** 3])
    poly = (poly - poly.mean(0)) / poly.std(0)
    design = np.column_stack([np.eye(4)[specimen], poly, species[:, None] * poly])
    return design, specimen, position, species, list(range(7, 10))


def test_log_cpm_rejects_bad_library():
    with pytest.raises(ValueError):
        log_cpm(np.ones((3, 2)), np.array([1., 0., 2.]))


def test_squeeze_var_recovers_prior_df():
    rng = np.random.default_rng(1)
    d0, df = 8., 6.
    true = 1. / (rng.chisquare(d0, 20000) / d0)          # scaled inverse chi-square, s0^2 = 1
    s2 = true * rng.chisquare(df, true.size) / df
    posterior, d0_hat, s02_hat = squeeze_var(s2, df)
    assert 5 < d0_hat < 12 and 0.8 < s02_hat < 1.25
    assert posterior.var() < s2.var()


def test_squeeze_var_without_extra_variation_pools_fully():
    rng = np.random.default_rng(2)
    s2 = rng.chisquare(10, 5000) / 10
    posterior, d0, _ = squeeze_var(s2, 10)
    assert d0 > 50
    assert np.ptp(posterior) < np.ptp(s2) / 3


def test_voom_weights_favour_high_counts():
    rng = np.random.default_rng(3)
    design, *_ = _layout()
    means = np.exp(rng.uniform(1, 8, 400))
    counts = rng.poisson(means, size=(design.shape[0], 400))
    library = counts.sum(1) + 1e5
    weights = voom_weights(log_cpm(counts, library), library, design)
    low, high = np.argsort(means)[:50], np.argsort(means)[-50:]
    assert weights[:, high].mean() > 3 * weights[:, low].mean()


def test_moderated_f_is_calibrated_under_null_and_detects_interaction():
    rng = np.random.default_rng(4)
    design, specimen, position, species, test = _layout()
    genes = 3000
    scale = rng.uniform(.5, 2, genes)                      # gene-specific noise for moderation
    y = rng.normal(0, 1, (design.shape[0], genes)) * scale
    for k in range(4):                                     # specimen offsets are absorbed
        y[specimen == k] += rng.normal(0, 2, genes)
    result = moderated_f(y, np.ones_like(y), design, test)
    assert abs(np.mean(result['p'] <= .05) - .05) < .015
    assert stats.kstest(result['p'], 'uniform').pvalue > 1e-3
    signal = y.copy()
    signal[:, :300] += 8 * scale[:300] * species[:, None] * (position[:, None] - .5)
    detected = moderated_f(signal, np.ones_like(y), design, test)
    assert np.mean(detected['p'][:300] <= .01) > .8


def test_moderated_f_rejects_rank_deficient_design():
    design, *_ = _layout()
    bad = np.column_stack([design, design[:, 0]])
    with pytest.raises(ValueError):
        moderated_f(np.zeros((design.shape[0], 5)), np.ones((design.shape[0], 5)), bad, [1])
