import numpy as np
import pytest
from scipy import sparse

from pseudospace import glucose_handling as gh
from pseudospace.positional_programs import bin_profiles, reference_profile


@pytest.mark.filterwarnings('ignore::RuntimeWarning')   # the emptied bin has no library
def test_bootstrap_deviation_centres_like_the_unresampled_profile():
    rng = np.random.default_rng(0)
    n_units, n_bins = 600, 5
    bins = np.repeat(np.arange(n_bins), n_units // n_bins)
    rate = np.vstack([np.full(n_bins, 2.), np.linspace(1, 4, n_bins), np.full(n_bins, 3.)]).T   # bins x genes
    lib = np.full(n_units, 100.)
    counts = sparse.csr_matrix(rng.poisson(rate[bins] * 5))
    usable = np.array([True, False, True])
    boot = gh.bootstrap_deviation(counts, lib, bins, n_bins, [1], usable, n_boot=40, seed=1)
    assert boot.shape == (40, 1, n_bins)
    prof = bin_profiles(counts, lib, bins, n_bins, np.zeros(n_units, int), [0])[0]
    plain = prof[1] - reference_profile(prof, usable)
    assert np.allclose(boot.mean(axis=0)[0], plain, atol=.05)
    assert np.all(np.diff(boot.mean(axis=0)[0]) > 0)
    few = gh.bootstrap_deviation(counts, lib, np.where(bins == 4, -1, bins), n_bins, [1], usable, n_boot=2, seed=1)
    assert np.isnan(few[:, :, 4]).all()


def test_balanced_relabelings():
    assert gh.balanced_relabelings(2) == (2, 1 / 3)
    others, p = gh.balanced_relabelings(3)
    assert others == 9 and np.isclose(p, .1)
    assert gh.balanced_relabelings(4)[0] == 34
    with pytest.raises(ValueError):
        gh.balanced_relabelings(0)


def test_within_group_sd_recovers_the_specimen_sd():
    rng = np.random.default_rng(2)
    n_genes = 20000
    strata = np.repeat([0, 1], n_genes // 2)
    sd_true = np.where(strata == 0, .2, .5)
    values = rng.normal(size=(n_genes, 4)) * sd_true[:, None] + rng.normal(size=(n_genes, 1))
    values[:50, 0] += 5                        # a few genes with a real specimen difference
    var, df, sd = gh.within_group_sd(values, [[0, 1], [2, 3]], strata)
    assert df == 2 and var.shape == (n_genes,)
    assert abs(sd[0] - .2) < .01 and abs(sd[1] - .5) < .02
    with pytest.raises(ValueError):
        gh.within_group_sd(values, [[0], [1]], strata)


def test_minimum_detectable_effect_matches_simulated_power():
    sd = .3
    for n in (3, 5):
        mde = gh.minimum_detectable_effect(sd, n)
        rng = np.random.default_rng(n)
        a = rng.normal(0, sd, (4000, n))
        b = rng.normal(mde, sd, (4000, n))
        from scipy import stats
        p = stats.ttest_ind(a, b, axis=1).pvalue
        assert abs((p < .05).mean() - .8) < .03     # the noncentral t makes the approximation slightly generous
    assert gh.minimum_detectable_effect(sd, 2) > gh.minimum_detectable_effect(sd, 4) > gh.minimum_detectable_effect(sd, 8)
    with pytest.raises(ValueError):
        gh.minimum_detectable_effect(sd, 1)


def test_luminal_glucose_balances_and_calibrates():
    n = 25
    s2 = np.r_[np.ones(15), np.zeros(10)]
    s1 = np.r_[np.zeros(15), np.ones(10)]
    F, up2, up1 = gh.luminal_glucose(s2, s1, 20., 1., c0=8.)
    assert F[0] == 1 and np.all(np.diff(F) <= 1e-12)
    assert np.isclose(1 - F[-1], up2.sum() + up1.sum(), atol=1e-6)
    assert up1[:15].sum() == 0 and up2[15:].sum() == 0
    a2, a1 = gh.calibrate_capacities(s2, s1, block_reabsorbed=.4, threshold=16.)
    blocked = gh.luminal_glucose(s2, s1, 0., a1, c0=8.)[0]
    assert np.isclose(1 - blocked[-1], .4, atol=1e-6)
    at_threshold = gh.luminal_glucose(s2, s1, a2, a1, c0=16.)[0]
    assert np.isclose(at_threshold[-1], .01, atol=1e-6)
    normal = gh.luminal_glucose(s2, s1, a2, a1, c0=8.)[0]
    assert normal[-1] < .01                     # below threshold almost nothing escapes
    with pytest.raises(ValueError):
        gh.luminal_glucose(s2, s1[:-1], 1., 1., c0=8.)


def test_switch_fit_matches_the_ladder_switch_model():
    from pseudospace import continuous_axis as ca
    rng = np.random.default_rng(3)
    centres = np.arange(.05, 2.5, .1)
    truth = [(.8, .5, -2.), (1.9, 1.2, 1.5), (1.3, .3, 1.)]
    Y = np.vstack([a / (1 + np.exp(-(centres - m_) * np.log(81) / w_)) for m_, w_, a in truth])
    Y += rng.normal(0, .05, Y.shape)
    weights = rng.integers(30, 200, len(centres)).astype(float)
    frac = np.column_stack([centres < 1, (centres >= 1) & (centres < 2), centres >= 2]).astype(float)
    lad = ca.fit_ladder([Y], [weights], centres, [frac])
    pos, width, amp = gh.switch_fit(Y, weights, centres)
    assert np.allclose(pos, lad['S position']) and np.allclose(width, lad['S width'])
    assert np.allclose(amp, lad['S amp'])
    assert np.allclose(pos, [m_ for m_, _, _ in truth], atol=.1)
    Y[:, 3] = np.nan
    assert np.isfinite(gh.switch_fit(Y, weights, centres)[0]).all()
