import numpy as np
import pytest
from scipy import sparse
from scipy.stats import spearmanr

from pseudospace.coordinate_evaluation import (absorption_ratio, binomial_thin_to_depth, concordance,
                                               fold_noise_variance, group_rank_correlations,
                                               log_normalize, poisson_count_split, reflect_unit,
                                               scale_counts, sha256_folds)


def _counts(rng, n=200, g=12, lam=3.0):
    return sparse.csr_matrix(rng.poisson(lam, size=(n, g)))


def test_group_rank_correlations_match_scipy_and_ignore_monotone_warps():
    rng = np.random.default_rng(0)
    x = _counts(rng, n=120, g=6)
    z = rng.uniform(size=120)
    groups = np.repeat(['a', 'b'], 60)
    table = group_rank_correlations(x, z, groups)
    dense = x.toarray()
    for j, label in enumerate(['a', 'b']):
        rows = groups == label
        for g in range(6):
            assert table.iloc[g, j] == pytest.approx(spearmanr(dense[rows, g], z[rows]).statistic)
    warped = group_rank_correlations(x, np.exp(5 * z), groups)
    np.testing.assert_allclose(warped.to_numpy(), table.to_numpy())


def test_group_rank_correlations_flag_constant_genes():
    x = sparse.csr_matrix(np.column_stack([np.zeros(30), np.arange(30)]))
    table = group_rank_correlations(x, np.arange(30.0), np.repeat('a', 30))
    assert np.isnan(table.iloc[0, 0]) and table.iloc[1, 0] == pytest.approx(1.0)


def test_concordance_skips_missing_genes():
    rho, n = concordance([1, 2, np.nan, 4, 5], [2, 4, 1, 8, 10])
    assert n == 4 and rho == pytest.approx(1.0)


def test_absorption_ratio_limits():
    original = np.zeros(5)
    frozen = np.array([.5, -.4, .3, .01, .2])
    full = absorption_ratio(original, frozen, frozen)
    assert full['median_ratio'] == pytest.approx(1.0) and full['slope'] == pytest.approx(1.0)
    assert full['n_ratio_genes'] == 4  # the 0.01 change is below the floor
    absorbed = absorption_ratio(original, frozen, original)
    assert absorbed['median_ratio'] == pytest.approx(0.0) and absorbed['slope'] == pytest.approx(0.0)


def test_thinning_hits_target_in_expectation_and_spares_shallow_rows():
    rng = np.random.default_rng(1)
    x = _counts(rng, n=50, g=40, lam=5.0)
    library = np.asarray(x.sum(axis=1)).ravel().astype(float)
    target = float(np.quantile(library, .3))
    draws = np.array([np.asarray(binomial_thin_to_depth(x, library, target, rng).sum(axis=1)).ravel()
                      for _ in range(300)])
    deep = library > target
    np.testing.assert_allclose(draws.mean(axis=0)[deep], target, rtol=.03)
    np.testing.assert_array_equal(draws[:, ~deep], np.broadcast_to(library[~deep], draws[:, ~deep].shape))


def test_count_split_adds_up_and_halves_in_expectation():
    rng = np.random.default_rng(2)
    x = _counts(rng, lam=8.0)
    a, b = poisson_count_split(x, .5, rng)
    np.testing.assert_array_equal((a + b).toarray(), x.toarray())
    assert (a.toarray() >= 0).all() and (b.toarray() >= 0).all()
    assert a.sum() / x.sum() == pytest.approx(.5, abs=.02)


def test_scale_counts_has_the_requested_expectation_and_leaves_other_entries():
    rng = np.random.default_rng(3)
    x = sparse.csr_matrix(np.full((40, 5), 20))
    rows, cols = np.arange(10), np.array([1, 3])
    log_factor = np.column_stack([np.full(10, np.log(2)), np.full(10, np.log(.5))])
    draws = np.stack([scale_counts(x, rows, cols, log_factor, rng).toarray() for _ in range(400)])
    assert draws[:, :10, 1].mean() == pytest.approx(40, rel=.03)
    assert draws[:, :10, 3].mean() == pytest.approx(10, rel=.03)
    np.testing.assert_array_equal(draws[:, 10:, :], 20)
    np.testing.assert_array_equal(draws[:, :, [0, 2, 4]], 20)


def test_fold_noise_variance_recovers_known_panel_noise():
    rng = np.random.default_rng(4)
    n = 4000
    truth = rng.uniform(size=n)
    species = np.repeat(['mouse', 'human'], n // 2)
    sd = np.where(species == 'human', .08, .03)
    coords = {k: truth + rng.normal(0, sd) for k in ('f0', 'f1', 'f2')}
    sigma2, _ = fold_noise_variance(coords, np.repeat(['m1', 'h1'], n // 2), np.tile(['S1', 'S2'], n // 2), species)
    assert sigma2['human'] == pytest.approx(.08 ** 2, rel=.1)
    assert sigma2['mouse'] == pytest.approx(.03 ** 2, rel=.1)


def test_reflect_unit_and_helpers():
    np.testing.assert_allclose(reflect_unit([-.2, .3, 1.25, 2.1]), [.2, .3, .75, .1])
    assert len(set(sha256_folds(['a', 'b', 'c'] * 10))) <= 3
    np.testing.assert_array_equal(sha256_folds(['Gene1', 'Gene2']), sha256_folds(['Gene1', 'Gene2']))
    out = log_normalize(sparse.csr_matrix([[1, 3], [0, 4]]), [4, 4], target_sum=4)
    np.testing.assert_allclose(out.toarray(), np.log1p([[1, 3], [0, 4]]))
