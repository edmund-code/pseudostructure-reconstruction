"""Synthetic checks for the frozen conditional count-rate atlas."""
import numpy as np
import pytest
pytest.importorskip('scipy')
pytest.importorskip('sklearn')
from scipy.special import logsumexp
from scipy.stats import nbinom, spearmanr

from pseudospace.count_rate_atlas import (
    _nb_log_density, fit_count_rate_atlas, project_count_rate_atlas,
)
from pseudospace.canonical_mean_atlas import _coordinates


def _draw_nb(rng, library, rates, theta=100.):
    mu = library[:, None] * rates
    p = theta / (theta + mu)
    return rng.negative_binomial(theta, p)


def test_smooth_count_rate_atlas_projects_heldout_spatial_order():
    rng = np.random.default_rng(70)
    n = 1400
    z = rng.uniform(0, 1, n)
    specimens = np.repeat(['A', 'B'], n//2)
    library = rng.integers(4500, 6000, size=n).astype(float)
    rates = np.column_stack([.04 + .035*z, .08 - .035*z, np.full(n, .12)])
    counts = _draw_nb(rng, library, rates)
    atlas = fit_count_rate_atlas(counts, library, z, specimens,
                                 bandwidth=.10, min_effective=12.)
    zq = np.linspace(.05, .95, 160)
    library_q = np.full(len(zq), 5000.)
    rate_q = np.column_stack([.04 + .035*zq, .08 - .035*zq, np.full(len(zq), .12)])
    query = _draw_nb(rng, library_q, rate_q)
    projected = project_count_rate_atlas(atlas, query, library_q, block_size=17)
    assert spearmanr(projected['z_MAP'], zq).statistic > .8
    np.testing.assert_allclose(projected['posterior'].sum(axis=1), 1.)
    assert np.isfinite(projected['log_density']).all()


def test_per_position_rates_are_equal_specimen_weighted_under_replication():
    z0 = np.linspace(0, 1, 100)
    z = np.r_[z0, np.tile(z0, 2)]
    specimens = np.r_[np.repeat('A', 100), np.repeat('B', 200)]
    library = np.full(len(z), 1000.)
    rates = np.column_stack([.10 + .05*z, .20 - .05*z])
    counts = np.rint(library[:, None] * rates).astype(int)
    grid = np.linspace(0, 1, 11)
    atlas = fit_count_rate_atlas(counts, library, z, specimens, grid=grid,
                                 bandwidth=.2, min_effective=4.)
    duplicate = np.flatnonzero(specimens == 'B')
    atlas_duplicated = fit_count_rate_atlas(
        np.vstack([counts, counts[duplicate]]), np.r_[library, library[duplicate]],
        np.r_[z, z[duplicate]], np.r_[specimens, specimens[duplicate]],
        grid=grid, bandwidth=.2, min_effective=4.)
    np.testing.assert_allclose(atlas['rates'], atlas_duplicated['rates'], atol=1e-12)


def test_exposure_scaling_preserves_poisson_map_and_more_counts_sharpen_posterior():
    rng = np.random.default_rng(8)
    z = np.tile(np.linspace(0, 1, 60), 2)
    specimens = np.repeat(['A', 'B'], 60)
    library = np.full(len(z), 10000.)
    rates = np.column_stack([.03 + .04*z, .09 - .04*z, np.full(len(z), .08)])
    counts = np.rint(library[:, None] * rates).astype(int)
    atlas = fit_count_rate_atlas(counts, library, z, specimens, theta=np.inf,
                                 grid=np.linspace(0, 1, 51), bandwidth=.15,
                                 min_effective=8.)
    zq = np.linspace(.15, .85, 35)
    L = np.full(len(zq), 5000.)
    rq = np.column_stack([.03 + .04*zq, .09 - .04*zq, np.full(len(zq), .08)])
    low = np.rint(L[:, None] * rq).astype(int)
    small = project_count_rate_atlas(atlas, low, L)
    large = project_count_rate_atlas(atlas, low*10, L*10)
    np.testing.assert_array_equal(small['z_MAP'], large['z_MAP'])
    assert np.mean(small['entropy']) > np.mean(large['entropy'])


def test_nb_log_density_matches_scipy_full_profile_likelihood():
    rng = np.random.default_rng(91)
    z = np.linspace(0, 1, 80)
    library = np.full(len(z), 1000.)
    rates = np.column_stack([.1 + .03*z, .2 - .03*z])
    counts = _draw_nb(rng, library, rates)
    atlas = fit_count_rate_atlas(counts, library, z, np.repeat('s', len(z)),
                                 grid=[0., .5, 1.], bandwidth=.3,
                                 min_effective=5., theta=40.)
    query = np.array([[12, 21]])
    exposure = np.array([300.])
    result = project_count_rate_atlas(atlas, query, exposure)
    theta = atlas['theta']
    point_log_likelihood = []
    for r in atlas['rates']:
        mu = exposure[0] * r
        p = theta / (theta + mu)
        point_log_likelihood.append(nbinom.logpmf(query[0], theta, p).sum())
    expected = logsumexp(point_log_likelihood) - np.log(len(point_log_likelihood))
    np.testing.assert_allclose(result['log_density'][0], expected, rtol=1e-11, atol=1e-11)


def test_atlas_refuses_sparse_support_and_invalid_exposure():
    rng = np.random.default_rng(15)
    z = rng.uniform(0, 1, 60)
    counts = rng.poisson(20, size=(60, 2))
    library = np.full(60, 1000.)
    with pytest.raises(ValueError, match='effective support'):
        fit_count_rate_atlas(counts, library, z, np.repeat('s', 60),
                             grid=np.linspace(0, 1, 11), bandwidth=.01,
                             min_effective=12.)
    with pytest.raises(ValueError, match='positive finite exposure'):
        fit_count_rate_atlas(counts, np.zeros(60), z, np.repeat('s', 60))


def test_high_magnitude_endpoint_likelihood_cannot_push_mean_past_grid():
    grid = np.linspace(0, 1, 101)
    count, exposure = 1_000_000, 10_000_000.
    endpoint_rate = count / exposure
    target_gap = 32.4
    adjacent_rate = endpoint_rate * (1 + np.sqrt(2 * target_gap / count))
    rates = np.full((len(grid), 1), .5)
    rates[-2, 0] = adjacent_rate
    rates[-1, 0] = endpoint_rate
    atlas = {'grid': grid, 'rates': rates, 'theta': np.inf,
             'feature_count': 1}
    query = np.array([[count]], dtype=np.int64)
    library = np.array([exposure])

    # Preserve a deterministic reproduction of the former normalization path:
    # the terminal mass rounds to one and the adjacent endpoint keeps ~1e-14.
    old_logp = _nb_log_density(query, library, rates, np.inf)[0]
    old_p = np.exp(old_logp - logsumexp(old_logp))
    old_mean = old_p @ grid
    assert old_p[-2] == pytest.approx(1e-14, rel=.1)
    assert old_p.sum() > 1.
    assert old_mean > 1.

    result = project_count_rate_atlas(atlas, query, library)
    assert np.isfinite(result['posterior']).all()
    np.testing.assert_allclose(result['posterior'].sum(axis=1), 1., atol=1e-14)
    assert np.all((result['z_mean'] >= 0.) & (result['z_mean'] <= 1.))
    assert np.isfinite(result['entropy']).all() and np.all(result['entropy'] >= 0.)
    _coordinates(result['z_mean'], len(query))
