"""Synthetic checks for the phase-3 landmark reconstruction module (no private data)."""
import numpy as np
import pytest
from scipy.stats import spearmanr

from pseudospace.landmark_reconstruction import (PathBasis, combine_half_posteriors, design_effect_temperature,
                                                 deviance_gain_repaired, fit_landmark_model, lrs_bootstrap_sd,
                                                 project_landmark_model, standardized_difference,
                                                 wsr_correlations, wsr_statistic, zone_f_statistic,
                                                 zone_log_prior, zone_order)


def simulate(n_per=(400, 250), genes=30, seed=0, dispersion=20.0, warp=None):
    rng = np.random.default_rng(seed)
    z = np.concatenate([rng.uniform(0, 1, n) for n in n_per])
    if warp is not None:
        z = warp(z)
    specimen = np.concatenate([[f's{i}'] * n for i, n in enumerate(n_per)])
    library = np.exp(rng.normal(np.log(3000), .5, len(z)))
    slopes = rng.choice([-1, 1], genes) * rng.uniform(1.5, 3.0, genes)
    level = np.log(rng.uniform(5e-4, 3e-3, genes))
    mu = library[:, None] * np.exp(level[None, :] + slopes[None, :] * (z[:, None] - .5))
    counts = rng.negative_binomial(dispersion, dispersion / (dispersion + mu))
    pts = np.array([1 / 6, 1 / 2, 5 / 6])
    block = {'counts': counts, 'library': library, 'specimen': specimen, 'shape_points': pts,
             'shape_values': slopes[:, None] * (pts - .5)[None, :], 'temperature': 1.0}
    return z, block, slopes


def test_basis_and_prior():
    b = PathBasis('spline', 6)
    assert b.at(np.linspace(0, 1, 11)).shape == (11, 6)
    assert PathBasis('linear').design(np.array([0., 1.])).shape == (2, 2)
    grid = np.linspace(0, 1, 20)
    lp = zone_log_prior(np.array([0, 1, -1]), np.array([0, .5, 1.]), grid, eps=.1)
    assert np.allclose(np.exp(lp).sum(axis=1), 1)
    assert np.exp(lp[0, grid < .5]).sum() == pytest.approx(.9)
    with pytest.raises(ValueError):
        zone_log_prior(np.array([0]), np.array([0, .5]), grid)
    with pytest.raises(ValueError):
        PathBasis('spline', 3)


def test_anchored_fit_converges_and_orders():
    z, block, _ = simulate()
    (res,), diag = fit_landmark_model([block], fine=80, coarse=30, fine_max=150)
    assert spearmanr(res['mean'], z).statistic > .9
    assert diag['converged'] and diag['iterations_fine'] < 150
    proj = project_landmark_model(res, block['counts'], block['library'])
    assert np.allclose(proj['mean'], res['mean'], atol=5e-3)
    assert np.all(res['sd'] >= 0)


def test_linear_basis_and_pooled_variant_run():
    z, block, _ = simulate(seed=2)
    (res,), diag = fit_landmark_model([block], basis='linear', balanced=False, fine=60, coarse=30, fine_max=80)
    assert spearmanr(res['mean'], z).statistic > .85


def test_shared_shapes_two_blocks():
    z0, b0, s0 = simulate(seed=3)
    rng = np.random.default_rng(9)
    z1, b1, _ = simulate(seed=4)
    pairs = [(g, g) for g in range(10)]
    out, diag = fit_landmark_model([b0, b1], shared=pairs, fine=60, coarse=30, fine_max=100)
    assert spearmanr(out[0]['mean'], z0).statistic > .85
    assert len(out) == 2 and np.isfinite(diag['objective']).all()
    assert rng is not None


def test_zone_prior_keeps_labels_and_combination():
    z, block, _ = simulate(seed=5)
    labels = np.digitize(z, [1 / 3, 2 / 3])
    block = dict(block, labels=labels, zone_edges=np.array([0, 1 / 3, 2 / 3, 1]), eps=.1)
    (res,), _ = fit_landmark_model([block], fine=60, coarse=30, fine_max=100)
    assert spearmanr(res['mean'], z).statistic > .9
    post = combine_half_posteriors(res['posterior'], res['posterior'])
    assert np.allclose(post.sum(axis=1), 1) and np.mean(post.max(axis=1)) >= np.mean(res['posterior'].max(axis=1))


def test_temperature_and_standardized_difference():
    rng = np.random.default_rng(0)
    m1, m2 = rng.normal(0, 2, 1000), rng.normal(0, 2, 1000)
    T, d = design_effect_temperature(m1, np.ones(1000), m2, np.ones(1000))
    assert T == pytest.approx(4.0, rel=.25)  # true scale of the difference is 2x the model's
    assert standardized_difference([1.], [1.], [0.], [0.])[0] == pytest.approx(1.0, rel=1e-6)


def test_lrs_bootstrap_scales_with_counts():
    rng = np.random.default_rng(1)
    small = lrs_bootstrap_sd(rng.poisson(2, (200, 20)), rng.poisson(2, (200, 20)), n_boot=50)
    large = lrs_bootstrap_sd(rng.poisson(50, (200, 20)), rng.poisson(50, (200, 20)), n_boot=50)
    assert np.median(small) > np.median(large) > 0


def test_zone_order_and_f_statistic():
    profiles = np.array([[0., 0.], [1., 0.], [2., 0.1], [3., 0.]])
    order, _ = zone_order(profiles[[0, 2, 1, 3]], start=0)
    assert order == [0, 2, 1, 3]
    rng = np.random.default_rng(2)
    zone = np.repeat([0, 1, 2], 100)
    library = np.full(300, 1000.)
    rate = np.column_stack([np.repeat([1e-3, 5e-3, 1e-2], 100), np.full(300, 5e-3)])
    counts = rng.poisson(library[:, None] * rate)
    F, rates = zone_f_statistic(counts, library, zone, np.repeat(['a', 'b'], 150), 3)
    assert F[0] > 20 * max(F[1], 1) and rates.shape == (2, 3)


def test_wsr_detects_transfer_and_not_noise():
    rng = np.random.default_rng(3)
    n, G = 300, 50
    seg = np.repeat(['S1', 'S2'], n // 2)
    z_tr, z_te = rng.uniform(size=n), rng.uniform(size=n)
    w = rng.normal(size=G)
    R_tr = np.outer(np.sin(3 * z_tr), w) + rng.normal(0, 2, (n, G))
    R_te = np.outer(np.sin(3 * z_te), w) + rng.normal(0, 2, (n, G))
    perms = {s: np.array([rng.permutation(n // 2) for _ in range(100)]) for s in ('S1', 'S2')}
    out = wsr_correlations(R_tr, z_tr, seg, R_te, z_te, seg, perms)
    obs, null = out['S1']
    assert wsr_statistic(obs) > np.max(wsr_statistic(null))
    noise = wsr_correlations(R_tr, z_tr, seg, rng.normal(size=(n, G)), z_te, seg, perms)['S1']
    assert (wsr_statistic(noise[1]) >= wsr_statistic(noise[0])).mean() > .01


def test_repaired_deviance_is_finite_and_excludes_sparse():
    rng = np.random.default_rng(4)
    n = 400
    seg = np.repeat(['S1', 'S2'], n // 2)
    z = rng.uniform(size=n)
    lib = np.full(n, 2000.)
    rate = np.column_stack([1e-3 * np.exp(2 * z), np.full(n, 1e-7)])
    y = rng.poisson(lib[:, None] * rate)
    tr, te = np.arange(n) % 2 == 0, np.arange(n) % 2 == 1
    out = deviance_gain_repaired(y[tr], lib[tr], z[tr], seg[tr], y[te], lib[te], z[te], seg[te])
    assert np.isfinite(out['coordinate_deviance']).all() and out['segment_deviance'][:, 1].sum() == 0
    assert out['coordinate_deviance'][:, 0].sum() < out['segment_deviance'][:, 0].sum()
