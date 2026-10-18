import numpy as np
import pytest

from pseudospace import tissue_location as tl

pytest.importorskip('patsy')


def test_program_score_and_cross_fitted_expectation():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(200, 3))
    ref = np.r_[np.ones(100, bool), np.zeros(100, bool)]
    X[~ref] += 2.
    s = tl.program_score(X, ref, floor=.01)
    assert abs(s[ref].mean()) < .1 and s[~ref].mean() > 1.5
    x = np.tile(np.linspace(0, 1, 100), 2)
    score = 2 * x + np.r_[np.zeros(100), np.ones(100)]          # case = reference + 1 everywhere
    spec = np.r_[np.repeat(['c1', 'c2'], 50), np.repeat(['k'], 100)]
    m = tl.position_expectation(x, score, spec, ['c1', 'c2'], width=.1)
    assert np.allclose((score - m)[100:], 1., atol=.15)
    # cross-fitting: c1's expectation comes from c2 only
    x2 = np.r_[np.linspace(0, 1, 50), np.linspace(0, 1, 50)]
    y2 = np.r_[np.zeros(50), np.ones(50)]
    m2 = tl.position_expectation(x2, y2, np.repeat(['c1', 'c2'], 50), ['c1', 'c2'], width=.2)
    assert np.allclose(m2[:50], 1.) and np.allclose(m2[50:], 0.)


def test_depth_effect_finds_depth_at_fixed_position_and_null_is_flat():
    rng = np.random.default_rng(1)
    n = 1500
    x = rng.uniform(0, 1, n)
    depth = np.clip(.6 * x + rng.uniform(0, .4, n), 0, 1)      # correlated with position, varies at fixed x
    strata = tl.equal_count_bins(x, 20)
    e = np.sin(3 * x) + .8 * depth + rng.normal(0, .3, n)
    out = tl.depth_effect(e, x, depth, strata, n_perm=300, seed=2)
    assert out['beta'] == pytest.approx(.8, abs=.12) and out['p'] < .01 and out['partial R2'] > .05
    flat = tl.depth_effect(np.sin(3 * x) + rng.normal(0, .3, n), x, depth, strata, n_perm=300, seed=3)
    assert flat['p'] > .01 and abs(flat['beta']) < .15
    assert np.bincount(strata).min() == np.bincount(strata).max() == 75


def test_morans_i_and_null_and_depth_and_mixture():
    rng = np.random.default_rng(4)
    xy = rng.uniform(0, 100, (800, 2))
    W = tl.knn_weights(xy, k=6)
    assert np.allclose(np.asarray(W.sum(axis=1)).ravel(), 1.) and W.diagonal().sum() == 0
    smooth = np.sin(xy[:, 0] / 10) + rng.normal(0, .2, 800)
    noise = rng.normal(size=800)
    assert tl.morans_i(smooth, W) > .5 and abs(tl.morans_i(noise, W)) < .1
    null = tl.moran_null(smooth, W, np.zeros(800, int), n_perm=200, seed=5)
    assert null.max() < .2
    d, pct = tl.glomerular_depth([[0, 0], [10, 0], [50, 0]], [[0, 1], [0, -1], [1, 0], [100, 100]], k=3)
    assert d[0] < d[1] < d[2] and pct.tolist() == pytest.approx([1 / 6, .5, 5 / 6])
    two = np.r_[rng.normal(0, .3, 300), rng.normal(3, .3, 100)]
    b = tl.mixture_bic(two)
    assert b['bic 2'] < b['bic 1'] - 10 and b['weights'][1] == pytest.approx(.25, abs=.05)
