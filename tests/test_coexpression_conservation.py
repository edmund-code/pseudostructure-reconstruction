import numpy as np
import pytest

from pseudospace import coexpression_conservation as cc


def _modular(rng, n_samples, n_modules=6, size=10, noise=.5):
    factors = rng.normal(size=(n_samples, n_modules))
    X = np.repeat(factors, size, axis=1) + noise * rng.normal(size=(n_samples, n_modules * size))
    return X


def test_spearman_network_and_rank_standardise():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(50, 4))
    X[:, 1] = X[:, 0] ** 3                                     # monotone: Spearman 1
    X[:, 3] = 2.                                                # constant
    C = cc.spearman_network(X)
    assert C[0, 1] == pytest.approx(1.) and C[3, 0] == 0 and np.allclose(np.diag(C), 1)
    R = cc.rank_standardise(C)
    assert R.max() == 1 and np.allclose(R, R.T) and R[0, 1] == 1
    assert np.allclose(cc.aggregate([C, C]), R)


def test_partner_auroc_identical_random_and_conservation():
    rng = np.random.default_rng(1)
    A = cc.spearman_network(_modular(rng, 300))
    B = cc.spearman_network(_modular(np.random.default_rng(2), 300))    # same module structure, new samples
    same = cc.conservation(A, A, k=5)
    assert np.allclose(same, 1.)
    rep = cc.conservation(A, B, k=5)
    assert np.median(rep) > .95
    perm = np.random.default_rng(3).permutation(A.shape[0])
    shuffled = cc.conservation(A, B[np.ix_(perm, perm)], k=5)  # gene labels scrambled in B
    assert abs(np.median(shuffled) - .5) < .1
    with pytest.raises(ValueError):
        cc.partner_auroc(A[:5, :5], B[:5, :5], k=10)


def test_binned_residuals_remove_position():
    x = np.linspace(0, 1, 400)
    X = np.c_[np.sin(3 * x), x]
    R = cc.binned_residuals(X, x, n_bins=20)
    assert np.abs(R).max() < .1
    assert np.allclose(R.reshape(20, 20, 2).mean(axis=1), 0)
