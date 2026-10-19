import numpy as np
import pytest
from scipy import sparse

from pseudospace import heterogeneity as ht


def test_moment_dispersion_recovers_nb_and_poisson():
    rng = np.random.default_rng(0)
    n_units = 4000
    lib = rng.uniform(500, 2000, n_units)
    rate = np.array([.002, .01, .01])
    alpha = np.array([.5, .05, 0.])
    cols = []
    for r, a in zip(rate, alpha):
        mu = lib * r
        lam = rng.gamma(1 / a, mu * a) if a > 0 else mu
        cols.append(rng.poisson(lam))
    y = sparse.csr_matrix(np.column_stack(cols))
    est, mean = ht.moment_dispersion(y, lib)
    assert est[0] == pytest.approx(.5, rel=.15) and est[1] == pytest.approx(.05, rel=.25)
    assert est[2] <= .01                                    # Poisson gene sits at the floor or near 0
    assert mean == pytest.approx(np.asarray(y.mean(axis=0)).ravel())


def test_residual_and_detection():
    mean = np.exp(np.linspace(-1, 3, 50))
    alpha = np.exp(-np.log(mean) + np.random.default_rng(1).normal(0, .05, 50))   # trend with a little noise
    alpha[10] *= np.e                                       # one gene one log unit above the trend
    eps, trend = ht.residual_dispersion(alpha, mean)
    assert eps[10] == pytest.approx(1, abs=.25) and np.nanmax(np.abs(np.delete(eps, 10))) < .3
    y = sparse.csr_matrix(np.array([[0, 1], [2, 0], [0, 0], [1, 1]]))
    assert ht.detection(y).tolist() == [.5, .5]
