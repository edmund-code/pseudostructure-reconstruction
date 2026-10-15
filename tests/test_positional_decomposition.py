import numpy as np
from scipy import sparse

from pseudospace import positional_decomposition as pd_
from pseudospace.positional_programs import strata_permutations

C = np.arange(.05, 2.5, .1)


def test_level_location_splits_a_shift_from_a_gradient():
    D = np.vstack([np.full(len(C), .8), C - C.mean(), .8 + (C - C.mean())])
    D[:, :2] = np.nan
    L, P, M = pd_.level_location(D)
    assert np.allclose(L, [.8, .1, .9], atol=1e-9)   # bins 0.05 and 0.15 are invalid
    assert M[0] < 1e-12 and np.isclose(M[1], M[2]) and np.isnan(P[0, 0])


def test_row_corr():
    a = np.vstack([C, C, np.ones(len(C))])
    b = np.vstack([2 * C + 1, -C, C])
    r = pd_.row_corr(a, b)
    assert np.isclose(r[0], 1) and np.isclose(r[1], -1) and np.isnan(r[2])


def test_cross_fitted_survival_is_unbiased_where_naive_flattens():
    rng = np.random.default_rng(0)
    g = np.vstack([np.sin(C * k) for k in np.linspace(1, 3, 400)])
    noise = 1.
    refs = [g + noise * rng.normal(size=g.shape) for _ in range(2)]
    cases = [.5 * g + noise * rng.normal(size=g.shape) for _ in range(2)]
    a, num, den = pd_.cross_fitted_survival(cases, refs)
    assert abs(pd_.pooled_ratio(num, den) - .5) < .05
    # naive: template from the mean of both refs, denominator from the same refs -> biased low
    t = (refs[0] + refs[1]) / 2
    naive = (np.mean(cases, axis=0) * t).sum() / (t * t).sum()
    assert naive < .45
    with np.testing.assert_raises(ValueError):
        pd_.cross_fitted_survival(cases, refs[:1])


def test_survival_class_edges():
    out = pd_.survival_class([-1, -.25, 0, .25, .5, .6, 1, 1.5, 2, np.nan])
    assert list(out) == ['reversed', 'erased', 'erased', 'erased', 'attenuated', 'attenuated', 'kept', 'kept',
                         'amplified', '']


def test_set_statistics_find_a_shifted_set():
    rng = np.random.default_rng(1)
    n = 600
    strata = np.repeat(np.arange(6), 100)
    perms = strata_permutations(strata, 500, rng)
    values = rng.normal(size=n)
    members = np.arange(0, 600, 30)
    values[members] += 2
    obs, z, p = pd_.set_mean_z(values, members, perms)
    assert obs > 1.5 and z > 4 and p < 1e-4
    prof = rng.normal(size=(n, len(C)))
    prof[members] += np.sin(C * 2)
    full, size, z2, p2 = pd_.set_profile_z(prof, members, perms)
    assert z2 > 4 and np.corrcoef(full, np.sin(C * 2))[0, 1] > .8
    num = np.abs(rng.normal(1, .1, (n, 4)))
    den = num.copy()
    num[members] *= -.5
    r, z3, _ = pd_.set_ratio_z(num, den, members, perms)
    assert np.isclose(r, -.5) and z3 < -4


def test_cosine_and_rms_ratio():
    a = np.array([1., -1, 2, 0])
    assert np.isclose(pd_.gradient_cosine(a, -a), -1) and np.isclose(pd_.gradient_cosine(a, 3 * a), 1)
    rng = np.random.default_rng(2)
    n = 500
    perms = strata_permutations(np.zeros(n, int), 300, rng)
    x, y = rng.normal(size=n), rng.normal(size=n)
    members = np.arange(20)
    y[members] = -x[members]
    obs, z, _ = pd_.set_cosine_z(x, y, members, perms)
    assert obs < -.99 and z < -4
    y2 = rng.normal(size=n)
    y2[members] = 4 * x[members]
    q, zq, _ = pd_.set_log_rms_ratio_z(x, y2, members, perms)
    assert abs(q - 2) < 1e-9 and zq > 4


def test_interior_dip():
    dip = 1 - .5 * np.exp(-(C - .8) ** 2 / .05)
    ok, pos, ex = pd_.interior_dip(dip, C)
    assert ok and abs(pos - .75) < .11 and ex > .5
    mono = C.copy()
    assert not pd_.interior_dip(mono, C)[0]
    shallow = 1 - .02 * np.exp(-(C - .8) ** 2 / .05) + .5 * (C > 2)
    assert not pd_.interior_dip(shallow, C)[0]


def test_pseudobulk_log_cpm():
    X = sparse.csr_matrix(np.array([[1, 3], [1, 1], [0, 8]]))
    out = pd_.pseudobulk_log_cpm(X, np.array([4., 2., 8.]), np.array(['a', 'a', 'b']), ['a', 'b'], pseudo=0.)
    with np.errstate(divide='ignore'):
        assert np.allclose(out, np.log2(np.array([[2 / 6, 4 / 6], [0, 1]]) * 1e6))
