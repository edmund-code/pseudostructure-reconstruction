import numpy as np
import pandas as pd
from scipy import sparse

from pseudospace import continuous_axis as ca

C = np.arange(.05, 2.5, .1)


def _label_frac(width=.15):
    """Label shares of a mouse whose label switches are logistic of 10–90% width ``width`` at 1 and 2."""
    s12 = 1 / (1 + np.exp(-(C - 1) * np.log(81) / width))
    s23 = 1 / (1 + np.exp(-(C - 2) * np.log(81) / width))
    return np.column_stack([1 - s12, s12 - s23, s23])


def _fit(profiles, noise=0., seed=0):
    rng = np.random.default_rng(seed)
    f = _label_frac()
    P = [np.asarray(profiles) + noise * rng.normal(size=np.shape(profiles)) for _ in range(2)]
    return ca.ladder_class(ca.fit_ladder(P, [np.ones(len(C))] * 2, C, [f, f]))


def test_label_fractions_rows_sum_to_one():
    bins = np.array([0, 0, 1, 1, 1, -1])
    labels = np.array(['a', 'b', 'b', 'b', 'c', 'a'])
    f = ca.label_fractions(bins, labels, 3, ['a', 'b', 'c'])
    assert np.allclose(f[0], [.5, .5, 0]) and np.allclose(f[1], [0, 2 / 3, 1 / 3]) and np.isnan(f[2]).all()


def test_ladder_recovers_each_shape():
    f = _label_frac()
    shapes = np.vstack([
        np.zeros(len(C)),                                   # flat
        .6 * (C - C.mean()),                                # gradient
        1.2 / (1 + np.exp(-(C - .6) * np.log(81) / .3)),    # switch inside S1
        np.exp(-4 * np.log(2) * (C - 1.5) ** 2 / .4 ** 2),  # bump in S2
        f[:, 1] * .8 + f[:, 2] * 1.6,                       # constant within labels
    ])
    t = _fit(shapes, noise=.02)
    assert list(t['class'][:4]) == ['flat', 'gradient', 'switch', 'bump']
    assert t.loc[4, 'class'] in ('discrete', 'discrete + within')   # DS wins ~2% of noise draws of a pure step
    assert abs(t.loc[2, 'S position'] - .6) < .06 and t.loc[2, 'location'] == 'inside S1'
    assert abs(t.loc[3, 'B position'] - 1.5) < .06 and t.loc[3, 'location'] == 'inside S2'
    assert t.loc[4, 'dBIC continuous'] < -6 and t.loc[2, 'dBIC continuous'] > 6
    assert t.loc[0, 'continuous vs discrete'] == 'flat'


def test_ladder_curve_and_half_point():
    f = _label_frac()
    row = pd.Series({'S amp': 1., 'S position': 1.2, 'S width': .5, 'D S2-S1': 1., 'D S3-S1': 1., 'L slope': 1.})
    s = ca.ladder_curve('S', row, C)
    assert abs(ca.half_point(C, s) - 1.2) < .02
    assert abs(ca.half_point(C, ca.ladder_curve('D', row, C, f)) - 1.) < .05
    bump = np.exp(-(C - 1.3) ** 2 / .05)
    assert abs(ca.half_point(C, bump) - 1.25) < .11
    assert np.isnan(ca.half_point(C, np.zeros(len(C))))


def test_centre_of_mass():
    D = np.vstack([np.zeros(len(C)), np.where(C > 2, 3., 0.)])
    com = ca.centre_of_mass(D, C)
    assert np.isclose(com[0], C.mean()) and com[1] > 1.8


def test_aligned_amplitude_finds_opposite_members():
    rng = np.random.default_rng(1)
    up = np.vstack([C - C.mean()] * 6)
    members = np.vstack([up, -up]) + .05 * rng.normal(size=(12, len(C)))
    test = np.vstack([up, -up]) + .05 * rng.normal(size=(12, len(C)))
    amp, signs = ca.aligned_amplitude(members, test)
    assert amp > 2 and len(set(signs[:6])) == 1 and signs[0] == -signs[-1]
    assert np.ptp(np.median(test, axis=0)) < .3


def test_decoder_recovers_bins():
    rng = np.random.default_rng(2)
    n_bins, n_genes = 10, 40
    rates = np.exp(rng.normal(-3, 1, (n_genes, 1)) + rng.normal(0, 1, (n_genes, n_bins)))
    true = rng.integers(0, n_bins, 500)
    lib = np.full(500, 2000.)
    counts = rng.poisson(lib[:, None] * rates[:, true].T)
    ll = ca.nb_bin_loglik(sparse.csr_matrix(counts), lib, rates, np.full(n_genes, .01))
    post = ca.log_posterior(ll, np.log(np.full(n_bins, 1 / n_bins)))
    assert np.allclose(np.exp(post).sum(axis=1), 1)
    assert (post.argmax(axis=1) == true).mean() > .8


def test_mutual_information():
    rng = np.random.default_rng(3)
    x = rng.integers(0, 4, 4000)
    assert abs(ca.plugin_mi(x, x, 4, 4) - 2) < .02
    assert abs(ca.corrected_mi(x, rng.integers(0, 4, 4000), 4, 4, rng)) < .01
    lab = (x >= 2).astype(int)
    # x is fully known from itself; within each half of x one bit remains
    assert abs(ca.within_label_mi(x, x, lab, 4, 4, rng) - 1) < .03


def test_end_contrast_and_reliability():
    D = np.vstack([C - C.mean(), np.zeros(len(C))])
    e = ca.end_contrast(D, n_end=8)
    assert e[0] > 1.5 and e[1] == 0
    assert np.isclose(ca.reliability(.5), 2 / 3)
    assert np.isclose(ca.disattenuated_r(.5, .5, .5), .75)


def test_graded_categories():
    out = ca.graded_categories([1, 1, 1, .1, 1], [1, 1, 1, -.1, 1], [1, -1, .1, 1, .3],
                               [1, -1, -.1, 1, .3])
    assert list(out) == ['conserved', 'reversed', 'first only', 'second only', '']


def test_local_change_of_a_step():
    step = np.where(C > 1, 1., 0.)[None]
    d = ca.local_change(step, half=2)
    pos = np.arange(2, len(C) - 1) * .1
    assert d.shape[1] == len(pos) and np.isclose(d[0, np.argmin(np.abs(pos - 1))], 1.)
    assert np.allclose(d[0, np.abs(pos - 1) > .25], 0)


def test_discrete_profiles_step_at_labels():
    rng = np.random.default_rng(4)
    n = 4000
    u = rng.uniform(0, 2.5, n)
    labels = np.where(u < 1, 'S1', np.where(u < 2, 'S2', 'S3'))
    bins = np.floor(u / .1).astype(int)
    rate = np.where(labels == 'S1', 1., np.where(labels == 'S2', 2., 4.))
    counts = sparse.csr_matrix(rng.poisson(rate * 50)[:, None])
    out = ca.discrete_profiles(counts, np.full(n, 1000.), labels, bins, 25, np.array(['m'] * n), ['m'],
                               ['S1', 'S2', 'S3'])
    y = out['m'][0]
    assert np.isclose(y[15] - y[5], 1, atol=.1) and np.isclose(y[22] - y[15], 1, atol=.1)


def test_matched_fold_changes_remove_composition():
    rng = np.random.default_rng(5)
    n = 6000
    cond = np.array(['case', 'ref'] * (n // 2))
    u = np.where(cond == 'case', rng.uniform(.5, 1, n), rng.uniform(0, 1, n))   # case sits late in the label
    bins = np.floor(u / .1).astype(int)
    rate = np.exp2(2 * u)                                                         # graded, no condition effect
    counts = sparse.csr_matrix(rng.poisson(rate * 200)[:, None])
    pooled, matched, shift = ca.matched_fold_changes(counts, np.full(n, 1000.), bins, np.array(['S1'] * n), cond, 10,
                                                     'S1')
    assert pooled[0] > .35 and abs(matched[0]) < .05 and shift > 2
