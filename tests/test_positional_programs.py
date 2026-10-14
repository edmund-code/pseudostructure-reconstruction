import numpy as np
import pandas as pd
from scipy import sparse

from pseudospace import positional_programs as pp


def test_fixed_bins_and_profiles():
    rng = np.random.default_rng(0)
    units = np.concatenate([rng.uniform(0, 2.5, 400), [-0.1, 2.7, np.nan]])
    idx, centres = pp.fixed_bins(units)
    assert len(centres) == 25 and np.allclose(centres[:2], [.05, .15])
    assert (idx[-3:] == -1).all() and (idx[:400] >= 0).all()
    spec = np.array(['a'] * 200 + ['b'] * 203)
    counts = sparse.csr_matrix(rng.poisson(5, (403, 3)))
    prof = pp.bin_profiles(counts, np.full(403, 100.), idx, 25, spec, ['a', 'b'], min_units=1)
    assert prof['a'].shape == (3, 25)
    assert np.allclose(np.nanmean(prof['a'], axis=1), 0)


def test_boundary_excess_bump_step_ramp():
    c = np.arange(.05, 2.5, .1)
    bump = pp.bump(c, 1., .4)
    step = (c >= 1).astype(float)
    ramp = pp.ramp(c, 1., 1.)
    E, T = pp.boundary_excess(np.vstack([bump, step, ramp, -bump]), c, 1)
    assert E[0] > .5 and E[1] <= 1e-12 and E[2] <= 1e-12
    assert T[3] > .5 and T[1] <= 1e-12
    E2, _ = pp.boundary_excess(pp.bump(c, 2., .4)[None], c, 2)
    assert E2[0] > .5


def test_program_profiles_and_reference():
    D = np.array([[1., 2.], [3., 4.], [5., np.nan]])
    out = pp.program_profiles(D, [np.array([0, 1, 2]), np.array([0])], min_members=2)
    assert np.allclose(out[0], [3., 3.]) and np.isnan(out[1]).all()
    assert np.allclose(pp.reference_profile(D, [True, True, False]), [2., 3.])


def test_strata_permutations_stay_in_strata():
    strata = np.array([0, 0, 1, 1, 1, 2])
    perms = pp.strata_permutations(strata, 20, np.random.default_rng(1))
    assert perms.shape == (20, 6)
    for p in perms:
        assert sorted(p) == list(range(6)) and (strata[p] == strata).all()


def test_redundancy_and_clusters():
    groups = pp.redundancy_groups({'a': 'xyz', 'b': 'xyw', 'c': 'pq', 'd': 'xy'}, threshold=.5)
    assert groups['a'] == groups['b'] == groups['d'] and groups['c'] != groups['a']
    c = np.arange(.05, 2.5, .1)
    rng = np.random.default_rng(2)
    rows = [pp.bump(c, m, .5) + rng.normal(0, .02, len(c)) for m in (0.5, 1.5, 2.2) for _ in range(5)]
    labels, k, scores = pp.cluster_profiles(np.vstack(rows), k_range=(2, 5))
    assert k == 3 and len(set(labels[:5])) == 1 and len(set(labels)) == 3


def test_matched_step_keeps_label_means():
    rng = np.random.default_rng(3)
    v = rng.normal(size=60)
    label = np.repeat(['S1', 'S2', 'S3'], 20)
    spec = np.tile(['a', 'b'], 30)
    s = pp.matched_step(v, label, spec)
    f = pd.DataFrame({'v': v, 's': s, 'l': label, 'p': spec}).groupby(['p', 'l'])[['v', 's']].mean()
    assert np.allclose(f.v, f.s)


def test_simulate_counts_mean_and_effect():
    rng = np.random.default_rng(4)
    lib = np.full(4000, 1000.)
    u = np.linspace(0, 2.5, 4000)
    eff = np.log(2) * 2 * pp.bump(u, 1., .5)[None]
    x = pp.simulate_counts(lib, np.array([.01]), np.zeros(4000), eff, np.array([.05]), rng)
    assert x.shape == (4000, 1)
    assert abs(x.sum() / (lib.sum() * .01) - 1) < .05
    assert x[np.abs(u - 1) < .1].mean() > 3 * x[u > 2].mean()


def test_segment_call_and_r2():
    assert pp.segment_call([[.3, .3, -.6], [.4, .26, -.7]]) == 'S1+S2'
    assert pp.segment_call([[.3, .1, -.6], [.1, .3, -.7]]) == 'none'
    assert pp.pair_r2([1, 2, 3, 4], [2, 4, 6, 8]) == 1.
    assert np.isnan(pp.pair_r2([1, 2], [1, 2]))
    rng = np.random.default_rng(5)
    assert pp.thin_to(np.array([1000]), 1000., 1000., rng)[0] == 1000
