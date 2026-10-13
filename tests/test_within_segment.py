import numpy as np
import pandas as pd
import pytest

pytest.importorskip('rpy2')
from pseudospace import within_segment as ws


def _needs_limma():
    try:
        from rpy2.robjects.packages import importr
        importr('limma')
    except Exception:                                   # noqa: BLE001 - any R failure means no limma here
        pytest.skip('R limma is not available')


def test_microdissection_segment_t_finds_the_segment_gene():
    _needs_limma()
    rng = np.random.default_rng(0)
    groups = np.array(['S1'] * 3 + ['S2'] * 5 + ['S3'] * 3)
    x = pd.DataFrame(rng.normal(5, .2, (300, 11)), index=[f'g{i}' for i in range(300)])
    x.loc['g0', groups == 'S3'] += 2.
    out = ws.microdissection_segment_t(x, groups, 'S3', covariate=x.mean(axis=1))
    assert out.loc['g0', 'log2_effect'] == pytest.approx(2., abs=.3)
    assert out.loc['g0', 't'] > 8 and np.abs(out.t.drop('g0')).median() < 1.5
    with pytest.raises(ValueError):
        ws.microdissection_segment_t(x, groups, 'S4')


def test_greedy_representatives_absorbs_by_overlap_of_the_smaller_core():
    cores = {'a': list('abcdef'), 'b': list('abcd'), 'c': list('uvwxyz'), 'd': list('aux')}
    reps = ws.greedy_representatives(cores, ['a', 'b', 'c', 'd'])
    assert reps == {'a': ['b'], 'c': ['d']}                   # d overlaps a by 1/3 and c by 2/3 -> absorbed by c
    reps = ws.greedy_representatives(cores, ['a', 'b', 'c', 'd'], overlap=.7)
    assert list(reps) == ['a', 'c', 'd']


def test_thirds_profile_places_the_peak_and_flags_flat():
    u = np.linspace(0, 2.95, 151)
    peak = np.exp(-((u - 2.5) / .1) ** 2)
    out = ws.thirds_profile(peak, u, 2)
    assert out['peak part'] == 1 and out['peak position'] == pytest.approx(2.5, abs=.02) and not out['flat']
    assert np.nansum(out['profile']) == pytest.approx(0, abs=1e-12)
    late = ws.thirds_profile(u, u, 2)
    assert late['peak part'] == 2
    assert ws.thirds_profile(.1 * u, u, 0)['flat']
    assert ws.thirds_profile(np.full(151, np.nan), u, 1)['peak part'] == -1


def test_replicated_spread_counts_only_agreement():
    rng = np.random.default_rng(1)
    a = rng.normal(size=(40, 3))
    assert ws.replicated_spread(a, a) == pytest.approx(((a - a.mean(0)) ** 2).sum())
    assert abs(ws.replicated_spread(a, rng.normal(size=(40, 3)))) < .5 * ws.replicated_spread(a, a)
    a[0, 0] = np.nan
    assert np.isfinite(ws.replicated_spread(a, a))


def test_stratified_permutation_stays_within_strata():
    strata = np.array([0, 0, 1, 1, 1, 2])
    perm = ws.stratified_permutation(strata, np.random.default_rng(2))
    assert sorted(perm) == list(range(6)) and (strata[perm] == strata).all()


def test_position_test_separates_position_grouped_pathways_from_random_ones():
    rng = np.random.default_rng(3)
    u = np.linspace(0, 2.95, 151)
    centres = rng.uniform(2.05, 2.9, 120)
    genes = np.exp(-((u[None, :] - centres[:, None]) / .15) ** 2)
    a = genes + rng.normal(0, .3, genes.shape)
    b = genes + rng.normal(0, .3, genes.shape)
    order = np.argsort(centres)
    grouped = [order[i:i + 10] for i in range(0, 120, 20)]           # pathways = genes with nearby peaks
    random_sets = [rng.choice(120, 10, replace=False) for _ in range(6)]
    strata = np.zeros(120, int)
    real = ws.position_test(a, b, u, grouped, 2, strata, n_perm=200, seed=0)
    noise = ws.position_test(a, b, u, random_sets, 2, strata, n_perm=200, seed=0)
    assert real['p'] < .01 and noise['p'] > .05
