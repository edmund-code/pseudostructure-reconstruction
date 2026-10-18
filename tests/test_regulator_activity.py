import numpy as np
import pandas as pd
import pytest
from scipy import stats

from pseudospace import regulator_activity as ra


def test_collectri_network_collapses_dimers_signs_and_maps():
    rows = pd.DataFrame({
        'source_genesymbol': ['HNF4A', 'HNF4A', 'NFKB1_RELA', 'REL_RELB', 'JUN_FOS', 'SMAD3_SMAD4', 'MYC'],
        'target_genesymbol': ['SLC5A2', 'CYP7B1', 'CCL2', 'CCL2', 'MMP1', 'SERPINE1', 'NOPE'],
        'is_stimulation': [True, False, True, False, True, True, True],
        'is_inhibition': [False, True, False, True, False, False, False]})
    h2m = {'SLC5A2': 'Slc5a2', 'CYP7B1': 'Cyp7b1', 'CCL2': 'Ccl2', 'MMP1': 'Mmp1a', 'SERPINE1': 'Serpine1'}
    net = ra.collectri_network(rows, h2m).set_index(['source', 'target']).weight
    assert net[('HNF4A', 'Slc5a2')] == 1 and net[('HNF4A', 'Cyp7b1')] == -1
    assert net[('NFKB', 'Ccl2')] == 1                       # one activating copy among the dimers keeps +1
    assert net[('AP1', 'Mmp1a')] == 1
    assert ('SMAD3_SMAD4', 'Serpine1') not in net.index and 'MYC' not in net.index.get_level_values(0)


def test_progeny_network_keeps_top_genes():
    recs = []
    for i, (path, gene, w, p) in enumerate([('TNFa', 'A', 2., .01), ('TNFa', 'B', -1., .02), ('TNFa', 'C', 5., .9),
                                            ('p53', 'A', 1., .5)]):
        recs += [(i, gene, 'pathway', path), (i, gene, 'weight', w), (i, gene, 'p_value', p)]
    rows = pd.DataFrame(recs, columns=['record_id', 'genesymbol', 'label', 'value'])
    net = ra.progeny_network(rows, {'A': 'a', 'B': 'b', 'C': 'c'}, top=2)
    assert sorted(map(tuple, net[net.source == 'TNFa'][['target', 'weight']].to_numpy())) == [('a', 2.), ('b', -1.)]
    assert net[net.source == 'p53'].target.tolist() == ['a']


def test_ulm_equals_ols_t_and_null_is_centred():
    rng = np.random.default_rng(0)
    genes = [f'g{i}' for i in range(300)]
    net = pd.DataFrame({'source': ['s1'] * 30 + ['s2'] * 5, 'target': genes[:30] + genes[30:35],
                        'weight': np.r_[np.ones(20), -np.ones(10), np.ones(5)]})
    W, sources = ra.weight_matrix(net, genes, min_targets=10)
    assert sources == ['s1'] and W.shape == (300, 1)
    y = rng.normal(size=300)
    y[:20] += 1.
    y[20:30] -= 1.
    t = ra.ulm(W, y)[0, 0]
    fit = stats.linregress(W.toarray()[:, 0], y)
    assert t == pytest.approx(fit.slope / fit.stderr)
    strata = np.repeat(np.arange(3), 100)
    perms = ra.stratified_permutations(strata, 400, seed=1)
    assert all(np.array_equal(np.sort(p[strata == s]), np.flatnonzero(strata == s)) for p in perms[:5] for s in range(3))
    obs, z = ra.ulm_null(W, y, perms)
    assert obs[0] == pytest.approx(t) and z[0] > 4
    _, z0 = ra.ulm_null(W, rng.normal(size=300), perms)
    assert abs(z0[0]) < 3.5
