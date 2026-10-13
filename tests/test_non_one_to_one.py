import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from pseudospace import non_one_to_one as no


def test_family_clean_allows_cross_hits_only_inside_the_family():
    inv = pd.DataFrame({'gene': ['A1', 'A2', 'B1', 'C1', 'D1'], 'component': [1, 1, 2, 3, -1],
                        'probe_status': ['ambiguous', 'near', 'ambiguous', 'clean', 'ambiguous'],
                        'offtarget_genes': ['A2', np.nan, 'A1', np.nan, 'A1']})
    assert no.family_clean(inv).tolist() == [True, True, False, True, False]


def test_family_sums_and_pairs():
    counts = sparse.csr_matrix(np.array([[1, 2, 3], [0, 5, 1]]))
    out = no.family_sums(counts, ['a', 'b', 'c'], {'F1': ['a', 'c'], 'F2': ['b']})
    assert out.toarray().tolist() == [[4, 2], [1, 5]]
    with pytest.raises(ValueError):
        no.family_sums(counts, ['a', 'b', 'c'], {'F': ['a', 'z']})
    pairs = no.family_pairs({'F1': ['c', 'a', 'b'], 'F2': ['x']})
    assert pairs[['a', 'b']].values.tolist() == [['a', 'b'], ['a', 'c'], ['b', 'c']]


def test_species_zonation_tests_against_the_typical_shape():
    pytest.importorskip('patsy')
    rng = np.random.default_rng(3)
    n = 700
    position = np.concatenate([rng.uniform(0, 1, n), rng.uniform(0, 1, n)])
    specimen = np.repeat(['m1', 'm2'], n)
    library = rng.uniform(800, 1200, 2 * n)
    shared = np.exp(.8 * position)                       # every gene shares this shape (composition drift)
    base_rate = np.exp(rng.uniform(np.log(2e-3), np.log(2e-2), 40))
    rates = np.outer(shared, base_rate)
    rates[:, 0] *= np.exp(3 * (position - .5))           # one zonated gene, rising beyond the shared shape
    counts = sparse.csr_matrix(rng.poisson(rates * library[:, None]))
    knots = np.quantile(position, [.25, .5, .75])
    grid = np.linspace(.02, .98, 41)
    z = no.species_zonation(counts, library, position, specimen, knots, grid)
    genes = [f'g{i}' for i in range(40)]
    calls = no.zonation_calls(z, genes, grid * 3)
    assert calls.loc['g0', 'zonated'] and calls.loc['g0', 'peak'] > 2.5
    assert calls.zonated.sum() <= 2                      # typical genes follow the shape: not zonated
    assert np.nanmedian(np.abs(z['centred']['m1'][1:])) < .3
    shape = z['shape']['m1']
    inside = np.flatnonzero(np.isfinite(shape))
    span = grid[inside[-1]] - grid[inside[0]]
    assert shape[inside[-1]] - shape[inside[0]] == pytest.approx(.8 * span / np.log(2), abs=.25)
    again = no.species_zonation(counts[:, :3], library, position, specimen, knots, grid, shapes=z['shape'])
    np.testing.assert_allclose(again['centred']['m2'][:3], z['centred']['m2'][:3], atol=.05, equal_nan=True)
    np.testing.assert_allclose(again['lr']['m1'], z['lr']['m1'][:3], rtol=1e-6)
