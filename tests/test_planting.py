import numpy as np
import pytest

from pseudospace import planting as pl


def test_choose_sets_respects_size_overlap_and_unique_genes():
    rng = np.random.default_rng(0)
    genes = [f'g{i}' for i in range(2000)]
    sets = {f's{j}': genes[j * 40:j * 40 + 50] for j in range(45)}       # neighbours overlap by 10 genes
    sets['tiny'] = genes[:5]
    chosen = pl.choose_sets(sets, rng, k=15, sizes=(20, 150), max_jaccard=.1, frac=.6, min_members=10)
    assert len(chosen) == 15 and 'tiny' not in chosen
    planted = [g for v in chosen.values() for g in v]
    assert len(planted) == len(set(planted))
    for a in chosen:
        for b in chosen:
            if a < b:
                A, B = set(sets[a]), set(sets[b])
                assert len(A & B) / len(A | B) < .1
    assert all(len(v) == 30 for v in chosen.values())
    with pytest.raises(ValueError):
        pl.choose_sets({'a': genes[:30]}, rng, k=2)


def test_effect_shapes():
    x = np.linspace(0, 1, 1001)
    planted = x > -1
    assert np.allclose(pl.effect_shape(x, 'level', planted), 1)
    loc = pl.effect_shape(x, 'located', planted)
    assert loc[0] < .01 and loc[-1] > .99 and abs(loc[650] - .5) < .01
    for kind in ('shape', 'bump'):
        assert abs(pl.effect_shape(x, kind, planted).mean()) < 1e-12
    b = pl.effect_shape(x, 'bump', planted)
    assert np.argmax(b) == 500
    with pytest.raises(ValueError):
        pl.effect_shape(x, 'nope', planted)


def test_plant_counts_scales_expectation():
    rng = np.random.default_rng(1)
    c = rng.poisson(20, size=(20000, 3)).astype(float)
    m = np.tile([.5, 1., 2.], (20000, 1))
    out = pl.plant_counts(c, m, rng)
    assert np.allclose(out.mean(axis=0) / c.mean(axis=0), [.5, 1., 2.], atol=.02)
    assert np.array_equal(out[:, 1], c[:, 1]) and np.all(out == np.round(out))
    assert (out[:, 0] <= c[:, 0]).all() and (out[:, 2] >= c[:, 2]).all()
