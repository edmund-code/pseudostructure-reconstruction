import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scipy.stats import spearmanr

from pseudospace.species_level import (LN2, band, bootstrap_spearman, clears_band, equal_segment_log_cpm,
                                       family_sums, grid_centred_level, limma_trend_species, log_cpm, pair_levels,
                                       pair_summary, segment_pseudobulk, species_difference, wilson_interval)


def test_grid_centred_level_removes_a_shared_position_dependent_shift():
    rng = np.random.default_rng(0)
    grid = np.linspace(0, 1, 21)
    effect = rng.normal(0, 1, 301)
    effect[0] = 3.0
    shared = .8 * np.sin(3 * grid) - .4              # depth / composition shift, the same for every gene
    delta = effect[:, None] + shared[None, :]
    level = grid_centred_level(pd.DataFrame(delta, index=[f'g{i}' for i in range(301)]))
    np.testing.assert_allclose(level.to_numpy(), (effect - np.median(effect)) / LN2, atol=1e-12)
    assert level.index[0] == 'g0'


def test_grid_centred_level_reference_weights_and_nonfinite_rows():
    grid = np.linspace(0, 1, 5)
    delta = np.array([[0., 0, 0, 0, 0], [1, 1, 1, 1, 1], [2, 2, 2, 2, 2], [50, 50, 50, 50, 50],
                      [np.nan, 1, 1, 1, 1]])
    # Reference excludes the outlier gene 3: median of genes 0-2 is 1 at every point.
    level = grid_centred_level(delta, reference=[True, True, True, False, True])
    np.testing.assert_allclose(level[:4], np.array([-1, 0, 1, 49]) / LN2)
    assert np.isnan(level[4])
    # Weights: a gene that differs only at the last grid point, weighted fully there.
    delta2 = np.zeros((3, 5))
    delta2[0, -1] = 2.
    w = np.array([0, 0, 0, 0, 1.])
    np.testing.assert_allclose(grid_centred_level(delta2, weights=w)[0], 2 / LN2)
    np.testing.assert_allclose(grid_centred_level(delta2)[0], 2 / 5 / LN2)
    with pytest.raises(ValueError):
        grid_centred_level(delta2, weights=np.ones(4))
    assert len(grid) == 5


def test_pair_levels_centre_each_pair_and_mask_invalid_specimens():
    rng = np.random.default_rng(1)
    base = rng.normal(-8, 1, (50, 11))
    effect = np.zeros(50)
    effect[:3] = [2., -1.5, 1.]
    curves = {'h1': base + effect[:, None] + .3, 'h2': base + effect[:, None] - .2 + .05 * np.arange(11),
              'm1': base, 'm2': base + .7}
    valid = {k: np.ones(50, bool) for k in curves}
    valid['m2'][1] = False
    pairs = [(h, m) for h in ('h1', 'h2') for m in ('m1', 'm2')]
    levels = pair_levels(curves, pairs, valid=valid, index=[f'g{i}' for i in range(50)])
    assert list(levels.columns) == ['h1 - m1', 'h1 - m2', 'h2 - m1', 'h2 - m2']
    np.testing.assert_allclose(levels.loc['g0', 'h1 - m1'], 2 / LN2, atol=1e-12)
    np.testing.assert_allclose(levels.loc['g2'], 1 / LN2, atol=1e-12)          # every pair, own centring
    assert levels.loc['g1', ['h1 - m2', 'h2 - m2']].isna().all()
    summary = pair_summary(levels)
    assert summary.loc['g0', 'pairs_same_sign'] and summary.loc['g0', 'pairs'] == 4
    assert not summary.loc['g1', 'pairs_same_sign'] and summary.loc['g1', 'pairs'] == 2
    np.testing.assert_allclose(summary.loc['g1', ['pair_min', 'pair_max']].astype(float), -1.5 / LN2, atol=1e-12)
    with pytest.raises(ValueError):
        pair_levels(curves, [('h1', 'm3')])


def test_segment_pseudobulk_and_equal_segment_log_cpm_by_hand():
    counts = np.array([[1, 0], [3, 2], [0, 5], [2, 2], [4, 0], [1, 1], [0, 3]], float)
    library = np.array([10, 20, 30, 10, 20, 30, 40.])
    unit = np.array(['a', 'a', 'a', 'b', 'b', 'b', 'b'])
    segment = np.array(['S1', 'S1', 'S2', 'S1', 'S2', 'S2', 'S1'])
    units, sums, totals = segment_pseudobulk(sparse.csr_matrix(counts), library, unit, segment, segments=('S1', 'S2'))
    assert units == ['a', 'b']
    np.testing.assert_array_equal(sums[0], [[4, 2], [0, 5]])
    np.testing.assert_array_equal(sums[1], [[2, 5], [5, 1]])
    np.testing.assert_array_equal(totals, [[30, 30], [50, 50]])
    values = equal_segment_log_cpm(sums, totals)
    expected_a0 = np.mean([np.log2(4.5 / 31 * 1e6), np.log2(.5 / 31 * 1e6)])
    np.testing.assert_allclose(values[0, 0], expected_a0)
    weighted = equal_segment_log_cpm(sums, totals, weights=[2, 1])
    np.testing.assert_allclose(weighted, (2 * log_cpm(sums, totals)[:, 0] + log_cpm(sums, totals)[:, 1]) / 3)
    with pytest.raises(ValueError):
        segment_pseudobulk(counts, library, unit, segment, segments=('S1', 'S2', 'S3'))
    np.testing.assert_allclose(species_difference(values, [True, False]), values[0] - values[1])
    with pytest.raises(ValueError):
        species_difference(values, [True, True])


def test_family_sums_add_member_columns_and_reject_unknown_members():
    counts = sparse.csr_matrix(np.array([[1, 2, 3, 4], [0, 1, 0, 5]], float))
    sums, names = family_sums(counts, ['A1', 'A2', 'B1', 'C1'], {'A': ['A1', 'A2'], 'B': ['B1'], 'C': ['C1', 'A1']})
    assert names == ['A', 'B', 'C']
    np.testing.assert_array_equal(sums.toarray(), [[3, 3, 5], [1, 0, 5]])
    with pytest.raises(ValueError):
        family_sums(counts, ['A1', 'A2', 'B1', 'C1'], {'A': ['A1', 'Z9']})


def test_yardsticks():
    lo, hi = band(np.r_[np.arange(101.), np.nan])
    assert (lo, hi) == (2.5, 97.5)
    np.testing.assert_array_equal(clears_band([99, 99, 1, 50], [1, -1, -1, 0], (lo, hi)), [True, False, True, False])
    low, high = wilson_interval(8, 10)
    assert 0.49 < low < 0.5 and 0.94 < high < 0.95             # 0.490, 0.943
    assert np.isnan(wilson_interval(0, 0)[0])
    rng = np.random.default_rng(3)
    x = rng.normal(size=400)
    y = x + rng.normal(size=400)
    y[:5] = np.nan
    out = bootstrap_spearman(x, y, n_boot=500, seed=1)
    keep = np.isfinite(y)
    assert out['n'] == 395
    np.testing.assert_allclose(out['rho'], spearmanr(x[keep], y[keep]).correlation)
    assert out['low'] < out['rho'] < out['high'] and out['high'] - out['low'] < .2
    assert bootstrap_spearman(x, y, n_boot=500, seed=1) == out


def _limma_available():
    pytest.importorskip('rpy2')
    from rpy2.robjects.packages import isinstalled
    if not isinstalled('limma'):
        pytest.skip('R package limma is not installed')


def _simulated_donors(seed=4):
    rng = np.random.default_rng(seed)
    n_genes, case = 2000, np.array([True] * 7 + [False] * 6)
    mean = rng.uniform(1, 12, n_genes)
    effect = np.zeros(n_genes)
    effect[:40] = 3.0
    effect[40:80] = -3.0
    noise_sd = .3 + 1.5 / (1 + mean)                     # a mean-variance trend, as on the log-CPM scale
    values = mean + .4 + effect * case[:, None] + rng.normal(0, 1, (13, n_genes)) * noise_sd
    values[case] += .4                                   # a constant offset that centring must remove
    return values, case, effect


def test_limma_trend_species_recovers_planted_genes_and_centres_the_coefficient():
    _limma_available()
    values, case, effect = _simulated_donors()
    genes = [f'g{i}' for i in range(values.shape[1])]
    out = limma_trend_species(values, case, genes=genes)
    np.testing.assert_allclose(out.coefficient, values[case].mean(axis=0) - values[~case].mean(axis=0), atol=1e-10)
    assert abs(np.median(out.centred)) < 1e-12
    strong = out.q.le(.05) & out.centred.abs().ge(1)
    assert strong.iloc[:80].mean() > .95 and strong.iloc[80:].mean() < .01
    assert (np.sign(out.centred.iloc[:80]) == np.sign(effect[:80])).all()
    np.testing.assert_allclose(out.se, out.stdev_unscaled * np.sqrt(out.s2_post))
    np.testing.assert_allclose(out.stdev_unscaled, np.sqrt(1 / 7 + 1 / 6))
    assert out.index[0] == 'g0'


def test_limma_without_trend_matches_squeeze_var_moderation():
    _limma_available()
    from pseudospace.pathway_pipelines import squeeze_var
    values, case, _ = _simulated_donors(seed=5)
    out = limma_trend_species(values, case, trend=False, robust=False)
    fitted = np.where(case[:, None], values[case].mean(axis=0), values[~case].mean(axis=0))
    s2 = ((values - fitted) ** 2).sum(axis=0) / (len(case) - 2)
    np.testing.assert_allclose(out.sigma ** 2, s2, rtol=1e-10)
    moderated = squeeze_var(s2, len(case) - 2)
    np.testing.assert_allclose(out.s2_post, moderated.var_post, rtol=1e-8)
    with pytest.raises(ValueError):
        limma_trend_species(values[:3], case[:3])
