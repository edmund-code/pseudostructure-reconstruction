import numpy as np
import pytest

from pseudospace.pathway_calibration import contrast_designs
from pseudospace.stats_gam import gam_internal_knots
from pseudospace.within_segment_nb import (beyond_step_designs, nb_segment_slopes, nb_specimen_curves, nested_lr,
                                           permuted_rows, pseudobulk_log_rates, step_ladder_designs,
                                           within_species_designs)

rank = np.linalg.matrix_rank


def _cohort(n=300, seed=0, noise=.1):
    rng = np.random.default_rng(seed)
    specimen = np.repeat(['m1', 'm2', 'h1', 'h2'], n)
    human = np.isin(specimen, ['h1', 'h2']).astype(float)
    s = rng.uniform(0, 1, len(specimen))
    segment = np.where(s < 1 / 3, 'PT-S1', np.where(s < 2 / 3, 'PT-S2', 'PT-S3'))
    flip = rng.uniform(size=len(s)) < noise
    segment = np.where(flip, rng.choice(['PT-S1', 'PT-S2', 'PT-S3'], len(s)), segment)
    return rng, specimen, human, s, segment


def _df(designs, comparisons):
    return {k: rank(designs[b]) - rank(designs[a]) for k, (a, b) in comparisons.items()}


def test_beyond_step_df_and_nesting():
    _, specimen, human, s, segment = _cohort()
    knots = gam_internal_knots(s, basis_df=6)
    relabel = np.isin(specimen, ['m1', 'h1']).astype(float)
    for group, nuisance, df2 in ((human, None, 12), (relabel, human, 6)):
        designs, comparisons, weights = beyond_step_designs(s, group, specimen, knots, segment, nuisance=nuisance)
        assert _df(designs, comparisons) == {'numerator': 6, 'denominator': df2}
        assert rank(np.column_stack([designs['reduced'], designs['full']])) == rank(designs['full'])
        assert np.all(weights > 0)


def test_case_only_steps_reproduce_contrast_designs():
    _, specimen, human, s, segment = _cohort()
    knots = gam_internal_knots(s, basis_df=6)
    designs, _, _ = beyond_step_designs(s, human, specimen, knots, segment, shared_steps=False)
    old, _ = contrast_designs(s, human, specimen, knots, segments=segment)
    np.testing.assert_array_equal(designs['reduced'], old['discrete'])
    np.testing.assert_array_equal(designs['full'], old['union'])
    shared, _, _ = beyond_step_designs(s, human, specimen, knots, segment)
    assert rank(shared['reduced']) == rank(old['discrete']) + 2


def test_ladder_df():
    _, specimen, human, s, segment = _cohort()
    knots = gam_internal_knots(s, basis_df=6)
    relabel = np.isin(specimen, ['m1', 'h1']).astype(float)
    bins = np.clip((s * 6).astype(int), 0, 5).astype(str)
    for labels, df1 in ((segment, 2), (bins, 5)):
        for group, nuisance, factor in ((human, None, 2), (relabel, human, 1)):
            designs, comparisons, _ = step_ladder_designs(s, group, specimen, knots, labels, nuisance=nuisance)
            assert _df(designs, comparisons) == {'numerator': df1, 'denominator': factor * df1}


def test_within_species_df():
    _, specimen, human, s, segment = _cohort()
    m = human == 0
    knots = gam_internal_knots(s[m], basis_df=6)
    for segments in (None, segment[m]):
        designs, comparisons, weights = within_species_designs(s[m], specimen[m], knots, segments)
        assert _df(designs, comparisons) == {'numerator': 6, 'denominator': 6}
        np.testing.assert_allclose(weights.sum(), m.sum())


def _counts(rng, log_rate, library, alpha=.05):
    mu = library * np.exp(log_rate)
    lam = rng.gamma(1 / alpha, alpha * mu)
    return rng.poisson(lam)


def test_beyond_steps_lr_small_for_pure_steps():
    rng, specimen, human, s, segment = _cohort(n=400, noise=0)
    knots = gam_internal_knots(s, basis_df=6)
    library = rng.uniform(2000, 6000, len(s))
    code = np.searchsorted(['PT-S1', 'PT-S2', 'PT-S3'], segment)
    steps_only = np.column_stack([_counts(rng, -6 + human * code * .8 + code * .3, library) for _ in range(20)])
    gradient = np.column_stack([_counts(rng, -6 + human * 2.5 * s + code * .3, library) for _ in range(20)])
    designs, comparisons, weights = beyond_step_designs(s, human, specimen, knots, segment)
    alpha = np.full(20, .05)
    a = nested_lr(steps_only, library, designs, comparisons, weights, alpha)
    b = nested_lr(gradient, library, designs, comparisons, weights, alpha)
    assert a['df1'] == 6 and a['df2'] == 12
    assert np.median(b['numerator']) > 5 * np.median(a['numerator'])


def test_segment_slopes_sign_and_size():
    rng = np.random.default_rng(3)
    n = 600
    specimen = np.repeat(['a', 'b'], n)
    s = np.tile(np.linspace(0, 1, n), 2)
    segment = np.where(s < .5, 'PT-S1', 'PT-S2')
    library = np.full(len(s), 4000.)
    up = _counts(rng, -5 + 1.5 * s, library)
    flat = _counts(rng, -5 + 0 * s, library)
    down = _counts(rng, -5 - 1.5 * s, library)
    out = nb_segment_slopes(np.column_stack([up, flat, down]), library, s, specimen, segment, np.full(3, .05),
                            order=('PT-S1', 'PT-S2', 'PT-S3'))
    for name in ('a', 'b'):
        t = out[name]
        assert t.loc[0, 'PT-S1'] > 3 and t.loc[2, 'PT-S1'] < -3 and abs(t.loc[1, 'PT-S1']) < 3
        assert t['PT-S3'].isna().all()


def test_specimen_curves_and_pseudobulk():
    rng = np.random.default_rng(4)
    n = 2000
    specimen = np.repeat(['a', 'b'], n // 2)
    s = rng.uniform(0, 1, n)
    library = rng.uniform(3000, 5000, n)
    truth = -4 + np.sin(3 * s)
    y = _counts(rng, truth, library, alpha=.01)[:, None]
    knots = gam_internal_knots(s, basis_df=6)
    grid = np.linspace(.05, .95, 31)
    curves = nb_specimen_curves(y, library, s, specimen, grid, knots, np.array([.01]))
    for name in ('a', 'b'):
        np.testing.assert_allclose(curves[name][0], -4 + np.sin(3 * grid), atol=.1)
    segment = np.where(s < .5, 'PT-S1', 'PT-S2')
    rates = pseudobulk_log_rates(y, library, specimen, segment)
    m = (specimen == 'a') & (segment == 'PT-S1')
    assert rates['a'][0, 0] == pytest.approx(np.log(y[m].sum() / library[m].sum()))
    assert np.isnan(rates['a'][0, 2])


def test_permuted_rows_stay_in_block():
    blocks = np.array(list('aabbbcc'))
    index = permuted_rows(blocks, np.random.default_rng(0))
    assert (blocks[index] == blocks).all() and sorted(index) == list(range(7))
