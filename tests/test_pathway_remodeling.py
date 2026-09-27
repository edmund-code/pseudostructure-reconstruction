"""Synthetic checks for notebook 12; no private data or donor-level claims."""
import pytest

np = pytest.importorskip('numpy')
pd = pytest.importorskip('pandas')
pytest.importorskip('scipy')
pytest.importorskip('patsy')
pytest.importorskip('statsmodels')
from scipy.stats import mannwhitneyu
from statsmodels.regression.linear_model import WLS

from pseudospace.pathway_remodeling import (
    fit_nested_trajectories, rank_auc, matched_pathway_tests, matching_strata,
    local_pathway_curves, spatial_descriptions, collapse_programs, gene_covariates,
)
from pseudospace.levelshape import build_ls_designs


def test_nested_effects_hc3_and_specimen_contrasts():
    rng = np.random.default_rng(12)
    position = np.tile(np.linspace(0, 1, 100), 4)
    human = np.repeat([0, 0, 1, 1], 100)
    specimen = np.repeat(['m1', 'm2', 'h1', 'h2'], 100)
    noise = rng.normal(0, .15, (400, 4))
    y = 3 + position[:, None] + noise
    y[:, 0] += 2 * human  # constant level
    y[:, 1] += 4 * human * (position - .5)  # crossover, zero average level
    y[:, 2] += np.repeat([1, -1, .7, -.7], 100)  # specimen nuisance only
    y[:, 3] = 0  # degenerate gene stays finite
    grid = np.linspace(.05, .95, 31)
    fit = fit_nested_trajectories(y, position, human, specimen, grid)
    assert fit['T_level'][0] > 100 * fit['T_spatial'][0]
    assert fit['T_spatial'][1] > 100 * fit['T_level'][1]
    assert fit['T_total'][2] < 5
    assert np.allclose(fit['delta'][0], 2, atol=.12)
    assert fit['delta'][1, 0] < -1 and fit['delta'][1, -1] > 1
    assert all(np.isfinite(v).all() for v in fit.values())
    assert np.allclose(fit['z'], fit['delta'] / fit['se'])
    # Independent implementation: local HC3 must equal statsmodels' sandwich.
    full = build_ls_designs(position, human, fit['knots'])[2]
    contrast_h = (specimen == 'h1').astype(float) - (specimen == 'h2').astype(float)
    contrast_m = (specimen == 'm1').astype(float) - (specimen == 'm2').astype(float)
    x = np.column_stack([full, contrast_m, contrast_h])
    reference = WLS(y[:, 1], x).fit(cov_type='HC3')
    from pseudospace.stats_gam import make_gam_design
    base = make_gam_design(grid, fit['knots'])
    p = base.shape[1]
    contrast = np.zeros((len(grid), x.shape[1]))
    contrast[:, p] = 1
    contrast[:, p + 1:2 * p] = base[:, 1:]
    expected = np.sqrt(np.einsum('ij,jk,ik->i', contrast, reference.cov_params(), contrast))
    assert np.allclose(fit['se'][1], expected)
    # Unbalanced row counts: weights affect both fit and robust uncertainty.
    keep = np.arange(400) >= 50
    uneven = fit_nested_trajectories(y[keep], position[keep], human[keep], specimen[keep], grid)
    full = build_ls_designs(position[keep], human[keep], uneven['knots'])[2]
    x = np.column_stack([full, contrast_m[keep], contrast_h[keep]])
    weights = np.where(specimen[keep] == 'm1', 1.75, .875)
    reference = WLS(y[keep, 1], x, weights=weights).fit(cov_type='HC3')
    base = make_gam_design(grid, uneven['knots'])
    contrast[:, p + 1:2 * p] = base[:, 1:]
    expected = np.sqrt(np.einsum('ij,jk,ik->i', contrast, reference.cov_params(), contrast))
    assert np.allclose(uneven['se'][1], expected)


def test_auc_matching_ties_and_relative_direction():
    values = np.array([1., 2., 2., 4., 8., 9.])
    assert rank_auc(values, [0, 2]) == mannwhitneyu(values[[0, 2]], values[[1, 3, 4, 5]]).statistic / 8
    scores = pd.DataFrame({'T_spatial': values, 'T_level': values[::-1]}, index=list('abcdef'))
    sets = {'top': ['e', 'f'], 'bottom': ['a', 'b']}
    result = matched_pathway_tests(scores, sets, np.zeros(6), n_null=999, seed=4)
    repeat = matched_pathway_tests(scores, sets, np.zeros(6), n_null=999, seed=4)
    pd.testing.assert_frame_equal(result, repeat)
    assert result.p_empirical.between(.001, 1).all()
    assert result.q_empirical.ge(result.p_empirical - 1e-12).all()
    # A completely fixed matched stratum has p=1, never a zero-variance false discovery.
    fixed = matched_pathway_tests(scores, {'fixed': ['a', 'b']}, [0, 0, 1, 1, 1, 1], n_null=99)
    assert (fixed.p_empirical == 1).all() and (fixed.fixed_member_fraction == 1).all()
    z = np.tile(np.array([-1., -2., -5., -6., -7., -8.])[:, None], (1, 5))
    curves = local_pathway_curves(z, list('abcdef'), {'relative': ['a', 'b']}, np.linspace(0, 1, 5))
    assert (curves.direction > 0).all() and (curves.median_member_z < 0).all()
    assert (curves.human_high_fraction == 0).all()
    descriptions = spatial_descriptions(curves)
    assert descriptions.affected_width.iloc[0] == 0
    assert descriptions.direction_sign_changes.iloc[0] == 0
    assert collapse_programs([], {}, {}, curves).empty
    assert collapse_programs(['relative'], {'relative': ['a', 'b']}, {}, curves).program.iloc[0] == 1


def test_validation_and_covariates():
    x = np.array([[0, 1], [1, 1], [0, 0], [1, 1.]])
    frame = gene_covariates(x, [0, 1, 0, 1], [0, 0, 1, 1], ['m', 'm', 'h', 'h'], ['a', 'b'])
    assert np.allclose(frame.detection, [.5, .75])
    assert np.array_equal(matching_strata(frame), matching_strata(frame.copy()))
    with pytest.raises(ValueError, match='both species'):
        fit_nested_trajectories(x, [0, 1, 0, 1], [0] * 4, ['m'] * 4, [0, 1])
    with pytest.raises(ValueError, match='unique'):
        rank_auc([1, 2, 3], [1, 1])
    with pytest.raises(ValueError, match='membership'):
        matched_pathway_tests(pd.DataFrame({'x': [1, 2]}, index=['a', 'b']), {'bad': ['c']}, [0, 0])


def test_saved_full_residuals_do_not_change_the_fit():
    rng = np.random.default_rng(8)
    s = np.tile(np.linspace(0, 1, 60), 4)
    human = np.repeat([0, 0, 1, 1], 60)
    sample = np.repeat(['m1', 'm2', 'h1', 'h2'], 60)
    y = 2 + rng.normal(size=(240, 3)) + s[:, None] * human[:, None]
    grid = np.linspace(.1, .9, 10)
    old = fit_nested_trajectories(y, s, human, sample, grid)
    fit = fit_nested_trajectories(y, s, human, sample, grid, return_residuals=True)
    for key in old:
        assert np.array_equal(old[key], fit[key])
    x = build_ls_designs(s, human, fit['knots'])[2]
    x = np.column_stack([x, (sample == 'm1').astype(float) - (sample == 'm2'),
                        (sample == 'h1').astype(float) - (sample == 'h2')])
    assert np.allclose(fit['residuals'], y - x @ np.linalg.lstsq(x, y, rcond=None)[0])


def test_residual_correlations_are_specimen_balanced_and_exact():
    from pseudospace.pathway_remodeling import residual_pathway_correlations, correlation_adjusted_rank_tests
    rng = np.random.default_rng(91)
    blocks = [rng.normal(size=(300, 4)), rng.normal(size=(40, 4))]
    blocks[0][:, 1] = blocks[0][:, 0] + .1 * blocks[0][:, 1]
    blocks[1][:, 1] = -blocks[1][:, 0] + .5 * blocks[1][:, 1]
    residuals = np.concatenate(blocks)
    sample = np.repeat(['large', 'small'], [300, 40])
    genes, sets = list('abcd'), {'pair': ['a', 'b'], 'triple': ['a', 'b', 'c']}
    audit = residual_pathway_correlations(residuals, sample, genes, sets)
    expected = [np.corrcoef(block[:, :2], rowvar=False)[0, 1] for block in blocks]
    assert np.allclose(audit.loc[audit.pathway_id.eq('pair'), 'rho_specimen'], expected)
    triple_expected = [np.corrcoef(block[:, :3], rowvar=False)[np.triu_indices(3, 1)].mean()
                       for block in blocks]
    assert np.allclose(audit.loc[audit.pathway_id.eq('triple'), 'rho_specimen'], triple_expected)
    scores = pd.Series([4., 3., 2., 1.], index=genes)
    test = correlation_adjusted_rank_tests(scores, sets, audit).iloc[0]
    balanced = np.tanh(np.arctanh(expected).mean())
    assert np.isclose(test.rho_residual, balanced)
    assert not np.isclose(test.rho_residual, np.corrcoef(residuals[:, :2], rowvar=False)[0, 1])
    # Undefined correlations stay explicit; original pathway size is unchanged.
    residuals[:, 1] = 0
    audit = residual_pathway_correlations(residuals, sample, genes, sets)
    tests = correlation_adjusted_rank_tests(scores, sets, audit).set_index('pathway_id')
    missing = tests.loc['pair']
    pair_audit = audit[audit.pathway_id.eq('pair')]
    assert pair_audit.rho_specimen.isna().all() and (pair_audit.n_variable_genes == 1).all()
    assert np.isnan(missing.p_corr) and np.isnan(missing.q_corr) and missing.n_genes == 2
    assert missing.auc == 1
    assert np.isclose(tests.loc['triple', 'q_corr'], min(2 * tests.loc['triple', 'p_corr'], 1))


def test_camera_rank_variance_ties_and_correlation_floor():
    from pseudospace.pathway_remodeling import correlation_adjusted_rank_tests
    genes = list('abcdef')
    sets = {'top': ['e', 'f'], 'bottom': ['a', 'b']}
    def audit(rho):
        return pd.DataFrame([{'pathway_id': p, 'specimen': s, 'n_genes': 2,
                              'n_variable_genes': 2, 'pair_fraction': 1., 'rho_specimen': rho}
                             for p in sets for s in ['mouse', 'human']])
    for values in ([1., 2., 3., 4., 5., 6.], [1., 2., 2., 3., 5., 5.]):
        stats = pd.Series(values, index=genes)
        zero = correlation_adjusted_rank_tests(stats, sets, audit(0)).set_index('pathway_id')
        positive = correlation_adjusted_rank_tests(stats, sets, audit(.4)).set_index('pathway_id')
        negative = correlation_adjusted_rank_tests(stats, sets, audit(-.4)).set_index('pathway_id')
        for pathway, members in sets.items():
            expected = mannwhitneyu(stats[members], stats.drop(members), alternative='greater',
                                    method='asymptotic', use_continuity=True)
            assert np.isclose(zero.loc[pathway, 'p_corr'], expected.pvalue)
            assert np.isclose(zero.loc[pathway, 'U'], expected.statistic)
        assert np.array_equal(zero.auc, positive.auc)
        assert np.allclose(zero.p_corr, negative.p_corr)
        assert positive.loc['top', 'p_corr'] > zero.loc['top', 'p_corr']
        assert (positive.variance_inflation > 1).all()
        assert np.allclose(positive.q_corr, np.minimum(positive.p_corr * [2, 1], 1))
    # limma RELEASE_3_22 rankSumTestWithCorrelation, df=Inf, tied scores, rho=.4.
    assert np.isclose(positive.loc['top', 'p_corr'], 0.066503999329956195, rtol=1e-13)
    flat = correlation_adjusted_rank_tests(pd.Series(1., index=genes), sets, audit(.4))
    assert (flat.auc == .5).all() and (flat.p_corr == 1).all() and (flat.q_corr == 1).all()
    with pytest.raises(ValueError, match='Duplicate'):
        correlation_adjusted_rank_tests(stats, sets, pd.concat([audit(0), audit(0)]))
