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
