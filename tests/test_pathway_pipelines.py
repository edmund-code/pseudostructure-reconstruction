import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scipy.special import ndtri_exp
from scipy.stats import chi2

from pseudospace.nb_gam import nb_group_difference_curves, nb_species_trajectories, nested_nb_lr
from pseudospace.pathway_calibration import COMPARISONS, contrast_designs
from pseudospace.pathway_pipelines import (camera_pr, clustering_first_inputs, common_grid, log_f_sf, moderated_f,
                                           over_representation, partition_groups, rank_first_gsea,
                                           rank_normal_scores, single_split_design, specimen_ranges,
                                           specimen_rate_floor, specimen_shape_designs, specimen_shape_lr,
                                           squeeze_var, standard_gene_model)
from pseudospace.stats_gam import gam_internal_knots, make_gam_design

SPECIMENS = ['m1', 'm2', 'h1', 'h2']
SHAPE = 0                                   # human-only smooth shape: a species x position gene


def _truth(position, specimen):
    """Natural-log rate per specimen and gene: specimen offsets, a common shape, and SHAPE's human-only bump."""
    human = np.isin(specimen, ['h1', 'h2']).astype(float)
    offset = pd.Series({'m1': 0., 'm2': .15, 'h1': .4, 'h2': .25})[specimen].to_numpy()
    common = .6 * position - .5 * position ** 2
    rate = np.log(2e-2) + offset[:, None] + common[:, None] + np.zeros((1, 6))
    rate[:, SHAPE] += 1.2 * human * np.sin(2 * np.pi * position)
    rate[:, 1:] += np.array([0., .3, -.4, .5, -.2])          # null genes: species offset only
    return rate


@pytest.fixture(scope='module')
def cohort():
    rng = np.random.default_rng(7)
    specimen = np.repeat(SPECIMENS, 250)
    position = rng.uniform(0, 1, len(specimen))
    position[specimen == 'm2'] = rng.uniform(.25, 1, 250)    # one specimen covers only part of the range
    library = rng.lognormal(np.log(3000), .3, len(specimen))
    log_mu = np.log(library)[:, None] + _truth(position, specimen)
    alpha = .05
    counts = rng.poisson(rng.gamma(1 / alpha, np.exp(log_mu) * alpha)).astype(float)
    human = np.isin(specimen, ['h1', 'h2'])
    genes = [f'g{i}' for i in range(counts.shape[1])]
    species = standard_gene_model(sparse.csr_matrix(counts), library, position, human.astype(float), specimen,
                                  genes=genes)
    return dict(counts=counts, library=library, position=position, specimen=specimen, human=human, genes=genes,
                species=species)


def test_species_split_finds_the_shape_gene_and_matches_nb_species_trajectories(cohort):
    model = cohort['species']
    t = model.statistics
    assert t.converged.all() and not t.separated.any()
    nulls = t.T_spatial.drop('g0')
    assert t.T_spatial['g0'] > 50 * nulls.max() and nulls.max() < chi2.ppf(.999, 6)
    np.testing.assert_allclose(model.grid, common_grid(cohort['position'], cohort['specimen']))
    assert len(model.grid) == 101 and model.grid[0] >= specimen_ranges(cohort['position'], cohort['specimen']).lo.max()
    ref = nb_species_trajectories(cohort['counts'], cohort['library'], cohort['position'], cohort['human'],
                                  cohort['specimen'], model.grid, model.dispersion)
    for key in ('T_level', 'T_spatial', 'T_total'):
        np.testing.assert_allclose(t[key], ref[key], rtol=1e-6, atol=1e-6)
    for mine, theirs in ((model.reference, 'mouse'), (model.case, 'human'), (model.delta, 'delta'), (model.se, 'se')):
        np.testing.assert_allclose(mine, ref[theirs], rtol=1e-5, atol=1e-6)
    assert dict(model.groups) == {'h1': 1, 'h2': 1, 'm1': 0, 'm2': 0}


def test_specimen_curves_recover_the_planted_rates_and_are_nan_outside_each_range(cohort):
    model = cohort['species']
    for name, curve in model.specimen_curves.items():
        truth = _truth(model.grid, np.repeat(name, len(model.grid))).T
        assert np.isfinite(curve).all()                       # the common grid lies inside every specimen's range
        assert np.abs(curve - truth).max() < .3 and np.median(np.abs(curve - truth)) < .05
    assert model.specimen_converged.all().all() and not model.specimen_separated.any().any()
    wide = np.linspace(.02, .98, 49)
    model = standard_gene_model(cohort['counts'], cohort['library'], cohort['position'], cohort['human'],
                                cohort['specimen'], wide, dispersion=cohort['species'].dispersion)
    ranges = specimen_ranges(cohort['position'], cohort['specimen'])
    for name, curve in model.specimen_curves.items():
        outside = (wide < ranges.lo[name]) | (wide > ranges.hi[name])
        np.testing.assert_array_equal(np.isnan(curve).all(axis=0), outside)
        assert np.isfinite(curve[:, ~outside]).all()
    assert np.isnan(model.specimen_curves['m2'][:, wide < .2]).all()


def test_relabeled_partition_absorbs_the_species_shape_through_the_nuisance(cohort):
    parts = partition_groups(cohort['specimen'], cohort['human'])
    assert list(parts) == ['species', 'swap:h1+m1', 'swap:h1+m2']
    np.testing.assert_array_equal(parts['species'][0], cohort['human'])
    assert parts['species'][1] is None
    group, nuisance = parts['swap:h1+m1']
    np.testing.assert_array_equal(group, np.isin(cohort['specimen'], ['h1', 'm1']))
    np.testing.assert_array_equal(nuisance, cohort['human'])
    args = (cohort['counts'], cohort['library'], cohort['position'], group, cohort['specimen'])
    relabeled = standard_gene_model(*args, nuisance_shape=nuisance, dispersion=cohort['species'].dispersion)
    assert relabeled.statistics.T_spatial.max() < chi2.ppf(.999, 6)
    # Estimated afresh, the dispersion is the species split's, so nothing changes.
    np.testing.assert_allclose(standard_gene_model(*args, nuisance_shape=nuisance).statistics.T_spatial,
                               relabeled.statistics.T_spatial, rtol=1e-6)
    # Without the nuisance the species shape leaks into the relabeled contrast (m2's partial range unbalances it).
    leaky = standard_gene_model(*args, dispersion=cohort['species'].dispersion).statistics.T_spatial.iloc[SHAPE]
    assert leaky > chi2.ppf(.999, 6) and leaky > 5 * relabeled.statistics.T_spatial.iloc[SHAPE]
    # Each relabeled group holds one specimen of each species, so the nuisance cancels in delta: nb_gam's curves.
    ref = nb_group_difference_curves(*args, relabeled.knots, relabeled.grid, relabeled.dispersion,
                                     nuisance_shape=nuisance)
    np.testing.assert_allclose(relabeled.statistics.T_spatial, ref['T_spatial'], rtol=1e-6)
    # Different IRLS starts: R glm's deviance criterion pins coefficients to about 1e-6.
    np.testing.assert_allclose(relabeled.delta, ref['delta'], atol=1e-5)
    np.testing.assert_allclose(relabeled.se, ref['se'], rtol=1e-4)


def _ranking(seed=3):
    rng = np.random.default_rng(seed)
    genes = [f'G{i:04d}' for i in range(1500)]
    statistic = pd.Series(rng.chisquare(6, len(genes)), index=genes)
    order = statistic.sort_values(ascending=False).index
    sets = {'planted': list(order[:25]) + list(rng.choice(order[25:], 5, replace=False))}
    sets.update({f'random{k}': list(rng.choice(genes, 30, replace=False)) for k in range(30)})
    return statistic, sets


def test_rank_first_gsea_calls_a_planted_set_and_is_deterministic():
    pytest.importorskip('gseapy')
    statistic, sets = _ranking()
    sets['too_small'] = list(statistic.index[:5])
    out = rank_first_gsea(statistic, sets, permutation_num=200, seed=1, threads=1)
    assert list(out.columns) == ['pathway_id', 'n', 'ES', 'NES', 'p', 'q', 'fwer_p', 'leading_edge', 'p_one_sided',
                                 'call']
    assert list(out.pathway_id) == [p for p in sets if p != 'too_small']
    table = out.set_index('pathway_id')
    assert table.call['planted'] and table.n['planted'] == 30 and table.NES['planted'] > 0
    assert set(table.leading_edge['planted']) <= set(sets['planted'])
    assert not table.call.drop('planted').any()
    again = rank_first_gsea(statistic.sample(frac=1, random_state=0), sets, permutation_num=200, seed=1, threads=1)
    pd.testing.assert_frame_equal(out, again)
    normal = rank_first_gsea(statistic, sets, metric='rank_normal', permutation_num=200, seed=1, threads=1)
    assert normal.set_index('pathway_id').loc['planted', 'call']


def test_rank_normal_scores_are_ordered_normal_quantiles_with_shared_ties():
    scores = rank_normal_scores(pd.Series([3., 1., 2., 2.], index=list('abcd')))
    assert scores['a'] == pytest.approx(-scores['b']) and scores['c'] == scores['d'] == pytest.approx(0.)
    assert scores['a'] > scores['c'] > scores['b']
    with pytest.raises(ValueError, match='finite'):
        rank_normal_scores([1., np.nan])


def _limma_available():
    pytest.importorskip('rpy2')
    from rpy2.robjects.packages import isinstalled
    if not isinstalled('limma'):
        pytest.skip('R package limma is not installed')


def _camera_formula(statistic, members, cor):
    """limma's cameraPR (use.ranks = FALSE): two-sample t of the set mean, set variance inflated by 1 + (m - 1) cor."""
    from scipy.stats import t as t_dist
    s = statistic.to_numpy()
    idx = statistic.index.get_indexer(sorted(set(members) & set(statistic.index)))
    g, m = len(s), len(idx)
    delta = g / (g - m) * (s[idx].mean() - s.mean())
    pooled = ((g - 1) * s.var(ddof=1) - delta ** 2 * m * (g - m) / g) / (g - 2)
    t = delta / np.sqrt(pooled * ((1 + (m - 1) * cor) / m + 1 / (g - m)))
    up, down = t_dist.sf(t, g - 2), t_dist.cdf(t, g - 2)
    return 2 * min(up, down), up


def test_camera_pr_reproduces_limmas_formula_and_calls_a_planted_set():
    _limma_available()
    statistic, sets = _ranking()
    statistic = rank_normal_scores(statistic)
    sets['partial'] = list(statistic.index[200:212]) + ['not_a_gene']
    sets['one_member'] = [statistic.index[0], 'not_a_gene']
    out = camera_pr(statistic, sets)
    assert list(out.columns) == ['pathway_id', 'n', 'direction', 'p', 'fdr', 'p_one_sided', 'call']
    assert list(out.pathway_id) == [p for p in sets if p != 'one_member']
    table = out.set_index('pathway_id')
    assert table.n['partial'] == 12
    for pathway in ('planted', 'random0', 'partial'):
        two_sided, upper = _camera_formula(statistic, sets[pathway], .01)
        assert table.p[pathway] == pytest.approx(two_sided, rel=1e-8)
        assert table.p_one_sided[pathway] == pytest.approx(upper, rel=1e-8)
    assert table.call['planted'] and table.direction['planted'] == 'Up' and table.call.sum() == 1
    ranks = camera_pr(statistic, sets, use_ranks=True).set_index('pathway_id')
    assert ranks.call['planted'] and not ranks.p.equals(table.p)
    stronger = camera_pr(statistic, sets, inter_gene_cor=.2).set_index('pathway_id')
    assert stronger.p['planted'] > table.p['planted']           # more correlation, less evidence


@pytest.fixture(scope='module')
def wide_cohort():
    """40 genes: one human-only bump, one gene absent from h2, the rest species offsets only."""
    rng = np.random.default_rng(11)
    specimen = np.repeat(SPECIMENS, 200)
    position = rng.uniform(0, 1, len(specimen))
    human = np.isin(specimen, ['h1', 'h2'])
    library = rng.lognormal(np.log(3000), .3, len(specimen))
    base = np.log(rng.uniform(5e-4, 2e-2, 40))
    rate = base[None, :] + (.4 * human[:, None] * rng.normal(0, 1, 40)[None, :]) + .5 * position[:, None]
    rate[:, 0] += 1.5 * human * np.sin(2 * np.pi * position)
    log_mu = np.log(library)[:, None] + rate
    counts = rng.poisson(rng.gamma(20, np.exp(log_mu) / 20)).astype(float)
    counts[specimen == 'h2', 1] = 0.                               # absent from one specimen: a separated fit
    genes = [f'g{i}' for i in range(40)]
    model = standard_gene_model(sparse.csr_matrix(counts), library, position, human.astype(float), specimen,
                                genes=genes)
    return dict(counts=counts, library=library, position=position, specimen=specimen, human=human, model=model)


def test_clustering_first_inputs_carry_the_shared_curves_into_the_external_pipeline(wide_cohort):
    pytest.importorskip('pseudospace_reconstruction')
    from pseudospace_reconstruction.compare import CompareParams, Design, compare_groups

    c, model = wide_cohort, wide_cohort['model']
    axis, fit = clustering_first_inputs(model, c['position'], c['human'], c['specimen'], c['library'],
                                        zone_names=['A', 'B'], zone_cuts=[.5])
    floor = specimen_rate_floor(c['library'], c['specimen'])
    assert floor['h2'] == pytest.approx(np.log(.5 / c['library'][c['specimen'] == 'h2'].sum()))
    assert list(fit.specimen_names) == sorted(SPECIMENS) and fit.specimen.shape == (4, 40, len(model.grid))
    np.testing.assert_allclose(axis.grid, model.grid)
    assert axis.weights.sum() == pytest.approx(1.) and (axis.weights > 0).all()
    for k, name in enumerate(fit.specimen_names):
        expected = np.maximum(model.specimen_curves[name], floor[name]) / np.log(2.)
        np.testing.assert_allclose(fit.specimen[k], expected, rtol=1e-6)
    h2 = list(fit.specimen_names).index('h2')
    assert model.specimen_curves['h2'][1].max() < floor['h2']     # the absent gene sits on the floor ...
    np.testing.assert_allclose(fit.specimen[h2, 1], floor['h2'] / np.log(2.), rtol=1e-6)
    assert np.isfinite(fit.specimen).all()                         # ... so every log ratio stays finite
    design = Design.from_labels(c['specimen'], np.where(c['human'], 'human', 'mouse'), case_label='human',
                                reference_label='mouse')
    result = compare_groups(fit, axis, design, node_spacing=.05, params=CompareParams(prior_k=0.))
    table = result.gene_table.set_index('gene')
    assert table.selected['g0'] and table.shape_replicated['g0']
    # No prior: the pooled gap is the group-mean difference of the log2 curves, minus the typical gene.
    log2 = {n: fit.specimen[k] for k, n in enumerate(fit.specimen_names)}
    raw = (log2['h1'] + log2['h2']) / 2 - (log2['m1'] + log2['m2']) / 2
    np.testing.assert_allclose(result.contrasts.D, raw - np.median(raw, axis=0), rtol=1e-5, atol=1e-6)


def test_single_split_design_keeps_one_mixed_split():
    pytest.importorskip('pseudospace_reconstruction')
    from pseudospace_reconstruction.compare import Design

    design = Design.from_labels(['h1', 'h2', 'm1', 'm2'], ['human', 'human', 'mouse', 'mouse'], case_label='human',
                                reference_label='mouse')
    assert len(design.mixed_splits) == 2
    one = single_split_design(design, ['m2', 'h1'])
    assert len(one.mixed_splits) == 1 and set(one.mixed_splits[0].high) == {'h1', 'm2'}
    assert single_split_design(design, ['h2', 'm2']).null_names != one.null_names   # the other side names it too
    assert one.pairs == design.pairs and one.case_samples == design.case_samples
    with pytest.raises(ValueError, match='0 mixed splits'):
        single_split_design(design, ['h1', 'h2'])


# Notebook 67: T_spatial normalized by specimen shape variation.

def test_specimen_shape_design_adds_twelve_columns_for_the_split_and_six_for_a_relabeling(cohort):
    c = cohort
    knots = gam_internal_knots(c['position'])
    parts = partition_groups(c['specimen'], c['human'])
    for label, expected in (('species', 12), ('swap:h1+m1', 6), ('swap:h1+m2', 6)):
        group, nuisance = parts[label]
        designs, weights = specimen_shape_designs(c['position'], group, c['specimen'], knots, nuisance_shape=nuisance)
        reference, ref_weights = contrast_designs(c['position'], group, c['specimen'], knots, nuisance_shape=nuisance)
        np.testing.assert_array_equal(designs['full'], reference['full'])
        np.testing.assert_array_equal(weights, ref_weights)
        full, spec = (np.linalg.matrix_rank(designs[k]) for k in ('full', 'spec'))
        assert spec - full == expected and spec == designs['spec'].shape[1] - (0 if nuisance is None else 6)
        # M_spec is every specimen's own intercept + shape, whatever the partition.
        spline = make_gam_design(c['position'], knots)
        own = np.column_stack([(c['specimen'] == n)[:, None] * spline for n in SPECIMENS])
        assert np.linalg.matrix_rank(np.column_stack([designs['spec'], own])) == spec == own.shape[1]


def _shape_null(seed=5, n_genes=300, n_per=150, tau=.15, alpha=.05):
    """2 + 2 specimens; every gene has species offsets and its own random smooth shape in each specimen (sd ``tau``
    per basis coefficient), but no group shape. Gene 0 adds a human-only bump: a true species x position gene."""
    rng = np.random.default_rng(seed)
    specimen = np.repeat(SPECIMENS, n_per)
    position = rng.uniform(0, 1, len(specimen))
    human = np.isin(specimen, ['h1', 'h2'])
    library = rng.lognormal(np.log(3000), .3, len(specimen))
    basis = make_gam_design(position, gam_internal_knots(position))[:, 1:]
    rate = (np.log(rng.uniform(3e-4, 3e-2, n_genes))[None, :] + .5 * position[:, None]
            + .3 * human[:, None] * rng.normal(0, 1, n_genes)[None, :])
    for name in SPECIMENS:
        rows = specimen == name
        rate[rows] += basis[rows] @ rng.normal(0, tau, (basis.shape[1], n_genes)) + rng.normal(0, .1, n_genes)
    rate[:, 0] += 1.2 * human * np.sin(2 * np.pi * position)
    counts = rng.poisson(rng.gamma(1 / alpha, library[:, None] * np.exp(rate) * alpha)).astype(float)
    return dict(counts=counts, library=library, position=position, human=human, specimen=specimen,
                alpha=np.full(n_genes, alpha))


def test_specimen_normalized_f_is_calibrated_where_raw_t_spatial_is_inflated():
    _limma_available()
    c = _shape_null()
    knots = gam_internal_knots(c['position'])
    covariate = np.log(c['counts'].mean(axis=0))
    for label, (group, nuisance) in partition_groups(c['specimen'], c['human']).items():
        designs, weights = contrast_designs(c['position'], group, c['specimen'], knots, nuisance_shape=nuisance)
        numerator = nested_nb_lr(c['counts'], designs, weights, np.log(c['library']), c['alpha'],
                                 {'T_spatial': COMPARISONS['T_spatial']})['T_spatial']
        den = specimen_shape_lr(c['counts'], c['library'], c['position'], group, c['specimen'], c['alpha'],
                                nuisance_shape=nuisance)
        assert den['df2'] == (12 if nuisance is None else 6) and den['converged'].all() and not den['separated'].any()
        f = moderated_f(numerator, 6, den['LR_spec'], den['df2'], covariate=covariate)
        nulls = slice(1, None) if label == 'species' else slice(None)
        raw = chi2.sf(numerator[nulls], 6)
        assert np.mean(raw <= .05) > .4                                   # specimen shapes count as group shape
        assert .02 <= np.mean(f.p.to_numpy()[nulls] <= .05) <= .12        # judged against the specimens: ~5%
        assert np.isfinite(f.d0).all() and (f.d0 > 5).all() and np.isfinite(f.z).all()
        if label == 'species':
            assert f.p.iloc[0] < 1e-6 and f.z.iloc[0] == f.z.max()


def _fit_f_dist(x, df1):
    """limma's legacy fitFDist without a covariate, written out (Smyth 2004)."""
    from scipy.optimize import brentq
    from scipy.special import digamma, polygamma

    x = np.maximum(np.asarray(x, float), 0)
    x = np.maximum(x, 1e-5 * np.median(x))
    e = np.log(x) + np.log(df1 / 2) - digamma(df1 / 2)
    evar = e.var(ddof=1) - polygamma(1, df1 / 2)
    if evar <= 0:
        return np.inf, x.mean()
    df2 = 2 * brentq(lambda y: polygamma(1, y) - evar, 1e-8, 1e8, xtol=1e-14, rtol=1e-14)
    return df2, np.exp(e.mean() - np.log(df2 / 2) + digamma(df2 / 2))


def test_squeeze_var_matches_limmas_formula_and_handles_an_infinite_prior():
    _limma_available()
    rng = np.random.default_rng(4)
    variance = pd.Series(.7 * rng.chisquare(12, 50) / 12 * rng.lognormal(0, .4, 50), index=[f'g{i}' for i in range(50)])
    variance.iloc[3] = 0.                                     # a zero LR: offset away from zero for the prior only
    out = squeeze_var(variance, 12)
    d0, s0 = _fit_f_dist(variance, 12)
    assert list(out.columns) == ['var', 'var_prior', 'df_prior', 'var_post'] and out.index.equals(variance.index)
    np.testing.assert_allclose(out.df_prior, d0, rtol=1e-6)
    np.testing.assert_allclose(out.var_prior, s0, rtol=1e-6)
    np.testing.assert_allclose(out.var_post, (12 * variance + d0 * s0) / (12 + d0), rtol=1e-6)
    flat = squeeze_var(np.linspace(.999, 1.001, 30), 6)       # no variance heterogeneity: d0 = inf, var_post = prior
    assert np.isinf(flat.df_prior).all()
    np.testing.assert_allclose(flat.var_post, flat['var'].mean())
    covariate = np.linspace(0, 5, 400)                        # a planted trend of the prior variance
    trend = np.exp(-.5 * covariate) * rng.chisquare(12, 400) / 12
    fitted = squeeze_var(trend, 12, covariate=covariate)
    assert fitted.var_prior.iloc[0] > 5 * fitted.var_prior.iloc[-1]
    robust = squeeze_var(np.r_[trend, 50.], 12, covariate=np.r_[covariate, 2.5], robust=True)
    assert robust.df_prior.iloc[-1] < robust.df_prior.iloc[:-1].median()   # the outlier gets less shrinkage


def test_log_f_sf_is_exact_in_the_far_tail_and_moderated_f_handles_an_infinite_prior():
    from scipy.stats import f as f_dist
    values = np.array([0., .3, 1., 4., 20.])
    for df2 in (6., 14.5, 300.):
        np.testing.assert_allclose(log_f_sf(values, 6, df2), f_dist.logsf(values, 6, df2), rtol=1e-10, atol=1e-14)
    np.testing.assert_allclose(log_f_sf(values, 6, np.inf), chi2.logsf(6 * values, 6), rtol=1e-10, atol=1e-14)
    far = log_f_sf(np.array([1e30, 1e6]), 6, np.array([30., np.inf]))
    assert np.isfinite(far).all() and far[0] < -700 and far[1] < -1e6
    mpmath = pytest.importorskip('mpmath')
    mpmath.mp.dps = 40
    x = mpmath.mpf(30) / (30 + 6 * mpmath.mpf(10) ** 30)
    exact = [float(mpmath.log(mpmath.betainc(15, 3, 0, x, regularized=True))),
             float(mpmath.log(mpmath.gammainc(3, 3e6, mpmath.inf, regularized=True)))]
    np.testing.assert_allclose(far, exact, rtol=1e-10)
    near = log_f_sf(np.array([1e-6, 1e-9]), 6, np.array([17.5, np.inf]))     # p within 1e-20 of 1: still below 0
    exact = [float(mpmath.log(1 - mpmath.betainc(3, 8.75, 0, 6e-6 / (17.5 + 6e-6), regularized=True))),
             float(mpmath.log(1 - mpmath.gammainc(3, 0, 3e-9, regularized=True)))]
    np.testing.assert_allclose(near, exact, rtol=1e-10)
    assert (near < 0).all() and np.isfinite(-ndtri_exp(near)).all()
    _limma_available()
    flat = moderated_f(pd.Series([0., 3., 60., 3000.]), 6, np.linspace(5.994, 6.006, 4), 6)
    assert np.isinf(flat.d0).all()
    np.testing.assert_allclose(flat.p.iloc[1:], chi2.sf(flat.F.iloc[1:] * 6, 6), rtol=1e-10)
    assert flat.z.iloc[0] == flat.z.iloc[1:].min() and np.isfinite(flat.z).all()     # F = 0 ties at the bottom


def test_over_representation_reproduces_the_hypergeometric_test():
    pytest.importorskip('gseapy')
    from scipy.stats import hypergeom
    genes = [f'G{i:03d}' for i in range(200)]
    listed = genes[:20]
    sets = {'planted': genes[:15] + genes[100:110], 'random': genes[15:17] + genes[120:150], 'empty': genes[150:170]}
    out = over_representation(listed, sets, genes)
    assert list(out.columns) == ['pathway_id', 'n', 'overlap', 'p', 'fdr', 'call']
    table = out.set_index('pathway_id')
    assert table.overlap.to_dict() == {'planted': 15, 'random': 2, 'empty': 0}
    assert table.p['planted'] == pytest.approx(hypergeom.sf(14, 200, 25, 20), rel=1e-10)
    assert table.p['random'] == pytest.approx(hypergeom.sf(1, 200, 32, 20), rel=1e-10)
    assert table.p['empty'] == 1 and table.fdr['empty'] == 1 and table.call.to_dict() == {
        'planted': True, 'random': False, 'empty': False}
    nothing = over_representation([], sets, genes)
    assert (nothing.p == 1).all() and not nothing.call.any()
    with pytest.raises(ValueError, match='background'):
        over_representation(['not_a_gene'], sets, genes)
