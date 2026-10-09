import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scipy.stats import chi2

from pseudospace.nb_gam import nb_group_difference_curves, nb_species_trajectories
from pseudospace.pathway_pipelines import (common_grid, partition_groups, rank_first_gsea, rank_normal_scores,
                                           specimen_ranges, standard_gene_model)

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
