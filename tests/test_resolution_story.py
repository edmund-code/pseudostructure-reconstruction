"""Synthetic checks for saved pathway trajectory summaries."""
import pytest

np = pytest.importorskip('numpy')
pd = pytest.importorskip('pandas')
pytest.importorskip('scipy')
pytest.importorskip('patsy')
pytest.importorskip('statsmodels')

from pseudospace.pathway_remodeling import fit_nested_trajectories
from pseudospace.resolution_story import (
    aggregate_saved_trajectories,
    controlled_resolution_calls,
    pathway_method_membership,
    signed_gene_heatmap,
)


def _membership_fixture():
    return pd.DataFrame({
        'pathway_id': ['spatial', 'conventional', 'segment', 'level', 'none'],
        'pathway_name': ['Spatial', 'Bulk', 'Segment', 'Level', 'None'],
        'q_family_DESeq2': [.8, .01, .8, .8, .8],
        'DESeq2_hit': [False, True, False, False, False],
        'bulk_hit': [False, True, False, False, False],
        'cluster_q_family_PT-S1': [.8, .8, .01, .8, .8],
        'cluster_q_family_PT-S2': [.8] * 5,
        'cluster_q_family_PT-S3': [.8] * 5,
        'cluster_hit': [False, False, True, False, False],
        'common_effect_T_level': [0, 0, 0, 1, 0],
        'common_q_empirical_T_level': [.8, .8, .8, .01, .8],
        'level_hit': [False, False, False, True, False],
        'common_effect_T_spatial': [1, 0, 0, 0, 0],
        'common_q_empirical_T_spatial': [.01, .8, .8, .8, .8],
        'spatial_hit': [True, False, False, False, False],
        'common_effect_T_total': [0, 0, 0, 1, 0],
        'common_q_empirical_T_total': [.8, .8, .8, .01, .8],
        'total_hit': [False, False, False, True, False],
    })


def test_pathway_method_membership_uses_saved_calls_and_classifies_discoveries():
    source = _membership_fixture()
    got = pathway_method_membership(source)
    assert got.pathway_id.tolist() == source.pathway_id.tolist()
    assert got.pathway_name.tolist() == source.pathway_name.tolist()
    assert got['whole_PT_DESeq2_GSEA'].tolist() == [False, True, False, False, False]
    assert got['S1_GSEA'].tolist() == [False, False, True, False, False]
    assert got['T_spatial'].tolist() == [True, False, False, False, False]
    assert got['spatial_only'].tolist() == [True, False, False, False, False]
    assert got['conventional_only'].tolist() == [False, True, True, False, False]
    assert got['discovery_category'].tolist() == [
        'Spatial only', 'Conventional only', 'Conventional only', 'Other/shared', 'Other/shared']
    assert got.common_spatial_effect.iloc[0] == 1



def test_pathway_membership_retains_q_boundary_and_signed_negative_nes():
    source = _membership_fixture()
    source.loc[1, 'q_family_DESeq2'] = .05
    source['NES_DESeq2'] = -2.0  # Signed significance accepts either direction.
    source.loc[0, 'common_q_empirical_T_spatial'] = .05
    source.loc[3, 'common_effect_T_total'] = -1
    source.loc[3, 'total_hit'] = False
    got = pathway_method_membership(source)
    assert got.whole_PT_DESeq2_GSEA.iloc[1]
    assert got.T_spatial.iloc[0]
    assert not got.T_total.iloc[3]
    source.loc[0, 'common_effect_T_spatial'] = 0
    source.loc[0, 'spatial_hit'] = False
    assert not pathway_method_membership(source).T_spatial.iloc[0]

def test_pathway_method_membership_includes_optional_control_but_not_conventional_union():
    source = _membership_fixture()
    source['common_effect_T_whole_pt_deseq2_abs'] = [1, 0, 0, 0, 0]
    source['common_q_empirical_T_whole_pt_deseq2_abs'] = [.01, .8, .8, .8, .8]
    got = pathway_method_membership(source)
    assert got.whole_PT_matched_AUC.tolist() == [True, False, False, False, False]
    assert got.any_conventional.tolist() == [False, True, True, False, False]


@pytest.mark.parametrize('column,value', [
    ('DESeq2_hit', [False, False, False, False, False]),
    ('cluster_hit', [False, False, False, False, False]),
    ('spatial_hit', [False, False, False, False, False]),
    ('level_hit', [False, False, False, False, False]),
    ('bulk_hit', ['false'] * 5),
])
def test_pathway_method_membership_rejects_mismatched_or_invalid_saved_calls(column, value):
    with pytest.raises(ValueError):
        pathway_method_membership(_membership_fixture().assign(**{column: value}))


def test_pathway_method_membership_rejects_partial_optional_control_and_invalid_ids():
    source = _membership_fixture()
    source['common_effect_T_whole_pt_deseq2_abs'] = 1
    with pytest.raises(ValueError, match='together'):
        pathway_method_membership(source)
    with pytest.raises(ValueError, match='unique'):
        pathway_method_membership(_membership_fixture().assign(pathway_id=['p'] * 5))
    with pytest.raises(ValueError, match='present'):
        pathway_method_membership(_membership_fixture().assign(pathway_id=['p', None, 'a', 'b', 'c']))


def test_pathway_mean_curves_and_hc3_match_fit_of_mean_and_preserve_covariance():
    rng = np.random.default_rng(170)
    sizes = {'m1': 80, 'm2': 55, 'h1': 70, 'h2': 45}
    specimen = np.concatenate([np.repeat(name, count) for name, count in sizes.items()])
    human = np.isin(specimen, ['h1', 'h2']).astype(float)
    position = np.concatenate([np.linspace(.01, .99, count) for count in sizes.values()])
    shared_noise = rng.normal(0, .25, len(position))
    y = np.column_stack([
        2 + position + human + shared_noise + rng.normal(0, .005, len(position)),
        3 + .5 * position + .7 * human + shared_noise + rng.normal(0, .005, len(position)),
        1 + 1.5 * position + 1.2 * human + shared_noise + rng.normal(0, .005, len(position)),
    ])
    genes = ['a', 'b', 'c']
    grid = np.linspace(.04, .96, 13)
    fit = fit_nested_trajectories(y, position, human, specimen, grid, basis_df=3,
                                  return_residuals=True)
    fit['genes'] = np.asarray(genes)
    result = aggregate_saved_trajectories(
        fit, position, human, specimen, {'shared': genes})
    row = result.loc[result.pathway_id.eq('shared')]
    reference = fit_nested_trajectories(y.mean(axis=1, keepdims=True), position, human,
        specimen, grid, basis_df=3)
    np.testing.assert_allclose(row.mouse, reference['mouse'][0])
    np.testing.assert_allclose(row.human, reference['human'][0])
    np.testing.assert_allclose(row.delta, reference['delta'][0])
    np.testing.assert_allclose(row.se_delta, reference['se'][0])
    np.testing.assert_allclose(row.position, (grid - grid.min()) / np.ptp(grid))
    assert row.n_genes.eq(3).all()

    independent_gene_se = np.sqrt(np.mean(fit['se'] ** 2, axis=0) / len(genes))
    assert np.any(row.se_delta.to_numpy() > independent_gene_se)


def test_pathway_summary_rejects_unknown_or_duplicate_members():
    position = np.tile(np.linspace(.01, .99, 40), 4)
    human = np.repeat([0, 0, 1, 1], 40).astype(float)
    specimen = np.repeat(['m1', 'm2', 'h1', 'h2'], 40)
    fit = fit_nested_trajectories(np.ones((160, 2)), position, human, specimen,
        np.linspace(.05, .95, 7), basis_df=3, return_residuals=True)
    fit['genes'] = np.asarray(['a', 'b'])
    with pytest.raises(ValueError, match='unknown genes'):
        aggregate_saved_trajectories(fit, position, human, specimen, {'bad': ['c']})
    with pytest.raises(ValueError, match='unique'):
        aggregate_saved_trajectories(fit, position, human, specimen, {'bad': ['a', 'a']})


def test_controlled_resolution_calls_uses_positive_method_specific_q_and_checks_saved_calls():
    comparison = pd.DataFrame({
        'pathway_id': ['retained', 'added', 'discrete_only', 'negative', 'pooled_only'],
        'effect_T_discrete_total': [1, 0, 1, -1, 0],
        'q_method_T_discrete_total': [.01, .01, .01, .01, .9],
        'effect_T_continuous_total': [1, 1, 0, -1, 1],
        'q_method_T_continuous_total': [.02, .03, .01, .01, .9],
        # Pooled calls intentionally disagree: this helper must ignore them.
        'controlled_discrete_hit_pooled': [False] * 5,
        'controlled_continuous_hit_pooled': [False] * 5,
        'controlled_discrete_hit': [True, False, True, False, False],
        'controlled_continuous_hit': [True, True, False, False, False],
    })
    got = controlled_resolution_calls(comparison)
    assert got.pathway_id.tolist() == comparison.pathway_id.tolist()
    assert got.controlled_added.tolist() == [False, True, False, False, False]
    assert got.controlled_retained.tolist() == [True, False, False, False, False]
    assert 'controlled_added' not in comparison
    with pytest.raises(ValueError, match='disagrees'):
        controlled_resolution_calls(comparison.assign(controlled_continuous_hit=False))
    derived = controlled_resolution_calls(comparison.drop(
        columns=['controlled_discrete_hit', 'controlled_continuous_hit']))
    assert derived.controlled_discrete_hit.tolist() == [True, False, True, False, False]
    assert derived.controlled_continuous_hit.tolist() == [True, True, False, False, False]


@pytest.mark.parametrize('change, message', [
    ({'pathway_id': ['p', 'p']}, 'unique'),
    ({'pathway_id': ['p', None]}, 'present'),
    ({'effect_T_discrete_total': [1, np.nan]}, 'finite'),
    ({'q_method_T_continuous_total': [.1, 1.1]}, 'finite'),
])
def test_controlled_resolution_calls_rejects_invalid_rows(change, message):
    base = pd.DataFrame({
        'pathway_id': ['p1', 'p2'], 'effect_T_discrete_total': [1, 1],
        'q_method_T_discrete_total': [.01, .02],
        'effect_T_continuous_total': [1, 1], 'q_method_T_continuous_total': [.01, .02],
    })
    with pytest.raises(ValueError, match=message):
        controlled_resolution_calls(base.assign(**change))


def test_signed_gene_heatmap_prioritizes_spatial_statistic_then_orders_peak_position():
    fit = {
        'genes': np.array(['z_gene', 'b_gene', 'a_gene', 'low']),
        'grid': np.linspace(0, 1, 4),
        'T_spatial': np.array([3., 5., 5., 1.]),
        'z': np.array([[0, 4, 1, 0], [0, 0, -3, 0], [0, 3, 0, 0], [9, 0, 0, 0.]]),
    }
    genes, z = signed_gene_heatmap(fit, ['z_gene', 'b_gene', 'a_gene', 'low'], max_genes=3)
    assert genes.tolist() == ['a_gene', 'z_gene', 'b_gene']
    np.testing.assert_array_equal(z, fit['z'][[2, 0, 1]])
    with pytest.raises(ValueError, match='unique'):
        signed_gene_heatmap(fit, ['a_gene', 'a_gene'])
    with pytest.raises(ValueError, match='unknown'):
        signed_gene_heatmap(fit, ['missing'])


def test_signed_gene_heatmap_rejects_misaligned_or_nonfinite_fit():
    fit = {'genes': np.array(['a', 'b']), 'grid': np.array([0., 1.]),
           'T_spatial': np.array([1., 2.]), 'z': np.ones((2, 2))}
    with pytest.raises(ValueError, match='align'):
        signed_gene_heatmap({**fit, 'z': np.ones((1, 2))}, ['a'])
    with pytest.raises(ValueError, match='finite'):
        signed_gene_heatmap({**fit, 'T_spatial': np.array([1., np.nan])}, ['a'])
    with pytest.raises(ValueError, match='present'):
        signed_gene_heatmap({**fit, 'genes': np.array(['a', None], dtype=object)}, ['a'])
    with pytest.raises(ValueError, match='align'):
        signed_gene_heatmap({**fit, 'grid': np.array([0., 0.])}, ['a'])
    with pytest.raises(ValueError, match='align'):
        signed_gene_heatmap({**fit, 'grid': np.array([0.])}, ['a'])
