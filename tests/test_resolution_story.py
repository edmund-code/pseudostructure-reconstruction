"""Synthetic checks for saved pathway trajectory summaries."""
import pytest

np = pytest.importorskip('numpy')
pd = pytest.importorskip('pandas')
pytest.importorskip('scipy')
pytest.importorskip('patsy')
pytest.importorskip('statsmodels')

from pseudospace.pathway_remodeling import fit_nested_trajectories
from pseudospace.resolution_story import aggregate_saved_trajectories


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
