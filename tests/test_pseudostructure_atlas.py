"""Synthetic end-to-end checks for the prototype pseudostructure wrapper."""
import numpy as np
import pytest
pytest.importorskip('scipy')
pytest.importorskip('sklearn')
from scipy.stats import spearmanr
from scipy import sparse

from pseudospace.pseudostructure_atlas import (
    fit_pseudostructure_atlas, project_pseudostructure_atlas,
    summarize_panel_sensitivity,
)


def _synthetic(seed=48, n_per_specimen=150):
    rng = np.random.default_rng(seed)
    genes = np.array([f'g{i}' for i in range(10)])
    specimens, z_values, counts = [], [], []
    library = np.full(2*n_per_specimen, 2500, dtype=int)
    for specimen in ('A', 'B'):
        z = np.linspace(.01, .99, n_per_specimen)
        p = np.column_stack([
            .10 + .02*z, .12 - .02*z,
            .08 + .03*np.sin(np.pi*z), .10 - .03*np.sin(np.pi*z),
            np.full(len(z), .08), np.full(len(z), .08),
            np.full(len(z), .08), np.full(len(z), .08),
            .15 + .05*z, .13 - .05*z,
        ])
        counts.extend([rng.multinomial(int(lib), prob) for lib, prob in zip(library[:n_per_specimen], p)])
        specimens.extend([specimen]*n_per_specimen)
        z_values.extend(z)
    counts = np.asarray(counts, dtype=int)
    z_values = np.asarray(z_values)
    anatomy = np.where(z_values < .25, 0, np.where(z_values < .75, 1, 2))
    return counts, np.asarray(library), np.asarray(specimens), anatomy, genes, z_values


def _fit():
    counts, library, specimens, anatomy, genes, z = _synthetic()
    order = np.random.default_rng(9).permutation(len(counts))
    counts,library,specimens,anatomy,z = [a[order] for a in [counts,library,specimens,anatomy,z]]
    reference = fit_pseudostructure_atlas(
        sparse.csr_matrix(counts), library, specimens, anatomy, genes,
        np.array([True]*8 + [False]*2), n_components=5, seed=7,
        grid=np.linspace(0, 1, 31), max_iter=60, tol=1e-5,
        rate_bandwidth=.16, min_effective=5.)
    return reference, counts, library, specimens, anatomy, genes, z


def test_fit_projects_heldout_profiles_with_both_orders_and_frozen_batching():
    reference, counts, library, specimens, anatomy, genes, z = _fit()
    # Independent count draws, never profiles used to fit the reference.
    query, query_library, query_specimens, _, query_genes, query_z = _synthetic(seed=149,n_per_specimen=100)
    np.testing.assert_array_equal(genes,query_genes)
    query_specimens = np.repeat('unseen',len(query))
    ids = np.array([f'q{i}' for i in range(len(query))])
    full = project_pseudostructure_atlas(reference, query, query_library,
                                         genes, ids, specimens=query_specimens)
    assert spearmanr(full['gaussian']['z_mean'], query_z).statistic > .75
    assert spearmanr(full['count']['z_mean'], query_z).statistic > .75
    assert np.isfinite(full['coordinate_disagreement']).all()
    np.testing.assert_allclose(full['gaussian']['posterior'].sum(axis=1),1.)
    np.testing.assert_allclose(full['count']['posterior'].sum(axis=1),1.)
    first = project_pseudostructure_atlas(reference, sparse.csr_matrix(query[:20]), query_library[:20],
                                          genes, ids[:20], specimens=query_specimens[:20])
    second = project_pseudostructure_atlas(reference, query[20:], query_library[20:],
                                           genes, ids[20:], specimens=query_specimens[20:])
    np.testing.assert_allclose(np.r_[first['gaussian']['z_mean'], second['gaussian']['z_mean']],
                               full['gaussian']['z_mean'])
    np.testing.assert_allclose(np.r_[first['count']['z_mean'], second['count']['z_mean']],
                               full['count']['z_mean'])


def test_masked_count_column_swap_and_gene_contract_are_frozen():
    reference, counts, library, specimens, anatomy, genes, _ = _fit()
    query, exposure, query_specimens = counts[:12].copy(), library[:12], specimens[:12]
    ids = np.array([f's{i}' for i in range(len(query))])
    baseline = project_pseudostructure_atlas(reference, query, exposure,
                                             genes, ids, specimens=query_specimens)
    swapped = query.copy()
    swapped[:, [8, 9]] = swapped[:, [9, 8]]
    changed = project_pseudostructure_atlas(reference, swapped, exposure,
                                            genes, ids, specimens=query_specimens)
    np.testing.assert_allclose(changed['gaussian']['z_mean'], baseline['gaussian']['z_mean'])
    np.testing.assert_allclose(changed['count']['z_mean'], baseline['count']['z_mean'])
    with pytest.raises(ValueError, match='exactly match'):
        project_pseudostructure_atlas(reference, query, exposure, genes[::-1], ids,
                                      specimens=query_specimens)
    with pytest.raises(ValueError, match='exactly match'):
        project_pseudostructure_atlas(reference, query[:, :-1], exposure, genes[:-1], ids,
                                      specimens=query_specimens)
    with pytest.raises(ValueError, match='unique'):
        project_pseudostructure_atlas(reference, query, exposure, genes,
                                      np.repeat('duplicate', len(query)),
                                      specimens=query_specimens)


def test_panel_sensitivity_requires_identical_ids_and_reports_ranges_not_pooling():
    reference, counts, library, specimens, anatomy, genes, _ = _fit()
    query = counts[:10]
    ids = np.array([f'id{i}' for i in range(len(query))])
    projection = project_pseudostructure_atlas(reference, query, library[:10], genes, ids)
    shifted = {'structure_ids': ids.copy(), 'gaussian': dict(projection['gaussian']),
               'count': dict(projection['count'])}
    shifted['gaussian']['z_mean'] = np.clip(projection['gaussian']['z_mean'] + .04, 0, 1)
    shifted['count']['z_mean'] = np.clip(projection['count']['z_mean'] - .03, 0, 1)
    summary = summarize_panel_sensitivity([projection, shifted])
    np.testing.assert_allclose(summary['gaussian']['z_median'],
                               np.median(np.vstack([projection['gaussian']['z_mean'], shifted['gaussian']['z_mean']]), axis=0))
    assert np.all(summary['gaussian']['z_min'] <= summary['gaussian']['z_max'])
    assert summary['n_panels'] == 2
    with pytest.raises(ValueError, match='identical structure ID order'):
        summarize_panel_sensitivity([projection, {**shifted, 'structure_ids': ids[::-1]}])


def test_fit_rejects_nonboolean_mask_and_insufficient_endpoint_support():
    counts, library, specimens, anatomy, genes, _ = _synthetic(n_per_specimen=30)
    with pytest.raises(ValueError, match='boolean vector'):
        fit_pseudostructure_atlas(counts, library, specimens, anatomy, genes,
                                  np.array([1]*8+[0]*2, dtype=int))
    with pytest.raises(ValueError, match='casefolding'):
        fit_pseudostructure_atlas(counts, library, specimens, anatomy,
                                  np.array(['Gene', 'gene'] + list(genes[2:])),
                                  np.array([True]*8+[False]*2))
    bad_anatomy = np.ones_like(anatomy)
    with pytest.raises(ValueError, match='endpoint labels'):
        fit_pseudostructure_atlas(counts, library, specimens, bad_anatomy,
                                  genes, np.array([True]*8+[False]*2))


def test_panel_summary_rejects_nonfinite_positions():
    good={'structure_ids':np.array(['a','b']),
          'gaussian':{'z_mean':np.array([.2,.8]),'entropy':np.array([1.,1.])},
          'count':{'z_mean':np.array([.3,.7]),'entropy':np.array([1.,1.])}}
    bad={**good,'count':{'z_mean':np.array([np.nan,.8]),'entropy':np.array([1.,1.])}}
    with pytest.raises(ValueError,match='finite and valid'):
        summarize_panel_sensitivity([good,bad])
