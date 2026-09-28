"""Synthetic checks for the appended interpretation stage; no private inputs."""
import numpy as np
import pytest

pd = pytest.importorskip('pandas')
pytest.importorskip('scipy')
from pseudospace.pathway_programs import (
    active_programs, overlap_sensitivity, program_figures, program_gene_evidence, representative_terms,
)


def test_overlap_matches_pairwise_auc_with_ties_and_exposes_shared_gene_dependence():
    stats = pd.Series([8., 2., 2., 0.], index=list('abcd'))
    sets = {'target': ['a', 'b'], 'other1': ['a'], 'other2': ['a', 'd']}
    result, frequency = overlap_sensitivity(stats, sets)
    target = result.set_index('pathway_id').loc['target']
    # a wins both comparisons; b ties c and beats d. Annotated a has less weight.
    assert target.original_effect == pytest.approx((1 + .75) / 2 - .5)
    assert target.overlap_weighted_effect == pytest.approx(np.average([1, .75], weights=[1/np.sqrt(3), 1]) - .5)
    assert target.overlap_effect_change < 0
    assert frequency.set_index('gene').loc['a', 'annotation_frequency'] == 3
    unchanged, _ = overlap_sensitivity(stats, {'target': ['a', 'b', 'b']})
    assert unchanged.overlap_effect_change.iloc[0] == 0  # duplicate annotation counts once
    with pytest.raises(ValueError, match='universe'):
        overlap_sensitivity(stats, {'bad': ['missing']})


def test_active_genes_drive_grouping_and_group_ids_are_order_invariant():
    sets = {'A': list('abxy'), 'B': list('abuv'), 'C': list('abxy')}
    members = pd.DataFrame([(p, g, g in active) for p, active in
        [('A', 'ab'), ('B', 'ab'), ('C', 'xy')] for g in sets[p]],
        columns=['pathway_id', 'gene', 'leading_edge'])
    members['broad_supporter'] = False
    members['spatial_driver'] = False
    curves = pd.DataFrame([(p, i, d, s) for p in sets for i, d, s in
        [(0, 0., -.5), (1, .2, 0.), (2, .4, .5)]],
        columns=['pathway_id', 'position', 'divergence', 'direction'])
    groups, pairs, active = active_programs(list(sets), sets, members, curves, distance_cut=.4)
    label = groups.set_index('pathway_id').paper_program
    assert label['A'] == label['B'] != label['C']
    ac = pairs[(pairs.source == 'A') & (pairs.target == 'C')].iloc[0]
    assert ac.full_jaccard == 1 and ac.active_jaccard == 0 and ac.similarity == pytest.approx(.5)
    reverse, _, _ = active_programs(list(reversed(sets)), sets, members, curves, distance_cut=.4)
    pd.testing.assert_frame_equal(groups, reverse)
    empty, empty_edges, _ = active_programs([], sets, members, curves)
    assert empty.empty and empty_edges.empty
    single, _, _ = active_programs(['A'], sets, members, curves)
    assert len(single) == 1
    bad = curves.copy()
    bad['position'] = bad.position.astype(float)
    bad.loc[bad.pathway_id.eq('B'), 'position'] += .1
    with pytest.raises(ValueError, match='coordinate grid'):
        active_programs(list(sets), sets, members, bad)

    atlas = pd.DataFrame({'pathway_id': list(sets), 'effect_T_spatial': [.4, .2, .1],
        'tested_fraction': [1., 1., 1.], 'n_retained': [1, 4, 4], 'n_planned': [7, 4, 4],
        'correlation_support': ['Correlation-supported'] * 3})
    terms = representative_terms(groups, atlas, active).set_index('pathway_id')
    assert not terms.loc['A', 'paper_priority']  # one successful run cannot look like 100% stability
    assert terms.loc['B', 'paper_priority']
    assert terms.representative_order.gt(0).sum() == 2  # duplicate active sets get one representative


def test_bidirectional_genes_do_not_cancel_and_singletons_are_not_called_private():
    genes = list('abcd')
    grid = np.array([0., .5, 1.])
    z = np.array([[3, 3, -3], [4, 4, -4], [-3, 3, -3], [-4, 4, -4]])
    fit = {'z': z, 'delta': z * .1, 'human': 1 + z * .1, 'mouse': np.ones_like(z)}
    terms = pd.DataFrame({'pathway_id': ['A', 'B'], 'paper_program': ['P', 'P']})
    members = pd.DataFrame({'pathway_id': ['A'] * 4 + ['B'], 'gene': genes + ['a'],
        'leading_edge': True, 'broad_supporter': False, 'spatial_driver': False, 'bulk_logcpm_effect': 0.})
    evidence, direction = program_gene_evidence(terms, members, genes, grid, fit)
    assert direction.direction_description.tolist() == ['mixed', 'human-high dominated', 'mouse-high dominated']
    assert z[:, 0].mean() == 0  # cancellation would conceal the strong opposing signals
    assert evidence.set_index('gene').loc['a', 'gene_role'] == 'shared active driver'
    assert evidence.set_index('gene').loc['b', 'gene_role'] == 'term-specific active gene'
    assert direction.n_active_genes.eq(4).all()  # shared a is counted once
    single, _ = program_gene_evidence(terms.iloc[:1], members, genes, grid, fit)
    assert single.gene_role.eq('single-term program').all()
    empty, direction_empty = program_gene_evidence(terms.iloc[:0], members, genes, grid, fit)
    assert empty.empty and direction_empty.empty
    with pytest.raises(ValueError, match='align'):
        program_gene_evidence(terms, members, genes, grid, {**fit, 'z': z.T})


def test_program_packet_prioritizes_drivers_and_keeps_complete_gene_panel():
    plt = pytest.importorskip('matplotlib.pyplot')
    genes = [f'g{i:02d}' for i in range(16)]
    grid = np.array([0., .5, 1.])
    table = pd.DataFrame({'paper_program': 'P', 'gene': genes, 'gene_role': 'single-term program',
        'n_active_terms': 1, 'peak_position': np.linspace(0, 1, 16), 'peak_abs_z': np.arange(16)[::-1],
        'spatial_driver': [False] * 15 + [True], 'leading_edge': False, 'broad_supporter': False})
    terms = pd.DataFrame({'paper_program': ['P'], 'pathway_id': ['A'], 'representative_order': [1]})
    curves = pd.DataFrame({'pathway_id': ['A'] * 3, 'position': grid,
                           'divergence': [0., .1, .2], 'direction': [-.2, 0., .2]})
    z = np.tile(grid, (16, 1))
    fit = {'z': z, 'human': z, 'mouse': -z}
    figures = dict(program_figures('P', terms, table, curves, genes, grid, fit))
    try:
        top_labels = [label.get_text() for label in figures['top_genes'].axes[0].get_yticklabels()]
        full_labels = [label.get_text() for label in figures['active_genes_01'].axes[0].get_yticklabels()]
        assert len(top_labels) == 15 and 'g15' in top_labels and 'g14' not in top_labels
        assert len(full_labels) == 16 and 'g14' in full_labels
    finally:
        for figure in figures.values():
            plt.close(figure)
