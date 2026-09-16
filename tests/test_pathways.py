"""Pathway bookkeeping: ortholog-aware membership, coverage stages, redundancy, member evidence."""
from __future__ import annotations

import pandas as pd
import pytest

from pseudospace.pathways import (
    build_pathway_membership,
    member_gene_evidence,
    resolve_member_symbols,
    summarize_pathway_redundancy,
)

UNIVERSE = ['Slc12a1', 'Umod', 'Kcnj1', 'Cldn16', 'Atp6v0d2', 'Rhcg', 'Aqp2']


def _ortholog_map():
    return pd.DataFrame({
        # letter-identical across species: reachable without the table (matched case-insensitively)
        'human_symbol': ['SLC12A1', 'UMOD', 'KCNJ1', 'ATP6V0D2', 'RHCG',
        # genuinely renamed orthologs: only the table gets these into the universe
                         'HUMAN_ONLY_NAME', 'ORPHAN'],
        'mouse_symbol': ['Slc12a1', 'Umod', 'Kcnj1', 'Atp6v0d2', 'Rhcg',
                         'Aqp2', 'NotAssayed'],
    })


def test_members_resolve_by_orthology_not_only_by_case():
    mapped, audit = resolve_member_symbols(
        ['SLC12A1', 'UMOD', 'HUMAN_ONLY_NAME', 'NEVER_IN_LIBRARY'],
        UNIVERSE, ortholog_map=_ortholog_map(),
    )
    # Slc12a1/Umod match case-insensitively; HUMAN_ONLY_NAME needs the ortholog table to become
    # Aqp2; NEVER_IN_LIBRARY is in no library at all. Returned symbols keep the universe spelling.
    assert mapped == ['Aqp2', 'Slc12a1', 'Umod']
    assert audit['n_requested'] == 4
    assert audit['n_with_ortholog'] == 1
    assert audit['n_assayed'] == 3
    assert audit['n_unmapped'] == 1


def test_members_whose_ortholog_is_absent_are_reported_not_silently_dropped():
    mapped, audit = resolve_member_symbols(['DIFFERENTNAME'], UNIVERSE, ortholog_map=_ortholog_map())
    assert mapped == []
    assert audit['n_with_ortholog'] == 0
    assert audit['n_unmapped'] == 1


def test_target_species_symbols_pass_through_unchanged():
    mapped, audit = resolve_member_symbols(['Slc12a1', 'Umod'], UNIVERSE)
    assert mapped == ['Slc12a1', 'Umod']
    assert audit['n_assayed'] == 2


def test_membership_reports_every_coverage_stage_and_keeps_excluded_pathways():
    library = {
        'small': ['SLC12A1', 'UMOD'],
        'kept': ['SLC12A1', 'UMOD', 'KCNJ1', 'ATP6V0D2', 'RHCG'],
        'too_big': ['SLC12A1', 'UMOD', 'KCNJ1', 'CLDN16', 'ATP6V0D2', 'RHCG', 'Aqp2'],
    }
    membership = build_pathway_membership(
        library, UNIVERSE, ortholog_map=_ortholog_map(), library_name='Test_2024',
        min_genes=3, max_genes=5, tested=['SLC12A1', 'UMOD', 'KCNJ1', 'ATP6V0D2'],
    )
    by_pathway = membership.set_index('pathway')
    assert by_pathway.loc['kept', 'retained']
    assert by_pathway.loc['kept', 'n_requested'] == 5
    assert by_pathway.loc['kept', 'n_with_ortholog'] == 0
    assert by_pathway.loc['kept', 'n_assayed'] == 5
    assert by_pathway.loc['kept', 'n_tested'] == 4
    assert not by_pathway.loc['small', 'retained']
    assert 'assayed members' in by_pathway.loc['small', 'exclusion_reason']
    # A cap must not delete the pathway silently: it stays with its reason recorded.
    assert not by_pathway.loc['too_big', 'retained']
    assert 'more than 5' in by_pathway.loc['too_big', 'exclusion_reason']


def test_membership_exposes_the_tested_subset_separately_from_the_assayed_set():
    """Downstream scoring indexes a tested-gene table, so it must not receive untested symbols."""
    library = {'mixed': ['SLC12A1', 'UMOD', 'KCNJ1', 'ATP6V0D2', 'RHCG']}
    membership = build_pathway_membership(
        library, UNIVERSE, library_name='Test_2024', min_genes=3,
        tested=['SLC12A1', 'UMOD', 'KCNJ1'],          # two assayed members fail the filter
    )
    row = membership.iloc[0]
    assert row['n_assayed'] == 5
    assert row['n_tested'] == 3
    assert row['genes_present'] == ['Kcnj1', 'Slc12a1', 'Umod']   # universe spelling
    assert row['n_genes_present'] == 3
    assert row['genes_assayed'] == ['Atp6v0d2', 'Kcnj1', 'Rhcg', 'Slc12a1', 'Umod']
    # without a tested set the assayed members are the present ones (backwards compatible)
    without_filter = build_pathway_membership(library, UNIVERSE, library_name='Test_2024', min_genes=3)
    assert without_filter.iloc[0]['genes_present'] == ['Atp6v0d2', 'Kcnj1', 'Rhcg',
                                                       'Slc12a1', 'Umod']


def test_redundancy_groups_pathways_sharing_almost_all_members():
    membership = pd.DataFrame({
        'library': ['L'] * 3,
        'pathway': ['ros_rns_phagocytes', 'insulin_receptor_recycling', 'unrelated'],
        'genes_present': [
            [f'Atp6v0d2_{i}' for i in range(21)],
            [f'Atp6v0d2_{i}' for i in range(20)] + ['extra'],
            ['Aqp2', 'Rhcg', 'Umod', 'Kcnj1'],
        ],
        'n_assayed': [21, 21, 4],
        'retained': [True, True, True],
    })
    pairs, groups = summarize_pathway_redundancy(membership, overlap_threshold=0.6)
    assert len(pairs) == 1
    assert pairs.iloc[0]['n_shared'] == 20
    assert pairs.iloc[0]['overlap'] == pytest.approx(20 / 21)
    grouped = groups.set_index(['library', 'pathway'])['redundancy_group']
    assert grouped.loc[('L', 'ros_rns_phagocytes')] == grouped.loc[('L', 'insulin_receptor_recycling')]
    assert grouped.loc[('L', 'unrelated')] != grouped.loc[('L', 'ros_rns_phagocytes')]
    assert groups.set_index(['library', 'pathway']).loc[('L', 'unrelated'), 'group_size'] == 1


def test_member_gene_evidence_flags_a_one_gene_driver():
    effects = pd.DataFrame({
        'gene': ['Atp6v0d2', 'Rhcg', 'Aqp2', 'Umod'],
        'level_effect_human_minus_mouse': [-3.0, -0.2, -0.1, -0.1],
    })
    evidence = member_gene_evidence(
        ['Atp6v0d2', 'Rhcg', 'Aqp2'], effects, effect_column='level_effect_human_minus_mouse',
    )
    assert evidence['n_members_present'] == 3
    assert evidence['strongest_member'] == 'Atp6v0d2'
    assert evidence['n_members_agreeing'] == 3
    assert evidence['fraction_members_agreeing'] == pytest.approx(1.0)
    # Without the strongest member the aggregate collapses towards zero.
    assert abs(evidence['member_effect_mean_without_strongest']) < abs(evidence['member_effect_mean'])


def test_member_gene_evidence_detects_opposing_members():
    effects = pd.DataFrame({
        'gene': ['A', 'B', 'C'],
        'effect': [2.0, -2.0, 0.5],
    })
    evidence = member_gene_evidence(['A', 'B', 'C'], effects, effect_column='effect')
    assert evidence['n_members_present'] == 3
    # the mean is +0.167, so the two positive members agree and the opposing one does not
    assert evidence['n_members_agreeing'] == 2
    assert evidence['fraction_members_agreeing'] == pytest.approx(2 / 3)
