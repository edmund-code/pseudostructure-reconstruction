"""Pathway gene-set handling: ortholog-aware membership, per-pathway coverage, redundancy.

The notebooks matched pathway symbols to the expression universe by uppercasing both sides. That
is not ortholog mapping: a library written in human symbols loses every member whose mouse ortholog
has a different symbol, and the failure is invisible because the library-wide match rate cannot
distinguish

1. members never present in the library,
2. members the accepted ortholog table cannot map,
3. members mapped but absent from the assayed matrix,
4. members assayed but dropped by the expression filter.

Everything here is pure bookkeeping on symbols, so it is testable without expression data. The
caller supplies the expression universe (``adata.var_names``) and, for the cross-species workflow,
the accepted ortholog table already used to build that universe.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = [
    'resolve_member_symbols',
    'build_pathway_membership',
    'summarize_pathway_redundancy',
    'member_gene_evidence',
]


def _as_upper_set(values) -> set[str]:
    return {str(value).upper() for value in values}


def resolve_member_symbols(members, universe, ortholog_map=None,
                           source_column='human_symbol', target_column='mouse_symbol'):
    """Resolve one gene set onto the expression universe, by orthology when a table is given.

    ``universe`` is the matrix ``var_names`` (mouse-symbol space in the cross-species workflow).
    With an ``ortholog_map`` every member symbol is first looked up in the source-species column
    and replaced by its accepted target-species symbol; members that are already in the universe
    are kept as they are, so a mouse-symbol library still matches exactly.

    Returns ``(mapped, audit)``: ``mapped`` is the sorted list of universe symbols the set reaches,
    ``audit`` counts each stage so a caller can report where members were lost.
    """
    # Keep the caller's own symbol spelling: it is what indexes the matrix.
    universe_lookup = {str(gene).upper(): str(gene) for gene in universe}
    requested = _as_upper_set(members)

    translated: dict[str, str] = {}
    if ortholog_map is not None and not ortholog_map.empty:
        frame = ortholog_map.copy()
        frame.columns = [str(column) for column in frame.columns]
        if source_column in frame.columns and target_column in frame.columns:
            for source, target in zip(frame[source_column].astype(str), frame[target_column].astype(str)):
                translated.setdefault(source.upper(), target)

    resolved: set[str] = set()
    n_via_ortholog = 0
    for symbol in requested:
        if symbol in universe_lookup:                      # already in target-species form
            resolved.add(symbol)
            continue
        target = translated.get(symbol)
        if target is not None and target.upper() in universe_lookup:
            resolved.add(target.upper())
            n_via_ortholog += 1

    mapped = sorted(universe_lookup[symbol] for symbol in resolved)
    audit = {
        'n_requested': len(requested),
        'n_with_ortholog': n_via_ortholog,
        'n_assayed': len(mapped),
        'n_unmapped': len(requested) - len(mapped),
        'n_unmapped_assayed_out': len(requested) - len(mapped) - n_via_ortholog
        if ortholog_map is not None else 0,
    }
    return mapped, audit


def build_pathway_membership(library_sets, universe, *, ortholog_map=None, library_name=None,
                             min_genes=10, max_genes=None, tested=None,
                             source_column='human_symbol', target_column='mouse_symbol'):
    """Membership and coverage for every pathway in one library.

    ``max_genes`` is reported, not silently applied: a pathway above the cap is kept in the table
    with ``exclusion_reason`` set, so the caller can state which biology a cap removes instead of
    dropping it. ``tested`` (genes surviving the expression filter) adds the fourth coverage stage.

    Returns a DataFrame with one row per pathway and the coverage columns
    ``n_requested`` / ``n_with_ortholog`` / ``n_assayed`` / ``n_tested`` / ``n_genes_present``,
    plus ``genes_present`` (members surviving the expression filter, or every assayed member when no
    ``tested`` set is given), ``genes_assayed``, ``retained`` and ``exclusion_reason``. Retention is
    decided on the *assayed* count, so a pathway is still reported when the filter removes members.
    """
    tested_upper = _as_upper_set(tested) if tested is not None else None
    rows = []
    for pathway, members in library_sets.items():
        if isinstance(members, dict):
            members = members.get('genes', [])
        mapped, audit = resolve_member_symbols(
            members, universe, ortholog_map=ortholog_map,
            source_column=source_column, target_column=target_column,
        )
        tested_members = ([gene for gene in mapped if gene.upper() in tested_upper]
                          if tested_upper is not None else list(mapped))
        n_tested = len(tested_members) if tested_upper is not None else np.nan
        if audit['n_assayed'] < int(min_genes):
            retained, reason = False, f'fewer than {int(min_genes)} assayed members'
        elif max_genes is not None and audit['n_assayed'] > int(max_genes):
            retained, reason = False, f'more than {int(max_genes)} assayed members'
        else:
            retained, reason = True, ''
        rows.append({
            'library': library_name,
            'pathway': str(pathway),
            'n_requested': audit['n_requested'],
            'n_with_ortholog': audit['n_with_ortholog'],
            'n_assayed': audit['n_assayed'],
            'n_tested': n_tested,
            # `genes_present` is what a downstream model may actually use: the members that survive
            # the expression filter. `genes_assayed` keeps the wider assayed list for coverage
            # reporting, so a caller can no longer index a tested-gene table with an untested symbol.
            'n_genes_present': len(tested_members),
            'genes_present': tested_members,
            'genes_assayed': mapped,
            'retained': retained,
            'exclusion_reason': reason,
        })
    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    return frame.sort_values(['retained', 'n_assayed'], ascending=[False, False]).reset_index(drop=True)


def summarize_pathway_redundancy(membership, *, gene_column='genes_present',
                                 overlap_threshold=0.6, label_columns=('library', 'pathway'),
                                 only_retained=True):
    """Group pathways that measure largely the same genes.

    Overlap is ``|A & B| / min(|A|, |B|)``, so a small set contained in a large one is redundant
    even though a Jaccard index would be small. Pathways above the threshold are joined by
    connected components; the returned tables are

    * ``pairs``: one row per overlapping pair with its shared genes and overlap value,
    * ``groups``: one row per pathway with ``redundancy_group`` and ``group_size``.
    """
    frame = membership.copy()
    if only_retained and 'retained' in frame.columns:
        frame = frame[frame['retained'].astype(bool)]
    frame = frame.reset_index(drop=True)
    labels = [f'{row[label_columns[0]]}: {row[label_columns[1]]}' if label_columns[0] in frame.columns
              else str(row[label_columns[1]]) for _, row in frame.iterrows()]
    gene_sets = [set(map(str, genes)) for genes in frame[gene_column]]

    pairs = []
    for i in range(len(frame)):
        for j in range(i + 1, len(frame)):
            smaller = min(len(gene_sets[i]), len(gene_sets[j]))
            if smaller == 0:
                continue
            shared = sorted(gene_sets[i] & gene_sets[j])
            overlap = len(shared) / smaller
            if overlap >= overlap_threshold:
                pairs.append({
                    'pathway_a': labels[i], 'pathway_b': labels[j],
                    'n_shared': len(shared), 'n_smaller': smaller,
                    'overlap': overlap, 'shared_genes': '; '.join(shared),
                })

    parent = list(range(len(frame)))

    def _find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def _union(a, b):
        ra, rb = _find(a), _find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)

    index_by_label = {label: index for index, label in enumerate(labels)}
    for pair in pairs:
        _union(index_by_label[pair['pathway_a']], index_by_label[pair['pathway_b']])

    roots = {}
    group_ids = []
    for index in range(len(frame)):
        root = _find(index)
        if root not in roots:
            roots[root] = len(roots) + 1
        group_ids.append(roots[root])
    sizes = pd.Series(group_ids).value_counts().to_dict()
    groups = frame[list(label_columns)].copy() if all(c in frame.columns for c in label_columns) else pd.DataFrame()
    groups['redundancy_group'] = group_ids
    groups['group_size'] = [sizes[group] for group in group_ids]
    if 'n_assayed' in frame.columns:
        groups['n_assayed'] = frame['n_assayed'].to_numpy()
    return pd.DataFrame(pairs), groups


def member_gene_evidence(members, gene_effects, *, effect_column,
                         present_column='gene', tolerance=0.0):
    """Summarise how the member genes of one pathway behave relative to its aggregate trend.

    ``gene_effects`` is a per-gene table (``present_column`` names the gene column,
    ``effect_column`` the signed effect), e.g. the fitted level effect of each gene. Returns a
    Series with the member count, the count and fraction of members agreeing in sign with the
    member mean, and the aggregate effect with the strongest contributor removed -- the check for a
    pathway claim resting on one gene.
    """
    frame = gene_effects[[present_column, effect_column]].dropna()
    frame = frame.assign(_gene=frame[present_column].astype(str).str.upper())
    wanted = [str(gene).upper() for gene in members]
    subset = frame[frame['_gene'].isin(wanted)].copy()
    subset['_abs'] = subset[effect_column].abs()
    subset = subset.sort_values('_abs', ascending=False)
    aggregate = float(subset[effect_column].mean()) if len(subset) else np.nan
    if len(subset):
        sign = np.sign(aggregate) if aggregate != 0 else 0.0
        agree = int((np.sign(subset[effect_column].to_numpy()) == sign).sum()) if sign else 0
    else:
        agree = 0
    without = float(subset[effect_column].iloc[1:].mean()) if len(subset) > 1 else np.nan
    return pd.Series({
        'n_members_present': int(len(subset)),
        'member_effect_mean': aggregate,
        'n_members_agreeing': agree,
        'fraction_members_agreeing': (agree / len(subset)) if len(subset) else np.nan,
        'strongest_member': (str(subset[present_column].iloc[0]) if len(subset) else ''),
        'strongest_member_effect': (float(subset[effect_column].iloc[0]) if len(subset) else np.nan),
        'member_effect_mean_without_strongest': without,
        'sign_flips_without_strongest': bool(
            len(subset) > 1 and np.sign(without) != np.sign(aggregate)
        ),
        'tolerance': float(tolerance),
    })
