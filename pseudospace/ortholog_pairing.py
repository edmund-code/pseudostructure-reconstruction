"""Alternative pairings of the shared ortholog space inside paralog families (notebook 64).

The paper's shared space keeps reciprocal-best HCOP pairs, so inside a multi-member orthogroup it
keeps whichever pair scored best on database support, not the paralog that carries the PT signal.
These helpers rebuild the pair table two ways:

* ``swap_pairs``: a family's main PT paralog replaces the least-expressed kept member when the
  reciprocal-best pick missed it;
* ``drop_pairs``: every kept pair that touches a multi-member family is removed;

and recompute notebook 38's probe-balance flag for any pair table.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

FAMILY_CLASSES = ('1:many', 'many:1', 'many:many')


def _main_member(table, component, candidates=None):
    """Highest-expressed PT-expressed, probe-clean member of ``component`` (None if there is none).

    ``table`` is indexed by symbol with columns component, pt_expressed, probe_clean and log2_cp10k.
    With ``candidates`` the choice is restricted to those symbols and the PT and probe filters are
    not applied (fallback to the kept members).
    """
    if candidates is None:
        members = table[table.component.eq(component) & table.pt_expressed.astype(bool) & table.probe_clean.astype(bool)]
    else:
        members = table.reindex([c for c in candidates if c in table.index])
    if members.empty:
        return None
    best = members.log2_cp10k.astype(float)
    return str(best.sort_values(ascending=False, kind='stable').index[0])


def family_of_pairs(pairs, mouse, human):
    """Component id of each kept pair (from its mouse member, else its human member); -1 if neither
    member lies in a multi-member family."""
    def component(table, symbol):
        if symbol in table.index and table.loc[symbol, 'ortholog_class'] in FAMILY_CLASSES:
            return int(table.loc[symbol, 'component'])
        return -1
    comp_m = np.array([component(mouse, m) for m in pairs.mouse_symbol])
    comp_h = np.array([component(human, h) for h in pairs.human_symbol])
    if ((comp_m >= 0) & (comp_h >= 0) & (comp_m != comp_h)).any():
        raise ValueError('A kept pair joins two different families.')
    return np.where(comp_m >= 0, comp_m, comp_h)


def swap_pairs(pairs, mouse, human):
    """Put each family's main PT paralog into the shared space when the kept pair(s) miss it.

    ``pairs`` has human_symbol, mouse_symbol and mapping_status. ``mouse`` and ``human`` are notebook
    62 gene tables (index = symbol; component, ortholog_class, pt_expressed, probe_clean,
    log2_cp10k). For each family and each species side separately: if the main paralog (highest
    pooled CP10k among PT-expressed, probe-clean members) is not already a kept member, it replaces
    the kept member of that side with the lowest pooled CP10k. Families whose main paralog is
    already kept are unchanged, so true sub-pairs (e.g. ACTB-Actb inside an actin family) are kept and
    the number of pairs does not change. Returns (new pair table, log of changed pairs).
    """
    pairs = pairs.reset_index(drop=True).copy()
    family = family_of_pairs(pairs, mouse, human)
    original = pairs.copy()
    for comp in sorted(set(family[family >= 0])):
        rows = np.flatnonzero(family == comp)
        for column, table in (('human_symbol', human), ('mouse_symbol', mouse)):
            main = _main_member(table, comp)
            kept = pairs.loc[rows, column]
            if main is None or main in set(kept):
                continue
            level = kept.map(table.log2_cp10k).astype(float).fillna(-np.inf)
            target = level.sort_values(kind='stable').index[0]
            pairs.loc[target, column] = main
            pairs.loc[target, 'mapping_status'] = 'family_main_paralog_swap'
    changed = (pairs.human_symbol != original.human_symbol) | (pairs.mouse_symbol != original.mouse_symbol)
    log = pd.DataFrame({'component': family[changed.to_numpy()],
                        'kept pair': (original.human_symbol + '–' + original.mouse_symbol)[changed].to_numpy(),
                        'new pair': (pairs.human_symbol + '–' + pairs.mouse_symbol)[changed].to_numpy(),
                        'human changed': (pairs.human_symbol != original.human_symbol)[changed].to_numpy(),
                        'mouse changed': (pairs.mouse_symbol != original.mouse_symbol)[changed].to_numpy()})
    for column in ('human_symbol', 'mouse_symbol'):
        if pairs[column].str.upper().duplicated().any():
            raise ValueError(f'Swapped map is not bijective in {column}.')
    return pairs.sort_values('human_symbol', key=lambda s: s.str.upper()).reset_index(drop=True), log


def drop_pairs(pairs, mouse, human):
    """Remove every kept pair whose human or mouse member lies in a multi-member family."""
    pairs = pairs.reset_index(drop=True)
    family = family_of_pairs(pairs, mouse, human)
    return pairs[family < 0].reset_index(drop=True), pairs[family >= 0].reset_index(drop=True)


def probe_balance(pairs, mouse_probe_counts, human_probe_counts, mouse_ids, human_ids):
    """Notebook 38's probe table for a pair table.

    ``*_probe_counts`` map Ensembl gene id -> included probes; ``*_ids`` map symbol -> Ensembl id.
    Returns mouse_symbol, human_symbol, mouse_probes, human_probes and probe_balance
    ('balanced' when equal, else 'mouse_fewer' or 'human_fewer').
    """
    table = pairs[['mouse_symbol', 'human_symbol']].copy()
    table['mouse_probes'] = table.mouse_symbol.map(mouse_ids).map(mouse_probe_counts).fillna(0).astype(int)
    table['human_probes'] = table.human_symbol.map(human_ids).map(human_probe_counts).fillna(0).astype(int)
    table['probe_balance'] = np.select([table.mouse_probes.eq(table.human_probes), table.mouse_probes.lt(table.human_probes)],
                                       ['balanced', 'mouse_fewer'], 'human_fewer')
    return table
