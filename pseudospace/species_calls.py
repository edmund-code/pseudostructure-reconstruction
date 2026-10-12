"""Descriptions of a list of positional human vs mouse pathway calls (notebook 77).

The rules come unchanged from notebooks 37 (E3), 38 (the "ours" arm), 43 §7 (categories), 47 (e) (peaks in segment
units) and 50 (overlap-preserving replication, conventional-only sets). They live here so that one code path serves
the old inputs (notebook 77's code-equivalence check) and the package inputs.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import norm, spearmanr
from statsmodels.stats.multitest import multipletests

from .continuous_descriptors import segment_units
from .pathway_decoys import welch_t
from .pathway_remodeling import local_pathway_curves, matching_strata, spatial_descriptions
from .replication_rules import arm_passes
from .set_collection_null import JointReplicationZ, collection_null, label_permutation_p, upper_p

CONTRASTS = ('S2-S1', 'S3c-early')
ZONATED_CALLS = ('up', 'down')
LED = ('mouse-led', 'human-led', 'shared', 'mixed')


# --- Categories (notebook 43 §7) --------------------------------------------------------------------------------

def member_category(row, contrasts=CONTRASTS):
    """One gene's category from its symmetric-class row (None when the gene is not classified).

    shared: conserved or reversal in some contrast; mouse-led / human-led: zonated in one species only there;
    mixed: mouse-led in one contrast and human-led in another; otherwise not zonated.
    """
    if row is None:
        return 'not classifiable'
    categories = []
    for contrast in contrasts:
        klass = row.get(contrast)
        if not isinstance(klass, str):
            continue
        mouse_zonated = row[contrast + ' mouse call'] in ZONATED_CALLS
        human_zonated = row[contrast + ' human call'] in ZONATED_CALLS
        if klass in ('conserved', 'reversal'):
            return 'shared'
        if mouse_zonated and not human_zonated:
            categories.append('mouse-led')
        elif human_zonated and not mouse_zonated:
            categories.append('human-led')
    if 'mouse-led' in categories and 'human-led' in categories:
        return 'mixed'
    return categories[0] if categories else 'not zonated'


def call_categories(gene_sets, ids, contributing, classes, *, min_zonated=3, share=.5):
    """Notebook 43 §7's pathway label from its contributing members (no fallback to all members).

    Fewer than ``min_zonated`` zonated contributing members: unresolved. Otherwise a summary of mouse zonation
    (>= ``share`` mouse-led), shared, human-led (each >= ``share``) or mixed.
    """
    contributing = set(contributing)
    rows = []
    for p in ids:
        members = [g for g in gene_sets[p] if g in contributing]
        counts = pd.Series([member_category(classes.loc[g] if g in classes.index else None) for g in members],
                           dtype=object).value_counts()
        zonated = int(counts.reindex(list(LED)).fillna(0).sum())
        shares = counts / zonated if zonated else counts * 0
        if zonated < min_zonated:
            label = 'unresolved'
        elif shares.get('mouse-led', 0) >= share:
            label = 'summary of mouse zonation'
        elif shares.get('shared', 0) >= share:
            label = 'shared'
        elif shares.get('human-led', 0) >= share:
            label = 'human-led'
        else:
            label = 'mixed'
        rows.append({'pathway_id': p, 'contributing_members': len(members), 'zonated_contributing': zonated,
                     **{f'n_{k}': int(counts.get(k, 0)) for k in (*LED, 'not zonated', 'not classifiable')},
                     'category': label})
    return pd.DataFrame(rows).set_index('pathway_id')


# --- Location and direction (notebooks 37 R3 and 47 e) -----------------------------------------------------------

def call_peaks(z, gene_sets, ids, grid):
    """Peak of each set's local rank divergence (notebook 37's descriptor) and the direction there.

    ``z`` is a genes x grid frame of pointwise human − mouse z. Genes with a non-finite value anywhere leave the
    ranking, and sets keep their remaining members. Returns peak_position, peak_divergence,
    median_member_z_at_peak, human_high and members_used per set.
    """
    z = pd.DataFrame(z)
    usable = np.isfinite(z.to_numpy(float)).all(axis=1)
    genes = list(z.index[usable])
    keep = set(genes)
    sets = {p: [g for g in gene_sets[p] if g in keep] for p in ids}
    sets = {p: v for p, v in sets.items() if 0 < len(v) < len(genes)}
    curves = local_pathway_curves(z.to_numpy(float)[usable], genes, sets, np.asarray(grid, float))
    out = spatial_descriptions(curves).set_index('pathway_id')[['peak_position', 'peak_divergence',
                                                                'median_member_z_at_peak']]
    out['human_high'] = out.median_member_z_at_peak.gt(0)
    out['members_used'] = pd.Series({p: len(v) for p, v in sets.items()})
    return out.reindex(list(ids))


def mean_transitions(transitions, specimens):
    """Mean S1→S2 and S2→S3 over ``specimens`` (columns specimen, transition, position); None if not increasing."""
    block = transitions[transitions.specimen.isin(list(specimens))]
    t = block.groupby('transition').position.mean()
    t12, t23 = float(t['S1→S2']), float(t['S2→S3'])
    return (t12, t23) if 0 < t12 < t23 < 1 else None


def to_segment_units(position, cuts, lo=0., hi=1.):
    """Notebook 47 (e)'s segment units (NaN everywhere when ``cuts`` is None)."""
    position = np.asarray(position, float)
    if cuts is None:
        return np.full(position.shape, np.nan)
    return segment_units(position, cuts[0], cuts[1], lo, hi)


def segment_index(units):
    """0, 1, 2 for S1, S2, S3 from segment units (NaN stays NaN)."""
    units = np.asarray(units, float)
    return np.where(np.isfinite(units), np.clip(np.floor(units), 0, 2), np.nan)


# --- Gene-removal checks (notebook 37's rule on a fixed gene statistic) ------------------------------------------

def set_effect(statistic, members):
    """Mean statistic of the members minus the mean of every other gene."""
    inside = statistic.index.isin(list(members))
    return float(statistic[inside].mean() - statistic[~inside].mean())


def subset_retention(statistic, gene_sets, ids, keep, *, camera, q=.10, ratio=.7, sizes=(10, 300)):
    """Notebook 37's gene-subset check, as a set test on a fixed statistic restricted to ``keep``.

    Sets keep their members in ``keep`` and are tested when 10–300 remain (at most |keep| − 1). ``camera`` is a
    set test returning pathway_id, direction and fdr (``camera_pr``). Retained: Up, FDR <= ``q`` over the
    subset's tested sets, and effect ratio (restricted / full) >= ``ratio``.
    """
    statistic = pd.Series(statistic, dtype=float)
    keep = set(keep)
    sub = statistic[statistic.index.isin(keep)]
    restricted = {p: [g for g in v if g in keep] for p, v in gene_sets.items()}
    restricted = {p: v for p, v in restricted.items() if sizes[0] <= len(v) <= min(sizes[1], len(sub) - 1)}
    tested = camera(sub, restricted).set_index('pathway_id')
    rows = []
    for p in ids:
        full = set_effect(statistic, gene_sets[p])
        if p not in tested.index:
            rows.append({'pathway_id': p, 'testable': False, 'members_kept': len([g for g in gene_sets[p] if g in keep]),
                         'direction': None, 'fdr': np.nan, 'effect_ratio': np.nan, 'retained': False})
            continue
        effect = set_effect(sub, restricted[p])
        r = effect / full if full != 0 else np.nan
        rows.append({'pathway_id': p, 'testable': True, 'members_kept': len(restricted[p]),
                     'direction': tested.loc[p, 'direction'], 'fdr': float(tested.loc[p, 'fdr']), 'effect_ratio': r,
                     'retained': bool(tested.loc[p, 'direction'] == 'Up' and tested.loc[p, 'fdr'] <= q
                                      and np.isfinite(r) and r >= ratio)})
    return pd.DataFrame(rows).set_index('pathway_id')


# --- External references (notebooks 37 E3 and 50) ----------------------------------------------------------------

def census_donors(table, meta, segment_map, early_parts, genes, *, min_cells=50):
    """Per-donor log2 CPM of late (S3), early and whole PT from Census segment pseudobulks (notebook 50)."""
    meta = meta[~meta.development_stage.str.contains('embryo|E16|Theiler', case=False)].copy()
    meta['segment'] = meta.cell_type.map(segment_map)
    meta = meta[meta.segment.notna() & meta.n_cells.ge(min_cells)]
    needed = set(segment_map.values())
    complete = meta.groupby('donor_id').segment.apply(lambda s: needed.issubset(set(s)))
    meta = meta[meta.donor_id.isin(complete[complete].index)]
    late, early, whole, info = {}, {}, {}, []

    def log_cpm_of(index):
        return np.log2(table.loc[index].sum() / meta.loc[index, 'total_counts_all_genes'].sum() * 1e6 + 1)

    for (donor, sex), block in meta.groupby(['donor_id', 'sex']):
        late_index, early_index = block.index[block.segment.eq('S3')], block.index[block.segment.isin(early_parts)]
        late[donor], early[donor] = log_cpm_of(late_index), log_cpm_of(early_index)
        whole[donor] = log_cpm_of(late_index.union(early_index))
        info.append({'donor': donor, 'sex': sex})
    frame = lambda values: pd.DataFrame(values).T.reindex(columns=list(genes))   # noqa: E731
    return frame(late), frame(early), frame(whole), pd.DataFrame(info).set_index('donor')


def lake_donors(sums, totals, groups, genes, *, min_nuclei=50, nuclei=None):
    """Per-donor log2 CPM (late = S3, early = S1 + S2, whole) from Lake donor|sex|segment sums on ``genes``."""
    meta = pd.DataFrame([g.split('|') for g in groups], columns=['donor', 'sex', 'segment'])
    meta['nuclei'] = nuclei
    table = meta.pivot_table(index='donor', columns='segment', values='nuclei', aggfunc='sum', fill_value=0)
    donors = table.index[(table.reindex(columns=['S1', 'S2', 'S3'], fill_value=0) >= min_nuclei).all(axis=1)]
    frame = pd.DataFrame(sums, columns=list(genes))
    totals = np.asarray(totals, float)

    def log_cpm(donor, parts):
        rows = meta.index[meta.donor.eq(donor) & meta.segment.isin(parts)]
        return np.log2(frame.loc[rows].sum() / totals[rows].sum() * 1e6 + 1)

    late = pd.DataFrame({d: log_cpm(d, ['S3']) for d in donors}).T
    early = pd.DataFrame({d: log_cpm(d, ['S1', 'S2']) for d in donors}).T
    whole = pd.DataFrame({d: log_cpm(d, ['S1', 'S2', 'S3']) for d in donors}).T
    return late, early, whole


def our_contrast(Y, genes, specimen, segment, names):
    """Per-specimen S3 − early (S1 + S2) of the mean per-unit expression (labels S1/S2/S3 without the 'PT-')."""
    out = {}
    for name in names:
        rows = specimen == name
        late = np.asarray(Y[rows & (segment == 'S3')].mean(axis=0)).ravel()
        early = np.asarray(Y[rows & np.isin(segment, ['S1', 'S2'])].mean(axis=0)).ravel()
        out[name] = late - early
    return pd.DataFrame(out, index=list(genes))


def replication_inputs(ours_contrast, human, mouse, covariates, gene_sets, human_late, human_early, mouse_late,
                       mouse_early, mouse_donors):
    """Notebook 50's ``reference_inputs`` with notebook 37's extras (centred contrasts, components)."""
    ours_raw = ours_contrast[human].mean(axis=1) - ours_contrast[mouse].mean(axis=1)
    mouse_gradient = ours_contrast[mouse].mean(axis=1)
    hd = human_late - human_early
    md = (mouse_late - mouse_early).loc[mouse_donors]
    level_h = pd.concat([human_late.mean(), human_early.mean()], axis=1).max(axis=1)
    level_m = pd.concat([mouse_late.loc[mouse_donors].mean(), mouse_early.loc[mouse_donors].mean()], axis=1).max(axis=1)
    keep = (covariates[['detection_mouse', 'detection_human']].min(axis=1).ge(.10) & level_h.ge(1) & level_m.ge(1)
            & hd.notna().all() & md.notna().all())
    g = keep.index[keep.to_numpy()]
    centred = ours_contrast - ours_contrast.loc[g].median()
    ours_c = centred[human].mean(axis=1) - centred[mouse].mean(axis=1)
    hc, mc = hd.sub(hd[g].median(axis=1), axis=0), md.sub(md[g].median(axis=1), axis=0)
    ext_raw, ext_c = hd.mean() - md.mean(), hc.mean() - mc.mean()
    t = pd.Series(welch_t(hc[g].to_numpy(), mc[g].to_numpy())[0], index=g)
    slope_ours = np.polyfit(mouse_gradient.loc[g], ours_raw.loc[g], 1)
    ours_slope_free = ours_raw - np.polyval(slope_ours, mouse_gradient)
    beta_ext = np.polyfit(md[g].mean(), ext_raw.loc[g], 1)[0]
    t_free = pd.Series(welch_t(hc[g].to_numpy(), (1 + beta_ext) * mc[g].to_numpy())[0], index=g)
    strata = matching_strata(pd.DataFrame({
        'mean_expression': ((level_h + level_m) / 2).loc[g].to_numpy(),
        'detection': covariates.loc[g, ['detection_mouse', 'detection_human']].min(axis=1).to_numpy(),
        'coverage': covariates.loc[g, 'coverage'].to_numpy()}), bins=3)
    gset = set(g)
    sets = {p: [x for x in v if x in gset] for p, v in gene_sets.items()}
    sets = {p: v for p, v in sets.items() if 5 <= len(v) < len(g)}
    residuals = pd.concat([hc[g] - hc[g].mean(), mc[g] - mc[g].mean()])
    component = lambda inter, grad: {  # noqa: E731
        'median interaction': float(inter.loc[g].median()), 'fraction of genes negative': float((inter.loc[g] < 0).mean()),
        'slope of interaction on mouse gradient': float(np.polyfit(grad.loc[g], inter.loc[g], 1)[0]),
        'R² of that slope': float(np.corrcoef(grad.loc[g], inter.loc[g])[0, 1] ** 2)}
    return {'genes': g, 'signed': np.sign(ours_c.loc[g]) * t, 'slope-free': np.sign(ours_slope_free.loc[g]) * t_free,
            'strata': strata, 'sets': sets, 'residuals': residuals.to_numpy(), 'ours_c': ours_c, 'ext_c': ext_c,
            'ours_component': {**component(ours_raw, mouse_gradient),
                               'median S3 − early, human': float(ours_contrast.loc[g, human].mean(axis=1).median()),
                               'median S3 − early, mouse': float(mouse_gradient.loc[g].median())},
            'component': {**component(ext_raw, md.mean()), 'median S3 − early, human': float(hd[g].mean().median()),
                          'median S3 − early, mouse': float(md[g].mean().median())}}


def scored_sets(ref, pathway_ids, programs):
    """The pathways (in the reference's testable sets) plus program unions with >= 5 informative members.

    ``programs``: dict kind -> Series pathway_id -> program number.
    """
    sets_g = ref['sets']
    out = {p: sets_g[p] for p in sorted(set(pathway_ids)) if p in sets_g}
    for kind, assignment in programs.items():
        for program, ids in assignment.groupby(assignment).groups.items():
            union = sorted({g for p in ids if p in sets_g for g in sets_g[p]})
            if len(union) >= 5:
                out[f'program:{kind}:{int(program)}'] = union
    return out


def selection_strata(ref, score):
    """Notebook 50's selection-matched strata: matching strata x the quintile of ``score`` among informative genes."""
    quintile = pd.qcut(pd.Series(score).loc[ref['genes']].rank(method='first'), 5, labels=False).to_numpy()
    return ref['strata'] * 5 + quintile


def collection_run(ref, statistic, sets, collections, strata, *, n_perm, seed):
    """Observed and permuted collection means, and every set's joint z (notebook 50's ``run_null``)."""
    scorer = JointReplicationZ(ref[statistic], strata, sets, ref['residuals'])
    observed, draws = collection_null(scorer, collections, strata, n_perm, seed)
    return {'observed': observed.to_numpy(), 'draws': draws.to_numpy(),
            'collections': np.array(list(observed.index), dtype=str), 'set_ids': np.array(scorer.ids, dtype=str),
            'set_z': scorer.z()}


def collection_summary(run, firsts, comparators, program_kinds, *, seed):
    """Notebook 50's summary row: means, p_random, D and p against each comparator, and program-level tests.

    ``program_kinds`` = (first kind, comparator kind) of the program collections, named '<kind> programs'.
    """
    collections = list(run['collections'])
    observed = pd.Series(run['observed'], index=collections)
    draws = pd.DataFrame(run['draws'], columns=collections)
    set_z = pd.Series(run['set_z'], index=run['set_ids'])
    first_kind, comparator_kind = program_kinds
    row = {f'mean z, {name}': observed[name] for name in collections}
    for name in (*firsts, f'{first_kind} programs', *comparators):
        row[f'p_random, {name}'] = upper_p(observed[name], draws[name])
    for first in firsts:
        for comparator in comparators:
            diff = observed[first] - observed[comparator]
            row[f'D {first} vs {comparator}'] = diff
            row[f'p {first} vs {comparator}'] = upper_p(diff, draws[first] - draws[comparator])
    first_z = set_z[[i for i in set_z.index if i.startswith(f'program:{first_kind}:')]]
    other_z = set_z[[i for i in set_z.index if i.startswith(f'program:{comparator_kind}:')]]
    row['programs (first / comparator)'] = f'{len(first_z)} / {len(other_z)}'
    row['p programs (label permutation)'] = (label_permutation_p(first_z, other_z, seed=seed)[0]
                                             if len(first_z) and len(other_z) else np.nan)
    row['p programs (gene relabeling)'] = upper_p(observed[f'{first_kind} programs'] - observed[f'{comparator_kind} programs'],
                                                  draws[f'{first_kind} programs'] - draws[f'{comparator_kind} programs'])
    return row


def replication_decision(row, first, comparator):
    """Notebook 48/50's rule: specific, better than random only, or not replicated."""
    beats_random = row[f'p_random, {first}'] <= .01
    beats_comparator = row[f'p {first} vs {comparator}'] <= .05 and row['p programs (label permutation)'] <= .05
    return 'specific' if beats_random and beats_comparator else 'better than random only' if beats_random else 'not replicated'


def per_call_replication(set_z, ids, *, q=.10):
    """Notebook 37 E3's count: BH <= ``q`` within the list on the one-sided p of the joint z, with z > 0."""
    z = pd.Series(set_z).reindex([p for p in ids if p in set(pd.Series(set_z).index)])
    p = norm.sf(z.to_numpy(float))
    fdr = pd.Series(multipletests(p, method='fdr_bh')[1], index=z.index)
    return pd.DataFrame({'z_joint': z, 'p_joint': p, 'q_within_list': fdr, 'replicated': fdr.le(q) & z.gt(0)})


def gene_agreement(ref, score, others=None, *, n_draws=1000, seed=12):
    """Notebook 37 E3's gene-level agreement on the top 10% of ``score`` among informative genes."""
    g = ref['genes']
    ours_c, ext_c = ref['ours_c'], ref['ext_c']
    score = pd.Series(score)
    top = g[score.loc[g].ge(score.loc[g].quantile(.9)).to_numpy()]
    rng = np.random.default_rng(seed)
    strata = pd.Series(ref['strata'], index=g)
    pools = {k: strata.index[strata.eq(k)] for k in strata.unique()}
    need = strata.loc[top].value_counts()
    draws = []
    for _ in range(n_draws):
        chosen = np.concatenate([rng.choice(pools[k], n, replace=False) for k, n in need.items()])
        draws.append(spearmanr(ours_c.loc[chosen], ext_c.loc[chosen]).statistic)

    def by(values):
        values = pd.Series(values)
        pick = g[values.loc[g].ge(values.loc[g].quantile(.9)).to_numpy()]
        return spearmanr(ours_c.loc[pick], ext_c.loc[pick]).statistic

    out = {'informative genes': len(g), 'gene ρ, all informative': spearmanr(ours_c.loc[g], ext_c.loc[g]).statistic,
           'gene ρ, top 10% of the score': spearmanr(ours_c.loc[top], ext_c.loc[top]).statistic,
           'gene ρ, expression-matched random (median)': float(np.median(draws)),
           'gene ρ, expression-matched random (97.5th pct)': float(np.percentile(draws, 97.5)),
           'sign agreement, top 10% of the score': float(np.mean(np.sign(ours_c.loc[top]) == np.sign(ext_c.loc[top])))}
    for label, values in (others or {}).items():
        out[f'gene ρ, top 10% {label}'] = by(values)
    return out


def average_level_inputs(genes, gene_sets, covariates, human_whole, mouse_whole, mouse_donors, whole_lfc):
    """Notebook 50's average-level test inputs for conventional-only sets (whole-PT human − mouse donors)."""
    info = genes
    hw = human_whole.sub(human_whole[info].median(axis=1), axis=0)
    mw = mouse_whole.loc[mouse_donors].sub(mouse_whole.loc[mouse_donors, info].median(axis=1), axis=0)
    t_main = pd.Series(welch_t(hw[info].to_numpy(), mw[info].to_numpy())[0], index=info)
    infoset = set(info)
    sets = {p: [x for x in v if x in infoset] for p, v in gene_sets.items()}
    sets = {p: v for p, v in sets.items() if 5 <= len(v) < len(info)}
    residuals = pd.concat([hw[info] - hw[info].mean(), mw[info] - mw[info].mean()])
    statistic = np.sign(whole_lfc.reindex(info).fillna(0)) * t_main
    strata = matching_strata(pd.DataFrame({
        'mean_expression': ((human_whole[info].mean() + mouse_whole.loc[mouse_donors, info].mean()) / 2).to_numpy(),
        'detection': covariates.loc[info, ['detection_mouse', 'detection_human']].min(axis=1).to_numpy(),
        'coverage': covariates.loc[info, 'coverage'].to_numpy()}), bins=3)
    return {'statistic': statistic, 'sets': sets, 'residuals': residuals.to_numpy(), 'strata': strata}


# --- Lead genes (notebook 38's "ours" arm) -------------------------------------------------------------------------

def ours_test(claim, ours_contrast, human, mouse):
    """Notebook 38's ``ours`` test for one claim row (mouse, human, interaction arms on S3 − S1 per specimen).

    Every claimed arm must hold in the species mean, in each specimen (signed arms) and in all four human–mouse
    pairings (interaction). Returns True, False, or None when nothing is claimed or measurable.
    """
    h = float(ours_contrast.loc[human].mean())
    m = float(ours_contrast.loc[mouse].mean())
    parts = [arm_passes(claim['mouse'], m, h), arm_passes(claim['human'], h, m), arm_passes(claim['interaction'], h - m)]
    for species, names in (('mouse', mouse), ('human', human)):
        if claim[species] not in ('', None, '0'):
            parts.append(all(arm_passes(claim[species], ours_contrast.loc[n]) for n in names))
    if claim['interaction'] not in ('', None):
        parts.append(all(arm_passes(claim['interaction'], ours_contrast.loc[a] - ours_contrast.loc[b])
                         for a in human for b in mouse))
    parts = [p for p in parts if p is not None]
    return all(parts) if parts else None
