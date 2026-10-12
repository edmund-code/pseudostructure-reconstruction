"""Published PT programs along the pseudospace (notebooks 79 and 80): replication, localization, added information.

Protocol: ``docs/results/pt-known-programs.md``. The programs are frozen lists from published tables
(``docs/results/known-programs/``). Per program the notebooks ask (a) whether our cohort reproduces the published
claim, (b) where along the PT the case - reference difference sits (segment level only), and (c) what that adds.

Everything here is arithmetic on specimen pseudobulks, on per-gene statistics the notebooks already have (the segment
and positional F-test z) and on fitted curves; no gene model is fitted here.

- ``read_programs``: a frozen program file, hash-checked, members inside the universe and none excluded.
- ``segment_pseudobulk``: log2 rate (sum counts + floor) / sum library per specimen and scope.
- ``centred_contrast``: a linear contrast over specimens minus its median over reference genes (notebook 76's
  centring, the positional analogue of median-of-ratios).
- ``partition_effects``: the group effect, its variance and residual df for the condition or a relabeling on the four
  pseudobulks of a scope (limma's two-group or species + relabel model in closed form).
- ``segment_contrast``: one species' segment-vs-other-segments contrast, its variance and df (zonation family).
- ``moderated_t``: limma's squeezeVar, then t on df + d0.
- ``program_test``: cameraPR of each program, with one-sided p in the claimed and the opposite direction.
- ``segment_summary``: program medians (or means) of per-gene contrasts by scope.
- ``replication_class``, ``transfer_class``, ``localization_class``: the protocol's decision rules (section 5).
- ``added_information``: the protocol's one-line template (section 4c).
- ``gene_span``, ``stratified_rank_normal``: the gene-length sensitivity (section 4a).
- ``bin_rates``: equal-count specimen bins for the member heatmaps.
"""
from __future__ import annotations

import gzip
import hashlib
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t as t_dist

from .aki_along_pt import empirical_p
from .pathway_pipelines import camera_pr, rank_normal_scores, squeeze_var

SEGMENTS = ('S1', 'S2', 'S3')
EFFECT_MIN, RATIO_MAX, ALPHA = .25, .25, .05


# --- Programs ---------------------------------------------------------------------------------------------------

def file_sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_programs(path, universe, exclude=(), sha256=None):
    """A frozen program file: ({program_id: members}, metadata frame indexed by program_id).

    Checks the sha256 (when given), that every member lies in ``universe`` and that no member is in ``exclude``.
    The members column is 'genes' (comma-separated mouse symbols, in the paper's ranking order).
    """
    if sha256 is not None and file_sha256(path) != sha256:
        raise ValueError(f'{path}: sha256 differs from the frozen value.')
    table = pd.read_csv(path).set_index('program_id')
    universe, exclude = set(map(str, universe)), set(map(str, exclude))
    programs = {}
    for pid, genes in table.genes.items():
        members = [g.strip() for g in str(genes).split(',') if g.strip()]
        outside = [g for g in members if g not in universe]
        banned = [g for g in members if g in exclude]
        if outside or banned:
            raise ValueError(f'{pid}: members outside the universe {outside[:5]} or excluded {banned[:5]}.')
        programs[str(pid)] = members
    return programs, table


# --- Pseudobulks and contrasts -------------------------------------------------------------------------------------

def segment_pseudobulk(counts, library, specimen, segment, scopes, *, genes=None, floor=.5):
    """{scope: genes x specimens frame of log2((sum counts + floor) / sum library)}.

    counts: units x genes (dense or sparse); scopes: {name: list of segment labels} (None = every unit).
    """
    specimen, segment = np.asarray(specimen).astype(str), np.asarray(segment).astype(str)
    library = np.asarray(library, float)
    n = counts.shape[0]
    if specimen.shape != (n,) or segment.shape != (n,) or library.shape != (n,) or (library <= 0).any():
        raise ValueError('Need aligned specimens, segments and positive libraries.')
    names = sorted(set(specimen))
    index = pd.RangeIndex(counts.shape[1]) if genes is None else pd.Index(genes, name='gene')
    out = {}
    for scope, labels in scopes.items():
        inside = np.ones(n, bool) if labels is None else np.isin(segment, list(labels))
        columns = {}
        for name in names:
            rows = inside & (specimen == name)
            if not rows.any():
                raise ValueError(f'No unit of {name} in scope {scope}.')
            total = np.asarray(counts[rows].sum(axis=0), float).ravel()
            columns[name] = np.log2((total + floor) / library[rows].sum())
        out[scope] = pd.DataFrame(columns, index=index)
    return out


def centred_contrast(rates, weights, reference=None):
    """sum_j weights[j] * rates[:, j], minus its median over the ``reference`` genes (all genes when None)."""
    weights = pd.Series(weights, dtype=float)
    value = rates[weights.index].to_numpy(float) @ weights.to_numpy()
    ref = np.ones(len(value), bool) if reference is None else np.asarray(reference, bool)
    return pd.Series(value - np.median(value[ref]), index=rates.index)


def partition_effects(rates, case, reference, partition='condition', ref_genes=None):
    """Effect, s2 (variance of the effect) and residual df on the four pseudobulks of one scope (protocol 3.3).

    case = (c1, c2), reference = (r1, r2) specimen names. condition: effect = mean(case) - mean(reference),
    s2 = (d_case^2 + d_ref^2) / 4 on 2 df (d_case = c1 - c2, d_ref = r1 - r2). 'relabeling 1' puts (c1, r1) in group 1:
    effect (d_case + d_ref) / 2, s2 = (d_case - d_ref)^2 / 4 on 1 df; 'relabeling 2' puts (c1, r2) in group 1: effect
    (d_case - d_ref) / 2, s2 = (d_case + d_ref)^2 / 4. Every contrast is centred by its own median over ``ref_genes``.
    """
    (c1, c2), (r1, r2) = case, reference
    centre = lambda w: centred_contrast(rates, w, ref_genes)          # noqa: E731
    if partition == 'condition':
        effect = centre({c1: .5, c2: .5, r1: -.5, r2: -.5})
        d_case, d_ref = centre({c1: 1., c2: -1.}), centre({r1: 1., r2: -1.})
        return pd.DataFrame({'effect': effect, 's2': (d_case ** 2 + d_ref ** 2) / 4, 'df': 2.})
    if partition == 'relabeling 1':
        effect = centre({c1: .5, c2: -.5, r1: .5, r2: -.5})
        residual = centre({c1: 1., c2: -1., r1: -1., r2: 1.})
    elif partition == 'relabeling 2':
        effect = centre({c1: .5, c2: -.5, r1: -.5, r2: .5})
        residual = centre({c1: 1., c2: -1., r1: 1., r2: -1.})
    else:
        raise ValueError("partition must be 'condition', 'relabeling 1' or 'relabeling 2'.")
    return pd.DataFrame({'effect': effect, 's2': residual ** 2 / 4, 'df': 1.})


def pair_contrasts(rates, case, reference, ref_genes=None):
    """{'c-r': centred c - r} for every case x reference pair."""
    return {f'{c}−{r}': centred_contrast(rates, {c: 1., r: -1.}, ref_genes) for c in case for r in reference}


def segment_contrast(rates, specimens, segment, segments=SEGMENTS, ref_genes=None):
    """One species' contrast of ``segment`` against the mean of the other segments (protocol 3.4).

    rates: {segment: genes x specimens}; specimens: the species' two specimens. Each specimen's contrast is centred by
    its median over ``ref_genes``; effect = their mean, s2 = (c1 - c2)^2 / 4 on 1 df. Returns effect, s2, df and the
    per-specimen contrasts.
    """
    others = [s for s in segments if s != segment]
    if segment not in segments or len(specimens) != 2:
        raise ValueError('Need a segment among ``segments`` and two specimens.')
    per = {}
    ref = None if ref_genes is None else np.asarray(ref_genes, bool)
    for name in specimens:
        value = rates[segment][name] - np.mean([rates[s][name] for s in others], axis=0)
        value = value.to_numpy(float)
        per[name] = value - np.median(value if ref is None else value[ref])
    a, b = (per[name] for name in specimens)
    out = pd.DataFrame({'effect': (a + b) / 2, 's2': (a - b) ** 2 / 4, 'df': 1.}, index=rates[segment].index)
    for name in specimens:
        out[name] = per[name]
    return out


def moderated_t(effect, s2, df, covariate=None, robust=False):
    """limma-style moderated t: squeezeVar(s2, df, covariate), t = effect / sqrt(s2_post) on df + d0, two-sided p."""
    effect = pd.Series(effect, dtype=float)
    squeezed = squeeze_var(pd.Series(np.asarray(s2, float), index=effect.index), df, covariate=covariate, robust=robust)
    if not (squeezed.var_post > 0).all():
        raise ValueError('A moderated variance is 0.')
    total = float(df) + squeezed.df_prior.to_numpy()
    t = effect.to_numpy() / np.sqrt(squeezed.var_post.to_numpy())
    p = 2 * t_dist.sf(np.abs(t), np.where(np.isinf(total), 1e12, total))
    return pd.DataFrame({'effect': effect.to_numpy(), 's2': squeezed['var'].to_numpy(), 's2_post': squeezed.var_post.to_numpy(),
                         'd0': squeezed.df_prior.to_numpy(), 'df_total': total, 't': t, 'p': p}, index=effect.index)


# --- Set tests -------------------------------------------------------------------------------------------------------

def program_test(statistic, programs, directions=None, *, inter_gene_cor=.01):
    """cameraPR of each program on one gene statistic, against all other genes.

    directions: {program: +1 or -1} (the claimed sign of the statistic; +1 when None). Returns per program n,
    direction (limma's), p (two-sided), p_claimed and p_opposite (one-sided p in the claimed and the opposite
    direction: p/2 when limma's direction matches, else 1 - p/2).
    """
    statistic = pd.Series(statistic, dtype=float)
    sets = {k: [g for g in v if g in statistic.index] for k, v in programs.items()}
    sets = {k: v for k, v in sets.items() if len(v) >= 2}
    if not sets:
        raise ValueError('No program has two members in the statistic.')
    table = camera_pr(statistic, sets, inter_gene_cor=inter_gene_cor).set_index('pathway_id')
    sign = pd.Series({k: (directions or {}).get(k, 1) for k in table.index}, dtype=float)
    up = table.direction.eq('Up').to_numpy()
    matches = np.where(sign.to_numpy() > 0, up, ~up)
    p = table.p.to_numpy(float)
    out = pd.DataFrame({'n': table.n, 'direction': table.direction, 'p': p,
                        'p_claimed': np.where(matches, p / 2, 1 - p / 2),
                        'p_opposite': np.where(matches, 1 - p / 2, p / 2)}, index=table.index)
    return out.reindex(list(programs))


def calibrated(p, null):
    """Empirical p of ``p`` against a pooled null sample of p (wraps ``aki_along_pt.empirical_p``)."""
    return empirical_p(p, null)


def gene_span(fasta_gz, symbols=None):
    """Genomic span (max end - min start + 1) per gene symbol (upper case) from Ensembl cDNA FASTA headers.

    Headers carry 'chromosome:<assembly>:<chr>:<start>:<end>:<strand> gene:<id> ... gene_symbol:<symbol>'. Only header
    lines are read. ``symbols`` restricts the output (case-insensitive).
    """
    pattern = re.compile(r':[^:\s]+:[^:\s]+:(\d+):(\d+):-?1 gene:\S+ .*?gene_symbol:(\S+)')
    lo, hi = {}, {}
    with gzip.open(fasta_gz, 'rt') as handle:
        for line in handle:
            if line.startswith('>'):
                m = pattern.search(line)
                if m:
                    start, end, symbol = int(m.group(1)), int(m.group(2)), m.group(3).upper()
                    lo[symbol] = min(lo.get(symbol, start), start)
                    hi[symbol] = max(hi.get(symbol, end), end)
    span = pd.Series({k: hi[k] - lo[k] + 1 for k in lo}, dtype=float)
    if symbols is not None:
        span = span.reindex(pd.Index(symbols).astype(str).str.upper())
    return span


def stratified_rank_normal(statistic, strata):
    """``rank_normal_scores`` within each stratum; genes with a missing stratum are dropped."""
    statistic, strata = pd.Series(statistic, dtype=float), pd.Series(strata, index=pd.Series(statistic).index)
    keep = strata.notna()
    out = [rank_normal_scores(statistic[keep & strata.eq(level)]) for level in strata[keep].unique()]
    return pd.concat(out).reindex(statistic.index[keep])


# --- Program summaries ---------------------------------------------------------------------------------------------

def segment_summary(contrasts, members, *, min_members=5, stat='median'):
    """Long table of program summaries: one row per program x contrast x scope.

    contrasts: {label: genes x scopes frame} (centred log2 contrasts); members: {program: genes}. A summary of fewer
    than ``min_members`` finite member values is NaN. Columns: program, contrast, scope, value, n, share_positive.
    """
    if stat not in ('median', 'mean'):
        raise ValueError("stat must be 'median' or 'mean'.")
    rows = []
    for label, frame in contrasts.items():
        for program, genes in members.items():
            present = [g for g in genes if g in frame.index]
            for scope in frame.columns:
                values = frame.loc[present, scope].to_numpy(float)
                values = values[np.isfinite(values)]
                ok = len(values) >= min_members
                rows.append({'program': program, 'contrast': label, 'scope': scope,
                             'value': (np.median(values) if stat == 'median' else np.mean(values)) if ok else np.nan,
                             'n': int(len(values)), 'share_positive': float((values > 0).mean()) if ok else np.nan})
    return pd.DataFrame(rows)


# --- Decision rules (protocol section 5) ---------------------------------------------------------------------------

def replication_class(q_claimed, specific, q_condition, q_pairs, q_opposite, *, effect_min=EFFECT_MIN, alpha=ALPHA):
    """Protocol 5.1. q_condition and q_pairs are program values signed in the claimed direction (Q = d * P) at the
    claim's scope; q_claimed and q_opposite are BH q of the one-sided program test in each direction."""
    pairs = np.asarray(q_pairs, float)
    if not np.isfinite(q_condition) or not np.isfinite(pairs).all():
        return 'not testable'
    if q_claimed <= alpha and specific and q_condition >= effect_min and (pairs > 0).all():
        return 'replicated'
    if q_opposite <= alpha and (pairs < 0).all():
        return 'opposite'
    if q_condition > 0 and (pairs > 0).all():
        return 'direction only'
    return 'not replicated'


def transfer_class(q, contrasts, *, effect_min=EFFECT_MIN, alpha=ALPHA):
    """Protocol 5.2 for one species: contrasts are the program's segment contrasts in its two specimens (claimed
    direction positive). Order: zonated, reversed, absent (signs differ or both within the margin), attenuated (both
    positive), else absent."""
    c = np.asarray(contrasts, float)
    if len(c) != 2 or not np.isfinite(c).all():
        return 'not testable'
    if q <= alpha and (c >= effect_min).all():
        return 'zonated'
    if (c <= -effect_min).all():
        return 'reversed'
    if np.sign(c[0]) != np.sign(c[1]) or (np.abs(c) < effect_min).all():
        return 'absent'
    if (c > 0).all():
        return 'attenuated'
    return 'absent'


def localization_class(q_between, q_condition, q_pairs, relabeled, *, effect_min=EFFECT_MIN, ratio_max=RATIO_MAX,
                       alpha=ALPHA, segments=SEGMENTS):
    """Protocol 5.3. q_condition: the program's signed values per segment (Q_k = d * P_k), length 3; q_pairs: pairs x
    segments; relabeled: relabelings x segments (signed the same way). Returns a dict: class, peak, margin, unanimous,
    specific, ratios, present (segments with Q >= effect_min in the condition and > 0 in every pair) and
    unresolved (segments with |Q| < effect_min)."""
    q = np.asarray(q_condition, float)
    pairs = np.atleast_2d(np.asarray(q_pairs, float))
    relab = np.atleast_2d(np.asarray(relabeled, float))
    if q.shape != (len(segments),) or pairs.shape[1] != len(segments) or relab.shape[1] != len(segments):
        raise ValueError('Need one value per segment for the condition, every pair and every relabeling.')
    if not np.isfinite(q).all() or not np.isfinite(pairs).all():
        return {'class': 'not testable', 'peak': None, 'margin': np.nan, 'unanimous': False, 'specific': False,
                'ratios': [np.nan] * len(relab), 'present': [], 'unresolved': []}
    peak = int(np.argmax(q))
    margin = float(q[peak] - np.delete(q, peak).max())
    unanimous = bool((pairs.argmax(axis=1) == peak).all())
    with np.errstate(divide='ignore', invalid='ignore'):
        ratios = np.abs(relab[:, peak]) / abs(q[peak]) if q[peak] != 0 else np.full(len(relab), np.inf)
    specific = bool(np.isfinite(ratios).all() and (ratios < ratio_max).all())
    present = [segments[k] for k in range(len(segments)) if q[k] >= effect_min and (pairs[:, k] > 0).all()]
    unresolved = [segments[k] for k in range(len(segments)) if abs(q[k]) < effect_min]
    signs = np.sign(np.r_[q, pairs.ravel()])
    if q_between <= alpha and specific and unanimous and margin >= effect_min:
        label = f'localized to {segments[peak]}'
    elif q_between <= alpha and specific:
        label = 'graded'
    elif (q_between > alpha and np.ptp(q) < effect_min and (signs == signs[0]).all() and signs[0] != 0
          and np.abs(q).max() >= effect_min):
        label = 'uniform within resolution'
    elif q_between > alpha and (np.abs(q) < effect_min).all():
        label = 'no difference'
    else:
        label = 'not resolved'
    return {'class': label, 'peak': segments[peak], 'margin': margin, 'unanimous': unanimous, 'specific': specific,
            'ratios': [float(r) for r in ratios], 'present': present, 'unresolved': unresolved}


def informative(replication, localization):
    """Protocol 5.4: 'position adds information', 'not informative' or 'replication only / other'."""
    if localization.startswith('localized') or localization in ('graded', 'uniform within resolution'):
        return 'position adds information'
    if replication in ('not replicated', 'direction only') and localization == 'not resolved':
        return 'not informative'
    return 'other'


def added_information(row):
    """Protocol 4c's template. ``row`` is a mapping with source, claim, replication, scope, value, pairs, q,
    localization, peak, unresolved (list), axis ('agrees' / 'differs' / 'n/a'), flags (list) and optionally
    transfer (zonation) and published_segment."""
    parts = [f"{row['source']}: {row['claim']}."]
    if row.get('transfer'):
        parts.append(f"Ours: {row['replication']} in the source species; in the other species {row['transfer']}.")
    else:
        parts.append(f"Ours: {row['replication']} at {row['scope']} (median {row['value']:+.2f} log2; "
                     f"{row['pairs']}/4 pairs; q {row['q']:.2g}).")
    loc = f"Along the PT: {row['localization']}"
    if row.get('peak') and not str(row['localization']).startswith(('localized', 'no difference')):
        loc += f", largest in {row['peak']}"
    unresolved = [s for s in row.get('unresolved', []) if s != row.get('peak')]
    if unresolved and not str(row['localization']).startswith('no difference'):
        loc += f"; {', '.join(unresolved)}: no difference resolved"
    parts.append(loc + '.')
    if row.get('axis') in ('agrees', 'differs'):
        parts.append(f"Axis {row['axis']}.")
    if row.get('flags'):
        parts.append('Flags: ' + '; '.join(row['flags']) + '.')
    return ' '.join(parts)


# --- Heatmap bins ------------------------------------------------------------------------------------------------------

def bin_rates(counts, library, position, specimen, segment, columns, *, n_bins=20, floor=.5):
    """Equal-count position bins per specimen: (log2 rate array specimens x bins x columns, majority segment frame
    specimens x bins, specimen order). ``columns`` are column indices of ``counts``."""
    position, specimen = np.asarray(position, float), np.asarray(specimen).astype(str)
    segment, library = np.asarray(segment).astype(str), np.asarray(library, float)
    names = sorted(set(specimen))
    rates = np.full((len(names), n_bins, len(columns)), np.nan)
    majority = pd.DataFrame(index=names, columns=range(n_bins), dtype=object)
    sub = counts[:, columns]
    for i, name in enumerate(names):
        rows = np.flatnonzero(specimen == name)
        order = rows[np.argsort(position[rows], kind='stable')]
        for b, chunk in enumerate(np.array_split(order, n_bins)):
            total = np.asarray(sub[chunk].sum(axis=0), float).ravel()
            rates[i, b] = np.log2((total + floor) / library[chunk].sum())
            majority.loc[name, b] = pd.Series(segment[chunk]).mode().iloc[0]
    return rates, majority, names
