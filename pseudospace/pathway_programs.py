"""Post-discovery interpretation for notebook 12; never refit or retest discovery."""
from collections import Counter
from hashlib import sha256
from itertools import combinations

import numpy as np
import pandas as pd


def overlap_sensitivity(statistics, gene_sets):
    """Downweight members by 1/sqrt(annotation frequency), keeping background fixed.

    Count each gene once per tested term across the full library collection.
    This is a descriptive weighted AUC, not PADOG or a new significance test.
    """
    if not statistics.index.is_unique or not np.isfinite(statistics).all():
        raise ValueError('Need unique genes and finite spatial statistics.')
    frequency = Counter(g for members in gene_sets.values() for g in set(members))
    if not set(frequency).issubset(statistics.index):
        raise ValueError('Every pathway member must belong to the tested universe.')
    annotation = pd.DataFrame({'gene': statistics.index})
    annotation['annotation_frequency'] = annotation.gene.map(frequency).fillna(0).astype(int)
    annotation['member_weight'] = 1 / np.sqrt(annotation.annotation_frequency.clip(lower=1))
    rows = []
    for pathway, members in gene_sets.items():
        members = sorted(set(members))
        background = np.sort(statistics.loc[~statistics.index.isin(members)].to_numpy())
        if not members or not len(background):
            raise ValueError('Need nonempty membership and complement.')
        values = statistics.loc[members].to_numpy()
        # Each member's probability of beating the SAME nonmember background; ties get half credit.
        wins = (np.searchsorted(background, values, 'left') +
                np.searchsorted(background, values, 'right')) / (2 * len(background))
        weights = np.array([1 / np.sqrt(frequency[g]) for g in members])
        original, weighted = wins.mean() - .5, np.average(wins, weights=weights) - .5
        rows.append((pathway, original, weighted, weighted - original))
    result = pd.DataFrame(rows, columns=['pathway_id', 'original_effect', 'overlap_weighted_effect',
                                         'overlap_effect_change'])
    result['original_effect_rank'] = result.original_effect.rank(ascending=False, method='min')
    result['overlap_effect_rank'] = result.overlap_weighted_effect.rank(ascending=False, method='min')
    result['overlap_rank_change'] = result.overlap_effect_rank - result.original_effect_rank
    return result, annotation


def active_programs(candidates, gene_sets, members, curves, *, distance_cut=.55):
    """Complete linkage: 50% active Jaccard, 25% D correlation, 25% S correlation.

    Negative/undefined correlations contribute zero. Full membership is reported
    but never used in clustering. Stable IDs hash the exact sorted member terms.
    """
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform

    ids = sorted(candidates)
    if len(ids) != len(set(ids)) or not 0 <= distance_cut <= 1:
        raise ValueError('Need unique candidates and a distance cut between zero and one.')
    roles = ['leading_edge', 'broad_supporter', 'spatial_driver']
    selected = members[members.pathway_id.isin(ids) & members[roles].any(axis=1)]
    active = {p: set(selected.loc[selected.pathway_id.eq(p), 'gene']) for p in ids}
    traces = {}
    for pathway in ids:
        trace = curves.loc[curves.pathway_id.eq(pathway)].sort_values('position')
        if (len(trace) < 2 or trace.position.duplicated().any() or not active[pathway] or
                not active[pathway].issubset(gene_sets[pathway]) or
                not np.isfinite(trace[['position', 'divergence', 'direction']]).all().all()):
            raise ValueError('Every candidate needs valid aligned curves and nonempty active genes.')
        traces[pathway] = trace
    def jaccard(left, right):
        return len(left & right) / len(left | right) if left | right else 0.
    distances = np.zeros((len(ids), len(ids)))
    rows = []
    # Ponytail: O(P²) pair comparisons suit hundreds of candidates; use sparse neighbors at atlas scale.
    for i, j in combinations(range(len(ids)), 2):
        left, right = ids[i], ids[j]
        a, b = traces[left], traces[right]
        if not np.array_equal(a.position, b.position):
            raise ValueError('Pathway curves must share the same coordinate grid.')
        correlations = []
        for column in ('divergence', 'direction'):
            correlation = (np.corrcoef(a[column], b[column])[0, 1]
                           if a[column].std() > 1e-10 and b[column].std() > 1e-10 else 0.)
            correlations.append(float(np.clip(correlation, 0, 1)))
        overlap = jaccard(active[left], active[right])
        similarity = .5 * overlap + .25 * sum(correlations)
        distances[i, j] = distances[j, i] = 1 - similarity
        rows.append((left, right, overlap, jaccard(set(gene_sets[left]), set(gene_sets[right])),
                     *correlations, similarity, similarity >= 1 - distance_cut))
    labels = (fcluster(linkage(squareform(distances), method='complete'), distance_cut,
                       criterion='distance') if len(ids) > 1 else np.ones(len(ids), int))
    stable_ids = {label: 'P_' + sha256('\n'.join(p for p, c in zip(ids, labels) if c == label)
                                     .encode()).hexdigest()[:10] for label in set(labels)}
    membership = pd.DataFrame({'pathway_id': ids, 'paper_program': [stable_ids[c] for c in labels]})
    edges = pd.DataFrame(rows, columns=['source', 'target', 'active_jaccard', 'full_jaccard',
                                       'divergence_correlation', 'direction_correlation',
                                       'similarity', 'network_edge'])
    return membership, edges, active


def representative_terms(membership, atlas, active):
    """Transparent draft representatives; biological interpretability needs review.

    Equal percentile weights: spatial effect, tested/requested coverage, retention
    over ALL planned runs, and coverage of the shared active core. A second term
    must add at least 10% of the program's active union. Missing robustness = 0.
    """
    terms = membership.merge(atlas, on='pathway_id', validate='one_to_one').copy()
    if len(terms) != len(membership):
        raise ValueError('Every candidate needs a row in the frozen evidence atlas.')
    terms['planned_retention'] = terms.n_retained.div(terms.n_planned.replace(0, np.nan)).fillna(0)
    terms['paper_priority'] = (terms.correlation_support.eq('Correlation-supported') &
                               terms.planned_retention.ge(.75))
    terms['active_core_coverage'] = 0.
    terms['representative_score'] = 0.
    terms['representative_order'] = 0
    for _, group in terms.groupby('paper_program', sort=True):
        counts = Counter(g for p in group.pathway_id for g in active[p])
        core = {g for g, n in counts.items() if n >= 2} or set(counts)
        terms.loc[group.index, 'active_core_coverage'] = [len(active[p] & core) / len(core)
                                                        for p in group.pathway_id]
        columns = ['effect_T_spatial', 'tested_fraction', 'planned_retention', 'active_core_coverage']
        scores = terms.loc[group.index, columns].fillna(0).rank(pct=True).mean(axis=1)
        terms.loc[group.index, 'representative_score'] = scores
        order = terms.loc[group.index].sort_values(['representative_score', 'pathway_id'],
                                                   ascending=[False, True])
        first = order.index[0]
        terms.loc[first, 'representative_order'] = 1
        first_genes = active[terms.loc[first, 'pathway_id']]
        novelty = order.pathway_id.map(lambda p: len(active[p] - first_genes) / len(counts))
        alternatives = order.loc[novelty.ge(.1)]
        if len(alternatives):
            terms.loc[alternatives.index[0], 'representative_order'] = 2
    return terms


def program_gene_evidence(terms, members, genes, grid, fit, *, z_threshold=2., fraction=.2):
    """Active genes and pointwise bidirectionality, without signed-mean cancellation."""
    index = pd.Index(genes)
    grid = np.asarray(grid)
    if (not index.is_unique or len(grid) < 2 or np.any(np.diff(grid) <= 0) or
            not np.isfinite(grid).all() or not 0 < fraction <= 1 or z_threshold <= 0):
        raise ValueError('Need unique genes, increasing finite grid, and positive direction thresholds.')
    for key in ('z', 'human', 'mouse', 'delta'):
        if np.shape(fit[key]) != (len(index), len(grid)) or not np.isfinite(fit[key]).all():
            raise ValueError('Gene trajectories must be finite and align with genes and grid.')
    roles = ['leading_edge', 'broad_supporter', 'spatial_driver']
    active_members = members[members[roles].any(axis=1)].drop_duplicates(['pathway_id', 'gene'])
    rows, directions = [], []
    for program, group in terms.groupby('paper_program', sort=True):
        detail = active_members[active_members.pathway_id.isin(group.pathway_id)]
        for gene, support in detail.groupby('gene', sort=True):
            if gene not in index:
                raise ValueError('An active gene is absent from the fitted universe.')
            idx = index.get_loc(gene)
            peak = int(np.argmax(abs(fit['z'][idx])))
            role = ('single-term program' if len(group) == 1 else
                    'shared active driver' if support.pathway_id.nunique() >= 2 else 'term-specific active gene')
            rows.append({'paper_program': program, 'gene': gene, 'gene_role': role,
                         'n_active_terms': support.pathway_id.nunique(),
                         'active_terms': ';'.join(sorted(support.pathway_id.unique())),
                         'peak_position': grid[peak], 'peak_abs_z': abs(fit['z'][idx, peak]),
                         'peak_delta': fit['delta'][idx, peak],
                         'bulk_logcpm_effect': support.bulk_logcpm_effect.iloc[0],
                         **{r: bool(support[r].any()) for r in roles}})
        active_genes = sorted(detail.gene.unique())
        if not active_genes:
            raise ValueError('Every program needs at least one active gene.')
        z = fit['z'][index.get_indexer(active_genes)]
        positive, negative = (z >= z_threshold).sum(axis=0), (z <= -z_threshold).sum(axis=0)
        threshold = max(2, int(np.ceil(fraction * len(z))))
        for k, position in enumerate(grid):
            up, down = positive[k] >= threshold, negative[k] >= threshold
            label = ('mixed' if up and down else 'human-high dominated' if up else
                     'mouse-high dominated' if down else 'weak / sparse')
            directions.append((program, position, len(z), positive[k], negative[k],
                               positive[k] / len(z), negative[k] / len(z), label))
    gene_table = pd.DataFrame(rows, columns=['paper_program', 'gene', 'gene_role', 'n_active_terms',
        'active_terms', 'peak_position', 'peak_abs_z', 'peak_delta', 'bulk_logcpm_effect', *roles])
    direction_table = pd.DataFrame(directions, columns=['paper_program', 'position', 'n_active_genes',
        'n_human_high', 'n_mouse_high', 'fraction_human_high', 'fraction_mouse_high', 'direction_description'])
    return gene_table, direction_table


def program_figures(program, terms, gene_table, curves, genes, grid, fit, *, heatmap_limit=6.):
    """Yield readable evidence panels; all active genes are paginated, never dropped.

    Fits are the original equal-specimen model means, not independent donor means.
    The four example genes balance shared and term-specific roles when available.
    """
    import matplotlib.pyplot as plt

    representatives = terms[(terms.paper_program == program) & (terms.representative_order > 0)]
    representatives = representatives.sort_values('representative_order')
    detail = gene_table[gene_table.paper_program.eq(program)].sort_values(['peak_position', 'gene'])
    index = pd.Index(genes)
    style = {'font.family': 'sans-serif', 'font.sans-serif': ['DejaVu Sans'],
             'font.size': 7, 'axes.titlesize': 8, 'axes.labelsize': 7, 'xtick.labelsize': 6,
             'ytick.labelsize': 6, 'legend.fontsize': 6, 'pdf.fonttype': 42, 'svg.fonttype': 'none',
             'axes.spines.top': False, 'axes.spines.right': False}
    with plt.rc_context(style):
        fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.7), layout='constrained')
        for row, color in zip(representatives.itertuples(), ['#0072B2', '#D55E00']):
            trace = curves[curves.pathway_id.eq(row.pathway_id)].sort_values('position')
            for ax, column in zip(axes, ['divergence', 'direction']):
                ax.plot(trace.position, trace[column], color=color, label=f'R{row.representative_order}')
        for ax, label, title in zip(axes, ['D(s): unsigned excess AUC', 'S(s): relative direction'],
                                    ['Divergence', 'Direction']):
            ax.axhline(0, color='.7', linewidth=.6, zorder=0)
            ax.set(xlabel='Established PT pseudospace', ylabel=label, title=title)
            ax.set_xlim(grid[0], grid[-1])
        axes[0].set_ylim(-.5, .5)
        axes[1].set_ylim(-1, 1)
        axes[0].legend(loc='upper left', bbox_to_anchor=(0, 1.02), frameon=False, ncol=2)
        fig.suptitle(program)
        yield 'pathway_curves', fig

        # Pagination preserves every active gene while keeping the text readable at export size.
        for start in range(0, len(detail), 36):
            page = detail.iloc[start:start + 36]
            idx = index.get_indexer(page.gene)
            fig, ax = plt.subplots(figsize=(7.2, max(2.4, .13 * len(page) + 1.1)), layout='constrained')
            step = (grid[-1] - grid[0]) / (len(grid) - 1)
            artist = ax.imshow(fit['z'][idx], aspect='auto', interpolation='none', cmap='RdBu_r',
                vmin=-heatmap_limit, vmax=heatmap_limit,
                extent=(grid[0] - step / 2, grid[-1] + step / 2, len(page) - .5, -.5))
            labels = [g + (' *' if n >= 2 else '') for g, n in zip(page.gene, page.n_active_terms)]
            ax.set(yticks=np.arange(len(page)), yticklabels=labels, xlabel='Established PT pseudospace',
                   title=f'{program} | active genes {start + 1}–{start + len(page)} / {len(detail)}')
            fig.colorbar(artist, ax=ax, label='Signed Z (human − mouse)', extend='both')
            yield f'active_genes_{start // 36 + 1:02d}', fig

        ranked = detail.sort_values(['peak_abs_z', 'gene'], ascending=[False, True])
        examples = ranked.groupby('gene_role', sort=True).head(2).head(4)
        examples = pd.concat([examples, ranked]).drop_duplicates('gene').head(4)
        fig, axes = plt.subplots(2, 2, figsize=(7.2, 4.4), layout='constrained')
        for ax, row in zip(axes.flat, examples.itertuples()):
            idx = index.get_loc(row.gene)
            ax.plot(grid, fit['human'][idx], color='#D55E00', label='Human', linewidth=1.2)
            ax.plot(grid, fit['mouse'][idx], color='#0072B2', label='Mouse', linestyle='--', linewidth=1.2)
            ax.set(title=f'{row.gene} | {row.gene_role}', xlabel='Established PT pseudospace',
                   ylabel='Fitted log-normalized expression', xlim=(grid[0], grid[-1]))
        for ax in axes.flat[len(examples):]:
            ax.set_visible(False)
        fig.legend(*axes.flat[0].get_legend_handles_labels(), loc='outside upper center', ncol=2, frameon=False)
        yield 'gene_fits', fig
