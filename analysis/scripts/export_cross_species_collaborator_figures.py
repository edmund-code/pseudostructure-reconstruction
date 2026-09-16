#!/usr/bin/env python3
"""Export a collaborator brief from completed notebook 03 results.

The figure contract and inferential limits are documented in
docs/workflows/cross-species-collaborator-figures.md. Inputs are read-only;
new artifacts and a copy of this script go to the private output directory.
All displayed data are real observations or explicitly described summaries.
"""
from __future__ import annotations

import argparse
import hashlib
import html
import json
import os
from pathlib import Path
import shutil
import sys
import textwrap
from types import SimpleNamespace

import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

PROJECT = next((p for p in [Path.cwd(), *Path(__file__).resolve().parents]
                if (p / 'pseudospace').is_dir()), Path.cwd())
sys.path.insert(0, str(PROJECT))
from pseudospace.levelshape import fit_single_condition_curves, run_level_shape
from pseudospace.specimen import specimen_balanced_curves
from pseudospace.stats_gam import gam_internal_knots

COLORS = {'mouse': '#0072B2', 'human': '#D55E00'}
SAMPLES = ('Ctrl1A2', 'Ctrl1A4', 'HUK1_COR1', 'HUK1_MED1')
SAMPLE_LABELS = {'Ctrl1A2': 'Mouse 1', 'Ctrl1A4': 'Mouse 2',
                 'HUK1_COR1': 'Human cortex 1', 'HUK1_MED1': 'Human cortex 2'}
PATHWAYS = (
    ('Reactome_2022', 'Citric Acid Cycle (TCA Cycle) R-HSA-71403', 'TCA cycle'),
    ('KEGG_2019_Mouse', 'Oxidative phosphorylation', 'Oxidative phosphorylation'),
    ('Reactome_2022', 'Peroxisomal Lipid Metabolism R-HSA-390918', 'Peroxisomal lipid metabolism'),
)
GENES = ('Mme', 'Slc27a2')
LAMBDA_GRID = np.logspace(-3, 3, 13)
FOOTER = 'Descriptive cohort: 2 mouse specimens; 2 human cortex slices from 1 donor. No confirmatory species test.'

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['DejaVu Sans'],
    'font.size': 11, 'axes.titlesize': 12, 'axes.labelsize': 11,
    'xtick.labelsize': 10, 'ytick.labelsize': 10,
    'legend.fontsize': 10, 'legend.frameon': False,
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.linewidth': 0.7, 'axes.edgecolor': '#65717C',
    'svg.fonttype': 'none', 'pdf.fonttype': 42,
    'savefig.facecolor': 'white', 'figure.facecolor': 'white',
})


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results-root', type=Path, default=Path(
        os.environ.get('PSEUDOSPACE_RESULTS_ROOT', PROJECT / 'results')))
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--qa-tools', type=Path, help='Directory containing audit_panel_alignment.py')
    return parser.parse_args()


def csv(root, name):
    return pd.read_csv(root / 'curves' / name)


def load_inputs(root):
    """Match saved fits to current tables, then reconstruct specimen summaries."""
    print('Reading current PT object and result tables.', flush=True)
    a = ad.read_h5ad(root / 'cross_species_pt_dpt.h5ad')
    keep = (a.obs['broad_tubule_marker_call'].astype(str).eq('PT')
            & a.obs['comparison_species'].astype(str).isin(('mouse', 'human'))
            & np.isfinite(a.obs['total_scanpy_dpt'].to_numpy(dtype=float)))
    a = a[keep.to_numpy()].copy()
    results = csv(root, 'gene_trajectory_comparison.csv').set_index('gene')
    names = a.var_names[a.var_names.isin(results.index)].to_numpy().astype(str)
    if len(names) != len(results):
        raise ValueError('Saved gene results do not match the current PT object.')
    if not a.var.loc[names, 'measured_in_both_inputs'].astype(bool).all():
        raise ValueError('A result gene is not measured in every input.')
    results = results.loc[names]
    expected = results['level_effect_human_minus_mouse'].to_numpy()
    fit = None
    cache_path = None
    for path in sorted((root / 'stage_cache').glob('section4_gene_fit__*.npz')):
        with np.load(path, allow_pickle=False) as z:
            if z['level_effect'].shape != expected.shape:
                continue
            if not np.allclose(z['level_effect'], expected, rtol=1e-7, atol=1e-8):
                continue
            if not np.allclose(z['shape_rms'], results['shape_rms'], rtol=1e-7, atol=1e-8):
                continue
            fit = {k: z[k] for k in ('lam_idx', 'curve_healthy', 'curve_aki')}
            cache_path = path
            break
    if fit is None:
        raise ValueError('No cached gene fit agrees with the complete current gene table.')
    Y = a[:, names].layers['lognorm'].tocsr().astype(np.float64)
    s = a.obs['total_scanpy_dpt'].to_numpy(dtype=float)
    species = a.obs['comparison_species'].astype(str).to_numpy()
    samples = a.obs['sample'].astype(str).to_numpy()
    if set(samples) != set(SAMPLES):
        raise ValueError('Unexpected specimen cohort; review the figure contract before exporting.')
    lo = max(np.percentile(s[species == k], 1) for k in COLORS)
    hi = min(np.percentile(s[species == k], 99) for k in COLORS)
    grid = np.linspace(lo, hi, fit['curve_healthy'].shape[1])
    knots = gam_internal_knots(s, basis_df=12)
    print(f'Reconstructing four specimen curves for {len(names):,} tested genes.', flush=True)
    specimens = {}
    for sample in SAMPLES:
        mask = samples == sample
        specimens[sample] = fit_single_condition_curves(
            Y[mask], s[mask], knots, grid, LAMBDA_GRID, fit['lam_idx'], support_pct=(1, 99))[0]
    balanced = {
        'mouse': specimen_balanced_curves({k: specimens[k] for k in SAMPLES[:2]}),
        'human': specimen_balanced_curves({k: specimens[k] for k in SAMPLES[2:]}),
    }
    saved_balanced = csv(root, 'specimen_balanced_gene_amplitudes.csv').set_index('gene').loc[names]
    np.testing.assert_allclose(
        np.nanmean(balanced['human'] - balanced['mouse'], axis=1),
        saved_balanced['balanced_level_effect_human_minus_mouse'], atol=1e-7, rtol=1e-7)

    pathway_table = csv(root, 'pathway_trajectory_comparison.csv')
    selection, pathway_scores = [], []
    mean = np.asarray(Y.mean(axis=0)).ravel()
    std = np.sqrt(np.maximum(np.asarray(Y.multiply(Y).mean(axis=0)).ravel() - mean ** 2, 1e-12))
    lookup = {g: i for i, g in enumerate(names)}
    for library, pathway, label in PATHWAYS:
        rows = pathway_table.loc[pathway_table.library.eq(library) & pathway_table.pathway.str.strip().eq(pathway)]
        if len(rows) != 1:
            raise ValueError(f'Expected one pathway entry: {library}: {pathway}')
        row = rows.iloc[0]
        members = [g.strip() for g in row.genes_present.split(';')]
        idx = np.array([lookup[g] for g in members])
        score = ((Y[:, idx].toarray() - mean[idx]) / std[idx]).mean(axis=1)
        pathway_scores.append(score)
        selection.append({'library': library, 'pathway': row.pathway,
                          'label': label, 'members': members, 'row': row})
    scores = np.column_stack(pathway_scores)
    pathway_fit = run_level_shape(scores, s, (species == 'human').astype(float), knots, grid, LAMBDA_GRID)
    saved_paths = csv(root, 'top_pathway_fitted_curves.csv')
    for i, entry in enumerate(selection):
        for kind, key in [('mouse', 'curve_healthy'), ('human', 'curve_aki')]:
            old = saved_paths.loc[saved_paths.library.eq(entry['library'])
                                  & saved_paths.pathway.eq(entry['pathway'])
                                  & saved_paths.species.eq(kind)].sort_values('pseudospace')
            if len(old) != len(grid):
                raise ValueError('Selected pathway is missing from saved curve outputs.')
            np.testing.assert_allclose(pathway_fit[key][i], old.fitted_module_z_score, rtol=1e-7, atol=1e-7)
    path_specimens = {}
    for sample in SAMPLES:
        mask = samples == sample
        path_specimens[sample] = fit_single_condition_curves(
            scores[mask], s[mask], knots, grid, LAMBDA_GRID, pathway_fit['lam_idx'])[0]
    print('Saved gene effects, balanced effects, and selected pathway curves verified.', flush=True)
    return SimpleNamespace(root=root, a=a, results=results, names=names, lookup=lookup, grid=grid,
                           fit=fit, specimens=specimens, balanced=balanced, pathways=selection,
                           pathway_fit=pathway_fit, path_specimens=path_specimens,
                           modules=csv(root, 'module_gene_assignment.csv').set_index('gene').loc[names],
                           enrichment=csv(root, 'module_pathway_enrichment.csv'), cache_path=cache_path)


def canvas(title, subtitle, rows=1, cols=1):
    fig, axes = plt.subplots(rows, cols, figsize=(13.33, 8.3), squeeze=False)
    fig.subplots_adjust(left=0.10, right=0.975, bottom=0.22, top=0.75, hspace=0.88, wspace=0.48)
    fig.text(0.06, 0.955, title, fontsize=21, weight='bold', va='top', color='#152B3C')
    fig.text(0.06, 0.90, subtitle, fontsize=12, va='top', color='#4B5965')
    fig.text(0.06, 0.035, FOOTER, fontsize=10, color='#4B5965')
    for i, ax in enumerate(axes.flat):
        ax.annotate(chr(97 + i), (0, 1), xycoords='axes fraction', xytext=(-34, 17),
                    textcoords='offset points', fontsize=14, weight='bold', annotation_clip=False)
        ax.set_axisbelow(True)
    return fig, axes


def trajectory_axis(ax, ylabel):
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.5, 1], ['Early PT', 'Middle', 'Late PT'])
    ax.set_xlabel('Position along inferred PT coordinate', labelpad=8)
    ax.set_ylabel(ylabel, labelpad=9)
    ax.grid(axis='y', color='#E7EBEE', linewidth=0.7)


def species_legend(fig, balanced=False):
    handles = [Line2D([], [], color=COLORS[k], lw=2.7, label=('Mouse' if k == 'mouse' else 'Human'))
               for k in COLORS]
    handles += [Line2D([], [], color='#6B747C', lw=1.2, ls=style, label=label)
                for style, label in [('--', 'Individual specimen / slice 1'), (':', 'Individual specimen / slice 2')]]
    fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(0.055, 0.854), ncol=4, columnspacing=2)


def save_sources(out, name, rows):
    frame = rows if isinstance(rows, pd.DataFrame) else pd.DataFrame(rows)
    frame.to_csv(out / 'source_data' / name, index=False)


def metabolic_figure(d, out):
    fig, axes = canvas('Metabolic gene sets show lower relative expression in human',
                       'Three selected metabolic themes; curves show expression and the lower panels show every member gene.',
                       rows=2, cols=3)
    species_legend(fig)
    curve_rows, effect_rows = [], []
    for j, entry in enumerate(d.pathways):
        ax, bx = axes[:, j]
        ax.set_title(entry['label'], pad=14)
        for n, sample in enumerate(SAMPLES):
            kind = 'mouse' if n < 2 else 'human'
            y = d.path_specimens[sample][j]
            ax.plot(d.grid, y, color=COLORS[kind], lw=1.1, alpha=0.65, ls='--' if n % 2 == 0 else ':')
            curve_rows.extend(dict(pathway=entry['label'], group=sample, species=kind,
                                   estimator='individual specimen fit', position=x, value=v)
                              for x, v in zip(d.grid, y))
        for kind, key in [('mouse', 'curve_healthy'), ('human', 'curve_aki')]:
            y = d.pathway_fit[key][j]
            ax.plot(d.grid, y, color=COLORS[kind], lw=2.6)
            curve_rows.extend(dict(pathway=entry['label'], group=kind, species=kind,
                                   estimator='pooled GAM', position=x, value=v) for x, v in zip(d.grid, y))
        trajectory_axis(ax, 'Relative pathway expression\n(mean gene z-score)' if j == 0 else '')
        effects = d.results.loc[entry['members'], 'level_effect_human_minus_mouse'].sort_values()
        lower = int((effects < 0).sum())
        # Every available member is plotted; vertical order is effect rank, not pseudospace.
        ys = np.arange(len(effects))
        bx.axvline(0, color='#8F969D', lw=0.9, ls='--')
        bx.scatter(effects, ys, c=[COLORS['mouse'] if x < 0 else COLORS['human'] for x in effects],
                   s=17, edgecolors='none', zorder=3)
        bx.set_title(f'{lower} / {len(effects)} members lower in human', pad=10)
        bx.set_yticks([])
        bx.set_ylabel('Each dot = one member gene' if j == 0 else '')
        bx.set_xlabel('Human − mouse level difference\n(log-normalized expression)', labelpad=8)
        bx.grid(axis='x', color='#E7EBEE', lw=0.7)
        for gene, effect in effects.items():
            effect_rows.append(dict(pathway=entry['label'], gene=gene, human_minus_mouse=effect))
    # Identical expression-score limits and gene-effect limits across the thematic panels.
    ymax = max(max(abs(v) for v in ax.get_ylim()) for ax in axes[0]) * 1.03
    for ax in axes[0]: ax.set_ylim(-ymax, ymax)
    xmin = min(ax.get_xlim()[0] for ax in axes[1]); xmax = max(ax.get_xlim()[1] for ax in axes[1])
    for ax in axes[1]: ax.set_xlim(xmin, xmax)
    fig.text(0.06, 0.093, 'Solid curves: pooled GAMs. Dashed/dotted: individual inputs. Member dots are effect sizes, not independent replicates.', fontsize=10)
    fig.text(0.06, 0.066, 'Interpretation: coordinated relative-expression offsets; these measurements do not establish lower absolute RNA or metabolic activity.', fontsize=10)
    save_sources(out, '01_pathway_curves.csv', curve_rows)
    save_sources(out, '01_member_gene_effects.csv', effect_rows)
    return fig


def gradient_figure(d, out):
    fig, axes = canvas('Overall expression gaps coexist with opposite PT gradients',
                       'Two illustrative genes identified in the current analysis. Mean-centering removes the offset while preserving gradient strength.',
                       rows=2, cols=2)
    species_legend(fig, balanced=True)
    rows, gradients = [], []
    early, late = d.grid <= np.quantile(d.grid, 1/3), d.grid >= np.quantile(d.grid, 2/3)
    for j, gene in enumerate(GENES):
        idx = d.lookup[gene]
        for r, centered in enumerate((False, True)):
            ax = axes[r, j]
            ax.set_title(gene + (' — expression level' if not centered else ' — spatial gradient'), pad=14)
            for n, sample in enumerate(SAMPLES):
                kind = 'mouse' if n < 2 else 'human'
                y = d.specimens[sample][idx].copy()
                if centered: y -= np.nanmean(y)
                ax.plot(d.grid, y, color=COLORS[kind], ls='--' if n % 2 == 0 else ':', lw=1.1, alpha=0.65)
                rows.extend(dict(gene=gene, group=sample, centered=centered, position=x, expression=v)
                            for x, v in zip(d.grid, y))
            for kind in COLORS:
                y = d.balanced[kind][idx].copy()
                if centered: y -= np.nanmean(y)
                ax.plot(d.grid, y, color=COLORS[kind], lw=2.7)
                rows.extend(dict(gene=gene, group=kind, centered=centered, position=x, expression=v)
                            for x, v in zip(d.grid, y))
            if centered: ax.axhline(0, color='#7E8790', lw=0.75, ls='--')
            trajectory_axis(ax, ('Expression\n(log-normalized)' if not centered else 'Centered expression\n(log-normalized)'))
        for sample in SAMPLES:
            y = d.specimens[sample][idx]
            gradients.append(dict(gene=gene, sample=sample, late_minus_early=float(np.nanmean(y[late]) - np.nanmean(y[early]))))
    extent = max(max(abs(v) for v in ax.get_ylim()) for ax in axes[1])
    for ax in axes[1]: ax.set_ylim(-extent, extent)
    fig.text(0.06, 0.093, 'Solid curves: equal-weight mean of the two available specimen curves per species. Thin lines expose specimen variation; no CI is claimed.', fontsize=10)
    fig.text(0.06, 0.066, 'Both examples rise toward late PT in mouse and fall in these human slices. Position is inferred expression order, not measured anatomical distance.', fontsize=10)
    save_sources(out, '02_gene_curves.csv', rows)
    save_sources(out, '02_individual_gradients.csv', gradients)
    return fig


def module_figure(d, out):
    fig, axes = canvas('Difference-curve modules connect spatial patterns to pathways',
                       'Modules group genes by the whole human − mouse difference profile, after removing each gene’s mean difference.',
                       rows=2, cols=2)
    fig.subplots_adjust(left=0.18, wspace=0.7)
    seeds = ('TNF-alpha Signaling via NF-kB', 'Oxidative Phosphorylation')
    titles = ('Late-rising difference module', 'Metabolic-associated difference module')
    terms = ((('TNF-alpha Signaling via NF-kB', 'TNF / NF-κB'), ('Hypoxia', 'Hypoxia'),
              ('IL-2/STAT5 Signaling', 'IL-2 / STAT5')),
             (('Oxidative Phosphorylation', 'Oxidative phosphorylation'), ('Citric Acid Cycle (TCA Cycle)', 'TCA cycle'),
              ('Peroxisomal Protein Import', 'Peroxisomal import')))
    delta = d.balanced['human'] - d.balanced['mouse']
    centered = delta - np.nanmean(delta, axis=1, keepdims=True)
    spread = np.nanstd(centered, axis=1, keepdims=True)
    standardized = np.divide(centered, spread, out=np.full_like(centered, np.nan), where=spread > 0)
    profile_rows, set_rows, membership_rows = [], [], []
    colors = ('#8064A2', '#398477')
    for j, seed in enumerate(seeds):
        hits = d.enrichment[d.enrichment.gene_set.str.contains(seed, regex=False)]
        module = hits.sort_values('p_value_adjusted').iloc[0].module
        mask = d.modules.response_module.eq(module).to_numpy()
        profiles = standardized[mask]
        q25, median, q75 = np.nanpercentile(profiles, [25, 50, 75], axis=0)
        ax, bx = axes[:, j]
        ax.fill_between(d.grid, q25, q75, color=colors[j], alpha=0.18, linewidth=0)
        ax.plot(d.grid, median, color=colors[j], lw=2.8)
        ax.axhline(0, color='#7E8790', lw=0.8, ls='--')
        ax.set_title(f'{titles[j]}\n{module}: {mask.sum():,} genes', pad=12)
        trajectory_axis(ax, 'Difference profile\n(own SD units)')
        for x, m, low, high in zip(d.grid, median, q25, q75):
            profile_rows.append(dict(module=module, position=x, median=m, q25=low, q75=high, n_genes=int(mask.sum())))
        membership_rows.extend(dict(gene=g, module=module) for g in d.names[mask])
        module_sets = d.enrichment[d.enrichment.module.eq(module)]
        plotted = []
        for term, label in terms[j]:
            candidates = module_sets[module_sets.gene_set.str.contains(term, regex=False)]
            row = candidates.sort_values('p_value_adjusted').iloc[0]
            plotted.append((row, label))
            set_rows.append(dict(row.to_dict(), display_label=label))
        maximum = max(max(r.n_overlap, r.expected_overlap) for r, _ in plotted)
        for y, (row, label) in enumerate(plotted):
            bx.plot([row.expected_overlap, row.n_overlap], [y, y], color='#AAB3BB', lw=1.3)
            bx.scatter(row.expected_overlap, y, facecolors='white', edgecolors='#4B5965', s=55, linewidths=1.2, zorder=3)
            bx.scatter(row.n_overlap, y, color=colors[j], s=55, zorder=4)
            bx.text(maximum * 1.20, y, f'q = {row.p_value_adjusted:.2g}', va='center', fontsize=10)
        bx.set_yticks(range(len(plotted)), [label for _, label in plotted])
        bx.set_ylim(len(plotted) - 0.5, -0.5)
        bx.set_xlim(0, maximum * 1.95)
        bx.set_xticks([0, 10, 20, 30, 40] if maximum > 40 else [0, 10, 20])
        bx.set_xlabel('Number of pathway members in module', labelpad=8)
        bx.set_title('Observed versus expected pathway membership', fontsize=10.5, pad=12)
        bx.grid(axis='x', color='#E7EBEE', lw=0.7)
    extent = max(max(abs(v) for v in ax.get_ylim()) for ax in axes[0])
    for ax in axes[0]: ax.set_ylim(-extent, extent)
    handles = [Line2D([], [], color='#596A7A', lw=2.5, label='Median shape; band = middle 50% of genes'),
               Line2D([], [], marker='o', linestyle='none', color='#4B5965', label='Observed membership'),
               Line2D([], [], marker='o', linestyle='none', color='#4B5965', markerfacecolor='white', label='Expected membership')]
    fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(0.055, 0.854), ncol=3)
    fig.text(0.06, 0.093, 'The band shows heterogeneity across genes, not uncertainty across donors. A rising difference can reflect a smaller human deficit toward late PT.', fontsize=10)
    fig.text(0.06, 0.066, f'Exploratory hypergeometric enrichment; BH across all {len(d.enrichment):,} module–pathway pairs. This does not test pathway activation.', fontsize=10)
    save_sources(out, '03_module_profiles.csv', profile_rows)
    save_sources(out, '03_module_enrichment.csv', set_rows)
    save_sources(out, '03_module_membership.csv', membership_rows)
    return fig


def robustness_figure(d, out):
    fig, axes = canvas('Level differences are stable; exact module boundaries are less stable',
                       'Sensitivity analyses ask what changes when we alter weighting or omit an input. They are not biological replication.',
                       rows=2, cols=2)
    fig.subplots_adjust(left=0.20, wspace=0.75)
    top = d.results[d.results.primary_eligible].nlargest(20, 'species_effect_rms')
    fraction = float(top.level_fraction.median())
    ax = axes[0, 0]
    values = np.sort(top.level_fraction.to_numpy())
    ax.scatter(values, np.arange(len(values)), s=25, color='#526F89')
    ax.set_xlim(0, 1.03); ax.set_yticks([])
    ax.set_xticks([0, .5, 1], ['0%', '50%', '100%'])
    ax.set_title(f'Top 20 genes\nMedian {fraction:.1%} level contribution', pad=14)
    ax.set_xlabel('Fraction of squared curve separation due to level')
    ax.set_ylabel('Ranked genes\n(one dot each)')
    agreement = csv(d.root, 'signed_ranking_pooled_vs_balanced.csv')
    ax = axes[0, 1]
    labels = ['Level', 'Amplitude', 'Redistribution', 'Early contrast', 'Late contrast']
    vals = agreement.spearman_pooled_vs_balanced.to_numpy()
    ax.scatter(vals, range(len(vals)), color='#398477', s=35)
    for i, v in enumerate(vals): ax.text(1.055, i, f'{v:.4f}', va='center', fontsize=10)
    ax.set_yticks(range(len(vals)), labels); ax.set_ylim(len(vals) - .5, -.5)
    ax.set_xlim(0, 1.28); ax.set_xticks([0, .5, 1])
    ax.set_title('Pooled versus equal-specimen weighting\nRank agreement', pad=14)
    ax.set_xlabel('Spearman rank correlation')
    coordinate = csv(d.root, 'coordinate_specimen_robustness.csv')
    stability = csv(d.root, 'module_response_stability.csv')
    stability = stability[stability.run.str.startswith('without_')].copy()
    for j, (frame, value, label) in enumerate([
        (coordinate, 'spearman', 'Coordinate ordering\nAfter specimen omission'),
        (stability, 'adjusted_rand_index_vs_reference', 'Module membership\nAfter specimen omission'),
    ]):
        ax = axes[1, j]
        sample = frame.dropped_sample if j == 0 else frame.run.str.replace('without_', '', regex=False)
        order = pd.Index(sample).get_indexer(SAMPLES)
        if (order < 0).any(): raise ValueError('Missing leave-one-specimen-out result.')
        vals = frame.iloc[order][value].to_numpy()
        ax.scatter(vals, range(4), s=35, color='#526F89' if j == 0 else '#8064A2')
        for i, v in enumerate(vals): ax.text(1.055, i, f'{v:.3f}', va='center', fontsize=10)
        ax.set_yticks(range(4), [f'Omit {SAMPLE_LABELS[s]}' for s in SAMPLES]); ax.set_ylim(3.5, -.5)
        ax.set_xlim(0, 1.28); ax.set_xticks([0, .5, 1]); ax.set_title(label, pad=14)
        ax.set_xlabel('Spearman rank correlation' if j == 0 else 'Adjusted Rand index (1 = identical partition)')
    for ax in axes.flat: ax.grid(axis='x', color='#E7EBEE', lw=0.7)
    normal = csv(d.root, 'normalisation_sensitivity.csv').iloc[0]
    fig.text(0.06, 0.093, f'Normalization check: level-effect correlation = {normal.level_effect_spearman_panel_vs_native:.6f}; top-20 sign agreement = {normal.top20_sign_agreement:.0%}.', fontsize=10)
    fig.text(0.06, 0.066, 'Coordinate checks reuse Harmony and an anchored root. Assay/capture differences and the single human donor remain unresolved.', fontsize=10)
    save_sources(out, '04_top_gene_level_fractions.csv', top.reset_index())
    save_sources(out, '04_weighting_agreement.csv', agreement)
    save_sources(out, '04_coordinate_omission.csv', coordinate)
    save_sources(out, '04_module_omission.csv', stability)
    save_sources(out, '04_normalization.csv', normal.to_frame().T)
    return fig


def text_page(title, subtitle, sections):
    fig = plt.figure(figsize=(13.33, 8.3))
    fig.text(0.06, 0.95, title, fontsize=23, weight='bold', va='top', color='#152B3C')
    fig.text(0.06, 0.885, subtitle, fontsize=13, va='top', color='#4B5965')
    y = 0.80
    for heading, body in sections:
        fig.text(0.06, y, heading, fontsize=15, weight='bold', va='top', color='#0072B2')
        y -= 0.045
        lines = '\n'.join(textwrap.fill(paragraph, width=112) for paragraph in body.split('\n'))
        fig.text(0.06, y, lines, fontsize=12, va='top', linespacing=1.5)
        y -= 0.030 * len(lines.splitlines()) + 0.026
    if y < .045:
        raise ValueError('Methods page would overflow; shorten or split the content.')
    fig.text(0.06, .035, FOOTER, fontsize=10, color='#4B5965')
    return fig


def score_page(d):
    return text_page('What does “fitted module z-score” mean?',
        'Clearer label: Relative pathway expression (mean gene z-score)', [
        ('1   Normalize each tubule’s gene counts',
         'Divide by the total counts over the retained gene panel, multiply by 10,000, then apply ln(1 + value). '
         'The input is relative expression, not absolute RNA abundance.'),
        ('2   Put each gene on its own reference scale',
         'For each gene: z = (log-normalized expression − pooled PT mean) / pooled PT standard deviation. '
         'Pool mouse and human structures together; do not standardize the species separately.'),
        ('3   Average the available member genes',
         'For each tubule, the pathway score is the arithmetic mean of its member-gene z-scores. '
         'Illustrative arithmetic: gene scores of +0.5, +1.0 and −0.5 give a pathway score of +0.33.'),
        ('4   Fit a smooth curve along PT',
         'A generalized additive model (GAM) estimates how that score varies with PT position in each species. '
         'The y-axis is this fitted expectation. The pathway score is not standardized again.'),
        ('How to read the number',
         '+1 means that members average one gene-specific SD above their pooled means. Zero is the pooled reference, '
         'not zero expression. The score is not a p-value, fold change, pathway SD or direct activity measurement.'),
    ])


def tests_page(d):
    return text_page('What was calculated, and what was actually tested?',
        'Three different kinds of evidence; their statistical units must remain distinct.', [
        ('GAM curves and effect sizes — descriptive',
         'Fit a shared smooth, then a species offset, then a species-by-position smooth. Choose smoothing by '
         'generalized cross-validation on the full model. Report level, amplitude and pattern differences. '
         'The cellwise F statistics are uncalibrated for species inference and are not used as species p-values.'),
        ('Module membership enrichment — exploratory overlap test',
         'Cluster centered, specimen-balanced human − mouse curves by whole-profile correlation. Ask whether '
         'a pathway contributes more members to a module than expected from eligible genes. Use a hypergeometric '
         f'tail and BH correction over all {len(d.enrichment):,} module–pathway pairs.'),
        ('Signed-ranking enrichment — exploratory competitive test',
         'For each question, compare the mean gene statistic in a pathway with the tested-gene background mean. '
         'Use the custom normal approximation with finite-population and squared-residual-correlation factors; '
         'BH correction spans five questions under both pooled and balanced weighting. This is not validated CAMERA.'),
        ('Robustness checks — sensitivity, not added replication',
         'Change the normalization denominator; weight specimens equally; omit one specimen and compare results. '
         'The two human slices share one donor. Neither many tubules nor many genes creates extra human donors.'),
    ])


def guide_html(d, out):
    cohort = csv(d.root, 'pt_cohort_summary.csv')
    enrich = csv(d.root, 'signed_enrichment_summary.csv')
    significant = d.enrichment[d.enrichment.p_value_adjusted < .05].copy()
    save_sources(out, 'all_significant_module_pathways.csv', significant)
    significant_display = significant[[
        'module', 'gene_set', 'n_overlap', 'expected_overlap', 'p_value_adjusted',
    ]].rename(columns={'module': 'Module', 'gene_set': 'Pathway', 'n_overlap': 'Observed members',
                       'expected_overlap': 'Expected members', 'p_value_adjusted': 'BH q'})
    significant_display['Expected members'] = significant_display['Expected members'].map(lambda v: f'{v:.1f}')
    significant_display['BH q'] = significant_display['BH q'].map(lambda v: f'{v:.3g}')
    enrich_display = enrich.rename(columns={
        'question': 'Question', 'n_genes_ranked': 'Genes ranked', 'n_sets_tested': 'Sets tested',
        'n_significant_against_background_after_BH': 'Sets with BH q < 0.05',
        'median_variance_inflation': 'Median variance factor', 'correlation_basis': 'Correlation basis',
    })
    descriptions = [
        ('01_metabolic_programs', 'Metabolic relative-expression differences',
         'The selected TCA, oxidative-phosphorylation and peroxisomal-lipid sets show lower human scores. '
         'The member panels display every tested member, so the collective direction can be assessed without relying on a pathway label. '
         'Thin curves expose specimen variation; they are not confidence intervals. These are selected metabolic themes, not independent discoveries.'),
        ('02_gene_gradients', 'Expression level and spatial gradient are different questions',
         'Mme and Slc27a2 are illustrative gradient candidates from the current analysis. Top panels retain expression level. '
         'Bottom panels subtract each curve’s mean, retaining gradient magnitude in log-normalized units. '
         'Equal-weight specimen means are thick; both individual inputs remain visible. These examples were selected after inspecting the data.'),
        ('03_response_modules', 'Pathway associations of shared difference shapes',
         'Profiles are the human-minus-mouse specimen-balanced difference, centered and divided by each difference curve’s own SD. '
         'The line is the median across all assigned genes; the band is the interquartile range across genes, not a donor confidence interval. '
         'The lower panels compare observed and expected membership for selected pathway annotations. '
         'A late-rising difference does not necessarily mean increasing absolute expression or activation in human.'),
        ('04_robustness', 'Robust conclusions and sensitive boundaries',
         'Level effects and signed rankings are stable under the checks shown. Module partitions change more when one input is omitted. '
         'Coordinate omission reuses the existing Harmony representation and anchored root. None of these checks removes assay or donor confounding.'),
    ]
    figure_sections = ''.join(f'<section><h2>{i+1}. {title}</h2><img src="{name}.png" alt="{html.escape(title)}">'
                              f'<p>{text}</p><p><a href="{name}.pdf">PDF</a> · <a href="{name}.svg">Editable SVG</a></p></section>'
                              for i, (name, title, text) in enumerate(descriptions))
    content = f'''<!doctype html><html lang="en"><meta charset="utf-8"><title>Human–mouse PT: collaborator brief</title>
<style>body{{font:17px/1.65 system-ui,sans-serif;color:#243443;background:#f3f6f8;margin:0}}
main{{max-width:1120px;margin:auto;padding:45px 35px;background:white}}h1,h2,h3{{line-height:1.25;color:#152b3c}}
h1{{font-size:38px}}h2{{margin-top:40px}}.lead{{font-size:21px}}.note{{padding:18px 22px;background:#eef4f8;border-left:4px solid #0072b2}}
img{{width:100%;height:auto}}table{{border-collapse:collapse;width:100%;font-size:15px}}td,th{{padding:9px;border-bottom:1px solid #dce4e9;text-align:left}}
code,pre{{background:#f3f6f8;font-size:15px}}pre{{padding:18px;white-space:pre-wrap}}a{{color:#006c9f}}section{{margin-top:42px}}li{{margin:9px 0}}th,td{{overflow-wrap:anywhere}}
@media print{{body{{background:white}}main{{padding:0}}section{{break-before:page}}}}</style><main>
<h1>Human–mouse proximal tubule expression</h1><p class="lead">A figure-led explanation of the strongest descriptive findings, how they were calculated, and what the tests mean.</p>
<p><a href="collaborator_brief.pdf">Download the six-page PDF brief</a> · <a href="source_data/">Figure source data</a></p>
<div class="note"><b>Design:</b> {len(d.a):,} PT structures and {len(d.names):,} tested genes; two mouse specimens versus one human donor represented by two cortex slices.
Each observation is a segmented tubule aggregate, not a single cell.
The inferred PT coordinate is specific to notebook 03. These results generate hypotheses; they do not establish population-level species differences.</div>
<h2>Read the pathway y-axis correctly</h2><p><b>“Fitted module z-score” is the GAM-smoothed mean of standardized member-gene expression.</b>
Here “module” means a predefined pathway gene set, not a discovered curve cluster. We use <b>Relative pathway expression (mean gene z-score)</b> in the figures.</p>
<pre>x[i,g] = ln(1 + 10,000 × counts[i,g] / retained-panel counts[i])
z[i,g] = (x[i,g] − mean_PT(x[g])) / SD_PT(x[g])
pathway_score[i,P] = mean of z[i,g] over tested members g in pathway P
plotted_value[P,species,s] = GAM-fitted pathway_score at PT position s</pre>
<p>The gene mean and population SD use all pooled mouse and human PT structures (so the reference is tubule-weighted).
Each gene receives equal weight after standardization. The pathway score is not standardized to unit variance afterwards.
For example, scores +0.5, +1.0 and −0.5 average to +0.33; this arithmetic example is illustrative, not a measured tubule.</p>
<ul><li><b>Zero:</b> the member genes average their pooled reference level; not no expression.</li>
<li><b>Positive:</b> members average above their own pooled means; negative means below.</li>
<li><b>Not encoded:</b> statistical significance, fold change, absolute RNA amount, or pathway activity.</li>
<li><b>Why standardize:</b> highly expressed genes do not dominate solely because of baseline abundance. Conversely, small or noisy changes can receive large relative weight when a gene’s SD is small.</li></ul>
{figure_sections}
<section><h2>Calculation and testing walkthrough</h2>
<h3>1. Define comparable measurements</h3><p>Human features are mapped to mouse symbols with an accepted reciprocal HCOP mapping, requiring at least three supporting databases.
The analysis gene set additionally requires a source feature in every input. This matters because a zero-filled ortholog target is not a measured zero.
After the earlier global expression/count filter (at least 5% of structures and 20 counts), the PT analysis requires detection in at least 2% of PT structures.
The current table has {len(d.names):,} tested genes. The broader integration object is retained to preserve reviewed clustering; the analysis gate is downstream.</p>
<h3>2. Normalize counts and define the coordinate</h3><p>Counts are normalized to 10,000 over the retained gene panel, followed by ln(1 + x).
Harmony is used for the representation underlying the trajectory, while expression models use the log-normalized count layer.
Section 4 uses the PT-specific recomputed DPT, not physical distance and not notebook 02’s coordinate. Curves are compared on the overlap of each species’ 1st–99th percentile support,
using {len(d.grid)} grid points. The fits use the eligible structures; the support restriction controls where predictions are compared.</p>
<h3>3. Fit smooth expression profiles</h3><pre>M0: expression = intercept + shared smooth(position)
M1: expression = intercept + shared smooth(position) + species offset
M2: expression = intercept + shared smooth(position) + species offset
                + species-specific change in smooth(position)</pre>
<p>These are Gaussian penalized cubic-spline fits on log-normalized gene expression, or on mean-z pathway scores.
Nine internal knots are used. For each response, generalized cross-validation selects one smoothing penalty from 13 candidates between 0.001 and 1,000 on M2;
the same penalty is reused for the nested models. GCV chooses smoothness, not biological significance.</p>
<h3>4. Decompose the fitted difference</h3><pre>delta(s) = human_curve(s) − mouse_curve(s)
level = mean_s delta(s)
shape_RMS = sqrt(mean_s [delta(s) − level]^2)
total_RMS^2 = level^2 + shape_RMS^2
level_fraction = level^2 / total_RMS^2
amplitude = maximum(curve) − minimum(curve)
pattern_RMS_z = RMS of the difference after standardizing each species curve separately</pre>
<p>Grid averages give equal weight to coordinate positions rather than tubule density. “Shape RMS” still mixes amplitude and profile-pattern changes.
The additional standardized pattern comparison separates profile geometry from level and amplitude.
Pattern status requires both curve amplitudes to be at least 0.15 and a standardized difference of at least 0.25 to be called “supported”.
These thresholds are reading rules, not significance tests. A level contribution of at least 50% determines the dominant level label; total effects below 0.05 are labeled “none”.</p>
<h3>5. Put specimens back into the summary</h3><p>Each specimen receives its own single-condition smooth, using the pooled fit’s selected penalty.
It is masked outside that specimen’s 1st–99th percentile support. Species summaries average the available specimen curves equally at each grid position.
At an endpoint supported by only one specimen, only that specimen contributes. The exported individual curves make this visible.
Primary response-module discovery uses these balanced curves; pooled and balanced signed rankings are both reported. The pathway curves in Figure 1 retain the notebook’s pooled estimator.</p>
<h3>6. Discover groups of similar curves</h3><p>Positional modules cluster balanced mouse curves. Response modules cluster balanced human-minus-mouse difference curves after mean-centering.
Curves are standardized over supported grid positions; average-linkage hierarchical clustering uses 1 − Pearson correlation.
The cut distance is 0.4, minimum curve amplitude is 0.05, and clusters smaller than 10 genes are unassigned.
Grouping uses the whole profile; peaks, widths and monotonicity describe the result. A constant offset is removed, but a gene with a large offset can still have a varying difference curve and enter a module.</p>
<h3>7. Distinguish the two enrichment tests</h3>
<table><tr><th>Analysis</th><th>Question and calculation</th><th>Interpretation</th></tr>
<tr><td>Module overlap</td><td>Within the discovery-eligible universe of N genes, a pathway has K members and a module n genes.
Expected overlap = nK/N. The hypergeometric upper-tail probability asks how unusual the observed overlap is.
BH correction covers all {len(d.enrichment):,} tested module–pathway pairs.</td><td>Exploratory annotation of a discovered module; not a test of activity or species generalization.</td></tr>
<tr><td>Signed rankings</td><td>For level, amplitude ratio, redistribution, early contrast and late contrast, compare each set’s mean gene statistic against the full tested background mean.
Use both pooled and balanced rankings, with BH over all ten question–weighting combinations.</td><td>Exploratory relative association with a particular expression question. The same underlying data support both weightings.</td></tr></table>
<p><b>Exact current competitive approximation:</b> the standard error is the background statistic SD divided by sqrt(K),
multiplied by sqrt((N − K)/(N − 1)) and sqrt(VIF). VIF = 1 + (m − 1) × mean squared off-diagonal residual correlation,
where m is the number of sampled pathway members. Residuals remove the unpenalized M2 design projection; calculations cap at 200 pathway genes and 2,000 structures.
A two-sided normal tail is used. This custom estimator is not validated CAMERA and is not calibrated for donor-level species inference.
The current notebook does not enable the optional gene-resampling diagnostic, so its reported p-values are not permutation p-values.</p>
<p><b>Redistribution statistic:</b> Spearman correlation between delta(s) and position. It describes direction, not the magnitude of change;
near-flat but nonconstant differences can still have large correlations. Inspect amplitude and individual curves before interpreting this ranking.</p>
<h3>8. Read robustness and significance honestly</h3><p>Gene/pathway cellwise F statistics are labeled uncalibrated and do not establish species significance.
Raw-count pseudobulk bins are exports for later modeling, not a completed count-based differential-expression test.
Normalization sensitivity compares the retained-panel denominator with the pre-filter ortholog-space denominator; it is not absolute-RNA calibration.
Specimen omission examines sensitivity of DPT ordering or module partitions. It does not add donors, and the coordinate check keeps Harmony fixed.</p>
<h3>Current input counts</h3>{cohort.to_html(index=False, border=0, float_format=lambda x: f'{x:.3f}')}
<h3>Exploratory competitive enrichment summary</h3>{enrich_display.to_html(index=False, border=0, float_format=lambda x: f'{x:.3f}')}
<p>Counts are overlapping pathway annotations, not independent mechanisms. Small or numerically zero saved p-values must not be described as zero probability or definitive biology.</p>
<details><summary>All {len(significant)} module–pathway associations with BH q &lt; 0.05</summary>
{significant_display.to_html(index=False, border=0)}</details>
</section><section><h2>How to present this to collaborators</h2>
<ol><li>Start with the one-donor human design and the meaning of the score.</li>
<li>Show the metabolic curves and their member-gene support as relative-expression findings.</li>
<li>Use the two gene examples to separate a vertical offset from a spatial-gradient difference.</li>
<li>Introduce response modules as groups of similar differences, then show the overlap evidence.</li>
<li>End with robust rankings, less stable module boundaries, and the need for independent human donors and assay-aware validation.</li></ol>
<p><b>What would strengthen the conclusions:</b> additional human donors and mouse specimens, validation using an independent measurement or anatomical coordinate,
and sensitivity to assay coverage, normalization, detection and clustering thresholds. Stronger pathway-activity claims require appropriate downstream signatures or functional evidence.</p>
<p>Figures were generated from saved outputs with reconstructed specimen curves checked numerically against the existing results.
No upstream clustering, coordinate or hypothesis test was changed for presentation. CSVs in <code>source_data</code> contain the plotted summaries;
the accompanying exporter and provenance manifest document reconstruction.</p></section></main></html>'''
    (out / 'collaborator_guide.html').write_text(content, encoding='utf-8')
    (out / 'figure_captions.md').write_text('\n\n'.join(f'## {title}\n\n{text}' for _, title, text in descriptions), encoding='utf-8')


def export_figure(fig, out, name, pdf, alignment):
    fig.canvas.draw()
    if alignment is not None and fig.axes:
        alignment(fig, json_out=out / 'qa' / f'{name}.alignment.json',
                  overlay_svg=out / 'qa' / f'{name}.alignment.svg',
                  tolerance_pt=1.5, gutter_tolerance_pt=1.5, strict=True)
    fig.savefig(out / f'{name}.pdf')
    fig.savefig(out / f'{name}.svg')
    fig.savefig(out / f'{name}.png', dpi=300)
    pdf.savefig(fig)
    plt.close(fig)
    print(f'Exported {name}', flush=True)


def main():
    args = parse_args()
    root = args.results_root.expanduser().resolve() / 'human_vs_healthy_mouse'
    out = args.output_dir or root / 'collaborator_figures'
    out.mkdir(parents=True, exist_ok=True)
    (out / 'source_data').mkdir(exist_ok=True)
    (out / 'qa').mkdir(exist_ok=True)
    alignment = None
    if args.qa_tools:
        sys.path.insert(0, str(args.qa_tools.expanduser().resolve()))
        from audit_panel_alignment import require_matplotlib_panel_alignment
        alignment = require_matplotlib_panel_alignment
    d = load_inputs(root)
    with PdfPages(out / 'collaborator_brief.pdf') as pdf:
        export_figure(score_page(d), out, '00_reading_the_score', pdf, alignment)
        for name, maker in [('01_metabolic_programs', metabolic_figure), ('02_gene_gradients', gradient_figure),
                            ('03_response_modules', module_figure), ('04_robustness', robustness_figure)]:
            export_figure(maker(d, out), out, name, pdf, alignment)
        export_figure(tests_page(d), out, '05_calculations_and_tests', pdf, alignment)
    guide_html(d, out)
    provenance = {
        'analysis': 'notebook 03 saved outputs; no upstream changes',
        'n_structures': len(d.a), 'n_tested_genes': len(d.names),
        'gene_fit_cache': d.cache_path.name,
        'gene_results_sha256': hashlib.sha256((root / 'curves/gene_trajectory_comparison.csv').read_bytes()).hexdigest(),
        'module_assignment_sha256': hashlib.sha256((root / 'curves/module_gene_assignment.csv').read_bytes()).hexdigest(),
        'exporter_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'validation': ['all cached gene level/shape effects match current table',
                       'all reconstructed balanced gene level effects match current table',
                       'all selected reconstructed pooled pathway curves match saved output'],
        'selection': {'metabolic_pathways': [p[2] for p in PATHWAYS], 'illustrative_genes': list(GENES),
                      'response_modules': 'best overlap q for TNF-alpha/NF-kB and oxidative phosphorylation'},
        'uncertainty': 'individual specimen curves shown; module bands describe gene IQR, not confidence intervals',
        'export_dimensions_inches': [13.33, 8.3], 'png_dpi': 300,
    }
    (out / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    shutil.copy2(__file__, out / Path(__file__).name)
    print('Figure pack and illustrated calculation guide complete.', flush=True)


if __name__ == '__main__':
    main()
