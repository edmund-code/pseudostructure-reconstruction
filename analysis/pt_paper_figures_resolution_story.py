# %% [markdown]
# # PT resolution story: broad signal and positional detail
#
# This notebook reuses saved results from notebooks 12 and 13. It asks what broad species differences are retained and what additional positional patterns remain visible in reconstructed PT pseudospace.
#
# `T_level` tests a broad species shift, `T_spatial` tests how that difference changes along PT, and `T_total` tests any trajectory difference. A constant species shift can be strong for level and weak for spatial; that is expected.

# %% [markdown]
# ## Load and validate saved results
#
# Set workflow paths and plotting defaults.

# %%
import argparse
import ast
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display

start = Path(__file__).resolve().parent if '__file__' in globals() else Path.cwd()
project = next(path for path in (start, *start.parents) if (path / 'pseudospace').is_dir())
sys.path.insert(0, str(project))
from pseudospace.pathway_remodeling import rank_auc

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--data-root', type=Path)
parser.add_argument('--results-root', type=Path)
args, _ = parser.parse_known_args()
data_root = Path(args.data_root or os.environ.get('PSEUDOSPACE_DATA_ROOT') or project / 'data').expanduser()
results_root = Path(args.results_root or os.environ.get('PSEUDOSPACE_RESULTS_ROOT') or project / 'results').expanduser()
source = results_root / 'pt_pathway_remodeling_scfates'
coordinate_path = results_root / 'minimal_pt_scfates' / 'cross_species_pt_scfates.h5ad'
output = results_root / 'pt_paper_figures_resolution_story'
output.mkdir(parents=True, exist_ok=True)
source_data = output / 'source_data'
source_data.mkdir(exist_ok=True)
ALPHA = .05
MIN_RETENTION = .75  # Retained / ALL planned runs; unavailable runs do not improve this fraction.
segments = ('PT-S1', 'PT-S2', 'PT-S3')
HUMAN_COLOR, MOUSE_COLOR = '#D55E00', '#0072B2'
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['DejaVu Sans'],
    'font.size': 10, 'axes.labelsize': 10, 'axes.titlesize': 11, 'legend.fontsize': 9,
    'axes.spines.top': False, 'axes.spines.right': False, 'lines.linewidth': 1.6})


def save_plot(fig, stem):
    """Save ordinary exploratory plot formats and display the notebook PNG."""
    for suffix in ('pdf', 'svg', 'png'):
        fig.savefig(output / f'exploratory_{stem}.{suffix}', dpi=200 if suffix == 'png' else None)
    plt.show()
    plt.close(fig)

print('Saved results:', source)
print('Exploratory plots and tables:', output)


# %% [markdown]
# ## Bind tables, fits, and metadata
#
# Validate that the saved statistics, common gene universe, structure positions, and reviewed segment labels refer to the same analysis run.

# %%
import anndata as ad

required_files = ('integrated_pathway_classification.csv', 'pathway_evidence_atlas.csv',
    'pathway_coverage.csv', 'gene_trajectories.npz', 'structure_expression_audit.csv',
    'conventional_gene_universe.csv', 'conventional_pt_ora.csv',
    'whole_pt_vs_gam_level_genes.csv', 'whole_pt_vs_gam_level_summary.csv',
    'cluster_signed_only_framework_review.csv', 'total_difference_summary.csv', 'run_manifest.json', 'gene_statistics.csv')
for filename in required_files:
    if not (source / filename).is_file():
        raise FileNotFoundError(f'Run notebook 12 first; missing {filename}')
if not coordinate_path.is_file():
    raise FileNotFoundError('Notebook 13 metadata artifact is required for reviewed segment positions.')
manifest12 = json.loads((source / 'run_manifest.json').read_text())
classification = pd.read_csv(source / 'integrated_pathway_classification.csv')
atlas = pd.read_csv(source / 'pathway_evidence_atlas.csv')
coverage = pd.read_csv(source / 'pathway_coverage.csv')
gene_comparison = pd.read_csv(source / 'whole_pt_vs_gam_level_genes.csv', index_col=0)
gene_comparison.index.name = 'gene'
global_summary = pd.read_csv(source / 'whole_pt_vs_gam_level_summary.csv')
ora = pd.read_csv(source / 'conventional_pt_ora.csv')
controlled_summary = pd.read_csv(source / 'total_difference_summary.csv')
cluster_review = pd.read_csv(source / 'cluster_signed_only_framework_review.csv')
with np.load(source / 'gene_trajectories.npz', allow_pickle=False) as archive:
    fit = {key: archive[key] for key in ('genes', 'structures', 'specimen', 'grid', 'knots',
        'mouse', 'human', 'delta', 'se', 'z', 'T_level', 'T_spatial', 'T_total')}
fit_gene_index = pd.Index(fit['genes'])
original_grid = fit['grid']
lo, hi = original_grid.min(), original_grid.max()
grid = (original_grid - lo) / (hi - lo)
if fit['mouse'].shape != (len(fit_gene_index), len(grid)) or not np.all(np.diff(grid) > 0):
    raise ValueError('Saved trajectory axes do not align.')
common_genes = pd.Index(pd.read_csv(source / 'conventional_gene_universe.csv')
    .query('common_count_support').gene)
if set(common_genes) != set(gene_comparison.index) or not np.isfinite(gene_comparison[['stat', 'Z_level']]).all().all():
    raise ValueError('Paired gene statistics must have exactly the common tested universe.')
all_sets = {row.pathway_id: ast.literal_eval(row.genes_present)
    for row in coverage.loc[coverage.tested_here].itertuples()}
gene_sets = {pathway: [gene for gene in all_sets[pathway] if gene in common_genes]
    for pathway in classification.pathway_id}
if not classification.pathway_id.is_unique or not pd.Series(gene_sets).map(len).between(10, 300).all():
    raise ValueError('Common pathway membership must be unique and contain 10–300 measured genes.')
backed = ad.read_h5ad(coordinate_path, backed='r')
obs = backed.obs.loc[fit['structures'], ['sample', 'comparison_species', 'segment_class', 'shared_pseudospace']].copy()
backed.file.close()
structure_audit = pd.read_csv(source / 'structure_expression_audit.csv', index_col=0).loc[obs.index]
np.testing.assert_array_equal(obs['sample'].astype(str), fit['specimen'])
np.testing.assert_allclose(obs.shared_pseudospace, structure_audit.shared_pseudospace, rtol=0, atol=1e-12)
if set(obs.segment_class.astype(str)) != set(segments):
    raise ValueError('Expected reviewed PT-S1/S2/S3 labels.')
if manifest12['n_genes'] != len(fit_gene_index) or manifest12['fit_input_fingerprint'] is None:
    raise ValueError('Notebook 12 run manifest is incompatible with saved fits.')
# Validate trajectory scores against the exact statistics used in the integrated table's source run.
saved_gene_stats = pd.read_csv(source / 'gene_statistics.csv', index_col=0).loc[fit_gene_index]
for statistic in ('T_level', 'T_spatial', 'T_total'):
    np.testing.assert_allclose(saved_gene_stats[statistic], fit[statistic], rtol=1e-9, atol=1e-9)
obs['position'] = (obs.shared_pseudospace - lo) / (hi - lo)
obs['within_display_support'] = obs.position.between(0, 1)
obs.to_csv(source_data / 'structures.csv', index_label='structure_id')
segment_medians = (obs.groupby(['sample', 'segment_class'], observed=True).position.median()
    .groupby('segment_class', observed=True).mean().reindex(segments))
if not np.all(np.diff(segment_medians) > 0):
    raise ValueError('Reviewed segment medians are not ordered along saved PT position.')
transition_guides = (segment_medians.to_numpy()[:-1] + segment_medians.to_numpy()[1:]) / 2
pd.DataFrame({'transition': ['S1/S2', 'S2/S3'], 'position': transition_guides,
    'definition': 'midpoint between equal-specimen means of segment median positions'}).to_csv(
        source_data / 'segment_guides.csv', index=False)
# Same nearest-grid histogram convention as notebook 12; tails use the end grid cells.
edges = np.r_[-np.inf, (original_grid[:-1] + original_grid[1:]) / 2, np.inf]
position_weights = {}
for label in ('whole_PT', *segments):
    weights = []
    for name in sorted(obs['sample'].unique()):
        group = obs.loc[obs['sample'].eq(name)]
        if label != 'whole_PT':
            group = group.loc[group.segment_class.astype(str).eq(label)]
        if len(group) == 0:
            raise ValueError(f'No positions for {name} {label}.')
        counts = np.histogram(group.shared_pseudospace, bins=edges)[0].astype(float)
        weights.append(counts / counts.sum())
    position_weights[label] = np.mean(weights, axis=0)
pd.DataFrame(position_weights, index=grid).rename_axis('position').reset_index().to_csv(
    source_data / 'position_distributions.csv', index=False)
print(f'{len(common_genes):,} common genes, {len(gene_sets):,} pathways, {len(obs):,} structures.')
print('Transition guides:', np.round(transition_guides, 3), '(labels overlap; descriptive only)')

# %% [markdown]
# ## Rank and review example pathways
#
# Apply the prespecified evidence and sensitivity gates, then check that the five named examples remain eligible and nonredundant.

# %%
atlas_columns = ['pathway_id', 'program', 'representative', 'n_testable', 'n_retained',
    'n_planned', 'q_corr', 'correlation_support', 'peak_position', 'peak_divergence',
    'affected_width', 'affected_grid_fraction', 'direction_sign_changes', 'spatial_driver_genes']
candidates = classification.merge(atlas[atlas_columns], on='pathway_id', how='left', validate='one_to_one')
candidates['robust_fraction_all_planned'] = candidates.n_retained / candidates.n_planned
ora_min = (ora.loc[ora.scope.eq('species_within_segment')].groupby('pathway_id').q_family.min())
candidates['conventional_ora_min_q'] = candidates.pathway_id.map(ora_min)
candidates['n_members_common'] = candidates.pathway_id.map({p: len(gs) for p, gs in gene_sets.items()})
candidates['native_bulk_hit'] = candidates.q_family_DESeq2.le(ALPHA)
candidates['native_segment_hit'] = candidates[[f'cluster_q_family_{s}' for s in segments]].le(ALPHA).any(axis=1)
candidates['unsigned_bulk_hit'] = candidates.common_effect_T_whole_pt_deseq2_abs.gt(0) & candidates.common_q_empirical_T_whole_pt_deseq2_abs.le(ALPHA)
candidates['original_spatial_supported'] = candidates.original_effect_T_spatial.gt(0) & candidates.original_q_empirical_T_spatial.le(ALPHA)
candidates['robust_spatial'] = (candidates.spatial_hit & candidates.original_spatial_supported
    & candidates.q_corr.le(ALPHA) & candidates.robust_fraction_all_planned.ge(MIN_RETENTION))
candidates['class_A_eligible'] = candidates.robust_spatial & candidates.native_segment_hit
candidates['class_B_eligible'] = (candidates.robust_spatial & ~candidates.native_bulk_hit
    & ~candidates.native_segment_hit & ~candidates.any_segment_magnitude_hit
    & ~candidates.unsigned_bulk_hit & candidates.conventional_ora_min_q.gt(ALPHA))
shape_records, candidate_curve_parts, driver_map = [], [], {}
for row in candidates.itertuples():
    idx = fit_gene_index.get_indexer(gene_sets[row.pathway_id])
    member_delta = fit['delta'][idx]
    mean_delta = member_delta.mean(axis=0)
    member_scores = saved_gene_stats.loc[gene_sets[row.pathway_id], 'T_spatial']
    driver_map[row.pathway_id] = member_scores.sort_values(ascending=False, kind='stable').head(3).index.tolist()
    averages = {label: float(mean_delta @ weights) for label, weights in position_weights.items()}
    whole = averages['whole_PT']
    departure = mean_delta - whole
    peak = int(np.argmax(np.abs(departure)))
    active = np.abs(departure) >= .5 * np.max(np.abs(departure))
    left = right = peak
    while left > 0 and active[left - 1]:
        left -= 1
    while right < len(grid) - 1 and active[right + 1]:
        right += 1
    # A coarse value cannot describe within-region fluctuation; use its own observed-position weights.
    within_segment_rms = np.mean([np.sqrt(np.sum(position_weights[s] * (mean_delta - averages[s]) ** 2)) for s in segments])
    absolute_mean = float(np.sum(position_weights['whole_PT'] * np.abs(mean_delta)))
    cancellation = 1 - abs(whole) / absolute_mean if absolute_mean > 1e-12 else 0.
    shape_records.append({'pathway_id': row.pathway_id,
        'mean_delta_min': mean_delta.min(), 'mean_delta_max': mean_delta.max(),
        'mean_delta_range': np.ptp(mean_delta), 'within_segment_rms': within_segment_rms,
        'whole_PT_fitted_mean_delta': whole, **{s + '_fitted_mean_delta': averages[s] for s in segments},
        'whole_PT_cancellation_fraction': cancellation,
        'peak_to_coarse_attenuation': 1 - max(abs(averages[s]) for s in segments) / max(np.max(np.abs(mean_delta)), 1e-12),
        'mean_curve_direction_reversal': mean_delta.min() < -.05 and mean_delta.max() > .05,
        'n_member_direction_reversals': int(((member_delta.min(axis=1) < -.05) & (member_delta.max(axis=1) > .05)).sum()),
        'spatial_departure_peak': grid[peak], 'spatial_departure_interval_left': grid[left],
        'spatial_departure_interval_right': grid[right], 'spatial_departure_width': (right-left)/max(len(grid)-1, 1),
        'peak_interval_crosses_guide': bool(np.any((transition_guides > grid[left]) & (transition_guides < grid[right]))),
        'major_driver_genes': ';'.join(driver_map[row.pathway_id]), 'n_major_driver_genes': len(driver_map[row.pathway_id])})
    candidate_curve_parts.append(pd.DataFrame({'pathway_id': row.pathway_id, 'position': grid,
        'original_position': original_grid, 'mean_delta': mean_delta, 'departure_from_whole_PT': departure}))
candidates = candidates.merge(pd.DataFrame(shape_records), on='pathway_id', validate='one_to_one')
# Exploratory ordering by observed common-universe spatial AUC effect; no new test is performed.
candidates['resolution_usefulness_score'] = candidates.common_effect_T_spatial.abs().rank(pct=True)
candidates['pt_relevance'] = 'Unreviewed: inspect measured member genes; pathway title alone is insufficient.'
candidates['example_class'] = np.select([candidates.class_A_eligible, candidates.class_B_eligible],
    ['Conventional evidence + spatial remodeling', 'Spatial evidence without native signed-GSEA significance'], default='not eligible for examples')
candidates = candidates.sort_values(['resolution_usefulness_score', 'pathway_id'], ascending=[False, True])
candidates['manuscript_rank_within_class'] = candidates.groupby('example_class').cumcount() + 1
candidate_curves = pd.concat(candidate_curve_parts, ignore_index=True)
candidate_curves.to_csv(source_data / 'candidate_mean_difference_curves.csv', index=False)
display(candidates.loc[candidates.class_A_eligible | candidates.class_B_eligible,
    ['pathway_name', 'example_class', 'resolution_usefulness_score', 'q_family_DESeq2',
     'cluster_best_q', 'common_q_empirical_T_spatial', 'robust_fraction_all_planned',
     'common_auc_T_spatial', 'common_effect_T_spatial', 'common_q_empirical_T_spatial',
     'robust_fraction_all_planned', 'major_driver_genes']].head(18))

# %% [markdown]
# ## Freeze five illustrative examples
#
# Two examples have conventional evidence; three are robust spatial calls without native conventional signed-GSEA calls.

# %%
selection_review = [
    ('A', 'KEGG_2019_Mouse::Drug metabolism', 'Drug metabolism',
     'Xenobiotic-handling annotation; Cyp2e1/Ces1d/Adh1 measured in PT.',
     'Native segment finding; the enrichment rank statistic independently supports position-dependent pathway remodeling.'),
    ('A', 'KEGG_2019_Mouse::Citrate cycle (TCA cycle)', 'Citrate cycle',
     'Mitochondrial carbon metabolism; Pck1/Idh1/Pcx measured in PT.',
     'Broad conventional finding with additional evidence for position-dependent remodeling.'),
    ('B', 'Reactome_2022::Transport Of Bile Salts And Organic Acids, Metal Ions And Amine Compounds R-HSA-425366',
     'Organic-acid / bile-salt transport', 'Measured Slc22a6/Slc13a3/Slc6a18 support a PT transport annotation.',
     'Robust spatial enrichment without native bulk or segment signed-GSEA significance.'),
    ('B', 'Reactome_2022::Mitochondrial Fatty Acid Beta-Oxidation R-HSA-77289',
     'Mitochondrial fatty-acid oxidation', 'PT metabolic annotation; Acsm3/Acadm/Acaa2 are measured members.',
     'Robust spatial enrichment without native bulk or segment signed-GSEA significance.'),
    ('B', 'Reactome_2022::Retinoid Metabolism And Transport R-HSA-975634',
     'Retinoid metabolism / transport', 'Measured Rbp4/Lpl/Akr1c18 members; hypothesis-generating retinoid/lipid annotation.',
     'Robust spatial enrichment without native bulk or segment signed-GSEA significance.')]
canonical_ids = {pathway.strip(): pathway for pathway in candidates.pathway_id}
if len(canonical_ids) != len(candidates):
    raise ValueError('Whitespace normalization collides across pathway IDs.')
selection = pd.DataFrame(selection_review, columns=['figure_class', 'requested_id', 'display_name',
    'pt_relevance_review', 'selection_rationale'])
selection['pathway_id'] = selection.requested_id.map(canonical_ids)
if selection.pathway_id.isna().any():
    raise ValueError('A reviewed example is absent from the current common universe.')
selection = selection.merge(candidates, on='pathway_id', validate='one_to_one')
for row in selection.itertuples():
    if not getattr(row, 'class_' + row.figure_class + '_eligible'):
        raise ValueError(f'Review the changed evidence before plotting {row.display_name}.')
if selection.program.isna().any() or selection.program.duplicated().any():
    raise ValueError('Selected examples must be assigned to distinct saved programs.')
redundancy_rows = []
for i, first in enumerate(selection.pathway_id):
    for second in selection.pathway_id.iloc[:i]:
        a, b = set(gene_sets[first]), set(gene_sets[second])
        overlap = len(a & b) / len(a | b)
        redundancy_rows.append({'first': first, 'second': second, 'jaccard': overlap})
        if overlap > .35:
            raise ValueError('Selected terms are redundant; revise the reviewed selection.')
for row in selection.itertuples():
    candidates.loc[candidates.pathway_id.eq(row.pathway_id), 'pt_relevance'] = row.pt_relevance_review
candidates['selected'] = candidates.pathway_id.isin(selection.pathway_id)
candidates.to_csv(output / 'manuscript_candidate_ranking.csv', index=False)
selection.to_csv(output / 'selected_examples.csv', index=False)
pd.DataFrame(redundancy_rows).to_csv(output / 'selected_example_redundancy.csv', index=False)
display(selection[['figure_class', 'display_name', 'q_family_DESeq2', 'cluster_best_q',
    'common_q_empirical_T_level', 'common_q_empirical_T_spatial', 'common_q_empirical_T_total',
    'robust_fraction_all_planned', 'q_corr', 'selection_rationale']])

# %% [markdown]
# ## Selected pathway evidence
#
# Signed GSEA NES summarizes directional gene-set association for the native whole-PT or segment analysis. `T_spatial` pathway enrichment is an unsigned rank-AUC test of gene-level spatial statistics against a covariate-matched null. Their effects and null distributions differ, so a fitted mean curve or nonsignificant signed NES does not explain a matched-null spatial result. The right-hand plot shows the actual gene ranks used by spatial enrichment; its q-value is the saved matched-null result, not a test derived from the ECDF.
#
# The separate member-gene curves use ±1.96 gene-level HC3 standard errors and are conditional on the observed structures and coordinate; the cohort has one human donor.

# %%
conventional_parts, gene_curve_parts = [], []
for _, row in selection.iterrows():
    conventional_parts.append(pd.DataFrame({
        'pathway_id': row['pathway_id'], 'scope': ['whole_PT', *segments],
        'NES': [row['NES_DESeq2'], *[row['cluster_NES_' + s] for s in segments]],
        'q': [row['q_family_DESeq2'], *[row['cluster_q_family_' + s] for s in segments]]}))
    for gene in gene_sets[row['pathway_id']]:
        i = fit_gene_index.get_loc(gene)
        gene_curve_parts.append(pd.DataFrame({'pathway_id': row['pathway_id'], 'gene': gene,
            'position': grid, 'original_position': original_grid,
            'mouse': fit['mouse'][i], 'human': fit['human'][i],
            'delta': fit['delta'][i], 'se_delta': fit['se'][i],
            'working_z': fit['z'][i], 'spatial_enrichment_driver': gene in driver_map[row['pathway_id']]}))
conventional_values = pd.concat(conventional_parts, ignore_index=True)
conventional_values['significant'] = conventional_values.q.le(ALPHA)
conventional_values.to_csv(source_data / 'selected_native_conventional_GSEA.csv', index=False)
gene_curves = pd.concat(gene_curve_parts, ignore_index=True)
gene_curves.to_csv(source_data / 'selected_member_gene_curves.csv', index=False)
gene_membership = pd.DataFrame([
    {'pathway_id': pathway, 'gene': gene, 'common_universe_member': True,
     'T_spatial': saved_gene_stats.loc[gene, 'T_spatial']}
    for pathway in selection.pathway_id for gene in gene_sets[pathway]])
gene_membership.to_csv(source_data / 'selected_full_gene_membership.csv', index=False)
gene_comparison.reset_index().to_csv(source_data / 'whole_PT_vs_GAM_level_genes.csv', index=False)
common_unsigned = classification[['pathway_id', 'common_effect_T_whole_pt_deseq2_abs',
    'common_effect_T_level', 'common_q_empirical_T_whole_pt_deseq2_abs', 'common_q_empirical_T_level']]
common_unsigned.to_csv(source_data / 'unsigned_whole_PT_vs_level.csv', index=False)


# %%
def plot_global_gene_agreement(ax):
    x = gene_comparison.stat.to_numpy()
    y = gene_comparison.Z_level.to_numpy()
    same_direction = np.sign(x) == np.sign(y)
    colors = np.where(same_direction, np.where(x > 0, HUMAN_COLOR, MOUSE_COLOR), '0.65')
    ax.scatter(x, y, s=5, c=colors, alpha=.35, rasterized=True)
    ax.axhline(0, color='0.5', lw=.8)
    ax.axvline(0, color='0.5', lw=.8)
    ax.set(xlabel='Whole-PT DESeq2 Wald statistic', ylabel='GAM-level signed Z',
           title='Broad species signal across common genes')
    rho = gene_comparison.stat.corr(gene_comparison.Z_level, method='spearman')
    agreement = same_direction.mean()
    ax.text(.02, .98, f'Spearman ρ = {rho:.3f}\nDirection agreement = {agreement:.1%}',
            transform=ax.transAxes, va='top')
    return rho, agreement


def plot_native_pathway_evidence(ax, pathway_id):
    rows = conventional_values.loc[conventional_values.pathway_id.eq(pathway_id)]
    colors = np.where(rows.NES > 0, HUMAN_COLOR, MOUSE_COLOR)
    ax.bar(rows.scope, rows.NES, color=colors, alpha=.8)
    ax.margins(y=.2)  # Leave room for the q-value labels above and below bars.
    ax.axhline(0, color='0.4', lw=.8)
    ax.set(ylabel='Signed GSEA NES', title='Native pathway evidence')
    for x, row in enumerate(rows.itertuples()):
        ax.text(x, row.NES + (.08 if row.NES >= 0 else -.12),
                f'q={row.q:.2g}', ha='center', va='bottom' if row.NES >= 0 else 'top', fontsize=8)
    ax.tick_params(axis='x', rotation=20)



# %% [markdown]
# ## Whole-PT broad signal
#
# Compare signed whole-PT DESeq2 statistics with the identifiable signed GAM-level coefficient statistic across the common genes.

# %%
fig, ax = plt.subplots(figsize=(10, 4), layout='constrained')
rho, direction_agreement = plot_global_gene_agreement(ax)
save_plot(fig, 'whole_pt_vs_gam_level')
print(f'Signed-statistic rank correlation: {rho:.3f}; direction agreement: {direction_agreement:.1%}.')

# %% [markdown]
# ## Where the coarse segment labels fall
#
# Show the equal-specimen segment position distributions; transition guides are descriptive and segment labels overlap.

# %%
fig, ax = plt.subplots(figsize=(10, 4), layout='constrained')
for label in segments:
    ax.step(grid, position_weights[label], where='mid', linewidth=1.8, label=label)
ax.set(xlabel='Normalized PT position', ylabel='Specimen-balanced fraction per grid cell',
       title='Reviewed coarse segment labels overlap along pseudospace')
ax.legend(frameon=False)
pd.DataFrame(position_weights, index=grid).rename_axis('position').reset_index().to_csv(
    source_data / 'exploratory_segment_position_weights.csv', index=False)
save_plot(fig, 'segment_position_distributions')

# %% [markdown]
# ## Pathway information profiles
#
# These evidence sets test different alternatives; significance is not expected to nest. The counts use notebook 12 pathway calls and their original FDR conventions.

# %%
information = classification[['pathway_id', 'pathway_name', 'bulk_hit', 'cluster_hit',
    'level_hit', 'spatial_hit', 'total_hit', 'GAM_level_hit', 'smooth_hit']].copy()
information['framework_hit_common'] = information[
    ['level_hit', 'spatial_hit', 'total_hit', 'GAM_level_hit', 'smooth_hit']].any(axis=1)
profiles = (information.groupby(['bulk_hit', 'cluster_hit', 'framework_hit_common'])
    .size().rename('pathways').reset_index().sort_values('pathways', ascending=False))
level_spatial_class = np.select(
    [information.level_hit & information.spatial_hit, information.level_hit, information.spatial_hit],
    ['Level + spatial', 'Level only', 'Spatial only'], default='Neither')
information['level_spatial_class'] = level_spatial_class
class_counts = (information.groupby(['level_spatial_class', 'total_hit']).size()
    .rename('pathways').reset_index())
information.to_csv(source_data / 'exploratory_pathway_calls.csv', index=False)
profiles.to_csv(source_data / 'exploratory_exact_information_profiles.csv', index=False)
class_counts.to_csv(source_data / 'exploratory_level_spatial_classes.csv', index=False)
print('Whole-PT / coarse-segment / continuous-framework evidence profiles:')
display(profiles)
print('Level and spatial classes, with total-test evidence shown separately:')
display(class_counts.pivot(index='level_spatial_class', columns='total_hit', values='pathways').fillna(0).astype(int))


# %%
fig, ax = plt.subplots(figsize=(10, 4), layout='constrained')
counts = information['level_spatial_class'].value_counts().reindex(['Level only', 'Level + spatial', 'Spatial only', 'Neither'], fill_value=0)
counts.plot.bar(ax=ax, color=['#4C78A8', '#72B7B2', '#F58518', '0.75'])
ax.set(xlabel='', ylabel='Pathways', title='Global level and position-dependent pathway calls')
ax.tick_params(axis='x', rotation=15)
save_plot(fig, 'level_spatial_pathway_classes')

# %% [markdown]
# ## Drug metabolism
#
# Native segment evidence and spatial rank evidence are both present; they use different statistics and nulls. The comparison table retains whole-PT, segment, and smooth segment q-values beside level, spatial, and total matched-AUC effects and q-values.

# %%
pathway_id = canonical_ids['KEGG_2019_Mouse::Drug metabolism']
example = selection.loc[selection.pathway_id.eq(pathway_id)].iloc[0]
evidence_rows = [
    {'evidence': 'GAM T_level', 'statistic': 'unsigned AUC − 0.5', 'effect': example.common_effect_T_level,
     'q': example.common_q_empirical_T_level},
    {'evidence': 'GAM T_spatial', 'statistic': 'unsigned AUC − 0.5', 'effect': example.common_effect_T_spatial,
     'q': example.common_q_empirical_T_spatial},
    {'evidence': 'GAM T_total', 'statistic': 'unsigned AUC − 0.5', 'effect': example.common_effect_T_total,
     'q': example.common_q_empirical_T_total},
    {'evidence': 'Whole-PT DESeq2', 'statistic': 'signed NES', 'effect': example.NES_DESeq2,
     'q': example.q_family_DESeq2},
]
for segment in segments:
    evidence_rows.append({'evidence': f'{segment} conventional', 'statistic': 'signed NES',
        'effect': example['cluster_NES_' + segment], 'q': example['cluster_q_family_' + segment]})
    evidence_rows.append({'evidence': f'{segment} smooth segment', 'statistic': 'signed NES',
        'effect': example['smooth_NES_' + segment], 'q': example['smooth_q_family_' + segment]})
print(example.display_name)
display(pd.DataFrame(evidence_rows))
spatial_scores = saved_gene_stats.loc[common_genes, 'T_spatial']
spatial_percentile = spatial_scores.rank(method='average') / len(spatial_scores)
is_member = common_genes.isin(gene_sets[pathway_id])
rank_data = pd.DataFrame({'gene': common_genes, 'T_spatial': spatial_scores.to_numpy(),
    'percentile_rank': spatial_percentile.to_numpy(), 'pathway_member': is_member})
member_indices = common_genes.get_indexer(gene_sets[pathway_id])
observed_auc = rank_auc(spatial_scores.to_numpy(), member_indices)
if not np.isclose(observed_auc, example.common_auc_T_spatial, atol=1e-12, rtol=0):
    raise ValueError('Member set and common T_spatial ranks do not reproduce the saved enrichment AUC.')
rank_data.to_csv(source_data / 'drug_metabolism_spatial_gene_ranks.csv', index=False)

fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout='constrained')
plot_native_pathway_evidence(axes[0], pathway_id)
for member_flag, label, color in [(False, 'Nonmembers', '0.55'), (True, 'Pathway members', '#6A51A3')]:
    values = np.sort(rank_data.loc[rank_data.pathway_member.eq(member_flag), 'percentile_rank'].to_numpy())
    axes[1].step(values, np.arange(1, len(values)+1) / len(values), where='post',
                 label=label, color=color)
axes[1].set(xlim=(0, 1), ylim=(0, 1), xlabel='Common-universe T_spatial percentile rank',
    ylabel='Empirical cumulative fraction', title='Gene ranks used in spatial enrichment')
axes[1].legend(frameon=False)
axes[1].text(.03, .97, f'AUC − 0.5 = {observed_auc - .5:.3f}\nSaved matched-null q = {example.common_q_empirical_T_spatial:.3g}',
    transform=axes[1].transAxes, va='top')
fig.suptitle('Drug metabolism')
save_plot(fig, 'pathway_drug_metabolism_native_and_spatial_ranks')

# %% [markdown]
# ### Member-gene trajectories: Drug metabolism
#
# Curves show three genes selected by their observed common-universe `T_spatial` scores. These fitted gene trajectories are illustrative conditional fits, not pathway enrichment evidence and not proof that aggregation caused a nonsignificant native GSEA result.

# %%
pathway_id = canonical_ids['KEGG_2019_Mouse::Drug metabolism']
fig, ax = plt.subplots(figsize=(10, 4), layout='constrained')
for gene in driver_map[pathway_id]:
    gene_curve = gene_curves.loc[gene_curves.pathway_id.eq(pathway_id) & gene_curves.gene.eq(gene)]
    ax.plot(gene_curve.position, gene_curve.delta, label=gene)
    ax.fill_between(gene_curve.position, gene_curve.delta - 1.96*gene_curve.se_delta,
                    gene_curve.delta + 1.96*gene_curve.se_delta, alpha=.12)
ax.axhline(0, color='0.5', lw=.8)
for guide in transition_guides:
    ax.axvline(guide, color='0.8', ls='--', lw=.8)
ax.set(xlim=(0, 1), xlabel='PT position', ylabel='Gene H − M difference',
       title='Top common-universe T_spatial member genes (illustrative)')
ax.legend(frameon=False)
save_plot(fig, 'drivers_drug_metabolism')

# %% [markdown]
# ## Citrate cycle
#
# Broad conventional finding with additional spatial rank evidence; the two tests ask different questions. The comparison table retains whole-PT, segment, and smooth segment q-values beside level, spatial, and total matched-AUC effects and q-values.

# %%
pathway_id = canonical_ids['KEGG_2019_Mouse::Citrate cycle (TCA cycle)']
example = selection.loc[selection.pathway_id.eq(pathway_id)].iloc[0]
evidence_rows = [
    {'evidence': 'GAM T_level', 'statistic': 'unsigned AUC − 0.5', 'effect': example.common_effect_T_level,
     'q': example.common_q_empirical_T_level},
    {'evidence': 'GAM T_spatial', 'statistic': 'unsigned AUC − 0.5', 'effect': example.common_effect_T_spatial,
     'q': example.common_q_empirical_T_spatial},
    {'evidence': 'GAM T_total', 'statistic': 'unsigned AUC − 0.5', 'effect': example.common_effect_T_total,
     'q': example.common_q_empirical_T_total},
    {'evidence': 'Whole-PT DESeq2', 'statistic': 'signed NES', 'effect': example.NES_DESeq2,
     'q': example.q_family_DESeq2},
]
for segment in segments:
    evidence_rows.append({'evidence': f'{segment} conventional', 'statistic': 'signed NES',
        'effect': example['cluster_NES_' + segment], 'q': example['cluster_q_family_' + segment]})
    evidence_rows.append({'evidence': f'{segment} smooth segment', 'statistic': 'signed NES',
        'effect': example['smooth_NES_' + segment], 'q': example['smooth_q_family_' + segment]})
print(example.display_name)
display(pd.DataFrame(evidence_rows))
spatial_scores = saved_gene_stats.loc[common_genes, 'T_spatial']
spatial_percentile = spatial_scores.rank(method='average') / len(spatial_scores)
is_member = common_genes.isin(gene_sets[pathway_id])
rank_data = pd.DataFrame({'gene': common_genes, 'T_spatial': spatial_scores.to_numpy(),
    'percentile_rank': spatial_percentile.to_numpy(), 'pathway_member': is_member})
member_indices = common_genes.get_indexer(gene_sets[pathway_id])
observed_auc = rank_auc(spatial_scores.to_numpy(), member_indices)
if not np.isclose(observed_auc, example.common_auc_T_spatial, atol=1e-12, rtol=0):
    raise ValueError('Member set and common T_spatial ranks do not reproduce the saved enrichment AUC.')
rank_data.to_csv(source_data / 'citrate_cycle_spatial_gene_ranks.csv', index=False)

fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout='constrained')
plot_native_pathway_evidence(axes[0], pathway_id)
for member_flag, label, color in [(False, 'Nonmembers', '0.55'), (True, 'Pathway members', '#6A51A3')]:
    values = np.sort(rank_data.loc[rank_data.pathway_member.eq(member_flag), 'percentile_rank'].to_numpy())
    axes[1].step(values, np.arange(1, len(values)+1) / len(values), where='post',
                 label=label, color=color)
axes[1].set(xlim=(0, 1), ylim=(0, 1), xlabel='Common-universe T_spatial percentile rank',
    ylabel='Empirical cumulative fraction', title='Gene ranks used in spatial enrichment')
axes[1].legend(frameon=False)
axes[1].text(.03, .97, f'AUC − 0.5 = {observed_auc - .5:.3f}\nSaved matched-null q = {example.common_q_empirical_T_spatial:.3g}',
    transform=axes[1].transAxes, va='top')
fig.suptitle('Citrate cycle')
save_plot(fig, 'pathway_citrate_cycle_native_and_spatial_ranks')

# %% [markdown]
# ### Member-gene trajectories: Citrate cycle
#
# Curves show three genes selected by their observed common-universe `T_spatial` scores. These fitted gene trajectories are illustrative conditional fits, not pathway enrichment evidence and not proof that aggregation caused a nonsignificant native GSEA result.

# %%
pathway_id = canonical_ids['KEGG_2019_Mouse::Citrate cycle (TCA cycle)']
fig, ax = plt.subplots(figsize=(10, 4), layout='constrained')
for gene in driver_map[pathway_id]:
    gene_curve = gene_curves.loc[gene_curves.pathway_id.eq(pathway_id) & gene_curves.gene.eq(gene)]
    ax.plot(gene_curve.position, gene_curve.delta, label=gene)
    ax.fill_between(gene_curve.position, gene_curve.delta - 1.96*gene_curve.se_delta,
                    gene_curve.delta + 1.96*gene_curve.se_delta, alpha=.12)
ax.axhline(0, color='0.5', lw=.8)
for guide in transition_guides:
    ax.axvline(guide, color='0.8', ls='--', lw=.8)
ax.set(xlim=(0, 1), xlabel='PT position', ylabel='Gene H − M difference',
       title='Top common-universe T_spatial member genes (illustrative)')
ax.legend(frameon=False)
save_plot(fig, 'drivers_citrate_cycle')

# %% [markdown]
# ## Organic-acid / bile-salt transport
#
# Robust spatial rank evidence without native whole-PT or segment signed-GSEA significance. The comparison table retains whole-PT, segment, and smooth segment q-values beside level, spatial, and total matched-AUC effects and q-values.

# %%
pathway_id = canonical_ids['Reactome_2022::Transport Of Bile Salts And Organic Acids, Metal Ions And Amine Compounds R-HSA-425366']
example = selection.loc[selection.pathway_id.eq(pathway_id)].iloc[0]
evidence_rows = [
    {'evidence': 'GAM T_level', 'statistic': 'unsigned AUC − 0.5', 'effect': example.common_effect_T_level,
     'q': example.common_q_empirical_T_level},
    {'evidence': 'GAM T_spatial', 'statistic': 'unsigned AUC − 0.5', 'effect': example.common_effect_T_spatial,
     'q': example.common_q_empirical_T_spatial},
    {'evidence': 'GAM T_total', 'statistic': 'unsigned AUC − 0.5', 'effect': example.common_effect_T_total,
     'q': example.common_q_empirical_T_total},
    {'evidence': 'Whole-PT DESeq2', 'statistic': 'signed NES', 'effect': example.NES_DESeq2,
     'q': example.q_family_DESeq2},
]
for segment in segments:
    evidence_rows.append({'evidence': f'{segment} conventional', 'statistic': 'signed NES',
        'effect': example['cluster_NES_' + segment], 'q': example['cluster_q_family_' + segment]})
    evidence_rows.append({'evidence': f'{segment} smooth segment', 'statistic': 'signed NES',
        'effect': example['smooth_NES_' + segment], 'q': example['smooth_q_family_' + segment]})
print(example.display_name)
display(pd.DataFrame(evidence_rows))
spatial_scores = saved_gene_stats.loc[common_genes, 'T_spatial']
spatial_percentile = spatial_scores.rank(method='average') / len(spatial_scores)
is_member = common_genes.isin(gene_sets[pathway_id])
rank_data = pd.DataFrame({'gene': common_genes, 'T_spatial': spatial_scores.to_numpy(),
    'percentile_rank': spatial_percentile.to_numpy(), 'pathway_member': is_member})
member_indices = common_genes.get_indexer(gene_sets[pathway_id])
observed_auc = rank_auc(spatial_scores.to_numpy(), member_indices)
if not np.isclose(observed_auc, example.common_auc_T_spatial, atol=1e-12, rtol=0):
    raise ValueError('Member set and common T_spatial ranks do not reproduce the saved enrichment AUC.')
rank_data.to_csv(source_data / 'organic_acid_transport_spatial_gene_ranks.csv', index=False)

fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout='constrained')
plot_native_pathway_evidence(axes[0], pathway_id)
for member_flag, label, color in [(False, 'Nonmembers', '0.55'), (True, 'Pathway members', '#6A51A3')]:
    values = np.sort(rank_data.loc[rank_data.pathway_member.eq(member_flag), 'percentile_rank'].to_numpy())
    axes[1].step(values, np.arange(1, len(values)+1) / len(values), where='post',
                 label=label, color=color)
axes[1].set(xlim=(0, 1), ylim=(0, 1), xlabel='Common-universe T_spatial percentile rank',
    ylabel='Empirical cumulative fraction', title='Gene ranks used in spatial enrichment')
axes[1].legend(frameon=False)
axes[1].text(.03, .97, f'AUC − 0.5 = {observed_auc - .5:.3f}\nSaved matched-null q = {example.common_q_empirical_T_spatial:.3g}',
    transform=axes[1].transAxes, va='top')
fig.suptitle('Organic-acid / bile-salt transport')
save_plot(fig, 'pathway_organic_acid_transport_native_and_spatial_ranks')

# %% [markdown]
# ### Member-gene trajectories: Organic-acid / bile-salt transport
#
# Curves show three genes selected by their observed common-universe `T_spatial` scores. These fitted gene trajectories are illustrative conditional fits, not pathway enrichment evidence and not proof that aggregation caused a nonsignificant native GSEA result.

# %%
pathway_id = canonical_ids['Reactome_2022::Transport Of Bile Salts And Organic Acids, Metal Ions And Amine Compounds R-HSA-425366']
fig, ax = plt.subplots(figsize=(10, 4), layout='constrained')
for gene in driver_map[pathway_id]:
    gene_curve = gene_curves.loc[gene_curves.pathway_id.eq(pathway_id) & gene_curves.gene.eq(gene)]
    ax.plot(gene_curve.position, gene_curve.delta, label=gene)
    ax.fill_between(gene_curve.position, gene_curve.delta - 1.96*gene_curve.se_delta,
                    gene_curve.delta + 1.96*gene_curve.se_delta, alpha=.12)
ax.axhline(0, color='0.5', lw=.8)
for guide in transition_guides:
    ax.axvline(guide, color='0.8', ls='--', lw=.8)
ax.set(xlim=(0, 1), xlabel='PT position', ylabel='Gene H − M difference',
       title='Top common-universe T_spatial member genes (illustrative)')
ax.legend(frameon=False)
save_plot(fig, 'drivers_organic-acid___bile-salt_transport')

# %% [markdown]
# ## Mitochondrial fatty-acid oxidation
#
# A robust position-dependent pathway result without native whole-PT or segment signed-GSEA significance. Native whole-PT and segment signed GSEA results are shown beside the actual common-universe gene ranks used in spatial enrichment. The comparison table retains whole-PT, segment, and smooth segment q-values beside level, spatial, and total matched-AUC effects and q-values.

# %%
pathway_id = canonical_ids['Reactome_2022::Mitochondrial Fatty Acid Beta-Oxidation R-HSA-77289']
example = selection.loc[selection.pathway_id.eq(pathway_id)].iloc[0]
evidence_rows = [
    {'evidence': 'GAM T_level', 'statistic': 'unsigned AUC − 0.5', 'effect': example.common_effect_T_level,
     'q': example.common_q_empirical_T_level},
    {'evidence': 'GAM T_spatial', 'statistic': 'unsigned AUC − 0.5', 'effect': example.common_effect_T_spatial,
     'q': example.common_q_empirical_T_spatial},
    {'evidence': 'GAM T_total', 'statistic': 'unsigned AUC − 0.5', 'effect': example.common_effect_T_total,
     'q': example.common_q_empirical_T_total},
    {'evidence': 'Whole-PT DESeq2', 'statistic': 'signed NES', 'effect': example.NES_DESeq2,
     'q': example.q_family_DESeq2},
]
for segment in segments:
    evidence_rows.append({'evidence': f'{segment} conventional', 'statistic': 'signed NES',
        'effect': example['cluster_NES_' + segment], 'q': example['cluster_q_family_' + segment]})
    evidence_rows.append({'evidence': f'{segment} smooth segment', 'statistic': 'signed NES',
        'effect': example['smooth_NES_' + segment], 'q': example['smooth_q_family_' + segment]})
print(example.display_name)
display(pd.DataFrame(evidence_rows))
spatial_scores = saved_gene_stats.loc[common_genes, 'T_spatial']
spatial_percentile = spatial_scores.rank(method='average') / len(spatial_scores)
is_member = common_genes.isin(gene_sets[pathway_id])
rank_data = pd.DataFrame({'gene': common_genes, 'T_spatial': spatial_scores.to_numpy(),
    'percentile_rank': spatial_percentile.to_numpy(), 'pathway_member': is_member})
member_indices = common_genes.get_indexer(gene_sets[pathway_id])
observed_auc = rank_auc(spatial_scores.to_numpy(), member_indices)
if not np.isclose(observed_auc, example.common_auc_T_spatial, atol=1e-12, rtol=0):
    raise ValueError('Member set and common T_spatial ranks do not reproduce the saved enrichment AUC.')
rank_data.to_csv(source_data / 'mitochondrial_fatty_acid_oxidation_spatial_gene_ranks.csv', index=False)

fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout='constrained')
plot_native_pathway_evidence(axes[0], pathway_id)
for member_flag, label, color in [(False, 'Nonmembers', '0.55'), (True, 'Pathway members', '#6A51A3')]:
    values = np.sort(rank_data.loc[rank_data.pathway_member.eq(member_flag), 'percentile_rank'].to_numpy())
    axes[1].step(values, np.arange(1, len(values)+1) / len(values), where='post',
                 label=label, color=color)
axes[1].set(xlim=(0, 1), ylim=(0, 1), xlabel='Common-universe T_spatial percentile rank',
    ylabel='Empirical cumulative fraction', title='Gene ranks used in spatial enrichment')
axes[1].legend(frameon=False)
axes[1].text(.03, .97, f'AUC − 0.5 = {observed_auc - .5:.3f}\nSaved matched-null q = {example.common_q_empirical_T_spatial:.3g}',
    transform=axes[1].transAxes, va='top')
fig.suptitle('Mitochondrial fatty-acid oxidation')
save_plot(fig, 'pathway_mitochondrial_fatty_acid_oxidation_native_and_spatial_ranks')

# %% [markdown]
# ### Member-gene trajectories: Mitochondrial fatty-acid oxidation
#
# Curves show three genes selected by their observed common-universe `T_spatial` scores. These fitted gene trajectories are illustrative conditional fits, not pathway enrichment evidence and not proof that aggregation caused a nonsignificant native GSEA result.

# %%
pathway_id = canonical_ids['Reactome_2022::Mitochondrial Fatty Acid Beta-Oxidation R-HSA-77289']
fig, ax = plt.subplots(figsize=(10, 4), layout='constrained')
for gene in driver_map[pathway_id]:
    gene_curve = gene_curves.loc[gene_curves.pathway_id.eq(pathway_id) & gene_curves.gene.eq(gene)]
    ax.plot(gene_curve.position, gene_curve.delta, label=gene)
    ax.fill_between(gene_curve.position, gene_curve.delta - 1.96*gene_curve.se_delta,
                    gene_curve.delta + 1.96*gene_curve.se_delta, alpha=.12)
ax.axhline(0, color='0.5', lw=.8)
for guide in transition_guides:
    ax.axvline(guide, color='0.8', ls='--', lw=.8)
ax.set(xlim=(0, 1), xlabel='PT position', ylabel='Gene H − M difference',
       title='Top common-universe T_spatial member genes (illustrative)')
ax.legend(frameon=False)
save_plot(fig, 'drivers_mitochondrial_fatty-acid_oxidation')

# %% [markdown]
# ## Retinoid metabolism / transport
#
# A robust position-dependent pathway result without native whole-PT or segment signed-GSEA significance. Native whole-PT and segment signed GSEA results are shown beside the actual common-universe gene ranks used in spatial enrichment. The comparison table retains whole-PT, segment, and smooth segment q-values beside level, spatial, and total matched-AUC effects and q-values.

# %%
pathway_id = canonical_ids['Reactome_2022::Retinoid Metabolism And Transport R-HSA-975634']
example = selection.loc[selection.pathway_id.eq(pathway_id)].iloc[0]
evidence_rows = [
    {'evidence': 'GAM T_level', 'statistic': 'unsigned AUC − 0.5', 'effect': example.common_effect_T_level,
     'q': example.common_q_empirical_T_level},
    {'evidence': 'GAM T_spatial', 'statistic': 'unsigned AUC − 0.5', 'effect': example.common_effect_T_spatial,
     'q': example.common_q_empirical_T_spatial},
    {'evidence': 'GAM T_total', 'statistic': 'unsigned AUC − 0.5', 'effect': example.common_effect_T_total,
     'q': example.common_q_empirical_T_total},
    {'evidence': 'Whole-PT DESeq2', 'statistic': 'signed NES', 'effect': example.NES_DESeq2,
     'q': example.q_family_DESeq2},
]
for segment in segments:
    evidence_rows.append({'evidence': f'{segment} conventional', 'statistic': 'signed NES',
        'effect': example['cluster_NES_' + segment], 'q': example['cluster_q_family_' + segment]})
    evidence_rows.append({'evidence': f'{segment} smooth segment', 'statistic': 'signed NES',
        'effect': example['smooth_NES_' + segment], 'q': example['smooth_q_family_' + segment]})
print(example.display_name)
display(pd.DataFrame(evidence_rows))
spatial_scores = saved_gene_stats.loc[common_genes, 'T_spatial']
spatial_percentile = spatial_scores.rank(method='average') / len(spatial_scores)
is_member = common_genes.isin(gene_sets[pathway_id])
rank_data = pd.DataFrame({'gene': common_genes, 'T_spatial': spatial_scores.to_numpy(),
    'percentile_rank': spatial_percentile.to_numpy(), 'pathway_member': is_member})
member_indices = common_genes.get_indexer(gene_sets[pathway_id])
observed_auc = rank_auc(spatial_scores.to_numpy(), member_indices)
if not np.isclose(observed_auc, example.common_auc_T_spatial, atol=1e-12, rtol=0):
    raise ValueError('Member set and common T_spatial ranks do not reproduce the saved enrichment AUC.')
rank_data.to_csv(source_data / 'retinoid_metabolism_transport_spatial_gene_ranks.csv', index=False)

fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout='constrained')
plot_native_pathway_evidence(axes[0], pathway_id)
for member_flag, label, color in [(False, 'Nonmembers', '0.55'), (True, 'Pathway members', '#6A51A3')]:
    values = np.sort(rank_data.loc[rank_data.pathway_member.eq(member_flag), 'percentile_rank'].to_numpy())
    axes[1].step(values, np.arange(1, len(values)+1) / len(values), where='post',
                 label=label, color=color)
axes[1].set(xlim=(0, 1), ylim=(0, 1), xlabel='Common-universe T_spatial percentile rank',
    ylabel='Empirical cumulative fraction', title='Gene ranks used in spatial enrichment')
axes[1].legend(frameon=False)
axes[1].text(.03, .97, f'AUC − 0.5 = {observed_auc - .5:.3f}\nSaved matched-null q = {example.common_q_empirical_T_spatial:.3g}',
    transform=axes[1].transAxes, va='top')
fig.suptitle('Retinoid metabolism / transport')
save_plot(fig, 'pathway_retinoid_metabolism_transport_native_and_spatial_ranks')

# %% [markdown]
# ### Member-gene trajectories: Retinoid metabolism / transport
#
# Curves show three genes selected by their observed common-universe `T_spatial` scores. These fitted gene trajectories are illustrative conditional fits, not pathway enrichment evidence and not proof that aggregation caused a nonsignificant native GSEA result.

# %%
pathway_id = canonical_ids['Reactome_2022::Retinoid Metabolism And Transport R-HSA-975634']
fig, ax = plt.subplots(figsize=(10, 4), layout='constrained')
for gene in driver_map[pathway_id]:
    gene_curve = gene_curves.loc[gene_curves.pathway_id.eq(pathway_id) & gene_curves.gene.eq(gene)]
    ax.plot(gene_curve.position, gene_curve.delta, label=gene)
    ax.fill_between(gene_curve.position, gene_curve.delta - 1.96*gene_curve.se_delta,
                    gene_curve.delta + 1.96*gene_curve.se_delta, alpha=.12)
ax.axhline(0, color='0.5', lw=.8)
for guide in transition_guides:
    ax.axvline(guide, color='0.8', ls='--', lw=.8)
ax.set(xlim=(0, 1), xlabel='PT position', ylabel='Gene H − M difference',
       title='Top common-universe T_spatial member genes (illustrative)')
ax.legend(frameon=False)
save_plot(fig, 'drivers_retinoid_metabolism___transport')

# %% [markdown]
# ## Cluster signed-GSEA-only diagnostic
#
# The diagnostic asks whether pathways called only in the conventional cluster-versus-smooth signed-GSEA comparison have support in any other continuous-framework test or are near its threshold.

# %%
fate_labels = ['Significant elsewhere', 'Near: q < 0.10', 'Near: 0.10 ≤ q < 0.15', 'No q < 0.15']
fate_counts = [int(cluster_review.any_framework_hit.sum()), int(cluster_review.near_threshold_010.sum()),
    int((~cluster_review.any_framework_hit & cluster_review.framework_q_lt_015 & ~cluster_review.framework_q_lt_010).sum()),
    int((~cluster_review.any_framework_hit & ~cluster_review.framework_q_lt_015).sum())]
fate_table = pd.DataFrame({'fate': fate_labels, 'pathways': fate_counts})
fate_table.to_csv(source_data / 'exploratory_cluster_only_fate.csv', index=False)
display(fate_table)
fig, ax = plt.subplots(figsize=(10, 4), layout='constrained')
ax.barh(fate_table.fate, fate_table.pathways, color=['#72B7B2', '#E5C494', '#F58518', '0.75'])
ax.set(xlabel='Pathways', title=f'{len(cluster_review)} cluster signed-GSEA-only pathways')
ax.invert_yaxis()
save_plot(fig, 'cluster_signed_only_fate')

# %% [markdown]
# ## Controlled discrete-versus-continuous benchmark
#
# Section 9c holds much of the statistical machinery constant. Treat it as a robustness/control analysis, separate from the native conventional-versus-pseudostructure comparison.

# %%
control = controlled_summary.set_index('comparison').loc['controlled_discrete_to_continuous']
control_table = pd.DataFrame([{'discrete_hits': int(control.baseline_hits), 'shared': int(control.overlap),
    'continuous_added': int(control.continuous_added)}])
control_table.to_csv(source_data / 'exploratory_controlled_benchmark.csv', index=False)
display(control_table)
fig, ax = plt.subplots(figsize=(10, 3.5), layout='constrained')
ax.barh(['Discrete hits'], [control.baseline_hits], color='0.7', label='Discrete significant')
ax.barh(['Continuous total-difference'], [control.overlap], color='#72B7B2', label='Shared')
ax.barh(['Continuous total-difference'], [control.continuous_added], left=[control.overlap], color='#F58518', label='Additional')
ax.set(xlabel='Pathways', title='Controlled matched-model comparison')
ax.legend(frameon=False)
save_plot(fig, 'controlled_discrete_continuous')

# %% [markdown]
# ## Manuscript-oriented interpretation
#
# Whole-PT and S1/S2/S3 analyses recover many strong species differences. The GAM/pseudostructure framework retains broad differences through level and total-difference tests. Its added spatial test uses gene-level `T_spatial` ranks and an unsigned, covariate-matched pathway null to test position-dependent remodeling. Native signed GSEA and spatial rank-AUC evidence target different alternatives and need not agree pathway by pathway. Example gene curves are descriptive fitted trajectories, not pathway enrichment results and not proof that aggregation caused a particular native GSEA result.

# %%
signed = global_summary.set_index('comparison').loc['signed GSEA']
unsigned = global_summary.set_index('comparison').loc['unsigned matched AUC']
print(f'Whole-PT DESeq2 versus GAM-level signed GSEA: {int(signed.shared)} shared, {int(signed.DESeq2_only)} DESeq2-only, {int(signed.GAM_level_only)} GAM-level-only; gene-statistic Spearman ρ={rho:.3f}.')
print(f'Unsigned matched pathway comparison: {int(unsigned.shared)}/{int(unsigned.DESeq2_hits)} bulk calls shared.')
print('Level/spatial classes:', counts.to_dict())
print(f'Cluster signed-GSEA-only fates: {dict(zip(fate_labels, fate_counts))}')
print(f'Controlled discrete-to-continuous comparison: {int(control.overlap)}/{int(control.baseline_hits)} shared; {int(control.continuous_added)} additional.')
summary = ('Conventional aggregation recovers broad species differences. The trajectory framework represents broad/global shifts through T_level and T_total, while T_spatial tests position-dependent remodeling using an unsigned matched-null rank statistic. Native signed GSEA and spatial rank-AUC evidence address different alternatives and need not nest.\n')
(output / 'exploratory_manuscript_interpretation.txt').write_text(summary)
print(summary)
