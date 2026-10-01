# %% [markdown]
# # 14 · The PT resolution story: manuscript figures from saved results
#
# **Question:** what biological information is retained, localized, and exposed as PT representation progresses from whole-PT aggregation to S1/S2/S3 and reconstructed pseudospace?
#
# This is a figure-generation notebook. It reads notebook 12's saved fits, pathway tests and robustness annotations, plus notebook 13's saved metadata. It does not estimate a coordinate, refit gene models, rerun DESeq2/GSEA, or introduce a new statistical benchmark.
#
# The paper contribution is **fine-grained pseudostructure**: conventional biology remains accessible while position-dependent modeling reveals organization that averages cannot display. `T_level`, `T_spatial`, and `T_total` jointly describe the framework. FDR-thresholded hit sets from different tests need not nest. All findings are exploratory: two healthy mouse specimens and two healthy **cortex** sections from one human donor; the medulla-labelled human section is also cortex.
#
# | Figure | Claim and evidence role | Destination |
# |---|---|---|
# | 1 | Information compression and preservation of broad gene signal | Main text; schematic-led comparison |
# | 2 | Shared conventional pathways acquire finer positional detail | Main text; representative examples |
# | 3 | Robust native-screen-negative pathways have structured fitted effects | Main text; decomposition of aggregation |
# | 4 | Information profiles overlap without requiring nested hit sets | Main text; full common-universe summary |
# | 5 | Recovery, screen discordances, and deliberately matched benchmark | Supplement; validation controls |
#
# Figures use Python/matplotlib, editable PDF/SVG and 600-dpi PNG, with tidy source data and legends. Selection is descriptive and transparent, **not an additional discovery test**.

# %% [markdown]
# ## 1 · Locate saved results and set the figure contract
#
# CLI flags override `PSEUDOSPACE_*` roots. Outputs remain under the configured results root, outside Git. The selected examples are reviewed pathway IDs, not an automatic list of the smallest q-values. Each run rechecks their native-screen status, robustness, and nonredundancy.
#
# The figure contract is 180 mm width, ≥5 pt glyphs, consistent human-orange/mouse-blue encoding, early→late PT orientation, and the saved common-support grid only. Pointwise HC3 bands describe conditional structure uncertainty, not independent-donor or simultaneous confidence intervals. The figure auditors use an optional configurable scripts directory; exports explicitly report when those auditors are unavailable.

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
from pseudospace.stage_cache import cached_frame, digest
from pseudospace import resolution_story

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
FIGURE_LOGIC_VERSION = '14-resolution-story-v1'
FIGURE_WIDTH_IN = 180 / 25.4
ALPHA = .05
MIN_RETENTION = .75  # Retained / ALL planned runs; unavailable runs do not improve this fraction.
segments = ('PT-S1', 'PT-S2', 'PT-S3')
HUMAN_COLOR, MOUSE_COLOR = '#D55E00', '#0072B2'
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['DejaVu Sans'],
    'font.size': 7, 'axes.labelsize': 7, 'axes.titlesize': 8, 'legend.fontsize': 6,
    'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': .7,
    'lines.linewidth': 1.3, 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
qa_scripts = Path(os.environ.get('PSEUDOSPACE_FIGURE_QA_SCRIPTS') or
    Path.home() / '.codex' / 'skills' / 'nature-figure' / 'scripts')
require_matplotlib_panel_alignment = None
qa_available = (qa_scripts / 'audit_panel_alignment.py').is_file()
if qa_available:
    sys.path.insert(0, str(qa_scripts))
    from audit_panel_alignment import require_matplotlib_panel_alignment
print('Read-only figure source:', source)
print('Figure output:', output)

# %% [markdown]
# ## 2 · Bind tables, fits, and metadata to the same run
#
# The common count-supported gene/pathway universe is notebook 12's conventional comparison universe. Pathway-mean and gene curves use its members consistently; notebook 12's local rank-divergence curve retains its original all-eligible-gene background and is labelled accordingly. Raw original and common q-values are preserved. A missing or mismatched artifact stops figure generation rather than silently filling in evidence.
#
# Normalize **only the display axis** so the saved common-support grid spans [0,1]. This does not change the coordinate or any fit. The reviewed labels overlap along position; midpoint guides between equal-specimen segment medians are descriptive transitions, not newly imposed bins or anatomical boundaries. Observed segment distributions are exported beside the guides.

# %%
import anndata as ad

required_files = ('integrated_pathway_classification.csv', 'pathway_evidence_atlas.csv',
    'pathway_coverage.csv', 'gene_trajectories.npz', 'structure_expression_audit.csv',
    'conventional_gene_universe.csv', 'conventional_pt_ora.csv', 'pathway_local_curves.csv',
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
local_curves = pd.read_csv(source / 'pathway_local_curves.csv')
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
# ## 3 · Rank manuscript candidates before choosing examples
#
# A shared example must pass native conventional **species-within-segment signed GSEA** and spatial modeling. A conventionally undetected example must fail whole-PT and all three native segment signed GSEA calls; we additionally exclude positive unsigned whole-PT/segment magnitude calls and native species ORA hits. Neither absence nor nonsignificance proves no effect. This conservative Class B rule uses the **native conventional pipeline**, never Section 9c's artificial matched benchmark.
#
# Both example classes require original and common positive spatial enrichment, residual-correlation support q ≤ 0.05, and retention in ≥75% of **all planned** sensitivity runs. The complete candidate table preserves failed gates, pathway size/program, driver identities, peak/width, sign reversals, boundary-spanning intervals, and native evidence. The manuscript-usefulness rank combines fitted positional variation, departure from segment averages, and reversals; it is a presentation aid, not a new statistical score. PT relevance is reviewed from actual member genes for selected examples and left explicitly unreviewed elsewhere.
#
# A **compression diagnostic** averages each fixed fitted difference over the observed specimen-balanced position distributions. Peak attenuation and cancellation quantify what these summaries discard while holding the curve fixed. They do not establish the cause of native GSEA's different q-values, which also depend on estimators and pathway nulls.

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
    member_variation = member_delta.std(axis=1)
    driver_indices = idx[np.argsort(-member_variation, kind='stable')[:3]]
    driver_map[row.pathway_id] = fit_gene_index[driver_indices].tolist()
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
        'major_driver_genes': ';'.join(driver_map[row.pathway_id]), 'n_major_driver_genes': len(driver_indices)})
    candidate_curve_parts.append(pd.DataFrame({'pathway_id': row.pathway_id, 'position': grid,
        'original_position': original_grid, 'mean_delta': mean_delta, 'departure_from_whole_PT': departure}))
candidates = candidates.merge(pd.DataFrame(shape_records), on='pathway_id', validate='one_to_one')
candidates['resolution_usefulness_score'] = (.5 * candidates.mean_delta_range.rank(pct=True)
    + .3 * candidates.within_segment_rms.rank(pct=True)
    + .2 * candidates.n_member_direction_reversals.rank(pct=True))
candidates['pt_relevance'] = 'Unreviewed: inspect measured member genes; pathway title alone is insufficient.'
candidates['example_class'] = np.select([candidates.class_A_eligible, candidates.class_B_eligible],
    ['A: shared; finer resolution', 'B: native-screen-negative; spatial'], default='not eligible for examples')
candidates = candidates.sort_values(['resolution_usefulness_score', 'pathway_id'], ascending=[False, True])
candidates['manuscript_rank_within_class'] = candidates.groupby('example_class').cumcount() + 1
candidate_curves = pd.concat(candidate_curve_parts, ignore_index=True)
candidate_curves.to_csv(source_data / 'candidate_mean_difference_curves.csv', index=False)
display(candidates.loc[candidates.class_A_eligible | candidates.class_B_eligible,
    ['pathway_name', 'example_class', 'resolution_usefulness_score', 'q_family_DESeq2',
     'cluster_best_q', 'common_q_empirical_T_spatial', 'robust_fraction_all_planned',
     'mean_curve_direction_reversal', 'within_segment_rms', 'major_driver_genes']].head(18))

# %% [markdown]
# ## 4 · Freeze a small, nonredundant example set
#
# **Class A:** drug metabolism and the citrate cycle. The former shows changing and reversing average effects and the latter resolves a broad mitochondrial metabolic difference into a position-dependent profile. These are distinct gene sets/programs.
#
# **Class B:** organic-acid/bile-salt transport, mitochondrial fatty-acid beta-oxidation, and retinoid metabolism/transport. These give a reversal/transport example, a gradually emerging metabolic difference, and heterogeneous lipid/retinoid member trajectories. Native conventional evidence is inspected for every chosen term; a title such as “bile salts” is not treated as proof of bile-acid flux or kidney mechanism.
#
# PT metabolic and organic-anion handling relevance is supported by primary experimental studies of [tubular fatty-acid oxidation](https://www.nature.com/articles/nm.3762) and [OAT1-linked transport/metabolism](https://pmc.ncbi.nlm.nih.gov/articles/PMC3173137/). These references support biological relevance, not validation of this cohort's species differences. The retinoid example is a hypothesis-generating annotation based on measured PT members and requires biological review.
#
# Selected IDs and rationale are editable below. Every selected term must still pass its class gate. A duplicate program or Jaccard overlap >0.35 stops the figure rather than silently accepting redundant examples.

# %%
selection_review = [
    ('A', 'KEGG_2019_Mouse::Drug metabolism', 'Drug metabolism',
     'Xenobiotic-handling annotation; Cyp2e1/Ces1d/Adh1 measured in PT.',
     'Native segment finding with a changing and sign-reversing member-average profile.'),
    ('A', 'KEGG_2019_Mouse::Citrate cycle (TCA cycle)', 'Citrate cycle',
     'Mitochondrial carbon metabolism; Pck1/Idh1/Pcx measured in PT.',
     'Broad conventional finding with a mid-PT trough and within-segment variation.'),
    ('B', 'Reactome_2022::Transport Of Bile Salts And Organic Acids, Metal Ions And Amine Compounds R-HSA-425366',
     'Organic-acid / bile-salt transport', 'Measured Slc22a6/Slc13a3/Slc6a18 support a PT transport annotation.',
     'Opposing fitted effects along position; averaging dilutes the signed mean.'),
    ('B', 'Reactome_2022::Mitochondrial Fatty Acid Beta-Oxidation R-HSA-77289',
     'Mitochondrial fatty-acid oxidation', 'PT metabolic annotation; Acsm3/Acadm/Acaa2 are measured members.',
     'Gradually emerging difference with heterogeneous member directions.'),
    ('B', 'Reactome_2022::Retinoid Metabolism And Transport R-HSA-975634',
     'Retinoid metabolism / transport', 'Measured Rbp4/Lpl/Akr1c18 members; hypothesis-generating retinoid/lipid annotation.',
     'A changing signed mean and mixed gene trajectories are hidden by a small number of summaries.')]
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
# ## 5 · Summarize saved pathway fits and export the plotted measurements
#
# Each pathway trajectory is the **equal-gene mean of saved log-normalized fits** over the common pathway membership. This abundance summary is not GSEA NES or local rank-divergence. The separate `D(s)` trace is notebook 12's competitive, unsigned local divergence in its original background.
#
# For pathway-mean uncertainty, average the saved residuals **before** forming the HC3 sandwich. This retains cross-gene residual covariance; averaging independent gene SEs would give an invalidly narrow band. The saved weighted full-model design is reconstructed algebraically, with equal species/specimen weights and within-species sum-to-zero specimen offsets. No expression model is refitted. Large residual arrays are read only on a cache miss and then released.
#
# The driver panels show the three member genes with largest fitted spatial SD, selected before drawing. All selected-member curves are exported, including undisplayed members. Driver roles explain the curve and are not extra significance calls.

# %%
def summarize_selected_pathways():
    with np.load(source / 'gene_trajectories.npz', allow_pickle=False) as archive:
        saved = {**fit, 'residuals': archive['residuals']}
        return resolution_story.aggregate_saved_trajectories(saved,
            obs.shared_pseudospace.to_numpy(float),
            obs.comparison_species.astype(str).eq('human').to_numpy(),
            obs['sample'].astype(str).to_numpy(),
            {p: gene_sets[p] for p in selection.pathway_id})

pathway_curves = cached_frame('selected_pathway_mean_hc3', summarize_selected_pathways,
    root=output / 'stage_cache', params={'figure_logic': FIGURE_LOGIC_VERSION},
    inputs={'fit': source / 'gene_trajectories.npz', 'obs': obs,
            'members': {p: gene_sets[p] for p in selection.pathway_id}},
    code=digest([Path(resolution_story.__file__), FIGURE_LOGIC_VERSION]))
pathway_curves.to_csv(source_data / 'selected_pathway_mean_curves.csv', index=False)
conventional_parts, gene_curve_parts, projection_parts = [], [], []
for row in selection.itertuples():
    conventional_parts.append(pd.DataFrame({
            'pathway_id': row.pathway_id, 'scope': ['whole_PT', *segments],
            'NES': [row.NES_DESeq2, *[selection.loc[selection.pathway_id.eq(row.pathway_id), 'cluster_NES_' + s].iloc[0] for s in segments]],
            'q': [row.q_family_DESeq2, *[selection.loc[selection.pathway_id.eq(row.pathway_id), 'cluster_q_family_' + s].iloc[0] for s in segments]]}))
    curve = pathway_curves.loc[pathway_curves.pathway_id.eq(row.pathway_id)]
    for label, weights in position_weights.items():
        projection_parts.append({'pathway_id': row.pathway_id, 'scope': label,
            'fitted_mean_delta': float(curve.delta.to_numpy() @ weights),
            'position_median': .5 if label == 'whole_PT' else segment_medians[label],
            'meaning': 'compression of same fitted curve, not DESeq2 estimate'})
    for gene in gene_sets[row.pathway_id]:
        i = fit_gene_index.get_loc(gene)
        gene_curve_parts.append(pd.DataFrame({'pathway_id': row.pathway_id, 'gene': gene,
            'position': grid, 'original_position': original_grid,
            'mouse': fit['mouse'][i], 'human': fit['human'][i],
            'delta': fit['delta'][i], 'se_delta': fit['se'][i],
            'working_z': fit['z'][i], 'plotted_driver': gene in driver_map[row.pathway_id]}))
conventional_values = pd.concat(conventional_parts, ignore_index=True)
conventional_values['significant'] = conventional_values.q.le(ALPHA)
conventional_values.to_csv(source_data / 'selected_native_conventional_GSEA.csv', index=False)
gene_curves = pd.concat(gene_curve_parts, ignore_index=True)
gene_curves.to_csv(source_data / 'selected_member_gene_curves.csv', index=False)
projections = pd.DataFrame(projection_parts)
projections.to_csv(source_data / 'selected_same_curve_compression.csv', index=False)
selected_local = local_curves.loc[local_curves.pathway_id.isin(selection.pathway_id)].copy()
selected_local['original_position'] = selected_local.position
selected_local['position'] = (selected_local.position - lo) / (hi - lo)
selected_local.to_csv(source_data / 'selected_local_rank_divergence.csv', index=False)
gene_comparison.reset_index().to_csv(source_data / 'whole_PT_vs_GAM_level_genes.csv', index=False)
common_unsigned = classification[['pathway_id', 'common_effect_T_whole_pt_deseq2_abs',
    'common_effect_T_level', 'common_q_empirical_T_whole_pt_deseq2_abs', 'common_q_empirical_T_level']]
common_unsigned.to_csv(source_data / 'unsigned_whole_PT_vs_level.csv', index=False)

# %% [markdown]
# ## 6 · Shared drawing/export helpers and final rendered QA
#
# Every figure has labelled panels, a legend file specifying the biological unit, score and uncertainty, and matching source-data files. Frame counts and set calls are descriptive summaries, so their bars have no artificial sampling error bars. Numerical scatter panels retain all paired observations and use rasterized marks in otherwise editable vector figures.
#
# The alignment gate checks final axes rectangles after layout. When auditor scripts are present, PDF glyph and collision audits run after every export. A failure stops the notebook; warnings remain visible for final-size review. Conceptual diagrams have no measured uncertainty. Reported screen unions are descriptive and are not a joint error-controlled claim across different test families.

# %%
import subprocess
from textwrap import fill

figure_records = []

def panel_label(ax, label):
    ax.annotate(label, xy=(0, 1), xycoords='axes fraction', xytext=(-12, 8),
        textcoords='offset points', fontsize=8, fontweight='bold', ha='left', va='bottom')

def export_figure(fig, stem, caption, *, exclude_axes=None, exemptions=None):
    if qa_available:
        require_matplotlib_panel_alignment(fig, json_out=str(output / (stem + '.alignment.json')),
            exclude_axes=exclude_axes or (), exemptions=exemptions or (), strict=True)
    else:
        (output / (stem + '.alignment.json')).write_text(json.dumps(
            {'status': 'NOT AUDITED', 'reason': 'figure auditor scripts unavailable'}))
    fig.savefig(output / (stem + '.pdf'))
    fig.savefig(output / (stem + '.svg'))
    fig.savefig(output / (stem + '.png'), dpi=600)
    (output / (stem + '_legend.md')).write_text(caption + '\n')
    audits = []
    if qa_available:
        for script, extra in [('audit_pdf_text.py', []), ('audit_figure_collisions.py',
            ['--json-out', str(output / (stem + '.collision.json'))])]:
            result = subprocess.run([sys.executable, str(qa_scripts / script), str(output / (stem + '.pdf')), *extra],
                capture_output=True, text=True)
            (output / (stem + '.' + script.removesuffix('.py') + '.txt')).write_text(result.stdout + result.stderr)
            if result.returncode:
                raise RuntimeError(f'Figure QA requires repair: {stem}, {script}. Read its saved report.')
            audits.append(script)
    figure_records.append({'figure': stem, 'formats': 'PDF; SVG; PNG 600 dpi',
        'audits': ';'.join(audits) if audits else 'not run', 'caption': caption})
    plt.show()
    plt.close(fig)

def broad_gene_scatter(ax, *, annotation=True):
    x, y = gene_comparison.stat.to_numpy(), gene_comparison.Z_level.to_numpy()
    agree = np.sign(x) == np.sign(y)
    colors = np.where(agree, np.where(x > 0, HUMAN_COLOR, MOUSE_COLOR), '0.6')
    ax.scatter(x, y, s=1.8, c=colors, alpha=.25, rasterized=True, linewidths=0)
    ax.axhline(0, color='0.6', lw=.5)
    ax.axvline(0, color='0.6', lw=.5)
    ax.set(xlabel='Whole-PT DESeq2 statistic', ylabel='Signed GAM level Z')
    rho = gene_comparison.stat.corr(gene_comparison.Z_level, method='spearman')
    direction = agree.mean()
    if annotation:
        ax.set_title(f'Spearman ρ = {rho:.3f}\nDirection agreement = {direction:.1%}', fontsize=7)
    return rho, direction

def conventional_panel(ax, pathway):
    values = conventional_values.loc[conventional_values.pathway_id.eq(pathway)]
    for x, row in enumerate(values.itertuples()):
        color = HUMAN_COLOR if row.NES > 0 else MOUSE_COLOR
        ax.plot([x, x], [0, row.NES], color=color, lw=1)
        ax.scatter([x], [row.NES], s=22, edgecolors=color,
            facecolors=color if row.significant else 'white', linewidths=.8, zorder=3)
    ax.axhline(0, color='0.6', lw=.5)
    ax.set(xticks=range(4), xticklabels=['PT', 'S1', 'S2', 'S3'], ylim=(-3, 3), ylabel='Signed GSEA NES')

def positional_guides(ax):
    for guide in transition_guides:
        ax.axvline(guide, color='0.75', lw=.6, ls='--', zorder=0)
    ax.set(xlim=(0, 1), xticks=[0, .5, 1])

def pathway_difference_panel(ax, pathway, *, projection=False):
    curve = pathway_curves.loc[pathway_curves.pathway_id.eq(pathway)]
    ax.fill_between(curve.position, curve.delta - 1.96 * curve.se_delta,
        curve.delta + 1.96 * curve.se_delta, color='0.8', linewidth=0, alpha=.8)
    ax.plot(curve.position, curve.delta, color='0.15')
    ax.axhline(0, color='0.55', lw=.6)
    local = selected_local.loc[selected_local.pathway_id.eq(pathway)].sort_values('position')
    peak = local.loc[local.divergence.idxmax(), 'position']
    ax.axvline(peak, color='#8C6BB1', lw=.8, ls=':', label='Peak local divergence')
    if projection:
        rows = projections.loc[projections.pathway_id.eq(pathway) & projections.scope.ne('whole_PT')]
        ax.scatter(rows.position_median, rows.fitted_mean_delta, marker='s', facecolors='white',
            edgecolors='0.35', s=20, zorder=4)
    positional_guides(ax)
    ax.set(ylabel='Mean H−M log expression')

def species_curve_panel(ax, pathway):
    curve = pathway_curves.loc[pathway_curves.pathway_id.eq(pathway)]
    for species, color in [('human', HUMAN_COLOR), ('mouse', MOUSE_COLOR)]:
        ax.plot(curve.position, curve[species], color=color, label=species.capitalize())
        ax.fill_between(curve.position, curve[species] - 1.96 * curve['se_' + species],
            curve[species] + 1.96 * curve['se_' + species], color=color, alpha=.15, linewidth=0)
    positional_guides(ax)
    ax.set(ylabel='Mean log expression')
    ax.legend(loc='lower left', bbox_to_anchor=(0, 1.16), ncol=2, frameon=False,
              fontsize=5, handlelength=.8, columnspacing=.6, handletextpad=.4, borderaxespad=0)

def gene_driver_panel(ax, pathway):
    driver_colors = ('#1B9E77', '#7570B3', '#A6761D')
    for gene, color in zip(driver_map[pathway], driver_colors):
        curve = gene_curves.loc[gene_curves.pathway_id.eq(pathway) & gene_curves.gene.eq(gene)]
        ax.plot(curve.position, curve.delta, color=color, label=gene)
        ax.fill_between(curve.position, curve.delta - 1.96 * curve.se_delta,
            curve.delta + 1.96 * curve.se_delta, color=color, alpha=.10, linewidth=0)
    ax.axhline(0, color='0.6', lw=.5)
    positional_guides(ax)
    ax.set(ylabel='Gene H−M difference')
    ax.legend(loc='lower right', bbox_to_anchor=(1, 1.26), ncol=3, frameon=False,
              fontsize=5, handlelength=.8, columnspacing=.6, handletextpad=.4, borderaxespad=0)

def position_density_panel(ax):
    for i, label in enumerate(segments):
        weights = position_weights[label]
        ax.fill_between(grid, i, i + .7 * weights / weights.max(), color='0.72', linewidth=0)
    positional_guides(ax)
    ax.set(yticks=[.3, 1.3, 2.3], yticklabels=['S1', 'S2', 'S3'], ylim=(-.1, 3),
           xlabel='Normalized PT position')
    ax.tick_params(axis='y', length=0, labelsize=5.5)
    ax.spines['left'].set_visible(False)


# %% [markdown]
# ## Figure 1 · Reconstructing position retains broad conventional signal
#
# The left panel shows within-specimen compression into one or three profiles versus retaining each structure's continuous coordinate. The right panel compares all common genes' signed whole-PT and level statistics. These different noise-scaled quantities are not expected to have equal numerical values; no identity line or cross-model superiority test is implied. The recovery result is computed from the saved data, not hardcoded.

# %%
fig, axes = plt.subplots(1, 2, figsize=(FIGURE_WIDTH_IN, 82/25.4), gridspec_kw={'width_ratios': [1.7, 1]})
fig.subplots_adjust(left=.04, right=.98, bottom=.18, top=.86, wspace=.42)
ax = axes[0]
ax.set(xlim=(0, 1), ylim=(0, 1))
ax.axis('off')
for y, label, question, analysis, count in (
    (.82, 'Whole PT', 'Overall species difference?', 'DESeq2 / pathway GSEA', 1),
    (.49, 'S1 / S2 / S3', 'Within a coarse PT region?', 'Segment DESeq2 / GSEA', 3),
    (.16, 'Pseudospace', 'How does difference vary with position?', 'GAM: level · spatial · total', 9)):
    xx = np.tile(np.linspace(.02, .15, 4), 3)
    yy = np.repeat([y-.045, y, y+.045], 4)
    ax.scatter(xx, yy, color='0.7', s=9, linewidths=0)
    ax.annotate('', xy=(.29, y), xytext=(.19, y), arrowprops={'arrowstyle': '->', 'lw': .8})
    target_x = np.linspace(.33, .47, count) if count > 1 else [.40]
    ax.scatter(target_x, np.full(count, y), color='0.35', s=24 if count < 4 else 8, linewidths=0)
    ax.text(.54, y+.07, label, fontsize=8, fontweight='bold', va='center')
    ax.text(.54, y-.01, analysis, fontsize=6.5, va='center')
    ax.text(.02, y-.12, question, fontsize=6.5)
ax.text(.02, .98, 'Structures within each specimen', fontsize=6.5)
ax.text(.32, -.06, 'Continuous coordinate retained →', fontsize=6)
for i, axis in enumerate(axes):
    panel_label(axis, chr(97+i))
rho, direction_agreement = broad_gene_scatter(axes[1])
unsigned_row = global_summary.set_index('comparison').loc['unsigned matched AUC']
fig.text(.56, .04, f'Unsigned pathway recovery: {int(unsigned_row.shared)}/{int(unsigned_row.DESeq2_hits)} bulk calls', fontsize=6.5)
export_figure(fig, 'figure_1_resolution_progression',
    'Figure 1 | Reconstructing PT position retains broad species signal. '
    'a, Conceptual within-specimen aggregation: all PT structures → one profile; reviewed labels → three profiles; '
    'pseudospace retains a continuous coordinate. Dot counts are illustrative. '
    'b, All common count-supported genes, whole-PT DESeq2 Wald statistic versus Mlevel human-minus-mouse coefficient / HC3 SE. '
    f'Spearman rho={rho:.3f}; direction agreement={direction_agreement:.1%}. Orange/blue denote human-high/mouse-high concordant directions; gray denotes disagreement. '
    'Statistics have different noise models and structure-level Z is conditional, not donor-calibrated. '
    f'Unsigned matched pathway comparison recovers {int(unsigned_row.shared)}/{int(unsigned_row.DESeq2_hits)} conventional magnitude calls. '
    'Two mouse specimens and two human cortex sections from one donor; exploratory.')


# %% [markdown]
# ## Figure 2 · Shared biological findings, finer positional resolution
#
# Each row is a shared pathway. Left: native whole-PT/S1/S2/S3 signed GSEA, with filled markers at the corresponding full-family q ≤ 0.05 and open markers otherwise. Center: the fixed member-average human-minus-mouse trajectory, with pointwise HC3 uncertainty; the purple dotted guide marks notebook 12's strongest local rank-divergence position. Right: major driver-gene differences, showing variation beyond a segment label. The small density strips show the actual overlap of reviewed segment position distributions.
#
# These are different summaries of the same biological example, not interchangeable estimators. NES has no valid confidence band in the saved results. Gene and pathway bands are conditional on observed structures and the saved coordinate.

# %%
def make_example_figure(example_class, stem):
    examples = selection.loc[selection.figure_class.eq(example_class)].reset_index(drop=True)
    n = len(examples)
    ncols = 3 if example_class == 'A' else 4
    fig, axes = plt.subplots(n+1, ncols,
        figsize=(FIGURE_WIDTH_IN, (140 if example_class == 'A' else 190)/25.4),
        gridspec_kw={'height_ratios': [1]*n + [.28]})
    fig.subplots_adjust(left=.09, right=.98, bottom=.11, top=.89, wspace=.78, hspace=.82)
    headings = ['Native aggregation', 'Continuous H−M', 'Driver genes'] if example_class == 'A' else [
        'Native aggregation', 'Human / mouse', 'Continuous H−M', 'Driver genes']
    for j, heading in enumerate(headings):
        fig.text((axes[0, j].get_position().x0 + axes[0, j].get_position().x1)/2,
                 .97, heading, ha='center', fontsize=8, fontweight='bold')
    letter = 0
    for i, row in enumerate(examples.itertuples()):
        conventional_panel(axes[i, 0], row.pathway_id)
        axes[i, 0].set_title(fill(row.display_name, width=21), fontsize=7, pad=12)
        if example_class == 'A':
            pathway_difference_panel(axes[i, 1], row.pathway_id)
            gene_driver_panel(axes[i, 2], row.pathway_id)
        else:
            species_curve_panel(axes[i, 1], row.pathway_id)
            pathway_difference_panel(axes[i, 2], row.pathway_id, projection=True)
            gene_driver_panel(axes[i, 3], row.pathway_id)
        for ax in axes[i]:
            panel_label(ax, chr(97+letter))
            letter += 1
    axes[-1, 0].axis('off')
    axes[-1, 0].text(0, .88, 'Filled: native q ≤ 0.05\nOpen: q > 0.05\nPT = whole PT', fontsize=6, va='top')
    for ax in axes[-1, 1:]:
        position_density_panel(ax)
    caption = ('Figure 2 | Shared pathways resolved continuously. ' if example_class == 'A' else
               'Figure 3 | Robust spatial pathway calls absent from the native conventional screens. ')
    caption += ('Rows: ' + '; '.join(examples.display_name) + '. '
        'Native signed NES: whole-PT BH across the full pathway family; segment BH pooled across all pathway×segment tests. '
        'Filled markers pass q≤0.05; open markers do not. Positive NES means human-high. '
        'Pathway curves average all common-membership gene fits, not GSEA NES or a pathway activation assay. '
        'Bands are mean ±1.96 pointwise structure-conditional HC3 SE, including cross-gene residual covariance for pathway means; '
        'they ignore independent-donor replication, spatial correlation and coordinate uncertainty. '
        'Driver panels show the three member genes with largest fitted spatial SD; gene colors identify genes, not species. '
        'Purple dotted lines mark maximum original-background D(s), not peak mean expression or a confidence interval. '
        'Gray dashed guides are midpoints of equal-specimen segment medians; reviewed labels overlap. '
        'Density strips are equal-specimen nearest-grid position distributions; tiny tails use end grid cells. '
        'Early→late display is normalized to saved common support; no extrapolation. ')
    if example_class == 'B':
        caption += ('Orange/blue species curves are human/mouse. Open squares in difference panels are the same fitted curve '
            'averaged over each observed segment distribution, plotted at the equal-specimen segment median. '
            'They are compression measurements, not DESeq2 estimates. '
            'All selected terms fail native whole-PT and segment signed GSEA, native segment ORA and the unsigned magnitude companions. '
            'Structured attenuation/cancellation is visible, but this does not establish that aggregation uniquely caused the native q-value discrepancy. ')
    caption += 'Two mouse specimens and two cortex sections from one human donor; exploratory.'
    export_figure(fig, stem, caption, exclude_axes=[axes[-1, 0]])

make_example_figure('A', 'figure_2_shared_pathways_finer_resolution')

# %% [markdown]
# ## Figure 3 · Native-screen-negative pathways with spatially structured signals
#
# Each example passes the robust spatial gates and fails the actual conventional pathway screens. Human/mouse curves, their difference, and major member trajectories expose signed cancellation, gradual emergence, and heterogeneous directions. The open squares summarize the **same fitted curve** over each segment's observed position distribution; the exported compression table quantifies attenuation without switching estimators.
#
# This figure supports “aggregation cannot describe the resolved pattern.” It does **not** prove that model differences played no role in the native screen discrepancy, or that every negative conventional result was caused by averaging. Retain the native q-values, correlation and sensitivity annotations alongside the examples.

# %%
make_example_figure('B', 'figure_3_spatial_native_screen_negative')

# %% [markdown]
# ## Figure 4 · Information profiles at progressively finer representation
#
# An UpSet-style view shows exact combinations of native whole-PT calls, native coarse calls, and **any common-universe framework call**. Framework support includes positive matched level/spatial/total evidence, signed level GSEA or smooth segment GSEA. These screen unions are descriptive and have different test opportunities; no joint FDR or perfect nesting is claimed.
#
# The separate level/spatial breakdown makes the added positional dimension visible while retaining omnibus total support. The large neither-level-nor-spatial background remains in a printed count, including total-only terms. All eight detection profiles and all eight level/spatial/total intersections are exported; no pathway is omitted from the underlying universe.

# %%
information = classification[['pathway_id', 'bulk_hit', 'cluster_hit', 'level_hit', 'spatial_hit',
    'total_hit', 'GAM_level_hit', 'smooth_hit']].copy()
information['framework_hit_common'] = information[['level_hit', 'spatial_hit', 'total_hit', 'GAM_level_hit', 'smooth_hit']].any(axis=1)
profile_names = ['bulk_hit', 'cluster_hit', 'framework_hit_common']
profiles = information.groupby(profile_names).size().rename('pathways').reindex(
    pd.MultiIndex.from_product([[False, True]]*3, names=profile_names), fill_value=0).reset_index()
profiles['n_representations'] = profiles[profile_names].sum(axis=1)
profiles = profiles.sort_values(['n_representations', 'pathways'], ascending=[False, False]).reset_index(drop=True)
profiles.to_csv(source_data / 'figure_4_exact_information_profiles.csv', index=False)
information.to_csv(source_data / 'figure_4_pathway_calls.csv', index=False)
intersections = information.groupby(['level_hit', 'spatial_hit', 'total_hit']).size().rename('pathways').reindex(
    pd.MultiIndex.from_product([[False, True]]*3, names=['level_hit', 'spatial_hit', 'total_hit']), fill_value=0)
intersections.to_csv(source_data / 'figure_4_framework_exact_intersections.csv')
level_class = np.select([information.level_hit & information.spatial_hit,
    information.level_hit, information.spatial_hit], ['Level + spatial', 'Level only', 'Spatial only'], default='Neither')
class_counts = pd.crosstab(level_class, information.total_hit).reindex(
    index=['Level only', 'Level + spatial', 'Spatial only', 'Neither'], columns=[False, True], fill_value=0)
class_counts.to_csv(source_data / 'figure_4_level_spatial_classes.csv')
fig = plt.figure(figsize=(FIGURE_WIDTH_IN, 105/25.4))
gs = fig.add_gridspec(2, 2, width_ratios=[1.65, 1], height_ratios=[2.3, 1],
    left=.14, right=.98, bottom=.17, top=.85, hspace=.3, wspace=.58)
counts_ax = fig.add_subplot(gs[0, 0])
matrix_ax = fig.add_subplot(gs[1, 0], sharex=counts_ax)
framework_ax = fig.add_subplot(gs[:, 1])
x = np.arange(len(profiles))
counts_ax.bar(x, profiles.pathways, color=np.where(profiles.framework_hit_common, '#8DA0CB', '0.75'))
counts_ax.set_yscale('symlog', linthresh=10)
counts_ax.tick_params(axis='y', labelsize=8)
counts_ax.set(ylabel='Pathway terms', title='Exact detection profiles')
counts_ax.tick_params(axis='x', bottom=False, labelbottom=False)
counts_ax.set_ylim(0, max(profiles.pathways.max()*2, 20))
for i, count in enumerate(profiles.pathways):
    counts_ax.annotate(str(count), (i, count), xytext=(0, 4), textcoords='offset points',
        ha='center', fontsize=6)
for i, row in enumerate(profiles.itertuples()):
    flags = [row.bulk_hit, row.cluster_hit, row.framework_hit_common]
    for j, flag in enumerate(flags):
        matrix_ax.scatter(i, 2-j, s=18, color='0.2' if flag else '0.85', linewidths=0)
    on = np.flatnonzero(flags)
    if len(on) > 1:
        matrix_ax.plot([i, i], [2-on.max(), 2-on.min()], color='0.3', lw=.7, zorder=0)
matrix_ax.set(yticks=[2,1,0], yticklabels=['Whole PT', 'S1/S2/S3', 'Framework'],
    xticks=[], ylim=(-.5, 2.5))
for spine in matrix_ax.spines.values():
    spine.set_visible(False)
for i, label in enumerate(['Level only', 'Level + spatial', 'Spatial only']):
    framework_ax.barh(i, class_counts.loc[label, True], color='#8DA0CB')
    framework_ax.barh(i, class_counts.loc[label, False], left=class_counts.loc[label, True], color='0.8')
    framework_ax.text(class_counts.loc[label].sum()+2, i, str(class_counts.loc[label].sum()), va='center', fontsize=6)
framework_ax.set(yticks=[0,1,2], yticklabels=['Level only', 'Level + spatial', 'Spatial only'],
    xlabel='Pathway terms', title='Global and spatial calls')
framework_ax.invert_yaxis()
framework_ax.set_xlim(0, max(class_counts.iloc[:3].sum(axis=1))*1.25)
framework_ax.text(.5, -.19, f'Neither: {int(class_counts.loc["Neither"].sum())}; total-only: {int(class_counts.loc["Neither", True])}',
    transform=framework_ax.transAxes, fontsize=6, ha='center')
fig.text(.66, .93, 'Blue: total supported   Gray: no total call', fontsize=6)
for ax, letter in [(counts_ax,'a'), (matrix_ax,'b'), (framework_ax,'c')]:
    panel_label(ax, letter)
export_figure(fig, 'figure_4_information_profiles',
    'Figure 4 | Progressively finer representation enables overlapping information profiles. '
    'a,b, Exact UpSet intersections for native whole-PT signed GSEA, any native segment signed GSEA, '
    'and any common-universe framework evidence (positive matched level/spatial/total, signed level GSEA or smooth segment GSEA). '
    'Counts include all common pathways; the count axis is symlog to retain both the background and small profiles. '
    'Each component retains notebook 12\'s own full-family q≤0.05 convention; this union is descriptive, not jointly FDR-controlled. '
    'c, Positive matched level/spatial classes, with total support shown independently; neither and total-only counts are printed. '
    'Detection profiles do not imply nested significance sets or universal power superiority. '
    'Pathway terms overlap in genes and are not independent programs. One human donor; exploratory.')

# %% [markdown]
# ## Figure 5 · Supplementary validation controls
#
# These panels reassure reviewers without replacing the biological narrative. Signed level agreement and unsigned pathway recovery assess retention of conventional signal. The cluster-only diagnostic separates support elsewhere from thresholds close to 0.05. The controlled benchmark is **deliberately matched**, with 66/66 discrete hits retained plus 23 additional calls in the saved run; it is not the native conventional comparator used to choose Figures 2–3. Numbers are recomputed from saved source tables.

# %%
fig, axes = plt.subplots(2, 2, figsize=(FIGURE_WIDTH_IN, 140/25.4))
fig.subplots_adjust(left=.22, right=.97, bottom=.11, top=.9, wspace=.9, hspace=.65)
broad_gene_scatter(axes[0,0])
axes[0,0].set_title(f'Signed broad-signal agreement\nρ = {rho:.3f}; directions = {direction_agreement:.1%}', fontsize=7, pad=14)
signed_row = global_summary.set_index('comparison').loc['signed GSEA']
x = classification.common_effect_T_whole_pt_deseq2_abs
y = classification.common_effect_T_level
bulk_unsigned = classification.common_effect_T_whole_pt_deseq2_abs.gt(0) & classification.common_q_empirical_T_whole_pt_deseq2_abs.le(ALPHA)
level_unsigned = classification.common_effect_T_level.gt(0) & classification.common_q_empirical_T_level.le(ALPHA)
axes[0,1].scatter(x, y, s=4, color='0.75', alpha=.4, rasterized=True, linewidths=0)
axes[0,1].scatter(x[bulk_unsigned], y[bulk_unsigned], s=16,
    color=np.where(level_unsigned[bulk_unsigned], '#0072B2', '#D55E00'), linewidths=0)
axes[0,1].set(xlabel='Bulk magnitude AUC effect', ylabel='GAM-level AUC effect',
    title=f'Unsigned recovery: {int((bulk_unsigned & level_unsigned).sum())}/{int(bulk_unsigned.sum())}\nEffect-rank ρ = {x.corr(y, method="spearman"):.3f}')
fate_labels = ['Significant elsewhere', 'Near: q < 0.10', 'Near: 0.10 ≤ q < 0.15', 'No q < 0.15']
fate_counts = [int(cluster_review.any_framework_hit.sum()), int(cluster_review.near_threshold_010.sum()),
    int((~cluster_review.any_framework_hit & cluster_review.framework_q_lt_015 & ~cluster_review.framework_q_lt_010).sum()),
    int((~cluster_review.any_framework_hit & ~cluster_review.framework_q_lt_015).sum())]
axes[1,0].barh(range(4), fate_counts, color=['#8DA0CB', '#E5C494', '#E69F00', '0.75'])
axes[1,0].set(yticks=range(4), yticklabels=fate_labels, xlabel='Pathway terms',
    title=f'{len(cluster_review)} cluster signed-GSEA-only terms')
axes[1,0].invert_yaxis()
axes[1,0].set_xlim(0, max(fate_counts)*1.2)
for i, count in enumerate(fate_counts):
    axes[1,0].text(count+.7, i, str(count), va='center', fontsize=6)
control = controlled_summary.set_index('comparison').loc['controlled_discrete_to_continuous']
axes[1,1].barh([0,1], [control.baseline_hits, control.overlap], color=['0.7', '#8DA0CB'])
axes[1,1].barh(1, control.continuous_added, left=control.overlap, color='#D55E00')
axes[1,1].set(yticks=[0,1], yticklabels=['Matched discrete', 'Matched continuous'], xlabel='Pathway terms',
    title='Controlled model benchmark')
axes[1,1].invert_yaxis()
axes[1,1].text(.5, .5, f'{int(control.overlap)}/{int(control.baseline_hits)} retained + {int(control.continuous_added)}',
    transform=axes[1,1].transAxes, ha='center', va='center', fontsize=6)
axes[1,1].text(.5, -.29, 'Matched control; not the native baseline',
    transform=axes[1,1].transAxes, ha='center', fontsize=6)
for ax, letter in zip(axes.flat, 'abcd'):
    panel_label(ax, letter)
pd.DataFrame({'fate':fate_labels, 'pathways':fate_counts}).to_csv(source_data / 'figure_5_cluster_only_fate.csv', index=False)
controlled_summary.to_csv(source_data / 'figure_5_controlled_benchmark.csv', index=False)
global_summary.to_csv(source_data / 'figure_5_global_agreement.csv', index=False)
export_figure(fig, 'figure_5_supporting_validation',
    'Figure 5 | Supporting validation, not the central biological comparison. '
    'a, Common-gene signed agreement (same source as Figure 1); HC3 Z remains structure-conditional. '
    f'Signed pathway GSEA overlap: {int(signed_row.shared)} shared, {int(signed_row.DESeq2_only)} DESeq2-only, {int(signed_row.GAM_level_only)} GAM-level-only. '
    'b, Whole-PT absolute DESeq2 Wald matched-AUC effects versus unsigned T_level effects, same common genes/sets/null; '
    'blue denotes conventional magnitude calls recovered by level, orange denotes conventional-only. '
    'c, Cluster-only in the specific cluster-versus-smooth signed segment screen; support elsewhere uses original/common matched tests, signed GAM level and smooth segment screens. '
    'Near q<0.10/q<0.15 categories are exclusive and require no significant framework call. '
    'd, Section 9c deliberately matched discrete/continuous WLS total-difference enrichment, with the same expression and statistical machinery. '
    'This benchmark is a robustness check, not the conventional baseline for biological examples. '
    'All q conventions are inherited unchanged from notebook 12; no intervals are added to descriptive term counts. One human donor; exploratory.')

# %% [markdown]
# ## 7 · Manuscript-ready interpretation and reproducibility record
#
# Report recovery, shared localization examples, and spatially structured native-screen-negative examples separately. “Missed” below means not passing the saved conventional pathway screens, not a demonstrated absence of a conventional gene effect. Compression measurements support attenuation/cancellation descriptions while differing estimators remain a possible source of screen discordance. The single human donor limits population inference.
#
# The record binds candidate decisions, exact inputs and figure code to exports. All plotted measurements are in tidy source-data tables, including undisplayed members and all pathway information profiles. Figure legends accompany each PDF/SVG/PNG. Before submission, biologically review the annotations, gene directions, support, cohort caveats and final-size panel readability.

# %%
figure_index = pd.DataFrame(figure_records)
figure_index.to_csv(output / 'figure_index.csv', index=False)
input_files = [source / name for name in required_files] + [coordinate_path]
record = {'figure_logic': FIGURE_LOGIC_VERSION, 'notebook12_logic': manifest12['logic_version'],
    'fit_input_fingerprint': manifest12['fit_input_fingerprint'],
    'inputs': {path.name: digest(path) for path in input_files},
    'helper_code': digest(Path(resolution_story.__file__)),
    'figure_notebook_code': digest(project / 'analysis' / 'notebooks' / '14_pt_paper_figures_resolution_story.ipynb'),
    'common_support': [float(lo), float(hi)], 'axis': 'saved support mapped to [0,1] for display only',
    'q_threshold': ALPHA, 'minimum_all_planned_retention': MIN_RETENTION,
    'selected': selection[['pathway_id', 'figure_class', 'selection_rationale']].to_dict('records'),
    'native_baseline': 'whole-PT and reviewed-segment DESeq2 signed GSEA; extra magnitude/ORA exclusions for Class B',
    'uncertainty': 'pointwise HC3; averaged residuals retain cross-gene covariance; conditional on structures and coordinate',
    'guides': 'midpoints of equal-specimen segment median positions; labels overlap',
    'qa_available': qa_available, 'n_common_genes': len(common_genes), 'n_common_pathways': len(classification)}
(output / 'figure_manifest.json').write_text(json.dumps(record, indent=2))
print('MANUSCRIPT SUMMARY — exploratory, one human donor')
print(f'Broad gene signal retained: Spearman ρ={rho:.3f}; direction agreement={direction_agreement:.1%}; '
      f'unsigned pathway recovery={int(unsigned_row.shared)}/{int(unsigned_row.DESeq2_hits)}.')
for row in selection.itertuples():
    print(f'Class {row.figure_class}: {row.display_name}; native bulk q={row.q_family_DESeq2:.3g}, '
          f'best segment q={row.cluster_best_q:.3g}, spatial q={row.common_q_empirical_T_spatial:.3g}, '
          f'all-planned retention={row.robust_fraction_all_planned:.0%}. {row.selection_rationale}')
print('Produced: Figure 1 a–b; Figure 2 a–f + densities; Figure 3 a–l + densities; '
      'Figure 4 a–c; Figure 5 a–d. PDF, SVG, PNG and tidy source data accompany every figure.')
manuscript_text = ('Whole-PT and S1/S2/S3 analyses captured broad and coarse regional species differences, '
    'while pseudostructure reconstruction retained broad signals and resolved their continuous spatial organization. '
    'Robust spatial pathway calls absent from the native conventional screens showed changing, heterogeneous or '
    'opposing fitted effects that coarse anatomical averages could not describe; estimator differences and limited '
    'donor replication remain relevant to interpretation.')
print(manuscript_text)
(output / 'manuscript_interpretation.txt').write_text(manuscript_text + '\n')
