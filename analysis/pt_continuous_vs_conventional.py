# %% [markdown]
# # Continuous versus conventional PT comparison
#
# **Question.** Which gene and pathway patterns in these mouse and human PT datasets are added by a continuous shared-coordinate comparison beyond bulk PT or three broad PT bins?
#
# All methods use notebook 03's PT structures and `lognorm` layer, notebook 09's exact eligible one-to-one ortholog list, and notebook 09's Reactome, Hallmark, and KEGG pathway membership. Fixed-coordinate Jeffreys mismatch from notebook 09 is the continuous comparator; G2G is a secondary alignment annotation. The two human sections, including the one named `HUK1_MED1`, are healthy **cortex from one donor**. All results are descriptive; pathway q-values test competitive gene rankings within these datasets and do not represent donor-level species inference.
#
# Run notebooks 03 and 09 first. Notebook 07 provides optional P1–P5 labels. This notebook does not recompute an embedding, ortholog mapping, G2G alignment, or pathway library.

# %% [markdown]
# ## 1. Setup, inputs, and comparison contract

# %%
import argparse
import ast
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.stats import mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests
from IPython.display import display

def project_dir():
    starts = [Path.cwd().resolve()]
    if '__file__' in globals():
        starts.insert(0, Path(__file__).resolve().parent)
    for start in starts:
        for candidate in (start, *start.parents):
            if (candidate / 'pseudospace').is_dir() and (candidate / 'data').is_dir():
                return candidate
    raise RuntimeError('Cannot locate project root.')

PROJECT_DIR = project_dir()
if str(PROJECT_DIR) not in sys.path:
    sys.path.insert(0, str(PROJECT_DIR))
parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--data-root', type=Path)
parser.add_argument('--results-root', type=Path)
args, _ = parser.parse_known_args()
DATA_ROOT = Path(args.data_root or os.environ.get('PSEUDOSPACE_DATA_ROOT') or PROJECT_DIR / 'data').expanduser().resolve()
RESULTS_ROOT = Path(args.results_root or os.environ.get('PSEUDOSPACE_RESULTS_ROOT') or PROJECT_DIR / 'results').expanduser().resolve()
UPSTREAM = RESULTS_ROOT / 'human_vs_healthy_mouse'
CONTINUOUS = RESULTS_ROOT / 'pt_genes2genes'
OUTPUT = RESULTS_ROOT / 'pt_continuous_vs_conventional'
SPECIMENS = {'mouse': ('Ctrl1A2', 'Ctrl1A4'), 'human': ('HUK1_COR1', 'HUK1_MED1')}
SEGMENTS = ('S1_like', 'S2_like', 'S3_like')
COLORS = {'mouse': '#0072B2', 'human': '#D55E00'}
OUTPUT.mkdir(parents=True, exist_ok=True)

def savefig(fig, name):
    fig.tight_layout()
    fig.savefig(OUTPUT / name, dpi=180, bbox_inches='tight')
    plt.show()
    plt.close(fig)

def read_required(name):
    path = CONTINUOUS / name
    if not path.exists():
        raise FileNotFoundError(f'{path} is required; finish notebook 09 before running notebook 10.')
    return pd.read_csv(path)

def unique_genes(frame, label):
    if 'gene' not in frame or frame.gene.isna().any() or frame.gene.duplicated().any():
        raise ValueError(f'{label} needs one non-null row per gene.')
    return frame


# %%
universe = unique_genes(read_required('g2g_gene_universe.csv'), 'eligible universe')
genes = universe.loc[universe.g2g_eligible.astype(bool), 'gene'].tolist()
fixed = unique_genes(read_required('fixed_pseudospace_gene_statistics.csv'), 'fixed statistics')
fixed_pairs = unique_genes(read_required('fixed_specimen_pair_consistency.csv'), 'fixed pairs')
g2g = unique_genes(read_required('g2g_vs_fixed_gene_comparison.csv'), 'G2G annotations')
required_fixed = {'fixed_overall_mismatch', 'fixed_early_mismatch', 'fixed_mid_mismatch',
                  'fixed_late_mismatch', 'min_specimen_detected_fraction',
                  'mouse_raw_decile_amplitude_lognorm', 'human_raw_decile_amplitude_lognorm'}
if not required_fixed.issubset(fixed):
    raise ValueError(f'Notebook 09 fixed statistic lacks {required_fixed - set(fixed)}.')
if set(genes) != set(fixed.gene) or set(genes) != set(fixed_pairs.gene) or set(genes) != set(g2g.gene):
    raise ValueError('Notebook 09 files disagree on the final eligible gene universe.')
adata = sc.read_h5ad(UPSTREAM / 'cross_species_pt_dpt.h5ad')
need_obs = {'comparison_species', 'sample', 'broad_tubule_marker_call', 'shared_pseudospace'}
if not need_obs.issubset(adata.obs) or 'lognorm' not in adata.layers:
    raise ValueError('Notebook 03 PT artifact lacks the expected labels, coordinate, or lognorm layer.')
adata = adata[adata.obs.broad_tubule_marker_call.astype(str).eq('PT'), genes].copy()
actual = {species: tuple(sorted(adata.obs.loc[adata.obs.comparison_species.astype(str).eq(species),
                                              'sample'].astype(str).unique())) for species in SPECIMENS}
if actual != SPECIMENS or not adata.obs_names.is_unique or not adata.var_names.is_unique:
    raise ValueError(f'PT structure/specimen contract failed: {actual}.')
X = adata.layers['lognorm']
X = X.tocsr() if sparse.issparse(X) else sparse.csr_matrix(X)
if not np.isfinite(X.data).all() or (X.data < 0).any():
    raise ValueError('lognorm contains non-finite or negative values.')
s = adata.obs.shared_pseudospace.to_numpy(float)
if not np.isfinite(s).all() or np.unique(s).size < 15:
    raise ValueError('Invalid shared PT coordinate.')
if s.min() < 0 or s.max() > 1:
    s = (s - s.min()) / (s.max() - s.min())
adata.obs['shared_s_09'] = s
fixed_grid = read_required('fixed_grid.csv').fixed_grid.to_numpy(float)
lo, hi = fixed_grid[0], fixed_grid[-1]
expected_support = (max(np.quantile(s[adata.obs.comparison_species.astype(str).eq('mouse')], .05),
                        np.quantile(s[adata.obs.comparison_species.astype(str).eq('human')], .05)),
                    min(np.quantile(s[adata.obs.comparison_species.astype(str).eq('mouse')], .95),
                        np.quantile(s[adata.obs.comparison_species.astype(str).eq('human')], .95)))
if not np.allclose((lo, hi), expected_support, atol=1e-6):
    raise ValueError('Notebook 09 fixed grid does not match this PT coordinate and structure set.')
print(f'{len(genes):,} eligible genes; {adata.n_obs:,} PT structures; fixed support [{lo:.3f}, {hi:.3f}].')
display(adata.obs.groupby(['comparison_species', 'sample'], observed=True).size().rename('structures').reset_index())

# %% [markdown]
# ## 2. Bulk and three-bin PT effects
#
# Each species mean is the equal average of its two specimen or section means; structure counts do not reweight the species. The discrete bins are fixed thirds of notebook 09's common-support interval, labeled **S1-like**, **S2-like**, **S3-like** because the human S2 marker validation is weak. Both methods retain all PT structures: structures beyond common support enter the nearest terminal bin. The continuous statistic evaluates the central common support using those structures' kernel weights. No boundaries are tuned to species differences.

# %%
edges = lo + (hi - lo) * np.array([1/3, 2/3])
segment_index = np.searchsorted(edges, s, side='right')
segment_table = pd.DataFrame({'structure': adata.obs_names, 'sample': adata.obs['sample'].astype(str).to_numpy(),
                              'species': adata.obs.comparison_species.astype(str).to_numpy(),
                              'shared_pseudospace': s, 'segment': np.array(SEGMENTS)[segment_index]})
segment_table.to_csv(OUTPUT / 'pt_structure_segment_assignments.csv', index=False)
segment_counts = segment_table.groupby(['sample', 'segment']).size().unstack(fill_value=0).reindex(columns=SEGMENTS)
if (segment_counts < 20).any().any():
    raise ValueError(f'At least one sample has fewer than 20 PT structures in a bin:\n{segment_counts}')
display(segment_counts)

def sample_summary(mask):
    selected = X[mask]
    return np.asarray(selected.mean(axis=0)).ravel(), np.asarray((selected > 0).mean(axis=0)).ravel()

sample_means, sample_detection = {}, {}
for species, samples in SPECIMENS.items():
    for sample in samples:
        mask = segment_table['sample'].eq(sample).to_numpy()
        sample_means[sample], sample_detection[sample] = sample_summary(mask)
bulk = pd.DataFrame({'gene': genes})
for sample in (*SPECIMENS['mouse'], *SPECIMENS['human']):
    bulk[f'{sample}_mean'] = sample_means[sample]
    bulk[f'{sample}_detection'] = sample_detection[sample]
for species, samples in SPECIMENS.items():
    bulk[f'bulk_{species}_mean'] = np.mean([sample_means[x] for x in samples], axis=0)
    bulk[f'detection_{species}'] = np.mean([sample_detection[x] for x in samples], axis=0)
bulk['bulk_effect'] = bulk.bulk_human_mean - bulk.bulk_mouse_mean
bulk['bulk_abs_effect'] = bulk.bulk_effect.abs()
bulk.to_csv(OUTPUT / 'bulk_pt_gene_statistics.csv', index=False)

segment = pd.DataFrame({'gene': genes})
for label in SEGMENTS:
    local = {}
    for species, samples in SPECIMENS.items():
        for sample in samples:
            mask = (segment_table['sample'].eq(sample) & segment_table.segment.eq(label)).to_numpy()
            local[sample], _ = sample_summary(mask)
            segment[f'{label}_{sample}_mean'] = local[sample]
    segment[f'{label}_mouse_mean'] = np.mean([local[x] for x in SPECIMENS['mouse']], axis=0)
    segment[f'{label}_human_mean'] = np.mean([local[x] for x in SPECIMENS['human']], axis=0)
    segment[f'{label}_effect'] = segment[f'{label}_human_mean'] - segment[f'{label}_mouse_mean']
effects = segment[[f'{x}_effect' for x in SEGMENTS]].to_numpy(float)
segment['segment_max_abs_effect'] = np.abs(effects).max(axis=1)
segment['segment_rms_effect'] = np.sqrt(np.square(effects).mean(axis=1))
segment['segment_of_max_effect'] = np.array(SEGMENTS)[np.abs(effects).argmax(axis=1)]
segment.to_csv(OUTPUT / 'segment_pt_gene_statistics.csv', index=False)
# Small executable contract: balanced effects equal arithmetic means of specimen effects.
assert np.allclose(bulk.bulk_effect, ((bulk.HUK1_COR1_mean + bulk.HUK1_MED1_mean) -
                                     (bulk.Ctrl1A2_mean + bulk.Ctrl1A4_mean)) / 2)
assert np.allclose(segment.segment_max_abs_effect, np.abs(effects).max(axis=1))

# %% [markdown]
# ## 3. Gene rankings, overlap, and descriptive candidates

# %%
annotation_cols = ['gene', 'alignment_similarity', 'g2g_warp_fraction', 'g2g_mismatch_fraction',
                   'g2g_cluster', 'mouse_positional_module', 'pairwise_median_similarity',
                   'pairwise_range_similarity']
gene_table = bulk.merge(segment, on='gene', validate='one_to_one').merge(
    fixed, on='gene', validate='one_to_one').merge(
    fixed_pairs.drop(columns='fixed_overall_mismatch'), on='gene', validate='one_to_one').merge(
    g2g[annotation_cols], on='gene', validate='one_to_one')
# Notebook 09's Jeffreys distance combines changes in local mean and variance.
# Decompose its saved local distributions so a variance effect is not called a shifted mean trajectory.
with np.load(CONTINUOUS / 'fixed_local_distributions.npz') as local:
    local_genes = local['genes'].astype(str).tolist()
    local_grid = local['fixed_grid'].copy()
    mu_mouse, mu_human = local['mu_mouse'].copy(), local['mu_human'].copy()
    sd_mouse, sd_human = local['sd_mouse'].copy(), local['sd_human'].copy()
if local_genes != genes or not np.allclose(local_grid, fixed_grid):
    raise ValueError('Notebook 09 local distributions disagree with this gene universe or grid.')
vm, vh = np.maximum(sd_mouse**2, .05**2), np.maximum(sd_human**2, .05**2)
mean_cost = .25 * (mu_mouse-mu_human)**2 * (1/vm+1/vh)
variance_cost = .25 * (vm/vh+vh/vm-2)
gene_table['fixed_mean_component'] = mean_cost.mean(axis=0)
gene_table['fixed_variance_component'] = variance_cost.mean(axis=0)
gene_table['fixed_mean_component_fraction'] = gene_table.fixed_mean_component / gene_table.fixed_overall_mismatch.clip(lower=1e-12)
if not np.allclose(gene_table.fixed_mean_component + gene_table.fixed_variance_component,
                   gene_table.fixed_overall_mismatch, rtol=1e-5, atol=1e-7):
    raise ValueError('Fixed-score decomposition does not reproduce notebook 09.')
local_difference = mu_human - mu_mouse
mouse_peak = fixed_grid[mu_mouse.argmax(axis=0)]
human_peak = fixed_grid[mu_human.argmax(axis=0)]
gene_table['max_local_mean_difference'] = np.abs(local_difference).max(axis=0)
gene_table['opposing_local_mean_directions'] = ((local_difference.min(axis=0) < -.05) &
    (local_difference.max(axis=0) > .05))
gene_table['within_same_bin_peak_shift'] = (
    (np.searchsorted(edges, mouse_peak) == np.searchsorted(edges, human_peak)) &
    (np.abs(mouse_peak-human_peak) >= np.diff(fixed_grid).mean()))
gene_table['mismatch_near_bin_boundary'] = np.min(np.abs(
    gene_table.fixed_max_mismatch_position.to_numpy()[:,None] - edges[None,:]),axis=1) <= np.diff(fixed_grid).mean()
score_cols = {'bulk': 'bulk_abs_effect', 'segment': 'segment_max_abs_effect',
              'continuous': 'fixed_overall_mismatch'}
for method, score in score_cols.items():
    gene_table[f'{method}_rank'] = gene_table[score].rank(ascending=False, method='average')
    gene_table[f'{method}_ordinal'] = gene_table[score].rank(ascending=False, method='first').astype(int)
gene_table.to_csv(OUTPUT / 'gene_method_rank_comparison.csv', index=False)
rank_rows = []
for a, b in (('bulk', 'segment'), ('bulk', 'continuous'), ('segment', 'continuous')):
    rho = spearmanr(gene_table[f'{a}_rank'], gene_table[f'{b}_rank']).statistic
    for n in (100, 250, 500, 1000):
        n = min(n, len(gene_table))
        left = set(gene_table.loc[gene_table[f'{a}_ordinal'].le(n), 'gene'])
        right = set(gene_table.loc[gene_table[f'{b}_ordinal'].le(n), 'gene'])
        overlap = len(left & right)
        rank_rows.append({'method_a': a, 'method_b': b, 'top_n': n,
                          'spearman': rho, 'overlap': overlap, 'jaccard': overlap / (2*n-overlap)})
rank_summary = pd.DataFrame(rank_rows)
rank_summary.to_csv(OUTPUT / 'gene_method_rank_summary.csv', index=False)
display(rank_summary)

# %%
# Equal-size top sets, with sensitivity. Ties are ordered by the stable eligible-gene order.
threshold_rows = []
for fraction in (.05, .10, .20):
    n = max(1, round(len(gene_table) * fraction))
    flags = {m: gene_table[f'{m}_ordinal'].le(n) for m in score_cols}
    for name, flag in flags.items():
        gene_table[f'top_{name}_{int(fraction*100)}pct'] = flag
    trajectory_only = flags['continuous'] & ~flags['bulk'] & ~flags['segment']
    threshold_rows.append({'top_fraction': fraction, 'top_n': n,
        'continuous_not_bulk': int((flags['continuous'] & ~flags['bulk']).sum()),
        'continuous_not_segment': int((flags['continuous'] & ~flags['segment']).sum()),
        'continuous_only': int(trajectory_only.sum()),
        'continuous_only_min_pair_q80': int((trajectory_only & gene_table.pairwise_min_divergence_percentile.ge(.8)).sum())})
threshold_sensitivity = pd.DataFrame(threshold_rows)
threshold_sensitivity.to_csv(OUTPUT / 'gene_threshold_sensitivity.csv', index=False)
strong = {m: gene_table[f'top_{m}_10pct'] for m in score_cols}
gene_table['discovery_category'] = np.select([
    strong['bulk'] & strong['segment'] & strong['continuous'],
    strong['bulk'] & ~strong['segment'] & strong['continuous'],
    ~strong['bulk'] & strong['segment'] & strong['continuous'],
    ~strong['bulk'] & ~strong['segment'] & strong['continuous'],
    strong['bulk'] & ~strong['segment'] & ~strong['continuous'],
    ~strong['bulk'] & strong['segment'] & ~strong['continuous'],
    strong['bulk'] & strong['segment'] & ~strong['continuous']],
    ['shared_all', 'bulk_continuous', 'segment_continuous', 'continuous_only',
     'bulk_only', 'segment_only', 'bulk_segment_only'], default='below_top_decile')
gene_table.to_csv(OUTPUT / 'gene_method_discovery_categories.csv', index=False)
# Inspect thresholds before interpreting candidate counts; all limits are rank/observed-QC guides.
thresholds = {'bulk_median': gene_table.bulk_abs_effect.median(),
              'segment_median': gene_table.segment_max_abs_effect.median(),
              'min_detection': .05, 'min_amplitude': .1, 'pair_min_percentile': .8}
display(pd.DataFrame({x: gene_table[x].quantile([.1,.25,.5,.75,.8,.9,.95]) for x in
    ['bulk_abs_effect','segment_max_abs_effect','fixed_overall_mismatch',
     'min_specimen_detected_fraction','pairwise_min_divergence_percentile',
     'mouse_raw_decile_amplitude_lognorm','human_raw_decile_amplitude_lognorm']}))
display(threshold_sensitivity)

# %%
quality = (gene_table.min_specimen_detected_fraction.ge(thresholds['min_detection']) &
           gene_table[['mouse_raw_decile_amplitude_lognorm',
                       'human_raw_decile_amplitude_lognorm']].max(axis=1).ge(thresholds['min_amplitude']))
robust = gene_table.pairwise_min_divergence_percentile.ge(thresholds['pair_min_percentile'])
trajectory = gene_table.loc[strong['continuous'] & ~strong['bulk'] & ~strong['segment'] &
    gene_table.bulk_abs_effect.le(thresholds['bulk_median']) &
    gene_table.segment_max_abs_effect.le(thresholds['segment_median']) & quality & robust].copy()
bulk_only = gene_table.loc[strong['bulk'] & ~strong['continuous'] &
    gene_table.continuous_ordinal.gt(.5*len(gene_table)) & quality].copy()
segment_only = gene_table.loc[strong['segment'] & ~strong['continuous'] &
    gene_table.continuous_ordinal.gt(.5*len(gene_table)) & quality].copy()
for frame, name, sort in [(trajectory, 'trajectory_only_candidate_genes.csv', 'fixed_overall_mismatch'),
                          (bulk_only, 'bulk_only_candidate_genes.csv', 'bulk_abs_effect'),
                          (segment_only, 'segment_only_candidate_genes.csv', 'segment_max_abs_effect')]:
    frame.sort_values(sort, ascending=False).to_csv(OUTPUT / name, index=False)
print({'trajectory_only': len(trajectory), 'bulk_only': len(bulk_only), 'segment_only': len(segment_only)})
# Binning misses require small *signed bin-mean* effects, not merely a low rank.
profile = read_required('fixed_pseudospace_mismatch_profiles.csv').set_index('gene').loc[genes]
profile_values = profile.to_numpy(float)
profile_slope = np.abs(np.diff(profile_values, axis=1)).mean(axis=1)
gene_table['mismatch_profile_step'] = profile_slope
# A continuous-only candidate is reported with its peak location and nearby boundary distance;
# these are diagnostics, not asserted biological mechanism labels.
gene_table['distance_to_nearest_bin_boundary'] = np.min(np.abs(
    gene_table.fixed_max_mismatch_position.to_numpy()[:, None] - edges[None, :]), axis=1)
missed_genes = gene_table.loc[gene_table.gene.isin(trajectory.gene)].copy()
missed_genes.to_csv(OUTPUT / 'discrete_binning_missed_genes.csv', index=False)

# %% [markdown]
# ## 4. Pathways on the same membership and gene universe
#
# The signed conventional statistics are human minus mouse lognorm effects. For each pathway, a two-sided Mann–Whitney test compares member gene scores with the remaining eligible genes; signed AUC−0.5 records direction. Each segment keeps its own test and BH correction covers all libraries and all three segment tests together. The fixed-coordinate AUC and q-values are imported unchanged from notebook 09. Since fixed mismatch is directionless, pathway ranking uses positive divergence effect; conventional ranking uses absolute signed effect. All q-values are competitive within-dataset summaries.

# %%
coverage = read_required('g2g_pathway_coverage.csv')
pathways = coverage.loc[coverage.retained.astype(bool), ['library','pathway','n_genes_present','genes_present']].rename(columns={'n_genes_present':'n_eligible_genes'}).copy()
pathways['members'] = pathways.genes_present.map(ast.literal_eval)
if pathways.duplicated(['library','pathway']).any() or pathways.empty:
    raise ValueError('Unexpected pathway membership table.')
if any(set(row.members) - set(genes) or len(row.members) != row.n_eligible_genes for row in pathways.itertuples()):
    raise ValueError('Pathway membership is inconsistent with notebook 09 eligible genes.')
fixed_ranked = read_required('fixed_ranked_pathway_enrichment.csv')
fixed_positional = read_required('fixed_positional_pathway_enrichment.csv')
fixed_coherence = read_required('fixed_pathway_coherence.csv')
if set(map(tuple, pathways[['library','pathway']].to_numpy())) != set(map(tuple, fixed_ranked[['library','pathway']].to_numpy())):
    raise ValueError('Continuous and conventional pathway test families differ.')

def ranked_pathway(frame, score, segment_name):
    stats = frame.set_index('gene')[score]
    rows = []
    for row in pathways.itertuples():
        inside = stats.loc[row.members].to_numpy(float)
        outside = stats.loc[~stats.index.isin(row.members)].to_numpy(float)
        u, p = mannwhitneyu(inside, outside, alternative='two-sided')
        auc = u / (len(inside)*len(outside))
        rows.append({'library': row.library, 'pathway': row.pathway, 'segment': segment_name,
                     'n_eligible_genes': len(inside), 'rank_auc': auc, 'rank_effect': auc-.5,
                     'median_gene_effect': np.median(inside), 'mean_gene_effect': np.mean(inside),
                     'p_value': p})
    return pd.DataFrame(rows)

bulk_path = ranked_pathway(gene_table, 'bulk_effect', 'bulk')
bulk_path['q_value'] = multipletests(bulk_path.p_value, method='fdr_bh')[1]
bulk_path.drop(columns='segment').to_csv(OUTPUT / 'bulk_ranked_pathway_enrichment.csv', index=False)
segment_path = pd.concat([ranked_pathway(gene_table, f'{label}_effect', label)
                          for label in SEGMENTS], ignore_index=True)
segment_path['q_value'] = multipletests(segment_path.p_value, method='fdr_bh')[1]
segment_path.to_csv(OUTPUT / 'segment_ranked_pathway_enrichment.csv', index=False)

# %%
keys = ['library', 'pathway', 'n_eligible_genes']
path_table = pathways[keys].merge(bulk_path.drop(columns='segment').add_prefix('bulk_').rename(columns={
    'bulk_library':'library','bulk_pathway':'pathway','bulk_n_eligible_genes':'n_eligible_genes'}), on=keys, validate='one_to_one')
for label in SEGMENTS:
    part = segment_path.loc[segment_path.segment.eq(label)].drop(columns='segment').add_prefix(f'{label}_').rename(columns={
        f'{label}_library':'library',f'{label}_pathway':'pathway',f'{label}_n_eligible_genes':'n_eligible_genes'})
    path_table = path_table.merge(part, on=keys, validate='one_to_one')
path_table = path_table.merge(fixed_ranked.add_prefix('continuous_').rename(columns={
    'continuous_library':'library','continuous_pathway':'pathway',
    'continuous_n_eligible_genes':'n_eligible_genes'}), on=keys, validate='one_to_one')
for region in ('early','mid','late'):
    part = fixed_positional.loc[fixed_positional.region.eq(region), ['library','pathway','rank_effect','q_value']].rename(
        columns={'rank_effect':f'continuous_{region}_effect','q_value':f'continuous_{region}_q_value'})
    path_table = path_table.merge(part, on=['library','pathway'], validate='one_to_one')
path_table = path_table.merge(fixed_coherence[['library','pathway','fixed_profile_entropy',
    'fixed_profile_entropy_null_z','fixed_mismatch_iqr','dominant_fixed_region',
    'normalized_cluster_entropy']], on=['library','pathway'], validate='one_to_one')
# AUC effects are comparable in rank space; lognorm means and Jeffreys distance remain separate scales.
path_table['segment_best_abs_effect'] = path_table[[f'{x}_rank_effect' for x in SEGMENTS]].abs().max(axis=1)
path_table['bulk_abs_rank_effect'] = path_table.bulk_rank_effect.abs()
path_table['continuous_divergence_effect'] = path_table.continuous_rank_effect
for method, col in [('bulk','bulk_abs_rank_effect'),('segment','segment_best_abs_effect'),
                    ('continuous','continuous_divergence_effect')]:
    path_table[f'{method}_rank'] = path_table[col].rank(ascending=False, method='average')
    path_table[f'{method}_ordinal'] = path_table[col].rank(ascending=False, method='first').astype(int)
path_table.to_csv(OUTPUT / 'pathway_method_comparison.csv', index=False)

# %% [markdown]
# ## 5. Pathway categories, cancellation, and binning loss

# %%
path_sensitivity = []
for fraction in (.05, .10, .20):
    n = max(1, round(len(path_table)*fraction))
    flags = {m: path_table[f'{m}_ordinal'].le(n) for m in ('bulk','segment','continuous')}
    for name, flag in flags.items():
        path_table[f'top_{name}_{int(fraction*100)}pct'] = flag
    path_sensitivity.append({'top_fraction':fraction,'top_n':n,
        'continuous_not_bulk':int((flags['continuous'] & ~flags['bulk']).sum()),
        'continuous_not_segment':int((flags['continuous'] & ~flags['segment']).sum()),
        'continuous_only':int((flags['continuous'] & ~flags['bulk'] & ~flags['segment']).sum())})
path_sensitivity = pd.DataFrame(path_sensitivity)
path_sensitivity.to_csv(OUTPUT / 'pathway_threshold_sensitivity.csv', index=False)
b, t, c = (path_table[f'top_{m}_10pct'] for m in ('bulk','segment','continuous'))
path_table['discovery_category'] = np.select([
    b & t & c, b & ~t & c, ~b & t & c, ~b & ~t & c,
    b & ~t & ~c, ~b & t & ~c, b & t & ~c],
    ['conventional_continuous','bulk_continuous','segment_continuous','continuous_only',
     'bulk_only','segment_only','bulk_segment_only'], default='below_top_decile')
# Keep a separate flag for discordant representation despite an exclusive display category.
path_table['method_discordant'] = ((b | t) ^ c)
path_table.to_csv(OUTPUT / 'pathway_method_discovery_categories.csv', index=False)

stats = gene_table.set_index('gene')
member_rows = []
for row in pathways.itertuples():
    members = stats.loc[row.members]
    signed = members.bulk_effect.to_numpy(float)
    member_rows.append({'library':row.library,'pathway':row.pathway,
        'median_gene_continuous_mismatch':members.fixed_overall_mismatch.median(),
        'abs_mean_bulk_gene_effect':abs(signed.mean()),
        'median_abs_bulk_gene_effect':np.median(abs(signed)),
        'fraction_positive_bulk':np.mean(signed > 0),
        'fraction_negative_bulk':np.mean(signed < 0),
        'mean_early_mismatch':members.fixed_early_mismatch.mean(),
        'mean_mid_mismatch':members.fixed_mid_mismatch.mean(),
        'mean_late_mismatch':members.fixed_late_mismatch.mean(),
        'fraction_technically_robust':np.mean(members.pairwise_min_divergence_percentile >= .8),
        'fraction_adequately_detected':np.mean(members.min_specimen_detected_fraction >= .05),
        'fraction_trajectory_only_genes':np.mean(members.index.isin(trajectory.gene)),
        'median_fixed_mean_component_fraction':members.fixed_mean_component_fraction.median(),
        'fraction_opposing_local_mean_directions':members.opposing_local_mean_directions.mean(),
        'fraction_within_same_bin_peak_shift':members.within_same_bin_peak_shift.mean(),
        'fraction_mismatch_near_bin_boundary':members.mismatch_near_bin_boundary.mean()})
path_table = path_table.merge(pd.DataFrame(member_rows), on=['library','pathway'], validate='one_to_one')
# Fixed divergence can include normalized level shifts. Require low conventional ranks,
# adequate member support, and coherence before calling a pathway a follow-up candidate.
coherence_cut = path_table.fixed_profile_entropy_null_z.median()
path_quality = (path_table.n_eligible_genes.ge(10) &
    path_table.fraction_adequately_detected.ge(.5) &
    path_table.fraction_technically_robust.ge(.25) &
    path_table.fixed_profile_entropy_null_z.le(coherence_cut))
continuous_path = path_table.loc[c & ~b & ~t & path_quality &
    path_table.bulk_abs_rank_effect.le(path_table.bulk_abs_rank_effect.median()) &
    path_table.segment_best_abs_effect.le(path_table.segment_best_abs_effect.median())].copy()
continuous_path.sort_values('continuous_rank').to_csv(OUTPUT / 'continuous_only_pathway_candidates.csv', index=False)
# Cancellation is measured both as signed-gene cancellation and spatial loss.
cancellation = path_table.loc[c & path_table.abs_mean_bulk_gene_effect.le(
    path_table.abs_mean_bulk_gene_effect.median())].copy()
cancellation['opposing_bulk_signs'] = cancellation[['fraction_positive_bulk','fraction_negative_bulk']].min(axis=1)
cancellation.sort_values(['continuous_rank','opposing_bulk_signs'], ascending=[True,False]).to_csv(
    OUTPUT / 'bulk_cancellation_pathways.csv', index=False)
binning = path_table.loc[c & ~t & path_table.segment_best_abs_effect.le(
    path_table.segment_best_abs_effect.median())].copy()
binning.sort_values('continuous_rank').to_csv(OUTPUT / 'discrete_binning_missed_pathways.csv', index=False)
path_table.to_csv(OUTPUT / 'pathway_method_comparison.csv', index=False)
path_table.to_csv(OUTPUT / 'pathway_method_discovery_categories.csv', index=False)
print({'continuous_only_pathways':len(continuous_path),'bulk_cancellation_pathways':len(cancellation),
       'binning_missed_pathways':len(binning),'coherence_null_z_median':coherence_cut})
display(path_sensitivity)

# %% [markdown]
# ## 6. Selected pathway gene evidence
#
# The three highest-ranked eligible continuous-only pathways and three concordant pathways are selected from observed ranks. Related pathway names can overlap heavily, so this is a follow-up shortlist, not independent discoveries. `fixed_profile_entropy_null_z` below zero means positional classes are more concentrated than size-matched random sets from notebook 09. It does not prove every member changes coherently.

# %%
selected = pd.concat([continuous_path.sort_values('continuous_rank').head(3),
    path_table.loc[path_table.discovery_category.eq('conventional_continuous')].sort_values('continuous_rank').head(3)],
    ignore_index=True).drop_duplicates(['library','pathway']).copy()
selected.to_csv(OUTPUT / 'selected_pathway_themes.csv', index=False)
if selected.empty:
    print('No pathway passed the displayed selection; inspect pathway_method_comparison.csv directly.')
else:
    selected_keys = set(map(tuple, selected[['library','pathway']].to_numpy()))
    detail_rows = []
    for row in pathways.itertuples():
        if (row.library,row.pathway) not in selected_keys:
            continue
        detail = gene_table.loc[gene_table.gene.isin(row.members)].copy()
        detail.insert(0, 'pathway', row.pathway)
        detail.insert(0, 'library', row.library)
        detail_rows.append(detail)
    selected_details = pd.concat(detail_rows, ignore_index=True)
    selected_details.to_csv(OUTPUT / 'selected_pathway_gene_details.csv', index=False)
    display(selected[['library','pathway','discovery_category','continuous_rank','bulk_rank',
                      'segment_rank','fixed_profile_entropy_null_z','fraction_technically_robust']])
if selected.empty:
    pd.DataFrame(columns=['library','pathway','gene']).to_csv(OUTPUT / 'selected_pathway_gene_details.csv', index=False)
interpretation_rows = []
for row in selected.itertuples():
    members = next(p.members for p in pathways.itertuples() if p.library == row.library and p.pathway == row.pathway)
    member_stats = gene_table.loc[gene_table.gene.isin(members)]
    modules = member_stats.mouse_positional_module.value_counts()
    modules = modules[modules.index.isin([f'P{i}' for i in range(1,6)])]
    best_segment = max(SEGMENTS, key=lambda x:abs(getattr(row,f'{x}_rank_effect')))
    strongest_region = max(('early','mid','late'),key=lambda x:getattr(row,f'mean_{x}_mismatch'))
    interpretation_rows.append({'library':row.library,'pathway':row.pathway,
        'A_mouse_position':modules.index[0] if len(modules) else 'P-module unavailable',
        'A_assigned_member_fraction':modules.sum()/len(members),
        'B_bulk_rank':row.bulk_rank,'B_bulk_mean_gene_effect':row.bulk_mean_gene_effect,
        'C_best_segment':best_segment,'C_segment_rank':row.segment_rank,
        'D_continuous_rank':row.continuous_rank,'D_largest_mismatch_region':strongest_region,
        'D_median_mean_component_fraction':row.median_fixed_mean_component_fraction,
        'D_opposing_mean_direction_fraction':row.fraction_opposing_local_mean_directions,
        'D_within_bin_peak_shift_fraction':row.fraction_within_same_bin_peak_shift,
        'D_boundary_peak_fraction':row.fraction_mismatch_near_bin_boundary,
        'E_entropy_null_z':row.fixed_profile_entropy_null_z,
        'E_robust_member_fraction':row.fraction_technically_robust,
        'E_detection_member_fraction':row.fraction_adequately_detected})
interpretations = pd.DataFrame(interpretation_rows)
interpretations.to_csv(OUTPUT / 'selected_pathway_interpretation.csv',index=False)
display(interpretations)

# %% [markdown]
# ## 7. Figures: design, rankings, overlap, and pathway effects

# %%
fig, axes = plt.subplots(1, 3, figsize=(12, 2.6), sharey=True)
for ax, title in zip(axes, ('Bulk PT','Three PT bins','Continuous PT')):
    ax.plot([0,1],[.5,.5], color='black', lw=1)
    ax.set(xlim=(0,1),ylim=(0,1),xlabel='shared PT pseudospace',title=title)
    ax.set_yticks([])
axes[0].fill_between([0,1], .35,.65, alpha=.25, color='gray')
for left,right,label in zip([0,1/3,2/3],[1/3,2/3,1],['S1-like','S2-like','S3-like']):
    axes[1].axvspan(left,right,alpha=.2)
    axes[1].text((left+right)/2,.7,label,ha='center',fontsize=8)
x = np.linspace(0,1,100)
axes[2].plot(x,.25+.5*x,color=COLORS['mouse'],label='mouse')
axes[2].plot(x,.5+.25*np.sin(2*np.pi*x),color=COLORS['human'],label='human')
axes[2].legend(fontsize=7)
savefig(fig,'fig01_analysis_design.png')

fig, axes = plt.subplots(1,3,figsize=(12,3.5))
for ax,(a,b) in zip(axes,(('bulk','segment'),('bulk','continuous'),('segment','continuous'))):
    ax.scatter(gene_table[f'{a}_rank'],gene_table[f'{b}_rank'],s=2,alpha=.1)
    ax.set(xlabel=f'{a} gene rank',ylabel=f'{b} gene rank',xscale='log',yscale='log')
    ax.text(.03,.95,f'ρ={spearmanr(gene_table[f"{a}_rank"],gene_table[f"{b}_rank"]).statistic:.2f}',
            transform=ax.transAxes,va='top')
savefig(fig,'fig02_gene_rank_comparison.png')

fig, axes = plt.subplots(1,2,figsize=(11,3.4))
for ax, frame, label in ((axes[0],gene_table,'genes'),(axes[1],path_table,'pathways')):
    counts = frame.discovery_category.value_counts()
    names = ['shared_all','conventional_continuous','bulk_continuous','segment_continuous',
             'continuous_only','bulk_only','segment_only','bulk_segment_only']
    names = [name for name in names if name in counts]
    ax.barh(names,[counts[name] for name in names],color=['#6B7280' if name != 'continuous_only' else '#C24A4A' for name in names])
    ax.set(xlabel=f'{label} in top-decile method sets',title=label)
    ax.invert_yaxis()
savefig(fig,'fig03_method_overlap.png')

fig, axes = plt.subplots(1,2,figsize=(11,4))
for ax,xcol,label in ((axes[0],'bulk_abs_rank_effect','|bulk signed AUC − 0.5|'),
                      (axes[1],'segment_best_abs_effect','best segment |signed AUC − 0.5|')):
    for category, group in path_table.groupby('discovery_category'):
        ax.scatter(group[xcol],group.continuous_divergence_effect,s=8,alpha=.3,
                   label=category if ax is axes[0] else None)
    ax.set(xlabel=label,ylabel='fixed divergence AUC − 0.5')
    for row in continuous_path.sort_values('continuous_rank').head(4).itertuples():
        ax.annotate(row.pathway[:25],(getattr(row,xcol),row.continuous_divergence_effect),fontsize=6)
axes[0].legend(fontsize=6,loc='best')
savefig(fig,'fig05_pathway_method_comparison.png')

# %% [markdown]
# ## 8. Representative genes and selected pathway curves
#
# Figures show raw structures, notebook 09's fixed-kernel local means, balanced bulk means, three-bin means, and notebook 09's continuous mismatch profile. Selection is based on the analysis tables; empty categories are omitted. The local curves are weighted descriptive summaries, not independently fitted biological replicates.

# %%
gene_index = {gene: i for i,gene in enumerate(genes)}
examples = []
for label, mask in [
    ('shared_all',gene_table.discovery_category.eq('shared_all')),
    ('bulk_only',gene_table.gene.isin(bulk_only.gene)),
    ('segment_continuous',gene_table.discovery_category.eq('segment_continuous')),
    ('trajectory_only',gene_table.gene.isin(trajectory.gene)),
    ('warp_rescuable',strong['continuous'] & gene_table.alignment_similarity.ge(.8) & gene_table.g2g_warp_fraction.ge(.3)),
    ('nonalignable',strong['continuous'] & gene_table.alignment_similarity.lt(.5) & robust & quality)]:
    subset = gene_table.loc[mask & quality].sort_values('fixed_overall_mismatch',ascending=False)
    if not subset.empty:
        examples.append((label,subset.iloc[0].gene))
pd.DataFrame(examples,columns=['category','gene']).to_csv(OUTPUT / 'representative_genes.csv',index=False)
if examples:
    fig, axes = plt.subplots(len(examples),5,figsize=(17,2.6*len(examples)),squeeze=False)
    for row,(label,gene) in enumerate(examples):
        i = gene_index[gene]
        for species in ('mouse','human'):
            mask = segment_table.species.eq(species).to_numpy()
            xx = X[mask,i].toarray().ravel()
            axes[row,0].scatter(s[mask],xx,s=2,alpha=.05,color=COLORS[species],label=species)
        for species, curve in [('mouse',mu_mouse[:,i]),('human',mu_human[:,i])]:
            axes[row,1].plot(fixed_grid,curve,color=COLORS[species],label=species)
        record = gene_table.loc[gene_table.gene.eq(gene)].iloc[0]
        axes[row,2].bar(['mouse','human'],[record.bulk_mouse_mean,record.bulk_human_mean],
                        color=[COLORS['mouse'],COLORS['human']])
        for species,offset in [('mouse',-.17),('human',.17)]:
            values=[record[f'{x}_{species}_mean'] for x in SEGMENTS]
            axes[row,3].bar(np.arange(3)+offset,values,width=.34,color=COLORS[species])
        axes[row,3].set_xticks(range(3),['S1','S2','S3'])
        axes[row,4].plot(fixed_grid,profile.loc[gene].to_numpy(float),color='#6A3D9A')
        for ax,title in zip(axes[row],('raw structures','fixed local means','bulk means','three-bin means','fixed mismatch')):
            ax.set_title(f'{gene}: {label}\n{title}',fontsize=8)
        axes[row,0].set_xlabel('shared s'); axes[row,1].set_xlabel('shared s')
        axes[row,4].set_xlabel('shared s')
        axes[row,0].set_ylabel('lognorm')
    axes[0,1].legend(fontsize=7)
    savefig(fig,'fig04_representative_genes.png')

# %%
if not selected.empty:
    pathway_gene_rows = []
    for row in selected.head(3).itertuples():
        members = next(p.members for p in pathways.itertuples() if p.library == row.library and p.pathway == row.pathway)
        top = gene_table.loc[gene_table.gene.isin(members)].sort_values('fixed_overall_mismatch',ascending=False).head(3)
        for gene in top.gene:
            pathway_gene_rows.append((row.pathway,gene))
    fig, axes = plt.subplots(len(pathway_gene_rows),2,figsize=(10,2.2*len(pathway_gene_rows)),squeeze=False)
    for j,(pathway,gene) in enumerate(pathway_gene_rows):
        i = gene_index[gene]
        axes[j,0].plot(fixed_grid,mu_mouse[:,i],color=COLORS['mouse'],label='mouse')
        axes[j,0].plot(fixed_grid,mu_human[:,i],color=COLORS['human'],label='human')
        axes[j,0].set(title=f'{pathway[:40]}: {gene}',xlabel='shared s',ylabel='local lognorm')
        record=gene_table.loc[gene_table.gene.eq(gene)].iloc[0]
        axes[j,1].bar(['bulk','S1','S2','S3'],[record.bulk_effect]+[record[f'{x}_effect'] for x in SEGMENTS],color='#7A8A99')
        axes[j,1].axhline(0,color='black',lw=.6)
        axes[j,1].set(ylabel='human − mouse lognorm')
    axes[0,0].legend(fontsize=7)
    savefig(fig,'fig06_selected_pathway_gene_curves.png')

# %% [markdown]
# ## 9. Data-driven synthesis
#
# A candidate is a follow-up priority, not a novel species effect by definition. Notebook 09's Jeffreys mismatch contains local mean and variance components, which the tables separate. A large variance component alone does not establish a shifted mean trajectory. In particular, a fixed divergence can reflect normalized level differences, and a low entropy pathway can still have opposing gene directions. The table below answers the requested comparison questions using observed ranks, technical pair checks, and member curves. Human sections remain one donor; these findings require independent donors and orthogonal validation.

# %%
summary = pd.DataFrame([
    {'question':'Bulk vs continuous gene rank agreement','answer':f"Spearman rho {rank_summary.loc[(rank_summary.method_a=='bulk') & (rank_summary.method_b=='continuous'),'spearman'].iloc[0]:.3f}"},
    {'question':'Three-bin vs continuous gene rank agreement','answer':f"Spearman rho {rank_summary.loc[(rank_summary.method_a=='segment') & (rank_summary.method_b=='continuous'),'spearman'].iloc[0]:.3f}"},
    {'question':'Top-decile continuous genes outside bulk top decile','answer':str(int((strong['continuous'] & ~strong['bulk']).sum()))},
    {'question':'Top-decile continuous genes outside segment top decile','answer':str(int((strong['continuous'] & ~strong['segment']).sum()))},
    {'question':'Strict trajectory-only genes','answer':str(len(trajectory))},
    {'question':'Strict genes with mean component >= half of fixed mismatch','answer':str(int(trajectory.fixed_mean_component_fraction.ge(.5).sum()))},
    {'question':'Strict trajectory-only fraction with G2G similarity < 0.5','answer':f"{trajectory.alignment_similarity.lt(.5).mean():.2f}" if len(trajectory) else 'none'},
    {'question':'Strict trajectory-only fraction warp-rich and G2G similarity >= 0.8','answer':f"{(trajectory.alignment_similarity.ge(.8) & trajectory.g2g_warp_fraction.ge(.3)).mean():.2f}" if len(trajectory) else 'none'},
    {'question':'Concordant top-decile pathways','answer':str(int(path_table.discovery_category.eq('conventional_continuous').sum()))},
    {'question':'Concordant pathway examples','answer':'; '.join(path_table.loc[path_table.discovery_category.eq('conventional_continuous')].nsmallest(3,'continuous_rank').pathway)},
    {'question':'Strict continuous-only pathway candidates','answer':str(len(continuous_path))},
    {'question':'Leading continuous-only pathways','answer':'; '.join(continuous_path.nsmallest(3,'continuous_rank').pathway)},
    {'question':'Observed binning-loss diagnostics','answer':f"{int(missed_genes.opposing_local_mean_directions.sum())} opposing local means; {int(missed_genes.within_same_bin_peak_shift.sum())} within-bin peak shifts; {int(missed_genes.mismatch_near_bin_boundary.sum())} near boundaries among strict gene candidates."},
    {'question':'Biological follow-up themes','answer':'; '.join(selected.head(5).pathway)},
    {'question':'Technical consistency guide','answer':'Every strict trajectory-only gene has four-pair minimum divergence percentile >= 0.8; pathway member fractions are tabulated separately.'},
    {'question':'Primary limitation','answer':'One human donor; sections are technical/spatial checks, not donor replicates.'}
])
summary.to_csv(OUTPUT / 'method_comparison_summary.csv',index=False)
display(summary)
display(selected[['library','pathway','discovery_category','bulk_rank','segment_rank',
                  'continuous_rank','dominant_fixed_region','fixed_profile_entropy_null_z']])
print('Inspect selected_pathway_gene_details.csv and Figures 4/6 before assigning mechanism labels.')
