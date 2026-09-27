# %% [markdown]
# # 12 · Where do human and mouse PT pathways differ?
#
# **Start with gene trajectories. Discover pathways from gene statistics. Return to pseudospace to explain the differences.**
#
# The PT coordinate is already established in notebook 03. We read it exactly as saved; we do not rebuild DPT, reverse it, or align species anew. This notebook is a separate analysis and leaves notebooks 06 and 09 intact.
#
# > **Reading guide**
# >
# > Each section introduces the question, the calculation, and how to read its output. Imports and settings appear where first needed. Numerical routines live in `pseudospace/pathway_remodeling.py` so the notebook stays focused on the analysis.
#
# | Object | Question |
# |---|---|
# | $T_{\mathrm{level},g}$ | Is there an approximately constant human–mouse difference? |
# | $T_{\mathrm{spatial},g}$ | Does the difference change with PT position? **Primary discovery.** |
# | $T_{\mathrm{total},g}$ | Is there any trajectory difference? |
# | $Z_g(s)$ | Where is the difference, and which species is higher? |
#
# > **What this cohort can tell us**
# >
# > The healthy comparison contains two mouse specimens and two human cortex sections from **one donor**. `HUK1_MED1` is a specimen name, not evidence of medullary tissue. Model uncertainty is conditional on observed structures; pathway nulls compare genes, not independent donors. All enrichment and stability results are exploratory.
#
# **Inputs:** notebook 03's `cross_species_pt_dpt.h5ad` and `ortholog_map_used.csv`; local Reactome, Hallmark, and KEGG JSON libraries. Notebook 09 is needed only for the final G2G/Jeffreys sensitivity comparison. Run in `kidney-pseudospace`. Tables, caches, and figures go to `results/pt_pathway_remodeling/` (or your configured results root), never into Git.

# %% [markdown]
# ## 1 · Find the input and output folders
#
# Set `PSEUDOSPACE_DATA_ROOT` and `PSEUDOSPACE_RESULTS_ROOT` for external data, or use `--data-root` and `--results-root` when running the script mirror. Command-line values take precedence. No scientific parameters are hidden in this setup cell.

# %%
import argparse
import os
import sys
from pathlib import Path

start = Path(__file__).resolve().parent if '__file__' in globals() else Path.cwd()
project = next(p for p in (start, *start.parents) if (p / 'pseudospace').is_dir())
sys.path.insert(0, str(project))
parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--data-root', type=Path)
parser.add_argument('--results-root', type=Path)
args, _ = parser.parse_known_args()
data_root = Path(args.data_root or os.environ.get('PSEUDOSPACE_DATA_ROOT') or project / 'data').expanduser()
results_root = Path(args.results_root or os.environ.get('PSEUDOSPACE_RESULTS_ROOT') or project / 'results').expanduser()
upstream = results_root / 'human_vs_healthy_mouse'
output = results_root / 'pt_pathway_remodeling'
output.mkdir(parents=True, exist_ok=True)
input_path = upstream / 'cross_species_pt_dpt.h5ad'
map_path = upstream / 'ortholog_map_used.csv'
for path in (input_path, map_path):
    if not path.is_file():
        raise FileNotFoundError(f'Run notebook 03 first: missing {path}')
print('Results folder:', output)

# %% [markdown]
# ## 2 · Check the cohort and preserve the coordinate
#
# We require the saved PT labels, shared pseudospace, raw counts, normalized expression, and measured-in-both-species flag. An ortholog column filled with structural zeros is not evidence that a gene was assayed.
#
# > **Support rule, fixed before fitting**
# >
# > Use the intersection of the four specimens' 1st–99th percentile coordinate intervals. The models, gene universe, and bulk benchmark use those same structures. This avoids comparing a well-observed trajectory with an extrapolated tail. Positions remain in notebook 03's original units.

# %%
import anndata as ad
import numpy as np
import pandas as pd
from IPython.display import display

adata = ad.read_h5ad(input_path)
required = {'comparison_species', 'sample', 'broad_tubule_marker_call', 'shared_pseudospace'}
if not required.issubset(adata.obs) or not {'counts', 'lognorm'}.issubset(adata.layers):
    raise ValueError('Notebook 03 artifact is missing required labels, coordinate, or expression layers.')
if 'measured_in_both_inputs' not in adata.var or adata.var.measured_in_both_inputs.dtype != bool:
    raise ValueError('A boolean measured_in_both_inputs audit is required.')
adata = adata[adata.obs.broad_tubule_marker_call.astype(str).eq('PT')].copy()
expected = {'mouse': ['Ctrl1A2', 'Ctrl1A4'], 'human': ['HUK1_COR1', 'HUK1_MED1']}
actual = adata.obs.groupby('comparison_species', observed=True)['sample'].apply(
    lambda v: sorted(v.astype(str).unique())).to_dict()
if actual != expected or not adata.obs_names.is_unique or not adata.var_names.is_unique:
    raise ValueError(f'Unexpected cohort or duplicate identifiers: {actual}')
position = adata.obs.shared_pseudospace.to_numpy(float)
if not np.isfinite(position).all() or position.min() < 0 or position.max() > 1:
    raise ValueError('Expected the already oriented [0, 1] coordinate from notebook 03.')
support = adata.obs.groupby('sample', observed=True).shared_pseudospace.quantile([.01, .99]).unstack()
lo, hi = support[.01].max(), support[.99].min()
if hi <= lo:
    raise ValueError('Specimens have no common coordinate support.')
adata = adata[(position >= lo) & (position <= hi)].copy()
position = adata.obs.shared_pseudospace.to_numpy(float)
human = adata.obs.comparison_species.astype(str).eq('human').to_numpy()
specimen = adata.obs['sample'].astype(str).to_numpy()
# Quantiles interpolate between structures: place grid endpoints inside retained observations.
lo = max(position[specimen == name].min() for name in np.unique(specimen))
hi = min(position[specimen == name].max() for name in np.unique(specimen))
support.to_csv(output / 'specimen_coordinate_support.csv')
display(adata.obs.groupby(['comparison_species', 'sample'], observed=True).size().rename('retained_PT_structures'))
print(f'Common observed interval: {lo:.3f}–{hi:.3f}; the coordinate is unchanged.')

# %% [markdown]
# ## 3 · Freeze the measurable gene universe
#
# Keep measured-in-both orthologs detected in at least **2% of structures on an equal-specimen basis**. We do not select DE genes or standardize gene expression. Per-species detection remains visible: measured-but-zero expression in one species is allowed.
#
# Mean expression, detection, and coverage are measured before looking at pathway results. Coverage is the fraction of occupied coordinate bins in which a gene is detected, averaged across specimens. These three quantities will define the matched pathway nulls.

# %%
from scipy import sparse
from pseudospace.pathway_remodeling import gene_covariates

min_detection = .02
expression = sparse.csr_matrix(adata.layers['lognorm'], dtype=float)
if not np.isfinite(expression.data).all() or (expression.data < 0).any():
    raise ValueError('Expected finite, nonnegative log-normalized expression.')
orthologs = pd.read_csv(map_path)
if not {'mouse_symbol', 'human_symbol'}.issubset(orthologs):
    raise ValueError('Accepted ortholog map lacks mouse/human symbols.')
if orthologs[['mouse_symbol', 'human_symbol']].isna().any().any() or orthologs.mouse_symbol.duplicated().any() or orthologs.human_symbol.duplicated().any():
    raise ValueError('Expected the accepted one-to-one ortholog map.')
universe = gene_covariates(expression, position, human, specimen, adata.var_names)
universe['measured_in_both'] = adata.var.measured_in_both_inputs.to_numpy()
universe['accepted_ortholog'] = universe.gene.isin(orthologs.mouse_symbol)
universe['eligible'] = universe.measured_in_both & universe.accepted_ortholog & universe.detection.ge(min_detection)
universe.to_csv(output / 'gene_universe.csv', index=False)
genes = universe.loc[universe.eligible, 'gene'].tolist()
if len(genes) < 30:
    raise ValueError('Too few eligible genes for competitive pathway testing.')
Y = expression[:, universe.eligible.to_numpy()]
covariates = universe.set_index('gene').loc[genes]
print(f'{len(genes):,} / {adata.n_vars:,} genes retained; no DE filtering.')
display(covariates[['mean_expression', 'detection', 'detection_mouse', 'detection_human', 'coverage']].describe().round(3))

# %% [markdown]
# ## 4 · Fit the three nested gene models
#
# For log-normalized gene expression, use a **Gaussian regression-spline GAM** with a small, fixed cubic basis:
#
# $$M_0: f(s)+\text{specimen contrasts},$$
# $$M_{\rm level}: f(s)+\text{species}+\text{specimen contrasts},$$
# $$M_{\rm full}: f(s)+\text{species}+\text{species}:f(s)+\text{specimen contrasts}.$$
#
# The specimen intercept contrasts sum to zero *within each species*. This makes the species comparison identifiable. Each species receives half the total weight; its specimens share that weight equally. The fitted curves therefore describe the equal-specimen comparison, with section offsets adjusted.
#
# > **Why a fixed basis?**
# >
# > Six spline basis columns, plus the intercept, provide a modest amount of flexibility. Using the same unpenalized basis in all models preserves exact nesting. This is a regression-spline GAM, not an NB count model or a smoothing-parameter-selected GAMM. We test smaller/larger bases later.
#
# For each comparison we divide the weighted residual improvement per added coefficient by the full comparison model's residual variance. These partial-F statistics account for residual noise, but are used as **conditional ranking statistics**, not donor-level p-values. $T_{\rm spatial}$ includes changes in amplitude as well as position; it is not a pure relocation statistic.
#
# The full model gives $\delta_g(s)=\hat\mu_{H,g}(s)-\hat\mu_{M,g}(s)$. Its local standard error uses a **heteroskedasticity-robust HC3 sandwich**, because balancing weights are not inverse-variance weights. It still treats structures as independent and does not account for donor replication, spatial correlation, or uncertainty in the established coordinate. $Z=\delta/SE$ is a working signal-to-noise measure, not a calibrated local biological test.
#
# [Weighted regression assumptions](https://www.statsmodels.org/stable/examples/notebooks/generated/wls.html).

# %%
from pseudospace.pathway_remodeling import fit_nested_trajectories
from pseudospace.stage_cache import cached_payload, cached_frame, digest
import pseudospace.pathway_remodeling as remodeling
import pseudospace.levelshape as levelshape
import pseudospace.stats_gam as stats_gam

basis_df = 6
NOTEBOOK_LOGIC_VERSION = '12-pathway-remodeling-v2'  # Bump when a cached calculation changes.
grid = np.linspace(lo, hi, 61)
# Require local data from every specimen, not only overlapping range endpoints.
local_counts = pd.DataFrame({name: [np.sum((specimen == name) & (abs(position - s) <= .08 * (hi - lo)))
                                  for s in grid] for name in np.unique(specimen)}, index=grid)
local_counts.to_csv(output / 'grid_support_counts.csv', index_label='position')
if local_counts.min().min() < 15:
    raise ValueError('A grid neighborhood has fewer than 15 structures in a specimen; inspect support before fitting.')
cache_code = digest([NOTEBOOK_LOGIC_VERSION, Path(remodeling.__file__), Path(levelshape.__file__), Path(stats_gam.__file__)])
fit_inputs = {'Y': digest(Y), 'position': position, 'human': human, 'specimen': specimen, 'genes': genes}
fit = cached_payload('nested_gene_models',
    lambda: fit_nested_trajectories(Y, position, human, specimen, grid, basis_df=basis_df),
    root=output / 'stage_cache', params={'basis_df': basis_df, 'grid': grid},
    inputs=fit_inputs, code=cache_code)
gene_stats = covariates.copy()
for name in ('T_level', 'T_spatial', 'T_total'):
    gene_stats[name] = fit[name]
gene_stats['mean_human_minus_mouse'] = fit['delta'].mean(axis=1)
gene_stats['spatial_rms'] = fit['delta'].std(axis=1)
gene_stats.to_csv(output / 'gene_statistics.csv')
np.savez_compressed(output / 'gene_trajectories.npz', genes=np.array(genes), **fit)
display(gene_stats.sort_values('T_spatial', ascending=False).head(10))

# %% [markdown]
# ## 5 · Load broad pathway libraries and audit coverage
#
# Map library members through notebook 03's accepted ortholog table. Retain pathways with **10–300 eligible genes**, counting after the measurement and detection filters. Test all eligible terms; choose biology only after discovery.
#
# > **Read the denominator**
# >
# > `measured_fraction` is the fraction of original library members reaching the measured-in-both universe. `tested_fraction` includes the detection filter too. A pathway name with low coverage is only a partial view of that pathway.

# %%
import json
from pseudospace.pathways import build_pathway_membership

libraries = ('Reactome_2022', 'MSigDB_Hallmark_2020', 'KEGG_2019_Mouse')
library_dir = data_root / 'mouse_vs_human' / 'pathway_gene_sets'
coverage_parts = []
for library in libraries:
    path = library_dir / f'{library}.json'
    if not path.is_file():
        raise FileNotFoundError(f'Missing local pathway library: {path}')
    coverage_parts.append(build_pathway_membership(json.loads(path.read_text()),
        universe.loc[universe.measured_in_both & universe.accepted_ortholog, 'gene'],
        ortholog_map=orthologs, library_name=library, tested=genes, min_genes=1))
coverage = pd.concat(coverage_parts, ignore_index=True)
coverage['pathway_id'] = coverage.library + '::' + coverage.pathway
coverage['measured_fraction'] = coverage.n_assayed / coverage.n_requested.clip(lower=1)
coverage['tested_fraction'] = coverage.n_tested / coverage.n_requested.clip(lower=1)
coverage['tested_here'] = coverage.n_tested.between(10, min(300, len(genes) - 1))
coverage.to_csv(output / 'pathway_coverage.csv', index=False)
gene_sets = coverage.loc[coverage.tested_here].set_index('pathway_id').genes_present.to_dict()
if not gene_sets:
    raise ValueError('No pathway has 10–300 eligible members.')
print(f'{len(gene_sets):,} eligible pathways across {len(libraries)} libraries.')
display(coverage.groupby('library').agg(terms=('pathway', 'size'), tested=('tested_here', 'sum')))

# %% [markdown]
# ## 6 · Primary discovery: matched rank-AUC enrichment
#
# For each statistic, compare all member genes with all eligible nonmembers:
#
# $$AUC_P=P(T_{g\in P}>T_{g\notin P})+\tfrac12P(\text{tie}),\qquad \text{effect}=AUC_P-0.5.$$
#
# An AUC of 0.70 means a random pathway member outranks a random background gene 70% of the time, counting ties as half. These statistics are unsigned: enrichment means **more remodeling**, not upregulation.
#
# > **Matched empirical null**
# >
# > Divide mean expression, detection, and coverage into coarse tertiles (ties stay together). Each random set has the pathway's exact size and stratum counts, sampled without replacement from the full eligible universe. Test the upper tail using $(1+\#\{AUC_{null}\ge AUC_{observed}\})/(B+1)$.
#
# BH correction spans **all pathways, libraries, and three model statistics together**. Nominal Mann–Whitney p-values are shown only as diagnostics. The matching addresses measurability; it does **not** preserve inter-gene correlation or create biological replication.
#
# `null_auc_sd` and `fixed_member_fraction` expose overly restricted nulls. With 9,999 draws the smallest empirical p-value is 0.0001; `mc_se` reports Monte Carlo uncertainty. Borderline discoveries need more draws and a rerun of the entire family, not selective extra testing of attractive terms.

# %%
from pseudospace.pathway_remodeling import matching_strata, matched_pathway_tests

n_null = 9999
seed = 12
strata = matching_strata(covariates, bins=3)
match_audit = covariates[['mean_expression', 'detection', 'coverage']].assign(stratum=strata)
match_audit.to_csv(output / 'matching_strata.csv')
rankings = gene_stats[['T_total', 'T_level', 'T_spatial']]
pathway_tests = cached_frame('matched_primary',
    lambda: matched_pathway_tests(rankings, gene_sets, strata, n_null=n_null, seed=seed),
    root=output / 'stage_cache', params={'n_null': n_null, 'seed': seed},
    inputs={'statistics': rankings, 'sets': gene_sets, 'strata': strata}, code=cache_code)
pathway_tests.to_csv(output / 'pathway_rank_auc.csv', index=False)
spatial = pathway_tests[pathway_tests.statistic.eq('T_spatial')].copy()
spatial = spatial.merge(coverage[['pathway_id', 'library', 'pathway', 'tested_fraction']], on='pathway_id')
spatial = spatial.sort_values(['q_empirical', 'effect', 'pathway_id'], ascending=[True, False, True])
display(spatial[['library', 'pathway', 'n_genes', 'effect', 'q_empirical', 'mc_se', 'tested_fraction']].head(15))

# %% [markdown]
# ### Freeze the candidates before inspecting local curves
#
# A discovery candidate has positive spatial enrichment and empirical **q ≤ 0.05**. This cutoff is a reproducible exploratory shortlist, not a confirmatory claim. If none pass, the notebook still exports the full screen and plots up to six highest-ranked terms clearly marked as illustrations. It never converts a top-ranked term into a discovery automatically.

# %%
candidate_ids = spatial.loc[spatial.q_empirical.le(.05) & spatial.effect.gt(0), 'pathway_id'].tolist()
illustration_ids = candidate_ids[:6] if candidate_ids else spatial.head(6).pathway_id.tolist()
print(f'{len(candidate_ids)} exploratory candidates.')
print('Plot selection:', 'discovery candidates' if candidate_ids else 'ranked illustrations; none passed the cutoff')

# %% [markdown]
# ## 7 · Complementary preranked GSEA and leading edges
#
# Use the **same numerical rankings and same pathway members** as the AUC screen. GSEA asks whether a subset concentrates near the top of the ranking; it need not agree with a pathway-wide AUC shift.
#
# Use the classic unweighted running sum (`weight=0`): the ordering of these F statistics matters, while their magnitudes are not expression effects. Ties are ordered reproducibly by gene name; their frequency is reported because large tied blocks can make leading edges unstable. Positive NES means enrichment at the **high-remodeling** end, not human-high expression.
#
# > **Different null, different role**
# >
# > GSEA uses 999 gene-set permutations, not the covariate-matched null above (9,999 draws). Its default p-value resolution is 0.001; increase `gsea_permutations` for a finer complementary screen. Its BH column spans all pathway/statistic pairs in this GSEA run. Leading edges support interpretation; GSEA agreement is not required for discovery. The package's own FDR is also retained with its original name.
#
# [GSEA guide: preranked interpretation and ties](https://docs.gsea-msigdb.org/GSEA/GSEA_User_Guide/).

# %%
import gseapy as gp
from importlib.metadata import version
from statsmodels.stats.multitest import multipletests

gsea_permutations = 999  # Complementary screen; increase for finer GSEA p-value resolution.

def run_gsea(statistic_frame):
    tables = []
    for statistic in statistic_frame:
        ordered = statistic_frame[statistic].rename_axis('gene').reset_index().sort_values(
            [statistic, 'gene'], ascending=[False, True], kind='stable')
        result = gp.prerank(rnk=ordered, gene_sets=gene_sets, min_size=10, max_size=300,
            permutation_num=gsea_permutations, weight=0, seed=seed, threads=2, outdir=None, no_plot=True)
        table = result.res2d.rename(columns={'Term': 'pathway_id'}).copy()
        table['statistic'] = statistic
        tables.append(table)
    frame = pd.concat(tables, ignore_index=True)
    # Package permutation p-values can be zero; expose the finite permutation floor.
    frame['p_resolution_limited'] = pd.to_numeric(frame['NOM p-val']).clip(lower=1 / (gsea_permutations + 1))
    frame['q_bh_family'] = multipletests(frame.p_resolution_limited, method='fdr_bh')[1]
    return frame

print('Tied-score fractions:', {c: round(1 - rankings[c].nunique() / len(rankings), 3) for c in rankings})
gsea = cached_frame('preranked_gsea', lambda: run_gsea(rankings), root=output / 'stage_cache',
    params={'seed': seed, 'permutations': gsea_permutations, 'weight': 0, 'gseapy': version('gseapy')},
    inputs={'rankings': rankings, 'sets': gene_sets}, code=cache_code)
gsea.to_csv(output / 'pathway_gsea.csv', index=False)
leading_edges = {r.pathway_id: str(r.Lead_genes).split(';') for r in
    gsea.loc[
        gsea.statistic.eq('T_spatial') & pd.to_numeric(gsea.NES).gt(0)].itertuples()}
display(gsea[gsea.statistic.eq('T_spatial')].sort_values('q_bh_family').head(10))

# %% [markdown]
# ## 8 · Return to pseudospace: magnitude and direction
#
# At each saved grid position, rank **all eligible genes**, then compute:
#
# $$D_P(s)=AUC_P(|Z(s)|)-0.5,\qquad S_P(s)=2\{AUC_P(Z(s))-0.5\}.$$
#
# - **Divergence $D$:** positive values indicate greater local separation than background. Opposite gene directions cannot cancel inside $|Z|$.
# - **Relative direction $S$:** positive values mean more human-directed Z values than background; negative values mean more mouse-directed values.
#
# > **Do not confuse relative with absolute direction**
# >
# > A positive $S$ can occur even when all members are negative but background genes are more negative. We therefore also report median member Z and the fraction of members with $\delta>0$. Use these and the actual fitted gene curves before calling a pathway human-high or mouse-high.
#
# These curves summarize ranks of gene differences. We never average pathway expression or gene expression z-scores. No gridwise significance is claimed, and discovery-selected curves are descriptive follow-up.

# %%
from pseudospace.pathway_remodeling import local_pathway_curves, spatial_descriptions

# Include illustrations only when needed; they remain explicitly separate from discoveries.
local_ids = list(dict.fromkeys(candidate_ids + illustration_ids))
local_sets = {p: gene_sets[p] for p in local_ids}
local_curves = local_pathway_curves(fit['z'], genes, local_sets, grid)
phenotypes = spatial_descriptions(local_curves)
phenotypes['discovery_candidate'] = phenotypes.pathway_id.isin(candidate_ids)
local_curves.to_csv(output / 'pathway_local_curves.csv', index=False)
phenotypes.to_csv(output / 'pathway_spatial_descriptions.csv', index=False)
display(phenotypes.head(10))

# %% [markdown]
# ### Read the spatial descriptions as measurements
#
# Peak position is the maximum of D on the supported grid. The affected width is the total grid fraction at or above half the positive peak, converted to coordinate units; it can include separated regions. A broad affected fraction describes a diffuse effect. Early/mid/late shares divide the supported interval into thirds **after discovery**; these are positional descriptions, not anatomical segment labels. Sign changes in S ignore a ±0.1 deadband to avoid counting tiny wiggles.

# %%
import matplotlib.pyplot as plt

plt.rcParams.update({'axes.spines.top': False, 'axes.spines.right': False, 'font.size': 10})
fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True, layout='constrained')
for pathway in illustration_ids:
    trace = local_curves[local_curves.pathway_id.eq(pathway)]
    label = pathway.replace('::', ': ', 1)
    axes[0].plot(trace.position, trace.divergence, label=label)
    axes[1].plot(trace.position, trace.direction)
axes[0].set(ylabel='D: divergence rank-AUC − 0.5', ylim=(-.5, .5),
    title='Exploratory candidates' if candidate_ids else 'Ranked illustrations — no discoveries at q ≤ 0.05')
axes[1].set(xlabel='Established shared PT pseudospace', ylabel='S: relative direction', ylim=(-1, 1))
for ax in axes:
    ax.axhline(0, color='0.6', lw=.7)
axes[1].legend(*axes[0].get_legend_handles_labels(), fontsize=8,
               loc='upper left', bbox_to_anchor=(0, -.22))
fig.savefig(output / 'pathway_local_curves.pdf', bbox_inches='tight')
plt.show()
plt.close(fig)

# %% [markdown]
# ## 9 · Compare with conventional PT pseudobulk
#
# The clean comparison is already inside the nested models: level versus spatial, using identical expression data and enrichment machinery.
#
# For a recognizable independent *method* benchmark, sum **raw counts per specimen**, normalize each pseudobulk library over all measured-in-both accepted orthologs (before detection filtering), and use log2(CPM + 1). The ranking is the absolute Welch t statistic across the two mouse specimens and two human sections; the signed difference in mean log2(CPM + 1) is retained as the bulk effect.
#
# > **This is a descriptive pseudobulk benchmark**
# >
# > The human sections are not independent donors. We therefore do not report biological DE p-values or pretend that a count-based donor-level DE model can solve the missing replication. The t statistic is a noise-scaled ranking of this dataset, not confirmatory DE. Summing log-normalized expression would not be pseudobulk and is not done here.
#
# Bulk uses the same supported structures, eligible genes, pathway members, matching strata, and rank-AUC/GSEA machinery. Its empirical BH family spans the bulk pathway tests. Primary model-test q-values remain frozen.

# %%
from pseudospace.specimen import pseudobulk_profiles

counts = sparse.csr_matrix(adata.layers['counts'], dtype=float)
if not np.isfinite(counts.data).all() or (counts.data < 0).any() or not np.allclose(counts.data, np.round(counts.data)):
    raise ValueError('Pseudobulk requires nonnegative raw integer counts.')
bulk_all = pseudobulk_profiles(counts, specimen, position, n_bins=1,
    gene_names=adata.var_names, min_structures_per_bin=1).set_index('specimen')
measured_genes = universe.loc[universe.measured_in_both & universe.accepted_ortholog, 'gene']
library_sizes = bulk_all[['count_' + g for g in measured_genes]].sum(axis=1)
if library_sizes.le(0).any():
    raise ValueError('An empty pseudobulk library cannot be normalized.')
bulk_expression = np.log2(bulk_all[['count_' + g for g in genes]].div(library_sizes, axis=0) * 1e6 + 1)
bulk_expression.columns = genes
bulk_expression.to_csv(output / 'specimen_pseudobulk_logcpm.csv')
bulk_mouse = bulk_expression.loc[expected['mouse']]
bulk_human = bulk_expression.loc[expected['human']]
bulk_effect = bulk_human.mean() - bulk_mouse.mean()
bulk_se = np.sqrt(bulk_human.var(ddof=1) / 2 + bulk_mouse.var(ddof=1) / 2).clip(lower=1e-6)
bulk_stats = pd.DataFrame({'bulk_logcpm_effect': bulk_effect, 'T_bulk': (bulk_effect / bulk_se).abs()})
bulk_stats.to_csv(output / 'bulk_gene_statistics.csv')
bulk_tests = cached_frame('matched_bulk',
    lambda: matched_pathway_tests(bulk_stats[['T_bulk']], gene_sets, strata, n_null=n_null, seed=seed),
    root=output / 'stage_cache', params={'n_null': n_null, 'seed': seed},
    inputs={'statistics': bulk_stats, 'sets': gene_sets, 'strata': strata}, code=cache_code)
bulk_tests.to_csv(output / 'bulk_pathway_rank_auc.csv', index=False)
bulk_gsea = cached_frame('bulk_gsea', lambda: run_gsea(bulk_stats[['T_bulk']]),
    root=output / 'stage_cache', params={'seed': seed, 'permutations': gsea_permutations, 'weight': 0, 'gseapy': version('gseapy')},
    inputs={'statistics': bulk_stats, 'sets': gene_sets}, code=cache_code)
bulk_gsea.to_csv(output / 'bulk_pathway_gsea.csv', index=False)

# %% [markdown]
# ### What does continuous analysis add?
#
# | Level/bulk evidence | Spatial evidence | Reading |
# |---|---|---|
# | High | Low | Level or bulk largely captures the observed difference |
# | Low | High | Candidate continuous-only discovery |
# | High | High | Bulk detects a difference; continuous analysis adds its spatial structure |
# | Low | Low | Weak evidence under this screen |
#
# “High” below uses positive effect and matched q ≤ 0.05. “Low” means not passing this exploratory cutoff, **not proof of no effect**. A pathway found by both methods remains valuable: report bulk rank effect, spatial rank effect, peak location, local direction, sign changes, and genes driving the peak.

# %%
all_tests = pd.concat([pathway_tests, bulk_tests], ignore_index=True)
comparison = all_tests.pivot(index='pathway_id', columns='statistic', values=['effect', 'q_empirical'])
comparison.columns = ['_'.join(c) for c in comparison.columns]
comparison = comparison.reset_index()
level_high = comparison.q_empirical_T_level.le(.05) & comparison.effect_T_level.gt(0)
bulk_high = comparison.q_empirical_T_bulk.le(.05) & comparison.effect_T_bulk.gt(0)
spatial_high = comparison.pathway_id.isin(candidate_ids)
comparison['comparison_group'] = np.select(
    [spatial_high & ~(level_high | bulk_high), spatial_high & (level_high | bulk_high),
     ~spatial_high & (level_high | bulk_high)],
    ['continuous-only candidate', 'bulk/level plus spatial structure', 'bulk/level captures difference'],
    default='weak evidence')
comparison = comparison.merge(phenotypes, on='pathway_id', how='left')
comparison['median_member_bulk_effect'] = comparison.pathway_id.map(
    {p: bulk_stats.loc[members, 'bulk_logcpm_effect'].median() for p, members in gene_sets.items()})
comparison.to_csv(output / 'pathway_continuous_vs_bulk.csv', index=False)
display(comparison.groupby('comparison_group').size().rename('pathways'))
display(comparison[comparison.pathway_id.isin(candidate_ids)].head(12))

# %% [markdown]
# ## 10 · Inspect the genes behind each pathway
#
# These roles can overlap; they are explanations, not new discovery filters:
#
# - **Broad supporters:** members in the top quarter of the full-universe spatial-statistic ranking.
# - **Leading edge:** genes in the positive spatial GSEA leading edge.
# - **Spatial drivers:** up to five members with the largest |Z| near the pathway's divergence peak (the peak grid point and its two neighbors on either side).
#
# The member table includes detection, bulk effect, and local direction. Inspect the actual human and mouse fitted curves for drivers, especially weakly detected genes. No pathway gets credit for an isolated large but uninterpretable curve.

# %%
gene_index = pd.Index(genes)
member_parts = []
for pathway in local_ids:
    detail = gene_stats.loc[gene_sets[pathway]].copy()
    peak = phenotypes.set_index('pathway_id').loc[pathway, 'peak_position']
    peak_index = int(np.argmin(abs(grid - peak)))
    window = slice(max(0, peak_index - 2), min(len(grid), peak_index + 3))
    idx = gene_index.get_indexer(detail.index)
    detail['pathway_id'] = pathway
    detail['broad_supporter'] = detail.T_spatial.ge(gene_stats.T_spatial.quantile(.75))
    detail['leading_edge'] = detail.index.isin(leading_edges.get(pathway, []))
    detail['peak_abs_z'] = np.abs(fit['z'][idx, window]).mean(axis=1)
    detail['peak_delta'] = fit['delta'][idx, window].mean(axis=1)
    detail['spatial_driver'] = detail.peak_abs_z.rank(ascending=False, method='first').le(5)
    detail['bulk_logcpm_effect'] = bulk_stats.loc[detail.index, 'bulk_logcpm_effect']
    member_parts.append(detail.reset_index())
member_evidence = pd.concat(member_parts, ignore_index=True)
member_evidence.to_csv(output / 'pathway_member_evidence.csv', index=False)
display(member_evidence[member_evidence.spatial_driver].head(15))

# %%
# Plot two drivers from each of the first three illustrated pathways.
for pathway in illustration_ids[:3]:
    drivers = member_evidence[member_evidence.pathway_id.eq(pathway)].nlargest(2, 'peak_abs_z')
    fig, axes = plt.subplots(1, 2, figsize=(10, 3), sharex=True, layout='constrained')
    for ax, row in zip(axes, drivers.itertuples()):
        i = gene_index.get_loc(row.gene)
        ax.plot(grid, fit['mouse'][i], color='#0072B2', label='mouse')
        ax.plot(grid, fit['human'][i], color='#D55E00', label='human')
        ax.set(title=f'{row.gene} · detection {row.detection:.0%}',
               xlabel='Shared PT pseudospace', ylabel='Fitted log-normalized expression')
        ax.legend(frameon=False)
    fig.suptitle(pathway.replace('::', ': ', 1))
    fig.savefig(output / f'drivers_{illustration_ids.index(pathway) + 1:02d}.pdf', bbox_inches='tight')
    plt.show()
    plt.close(fig)

# %% [markdown]
# ## 11 · Collapse redundant terms after discovery
#
# Cluster only the frozen discovery candidates using equal-weight similarity from member overlap, leading-edge overlap, D-curve correlation, and S-curve correlation. Complete linkage avoids merging a whole pathway hierarchy through a chain of weakly connected terms. The distance cutoff is a presentation choice, not a test; missing leading edges contribute no similarity.
#
# The representative is the member with the strongest primary spatial evidence. Every original term and its statistics remain in the exported tables. The data determine how many programs remain; we do not force 10–20 groups or assign biological program names automatically.

# %%
from pseudospace.pathway_remodeling import collapse_programs

programs = collapse_programs(candidate_ids, gene_sets, leading_edges, local_curves, distance_cut=.45)
programs = programs.merge(spatial[['pathway_id', 'effect', 'q_empirical']], on='pathway_id', how='left')
representatives = programs.sort_values(['q_empirical', 'effect'], ascending=[True, False]).drop_duplicates('program')
programs['representative'] = programs.pathway_id.isin(representatives.pathway_id)
programs.to_csv(output / 'pathway_program_membership.csv', index=False)
display(representatives.head(20))

# %% [markdown]
# ## 12 · Ask whether the spatial evidence is fragile
#
# Refit all genes, then retest the full eligible pathway family under:
#
# 1. Leaving out each mouse and each human section in turn.
# 2. Spline basis sizes 4 and 8 instead of 6.
# 3. Raising the detection threshold from 2% to 5%.
#
# The coordinate and grid stay fixed. Each refit rebalances the remaining specimens; the detection sensitivity updates membership and the matched gene background. Every sensitivity run corrects across **all three statistics and all eligible pathways**, so candidate-only BH cannot make retention easier.
#
# > **Stability of this dataset**
# >
# > A pathway is retained if its spatial effect stays positive and its matched q stays ≤ 0.05. Report retained / testable runs **and** testable / planned runs; unavailable pathways are not silently counted as failures or successes. Leave-section-out is a section sensitivity, not human donor resampling. Each run uses the same number of matched draws as discovery so the Monte Carlo floor does not change.
#
# Model fits and pathway tests are cached separately. The existing cache switch `PSEUDOSPACE_STAGE_CACHE=0` forces a rebuild. This is the slowest section on a first run. If no candidates pass discovery, this section records an empty stability table and skips refitting.

# %%
variants = [(f'omit_{name}', specimen != name, basis_df, min_detection) for name in np.unique(specimen)]
variants += [('basis_4', np.ones(len(position), bool), 4, min_detection),
             ('basis_8', np.ones(len(position), bool), 8, min_detection),
             ('detection_5pct', np.ones(len(position), bool), basis_df, .05)]
stability_parts = []
for label, keep_rows, df, detection_cut in (variants if candidate_ids else []):
    keep_genes = covariates.detection.ge(detection_cut).to_numpy()
    variant_genes = gene_index[keep_genes]
    variant_sets = {p: sorted(set(members) & set(variant_genes)) for p, members in gene_sets.items()}
    variant_sets = {p: members for p, members in variant_sets.items() if 10 <= len(members) <= min(300, len(variant_genes) - 1)}
    if not variant_sets:
        stability_parts.append(pd.DataFrame({'pathway_id': candidate_ids, 'variant': label,
            'testable': False, 'retained': False, 'effect': np.nan, 'q_empirical': np.nan}))
        continue
    variant_Y = Y[keep_rows][:, keep_genes]
    variant_inputs = {'Y': digest(variant_Y), 'position': position[keep_rows], 'human': human[keep_rows],
                      'specimen': specimen[keep_rows], 'genes': variant_genes.tolist()}
    variant_fit = cached_payload(label + '_fit',
        lambda: fit_nested_trajectories(variant_Y, position[keep_rows], human[keep_rows], specimen[keep_rows], grid, basis_df=df),
        root=output / 'stage_cache', params={'basis_df': df, 'grid': grid}, inputs=variant_inputs, code=cache_code)
    variant_scores = pd.DataFrame({name: variant_fit[name] for name in rankings}, index=variant_genes)
    variant_covariates = gene_covariates(variant_Y, position[keep_rows], human[keep_rows], specimen[keep_rows], variant_genes)
    variant_strata = matching_strata(variant_covariates, bins=3)
    variant_tests = cached_frame(label + '_pathways',
        lambda: matched_pathway_tests(variant_scores, variant_sets, variant_strata, n_null=n_null, seed=seed),
        root=output / 'stage_cache', params={'n_null': n_null, 'seed': seed},
        inputs={'statistics': variant_scores, 'sets': variant_sets, 'strata': variant_strata}, code=cache_code)
    variant_spatial = variant_tests[variant_tests.statistic.eq('T_spatial')].set_index('pathway_id')
    candidate_run = variant_spatial.reindex(candidate_ids).reset_index()
    candidate_run['variant'] = label
    candidate_run['testable'] = candidate_run.q_empirical.notna()
    candidate_run['retained'] = candidate_run.effect.gt(0) & candidate_run.q_empirical.le(.05)
    stability_parts.append(candidate_run)
    print(label, ':', int(candidate_run.retained.sum()), 'of', int(candidate_run.testable.sum()), 'testable candidates retained')
stability = (pd.concat(stability_parts, ignore_index=True) if stability_parts else
             pd.DataFrame(columns=['pathway_id', 'variant', 'testable', 'retained']))
stability.to_csv(output / 'pathway_stability_runs.csv', index=False)
stability_summary = stability.groupby('pathway_id').agg(n_testable=('testable', 'sum'), n_retained=('retained', 'sum'))
stability_summary['n_planned'] = len(variants)
stability_summary['retention_fraction'] = stability_summary.n_retained / stability_summary.n_testable.replace(0, np.nan)
stability_summary.to_csv(output / 'pathway_stability_summary.csv')
display(stability_summary.head(15))

# %% [markdown]
# ## 13 · Optional G2G / Jeffreys sensitivity from notebook 09
#
# Read the saved gene statistics only after freezing discovery. This comparison uses the intersection of eligible genes and reruns matched rank-AUC on the same intersected membership for all methods. High `1 − alignment_similarity` means poor G2G alignment; high fixed Jeffreys mismatch describes distributional differences without warping.
#
# These methods measure different things: GAM spatial statistics isolate nonconstant mean differences, G2G permits alignment, and Jeffreys mismatch includes local variance differences. Agreement is useful sensitivity evidence; disagreement can be informative. This is **not** independent validation or a replacement coordinate. A missing notebook-09 artifact is recorded as unavailable and does not stop primary analysis.

# %%
sensitivity_path = results_root / 'pt_genes2genes' / 'g2g_vs_fixed_gene_comparison.csv'
g2g_status = 'unavailable: run notebook 09 to add this sensitivity comparison'
g2g_tests = pd.DataFrame()
if sensitivity_path.is_file():
    saved = pd.read_csv(sensitivity_path).set_index('gene')
    columns = ['alignment_similarity', 'fixed_overall_mismatch']
    if not saved.index.is_unique or not set(columns).issubset(saved):
        raise ValueError('Notebook 09 sensitivity table has an unexpected schema.')
    shared = gene_index.intersection(saved.dropna(subset=columns).index, sort=False)
    if len(shared) < 30:
        raise ValueError('Too few genes shared with notebook 09.')
    sensitivity_scores = gene_stats.loc[shared, ['T_spatial']].assign(
        G2G_divergence=1 - saved.loc[shared, 'alignment_similarity'],
        Jeffreys_mismatch=saved.loc[shared, 'fixed_overall_mismatch'])
    sensitivity_sets = {p: sorted(set(members) & set(shared)) for p, members in gene_sets.items()}
    sensitivity_sets = {p: members for p, members in sensitivity_sets.items() if 10 <= len(members) <= min(300, len(shared) - 1)}
    sensitivity_strata = matching_strata(covariates.loc[shared], bins=3)
    g2g_tests = cached_frame('g2g_jeffreys_sensitivity',
        lambda: matched_pathway_tests(sensitivity_scores, sensitivity_sets, sensitivity_strata, n_null=n_null, seed=seed),
        root=output / 'stage_cache', params={'n_null': n_null, 'seed': seed},
        inputs={'statistics': sensitivity_scores, 'sets': sensitivity_sets, 'strata': sensitivity_strata}, code=cache_code)
    g2g_tests.to_csv(output / 'g2g_jeffreys_pathway_sensitivity.csv', index=False)
    concordance = sensitivity_scores.corr(method='spearman')
    concordance.to_csv(output / 'g2g_jeffreys_gene_rank_concordance.csv')
    display(concordance)
    g2g_status = f'completed on {len(shared)} shared genes; coordinate and support differ in notebook 09'
print(g2g_status)

# %% [markdown]
# ## 14 · Final evidence table and reproducibility record
#
# The final table retains **all original pathway statistics** and adds program membership, bulk comparison, spatial descriptions, and stability. It is ready for biological review, not a list of validated mechanisms.
#
# > **Before writing a biological conclusion**
# >
# > Check coverage and null resolution; inspect absolute gene directions and fitted curves; consider specimen and basis sensitivity; keep correlation and one-human-donor limitations beside the result. Neither a small pathway q-value nor a stable section-level pattern is population-level human inference.
#
# The table preserves both continuous-only candidates and pathways for which continuous analysis adds detail to a bulk-accessible difference. Program grouping and driver selection happen after discovery.

# %%
final_table = comparison.merge(programs[['pathway_id', 'program', 'representative']], on='pathway_id', how='left')
final_table = final_table.merge(stability_summary, on='pathway_id', how='left')
final_table = final_table.merge(coverage[['pathway_id', 'library', 'pathway', 'n_requested', 'n_assayed',
    'n_tested', 'measured_fraction', 'tested_fraction']], on='pathway_id', how='left')
for role in ('broad_supporter', 'leading_edge', 'spatial_driver'):
    role_genes = member_evidence[member_evidence[role]].groupby('pathway_id').gene.agg(';'.join)
    final_table[role + '_genes'] = final_table.pathway_id.map(role_genes).fillna('')
final_table['caveat'] = 'Exploratory gene-set evidence; two mouse specimens, two cortex sections from one human donor.'
final_table.to_csv(output / 'pathway_evidence_atlas.csv', index=False)
manifest = {'logic_version': NOTEBOOK_LOGIC_VERSION, 'implementation_fingerprint': cache_code,
    'fit_input_fingerprint': digest(fit_inputs), 'input_fingerprint': digest(input_path),
    'ortholog_fingerprint': digest(map_path), 'expression_fingerprint': fit_inputs['Y'],
    'coordinate': '03 shared_pseudospace; unchanged', 'common_support': [float(lo), float(hi)],
    'n_genes': len(genes), 'n_pathways': len(gene_sets), 'basis_df': basis_df,
    'detection_threshold': min_detection, 'matching_bins': 3, 'n_null': n_null, 'gsea_permutations': gsea_permutations, 'seed': seed,
    'libraries': {name: digest(library_dir / f'{name}.json') for name in libraries},
    'packages': {name: version(name) for name in ['numpy', 'scipy', 'pandas', 'anndata', 'patsy', 'statsmodels', 'gseapy']},
    'g2g_status': g2g_status, 'n_candidates': len(candidate_ids),
    'caveat': final_table.caveat.iloc[0]}
(output / 'run_manifest.json').write_text(json.dumps(manifest, indent=2))
display(final_table[final_table.pathway_id.isin(candidate_ids)].sort_values('q_empirical_T_spatial').head(20))
print('Saved evidence atlas and run manifest to', output)
