# %% [markdown]
# # 12 · Where do human and mouse PT pathways differ?
#
# **Start with gene trajectories. Discover pathways from gene statistics. Return to pseudospace to explain the differences.**
#
# The primary PT coordinate comes from notebook 13’s nonbranching scFates curve. We read it exactly as saved; we do not refit the curve or align species anew. Notebook 13 compares it with DPT on the same PT structures and embedding. Its refitted PT labels can differ from notebook 03, so changes from earlier DPT-based results are not coordinate-only effects. This notebook is a separate pathway analysis and leaves notebooks 06 and 09 intact.
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
# **Inputs:** notebook 13’s `cross_species_pt_scfates.h5ad` (structure identities and coordinates only), `ortholog_map_used.csv`, and the four original `tubule_by_gene/*_tubule_by_gene_caleb.h5ad` count matrices; local Reactome, Hallmark, and KEGG JSON libraries. Run in `kidney-pseudospace`. Tables, caches, and figures go to `results/pt_pathway_remodeling_scfates/` (or your configured results root), never into Git.
#

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
upstream = results_root / 'minimal_pt_scfates'
output = results_root / 'pt_pathway_remodeling_scfates'
output.mkdir(parents=True, exist_ok=True)
input_path = upstream / 'cross_species_pt_scfates.h5ad'
map_path = upstream / 'ortholog_map_used.csv'
for path in (input_path, map_path):
    if not path.is_file():
        raise FileNotFoundError(f'Run notebook 13 first: missing {path}')
print('Results folder:', output)

# %% [markdown]
# ## 2 · Check the cohort and preserve the coordinate
#
# We read only the saved PT labels, specimen identities, and shared pseudospace from notebook 13. Its artifact contains metadata only, so section 3 reloads expression from the original specimens. Structure selection and coordinate orientation remain unchanged.
#
# > **Support rule, fixed before fitting**
# >
# > Use the intersection of the four specimens' 1st–99th percentile coordinate intervals. The models, gene universe, and bulk benchmark use those same structures. This avoids comparing a well-observed trajectory with an extrapolated tail. Positions remain on notebook 13’s saved [0, 1] scFates scale.

# %%
import anndata as ad
import numpy as np
import pandas as pd
from IPython.display import display

# Read metadata only: notebook 13 owns the coordinate, not this analysis gene universe.
saved_pt = ad.read_h5ad(input_path, backed='r')
pt_obs = saved_pt.obs.copy()
saved_pt.file.close()
adata = ad.AnnData(obs=pt_obs)
required = {'comparison_species', 'sample', 'broad_tubule_marker_call', 'shared_pseudospace'}
if not required.issubset(adata.obs):
    raise ValueError('Notebook 13 artifact is missing required labels or coordinate.')
adata = adata[adata.obs.broad_tubule_marker_call.astype(str).eq('PT')].copy()
expected = {'mouse': ['Ctrl1A2', 'Ctrl1A4'], 'human': ['HUK1_COR1', 'HUK1_MED1']}
actual = adata.obs.groupby('comparison_species', observed=True)['sample'].apply(
    lambda v: sorted(v.astype(str).unique())).to_dict()
if actual != expected or not adata.obs_names.is_unique or not adata.var_names.is_unique:
    raise ValueError(f'Unexpected cohort or duplicate identifiers: {actual}')
position = adata.obs.shared_pseudospace.to_numpy(float)
if not np.isfinite(position).all() or position.min() < 0 or position.max() > 1:
    raise ValueError('Expected the already oriented [0, 1] coordinate from notebook 13.')
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
# ## 3 · Start from all accepted orthologs and apply two rules
#
# We rebuild expression for the **same retained PT structures**, starting from every pair in the accepted one-to-one ortholog map. We do not inherit notebook 13’s 5% detection or 20-total-count embedding filters: those were used to construct its embedding, whose saved coordinate we continue to use unchanged.
#
# | Step | Rule | Purpose |
# |---|---|---|
# | Starting panel | All accepted one-to-one human–mouse ortholog pairs | Define comparable gene identities |
# | 1. Measurement availability | The gene must be present in every specimen's original feature list | Exclude missing measurements masquerading as zeros |
# | 2. PT detection | Mean of the four specimen-specific detection fractions ≥2% | Remove genes with very little observed PT expression |
#
# > **Detection means expression greater than zero**
# >
# > Each specimen contributes equally, regardless of its number of structures. The 2% threshold is not required separately in each species or each specimen. A gene measured but unexpressed in one species can remain eligible. A feature missing from an input cannot.
#
# ### Normalize before applying the detection threshold
#
# For each structure, sum raw counts over **all orthologs measured in every specimen**, scale to 10,000 counts, and apply natural-log `log1p`. That measured panel defines the denominator before the 2% filter. Changing the detection cutoff therefore cannot change the normalization of genes that remain. Pseudobulk uses raw counts and the same measured panel for its library-size denominator.
#
# > **A new analysis universe requires a new run**
# >
# > Restoring genes and changing the normalization denominator can change fitted curves, ranks, candidate pathways, and program groups. The old 9,914-gene results are not assumed to persist. The coordinate and retained structure identities stay fixed. Existing cache keys are invalidated by the new input and logic version.
#
# Mean expression and positional coverage are **matching covariates, not extra eligibility filters**. Coverage is the fraction of occupied coordinate bins with detection, averaged across specimens. All accepted pairs appear in the audit, including unmeasured genes that cannot enter testing.
#

# %%
from scipy import sparse
from pseudospace.pathway_inputs import rebuild_pt_expression
from pseudospace.pathway_remodeling import gene_covariates

orthologs = pd.read_csv(map_path)
source_paths = {name: data_root / 'tubule_by_gene' / f'{name}_tubule_by_gene_caleb.h5ad'
                for names in expected.values() for name in names}
# Recover counts for exactly these structures; normalize on the full measured ortholog panel.
pt_obs = adata.obs.copy()
adata = rebuild_pt_expression(pt_obs, source_paths, orthologs, target_sum=1e4)
pd.testing.assert_frame_equal(adata.obs[pt_obs.columns], pt_obs)
adata.var.to_csv(output / 'ortholog_measurement_audit.csv')
adata.obs[['sample', 'shared_pseudospace', 'ortholog_library_size']].to_csv(
    output / 'structure_expression_audit.csv')

min_detection = .02
expression = sparse.csr_matrix(adata.layers['lognorm'], dtype=float)
universe = gene_covariates(expression, position, human, specimen, adata.var_names)
universe['measured_in_both'] = adata.var.measured_in_both_inputs.to_numpy()
universe['accepted_ortholog'] = True  # The starting columns are exactly the accepted pairs.
universe['eligible'] = universe.measured_in_both & universe.detection.ge(min_detection)
universe.to_csv(output / 'gene_universe.csv', index=False)
genes = universe.loc[universe.eligible, 'gene'].tolist()
if len(genes) < 30:
    raise ValueError('Too few eligible genes for competitive pathway testing.')
Y = expression[:, universe.eligible.to_numpy()]
covariates = universe.set_index('gene').loc[genes]

# A three-row flow table is the complete gene-selection story for the paper.
gene_filter_flow = pd.DataFrame({'stage': ['Accepted one-to-one orthologs',
    'Measured in every specimen', 'Equal-specimen PT detection ≥2%'],
    'n_genes': [len(universe), int(universe.measured_in_both.sum()), len(genes)]})
gene_filter_flow['removed_at_step'] = gene_filter_flow.n_genes.shift(fill_value=len(universe)) - gene_filter_flow.n_genes
gene_filter_flow.to_csv(output / 'gene_filter_flow.csv', index=False)
display(gene_filter_flow)
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
#
# The full-model residuals are saved with the trajectory export for the correlation check below. They are **raw observed minus fitted expression**, including the fitted specimen effects. `residuals` has structures as rows and genes as columns; the saved `structures`, `specimen`, and `genes` arrays identify those axes. Only the primary fit retains this large matrix. After updating, restart the kernel and run all cells so the updated helper and residuals are loaded.

# %%
from pseudospace.pathway_remodeling import fit_nested_trajectories
from pseudospace.stage_cache import cached_payload, cached_frame, digest
import pseudospace.pathway_remodeling as remodeling
import pseudospace.pathway_inputs as pathway_inputs
import pseudospace.cross_species as cross_species
import pseudospace.levelshape as levelshape
import pseudospace.stats_gam as stats_gam

basis_df = 6
NOTEBOOK_LOGIC_VERSION = '12-pathway-remodeling-v6-scfates-coordinate'  # Bump when a cached calculation changes.
grid = np.linspace(lo, hi, 61)
# Require local data from every specimen, not only overlapping range endpoints.
local_counts = pd.DataFrame({name: [np.sum((specimen == name) & (abs(position - s) <= .08 * (hi - lo)))
                                  for s in grid] for name in np.unique(specimen)}, index=grid)
local_counts.to_csv(output / 'grid_support_counts.csv', index_label='position')
if local_counts.min().min() < 15:
    raise ValueError('A grid neighborhood has fewer than 15 structures in a specimen; inspect support before fitting.')
cache_code = digest([NOTEBOOK_LOGIC_VERSION, Path(remodeling.__file__), Path(levelshape.__file__), Path(stats_gam.__file__),
                     Path(pathway_inputs.__file__), Path(cross_species.__file__)])
fit_inputs = {'Y': digest(Y), 'position': position, 'human': human, 'specimen': specimen, 'genes': genes}
fit = cached_payload('nested_gene_models',
    lambda: fit_nested_trajectories(Y, position, human, specimen, grid, basis_df=basis_df, return_residuals=True),
    root=output / 'stage_cache', params={'basis_df': basis_df, 'grid': grid, 'return_residuals': True},
    inputs=fit_inputs, code=cache_code)
gene_stats = covariates.copy()
for name in ('T_level', 'T_spatial', 'T_total'):
    gene_stats[name] = fit[name]
gene_stats['mean_human_minus_mouse'] = fit['delta'].mean(axis=1)
gene_stats['spatial_rms'] = fit['delta'].std(axis=1)
gene_stats.to_csv(output / 'gene_statistics.csv')
np.savez_compressed(output / 'gene_trajectories.npz', genes=np.array(genes),
    structures=adata.obs_names.to_numpy(dtype=str), specimen=np.asarray(specimen, str), **fit)
display(gene_stats.sort_values('T_spatial', ascending=False).head(10))

# %% [markdown]
# ### Visual check: what each nested model permits
#
# The example gene is chosen by the largest observed $T_{\rm spatial}$ **only to illustrate the fit**. Each panel shows the same retained structures and a one-gene refit with the primary model’s knots and balancing weights. The plotted full-model curves are checked against the saved primary fit. Gray is the shared curve in $M_0$; $M_{\rm level}$ allows a constant species offset; $M_{\rm full}$ also allows that difference to vary with position. Points are individual log-normalized structures, colored by species. They are not independent donor replicates. The curves cover the common observed interval and show the equal-specimen comparison with specimen intercepts set to their within-species mean. A constrained Gaussian model can predict values below zero; these are fitted values on the log-normalized scale, not observed counts. Compare the three constraints rather than treating the chosen gene as independent evidence.

# %%
import matplotlib.pyplot as plt
from pseudospace.levelshape import build_ls_designs

example_gene = int(np.argmax(fit['T_spatial']))
example_name = genes[example_gene]
example_values = Y[:, example_gene].toarray().ravel()
example_weights = np.zeros(len(position))
example_nuisance = []
for group in (False, True):
    names = np.unique(specimen[human == group])
    for name in names:
        mask = specimen == name
        example_weights[mask] = len(position) / (2 * len(names) * mask.sum())
    example_nuisance.extend((specimen == name).astype(float) - (specimen == names[-1]).astype(float)
                            for name in names[:-1])
model_designs = [np.column_stack([design, *example_nuisance])
                 for design in build_ls_designs(position, human.astype(float), fit['knots'])[:3]]
root_weights = np.sqrt(example_weights)
model_betas = [np.linalg.lstsq(root_weights[:, None] * design,
                              root_weights * example_values, rcond=None)[0]
               for design in model_designs]
example_curves = {}
for label, group in [('mouse', 0.), ('human', 1.)]:
    grid_designs = build_ls_designs(grid, np.full(len(grid), group), fit['knots'])[:3]
    example_curves[label] = [design @ beta[:design.shape[1]]
                             for design, beta in zip(grid_designs, model_betas)]
np.testing.assert_allclose(example_curves['mouse'][2], fit['mouse'][example_gene])
np.testing.assert_allclose(example_curves['human'][2], fit['human'][example_gene])

colors = {'mouse': '#0072B2', 'human': '#D55E00'}
fig, axes = plt.subplots(1, 3, figsize=(12, 3.6), sharex=True, sharey=True, layout='constrained')
for i, (ax, title) in enumerate(zip(axes, ['M0: shared curve', 'Mlevel: parallel curves',
                                            'Mfull: position-dependent difference'])):
    for label, group in [('mouse', 0), ('human', 1)]:
        mask = human == group
        ax.scatter(position[mask], example_values[mask], s=2, alpha=.07,
                   color=colors[label], rasterized=True)
    if i == 0:
        ax.plot(grid, example_curves['mouse'][i], color='0.15', lw=2, label='Shared')
    else:
        for label in ('mouse', 'human'):
            ax.plot(grid, example_curves[label][i], color=colors[label], lw=2,
                    label=label.capitalize())
    ax.set(title=title, xlabel='Established shared PT pseudospace', xlim=(position.min(), position.max()))
    ax.legend(frameon=False, fontsize=8)
axes[0].set_ylabel('log1p(counts per 10,000)')
fig.suptitle(f'Nested models for {example_name} (illustration selected by spatial score)')
fig.savefig(output / 'nested_gam_example.pdf', bbox_inches='tight')
plt.show()
plt.close(fig)

# %% [markdown]
# ## 5 · Load broad pathway libraries and audit coverage
#
# Map library members through notebook 13’s accepted ortholog table. Retain pathways with **10–300 eligible genes**, counting after the measurement and detection filters. Test all eligible terms; choose biology only after discovery.
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
# ## 5a · Audit whether the matched null matters before discovery
#
# A gene's spatial score may depend on how readily that gene is measured. First, plot $T_{\rm spatial}$ against mean log-normalized expression, detection fraction, and occupied-bin coverage across **all eligible genes**. The vertical axis uses $\log(1+T_{\rm spatial})$ only to display its long tail; the Spearman correlations use the original score. These relationships indicate possible measurement bias, but do not by themselves establish whether pathway results change.
#
# Then test the **same complete pathway family and all three gene statistics** with two random-set nulls. The uniform null samples any eligible genes while preserving pathway size. The matched null also preserves each pathway's counts in the expression/detection/coverage strata. Both use 9,999 draws and the same seed; for each statistic, BH correction spans every tested pathway across all three libraries. Their observed AUCs must agree exactly: only the null distribution and resulting p/q values may differ. A large change in the candidate list shows that the sampling rule matters in this cohort; it does not prove that either null is biologically calibrated, especially for terms near the q cutoff or the Monte Carlo resolution limit. The comparison is diagnostic; the predeclared matched discovery rule in section 6 remains fixed. Neither null accounts for donor replication or within-pathway gene correlation.

# %%
import matplotlib.pyplot as plt

covariate_names = {'mean_expression': 'Mean log-normalized expression',
                   'detection': 'Detection fraction',
                   'coverage': 'Pseudospace coverage'}
covariate_rho = gene_stats[list(covariate_names)].corrwith(gene_stats.T_spatial, method='spearman')
covariate_rho.rename('spearman_with_T_spatial').to_csv(output / 'spatial_statistic_covariate_correlations.csv')
fig, axes = plt.subplots(1, 3, figsize=(12, 3.5), sharey=True, layout='constrained')
hexes = []
for ax, (column, label) in zip(axes, covariate_names.items()):
    hexes.append(ax.hexbin(gene_stats[column], np.log1p(gene_stats.T_spatial),
                           gridsize=40, mincnt=1, bins='log', cmap='Blues'))
    ax.set(xlabel=label, title=f'Spearman ρ = {covariate_rho[column]:.2f}')
common_max = max(hb.get_array().max() for hb in hexes)
for hb in hexes:
    hb.set_clim(1, common_max)
axes[0].set_ylabel('log1p(spatial gene statistic)')
fig.colorbar(hexes[-1], ax=axes, label='Genes per hexagon (log color scale)')
fig.savefig(output / 'spatial_statistic_measurement_covariates.pdf', bbox_inches='tight')
plt.show()
plt.close(fig)

# %% [markdown]
# ### Read the actual matching buckets
#
# The null uses quantile cuts with ties kept together, so a covariate can have fewer than three buckets. Each matrix fixes one **coverage** bucket; rows are mean-expression buckets and columns are detection buckets. A cell gives the number of eligible genes in that joint stratum, including zero-count cells. The interval labels show the cut rules ($<$ below a cut; $\ge$ at or above it), and the exported boundary table retains the unrounded cut values. The plot's assignments are checked against the strata used by the pathway test. Small or empty joint buckets restrict the replacement genes available to a matched draw; inspect `fixed_member_fraction` in the later pathway audit.

# %%
from pseudospace.pathway_remodeling import matching_strata

matching_columns = ('mean_expression', 'detection', 'coverage')
actual_strata = matching_strata(covariates, bins=3)
bucket_codes, bucket_labels, boundary_rows = {}, {}, []
for name in matching_columns:
    values = covariates[name].to_numpy(float)
    edges = np.unique(np.quantile(values, np.linspace(0, 1, 4)[1:-1]))
    codes = np.searchsorted(edges, values, side='right')
    labels = (['all values'] if not len(edges) else
              [f'< {edges[0]:.4g}'] +
              [f'{low:.4g}–< {high:.4g}' for low, high in zip(edges[:-1], edges[1:])] +
              [f'≥ {edges[-1]:.4g}'])
    bucket_codes[name], bucket_labels[name] = codes, labels
    for bucket, label in enumerate(labels):
        boundary_rows.append({'covariate': name, 'bucket': bucket + 1, 'interval': label,
                              'lower_cut': edges[bucket - 1] if bucket else np.nan,
                              'upper_cut': edges[bucket] if bucket < len(edges) else np.nan,
                              'n_genes': int(np.sum(codes == bucket))})
joint_codes = np.column_stack([bucket_codes[name] for name in matching_columns])
np.testing.assert_array_equal(np.unique(joint_codes, axis=0, return_inverse=True)[1], actual_strata)
bucket_counts = np.zeros(tuple(len(bucket_labels[name]) for name in matching_columns), dtype=int)
np.add.at(bucket_counts, tuple(bucket_codes[name] for name in matching_columns), 1)
assert bucket_counts.sum() == len(covariates)
pd.DataFrame(boundary_rows).to_csv(output / 'matching_bucket_boundaries.csv', index=False)
pd.DataFrame([{'expression_bucket': i + 1, 'detection_bucket': j + 1,
               'coverage_bucket': k + 1, 'n_genes': int(bucket_counts[i, j, k])}
              for i, j, k in np.ndindex(bucket_counts.shape)]).to_csv(
    output / 'matching_joint_bucket_counts.csv', index=False)

coverage_labels = bucket_labels['coverage']
fig, axes = plt.subplots(1, len(coverage_labels), figsize=(5.5 * len(coverage_labels), 4.7),
                         sharey=True, layout='constrained')
axes = np.atleast_1d(axes)
for k, ax in enumerate(axes):
    counts = bucket_counts[:, :, k]
    ax.imshow(counts, origin='lower', cmap='Blues', vmin=0,
              vmax=max(1, int(bucket_counts.max())))
    for i, j in np.ndindex(counts.shape):
        ax.text(j, i, f'{counts[i, j]:,}', ha='center', va='center',
                color='white' if counts[i, j] > bucket_counts.max() / 2 else '0.15')
    ax.set(xticks=np.arange(counts.shape[1]), yticks=np.arange(counts.shape[0]),
           xticklabels=bucket_labels['detection'], yticklabels=bucket_labels['mean_expression'],
           xlabel='Detection fraction',
           title=f'Coverage {coverage_labels[k]}  (n={counts.sum():,})')
    ax.tick_params(axis='x', labelrotation=25)
axes[0].set_ylabel('Mean log-normalized expression')
fig.suptitle(f'Gene counts in matching strata (n={len(covariates):,}; ' +
             f'{np.count_nonzero(bucket_counts)}/{bucket_counts.size} occupied)')
fig.savefig(output / 'matching_strata_matrix.pdf', bbox_inches='tight')
plt.show()
plt.close(fig)

# %%
from pseudospace.pathway_remodeling import matching_strata, matched_pathway_tests

n_null = 9999
seed = 12
strata = matching_strata(covariates, bins=3)
rankings = gene_stats[['T_total', 'T_level', 'T_spatial']]
matched_null_audit = cached_frame('matched_primary',
    lambda: matched_pathway_tests(rankings, gene_sets, strata, n_null=n_null, seed=seed),
    root=output / 'stage_cache', params={'n_null': n_null, 'seed': seed},
    inputs={'statistics': rankings, 'sets': gene_sets, 'strata': strata}, code=cache_code)
uniform_strata = np.zeros(len(strata), dtype=int)
uniform_null_audit = cached_frame('uniform_size_null',
    lambda: matched_pathway_tests(rankings, gene_sets, uniform_strata, n_null=n_null, seed=seed),
    root=output / 'stage_cache', params={'n_null': n_null, 'seed': seed},
    inputs={'statistics': rankings, 'sets': gene_sets, 'strata': uniform_strata},
    code=digest([cache_code, 'uniform-size-null-v1']))
uniform_null_audit.to_csv(output / 'pathway_uniform_size_null.csv', index=False)
null_comparison = matched_null_audit.merge(
    uniform_null_audit, on=['pathway_id', 'statistic'], suffixes=('_matched', '_uniform'),
    validate='one_to_one')
assert len(null_comparison) == len(matched_null_audit) == len(uniform_null_audit)
for column in ('n_genes', 'auc', 'effect'):
    np.testing.assert_allclose(null_comparison[f'{column}_matched'],
                               null_comparison[f'{column}_uniform'], rtol=0, atol=1e-12)
null_comparison.to_csv(output / 'pathway_null_comparison.csv', index=False)

spatial_null = null_comparison.loc[null_comparison.statistic.eq('T_spatial')].copy()
positive = spatial_null.effect_matched.gt(0)
spatial_null['matched_hit'] = positive & spatial_null.q_empirical_matched.le(.05)
spatial_null['uniform_hit'] = positive & spatial_null.q_empirical_uniform.le(.05)
spatial_null['cutoff_discordant'] = spatial_null.matched_hit.ne(spatial_null.uniform_hit)
null_summary = pd.Series({
    'tested_spatial_pathways': len(spatial_null),
    'matched_candidates': int(spatial_null.matched_hit.sum()),
    'uniform_candidates': int(spatial_null.uniform_hit.sum()),
    'both': int((spatial_null.matched_hit & spatial_null.uniform_hit).sum()),
    'matched_only': int((spatial_null.matched_hit & ~spatial_null.uniform_hit).sum()),
    'uniform_only': int((~spatial_null.matched_hit & spatial_null.uniform_hit).sum()),
    'median_absolute_q_difference': float((spatial_null.q_empirical_matched -
                                            spatial_null.q_empirical_uniform).abs().median()),
})
null_summary.rename('value').to_csv(output / 'pathway_null_comparison_summary.csv')
display(null_summary.to_frame('value'))

x = -np.log10(spatial_null.q_empirical_uniform.clip(lower=1 / (n_null + 1)))
y = -np.log10(spatial_null.q_empirical_matched.clip(lower=1 / (n_null + 1)))
fig, ax = plt.subplots(figsize=(6, 5), layout='constrained')
same = ~spatial_null.cutoff_discordant
ax.scatter(x[same], y[same], s=11, color='0.55', alpha=.35,
           rasterized=True, label='Same candidate status')
ax.scatter(x[~same], y[~same], s=24, color='#D55E00', alpha=.85,
           rasterized=True, label='Cutoff changed')
limit = max(float(x.max()), float(y.max()), -np.log10(.05)) + .2
ax.plot([0, limit], [0, limit], color='0.65', lw=.8)
ax.axvline(-np.log10(.05), color='0.35', ls='--', lw=.8)
ax.axhline(-np.log10(.05), color='0.35', ls='--', lw=.8)
ax.set(xlim=(0, limit), ylim=(0, limit), xlabel='−log10(q), uniform size null',
       ylabel='−log10(q), matched null', title='Does covariate matching change pathway evidence?')
ax.legend(frameon=False, fontsize=8, loc='upper left')
fig.savefig(output / 'pathway_null_comparison.pdf', bbox_inches='tight')
plt.show()
plt.close(fig)

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
# BH correction is **separate for each statistic** and spans every tested pathway across KEGG, Reactome, and Hallmark within that statistic. $T_{\rm spatial}$ is the primary discovery family; $T_{\rm level}$ and $T_{\rm total}$ are secondary characterization families. At the pathway level, the null is competitive enrichment of each gene statistic against matched genes, not a donor-level test of a species effect. Nominal Mann–Whitney p-values are shown only as diagnostics. The matching addresses measurability; it does **not** preserve inter-gene correlation or create biological replication.
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
family_summary = (pathway_tests.assign(hit=lambda frame: frame.q_empirical.le(.05) & frame.effect.gt(0))
    .groupby('statistic').agg(tested=('pathway_id', 'size'), exploratory_hits=('hit', 'sum')))
family_summary.to_csv(output / 'pathway_bh_family_summary.csv')
display(family_summary)
spatial = pathway_tests[pathway_tests.statistic.eq('T_spatial')].copy()
spatial = spatial.merge(coverage[['pathway_id', 'library', 'pathway', 'tested_fraction']], on='pathway_id')
spatial = spatial.sort_values(['q_empirical', 'effect', 'pathway_id'], ascending=[True, False, True])
display(spatial[['library', 'pathway', 'n_genes', 'effect', 'q_empirical', 'mc_se', 'tested_fraction']].head(15))

# %% [markdown]
# ### Primary spatial pathway screen
#
# Each dot is one tested pathway. Orange dots pass the predeclared positive-effect, matched-q discovery rule; labels identify only six leading terms. The y-axis is pathway-level evidence for **T_spatial**, not a local gridwise test. These findings remain exploratory with one human donor.
#

# %%
from textwrap import fill
fig, ax = plt.subplots(figsize=(8, 5), layout='constrained')
hits = spatial.q_empirical.le(.05) & spatial.effect.gt(0)
for mask, color, label, size in [(~hits, '0.75', 'Other tested pathways', 13),
                                  (hits, '#D55E00', 'Spatial candidates', 23)]:
    ax.scatter(spatial.loc[mask, 'effect'],
               -np.log10(spatial.loc[mask, 'q_empirical'].clip(lower=1e-300)),
               s=size, color=color, alpha=.7, label=f'{label} (n={mask.sum()})')
leading = spatial.loc[hits].head(6)
ax.scatter(leading.effect, -np.log10(leading.q_empirical.clip(lower=1e-300)),
           s=42, facecolors='none', edgecolors='black', linewidth=.8)
labels = '\n'.join(fill(row.pathway, width=36) for row in leading.itertuples())
ax.text(1.02, .98, 'Leading matched-q terms\n\n' + labels,
        transform=ax.transAxes, va='top', fontsize=8)
ax.axhline(-np.log10(.05), color='0.4', ls='--', lw=.8)
ax.axvline(0, color='0.4', lw=.8)
ax.set(xlabel='Spatial pathway effect (AUC − 0.5)', ylabel='−log10(matched empirical q)',
       title='Primary discovery: spatial trajectory remodeling')
ax.legend(frameon=False, loc='upper left')
fig.savefig(output / 'pathway_spatial_discovery.pdf', bbox_inches='tight')
plt.show()
plt.close(fig)


# %% [markdown]
# ### Visual check: what a pathway AUC measures
#
# For the top-ranked spatial pathway, these curves show where its tested member genes and the other eligible genes fall in the $T_{\rm spatial}$ ranking. A pathway shifted toward high percentiles has an AUC above 0.5. The displayed matched $q$ comes from the separate covariate-matched random-set test; the curves themselves are descriptive and do not show that null. The top row is selected for illustration, not as an additional discovery rule.

# %%
from textwrap import fill

example_pathway = spatial.iloc[0]
member_mask = gene_stats.index.isin(gene_sets[example_pathway.pathway_id])
percentile = gene_stats.T_spatial.rank(pct=True).to_numpy()
fig, ax = plt.subplots(figsize=(7, 4), layout='constrained')
for selected, label, color in [(~member_mask, 'Other eligible genes', '0.55'),
                               (member_mask, 'Pathway members', '#D55E00')]:
    values = np.sort(percentile[selected])
    ax.step(np.r_[0, values], np.r_[0, np.arange(1, len(values) + 1) / len(values)],
            where='post', color=color, lw=2, label=f'{label} (n={len(values):,})')
ax.set(xlim=(0, 1), ylim=(0, 1), xlabel='Percentile of spatial gene statistic',
       ylabel='Cumulative fraction of genes', title=f'Pathway rank distribution\n{fill(example_pathway.pathway, width=60)}')
ax.legend(frameon=False, loc='upper left',
          title=f"AUC = {example_pathway.auc:.2f}; matched q = {example_pathway.q_empirical:.2g}")
fig.savefig(output / 'pathway_auc_example.pdf', bbox_inches='tight')
plt.show()
plt.close(fig)

# %% [markdown]
# ### Freeze the candidates before inspecting local curves
#
# A discovery candidate has positive spatial enrichment and empirical **q ≤ 0.05 within the spatial pathway family**. This cutoff is a reproducible exploratory shortlist, not a confirmatory claim. If none pass, the notebook still exports the full screen and plots up to six highest-ranked terms clearly marked as illustrations. It never converts a top-ranked term into a discovery automatically.

# %%
candidate_ids = spatial.loc[spatial.q_empirical.le(.05) & spatial.effect.gt(0), 'pathway_id'].tolist()
illustration_ids = candidate_ids[:6] if candidate_ids else spatial.head(6).pathway_id.tolist()
print(f'{len(candidate_ids)} exploratory candidates.')
print('Plot selection:', 'discovery candidates' if candidate_ids else 'ranked illustrations; none passed the cutoff')

# %% [markdown]
# ### 6a · Correlation-aware competitive enrichment
#
# **Does a pathway's rank-sum evidence survive accounting for its genes moving together?** We keep the same $T_{\mathrm{spatial}}$, AUC, matched null, and matched q-values. This section adds a second significance calculation, then annotates the existing candidates.
#
# > **Remove the modeled signal before estimating correlation**
# >
# > Use residuals from the **full model**: $r_g=y_g-\widehat y_g$. Baseline position, species, species-by-position, and specimen effects have already been fitted. Within each specimen, calculate the mean Pearson correlation across distinct pathway-gene pairs. Fisher-transform those specimen means, average them **equally**, and transform back:
# >
# > $\bar\rho_P=\tanh\{J^{-1}\sum_j\operatorname{atanh}(\rho_{P,j})\}$, then $\rho_P^*=\max(0,\bar\rho_P)$.
#
# There is no pooled-structure correlation and no random gene or structure cap. The helper computes the exact mean pair correlation using normalized residual-column sums, without building a large pairwise matrix. Constant residual columns have undefined Pearson correlation and are excluded **only from correlation estimation**; the original pathway membership and AUC stay intact. The exported per-specimen audit records variable-gene counts and pair coverage. The available-pair estimate is applied to the original pathway size, an approximation to inspect when coverage is low. If a specimen has fewer than two variable members or four structures, the pathway is marked **correlation unavailable**, never assigned zero correlation silently.

# %%
from pseudospace.pathway_remodeling import residual_pathway_correlations, correlation_adjusted_rank_tests

residual_correlations = cached_frame('spatial_residual_correlations',
    lambda: residual_pathway_correlations(fit['residuals'], specimen, genes, gene_sets),
    root=output / 'stage_cache', inputs={'residuals': digest(fit['residuals']),
        'specimen': specimen, 'genes': genes, 'sets': gene_sets}, code=cache_code)
residual_correlations.to_csv(output / 'pathway_residual_correlations_by_specimen.csv', index=False)
display(residual_correlations.groupby('specimen').agg(
    pathways=('pathway_id', 'size'), estimable=('rho_specimen', 'count'),
    median_rho=('rho_specimen', 'median'), min_pair_coverage=('pair_fraction', 'min')))

# %% [markdown]
# #### Same rank sum, correlation-adjusted variance
#
# For $m$ pathway genes and $n=G-m$ background genes, keep $U=mn\,AUC$. Use the CAMERA rank variance:
#
# $$V_\rho=\frac{\arcsin(1)mn+\arcsin(1/2)mn(n-1)+\arcsin(\rho/2)m(m-1)n(n-1)+\arcsin((1+\rho)/2)m(m-1)n}{2\pi}.$$
#
# At $\rho=0$, this reduces to $mn(G+1)/12$. Following [limma's rank-test implementation](https://github.com/bioc/limma/blob/RELEASE_3_22/R/rankSumTestWithCorrelation.R), multiply this variance by the tie factor $1-\sum_k(t_k^3-t_k)/(G^3-G)$ and use the upper-tail continuity correction:
#
# $$Z_{\mathrm{corr}}=\frac{U-mn/2-0.5}{\sqrt{V_\rho\,\mathrm{tie\ factor}}},\qquad p_{\mathrm{corr}}=1-\Phi(Z_{\mathrm{corr}}).$$
#
# BH correction spans **all tested spatial pathways across all three libraries**, not just the candidates. Missing correlation estimates retain a conservative p=1 slot in that family but are displayed as unavailable; entirely tied rankings give p=1. The matched spatial q-values use their own full spatial-pathway family across the three libraries.
#
# > **Two complementary checks**
# >
# > Matched q-values address measurability; correlation q-values account approximately for within-pathway dependence. This is a **CAMERA-style adaptation to unsigned F statistics**, with specimen-balanced residual correlation estimates and a normal tail. It is not a direct `camera()`/`cameraPR()` analysis, a joint correction for both biases, or donor-level inference. Residual correlation is a working approximation to dependence among these F statistics. The one-human-donor and spatial-dependence limitations remain.
#
# We retain the existing candidates. Label them **Correlation-supported** when $q_{\rm corr}\le0.05$, **Correlation-sensitive** when it exceeds 0.05, and **Correlation unavailable** when estimation fails. A correlation-sensitive pathway can still have a large, biologically coherent effect; its genes offer fewer independent pieces of evidence.

# %%
correlation_tests = correlation_adjusted_rank_tests(gene_stats.T_spatial, gene_sets, residual_correlations)
# Enforce the central contract: changing significance must not change the pathway score.
auc_check = spatial.set_index('pathway_id').auc.reindex(correlation_tests.pathway_id)
assert np.allclose(auc_check, correlation_tests.auc, rtol=0, atol=1e-12)
correlation_columns = ['pathway_id', 'rho_residual', 'rho_used', 'p_corr', 'q_corr',
                       'variance_inflation', 'min_pair_fraction', 'n_specimens_estimable']
spatial = spatial.merge(correlation_tests[correlation_columns], on='pathway_id', validate='one_to_one')
spatial['correlation_support'] = np.select(
    [spatial.q_corr.isna(), spatial.q_corr.le(.05)],
    ['Correlation unavailable', 'Correlation-supported'], default='Correlation-sensitive')
correlation_tests.to_csv(output / 'pathway_spatial_correlation_tests.csv', index=False)
spatial.to_csv(output / 'pathway_spatial_evidence.csv', index=False)
# Membership of candidate_ids is still determined only by matched q and positive AUC effect.
correlation_candidates = spatial[spatial.pathway_id.isin(candidate_ids)]
display(correlation_candidates.groupby('correlation_support').size().rename('existing_candidates'))
display(correlation_candidates[['pathway', 'auc', 'q_empirical', 'rho_residual',
                                'q_corr', 'correlation_support', 'min_pair_fraction']].head(15))

# %% [markdown]
# #### Diagnostic: matched evidence versus correlation-adjusted evidence
#
# Every dot is a tested spatial pathway. The dashed lines mark q=0.05; color shows the nonnegative residual correlation used in the variance calculation. The diagonal is an agreement guide, not a null expectation: these are different tests and BH families. Points right of the vertical cutoff with positive AUC effect remain candidates even when below the horizontal cutoff. Unavailable correlation estimates are counted separately. This plot diagnoses sensitivity to the working correlation correction; it does not establish population-level calibration.

# %%
import matplotlib.pyplot as plt

finite_correlation = spatial.dropna(subset=['q_corr'])
fig, ax = plt.subplots(figsize=(6.5, 5.5), layout='constrained')
points = ax.scatter(-np.log10(finite_correlation.q_empirical.clip(lower=1e-300)),
                    -np.log10(finite_correlation.q_corr.clip(lower=1e-300)),
                    c=finite_correlation.rho_used, cmap='viridis', vmin=0, s=18, alpha=.8)
for support, marker, edge in [('Correlation-supported', 'o', '#D55E00'),
                               ('Correlation-sensitive', 's', '#0072B2')]:
    marked = correlation_candidates[correlation_candidates.correlation_support.eq(support)]
    ax.scatter(-np.log10(marked.q_empirical.clip(lower=1e-300)),
               -np.log10(marked.q_corr.clip(lower=1e-300)),
               s=42, facecolors='none', edgecolors=edge, marker=marker,
               linewidth=1.1, label=f'{support} candidates (n={len(marked)})')
ax.legend(frameon=False, fontsize=8)
cutoff = -np.log10(.05)
ax.axvline(cutoff, color='0.4', ls='--', lw=.8)
ax.axhline(cutoff, color='0.4', ls='--', lw=.8)
limit = max(*ax.get_xlim(), *ax.get_ylim(), cutoff + .2)
ax.plot([0, limit], [0, limit], color='0.7', lw=.8)
ax.set(xlim=(0, limit), ylim=(0, limit), xlabel='−log10(matched empirical q)',
       ylabel='−log10(correlation-adjusted q)', title='Spatial pathway evidence: two complementary checks')
fig.colorbar(points, ax=ax, label='Residual correlation used (negative values floored at 0)')
fig.savefig(output / 'pathway_correlation_diagnostic.pdf', bbox_inches='tight')
plt.show()
plt.close(fig)
print('Correlation estimates unavailable:', int(spatial.q_corr.isna().sum()))


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
# ### AUC–GSEA concordance for spatial remodeling
#
# Both screens use the same $T_{\rm spatial}$ ranking and pathway memberships. Compare their effect summaries for every pathway returned by both methods. Spearman correlation uses pathways with **positive NES**; the scatterplot retains both NES directions. The 2×2 table counts positive-effect AUC discoveries under the matched spatial BH family and positive-NES GSEA findings under its existing BH family (all pathway–statistic pairs). These different nulls and correction families can give different significance calls. Agreement is descriptive and does not enter the primary candidate rule.
#

# %%
auc_gsea = spatial[['pathway_id', 'effect', 'q_empirical']].merge(
    gsea.loc[gsea.statistic.eq('T_spatial'), ['pathway_id', 'NES', 'q_bh_family']],
    on='pathway_id', how='inner', validate='one_to_one')
auc_gsea['NES'] = pd.to_numeric(auc_gsea.NES, errors='coerce')
auc_gsea = auc_gsea.loc[np.isfinite(auc_gsea.effect) & np.isfinite(auc_gsea.NES)].copy()
auc_gsea['auc_hit'] = auc_gsea.effect.gt(0) & auc_gsea.q_empirical.le(.05)
auc_gsea['gsea_hit'] = auc_gsea.NES.gt(0) & auc_gsea.q_bh_family.le(.05)
positive_nes = auc_gsea.loc[auc_gsea.NES.gt(0)]
rho = positive_nes.effect.corr(positive_nes.NES, method='spearman') if len(positive_nes) >= 2 else np.nan
print(f'Shared spatial pathways: {len(auc_gsea)} / {len(spatial)} AUC pathways; '
      f'positive-NES Spearman rho: {rho:.3f} (n={len(positive_nes)})')
overlap = pd.crosstab(auc_gsea.auc_hit, auc_gsea.gsea_hit).reindex(
    index=[False, True], columns=[False, True], fill_value=0)
overlap.index.name = 'AUC: positive effect and matched q ≤ 0.05'
overlap.columns.name = 'GSEA: positive NES and BH q ≤ 0.05'
assert overlap.to_numpy().sum() == len(auc_gsea)
display(overlap)
fig, ax = plt.subplots(figsize=(6, 4.5))
groups = [
    (~auc_gsea.auc_hit & ~auc_gsea.gsea_hit, 'Neither', '0.75'),
    (auc_gsea.auc_hit & ~auc_gsea.gsea_hit, 'AUC only', '#D55E00'),
    (~auc_gsea.auc_hit & auc_gsea.gsea_hit, 'GSEA only', '#009E73'),
    (auc_gsea.auc_hit & auc_gsea.gsea_hit, 'Both', '#0072B2'),
]
for mask, label, color in groups:
    ax.scatter(auc_gsea.loc[mask, 'effect'], auc_gsea.loc[mask, 'NES'],
               color=color, s=14 if label == 'Neither' else 25,
               alpha=.55 if label == 'Neither' else .85, label=f'{label} (n={mask.sum()})')
ax.axvline(0, color='0.35', ls='--', lw=.8)
ax.axhline(0, color='0.35', ls='--', lw=.8)
ax.set(xlabel='Spatial AUC − 0.5', ylabel='Unweighted GSEA NES',
       title='Spatial pathway AUC–GSEA concordance')
ax.legend(frameon=False)
fig.savefig(output / 'pathway_auc_gsea_concordance.pdf', bbox_inches='tight')
plt.show()
plt.close(fig)


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
# ### Where and how broadly do candidate pathways diverge?
#
# Peak position and width describe the existing D(s) curves. Color shows the median **absolute signed gene Z** at that peak: orange is human-high, blue is mouse-high. Early and late refer only to the shared PT coordinate, not anatomical S1/S2/S3 segments. Point size measures the fraction of the grid at or above half the positive peak, not significance.
#

# %%
if len(phenotypes):
    fig, ax = plt.subplots(figsize=(8, 5), layout='constrained')
    direction_limit = max(1, np.nanmax(np.abs(phenotypes.median_member_z_at_peak)))
    points = ax.scatter(phenotypes.peak_position, phenotypes.peak_divergence,
                        s=25 + 140 * phenotypes.affected_grid_fraction,
                        c=phenotypes.median_member_z_at_peak, cmap='RdBu_r',
                        vmin=-direction_limit, vmax=direction_limit,
                        edgecolor='0.3', linewidth=.3, alpha=.8)
    for pathway in illustration_ids[:5]:
        row = phenotypes.set_index('pathway_id').loc[pathway]
        ax.annotate(fill(pathway.split('::', 1)[-1], width=25),
                    (row.peak_position, row.peak_divergence),
                    xytext=(4, 4), textcoords='offset points', fontsize=7)
    ax.set(xlim=(grid[0], grid[-1]), xlabel='Established shared PT pseudospace',
           ylabel='Peak D(s): divergence rank-AUC − 0.5',
           title='Spatial location and breadth of pathway remodeling')
    fig.colorbar(points, ax=ax, label='Median member Z at peak (human − mouse)')
    fig.savefig(output / 'pathway_spatial_phenotypes.pdf', bbox_inches='tight')
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
#
# Library sizes sum counts over the **full measured ortholog panel before the PT detection filter**, matching section 3's denominator definition. The inherited notebook 13 embedding gene filters play no role in this benchmark.

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
# “High” below uses positive effect and matched q ≤ 0.05. In the exported atlas, `q_spatial`, `q_level`, and `q_total` are the empirical BH values from their separate pathway families. “Low” means not passing this exploratory cutoff, **not proof of no effect**. A pathway found by both methods remains valuable: report level/bulk and spatial AUC effects, peak location, local direction, sign changes, and genes driving the peak. Different q-value statuses do not test whether one effect is statistically stronger than another; compare the effect sizes directly.

# %%
all_tests = pd.concat([pathway_tests, bulk_tests], ignore_index=True)
comparison = all_tests.pivot(index='pathway_id', columns='statistic', values=['effect', 'q_empirical'])
comparison.columns = ['_'.join(c) for c in comparison.columns]
comparison = comparison.reset_index()
for name in ('spatial', 'level', 'total'):
    comparison[f'q_{name}'] = comparison[f'q_empirical_T_{name}']
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
# ### What does continuous position reveal beyond pseudobulk?
#
# The scatter compares the existing bulk and spatial AUC effects; color uses the notebook's frozen comparison groups, which also consider the level screen. The count panel keeps the large weak-evidence background visible without hiding the smaller groups. A bulk/level negative screen does not prove a zero overall species effect.
#

# %%
group_colors = {'continuous-only candidate': '#D55E00',
                'bulk/level plus spatial structure': '#0072B2',
                'bulk/level captures difference': '#009E73', 'weak evidence': '0.75'}
group_order = list(group_colors)
fig, axes = plt.subplots(1, 2, figsize=(11, 5), gridspec_kw={'width_ratios': [2, 1]},
                         layout='constrained')
for group in reversed(group_order):
    rows = comparison[comparison.comparison_group.eq(group)]
    axes[0].scatter(rows.effect_T_bulk, rows.effect_T_spatial, s=15 if group == 'weak evidence' else 25,
                    alpha=.55 if group == 'weak evidence' else .8,
                    color=group_colors[group], label=group)
axes[0].axhline(0, color='0.5', lw=.7)
axes[0].axvline(0, color='0.5', lw=.7)
axes[0].set(xlabel='Pseudobulk pathway effect (AUC − 0.5)',
            ylabel='Spatial pathway effect (AUC − 0.5)',
            title='Continuous versus conventional evidence')
axes[0].legend(frameon=False, fontsize=8)
group_counts = comparison.comparison_group.value_counts().reindex(group_order, fill_value=0)
axes[1].barh(np.arange(len(group_counts)), group_counts, color=[group_colors[g] for g in group_order])
axes[1].set(yticks=np.arange(len(group_counts)), yticklabels=group_order,
            xlabel='Tested pathways (log scale)', title='Screen categories')
axes[1].set_xscale('symlog', linthresh=1)
axes[1].invert_yaxis()
for i, count in enumerate(group_counts):
    axes[1].annotate(str(count), (max(count, 1), i), xytext=(4, 0),
                     textcoords='offset points', va='center', fontsize=8)
fig.savefig(output / 'pathway_continuous_vs_bulk_overview.pdf', bbox_inches='tight')
plt.show()
plt.close(fig)


# %% [markdown]
# ### Unsigned PT-S1, PT-S2, and PT-S3 cluster diagnostic
#
# Notebook 13 confirmed the reviewed Leiden clusters labeled `PT-S1`, `PT-S2`, and `PT-S3`. Use those saved `segment_class` labels on the **same common-support structures** as the continuous model and whole-PT pseudobulk. These are cluster labels, not equal-width pseudospace bins. Keep the saved coordinate and labels unchanged; report the number of structures per specimen and cluster.
#
# Within each specimen and cluster, sum raw counts and calculate log2(CPM + 1), using **all accepted orthologs measured in every specimen** for that cluster library's denominator. The same PT-wide 2% detection rule defines eligible genes, and the same pathway members and matching strata are used throughout. For each cluster, rank genes by the absolute Welch statistic between the two mouse specimens and two human cortex sections, then run the matched rank-AUC and unweighted preranked GSEA screens used for whole-PT pseudobulk. Signed mean log2(CPM + 1) differences retain direction. Each cluster's matched BH correction spans all tested pathways across the three libraries; GSEA BH spans all three cluster screens together.
#
# These cluster-level rankings are **descriptive**: two human sections come from one donor, and cluster membership was learned from the joint dataset. A cluster-specific hit does not prove a biological segment effect or test whether its effect differs from another cluster. Compare the pathway effects and gene directions with the existing continuous spatial screen; the primary discovery list is unchanged. Report nominal matched-test signals and GSEApy’s per-cluster FDR alongside the stricter matched BH and pooled GSEA BH results. The latter GSEA correction uses a 0.001 p-value floor from 999 permutations; near-zero counts can reflect that resolution. GSEApy FDR covers pathways within each cluster run, not the three clusters jointly.

# %%
segments = ('PT-S1', 'PT-S2', 'PT-S3')
if 'segment_class' not in adata.obs:
    raise ValueError('Notebook 13 PT artifact lacks reviewed segment_class labels.')
labels = adata.obs.segment_class.astype(str).to_numpy()
if set(labels) != set(segments):
    raise ValueError(f'Expected only reviewed PT-S1/S2/S3 labels, got {sorted(set(labels))}.')
segment_counts = (pd.DataFrame({'specimen': specimen, 'segment_class': labels})
    .groupby(['specimen', 'segment_class']).size().unstack(fill_value=0)
    .reindex(index=bulk_all.index, columns=segments, fill_value=0))
if segment_counts.lt(20).any().any():
    raise ValueError(f'Fewer than 20 structures in a reviewed cluster and specimen:\n{segment_counts}')
segment_counts.to_csv(output / 'reviewed_pt_cluster_counts.csv')
display(segment_counts)

segment_raw_parts = []
for label in segments:
    keep = labels == label
    part = pseudobulk_profiles(counts[keep], specimen[keep], position[keep], n_bins=1,
        gene_names=adata.var_names, min_structures_per_bin=20).set_index('specimen')
    segment_raw_parts.append(part.reindex(bulk_all.index))
segment_raw = pd.concat(segment_raw_parts, keys=segments, names=['segment_class', 'specimen'])
count_columns = ['count_' + gene for gene in adata.var_names]
np.testing.assert_array_equal(segment_raw[count_columns].groupby('specimen').sum().loc[bulk_all.index].to_numpy(),
                              bulk_all[count_columns].to_numpy())
segment_library = segment_raw[['count_' + gene for gene in measured_genes]].sum(axis=1)
if segment_library.le(0).any():
    raise ValueError('An empty reviewed-cluster pseudobulk library cannot be normalized.')
segment_expression = np.log2(segment_raw[['count_' + gene for gene in genes]].div(segment_library, axis=0) * 1e6 + 1)
segment_expression.columns = genes
segment_expression.to_csv(output / 'reviewed_pt_cluster_pseudobulk_logcpm.csv')

segment_gene_parts = []
segment_rankings = pd.DataFrame(index=pd.Index(genes, name='gene'))
for label in segments:
    block = segment_expression.loc[label]
    mouse_block = block.loc[expected['mouse']]
    human_block = block.loc[expected['human']]
    effect = human_block.mean() - mouse_block.mean()
    se = np.sqrt(human_block.var(ddof=1) / 2 + mouse_block.var(ddof=1) / 2).clip(lower=1e-6)
    statistic = 'T_' + label
    segment_rankings[statistic] = (effect / se).abs().reindex(genes)
    segment_gene_parts.append(pd.DataFrame({'segment_class': label, 'gene': genes,
        'human_minus_mouse_logcpm': effect.reindex(genes).to_numpy(),
        'T_welch_abs': segment_rankings[statistic].to_numpy()}))
segment_gene_stats = pd.concat(segment_gene_parts, ignore_index=True)
segment_gene_stats.to_csv(output / 'reviewed_pt_cluster_gene_statistics.csv', index=False)
segment_tests = cached_frame('matched_reviewed_pt_clusters',
    lambda: matched_pathway_tests(segment_rankings, gene_sets, strata, n_null=n_null, seed=seed),
    root=output / 'stage_cache', params={'n_null': n_null, 'seed': seed, 'segments': segments},
    inputs={'statistics': segment_rankings, 'sets': gene_sets, 'strata': strata}, code=cache_code)
segment_tests['segment_class'] = segment_tests.statistic.str.removeprefix('T_')
segment_tests.to_csv(output / 'reviewed_pt_cluster_pathway_rank_auc.csv', index=False)
segment_gsea = cached_frame('reviewed_pt_cluster_gsea', lambda: run_gsea(segment_rankings),
    root=output / 'stage_cache',
    params={'seed': seed, 'permutations': gsea_permutations, 'weight': 0,
            'gseapy': version('gseapy'), 'segments': segments},
    inputs={'statistics': segment_rankings, 'sets': gene_sets}, code=cache_code)
segment_gsea['segment_class'] = segment_gsea.statistic.str.removeprefix('T_')
segment_gsea.to_csv(output / 'reviewed_pt_cluster_pathway_gsea.csv', index=False)

segment_comparison = comparison[['pathway_id', 'effect_T_spatial', 'q_empirical_T_spatial',
                                 'effect_T_bulk', 'q_empirical_T_bulk']].copy()
for label in segments:
    rows = segment_tests.loc[segment_tests.segment_class.eq(label),
                             ['pathway_id', 'effect', 'p_empirical', 'q_empirical']]
    segment_comparison = segment_comparison.merge(rows.rename(columns={
        'effect': f'effect_{label}', 'p_empirical': f'p_{label}',
        'q_empirical': f'q_{label}'}), on='pathway_id', validate='one_to_one')
    gsea_rows = segment_gsea.loc[segment_gsea.segment_class.eq(label),
                                 ['pathway_id', 'NES', 'FDR q-val', 'q_bh_family']]
    segment_comparison = segment_comparison.merge(gsea_rows.rename(columns={
        'NES': f'NES_{label}', 'FDR q-val': f'gsea_fdr_{label}',
        'q_bh_family': f'gsea_family_q_{label}'}), on='pathway_id', validate='one_to_one')
segment_comparison['primary_spatial_candidate'] = segment_comparison.pathway_id.isin(candidate_ids)
segment_comparison.to_csv(output / 'pathway_continuous_vs_reviewed_pt_clusters.csv', index=False)
segment_screen_summary = (segment_tests.assign(
    matched_nominal=lambda frame: frame.p_empirical.le(.05) & frame.effect.gt(0),
    matched_hit=lambda frame: frame.q_empirical.le(.05) & frame.effect.gt(0))
    .groupby('segment_class').agg(tested=('pathway_id', 'size'),
        matched_nominal=('matched_nominal', 'sum'), matched_hits=('matched_hit', 'sum')))
gsea_screen_counts = (segment_gsea.assign(
    gsea_package_fdr_hit=lambda frame: pd.to_numeric(frame['FDR q-val']).le(.05) & pd.to_numeric(frame.NES).gt(0),
    gsea_family_hit=lambda frame: frame.q_bh_family.le(.05) & pd.to_numeric(frame.NES).gt(0))
    .groupby('segment_class').agg(gsea_package_fdr_hits=('gsea_package_fdr_hit', 'sum'),
                                  gsea_family_bh_hits=('gsea_family_hit', 'sum')))
segment_screen_summary = segment_screen_summary.join(gsea_screen_counts)
segment_screen_summary['spatial_candidate_overlap'] = 0
for label in segments:
    segment_screen_summary.loc[label, 'spatial_candidate_overlap'] = int((
        segment_comparison.primary_spatial_candidate & segment_comparison[f'q_{label}'].le(.05)
        & segment_comparison[f'effect_{label}'].gt(0)).sum())
segment_screen_summary.to_csv(output / 'reviewed_pt_cluster_screen_summary.csv')
display(segment_screen_summary)
display(segment_comparison.loc[segment_comparison.primary_spatial_candidate]
    .sort_values('q_empirical_T_spatial').head(12))


# %% [markdown]
# ### What does the unsigned cluster diagnostic miss?
#
# Keep the primary matched rank-AUC criterion visible: positive effect and empirical BH q ≤ 0.05. Also count positive GSEA results at GSEApy’s **per-cluster** FDR q ≤ 0.05. The named continuous-only export contains spatial candidates with neither adjusted cluster signal; the full comparison table retains the separate matched and GSEA calls. Nominal matched p-values are shown as exploratory rankings, not adjusted hits. The pooled three-cluster GSEA BH column remains available but is resolution-limited by 999 permutations.
#
# These tests have different nulls and correction families. A cluster screen missing either cutoff is not evidence that the pathway is absent from that segment, nor a test that the continuous and cluster effects differ. With two mice and two cortex sections from one human donor, inspect the effect sizes and member-gene curves before interpreting a continuous-only label.
#
# This unsigned comparison ranks the *magnitude* of human–mouse differences. Its continuous-only export is retained for provenance; the signed count-model benchmark below is the conventional comparison for the paper.
#

# %%
cluster_screen = segment_comparison.copy()
cluster_screen['continuous_hit'] = cluster_screen.pathway_id.isin(candidate_ids)
for label in segments:
    cluster_screen[f'hit_{label}'] = (cluster_screen[f'effect_{label}'].gt(0)
                                       & cluster_screen[f'q_{label}'].le(.05))
cluster_screen['any_cluster_hit'] = cluster_screen[[f'hit_{label}' for label in segments]].any(axis=1)
cluster_screen['min_cluster_q'] = cluster_screen[[f'q_{label}' for label in segments]].min(axis=1)
cluster_effects = cluster_screen[[f'effect_{label}' for label in segments]]
cluster_screen['max_cluster_effect'] = cluster_effects.max(axis=1)
cluster_screen['strongest_cluster'] = cluster_effects.idxmax(axis=1).str.removeprefix('effect_')
cluster_screen['any_cluster_gsea_hit'] = np.column_stack([
    cluster_screen[f'NES_{label}'].gt(0) & cluster_screen[f'gsea_family_q_{label}'].le(.05)
    for label in segments]).any(axis=1)
cluster_screen['any_cluster_gsea_package_hit'] = np.column_stack([
    cluster_screen[f'NES_{label}'].gt(0) & cluster_screen[f'gsea_fdr_{label}'].le(.05)
    for label in segments]).any(axis=1)
cluster_screen['any_cluster_adjusted_signal'] = (cluster_screen.any_cluster_hit
    | cluster_screen.any_cluster_gsea_package_hit)
cluster_screen['continuous_only_matched_screen'] = (cluster_screen.continuous_hit
    & ~cluster_screen.any_cluster_hit)
cluster_screen['comparison_group_vs_clusters'] = np.select(
    [cluster_screen.continuous_hit & ~cluster_screen.any_cluster_adjusted_signal,
     cluster_screen.continuous_hit & cluster_screen.any_cluster_adjusted_signal,
     ~cluster_screen.continuous_hit & cluster_screen.any_cluster_adjusted_signal],
    ['continuous-only across adjusted screens', 'continuous and cluster signal',
     'cluster-only signal'], default='neither adjusted screen')
assert cluster_screen.pathway_id.is_unique
assert cluster_screen.continuous_hit.equals(cluster_screen.primary_spatial_candidate)
assert len(cluster_screen) == len(segment_comparison)
cluster_screen.to_csv(output / 'pathway_continuous_vs_reviewed_pt_clusters.csv', index=False)
continuous_only = cluster_screen.loc[cluster_screen.comparison_group_vs_clusters.eq(
    'continuous-only across adjusted screens')].sort_values(
    ['q_empirical_T_spatial', 'effect_T_spatial'], ascending=[True, False])
assert (continuous_only.continuous_hit & ~continuous_only.any_cluster_adjusted_signal).all()
continuous_only.to_csv(output / 'pathway_continuous_only_vs_pt_clusters.csv', index=False)
cluster_screen_summary = (cluster_screen.groupby('comparison_group_vs_clusters', observed=True)
    .size().reindex(['continuous-only across adjusted screens', 'continuous and cluster signal',
                     'cluster-only signal', 'neither adjusted screen'], fill_value=0)
    .rename('pathways'))
cluster_screen_summary.to_csv(output / 'pathway_continuous_vs_pt_clusters_summary.csv')
display(cluster_screen_summary)
print('Spatial candidates with no matched cluster hit:', int(cluster_screen.continuous_only_matched_screen.sum()))
print('Spatial candidates also found by per-cluster GSEA FDR:', int((
    cluster_screen.continuous_hit & cluster_screen.any_cluster_gsea_package_hit).sum()))
display(continuous_only[['pathway_id', 'effect_T_spatial', 'q_empirical_T_spatial',
                         'strongest_cluster', 'max_cluster_effect', 'min_cluster_q',
                         'any_cluster_gsea_package_hit']].head(20).round(3))
cluster_gsea_top = (segment_gsea.loc[pd.to_numeric(segment_gsea.NES).gt(0)]
    .sort_values(['segment_class', 'FDR q-val', 'NES'], ascending=[True, True, False])
    .groupby('segment_class', sort=False).head(8)
    [['segment_class', 'pathway_id', 'NES', 'FDR q-val', 'q_bh_family']])
display(cluster_gsea_top)

# %% [markdown]
# ## 9a · Conventional PT analysis: signed pseudobulk DE
#
# The reviewed PT-S1/S2/S3 assignments define two study-style contrasts on the same common-support structures:
#
# 1. **Human versus mouse within each segment.** Sum raw counts by specimen × segment; fit a separate DESeq2 negative-binomial model (`~ species`) to each segment and rank every eligible gene by its **signed Wald statistic**. Positive log2 fold change and statistic mean human-high. This is the discrete comparator for continuous spatial remodeling.
# 2. **Segment versus the other PT segments.** Within each specimen, compare the target segment with the sum of the other two using a paired `~ specimen + segment_group` design. Positive values mean segment-enriched. These contrasts characterize the reviewed cluster labels; they are not human–mouse tests.
#
# A fixed filter requires at least 10 raw counts in at least two pseudobulks **within each of S1, S2, and S3**. This produces one common gene universe for all three species contrasts and for the continuous reanalysis below. Retain the original PT-wide universe and candidate list separately. Every fitted count is an original measured ortholog; no Harmony expression or normalized values enter DESeq2. Counts, sample depth, DEG tables, and filtering decisions are exported. The existing logCPM profiles also record whether both human sections lie above or below both mouse specimens for each model-selected gene; this is a section-level consistency check, not independent-donor replication. Reactome, Hallmark, and KEGG are the same pathway libraries used by the continuous screen.
#
# **Inference limit:** the two human sections come from one donor. DESeq2 estimates specimen/section-level variation but cannot estimate variation among human donors here. Its Wald p-values, DEG FDR values, GSEA results, and ORA results are descriptive model outputs, not population-level evidence for a species effect. Cluster assignment was learned on these same data. Interpret effects, replicate consistency, and gene membership before biological claims. [DESeq2 method](https://pmc.ncbi.nlm.nih.gov/articles/PMC4302049/); [pseudobulk rationale](https://www.nature.com/articles/s41467-021-25960-2); [PyDESeq2 workflow](https://pydeseq2.readthedocs.io/en/v0.5.4/auto_examples/plot_minimal_pydeseq2_pipeline.html).

# %%
import pseudospace.conventional_pt as conventional_pt
from pseudospace.conventional_pt import deseq2_contrast

raw_blocks = {}
for label in segments:
    block = segment_raw.loc[label, ['count_' + gene for gene in genes]].copy()
    block.columns = genes
    block = block.loc[bulk_all.index]
    if not np.equal(block.to_numpy(), np.floor(block.to_numpy())).all():
        raise ValueError(f'{label}: pseudobulk has noninteger counts.')
    raw_blocks[label] = block.astype(np.int64)

common_count_support = np.logical_and.reduce([
    raw_blocks[label].ge(10).sum(axis=0).ge(2).to_numpy() for label in segments])
conventional_genes = pd.Index(genes)[common_count_support]
if len(conventional_genes) < 1000:
    raise ValueError('Too few genes have count support in all reviewed PT segments.')
conventional_gene_sets = {pathway: [gene for gene in members if gene in conventional_genes]
    for pathway, members in gene_sets.items()}
conventional_gene_sets = {pathway: members for pathway, members in conventional_gene_sets.items()
                          if 10 <= len(members) <= 300}
if not conventional_gene_sets:
    raise ValueError('No pathways remain in the common conventional gene universe.')
conventional_coverage = pd.DataFrame({'pathway_id': list(gene_sets),
    'original_members': [len(gene_sets[pathway]) for pathway in gene_sets],
    'shared_members': [len(set(gene_sets[pathway]) & set(conventional_genes)) for pathway in gene_sets]})
conventional_coverage['tested'] = conventional_coverage.pathway_id.isin(conventional_gene_sets)
conventional_coverage.to_csv(output / 'conventional_pathway_coverage.csv', index=False)
(pd.DataFrame({'gene': genes, 'common_count_support': common_count_support})
    .to_csv(output / 'conventional_gene_universe.csv', index=False))

sample_qc = (segment_counts.stack().rename('n_structures').reset_index()
    .rename(columns={'level_1': 'segment_class'}))
sample_qc['species'] = sample_qc.specimen.map({name: species for species, names in expected.items()
                                                for name in names})
sample_qc['library_counts'] = [int(raw_blocks[row.segment_class].loc[row.specimen].sum())
                               for row in sample_qc.itertuples()]
sample_qc['detected_genes'] = [int(raw_blocks[row.segment_class].loc[row.specimen].gt(0).sum())
                               for row in sample_qc.itertuples()]
assert sample_qc.species.notna().all() and sample_qc.library_counts.gt(0).all()
sample_qc.to_csv(output / 'conventional_pseudobulk_qc.csv', index=False)
display(sample_qc)
print(f'Common DE universe: {len(conventional_genes):,} genes; '
      f'{len(conventional_gene_sets):,} comparable pathway sets.')

# %%
conventional_code = digest([cache_code, Path(conventional_pt.__file__),
                            'signed-deseq2-v1', version('pydeseq2')])
conventional_de_parts = []
for label in segments:
    model_counts = raw_blocks[label][conventional_genes]
    model_metadata = pd.DataFrame({'species': [
        'mouse' if name in expected['mouse'] else 'human' for name in model_counts.index]},
        index=model_counts.index)
    result = cached_frame('conventional_species_' + label,
        lambda counts=model_counts, metadata=model_metadata: deseq2_contrast(
            counts, metadata, design='~species', contrast=('species', 'human', 'mouse')),
        root=output / 'stage_cache', params={'design': '~species', 'contrast': 'human_vs_mouse',
                                            'version': version('pydeseq2')},
        inputs={'counts': model_counts, 'metadata': model_metadata}, code=conventional_code)
    result['scope'] = 'species_within_segment'
    result['segment_class'] = label
    conventional_de_parts.append(result)

    other = [part for part in segments if part != label]
    rest_counts = raw_blocks[other[0]][conventional_genes] + raw_blocks[other[1]][conventional_genes]
    paired_counts = pd.concat([model_counts, rest_counts], keys=['target', 'rest'],
                              names=['segment_group', 'specimen'])
    paired_metadata = paired_counts.index.to_frame(index=False)
    paired_metadata.index = [f'{specimen}|{group}' for group, specimen in paired_counts.index]
    paired_counts.index = paired_metadata.index
    marker_result = cached_frame('conventional_marker_' + label,
        lambda counts=paired_counts, metadata=paired_metadata: deseq2_contrast(
            counts, metadata, design='~specimen + segment_group',
            contrast=('segment_group', 'target', 'rest')),
        root=output / 'stage_cache', params={'design': '~specimen + segment_group',
            'contrast': 'target_vs_rest', 'version': version('pydeseq2')},
        inputs={'counts': paired_counts, 'metadata': paired_metadata}, code=conventional_code)
    marker_result['scope'] = 'segment_vs_rest'
    marker_result['segment_class'] = label
    conventional_de_parts.append(marker_result)

conventional_de = pd.concat(conventional_de_parts, ignore_index=True)
assert conventional_de.groupby(['scope', 'segment_class']).gene.nunique().eq(len(conventional_genes)).all()
conventional_de['logcpm_complete_pair_separation'] = pd.Series(pd.NA, index=conventional_de.index, dtype='boolean')
for label in segments:
    block = segment_expression.loc[label, conventional_genes]
    human_above = block.loc[expected['human']].min().gt(block.loc[expected['mouse']].max())
    mouse_above = block.loc[expected['mouse']].min().gt(block.loc[expected['human']].max())
    mask = conventional_de.scope.eq('species_within_segment') & conventional_de.segment_class.eq(label)
    rows = conventional_de.loc[mask]
    conventional_de.loc[mask, 'logcpm_complete_pair_separation'] = np.where(
        rows.log2FoldChange.gt(0), human_above.reindex(rows.gene).to_numpy(),
        mouse_above.reindex(rows.gene).to_numpy())
conventional_de.to_csv(output / 'conventional_pt_deseq2_gene_statistics.csv', index=False)
species_de = conventional_de.loc[conventional_de.scope.eq('species_within_segment')].copy()
species_de['selected_gene'] = species_de.padj.le(.05) & species_de.log2FoldChange.abs().ge(1)
separation_audit = species_de.loc[species_de.selected_gene].groupby('segment_class').agg(
    selected_genes=('gene', 'size'),
    complete_pair_separation=('logcpm_complete_pair_separation', 'sum'),
    fraction_complete_pair_separation=('logcpm_complete_pair_separation', 'mean'))
separation_audit.to_csv(output / 'conventional_species_pair_separation.csv')
display(separation_audit)
de_summary = conventional_de.assign(
    human_or_target_high=lambda frame: frame.padj.le(.05) & frame.log2FoldChange.ge(1),
    mouse_or_rest_high=lambda frame: frame.padj.le(.05) & frame.log2FoldChange.le(-1))
de_summary = de_summary.groupby(['scope', 'segment_class']).agg(
    genes=('gene', 'size'), tested=('stat', 'count'),
    higher=('human_or_target_high', 'sum'), lower=('mouse_or_rest_high', 'sum'))
de_summary.to_csv(output / 'conventional_pt_de_summary.csv')
display(de_summary)

# %% [markdown]
# ### Directional pathway enrichment
#
# For each signed Wald ranking, use weighted preranked GSEA (`weight=1`) so magnitude and direction both contribute. GSEApy's multilevel procedure resolves p-values below the 999-permutation floor of the earlier unsigned diagnostic. BH correction spans **all pathways and all three segment contrasts within each contrast type** (species or segment identity); retain GSEApy's per-run FDR as a separate column. Positive NES is human-high for the species test and segment-enriched for the marker test; negative NES reverses that direction. The same shared gene/pathway universe is used for all six rankings.
#
# A second, thresholded view runs hypergeometric ORA on DESeq2 genes with adjusted gene p ≤ 0.05 and |log2FC| ≥ 1, separately for each direction. The tested-gene universe is the background. ORA is supplementary because its DEG threshold and the present replication can strongly affect the list. Its BH family covers all pathway × segment × direction tests within a contrast type. [GSEApy preranked methods](https://gseapy.readthedocs.io/en/latest/gseapy_example.html).

# %%
from scipy.stats import hypergeom

def conventional_gsea(de_rows):
    ranking = de_rows[['gene', 'stat']].dropna().sort_values(
        ['stat', 'gene'], ascending=[False, True], kind='stable')
    if ranking.gene.duplicated().any() or len(ranking) < 1000:
        raise ValueError('Signed GSEA needs unique genes and adequate testable coverage.')
    fit_result = gp.prerank(rnk=ranking, gene_sets=conventional_gene_sets,
        min_size=10, max_size=300, weight=1, method='multilevel',
        seed=seed, threads=2, outdir=None, no_plot=True)
    table = fit_result.res2d.rename(columns={'Term': 'pathway_id'}).copy()
    table['scope'] = de_rows.scope.iloc[0]
    table['segment_class'] = de_rows.segment_class.iloc[0]
    return table

conventional_gsea_parts = []
for (scope, label), rows in conventional_de.groupby(['scope', 'segment_class'], sort=True):
    result = cached_frame('conventional_gsea_' + scope + '_' + label,
        lambda data=rows: conventional_gsea(data), root=output / 'stage_cache',
        params={'method': 'multilevel', 'weight': 1, 'seed': seed,
                'gseapy': version('gseapy')},
        inputs={'de': rows[['gene', 'stat']], 'sets': conventional_gene_sets},
        code=digest([conventional_code, 'multilevel-gsea-v1']))
    conventional_gsea_parts.append(result)
conventional_gsea_results = pd.concat(conventional_gsea_parts, ignore_index=True)
conventional_gsea_results['NES'] = pd.to_numeric(conventional_gsea_results.NES)
conventional_gsea_results['p_nominal'] = pd.to_numeric(conventional_gsea_results['NOM p-val'])
if conventional_gsea_results.p_nominal.isna().any():
    raise ValueError('GSEA returned missing nominal p-values; inspect the rankings.')
conventional_gsea_results['q_family'] = conventional_gsea_results.groupby('scope').p_nominal.transform(
    lambda values: multipletests(values, method='fdr_bh')[1])
conventional_gsea_results['direction'] = np.where(
    conventional_gsea_results.scope.eq('species_within_segment'),
    np.where(conventional_gsea_results.NES.gt(0), 'human_high', 'mouse_high'),
    np.where(conventional_gsea_results.NES.gt(0), 'segment_enriched', 'rest_enriched'))
conventional_gsea_results.to_csv(output / 'conventional_pt_signed_gsea.csv', index=False)

def conventional_ora(de_rows):
    tested = de_rows.loc[de_rows.stat.notna()].set_index('gene')
    universe = set(tested.index)
    rows = []
    for direction, selected in [('higher', tested.loc[tested.padj.le(.05) & tested.log2FoldChange.ge(1)]),
                                ('lower', tested.loc[tested.padj.le(.05) & tested.log2FoldChange.le(-1)])]:
        selected_genes = set(selected.index)
        for pathway, members in conventional_gene_sets.items():
            member_genes = set(members) & universe
            overlap = member_genes & selected_genes
            rows.append({'scope': de_rows.scope.iloc[0], 'segment_class': de_rows.segment_class.iloc[0],
                         'direction': direction, 'pathway_id': pathway,
                         'n_tested_genes': len(universe), 'n_deg': len(selected_genes),
                         'n_pathway_genes': len(member_genes), 'n_overlap': len(overlap),
                         'p_value': hypergeom.sf(len(overlap) - 1, len(universe),
                                                len(member_genes), len(selected_genes)) if member_genes else 1.,
                         'overlap_genes': ';'.join(sorted(overlap))})
    return pd.DataFrame(rows)

conventional_ora_results = pd.concat([conventional_ora(rows) for _, rows in
    conventional_de.groupby(['scope', 'segment_class'], sort=True)], ignore_index=True)
conventional_ora_results['q_family'] = conventional_ora_results.groupby('scope').p_value.transform(
    lambda values: multipletests(values, method='fdr_bh')[1])
conventional_ora_results.to_csv(output / 'conventional_pt_ora.csv', index=False)
pathway_summary = (conventional_gsea_results.assign(hit=lambda frame: frame.q_family.le(.05))
    .groupby(['scope', 'segment_class', 'direction']).hit.agg(['sum', 'count']))
pathway_summary.to_csv(output / 'conventional_pt_pathway_summary.csv')
display(pathway_summary)

# %% [markdown]
# ### What does continuous position add beyond S1/S2/S3?
#
# Retest the continuous spatial statistic on the **same genes and pathway members** used by the signed discrete benchmark; leave the notebook's original primary list unchanged. A robust continuous candidate passes the original and shared-universe matched screens. Compare it with directional segment GSEA after correcting across the full three-segment family. Also apply the same matched rank-AUC test to the **absolute DESeq2 Wald statistics** within S1/S2/S3, using the same shared gene universe and matching strata as the continuous reanalysis. Correct its empirical p-values across all pathways and all three segment tests for the comparison, while retaining the per-segment q-values as diagnostics. This magnitude companion prevents mixed human-high and mouse-high members from being mistaken for a gain of continuous position merely because signed GSEA cancels them. A continuous-only *screen result* requires that all three GSEA contrasts were testable and that neither discrete screen passed its adjusted cutoff. This does **not** prove a segment effect is zero or establish a significant difference between methods: the models and pathway nulls differ, and one human donor cannot support a population novelty claim. Export the full table, screen-only shortlist, gene-level contrasts, and heatmaps for biological review.

# %%
shared_strata = pd.Series(strata, index=pd.Index(genes)).loc[conventional_genes].to_numpy()
absolute_segment_scores = pd.DataFrame(index=conventional_genes)
for label in segments:
    rows = conventional_de.loc[(conventional_de.scope.eq('species_within_segment'))
                               & (conventional_de.segment_class.eq(label))].set_index('gene')
    absolute_segment_scores['abs_wald_' + label] = rows.loc[conventional_genes, 'stat'].abs()
if not np.isfinite(absolute_segment_scores.to_numpy()).all():
    raise ValueError('Magnitude benchmark contains untestable Wald statistics.')
absolute_segment_tests = cached_frame('conventional_abs_wald_auc',
    lambda: matched_pathway_tests(absolute_segment_scores, conventional_gene_sets,
                                  shared_strata, n_null=n_null, seed=seed),
    root=output / 'stage_cache', params={'n_null': n_null, 'seed': seed,
                                        'ranking': 'absolute_deseq2_wald'},
    inputs={'scores': absolute_segment_scores, 'sets': conventional_gene_sets,
            'strata': shared_strata},
    code=digest([conventional_code, Path(remodeling.__file__), 'abs-wald-auc-v1']))
absolute_segment_tests['q_three_segments'] = multipletests(
    absolute_segment_tests.p_empirical, method='fdr_bh')[1]
absolute_segment_tests.to_csv(output / 'conventional_pt_abs_deseq2_rank_auc.csv', index=False)
shared_scores = gene_stats.loc[conventional_genes, ['T_spatial']]
shared_spatial = cached_frame('conventional_shared_spatial_auc',
    lambda: matched_pathway_tests(shared_scores, conventional_gene_sets, shared_strata,
                                  n_null=n_null, seed=seed),
    root=output / 'stage_cache', params={'n_null': n_null, 'seed': seed,
                                        'universe': 'shared_de_genes'},
    inputs={'scores': shared_scores, 'sets': conventional_gene_sets, 'strata': shared_strata},
    code=digest([cache_code, Path(remodeling.__file__), 'shared-conventional-v1']))
shared_spatial.to_csv(output / 'conventional_shared_universe_spatial_auc.csv', index=False)

species_gsea = conventional_gsea_results.loc[
    conventional_gsea_results.scope.eq('species_within_segment')].copy()
benchmark = shared_spatial[['pathway_id', 'effect', 'q_empirical']].rename(columns={
    'effect': 'spatial_effect_shared', 'q_empirical': 'spatial_q_shared'})
benchmark = benchmark.merge(spatial[['pathway_id', 'effect', 'q_empirical']].rename(columns={
    'effect': 'spatial_effect_original', 'q_empirical': 'spatial_q_original'}),
    on='pathway_id', validate='one_to_one')
for label in segments:
    values = species_gsea.loc[species_gsea.segment_class.eq(label),
        ['pathway_id', 'NES', 'p_nominal', 'q_family']]
    benchmark = benchmark.merge(values.rename(columns={
        'NES': f'NES_{label}', 'p_nominal': f'gsea_p_{label}',
        'q_family': f'gsea_q_{label}'}), on='pathway_id', how='left', validate='one_to_one')
for label in segments:
    magnitude = absolute_segment_tests.loc[absolute_segment_tests.statistic.eq('abs_wald_' + label),
        ['pathway_id', 'effect', 'q_three_segments']]
    benchmark = benchmark.merge(magnitude.rename(columns={
        'effect': f'magnitude_effect_{label}', 'q_three_segments': f'magnitude_q_{label}'}),
        on='pathway_id', validate='one_to_one')
benchmark['all_segments_tested'] = benchmark[[f'NES_{label}' for label in segments]].notna().all(axis=1)
benchmark['continuous_candidate'] = (benchmark.spatial_effect_original.gt(0)
    & benchmark.spatial_q_original.le(.05) & benchmark.spatial_effect_shared.gt(0)
    & benchmark.spatial_q_shared.le(.05))
benchmark['any_segment_gsea_hit'] = benchmark[[f'gsea_q_{label}' for label in segments]].le(.05).any(axis=1)
benchmark['any_segment_magnitude_hit'] = np.column_stack([
    benchmark[f'magnitude_effect_{label}'].gt(0) & benchmark[f'magnitude_q_{label}'].le(.05)
    for label in segments]).any(axis=1)
benchmark['any_segment_signal'] = (benchmark.any_segment_gsea_hit
    | benchmark.any_segment_magnitude_hit)
benchmark['best_segment_q'] = benchmark[[f'gsea_q_{label}' for label in segments]].min(axis=1)
benchmark['best_magnitude_q'] = benchmark[[f'magnitude_q_{label}' for label in segments]].min(axis=1)
benchmark['max_abs_segment_NES'] = benchmark[[f'NES_{label}' for label in segments]].abs().max(axis=1)
benchmark['comparison_group'] = np.select(
    [~benchmark.all_segments_tested,
     benchmark.continuous_candidate & ~benchmark.any_segment_signal,
     benchmark.continuous_candidate & benchmark.any_segment_signal,
     ~benchmark.continuous_candidate & benchmark.any_segment_signal],
    ['not comparable: discrete coverage', 'continuous screen only',
     'continuous and discrete', 'discrete screen only'], default='neither screen')
benchmark = benchmark.merge(phenotypes[['pathway_id', 'peak_position', 'affected_width']],
                            on='pathway_id', how='left', validate='one_to_one')
assert benchmark.pathway_id.is_unique
assert benchmark.loc[benchmark.comparison_group.eq('continuous screen only'),
                     'all_segments_tested'].all()
benchmark['max_segment_magnitude_effect'] = benchmark[[
    f'magnitude_effect_{label}' for label in segments]].max(axis=1)
benchmark.to_csv(output / 'conventional_pt_vs_continuous_pathways.csv', index=False)
continuous_screen_only = benchmark.loc[benchmark.comparison_group.eq('continuous screen only')].sort_values(
    ['spatial_q_shared', 'spatial_effect_shared'], ascending=[True, False])
continuous_screen_only.to_csv(output / 'conventional_pt_continuous_screen_only.csv', index=False)
benchmark_counts = benchmark.comparison_group.value_counts().rename('pathways')
benchmark_counts.to_csv(output / 'conventional_pt_vs_continuous_summary.csv')
display(benchmark_counts)
colors = {'neither screen': '0.75', 'continuous screen only': '#D55E00',
          'continuous and discrete': '#0072B2', 'discrete screen only': '#009E73',
          'not comparable: discrete coverage': '#555555'}
fig, ax = plt.subplots(figsize=(7, 6), layout='constrained')
for group in colors:
    rows = benchmark.loc[benchmark.comparison_group.eq(group)]
    if rows.empty:
        continue
    ax.scatter(rows.max_segment_magnitude_effect, rows.spatial_effect_shared,
               s=15 if group == 'neither screen' else 24,
               alpha=.45 if group == 'neither screen' else .8,
               color=colors[group], label=f'{group} (n={len(rows)})')
ax.axhline(0, color='0.5', lw=.7)
ax.axvline(0, color='0.5', lw=.7)
ax.set(xlabel='Strongest discrete pathway magnitude effect (rank-AUC − 0.5)',
       ylabel='Continuous spatial pathway effect (rank-AUC − 0.5)',
       title='Continuous versus S1/S2/S3 magnitude screens')
ax.legend(frameon=False, fontsize=8, loc='best')
fig.savefig(output / 'conventional_pt_vs_continuous_effects.pdf', bbox_inches='tight')
plt.show()
plt.close(fig)
display(continuous_screen_only[['pathway_id', 'spatial_effect_shared', 'spatial_q_shared',
                                'best_segment_q', 'best_magnitude_q', 'max_abs_segment_NES',
                                'peak_position', 'affected_width']].head(20).round(3))

# %%
from textwrap import shorten

def pathway_segment_heatmap(scope, title, filename, top_per_segment=7):
    rows = conventional_gsea_results.loc[conventional_gsea_results.scope.eq(scope)]
    chosen = list(dict.fromkeys(rows.sort_values(['q_family', 'p_nominal'])
        .groupby('segment_class', sort=False).head(top_per_segment).pathway_id))
    if scope == 'species_within_segment':
        chosen = list(dict.fromkeys(chosen + continuous_screen_only.pathway_id.head(8).tolist()))
    if not chosen:
        return
    nes = rows.pivot(index='pathway_id', columns='segment_class', values='NES').reindex(
        index=chosen, columns=segments)
    q = rows.pivot(index='pathway_id', columns='segment_class', values='q_family').reindex(
        index=chosen, columns=segments)
    limit = max(2., float(np.nanquantile(abs(nes.to_numpy(float)), .98)))
    fig, ax = plt.subplots(figsize=(8, max(5, .30 * len(chosen) + 1.5)), layout='constrained')
    image = ax.imshow(nes.to_numpy(float), cmap='RdBu_r', vmin=-limit, vmax=limit,
                      aspect='auto', interpolation='none')
    ax.set(xticks=range(len(segments)), xticklabels=segments,
           yticks=range(len(chosen)),
           yticklabels=[shorten(item.replace('Reactome_2022::', 'Reactome: ')
                                .replace('MSigDB_Hallmark_2020::', 'Hallmark: ')
                                .replace('KEGG_2019_Mouse::', 'KEGG: '),
                                width=70, placeholder='…') for item in chosen],
           title=title, xlabel='Dots: pooled three-segment BH q ≤ 0.05; blank: not testable')
    for i in range(len(chosen)):
        for j in range(len(segments)):
            if pd.notna(q.iat[i, j]) and q.iat[i, j] <= .05:
                color = 'white' if abs(nes.iat[i, j]) > limit / 2 else 'black'
                ax.text(j, i, '•', ha='center', va='center', color=color, fontsize=12)
    fig.colorbar(image, ax=ax, label='Signed GSEA NES')
    fig.savefig(output / filename, bbox_inches='tight')
    plt.show()
    plt.close(fig)

pathway_segment_heatmap('species_within_segment',
    'Human-high (red) versus mouse-high (blue) pathways by PT segment',
    'conventional_species_pathway_heatmap.pdf')
pathway_segment_heatmap('segment_vs_rest',
    'Segment-enriched (red) versus other-PT-enriched (blue) pathways',
    'conventional_segment_marker_pathway_heatmap.pdf')

# %% [markdown]
# ## 9b · Match the pathway question and enrichment test
#
# The earlier “discrete screen only” label compares **signed segment GSEA for any species difference** with an **unsigned test of position-dependent shape**. Those are different hypotheses. Here we ask the same directional question in each reviewed segment on both sides: which pathways are human-high or mouse-high? We use the same count-supported ortholog universe, pathway membership, weighted multilevel GSEA, and BH family spanning all pathway × segment tests.
#
# For the smooth model, evaluate its full-model human–mouse working $Z_g(s)$ curve over the **observed positions of each specimen in that segment**. Average positions within each specimen, then average the four specimen distributions equally. The resulting gene ranking is an average local working $Z$, not a calibrated donor-level test statistic. Segment labels define where the smooth curve is queried in this *benchmark*; this is not a cluster-independent discovery test. The fitted smooth curves themselves are unchanged.
#
# The count-based DESeq2 Wald statistic and the structure-level smooth-model $Z$ still have different noise models. Neither screen is guaranteed to contain all hits from the other. This paired GSEA comparison isolates much of the earlier pathway-test mismatch; the original $T_{\rm spatial}$ screen continues to answer the distinct question of changing trajectory shape.

# %%
# Compare signed pathway enrichment after holding gene/pathway universe and GSEA method fixed.
# Histogram weights approximate each segment's observed coordinate distribution on the fit grid.
grid_edges = np.r_[-np.inf, (grid[:-1] + grid[1:]) / 2, np.inf]
shared_gene_index = pd.Index(genes).get_indexer(conventional_genes)
assert (shared_gene_index >= 0).all() and fit['z'].shape == (len(genes), len(grid))
continuous_local_gsea_parts = []
continuous_local_weights = {}
for label in segments:
    specimen_weights = []
    for name in expected['mouse'] + expected['human']:
        locations = adata.obs.loc[
            adata.obs['sample'].astype(str).eq(name)
            & adata.obs.segment_class.astype(str).eq(label),
            'shared_pseudospace'].to_numpy(float)
        if len(locations) < 15:
            raise ValueError(f'{name} {label}: too few structures for a segment-matched curve query.')
        histogram = np.histogram(locations, bins=grid_edges)[0].astype(float)
        specimen_weights.append(histogram / histogram.sum())
    weights = np.mean(specimen_weights, axis=0)
    assert np.isclose(weights.sum(), 1) and len(specimen_weights) == 4
    continuous_local_weights[label] = weights
    scores = fit['z'][shared_gene_index] @ weights
    if not np.isfinite(scores).all():
        raise ValueError(f'{label}: nonfinite smooth-model working Z ranking.')
    ranking = pd.DataFrame({'gene': conventional_genes, 'stat': scores,
                            'scope': 'smooth_species_within_segment',
                            'segment_class': label})
    result = cached_frame('smooth_signed_gsea_' + label,
        lambda data=ranking: conventional_gsea(data), root=output / 'stage_cache',
        params={'method': 'multilevel', 'weight': 1, 'seed': seed,
                'gseapy': version('gseapy'), 'coordinate_query': 'equal_specimen_segment_positions'},
        inputs={'ranking': ranking[['gene', 'stat']], 'sets': conventional_gene_sets},
        code=digest([conventional_code, Path(remodeling.__file__), 'smooth-segment-gsea-v1']))
    continuous_local_gsea_parts.append(result)
pd.DataFrame(continuous_local_weights, index=grid).to_csv(
    output / "continuous_segment_coordinate_weights.csv", index_label="position")
continuous_local_gsea = pd.concat(continuous_local_gsea_parts, ignore_index=True)
continuous_local_gsea['NES'] = pd.to_numeric(continuous_local_gsea.NES)
continuous_local_gsea['p_nominal'] = pd.to_numeric(continuous_local_gsea['NOM p-val'])
if continuous_local_gsea.p_nominal.isna().any():
    raise ValueError('Smooth segment GSEA returned missing nominal p-values.')
continuous_local_gsea['q_family'] = multipletests(
    continuous_local_gsea.p_nominal, method='fdr_bh')[1]
continuous_local_gsea['direction'] = np.where(
    continuous_local_gsea.NES.gt(0), 'human_high', 'mouse_high')
continuous_local_gsea.to_csv(output / 'continuous_segment_signed_gsea.csv', index=False)

discrete_signed_gsea = conventional_gsea_results.loc[
    conventional_gsea_results.scope.eq('species_within_segment')].copy()
for rows in (discrete_signed_gsea, continuous_local_gsea):
    assert rows.groupby('segment_class').pathway_id.nunique().reindex(segments).eq(len(benchmark)).all()
    assert set(rows.pathway_id) == set(benchmark.pathway_id)
signed_comparison = benchmark[['pathway_id', 'comparison_group', 'spatial_q_original',
                               'spatial_q_shared']].copy()
for prefix, rows in [('cluster', discrete_signed_gsea), ('smooth', continuous_local_gsea)]:
    by_pathway = rows.groupby('pathway_id').agg(
        **{f'{prefix}_best_q': ('q_family', 'min'),
           f'{prefix}_best_abs_NES': ('NES', lambda x: x.abs().max())})
    signed_comparison = signed_comparison.merge(by_pathway, on='pathway_id', validate='one_to_one')
    signed_comparison[f'{prefix}_signed_hit'] = signed_comparison[f'{prefix}_best_q'].le(.05)
signed_comparison['signed_gsea_group'] = np.select(
    [signed_comparison.cluster_signed_hit & signed_comparison.smooth_signed_hit,
     signed_comparison.cluster_signed_hit,
     signed_comparison.smooth_signed_hit],
    ['both signed GSEA screens', 'cluster signed GSEA only', 'smooth signed GSEA only'],
    default='neither signed GSEA screen')
signed_comparison.to_csv(output / 'conventional_vs_smooth_signed_gsea_pathways.csv', index=False)
signed_counts = signed_comparison.signed_gsea_group.value_counts().rename('pathways')
signed_counts.to_csv(output / 'conventional_vs_smooth_signed_gsea_summary.csv')
display(signed_counts)
prior_discrete = signed_comparison.comparison_group.eq('discrete screen only')
print('Prior spatial-shape-comparison “discrete screen only” pathways also detected by '
      'smooth segment GSEA:', int((prior_discrete & signed_comparison.smooth_signed_hit).sum()),
      '/', int(prior_discrete.sum()))


# %% [markdown]
# The overlap table is the **like-for-like signed pathway screen**. A pathway in “cluster signed GSEA only” is a screen discordance, not proof that the smooth model cannot represent its biology. Review its NES and q-values in all three segments, individual genes, and specimen-level effects. The shape test remains useful for asking what changes *within* and *between* segments, beyond broad level differences. A strict superset test would fit discrete segment effects and added smooth positional effects in one nested expression model; these two existing models are not nested.

# %%
# Show overlap and the directional signal behind it, with one color scale for both methods.
overlap_signed = pd.crosstab(
    signed_comparison.cluster_signed_hit, signed_comparison.smooth_signed_hit).reindex(
    index=[False, True], columns=[False, True], fill_value=0)
fig, ax = plt.subplots(figsize=(5, 4), layout='constrained')
ax.imshow(overlap_signed.to_numpy(), cmap='Blues', aspect='auto')
ax.set(xticks=[0, 1], xticklabels=['No', 'Yes'], yticks=[0, 1], yticklabels=['No', 'Yes'],
       xlabel='Smooth segment GSEA q ≤ 0.05', ylabel='Cluster segment GSEA q ≤ 0.05',
       title='Signed pathway screen overlap')
for (i, j), count in np.ndenumerate(overlap_signed.to_numpy()):
    ax.text(j, i, str(count), ha='center', va='center',
            color='white' if count > overlap_signed.to_numpy().max() / 2 else 'black')
fig.savefig(output / 'conventional_vs_smooth_signed_gsea_overlap.pdf', bbox_inches='tight')
plt.show()
plt.close(fig)

prior_ids = set(signed_comparison.loc[
    signed_comparison.comparison_group.eq('discrete screen only'), 'pathway_id'])
prior_top = (discrete_signed_gsea.loc[discrete_signed_gsea.pathway_id.isin(prior_ids)]
    .groupby('pathway_id').q_family.min().nsmallest(8).index.tolist())
cluster_top = signed_comparison.loc[
    signed_comparison.signed_gsea_group.eq('cluster signed GSEA only')].nsmallest(
    5, 'cluster_best_q').pathway_id.tolist()
smooth_top = signed_comparison.loc[
    signed_comparison.signed_gsea_group.eq('smooth signed GSEA only')].nsmallest(
    5, 'smooth_best_q').pathway_id.tolist()
chosen = list(dict.fromkeys(prior_top + cluster_top + smooth_top))
combined_signed = pd.concat([
    discrete_signed_gsea.assign(method='Cluster'),
    continuous_local_gsea.assign(method='Smooth')], ignore_index=True)
columns = pd.MultiIndex.from_product([['Cluster', 'Smooth'], segments],
                                      names=['method', 'segment_class'])
nes = combined_signed.pivot(index='pathway_id', columns=['method', 'segment_class'],
                            values='NES').reindex(index=chosen, columns=columns)
q = combined_signed.pivot(index='pathway_id', columns=['method', 'segment_class'],
                          values='q_family').reindex(index=chosen, columns=columns)
assert np.isfinite(nes.to_numpy(float)).all() and np.isfinite(q.to_numpy(float)).all()
limit = max(2., float(np.nanquantile(abs(nes.to_numpy(float)), .98)))
fig, ax = plt.subplots(figsize=(10, max(5, .35 * len(chosen) + 1.5)), layout='constrained')
image = ax.imshow(nes.to_numpy(float), cmap='RdBu_r', vmin=-limit, vmax=limit,
                  aspect='auto', interpolation='none')
ax.axvline(2.5, color='black', lw=1)
ax.set(xticks=range(len(columns)),
       xticklabels=[f'{method}\n{segment}' for method, segment in columns],
       yticks=range(len(chosen)),
       yticklabels=[shorten(item.replace('Reactome_2022::', 'Reactome: ')
                            .replace('MSigDB_Hallmark_2020::', 'Hallmark: ')
                            .replace('KEGG_2019_Mouse::', 'KEGG: '),
                            width=65, placeholder='…') for item in chosen],
       title='Human-high (red) and mouse-high (blue) signed pathway enrichment')
for i, j in np.ndindex(q.shape):
    if q.iat[i, j] <= .05:
        ax.text(j, i, '•', ha='center', va='center',
                color='white' if abs(nes.iat[i, j]) > limit / 2 else 'black', fontsize=12)
fig.colorbar(image, ax=ax, label='GSEA NES; dot: pooled three-segment BH q ≤ 0.05')
fig.savefig(output / 'conventional_vs_smooth_signed_gsea_heatmap.pdf', bbox_inches='tight')
plt.show()
plt.close(fig)
display(signed_comparison.loc[signed_comparison.signed_gsea_group.eq(
    'cluster signed GSEA only')].nsmallest(12, 'cluster_best_q')[
    ['pathway_id', 'cluster_best_q', 'smooth_best_q', 'spatial_q_shared']])


# %% [markdown]
# ## 9c · Benchmark total human–mouse differences with shared models
#
# This benchmark asks the **same biological question** in each model: do human and mouse differ anywhere among the measured structures? Every method gets one ranking opportunity per gene, and each pathway family is tested once per method. There is no threshold tuning or second pass over a selected pathway list.
#
# The controlled expression comparison uses shared expression and the same within-species, sum-to-zero specimen intercepts throughout. The discrete model adds human-by-segment steps to that common baseline; the smooth model uses the existing position-dependent species model. The augmented model contains both sets of terms, so its increment asks what position adds after segment-level species differences are represented. This controlled WLS benchmark is a descriptive sensitivity analysis, not the standard pseudobulk test. The absolute maximum DESeq2 Wald score is a magnitude companion to signed GSEA, not a replacement for it; selecting the largest of three segment Wald magnitudes is included in the matched competitive null. The count and structure models retain different noise assumptions.
#
# The matched null controls pathway size and the existing expression, detection, and coverage strata. The incremental position result is secondary explanatory support: it does not establish a direction or prove within-segment biology, and it is not a primary discovery filter. With two mouse specimens and one human donor, all results remain descriptive and structure uncertainty is substantial. A benchmark screen need not contain every hit from another model.
#
# **Read the comparisons separately.** The primary controlled screen holds the expression scale, Gaussian noise model, shared spline baseline, specimen adjustment, gene/pathway universe, enrichment test, and correction rule fixed; only the representation of species differences changes (segment steps versus smooth position). It tests competitive enrichment of conditional gene scores, not whether a pathway has a donor-level species effect. The native max-Wald comparison aligns magnitude and the pathway test but keeps DESeq2's different count noise model. The signed-GSEA reference retains its directional question and is reported separately; its recovery fraction cannot be substituted with controlled-model retention.
#
# Primary calls use positive AUC effect and method-specific BH q ≤ 0.05, each over the identical pathway family. Joint BH across both controlled methods is a reported sensitivity analysis. Incremental position support has its own full-pathway correction and can occur even without a total-screen call; do not add that count to the continuous discovery count. No original-universe confirmation is required here because both sides use the same shared universe.
#

# %%
from scipy import sparse
from pseudospace.levelshape import build_ls_designs
import inspect

BENCHMARK_LOGIC_VERSION = '12-total-difference-benchmark-v3'

def benchmark_nested_gene_scores(expression, base, discrete, smooth, weights, block_size=256):
    """Blockwise partial-F scores from explicitly nested weighted least-squares designs."""
    if getattr(expression, 'ndim', None) != 2:
        raise ValueError('Expression must be a two-dimensional matrix.')
    n = expression.shape[0]
    if not sparse.issparse(expression):
        expression = np.asarray(expression, dtype=float)
    matrices = [np.asarray(x, dtype=float) for x in (base, discrete, smooth)]
    base, discrete, smooth = matrices
    weights = np.asarray(weights, dtype=float)
    if any(x.ndim != 2 or x.shape[0] != n for x in matrices):
        raise ValueError('Every design must be a matrix with the expression row count.')
    if any(x.shape[1] == 0 or not np.isfinite(x).all() for x in matrices):
        raise ValueError('Design matrices must have columns and finite values.')
    if weights.shape != (n,) or not np.isfinite(weights).all() or (weights <= 0).any():
        raise ValueError('Weights must be finite, positive, and match expression rows.')
    if expression.shape[1] == 0 or block_size < 1:
        raise ValueError('Expression needs genes and block_size must be positive.')
    if sparse.issparse(expression):
        if not np.isfinite(expression.data).all():
            raise ValueError('Sparse expression contains nonfinite values.')
    elif not np.isfinite(expression).all():
        raise ValueError('Expression contains nonfinite values.')

    def orthobasis(design):
        weighted = np.sqrt(weights)[:, None] * design
        u, singular, _ = np.linalg.svd(weighted, full_matrices=False)
        tol = np.finfo(float).eps * max(weighted.shape) * (singular[0] if len(singular) else 0.)
        rank = int(np.sum(singular > tol))
        return u[:, :rank], rank, tol

    q0, r0, _ = orthobasis(base)
    qd, rd, _ = orthobasis(discrete)
    qs, rs, _ = orthobasis(smooth)
    union = np.column_stack([discrete, smooth])
    qu, ru, _ = orthobasis(union)
    nesting_tol = 1e-8
    for name, q in [('discrete', qd), ('smooth', qs), ('union', qu)]:
        if np.linalg.norm(q0 - q @ (q.T @ q0)) > nesting_tol * max(1., np.sqrt(r0)):
            raise ValueError(f'Base design is not nested in {name} design.')
    for name, q in [('discrete', qd), ('smooth', qs)]:
        if np.linalg.norm(q - qu @ (qu.T @ q)) > nesting_tol * max(1., np.sqrt(q.shape[1])):
            raise ValueError(f'{name} design is not nested in the augmented design.')
    if n <= ru + 2:
        raise ValueError('Augmented design leaves too few residual degrees of freedom.')
    df_d, df_s, df_pos = rd-r0, rs-r0, ru-rd
    if min(df_d, df_s, df_pos) <= 0:
        raise ValueError('Each requested comparison must add positive design degrees of freedom.')
    sse0 = np.zeros(expression.shape[1]); ssed = np.zeros_like(sse0)
    sses = np.zeros_like(sse0); sseu = np.zeros_like(sse0)
    sqrtw = np.sqrt(weights)
    for start in range(0, expression.shape[1], block_size):
        stop = min(start + block_size, expression.shape[1])
        block = expression[:, start:stop]
        block = block.toarray() if sparse.issparse(block) else np.asarray(block)
        if not np.isfinite(block).all():
            raise ValueError('Expression block contains nonfinite values.')
        yw = sqrtw[:, None] * block
        for target, q in [(sse0, q0), (ssed, qd), (sses, qs), (sseu, qu)]:
            resid = yw - q @ (q.T @ yw)
            target[start:stop] = np.einsum('ij,ij->j', resid, resid)
    def partial_f(reduced, full, df_added, rank_full):
        improvement = reduced - full
        tol = 1e-10 * np.maximum(1., reduced)
        if np.any(improvement < -tol):
            raise ValueError('Nested-model SSE materially increased in the larger model.')
        improvement = np.maximum(improvement, 0.)
        if np.any(full <= 0):
            raise ValueError('A full-model residual SSE is nonpositive.')
        return (improvement / df_added) / (full / (n-rank_full))
    return {'T_discrete_total': partial_f(sse0, ssed, df_d, rd),
            'T_continuous_total': partial_f(sse0, sses, df_s, rs),
            'T_position_given_segments': partial_f(ssed, sseu, df_pos, ru),
            'design_df': np.array([n, r0, rd, rs, ru, df_d, df_s, df_pos, n-ru], dtype=np.int64)}

def regression_check_benchmark_nested_gene_scores():
    """Small data-free check against independent weighted least-squares projections."""
    rng = np.random.default_rng(1209)
    n = 84
    position_check = np.linspace(0., 1., n)
    human_check = np.tile([0., 1.], n // 2)
    baseline = np.column_stack([np.ones(n), position_check])
    step = human_check[:, None]
    curve = (human_check * position_check)[:, None]
    base = baseline
    discrete = np.column_stack([base, step])
    smooth = np.column_stack([base, curve])
    values = np.column_stack([
        baseline @ np.array([1.2, -.3]) + .7 * human_check,
        baseline @ np.array([-.2, .9]) + .4 * human_check * position_check,
    ]) + rng.normal(0, .15, size=(n, 2))
    weights = rng.uniform(.5, 1.5, size=n)
    expression = sparse.csr_matrix(values)
    got = benchmark_nested_gene_scores(expression, base, discrete, smooth, weights)
    union = np.column_stack([discrete, smooth])
    def independent_f(reduced, full, added_df, rank_full):
        scores = []
        for j in range(values.shape[1]):
            def sse(design):
                root = np.sqrt(weights)
                beta = np.linalg.lstsq(root[:, None] * design,
                    root * values[:, j], rcond=None)[0]
                residual = root * (values[:, j] - design @ beta)
                return residual @ residual
            scores.append(((sse(reduced) - sse(full)) / added_df)
                          / (sse(full) / (n - rank_full)))
        return np.asarray(scores)
    rank_d = np.linalg.matrix_rank(np.sqrt(weights)[:, None] * discrete)
    rank_s = np.linalg.matrix_rank(np.sqrt(weights)[:, None] * smooth)
    rank_u = np.linalg.matrix_rank(np.sqrt(weights)[:, None] * union)
    np.testing.assert_allclose(got['T_discrete_total'],
        independent_f(base, discrete, rank_d - np.linalg.matrix_rank(np.sqrt(weights)[:, None] * base), rank_d))
    np.testing.assert_allclose(got['T_continuous_total'],
        independent_f(base, smooth, rank_s - np.linalg.matrix_rank(np.sqrt(weights)[:, None] * base), rank_s))
    np.testing.assert_allclose(got['T_position_given_segments'],
        independent_f(discrete, union, rank_u-rank_d, rank_u))
    redundant = benchmark_nested_gene_scores(expression,
        np.column_stack([base, base[:, 0]]), np.column_stack([discrete, base[:, 0]]),
        np.column_stack([smooth, base[:, 0]]), weights)
    np.testing.assert_allclose(redundant['T_discrete_total'], got['T_discrete_total'])
    np.testing.assert_allclose(redundant['T_continuous_total'], got['T_continuous_total'])
    try:
        benchmark_nested_gene_scores(expression, np.column_stack([base, curve]),
                                     discrete, smooth, weights)
    except ValueError:
        pass
    else:
        raise AssertionError('A nonnested reduced design must be rejected.')
    for bad_weights in (weights[:-1], np.r_[weights[:-1], 0.], np.r_[weights[:-1], np.nan]):
        try:
            benchmark_nested_gene_scores(expression, base, discrete, smooth, bad_weights)
        except ValueError:
            pass
        else:
            raise AssertionError('Mismatched or invalid weights must be rejected.')
    return 'projection regression checks passed'

assert regression_check_benchmark_nested_gene_scores() == 'projection regression checks passed'

# The baseline and nuisance terms are rebuilt from the same specimens and weights as the original fits.
model_weights = np.zeros(len(position), dtype=float)
model_nuisance = []
for group in (False, True):
    names = np.unique(specimen[human == group])
    if len(names) != 2:
        raise ValueError('The benchmark expects two specimens within each species.')
    for name in names:
        mask = specimen == name
        model_weights[mask] = len(position) / (2 * len(names) * mask.sum())
    model_nuisance.extend((specimen == name).astype(float) - (specimen == names[-1]).astype(float)
                          for name in names[:-1])
base_design, _, smooth_design = build_ls_designs(
    position, human.astype(float), fit['knots'])[:3]
base_design = np.column_stack([base_design, *model_nuisance])
smooth_design = np.column_stack([smooth_design, *model_nuisance])
segment_levels = ('PT-S1', 'PT-S2', 'PT-S3')
if set(labels) != set(segment_levels):
    raise ValueError('Unexpected reviewed PT segment labels in the benchmark.')
segment_steps = np.column_stack([
    human.astype(float) * (labels == segment).astype(float) for segment in segment_levels])
discrete_design = np.column_stack([base_design, segment_steps])
benchmark_gene_columns = pd.Index(genes).get_indexer(conventional_genes)
if (benchmark_gene_columns < 0).any():
    raise ValueError('Count-supported genes could not be aligned to expression.')
benchmark_Y = Y[:, benchmark_gene_columns]
if benchmark_Y.shape[1] != len(conventional_genes):
    raise ValueError('Count-supported gene columns could not be aligned to expression.')

benchmark_inputs = {'expression': digest(benchmark_Y), 'designs': [base_design, discrete_design,
    smooth_design], 'weights': model_weights, 'genes': conventional_genes.tolist()}
benchmark_code = digest([BENCHMARK_LOGIC_VERSION, inspect.getsource(benchmark_nested_gene_scores)])
benchmark_scores = cached_payload('total_difference_gene_scores',
    lambda: benchmark_nested_gene_scores(benchmark_Y, base_design, discrete_design,
                                         smooth_design, model_weights),
    root=output / 'stage_cache', params={'logic': BENCHMARK_LOGIC_VERSION},
    inputs=benchmark_inputs, code=benchmark_code)
# Validate the cached or freshly calculated result against the original total score on every run.
np.testing.assert_allclose(benchmark_scores['T_continuous_total'],
    gene_stats.loc[conventional_genes, 'T_total'].to_numpy(float), rtol=1e-6, atol=1e-8)
design_audit = pd.Series(benchmark_scores['design_df'], name='value', index=[
    'n', 'base_rank', 'discrete_rank', 'smooth_rank', 'augmented_rank',
    'discrete_added', 'smooth_added', 'position_given_segments_added', 'augmented_residual'])
design_audit.to_csv(output / 'total_difference_design_audit.csv')


# %%
# One ranking per method, with the maximum absolute segment Wald included as its own method.
native_max_wald = absolute_segment_scores.max(axis=1).rename('T_deseq2_any_segment')
score_frame = pd.DataFrame({name: benchmark_scores[name] for name in (
    'T_discrete_total', 'T_continuous_total', 'T_position_given_segments')}, index=conventional_genes)
score_frame.index.name = 'gene'
all_benchmark_rankings = score_frame.join(native_max_wald, validate='one_to_one')
all_benchmark_rankings.to_csv(output / 'total_difference_gene_rankings.csv')

benchmark_pathway_tests = cached_frame('total_difference_matched_pathways',
    lambda: matched_pathway_tests(all_benchmark_rankings, conventional_gene_sets, shared_strata,
                                  n_null=n_null, seed=seed),
    root=output / 'stage_cache', params={'n_null': n_null, 'seed': seed,
        'methods': list(all_benchmark_rankings.columns)},
    inputs={'rankings': all_benchmark_rankings, 'sets': conventional_gene_sets,
            'strata': shared_strata}, code=digest([BENCHMARK_LOGIC_VERSION,
                Path(remodeling.__file__), 'matched-pathways-v1']))
benchmark_pathway_tests['q_method'] = benchmark_pathway_tests.groupby(
    'statistic').p_empirical.transform(lambda values: multipletests(values, method='fdr_bh')[1])
controlled_mask = benchmark_pathway_tests.statistic.isin(
    ['T_discrete_total', 'T_continuous_total'])
benchmark_pathway_tests['q_controlled_pooled'] = np.nan
benchmark_pathway_tests.loc[controlled_mask, 'q_controlled_pooled'] = multipletests(
    benchmark_pathway_tests.loc[controlled_mask, 'p_empirical'], method='fdr_bh')[1]
benchmark_pathway_tests.to_csv(output / 'total_difference_all_pathway_tests.csv', index=False)
benchmark_manifest = {
    'logic_version': BENCHMARK_LOGIC_VERSION,
    'n_null': int(n_null), 'seed': int(seed), 'alpha': .05,
    'pathway_family': 'all common-count-supported pathways tested once per ranking method',
    'method_q_family': 'BH across all tested pathways separately within each method',
    'controlled_pooled_q_family': 'BH across discrete-total and continuous-total tests jointly',
    'incremental_q_family': 'BH across all pathways for position-given-segments; secondary support',
    'ranking_methods': list(all_benchmark_rankings.columns),
    'native_magnitude_method': 'maximum absolute DESeq2 Wald statistic across PT-S1/S2/S3 per gene',
    'expression_shape': [int(benchmark_Y.shape[0]), int(benchmark_Y.shape[1])],
    'n_genes': int(len(conventional_genes)), 'n_pathways': int(len(conventional_gene_sets)),
    'design_df': {key: int(value) for key, value in design_audit.items()},
    'gene_score_digest': digest(all_benchmark_rankings),
    'matched_test_cache_key': digest([BENCHMARK_LOGIC_VERSION, n_null, seed,
        all_benchmark_rankings, conventional_gene_sets, shared_strata]),
    'gene_score_cache_code': benchmark_code,
}
with (output / 'total_difference_benchmark_manifest.json').open('w') as handle:
    json.dump(benchmark_manifest, handle, indent=2, sort_keys=True)

# Wide pathway comparison preserves effect and both method-specific and pooled controlled q-values.
benchmark_comparison = benchmark_pathway_tests.pivot(index='pathway_id', columns='statistic',
    values=['effect', 'p_empirical', 'q_method', 'q_controlled_pooled'])
benchmark_comparison.columns = ['_'.join(map(str, column)).rstrip('_')
                                for column in benchmark_comparison.columns]
benchmark_comparison = benchmark_comparison.reset_index()
native_hits = signed_comparison[['pathway_id', 'cluster_signed_hit', 'cluster_best_q']].rename(
    columns={'cluster_signed_hit': 'native_signed_gsea_hit',
             'cluster_best_q': 'native_signed_gsea_best_q'})
benchmark_comparison = benchmark_comparison.merge(native_hits, on='pathway_id',
                                                    validate='one_to_one')
benchmark_comparison['controlled_discrete_hit'] = (
    benchmark_comparison.effect_T_discrete_total.gt(0)
    & benchmark_comparison.q_method_T_discrete_total.le(.05))
benchmark_comparison['controlled_continuous_hit'] = (
    benchmark_comparison.effect_T_continuous_total.gt(0)
    & benchmark_comparison.q_method_T_continuous_total.le(.05))
benchmark_comparison['controlled_discrete_hit_pooled'] = (
    benchmark_comparison.effect_T_discrete_total.gt(0)
    & benchmark_comparison.q_controlled_pooled_T_discrete_total.le(.05))
benchmark_comparison['controlled_continuous_hit_pooled'] = (
    benchmark_comparison.effect_T_continuous_total.gt(0)
    & benchmark_comparison.q_controlled_pooled_T_continuous_total.le(.05))
benchmark_comparison['incremental_position_support'] = (
    benchmark_comparison.effect_T_position_given_segments.gt(0)
    & benchmark_comparison.q_method_T_position_given_segments.le(.05))
benchmark_comparison['native_magnitude_hit'] = (
    benchmark_comparison.effect_T_deseq2_any_segment.gt(0)
    & benchmark_comparison.q_method_T_deseq2_any_segment.le(.05))
benchmark_comparison['controlled_group'] = np.select(
    [benchmark_comparison.controlled_discrete_hit & benchmark_comparison.controlled_continuous_hit,
     benchmark_comparison.controlled_discrete_hit, benchmark_comparison.controlled_continuous_hit],
    ['both', 'discrete only', 'continuous only'], default='neither')
benchmark_comparison['native_magnitude_group'] = np.select(
    [benchmark_comparison.native_magnitude_hit & benchmark_comparison.controlled_continuous_hit,
     benchmark_comparison.native_magnitude_hit, benchmark_comparison.controlled_continuous_hit],
    ['both', 'native magnitude only', 'continuous only'], default='neither')
benchmark_comparison.to_csv(output / 'total_difference_pathway_comparison.csv', index=False)

controlled_discrete = benchmark_comparison.controlled_discrete_hit
controlled_continuous = benchmark_comparison.controlled_continuous_hit
native_magnitude = benchmark_comparison.native_magnitude_hit
native_signed = benchmark_comparison.native_signed_gsea_hit.fillna(False)
controlled_overlap = pd.crosstab(controlled_discrete, controlled_continuous).reindex(
    index=[False, True], columns=[False, True], fill_value=0)
native_overlap = pd.crosstab(native_magnitude, benchmark_comparison.controlled_continuous_hit).reindex(
    index=[False, True], columns=[False, True], fill_value=0)
summary_rows = []
for name, baseline, added in [
    ('controlled_discrete_to_continuous', controlled_discrete, controlled_continuous),
    ('native_signed_to_controlled_continuous', native_signed, controlled_continuous),
    ('controlled_pooled_BH', benchmark_comparison.controlled_discrete_hit_pooled,
     benchmark_comparison.controlled_continuous_hit_pooled)]:
    nbase = int(baseline.sum())
    retained = int((baseline & added).sum())
    summary_rows.append({'comparison': name, 'baseline_hits': nbase,
        'continuous_hits': int(added.sum()), 'overlap': retained,
        'baseline_retention_fraction': retained / nbase if nbase else np.nan,
        'continuous_added': int((added & ~baseline).sum()),
        'incremental_position_support': int(benchmark_comparison.incremental_position_support.sum())})
summary_rows.append({'comparison': 'native_magnitude_vs_controlled_continuous',
    'baseline_hits': int(native_magnitude.sum()), 'continuous_hits': int(controlled_continuous.sum()),
    'overlap': int((native_magnitude & controlled_continuous).sum()),
    'baseline_retention_fraction': ((native_magnitude & controlled_continuous).sum()
        / native_magnitude.sum() if native_magnitude.sum() else np.nan),
    'continuous_added': int((controlled_continuous & ~native_magnitude).sum()),
    'incremental_position_support': int(benchmark_comparison.incremental_position_support.sum())})
incremental_by_category = (benchmark_comparison.groupby('controlled_group')
    .incremental_position_support.sum().rename('incremental_position_support').reset_index())
benchmark_summary = pd.DataFrame(summary_rows)
benchmark_summary.to_csv(output / 'total_difference_summary.csv', index=False)
incremental_by_category.to_csv(output / 'total_difference_incremental_by_category.csv', index=False)
controlled_overlap.to_csv(output / 'total_difference_controlled_overlap.csv')
native_overlap.to_csv(output / 'total_difference_native_magnitude_overlap.csv')


# %%
# Two overlap views keep controlled and native count-model comparisons distinct.
fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), layout='constrained')
for ax, table, xlabel, ylabel, title in [
    (axes[0], controlled_overlap, 'Controlled continuous hit', 'Controlled discrete hit',
     'Controlled model overlap'),
    (axes[1], native_overlap, 'Controlled continuous hit', 'Native max-Wald hit',
     'Native magnitude overlap')]:
    values = table.to_numpy()
    ax.imshow(values, cmap='Blues', aspect='equal')
    ax.set(xticks=[0, 1], xticklabels=['No', 'Yes'], yticks=[0, 1], yticklabels=['No', 'Yes'],
           xlabel=xlabel, ylabel=ylabel, title=title)
    for (i, j), count in np.ndenumerate(values):
        ax.text(j, i, str(count), ha='center', va='center',
                color='white' if count > values.max() / 2 else 'black')
fig.savefig(output / 'total_difference_overlap_heatmaps.pdf', bbox_inches='tight')
plt.show()
plt.close(fig)

benchmark_colors = {'neither': '0.75', 'discrete only': '#009E73',
                    'continuous only': '#D55E00', 'both': '#0072B2'}
fig, ax = plt.subplots(figsize=(7, 5.5), layout='constrained')
for group, color in benchmark_colors.items():
    rows = benchmark_comparison.loc[benchmark_comparison.controlled_group.eq(group)]
    ax.scatter(rows.effect_T_discrete_total, rows.effect_T_continuous_total,
        s=22 if group != 'neither' else 14, alpha=.8 if group != 'neither' else .45,
        color=color, label=f'{group} (n={len(rows)})')
ax.axhline(0, color='0.5', lw=.7); ax.axvline(0, color='0.5', lw=.7)
ax.set(xlabel='Controlled discrete total-difference effect (AUC − 0.5)',
       ylabel='Controlled continuous total-difference effect (AUC − 0.5)',
       title='Matched pathway effects: discrete and continuous models')
ax.legend(frameon=False, fontsize=8)
fig.savefig(output / 'total_difference_controlled_effects.pdf', bbox_inches='tight')
plt.show()
plt.close(fig)

# Concise discrepancy table: strongest q/effect evidence from each total-difference ranking.
show_cols = ['pathway_id', 'effect_T_discrete_total', 'q_method_T_discrete_total',
    'effect_T_continuous_total', 'q_method_T_continuous_total',
    'effect_T_position_given_segments', 'q_method_T_position_given_segments',
    'effect_T_deseq2_any_segment', 'q_method_T_deseq2_any_segment',
    'native_signed_gsea_hit', 'native_signed_gsea_best_q']
top_discrete = benchmark_comparison.loc[benchmark_comparison.controlled_group.eq('discrete only')].sort_values(
    ['q_method_T_discrete_total', 'effect_T_discrete_total', 'pathway_id'],
    ascending=[True, False, True]).head(10)
top_continuous = benchmark_comparison.loc[benchmark_comparison.controlled_group.eq('continuous only')].sort_values(
    ['q_method_T_continuous_total', 'effect_T_continuous_total', 'pathway_id'],
    ascending=[True, False, True]).head(10)
top_discrepancy = pd.concat([top_discrete.assign(top_method='controlled discrete only'),
    top_continuous.assign(top_method='controlled continuous only')], ignore_index=True)
top_discrepancy = top_discrepancy[['top_method', *show_cols]]
top_discrepancy.to_csv(output / 'total_difference_top_discrepancies.csv', index=False)
for group in ('both', 'continuous only', 'discrete only'):
    benchmark_comparison.loc[benchmark_comparison.controlled_group.eq(group)].to_csv(
        output / ('total_difference_' + group.replace(' ', '_') + '.csv'), index=False)
display(benchmark_summary)
display(top_discrepancy.round(4))

# Report retention without substituting one comparator for another.
from IPython.display import Markdown
controlled_report = benchmark_summary.set_index('comparison').loc['controlled_discrete_to_continuous']
native_report = benchmark_summary.set_index('comparison').loc['native_magnitude_vs_controlled_continuous']
signed_report = benchmark_summary.set_index('comparison').loc['native_signed_to_controlled_continuous']
display(Markdown(
    f"**Controlled detection:** continuous retains {int(controlled_report.overlap)} of "
    f"{int(controlled_report.baseline_hits)} step-model hits and adds "
    f"{int(controlled_report.continuous_added)} pathways. "
    f"**Native magnitude companion:** it retains {int(native_report.overlap)} of "
    f"{int(native_report.baseline_hits)} max-Wald hits and adds "
    f"{int(native_report.continuous_added)} pathways. "
    f"**Conventional signed GSEA reference:** it detects {int(signed_report.overlap)} of "
    f"{int(signed_report.baseline_hits)} directional hits. These are different comparisons; "
    "the controlled overlap does not establish complete recovery of conventional signed GSEA. "
    "Additional screen calls need member-gene and specimen-level review before claiming "
    "that within-segment or boundary-spanning patterns caused their detection."))


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
from textwrap import fill

def plot_pathway_evidence(pathway, detail):
    """Show existing pathway curves, prioritized signed gene evidence, and fitted examples."""
    selected = (detail.sort_values(
        ['spatial_driver', 'leading_edge', 'broad_supporter', 'peak_abs_z', 'gene'],
        ascending=[False, False, False, False, True]).head(15).copy())
    indices = gene_index.get_indexer(selected.gene)
    selected = selected.iloc[np.argsort(np.argmax(np.abs(fit['z'][indices]), axis=1), kind='stable')]
    indices = gene_index.get_indexer(selected.gene)
    trace = local_curves[local_curves.pathway_id.eq(pathway)].sort_values('position')
    fig = plt.figure(figsize=(11, 8), layout='constrained')
    layout = fig.add_gridspec(3, 2, height_ratios=[1, 2.4, 1.2])
    for column, label, axis in [('divergence', 'D(s): unsigned excess AUC', fig.add_subplot(layout[0, 0])),
                                ('direction', 'S(s): relative direction', fig.add_subplot(layout[0, 1]))]:
        axis.plot(trace.position, trace[column], color='#0072B2')
        axis.axhline(0, color='0.6', lw=.7)
        axis.set(xlim=(grid[0], grid[-1]), xlabel='Shared PT pseudospace', ylabel=label)
    axis = fig.add_subplot(layout[1, :])
    step = (grid[-1] - grid[0]) / (len(grid) - 1)
    heat = axis.imshow(fit['z'][indices], cmap='RdBu_r', vmin=-6, vmax=6,
        aspect='auto', interpolation='none',
        extent=(grid[0] - step / 2, grid[-1] + step / 2, len(selected) - .5, -.5))
    labels = []
    for row in selected.itertuples():
        roles = ''.join(letter for field, letter in [('spatial_driver', 'D'),
            ('leading_edge', 'E'), ('broad_supporter', 'B')] if getattr(row, field))
        bulk_sign = '+' if row.bulk_logcpm_effect > 0 else '−' if row.bulk_logcpm_effect < 0 else '0'
        labels.append(f'{row.gene}  [{roles or "–"}]  det {row.detection:.0%}  bulk {bulk_sign}')
    axis.set(yticks=np.arange(len(selected)), yticklabels=labels,
             title='Top gene evidence, ordered by peak |Z| position')
    fig.colorbar(heat, ax=axis, label='Signed Z (human − mouse)', extend='both')
    drivers = detail.sort_values(['peak_abs_z', 'gene'], ascending=[False, True]).head(2)
    for j, row in enumerate(drivers.itertuples()):
        axis = fig.add_subplot(layout[2, j])
        index = gene_index.get_loc(row.gene)
        axis.plot(grid, fit['human'][index], color='#D55E00', label='Human')
        axis.plot(grid, fit['mouse'][index], color='#0072B2', ls='--', label='Mouse')
        axis.set(xlim=(grid[0], grid[-1]), title=f'{row.gene} · detection {row.detection:.0%}',
                 xlabel='Shared PT pseudospace', ylabel='Fitted log-normalized expression')
        axis.legend(frameon=False, fontsize=8)
    fig.suptitle(fill(pathway.replace('::', ': ', 1), width=65))
    return fig

for number, pathway in enumerate(illustration_ids[:3], 1):
    detail = member_evidence[member_evidence.pathway_id.eq(pathway)]
    fig = plot_pathway_evidence(pathway, detail)
    fig.savefig(output / f'pathway_gene_evidence_{number:02d}.pdf', bbox_inches='tight')
    plt.show()
    plt.close(fig)


# %% [markdown]
# The heatmap shows **signed** gene evidence: red is human-high and blue is mouse-high. D/E/B mark the existing driver, leading-edge, and broad-support roles; detection and bulk sign give context for each row. Inspect the fitted curves below before treating an extreme Z as an interpretable expression change. Color saturation is for display, not a local significance call.
#

# %% [markdown]
# ## 11 · Collapse redundant terms after discovery
#
# Cluster only the frozen discovery candidates using equal-weight similarity from member overlap, leading-edge overlap, D-curve correlation, and S-curve correlation. Complete linkage avoids merging a whole pathway hierarchy through a chain of weakly connected terms. The distance cutoff is a presentation choice, not a test; missing leading edges contribute no similarity.
#
# The representative is the member with the strongest primary spatial evidence. Every original term and its statistics remain in the exported tables. The data determine how many programs remain; we do not force 10–20 groups or assign biological program names automatically.

# %%
from pseudospace.pathway_remodeling import collapse_programs

programs = collapse_programs(candidate_ids, gene_sets, leading_edges, local_curves, distance_cut=.45)
programs = programs.merge(spatial[['pathway_id', 'effect', 'q_empirical', 'q_corr', 'correlation_support']], on='pathway_id', how='left')
representatives = programs.sort_values(['q_empirical', 'effect'], ascending=[True, False]).drop_duplicates('program')
programs['representative'] = programs.pathway_id.isin(representatives.pathway_id)
programs.to_csv(output / 'pathway_program_membership.csv', index=False)
display(representatives.head(20))

# %% [markdown]
# ### How do the candidate pathways cluster?
#
# This dendrogram uses the **same** four equal-weight similarities and complete linkage as the program assignments above. The dashed line is the existing distance cutoff of 0.45. Branch height describes pathway similarity, not statistical significance; labels are retained in the saved PDF for close inspection. The grouping is checked against the membership table before plotting.
#

# %%
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.spatial.distance import squareform

if len(candidate_ids) > 1:
    traces = {p: local_curves[local_curves.pathway_id.eq(p)].sort_values('position') for p in candidate_ids}
    distances = np.zeros((len(candidate_ids), len(candidate_ids)))
    def jaccard(left, right):
        left, right = set(left), set(right)
        return len(left & right) / len(left | right) if left | right else 0.
    for i, left in enumerate(candidate_ids):
        for j in range(i):
            right = candidate_ids[j]
            similarities = [jaccard(gene_sets[left], gene_sets[right]),
                            jaccard(leading_edges.get(left, []), leading_edges.get(right, []))]
            for column in ('divergence', 'direction'):
                a, b = traces[left][column].to_numpy(), traces[right][column].to_numpy()
                correlation = np.corrcoef(a, b)[0, 1] if np.std(a) > 1e-10 and np.std(b) > 1e-10 else 0.
                similarities.append(max(float(correlation), 0.))
            distances[i, j] = distances[j, i] = 1 - np.mean(similarities)
    tree = linkage(squareform(distances), method='complete')
    assigned = programs.set_index('pathway_id').loc[candidate_ids, 'program'].to_numpy()
    assert np.array_equal(fcluster(tree, .45, criterion='distance'), assigned)
    labels = [p.replace('::', ': ', 1) for p in candidate_ids]
    fig, ax = plt.subplots(figsize=(11, max(5, .23 * len(labels))), layout='constrained')
    dendrogram(tree, labels=labels, orientation='left', leaf_font_size=6,
               color_threshold=.45, above_threshold_color='0.7', ax=ax)
    ax.axvline(.45, color='0.3', ls='--', lw=.8)
    ax.set(xlabel='Complete-linkage distance (1 − mean similarity)',
           title='Exploratory candidate pathway clustering')
    fig.savefig(output / 'pathway_program_dendrogram.pdf', bbox_inches='tight')
    plt.show()
    plt.close(fig)
elif candidate_ids:
    print('One candidate: no pairwise pathway tree to display.')
else:
    print('No spatial candidates: pathway clustering was skipped.')


# %% [markdown]
# ## 12 · Ask whether the spatial evidence is fragile
#
# Refit all genes, then retest the full eligible pathway family under:
#
# 1. Leaving out each mouse and each human section in turn.
# 2. Spline basis sizes 4 and 8 instead of 6.
# 3. Raising the detection threshold from 2% to 5%.
#
# The coordinate and grid stay fixed. Each refit rebalances the remaining specimens; the detection sensitivity updates membership and the matched gene background. Every sensitivity run corrects **within each statistic across all eligible pathways from all three libraries**, so candidate-only BH cannot make retention easier.
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
# ### Which candidates survive the planned sensitivity runs?
#
# Green means the original spatial candidate remains positive with matched q ≤ 0.05 in that run; pale gray means testable but not retained; dark gray means unavailable. Rows are ordered by retained/testable runs, then the original matched q. The notebook samples 24 rows across the most, middle, and least retained candidates and saves the complete candidate matrix separately. These are section and model sensitivities, not donor replication.
#

# %%
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch

if candidate_ids:
    q_order = spatial.set_index('pathway_id').q_empirical
    ordered_ids = (stability_summary.assign(q=q_order)
        .sort_values(['retention_fraction', 'n_retained', 'q'], ascending=[False, False, True]).index.tolist())
    variant_order = [item[0] for item in variants]
    stability_lookup = stability.set_index(['pathway_id', 'variant'])
    def draw_stability(ids):
        values = np.array([[int(stability_lookup.loc[(pathway, variant), 'retained'])
                            if stability_lookup.loc[(pathway, variant), 'testable'] else -1
                            for variant in variant_order] for pathway in ids])
        fig, ax = plt.subplots(figsize=(10, max(4, .27 * len(ids) + 1.3)), layout='constrained')
        ax.imshow(values, aspect='auto', interpolation='none',
                  cmap=ListedColormap(['#555555', '#e5e5e5', '#009E73']),
                  norm=BoundaryNorm([-1.5, -.5, .5, 1.5], 3))
        labels = [spatial.set_index('pathway_id').loc[pathway, 'pathway'] for pathway in ids]
        ax.set(xticks=np.arange(len(variant_order)), xticklabels=variant_order,
               yticks=np.arange(len(ids)), yticklabels=labels,
               title='Candidate retention across planned sensitivities')
        ax.tick_params(axis='x', labelrotation=45)
        ax.legend(handles=[Patch(color=c, label=label) for c, label in
            [('#009E73', 'Retained'), ('#e5e5e5', 'Not retained'), ('#555555', 'Unavailable')]],
            loc='upper left', bbox_to_anchor=(1, 1), frameon=False)
        return fig
    middle = len(ordered_ids) // 2
    display_ids = list(dict.fromkeys(ordered_ids[:8] + ordered_ids[max(0, middle - 4):middle + 4] + ordered_ids[-8:]))
    display_ids.sort(key=ordered_ids.index)
    fig = draw_stability(display_ids)
    fig.savefig(output / 'pathway_stability_overview.pdf', bbox_inches='tight')
    plt.show()
    plt.close(fig)
    if len(ordered_ids) > 24:
        fig = draw_stability(ordered_ids)
        fig.savefig(output / 'pathway_stability_all_candidates.pdf', bbox_inches='tight')
        plt.close(fig)


# %% [markdown]
# ## 14 · Final evidence table and reproducibility record
#
# The final table retains **all original pathway statistics** and adds program membership, bulk comparison, spatial descriptions, stability, and correlation-adjusted spatial evidence. It is ready for biological review, not a list of validated mechanisms.
#
# > **Before writing a biological conclusion**
# >
# > Check coverage and null resolution; inspect absolute gene directions and fitted curves; consider specimen and basis sensitivity; keep correlation and one-human-donor limitations beside the result. Neither a small pathway q-value nor a stable section-level pattern is population-level human inference.
#
# The table preserves both continuous-only candidates and pathways for which continuous analysis adds detail to a bulk-accessible difference. Program grouping and driver selection happen after discovery.
#
# Correlation support annotates the primary spatial screen only. The primary matched-null candidate rule, bulk comparison, program representatives, and stability definition remain unchanged; `q_corr` is not an extra discovery or retention filter.

# %%
final_table = comparison.merge(programs[['pathway_id', 'program', 'representative']], on='pathway_id', how='left')
final_table = final_table.merge(stability_summary, on='pathway_id', how='left')
final_table = final_table.merge(coverage[['pathway_id', 'library', 'pathway', 'n_requested', 'n_assayed',
    'n_tested', 'measured_fraction', 'tested_fraction']], on='pathway_id', how='left')
for role in ('broad_supporter', 'leading_edge', 'spatial_driver'):
    role_genes = member_evidence[member_evidence[role]].groupby('pathway_id').gene.agg(';'.join)
    final_table[role + '_genes'] = final_table.pathway_id.map(role_genes).fillna('')
final_table = final_table.merge(spatial[correlation_columns + ['correlation_support']],
    on='pathway_id', how='left', validate='one_to_one')
final_table['caveat'] = 'Exploratory gene-set evidence; two mouse specimens, two cortex sections from one human donor.'
final_table.to_csv(output / 'pathway_evidence_atlas.csv', index=False)

conventional_review = benchmark.merge(final_table[['pathway_id', 'program', 'q_corr',
    'correlation_support', 'n_testable', 'n_retained', 'retention_fraction']],
    on='pathway_id', how='left', validate='one_to_one')
conventional_review.to_csv(output / 'conventional_pt_vs_continuous_evidence_review.csv', index=False)
display(conventional_review.loc[conventional_review.comparison_group.eq('continuous screen only')]
    .sort_values('spatial_q_shared')[['pathway_id', 'spatial_q_shared', 'best_segment_q', 'best_magnitude_q',
        'peak_position', 'q_corr', 'correlation_support', 'n_retained', 'n_testable']].head(20))

manifest = {'logic_version': NOTEBOOK_LOGIC_VERSION, 'implementation_fingerprint': cache_code,
    'fit_input_fingerprint': digest(fit_inputs), 'input_fingerprint': digest(input_path),
    'ortholog_fingerprint': digest(map_path),
    'source_count_fingerprints': {name: digest(path) for name, path in source_paths.items()},
    'expression_reconstruction': dict(adata.uns['expression_reconstruction']),
    'retained_structures_fingerprint': digest(adata.obs_names.to_numpy()),
    'gene_filter_flow': gene_filter_flow.to_dict(orient='records'), 'expression_fingerprint': fit_inputs['Y'],
    'coordinate': '13 scFates shared_pseudospace; unchanged', 'common_support': [float(lo), float(hi)],
    'n_genes': len(genes), 'n_pathways': len(gene_sets), 'basis_df': basis_df,
    'detection_threshold': min_detection, 'matching_bins': 3, 'n_null': n_null, 'gsea_permutations': gsea_permutations, 'seed': seed,
    'libraries': {name: digest(library_dir / f'{name}.json') for name in libraries},
    'packages': {name: version(name) for name in ['numpy', 'scipy', 'pandas', 'anndata', 'patsy', 'statsmodels', 'gseapy']},
    'n_candidates': len(candidate_ids),
    'correlation_test': 'spatial CAMERA-style rank variance; equal-specimen Fisher mean; normal tail',
    'empirical_bh_families': 'separate T_spatial, T_level, T_total; all tested pathways across libraries',
    'correlation_bh_family': 'all tested spatial pathways across libraries',
    'correlation_candidate_counts': correlation_candidates.correlation_support.value_counts().to_dict(),
    'caveat': final_table.caveat.iloc[0]}
(output / 'run_manifest.json').write_text(json.dumps(manifest, indent=2))
display(final_table[final_table.pathway_id.isin(candidate_ids)].sort_values('q_empirical_T_spatial').head(20))
print('Saved evidence atlas and run manifest to', output)


# %% [markdown]
# ## 15 · Freeze discovery and begin paper interpretation
#
# > **A new stage, with the original analysis preserved**
# >
# > Sections 1–14 establish the coordinate, gene statistics, candidate pathways, correlation annotation, and robustness. Everything below interprets those results. It does not refit the coordinate, change the spatial AUC, recalculate discovery q-values, or replace the earlier program table.
#
# The paper question is: **which recurring groups of genes explain the human–mouse differences, where do they differ along PT pseudospace, and what does this reveal beyond an overall species contrast?**
#
# We retain **every original candidate**. Correlation support and retention across sensitivity runs are evidence annotations. For reading order only, a `paper_priority` flag later marks correlation-supported terms retained in at least 75% of **all planned** runs; unavailable runs therefore cannot make robustness look better. This is not a second discovery rule.
#
# > **What these data can support**
# >
# > There are two healthy mouse specimens and two healthy **cortex** sections from **one human donor**. The medulla-labelled section is also cortex. Programs, driver genes, and direction descriptions are exploratory observations in this cohort, not validated mechanisms or population-level species effects.
#
# Run the earlier sections first. The new files live in a fingerprinted `paper_interpretation/` folder, leaving the discovery exports in place. The fingerprint covers the frozen inputs and interpretation implementation; rerunning the same interpretation preserves a manually edited review worksheet.

# %%
from pseudospace import pathway_programs as paper_helpers

# Copy the candidate rows so paper annotations cannot change the discovery table.
frozen_candidates = final_table[final_table.pathway_id.isin(candidate_ids)].copy()
assert set(frozen_candidates.pathway_id) == set(candidate_ids)
paper_input_hash = digest({'atlas': final_table, 'members': member_evidence,
    'curves': local_curves, 'gene_statistics': gene_stats, 'gene_sets': gene_sets,
    'fits': {key: fit[key] for key in ('grid', 'z', 'human', 'mouse', 'delta')}, 'manifest': manifest})
# Bump this version if changing the appended notebook code or its settings.
paper_version = '12-paper-interpretation-v1'
paper_run = digest({'inputs': paper_input_hash, 'version': paper_version,
    'helper': project / 'pseudospace' / 'pathway_programs.py'})
paper_output = output / 'paper_interpretation' / paper_run
paper_output.mkdir(parents=True, exist_ok=True)
frozen_candidates.to_csv(paper_output / 'frozen_candidates.csv', index=False)
print(f'{len(frozen_candidates)} frozen candidates. Interpretation folder: {paper_output}')

# %% [markdown]
# ## 16 · Ask whether frequently annotated genes dominate the effect
#
# A gene appearing in many pathway definitions can make several terms look like separate discoveries. We add **one overlap sensitivity**, inspired by the downweighting idea in [PADOG](https://tarcalab.med.wayne.edu/software/padog/), without implementing PADOG or adding another enrichment test.
#
# For each gene, count the number of **all tested terms across the loaded libraries** containing it, and give it weight $w_g=1/\sqrt{f_g}$. For each pathway, compute each member's probability of outranking the same nonmember background on $T_{\mathrm{spatial}}$, with half credit for ties. Average those probabilities using the weights, then subtract 0.5.
#
# > **How to read the sensitivity**
# >
# > The original AUC effect is the unweighted mean of the same member contributions. A negative `overlap_effect_change` means downweighting broadly annotated members weakens that effect. A positive rank change means the pathway falls in the **effect ranking across all tested pathways**. These are descriptive changes, not new p-values; the original matched-null q-value does not apply to the weighted effect.
#
# Library redundancy affects annotation frequency, so this check is conditional on this exact library collection. We save the frequencies and weights. No pathway is removed because of this sensitivity.

# %%
overlap_check, annotation_weights = paper_helpers.overlap_sensitivity(gene_stats.T_spatial, gene_sets)
# Confirm that the sensitivity starts from the original spatial effect.
original_check = overlap_check.merge(final_table[['pathway_id', 'effect_T_spatial']], on='pathway_id')
assert np.allclose(original_check.original_effect, original_check.effect_T_spatial)
overlap_check.to_csv(paper_output / 'overlap_sensitivity.csv', index=False)
annotation_weights.to_csv(paper_output / 'annotation_weights.csv', index=False)
display(overlap_check[overlap_check.pathway_id.isin(candidate_ids)]
    .sort_values('overlap_effect_change').head(12))

# %% [markdown]
# ## 17 · Group pathways by their active genes and spatial behavior
#
# [EnrichmentMap](https://doi.org/10.1371/journal.pone.0013984) motivates reading related terms together. Here, the grouping uses the evidence that drove the spatial signal:
#
# | Ingredient | Definition | Weight in grouping |
# |---|---|---:|
# | Active genes | Union of the existing positive GSEA leading edge, broad supporters, and local spatial drivers | 50%: Jaccard overlap |
# | Divergence $D(s)$ | Original unsigned local enrichment curve | 25%: positive Pearson correlation |
# | Direction $S(s)$ | Original relative signed enrichment curve | 25%: positive Pearson correlation |
# | Full membership | All tested genes assigned to each term | Supplemental only; 0% |
#
# A constant curve supplies no correlation evidence; missing curves or active genes stop the analysis; negative correlations contribute zero. The active-gene definitions come directly from section 10, so we introduce no new driver-selection thresholds.
#
# Complete linkage at distance **0.55** requires every pair within a group to have combined similarity of at least **0.45**. This is a transparent, descriptive resolution choice, not an optimized biological boundary. We do **not** tune it to obtain 8–20 groups. Singletons remain visible; broader biological themes can be proposed during review.
#
# > **A grouping is not yet a named biological program**
# >
# > These are draft program groups. Similarity can reflect shared genes, similar curves, or both. The pairwise table exposes each component. A network edge passes the same similarity threshold, but an edge connecting two groups does not merge them: complete linkage is stricter than connected components.
#
# Every candidate becomes a network node. The exported edge table can be read by Cytoscape; the complete pairwise table also retains absent edges. Program IDs depend on their exact member terms, allowing annotations to be traced to a specific grouping.

# %%
paper_membership, paper_pairs, active_sets = paper_helpers.active_programs(
    candidate_ids, gene_sets, member_evidence, local_curves, distance_cut=.55)
paper_membership.to_csv(paper_output / 'program_membership.csv', index=False)
paper_pairs.to_csv(paper_output / 'pathway_pairwise_similarity.csv', index=False)
paper_pairs[paper_pairs.network_edge].to_csv(paper_output / 'network_edges.csv', index=False)
assert set(paper_membership.pathway_id) == set(candidate_ids)
print(f'{len(candidate_ids)} terms → {paper_membership.paper_program.nunique()} draft groups.')
display(paper_membership.groupby('paper_program').size().rename('n_terms').value_counts().sort_index())

# %% [markdown]
# ### Why do related terms form draft programs?
#
# The matrix uses the **same combined similarity** as the descriptive grouping. It shows terms from the largest multi-term groups; all pairwise values remain in the exported table. Blocks within a group help explain the shared program label, while similarity across groups need not merge them under complete linkage.
#

# %%
multi = paper_membership.groupby('paper_program').filter(lambda rows: len(rows) > 1)
if len(multi) > 1:
    selected_programs = multi.paper_program.value_counts().head(4).index
    selected_ids = multi[multi.paper_program.isin(selected_programs)].pathway_id.head(24).tolist()
    selected_ids = sorted(selected_ids, key=lambda pathway: (
        paper_membership.set_index('pathway_id').loc[pathway, 'paper_program'],
        spatial.set_index('pathway_id').loc[pathway, 'q_empirical'], pathway))
    similarity = pd.DataFrame(np.eye(len(selected_ids)), index=selected_ids, columns=selected_ids)
    for row in paper_pairs.itertuples():
        if row.source in similarity.index and row.target in similarity.columns:
            similarity.loc[row.source, row.target] = similarity.loc[row.target, row.source] = row.similarity
    fig, ax = plt.subplots(figsize=(max(7, .4 * len(selected_ids) + 3),
                                    max(6, .4 * len(selected_ids) + 2)), layout='constrained')
    artist = ax.imshow(similarity, cmap='viridis', vmin=0, vmax=1, interpolation='none')
    names = [spatial.set_index('pathway_id').loc[pathway, 'pathway'][:38] for pathway in selected_ids]
    ax.set(xticks=np.arange(len(names)), xticklabels=names,
           yticks=np.arange(len(names)), yticklabels=names,
           title='Similarity used for draft program grouping')
    ax.tick_params(axis='x', labelrotation=90, labelsize=7)
    ax.tick_params(axis='y', labelsize=7)
    fig.colorbar(artist, ax=ax, label='Combined active-gene and curve similarity')
    fig.savefig(paper_output / 'program_similarity_overview.pdf', bbox_inches='tight')
    plt.show()
    plt.close(fig)


# %% [markdown]
# ## 18 · Choose draft representatives and retain bulk context
#
# We choose one or two terms to make each packet readable, while retaining **all terms** in the supplementary table. The first representative maximizes the equally weighted percentile ranks of spatial effect, tested/requested coverage, retention over all planned sensitivity runs, and coverage of the group's shared active genes. The shared core contains genes active in at least two terms; if there is no shared core, use the active union.
#
# A second representative must add at least 10% of the group's active-gene union beyond the first. Among eligible second terms, the same score determines selection. Ties resolve by pathway ID. These are **draft computational representatives**: biological interpretability still requires reviewing the term names and genes. No automatic label is presented as a validated mechanism.
#
# > **What did bulk already know?**
# >
# > Each representative retains its original bulk, level, and spatial effects and q-values. `continuous-only candidate` means the spatial screen found a candidate while the chosen bulk/level screens did not; it does not prove absence of a bulk effect. `bulk/level plus spatial structure` means pseudospace adds localization or restructuring to a difference already accessible to bulk/level analysis.
# >
# > The original `bulk/level captures difference` terms remain in the complete tested atlas. They are not silently promoted into spatial programs. A program whose terms span several comparison classes is labelled **mixed context**, rather than being assigned the most favorable class.

# %%
paper_terms = paper_helpers.representative_terms(paper_membership, frozen_candidates, active_sets)
paper_terms = paper_terms.merge(overlap_check, on='pathway_id', validate='one_to_one')
paper_terms.to_csv(paper_output / 'all_candidate_terms.csv', index=False)
paper_terms.to_csv(paper_output / 'network_nodes.csv', index=False)
final_table.to_csv(paper_output / 'complete_discovery_atlas.csv', index=False)

paper_summary = paper_terms.groupby('paper_program').agg(
    n_terms=('pathway_id', 'size'), n_priority_terms=('paper_priority', 'sum'),
    n_correlation_supported=('correlation_support', lambda x: x.eq('Correlation-supported').sum()),
    median_planned_retention=('planned_retention', 'median'),
    minimum_overlap_effect_change=('overlap_effect_change', 'min'),
    bulk_context=('comparison_group', lambda x: x.iloc[0] if x.nunique() == 1 else 'mixed context'))
paper_representatives = paper_terms[paper_terms.representative_order.gt(0)].sort_values(
    ['paper_program', 'representative_order'])
paper_summary['draft_label'] = paper_representatives.groupby('paper_program').pathway.first()
paper_summary['representative_terms'] = paper_representatives.groupby('paper_program').pathway_id.agg(';'.join)
paper_summary = paper_summary.sort_values(['n_priority_terms', 'n_terms'], ascending=False)
display(paper_summary.head(15))

# %% [markdown]
# ## 19 · Open each program: shared drivers, private genes, and mixed directions
#
# A pathway's signed average can hide simultaneous human-high and mouse-high genes. Following the **directionality distinction**, rather than the tests, in [piano](https://doi.org/10.1093/nar/gkt111), we describe both sides separately.
#
# At each grid point, count active genes with $Z\geq2$ and $Z\leq-2$. A side is descriptively prominent when it includes at least **two genes and 20% of active genes**. Both sides prominent means `mixed`; one side prominent means `human-high dominated` or `mouse-high dominated`; otherwise the label is `weak / sparse`. These are display rules for working statistics, **not local significance thresholds**. We export the actual counts and fractions so the labels can be checked.
#
# > **Two different patterns**
# >
# > Simultaneous opposing genes at the same position produce a mixed pattern. Human-high dominance in one region and mouse-high dominance elsewhere produces a positional direction switch. Neither pattern should be compressed into a single signed mean. $S(s)$ itself is relative to background genes and is not proof that most pathway members are absolutely human-high.
#
# A **shared active driver** supports at least two terms within the draft group. A **term-specific active gene** supports exactly one term within a multi-term group; it need not be unique elsewhere in the atlas. Single-term groups are labelled separately. Large local $|Z|$ identifies genes to inspect, not necessarily a large expression effect: check the saved $\delta$, detection information in the member table, and actual fitted curves.

# %%
paper_genes, paper_directions = paper_helpers.program_gene_evidence(
    paper_terms, member_evidence, genes, grid, fit, z_threshold=2., fraction=.2)
paper_genes.to_csv(paper_output / 'program_gene_evidence.csv', index=False)
paper_directions.to_csv(paper_output / 'program_direction_by_position.csv', index=False)
member_evidence[member_evidence.pathway_id.isin(candidate_ids)].to_csv(
    paper_output / 'candidate_member_evidence.csv', index=False)
paper_summary['n_active_genes'] = paper_genes.groupby('paper_program').gene.nunique()
paper_summary['n_shared_active_genes'] = paper_genes[paper_genes.n_active_terms.ge(2)].groupby('paper_program').size()
paper_summary['n_shared_active_genes'] = paper_summary.n_shared_active_genes.fillna(0).astype(int)
paper_summary['mixed_grid_fraction'] = paper_directions.assign(
    mixed=paper_directions.direction_description.eq('mixed')).groupby('paper_program').mixed.mean()
paper_summary['positional_direction_switch'] = paper_directions.groupby('paper_program').direction_description.agg(
    lambda x: {'human-high dominated', 'mouse-high dominated'}.issubset(set(x)))
paper_summary.to_csv(paper_output / 'program_summary.csv')
display(paper_summary.head(15))

# %% [markdown]
# ### Preserve the numerical source of every panel
#
# The figure source contains all active genes at every fitted grid position, including unclipped signed $Z$, human and mouse fitted expression, and their difference. The pathway source contains the original candidate $D(s)$ and $S(s)$ curves. Program membership and representative tables link these sources to every packet. No residual matrix or new fit is needed here.

# %%
paper_gene_order = sorted(paper_genes.gene.unique())
paper_gene_rows = gene_index.get_indexer(paper_gene_order)
assert (paper_gene_rows >= 0).all()
paper_gene_source = pd.DataFrame({'gene': np.repeat(paper_gene_order, len(grid)),
    'position': np.tile(grid, len(paper_gene_order))})
for measure in ('z', 'human', 'mouse', 'delta'):
    paper_gene_source[measure] = fit[measure][paper_gene_rows].ravel()
paper_gene_source.to_csv(paper_output / 'gene_curve_source.csv.gz', index=False)
local_curves[local_curves.pathway_id.isin(candidate_ids)].to_csv(
    paper_output / 'pathway_curve_source.csv', index=False)

# %% [markdown]
# ## 20 · Generate an evidence packet for every draft program
#
# **Figure question:** which genes produce this program's spatial difference, and what does position add to its bulk comparison?
#
# | Panel | Evidence it contributes | Reading rule |
# |---|---|---|
# | Representative $D(s)$ and $S(s)$ | Where unsigned divergence and relative direction occur | R1/R2 map to full term names and bulk context in the packet caption |
# | Active-gene heatmap | Which genes differ at which positions, including opposing signs | Order by each gene's maximum $|Z|$ position; `*` marks shared active drivers |
# | Human and mouse gene fits | Whether the local statistic corresponds to a plausible expression pattern | Up to four examples, balancing shared and term-specific roles when available |
#
# A prioritized 15-gene heatmap is shown first; supplementary heatmaps include **every active gene**, paginated at 36 rows for readability. A shared symmetric color scale saturates at $Z=\pm6$; source tables retain the full values. D and S use common axes across programs. The curve examples take the two strongest peak-$|Z|$ genes per role, then fill remaining slots up to four by peak $|Z|$; all remaining genes are retained in the source table and heatmaps.
#
# > **Figure scope and limitations**
# >
# > These are paper-planning evidence panels, not a final manuscript figure selection. Human and mouse curves are the existing equal-specimen fitted means with zero specimen contrasts. No donor-level confidence interval can be justified from one human donor, so these descriptive curves have no uncertainty ribbons. HC3 $Z$ does not account for spatial dependence, donor replication, or coordinate uncertainty. Differences in sampled structures or alignment remain alternative explanations.
#
# Exports are editable PDF/SVG panels at 183 mm width, with 300-dpi PNG previews. Each program receives a caption with its representative terms, bulk context, robustness, overlap sensitivity, and cohort caveat. All packets are saved; only the first two groups in the evidence-based reading order are displayed inline to keep the notebook manageable. Change `preview_programs` here to inspect others.
#

# %%
from IPython.display import Markdown

preview_programs = paper_summary.head(2).index.tolist()
for paper_program in paper_summary.index:
    packet = paper_output / paper_program
    packet.mkdir(exist_ok=True)
    packet_terms = paper_representatives[paper_representatives.paper_program.eq(paper_program)]
    caption = [f'# {paper_program}: draft program evidence',
        f"Bulk context: {paper_summary.loc[paper_program, 'bulk_context']}.",
        'Representatives are computational suggestions; biological names and interpretation require review.']
    for row in packet_terms.itertuples():
        caption.append(f'R{row.representative_order}: {row.pathway}. {row.comparison_group}; '
            f'spatial AUC effect {row.effect_T_spatial:.3f}, matched q {row.q_empirical_T_spatial:.3g}; '
            f'bulk q {row.q_empirical_T_bulk:.3g}, level q {row.q_empirical_T_level:.3g}; '
            f'{row.correlation_support}, correlation q {row.q_corr:.3g}; '
            f'retained {int(row.n_retained)}/{int(row.n_planned)} planned runs '
            f'({int(row.n_testable)} testable); overlap effect change {row.overlap_effect_change:.3f}.')
    caption.extend(['D and S are original pathway curves, not a new program score. '
        'S is relative to the gene background; consult absolute signed Z as well.',
        'Heatmaps: the primary panel shows up to 15 prioritized genes; supplementary pages show every active gene. '
    'Rows follow peak |Z| position; * = active in two or more group terms. '
        'Colors saturate at ±6; raw values are in gene_curve_source.csv.gz. '
        'Each gene occurs once per program even if it supports several terms.',
        'Gene fits: up to four examples chosen by peak |Z| within shared/term-specific roles; '
        'human solid orange, mouse dashed blue. Original balanced model means, no uncertainty ribbons. '
        'Expression is log-normalized; early/late refer to the established coordinate, not validated anatomical boundaries.',
        'Exploratory: two healthy mouse specimens; two healthy cortex sections from ONE human donor. '
        'Structural HC3 statistics are not donor-level tests. No new p-values or claims of pathway activation.',
        'Source: pathway_curve_source.csv, gene_curve_source.csv.gz, program_gene_evidence.csv, '
        'candidate_member_evidence.csv, all_candidate_terms.csv in the parent folder.'])
    (packet / 'caption.md').write_text('\n\n'.join(caption))
    if paper_program in preview_programs:
        display(Markdown('\n\n'.join(caption)))
    for panel_name, paper_figure in paper_helpers.program_figures(
            paper_program, paper_terms, paper_genes, local_curves, genes, grid, fit):
        paper_figure.savefig(packet / f'{panel_name}.pdf')
        paper_figure.savefig(packet / f'{panel_name}.svg')
        paper_figure.savefig(packet / f'{panel_name}.png', dpi=300)
        if paper_program in preview_programs and panel_name in {'pathway_curves', 'top_genes', 'gene_fits'}:
            plt.show()
        plt.close(paper_figure)
print(f'Saved evidence packets for {len(paper_summary)} draft groups.')


# %% [markdown]
# ## 21 · Turn the packets into a paper outline through biological review
#
# The worksheet is deliberately unfinished: the data can nominate a program but cannot supply a trustworthy biological label or mechanism by itself. Review the **genes and curves**, not just the pathway title. Edit the CSV outside the notebook; reruns of the same fingerprint keep it intact. Changing inputs, helper code, or `paper_version` creates a new folder rather than applying old annotations to new groups.
#
# | Review question | What to record |
# |---|---|
# | Is this one coherent biological theme? | Reviewed program label, or a proposed merge/split with a reason |
# | Are the representatives interpretable? | One or two reviewed pathway IDs and the rationale; defaults are suggestions |
# | What does pseudospace add? | Localization, directional reorganization, or little added information over bulk/level |
# | Which genes carry the claim? | Shared versus term-specific genes; identify domination by a small gene family |
# | Could the pattern reflect the model or sampling? | Detection, support, specimen sensitivity, overlap sensitivity, and coordinate concerns |
# | Is there independent support? | Verified literature references and a specific independent-donor or experimental validation plan |
#
# > **Suggested paper structure**
# >
# > 1. Establish the fixed coordinate, shared support, cohort, and analysis question using the existing workflow.
# > 2. Show the frozen pathway screen with effect sizes, bulk comparison, correlation annotations, and robustness; retain the complete tested family in supplements.
# > 3. Use a small, biologically reviewed set of program packets to demonstrate what continuous position adds. Include mixed-direction programs where appropriate, rather than selecting only clean signed changes.
# > 4. Separate observed expression patterns from proposed biological explanations. Independent human donors and external validation are subsequent evidence, not something this notebook has already supplied.
#
# The number of programs is an outcome of the documented grouping rule. A smaller main-figure selection should follow biological review and complementary evidence, not a target count or the smallest q-values. Correlation-sensitive and unstable terms remain in the supplements and their caveats travel with them.

# %%
review_path = paper_output / 'biological_review.csv'
if not review_path.exists():
    review_sheet = paper_summary.reset_index()[['paper_program', 'draft_label', 'representative_terms', 'bulk_context']]
    for review_column in ('reviewed_label', 'reviewed_representative_ids', 'merge_split_proposal',
                         'gene_based_interpretation', 'what_position_adds', 'alternative_explanations',
                         'literature_references', 'external_validation_plan', 'main_figure_rationale', 'review_status'):
        review_sheet[review_column] = ''
    review_sheet.to_csv(review_path, index=False)
paper_review = pd.read_csv(review_path, keep_default_na=False)
assert paper_review.paper_program.is_unique and set(paper_review.paper_program) == set(paper_summary.index)

# Record interpretation choices separately from the untouched discovery manifest.
paper_manifest = {'version': paper_version, 'run': paper_run, 'frozen_input_hash': paper_input_hash,
    'discovery_manifest': manifest, 'interpretation_code': digest(project / 'pseudospace' / 'pathway_programs.py'),
    'overlap': 'member weights 1/sqrt(frequency in all tested terms); unchanged nonmember background; no p-values',
    'similarity_weights': {'active_jaccard': .5, 'positive_D_correlation': .25, 'positive_S_correlation': .25},
    'linkage': 'complete', 'distance_cut': .55, 'priority_planned_retention': .75,
    'direction_description': {'absolute_z': 2, 'minimum_fraction': .2, 'minimum_genes': 2},
    'heatmap_z_saturation': 6, 'n_frozen_candidates': len(candidate_ids),
    'n_draft_groups': len(paper_summary), 'n_priority_terms': int(paper_terms.paper_priority.sum()),
    'review_status': 'Draft groups and representatives; worksheet is not automatically applied to figures.'}
(paper_output / 'interpretation_manifest.json').write_text(json.dumps(paper_manifest, indent=2))
# Final guard: no discovery input was mutated by interpretation.
assert paper_input_hash == digest({'atlas': final_table, 'members': member_evidence,
    'curves': local_curves, 'gene_statistics': gene_stats, 'gene_sets': gene_sets,
    'fits': {key: fit[key] for key in ('grid', 'z', 'human', 'mouse', 'delta')}, 'manifest': manifest})
display(paper_review.head(10))
print('Biological review worksheet:', review_path)
