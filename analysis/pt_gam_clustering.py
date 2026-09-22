# %% [markdown]
# # Gene-level GAM curve clustering of the human vs healthy-mouse PT pseudospace
#
# Unsupervised, **gene-level** discovery of recurring continuous PT expression programs, from the GAM
# curves notebook 03 already validated. Two catalogs are built, plus the map between them:
#
# ```
# normal mouse PT organisation     mouse GAMs -> positional PCA -> positional modules
# cross-species remodelling        human GAM - mouse GAM -> remove the mean species offset
#                                  -> response PCA -> response modules
# biological interpretation        positional module x response module -> pathway enrichment
# ```
#
# The conceptual outputs are `gene -> normal mouse positional module` and
# `gene -> human-vs-mouse response module`. Together they say what spatial PT program a gene normally
# belongs to **and** how that program changes in human.
#
# This is deliberately the opposite direction of travel from a pathway-first analysis: here the
# clusters are discovered first and pathway annotations are loaded only **after** every cluster
# assignment and every cluster-number decision is frozen. No pathway is ever aggregated into an
# expression curve, and pathway membership never influences which genes land in which module.
#
# **Input.** This notebook consumes notebook 03's PT-only artifact (`cross_species_pt_dpt.h5ad`) and
# its accepted ortholog map. It never re-does integration, cell typing, the diffusion map, DPT, the
# PT trajectory or the ortholog definition, and it never includes a nephron segment outside PT.
#
# **Inference unit.** The comparison is 4 specimens: 2 mouse specimens and 2 human *sections of one
# donor*. Averaging the two human sections with equal weight stops one section dominating because it
# contains more PT structures; it does **not** create independent human biological replication. Every
# number below is a descriptive summary of a 2-vs-1-donor comparison, not a species-level test.
#
# Run notebook 03 first: this notebook stops with a clear error if its artifacts are missing.
#

# %% [markdown]
# ## Setup
#
# Configuration, the upstream dependencies, and the parameters every analysis below shares.

# %%
# Purpose: bootstrap the run: cache redirects, repository root, roots, and output directories.
# No IPython magics here on purpose: notebooks 02, 04, 05 and 06 carry none, and a magic is a live
# line in the .ipynb that the synthetic dry-run harness cannot exec.

import argparse
import json
import os
import sys
import warnings
from pathlib import Path

# These must be set before importing matplotlib, numba, Scanpy, or modules that import them.
CACHE_ROOT = Path(os.environ.get('PSEUDOSPACE_CACHE_ROOT', '/tmp/pseudospace_pt_gam_clustering_cache'))
for _subdir in ('matplotlib', 'numba', 'xdg'):
    (CACHE_ROOT / _subdir).mkdir(parents=True, exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR', str(CACHE_ROOT / 'matplotlib'))
os.environ.setdefault('NUMBA_CACHE_DIR', str(CACHE_ROOT / 'numba'))
os.environ.setdefault('XDG_CACHE_HOME', str(CACHE_ROOT / 'xdg'))


def _find_project_dir() -> Path:
    starts = [Path.cwd().resolve()]
    if '__file__' in globals():
        starts.insert(0, Path(__file__).resolve().parent)
    for start in starts:
        for candidate in (start, *start.parents):
            if (candidate / 'pseudospace').is_dir() and (candidate / 'data').is_dir():
                return candidate
    raise RuntimeError('Could not locate the repository root containing pseudospace/ and data/')


def _workflow_roots(project_dir: Path) -> tuple[Path, Path]:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--data-root', type=Path)
    parser.add_argument('--results-root', type=Path)
    args, _ = parser.parse_known_args()
    data_root = args.data_root or os.environ.get('PSEUDOSPACE_DATA_ROOT') or project_dir / 'data'
    results_root = (args.results_root or os.environ.get('PSEUDOSPACE_RESULTS_ROOT')
                    or project_dir / 'results')
    return Path(data_root).expanduser().resolve(), Path(results_root).expanduser().resolve()


PROJECT_DIR = _find_project_dir()
if str(PROJECT_DIR) not in sys.path:
    sys.path.insert(0, str(PROJECT_DIR))
DATA_ROOT, RESULTS_ROOT = _workflow_roots(PROJECT_DIR)


def _rel(path) -> str:
    """Path relative to the repository when it is inside it, absolute otherwise.

    A run may legitimately point the results root outside the repository (the synthetic dry-run
    harness does), and an unconditional ``relative_to`` would then raise instead of printing.
    """
    path = Path(path)
    try:
        return str(path.relative_to(PROJECT_DIR))
    except ValueError:
        return str(path)


# Notebook 03 owns the pseudospace, the PT subset, the reviewed labels and the ortholog map. This
# notebook reads them and never writes into that directory.
UPSTREAM_DIR = RESULTS_ROOT / 'human_vs_healthy_mouse'
OUTPUT_DIR = RESULTS_ROOT / 'pt_gam_clustering'
DIAGNOSTIC_DIR = OUTPUT_DIR / 'diagnostics'
for _directory in (OUTPUT_DIR, DIAGNOSTIC_DIR):
    _directory.mkdir(parents=True, exist_ok=True)

DPT_OUTPUT_PATH = UPSTREAM_DIR / 'cross_species_pt_dpt.h5ad'
ORTHOLOG_MAP_PATH = UPSTREAM_DIR / 'ortholog_map_used.csv'

STAGE_CACHE_ENABLED = os.environ.get('PSEUDOSPACE_STAGE_CACHE', '1').strip().lower() not in ('0', 'false', 'no', '')
# 03's cache already holds the pooled nested per-gene fit; this notebook takes only its GCV lambda from
# that entry (so it must be a hit, not a recomputation) and stores its own per-specimen fits beside it.
STAGE_CACHE_DIR = UPSTREAM_DIR / 'stage_cache'
STAGE_CACHE_DIR.mkdir(parents=True, exist_ok=True)

print(f'Project: {PROJECT_DIR.name}')
print(f'Reading the PT artifact from: {_rel(UPSTREAM_DIR)}')
print(f'Results will be written to: {_rel(OUTPUT_DIR)}')
print(f'Stage cache: {"on" if STAGE_CACHE_ENABLED else "off"} ({_rel(STAGE_CACHE_DIR)})')


# %%
# Purpose: import the analysis stack and the project helpers this notebook reuses.
# The GAM, curve-module, pathway and cache machinery all already exist in `pseudospace`; nothing here
# re-implements a fit, a membership resolver or a caching scheme.
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from IPython.display import display
from scipy.cluster.hierarchy import fcluster, linkage
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score, silhouette_score

from pseudospace.levelshape import fit_single_condition_curves, summarize_curve_effects
from pseudospace.modules import (
    curve_descriptors,
    difference_curves,
    enrich_modules,
    module_stability,
)
from pseudospace.pathways import build_pathway_membership, summarize_pathway_redundancy
from pseudospace.specimen import specimen_balanced_curves
from pseudospace.stats_gam import (
    as_csr,
    gam_internal_knots,
    safe_spearman,
    zscore_rows,
)
from pseudospace.stage_cache import (
    cached_frame,
    cached_payload,
    cached_run_level_shape,
    code_digest,
    digest,
)

import matplotlib as _mpl
_mpl_cfg = Path(_mpl.get_configdir()).resolve()
if not _mpl_cfg.is_relative_to(CACHE_ROOT.resolve()):
    warnings.warn(
        f'matplotlib is caching to {_mpl_cfg}, not under {CACHE_ROOT.resolve()}. '
        'Run the bootstrap cell before importing Scanpy/matplotlib.', stacklevel=2)
else:
    print(f'Cache redirects active ({CACHE_ROOT}).')


# %%
# Purpose: every threshold this notebook applies, as one named knob each.
# Fitting conventions are byte-identical to notebook 03's Section 4. That is what makes the pooled fit
# a stage-cache hit rather than a re-fit, and what makes the curves here the same curves 03 and 05
# describe. Changing them here would silently produce different curves from the rest of the project.
GAM_CONFIG = {
    'n_internal_knots': 9,
    'lambda_grid': np.logspace(-3, 3, 13),
    'common_support_pct': (1, 99),
    'grid_points': 101,
    'min_detected_fraction': 0.02,
}

# Gene-universe and module-eligibility gates. None of them is a p-value, none is a differential
# expression call, and none may be tuned on pathway results.
UNIVERSE_CONFIG = {
    'measurement': 'lognorm',
    'min_detected_fraction': GAM_CONFIG['min_detected_fraction'],
    # A nearly flat curve must not be standardised and clustered: dividing it by its own tiny spread
    # turns numerical noise into an apparently strong shape. These are the absolute-magnitude gates
    # applied BEFORE scaling, in log-normalised expression units.
    'positional_min_peak_to_peak': 0.05,
    'response_min_peak_to_peak': 0.05,
}

# Clustering is hierarchical (deterministic: one Ward linkage, cut at k) on the retained PC scores.
CLUSTER_CONFIG = {
    'k_values': tuple(range(3, 13)),
    'pc_variance_target': 0.90,
    'pc_max_components': 15,
    # On raw scores PC1 carries ~64% of the variance and would decide membership alone; scaling each
    # retained component to unit variance gives the minor shape modes equal weight in the distance.
    'pc_whiten': True,
    'min_module_size': 30,
    'silhouette_tolerance': 0.01,
    # silhouette_score is O(n^2) in the number of features; a fixed subsample keeps the k sweep cheap
    # and the diagnostic deterministic.
    'silhouette_sample_size': 2000,
    'ward_method': 'ward',
}

PATHWAY_CONFIG = {
    'libraries': ['Reactome_2022', 'MSigDB_Hallmark_2020', 'KEGG_2019_Mouse'],
    'min_genes': 10,
    'overlap_threshold': 0.6,
    'top_per_module': 5,
}
PATHWAY_LIBRARY_DIR = DATA_ROOT / 'mouse_vs_human' / 'pathway_gene_sets'

SPECIES_COLORS = {'mouse': '#0072B2', 'human': '#D55E00'}
MODULE_CMAP = plt.get_cmap('tab10')
# The centroid-derived descriptor vocabulary. The label is read off the module's own centroid AFTER
# clustering; it never decides membership.
DIRECTION_DESCRIPTOR = {
    'early-declining': 'early-high, monotonic decline',
    'early-peaked': 'early-peaked',
    'mid-peaked': 'mid-PT peak',
    'mid-plateau': 'mid-PT plateau',
    'late-rising': 'monotonic rise to late',
    'late-peaked': 'late-high',
    'flat': 'flat',
}
figure_index = []


# %%
# Purpose: small plotting and numeric helpers shared by every section below.
def _save_figure(fig, filename: str):
    """Write one figure into this notebook's own output directory and record it."""
    path = OUTPUT_DIR / filename
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.show()
    plt.close(fig)
    figure_index.append(filename)
    print('figure:', filename)


def _peak_to_peak(values):
    values = np.asarray(values, dtype=float)
    return np.nanmax(values, axis=1) - np.nanmin(values, axis=1)


def _grid_thirds(grid):
    """Early / middle / late masks over the common grid: three equal blocks of grid points."""
    index = np.arange(len(grid))
    third = len(grid) / 3.0
    return index < third, (index >= third) & (index < 2 * third), index >= 2 * third


def _curve_eligibility(curves, gate, min_peak_to_peak, min_finite_fraction=0.5):
    """Which features may enter a clustering, and which grid columns that run can actually see.

    A run whose curves come from a single specimen is undefined outside that specimen's own
    pseudospace support. Requiring every grid point to be finite would then reject EVERY feature and
    an omission run would report zero modules, which is exactly how a leave-one-out probe comes back
    empty. The support is therefore the columns for which at least ``min_finite_fraction`` of the
    features are finite (the same rule ``pseudospace.modules`` standardises on), a feature must be
    finite across all of them, and the caller's absolute-magnitude gate is applied on top.
    Returns ``(eligible, column_ok)``.
    """
    values = np.asarray(curves, dtype=float)
    finite = np.isfinite(values)
    column_ok = finite.mean(axis=0) >= float(min_finite_fraction)
    if column_ok.sum() < 3:                    # too little overlap left to describe any shape
        column_ok = np.ones(values.shape[1], dtype=bool)
    eligible = (finite[:, column_ok].all(axis=1)
                & (np.asarray(gate, dtype=float) >= min_peak_to_peak))
    return eligible, column_ok


def _row_spearman(matrix, reference):
    return np.array([safe_spearman(row, reference) for row in np.asarray(matrix)], dtype=float)


def _curve_metric_table(curves_mouse, curves_human, grid, feature_names):
    """Per-gene metrics of the two fitted species curves, the raw difference and the centred response.

    ``level_shift`` is the mean difference over the grid - a global human-vs-mouse offset.
    ``response_rms`` and ``response_peak_to_peak`` are computed on the **centred** response
    ``R = D - level_shift``, so they describe pseudospace-dependent remodelling on their own. That
    separation is the point: a gene can differ hugely between species overall and not be spatially
    remodelled at all, and it is then described here rather than forced into a response module.
    """
    mouse = np.asarray(curves_mouse, dtype=float)
    human = np.asarray(curves_human, dtype=float)
    difference = human - mouse
    level_shift = np.nanmean(difference, axis=1)
    centred = difference - level_shift[:, None]
    early, middle, late = _grid_thirds(grid)

    table = pd.DataFrame({
        'gene': np.asarray(feature_names),
        'mouse_mean_expression': np.nanmean(mouse, axis=1),
        'human_mean_expression': np.nanmean(human, axis=1),
        'mouse_peak_to_peak': _peak_to_peak(mouse),
        'human_peak_to_peak': _peak_to_peak(human),
        'amplitude_change': _peak_to_peak(human) - _peak_to_peak(mouse),
        'level_shift': level_shift,
        'response_rms': np.sqrt(np.nanmean(centred ** 2, axis=1)),
        'response_peak_to_peak': _peak_to_peak(centred),
        'early_mean_difference': np.nanmean(difference[:, early], axis=1),
        'middle_mean_difference': np.nanmean(difference[:, middle], axis=1),
        'late_mean_difference': np.nanmean(difference[:, late], axis=1),
        'max_abs_difference': np.nanmax(np.abs(difference), axis=1),
        'mouse_human_curve_spearman': [
            safe_spearman(left, right) for left, right in zip(mouse, human)
        ],
    })
    # The level/amplitude/pattern split from notebook 03's Section 4: `pattern_rms_z` is invariant to
    # level and amplitude, so it moves only when peak position, width or monotonicity changes.
    effects = summarize_curve_effects(mouse, human, feature_names=np.asarray(feature_names))
    effects = effects.rename(columns={'feature': 'gene'})
    extra = ['level_fraction', 'shape_fraction', 'pattern_rms_z', 'amplitude_ratio',
             'amplitude_log2_ratio', 'pattern_status', 'difference_type']
    return table.merge(effects[['gene', *extra]], on='gene', how='left')


def _curve_pca(standardized, config):
    """PCA through the retained-component rule, applied identically wherever curves are clustered.

    ``pc_variance_target`` caps the component count at the first component reaching that cumulative
    share, and ``pc_max_components`` stops a diffuse tail being kept: dozens of tiny components would
    each contribute noise to the distance.
    """
    pca = PCA(n_components=min(standardized.shape), random_state=0).fit(standardized)
    cumulative = np.cumsum(pca.explained_variance_ratio_)
    n_pc = int(np.clip(np.searchsorted(cumulative, config['pc_variance_target']) + 1, 2,
                       min(config['pc_max_components'], standardized.shape[1])))
    scores = pca.transform(standardized)[:, :n_pc]
    if config['pc_whiten']:
        spread = scores.std(axis=0, ddof=1)
        scores = scores / np.where(spread > 0, spread, 1.0)
    return pca, scores, n_pc, cumulative


def _k_diagnostics(tree, scores, config):
    """Silhouette, smallest module and largest module fraction for every candidate k.

    One linkage serves every k, so the sweep costs one O(n^2) build rather than one per k. The
    largest-module fraction is reported as a diagnostic only: it does not select k.
    """
    rows = []
    for k in config['k_values']:
        labels = fcluster(tree, int(k), criterion='maxclust')
        sizes = np.bincount(labels)[1:]
        rows.append({
            'k': int(k),
            'silhouette': float(silhouette_score(
                scores, labels,
                sample_size=min(config['silhouette_sample_size'], labels.size),
                random_state=0)),
            'smallest_module': int(sizes.min()),
            'largest_module_fraction': float(sizes.max() / sizes.sum()),
        })
    return pd.DataFrame(rows)


def _select_k(diagnostics, config):
    """Parsimonious k: best silhouette among the ks that avoid a tiny module, then the smallest k
    within ``silhouette_tolerance`` of that best.

    The guard is "do not fragment the transcriptome into many tiny modules"; the tolerance is the
    parsimony rule, so an equally good larger split is not preferred. Pathway enrichment never enters
    this choice.
    """
    usable = diagnostics[diagnostics['smallest_module'] >= config['min_module_size']]
    if usable.empty:
        usable = diagnostics
    best = float(usable['silhouette'].max())
    within = usable[usable['silhouette'] >= best - config['silhouette_tolerance']]
    return int(within['k'].min()), best


def _module_order(raw_labels, curves, grid, prefix):
    """Rename raw cluster ids to ``prefix1..prefixK``, ordered early -> late by centroid peak.

    Module names are arbitrary; ordering them by where their centroid peaks makes the catalogs and
    the cross-tab readable without implying any ranking of importance.
    """
    peaks = {}
    for cluster in np.unique(raw_labels):
        centroid = np.nanmean(curves[raw_labels == cluster], axis=0)
        finite = np.isfinite(centroid)
        peaks[int(cluster)] = (float(grid[finite][int(np.argmax(centroid[finite]))])
                               if finite.sum() >= 3 else np.inf)
    ordered = sorted(peaks, key=lambda cluster: (peaks[cluster], -int((raw_labels == cluster).sum())))
    mapping = {cluster: f'{prefix}{position + 1}' for position, cluster in enumerate(ordered)}
    return np.array([mapping[int(cluster)] for cluster in raw_labels], dtype=object)


def _cluster_curves(curves, grid, feature_names, config, prefix, representativeness):
    """Standardise -> PCA -> whitened Ward clustering -> k diagnostics -> labels, centroids, catalog.

    ``curves`` is the unstandardised feature x grid matrix of curves that have already passed their
    absolute-magnitude gate. Membership is decided on curve **shape** only: each curve is z-scored
    over the grid first, so level and amplitude cannot decide a module.
    """
    curves = np.asarray(curves, dtype=float)
    feature_names = np.asarray(feature_names)
    standardized = zscore_rows(curves)
    pca, scores, n_pc, cumulative = _curve_pca(standardized, config)
    tree = linkage(scores, method=config['ward_method'])
    diagnostics = _k_diagnostics(tree, scores, config)
    k, best_silhouette = _select_k(diagnostics, config)
    labels = _module_order(fcluster(tree, int(k), criterion='maxclust'), curves, grid, prefix)

    rows, assignment, centroid_rows = [], [], []
    for module in sorted(set(labels), key=lambda name: int(name[len(prefix):])):
        members = labels == module
        raw_centroid = np.nanmean(curves[members], axis=0)
        std_centroid = np.nanmean(standardized[members], axis=0)
        correlations = _row_spearman(standardized[members], std_centroid)
        amplitude = _peak_to_peak(curves[members])
        # Representative genes must BOTH follow the module centroid and carry a real effect, so the
        # ranking statistic multiplies the two rather than taking either alone.
        rank = correlations * np.asarray(representativeness)[members]
        best = np.argsort(-np.nan_to_num(rank, nan=-np.inf))[:5]
        descriptors = curve_descriptors(raw_centroid[None, :], grid).iloc[0].to_dict()
        direction = descriptors['direction']
        rows.append({
            'module': module,
            'n_genes': int(members.sum()),
            'peak_position': descriptors['peak_position'],
            'half_max_width': descriptors['half_max_width'],
            'monotonicity': descriptors['monotonicity'],
            'direction': direction,
            'descriptor': DIRECTION_DESCRIPTOR.get(direction, direction),
            'median_amplitude': float(np.nanmedian(amplitude)),
            'median_member_to_centroid_correlation': float(np.nanmedian(correlations)),
            'representative_genes': '; '.join(feature_names[members][best]),
        })
        assignment.append(pd.DataFrame({
            'gene': feature_names[members],
            'module': module,
            'member_to_centroid_correlation': correlations,
            'representativeness': rank,
        }))
        for kind, values in (('standardized_centroid', std_centroid),
                            ('unscaled_centroid', raw_centroid)):
            centroid_rows.append(pd.DataFrame({
                'module': module, 'centroid': kind,
                'pseudospace': np.asarray(grid, dtype=float), 'value': values,
            }))

    return {
        'standardized': standardized,
        'pca': pca,
        'scores': scores,
        'n_pc': n_pc,
        'cumulative_variance': cumulative,
        'diagnostics': diagnostics,
        'k': int(k),
        'best_silhouette': best_silhouette,
        'labels': pd.Series(labels, index=feature_names, name='module'),
        'catalog': pd.DataFrame(rows).sort_values('module').reset_index(drop=True),
        'assignments': pd.concat(assignment, ignore_index=True),
        'centroids': pd.concat(centroid_rows, ignore_index=True),
        'loadings': pd.DataFrame({
            'pseudospace': np.asarray(grid, dtype=float),
            **{f'pc{component + 1}_loading': pca.components_[component]
               for component in range(n_pc)},
        }),
    }


def _pca_scores_frame(result, feature_names):
    frame = pd.DataFrame({'gene': np.asarray(feature_names)})
    for component in range(result['n_pc']):
        frame[f'pc{component + 1}_score'] = result['scores'][:, component]
    return frame


def _k_sweep_figure(diagnostics, title, filename):
    """Figure: explained variance plus per-k objective diagnostics, with the chosen k marked."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))
    axes[0].plot(np.arange(1, len(diagnostics) + 1), diagnostics['silhouette'], marker='o', ms=3)
    axes[0].set_xlabel('k (candidate clusters)')
    axes[0].set_ylabel('silhouette (deterministic subsample)')
    axes[0].set_title('Objective k diagnostics')
    twin = axes[0].twinx()
    twin.plot(np.arange(1, len(diagnostics) + 1),
              diagnostics['largest_module_fraction'], color='#E15759', ls='--', marker='.')
    twin.set_ylabel('largest module fraction', color='#E15759')
    twin.tick_params(axis='y', colors='#E15759')
    axes[1].plot(np.arange(1, len(diagnostics) + 1), diagnostics['smallest_module'], marker='o',
                 ms=3, color='#59A14F')
    axes[1].axhline(CLUSTER_CONFIG['min_module_size'], color='k', ls=':', lw=1)
    axes[1].set_yscale('log')
    axes[1].set_xlabel('k (candidate clusters)')
    axes[1].set_ylabel('smallest module (genes)')
    axes[1].set_title(f"Fragmentation guard (>= {CLUSTER_CONFIG['min_module_size']} genes)")
    fig.suptitle(title, fontsize=11)
    _save_figure(fig, filename)


def _format_k_diagnostics(diagnostics, k, best):
    frame = diagnostics.copy()
    frame['selected'] = frame['k'].eq(int(k))
    print(f'selected k = {k} (best silhouette {best:.3f}); '
          f'smallest module {int(frame.loc[frame["selected"], "smallest_module"].iloc[0])} genes, '
          f'largest module {frame.loc[frame["selected"], "largest_module_fraction"].iloc[0]:.1%}')
    return frame


print('Setup complete.')


# %% [markdown]
# ## 1 - Input and PT-only scope
#
# Start from the PT-specific artifact notebook 03 wrote. It already contains reviewed PT structures
# only, both species' comparison samples, the PT S1-to-S3 DPT coordinate, the accepted ortholog
# mapping, log-normalised expression and the specimen identifiers. Everything this notebook needs is
# verified here and an absence fails loudly rather than propagating a silently wrong analysis.
#
# Nothing in this notebook re-does integration, cell typing, DPT, the PT trajectory or the ortholog
# map, and no nephron segment outside PT enters the object.

# %%
# Purpose: load the PT-only artifact and verify the exact contract this notebook relies on.
from pseudospace.vocabulary import KEEP_TUBULE_CLASSES

for required_path in (DPT_OUTPUT_PATH, ORTHOLOG_MAP_PATH):
    if not required_path.exists():
        raise FileNotFoundError(
            f'{required_path} is missing. Run analysis/notebooks/03_human_vs_healthy_mouse.ipynb '
            'first: this notebook consumes its PT pseudospace coordinate and accepted ortholog map.'
        )

adata_full = sc.read_h5ad(DPT_OUTPUT_PATH)
if UNIVERSE_CONFIG['measurement'] not in adata_full.layers:
    raise KeyError(f'{DPT_OUTPUT_PATH.name} is missing layers["{UNIVERSE_CONFIG["measurement"]}"].')
adata_full.X = adata_full.layers[UNIVERSE_CONFIG['measurement']].copy()

missing_obs = [column for column in ('comparison_species', 'sample', 'total_scanpy_dpt')
               if column not in adata_full.obs.columns]
if missing_obs:
    raise KeyError(f'{DPT_OUTPUT_PATH.name} is missing obs columns: {missing_obs}')
if 'measured_in_both_inputs' not in adata_full.var.columns:
    raise KeyError(f'{DPT_OUTPUT_PATH.name} is missing var["measured_in_both_inputs"].')

# The artifact is documented as PT-only. Verify that where it is self-describing, and fall back to the
# reviewed class column for the PT restriction so an older object cannot smuggle in other segments.
pt_restriction = np.ones(adata_full.n_obs, dtype=bool)
if 'broad_tubule_marker_call' in adata_full.obs.columns:
    broad = adata_full.obs['broad_tubule_marker_call'].astype(str)
    unexpected = sorted(set(broad) - set(KEEP_TUBULE_CLASSES))
    if unexpected:
        raise ValueError(f'{DPT_OUTPUT_PATH.name} contains non-nephron classes: {unexpected}')
    non_pt = sorted(set(broad) - {'PT'})
    if non_pt:
        print(f'NOTE: the artifact is not PT-only (also carries {non_pt}); restricting to PT here.')
        pt_restriction = broad.eq('PT').to_numpy()
temp_mask = (pt_restriction
             & adata_full.obs['comparison_species'].astype(str).isin(SPECIES_COLORS)
             & np.isfinite(adata_full.obs['total_scanpy_dpt'].to_numpy(dtype=float))).to_numpy()
adata_pt = adata_full[temp_mask].copy()
del temp_mask

if adata_pt.n_obs < 100 or adata_pt.obs['comparison_species'].nunique() != 2:
    raise ValueError('Insufficient PT structures from both species.')

s = adata_pt.obs['total_scanpy_dpt'].to_numpy(dtype=float)
species = adata_pt.obs['comparison_species'].astype(str).to_numpy()
comparison_indicator = (species == 'human').astype(float)
samples = adata_pt.obs['sample'].astype(str).to_numpy()

# Species comes from the data, not from a hardcoded sample list, so a renamed section cannot silently
# be treated as the wrong species.
sample_species = (adata_pt.obs[['sample', 'comparison_species']].astype(str)
                  .drop_duplicates().set_index('sample')['comparison_species'].to_dict())
species_by_sample = dict(sample_species)
if len(set(species_by_sample.values())) != 2:
    raise ValueError(f'Samples do not map onto two species: {species_by_sample}')
for name in sorted(species_by_sample):
    observed = set(species[samples == name])
    if observed != {species_by_sample[name]}:
        raise ValueError(f'Sample {name} carries more than one species: {observed}')
mouse_specimens = sorted(name for name, kind in species_by_sample.items() if kind == 'mouse')
human_specimens = sorted(name for name, kind in species_by_sample.items() if kind == 'human')
if not mouse_specimens or not human_specimens:
    raise ValueError('Both species need at least one specimen.')

pt_cohort_summary = (
    adata_pt.obs.groupby(['comparison_species', 'sample'], observed=True)
    .agg(n_pt_structures=('sample', 'size'),
         dpt_min=('total_scanpy_dpt', 'min'),
         dpt_max=('total_scanpy_dpt', 'max'))
    .reset_index()
)
pt_cohort_summary.to_csv(DIAGNOSTIC_DIR / 'pt_cohort_summary.csv', index=False)
display(pt_cohort_summary.round(3))
print(f'PT structures: {adata_pt.n_obs:,}; genes assayed: {adata_pt.n_vars:,}')
print(f'Mouse specimens: {mouse_specimens}; human sections: {human_specimens} '
      '(two sections of ONE donor, biologically one replicate)')


# %% [markdown]
# ## 2 - Gene universe
#
# Filtering happens **here**, after loading the PT-only object, and it is deliberately minimal.
#
# A gene is eligible when it was measured in **both** species and is detected in at least 2% of the PT
# structures in this notebook. `measured_in_both_inputs` is the availability rule from notebook 03:
# the accepted ortholog map creates a target column for every pair, so a gene absent from one input's
# feature list is a structural zero there and must not enter a fit or an enrichment background.
#
# Mitochondrial genes, ribosomal genes and PT-axis genes are **not** excluded, and no flag is created
# for any of them. Nothing is filtered on differential-expression significance, a GAM p-value, a
# human-mouse effect size or pathway membership.

# %%
# Purpose: the gene universe, its audit, and the matrix every fit below shares.
Y_all = as_csr(adata_pt.layers[UNIVERSE_CONFIG['measurement']])
detected_structures = np.asarray((Y_all > 0).sum(axis=0)).ravel()
gene_mean_all = np.asarray(Y_all.mean(axis=0)).ravel()
min_detected = int(np.ceil(UNIVERSE_CONFIG['min_detected_fraction'] * adata_pt.n_obs))
measured_in_both = adata_pt.var['measured_in_both_inputs'].to_numpy(dtype=bool)

eligible_genes = measured_in_both & (detected_structures >= min_detected)
gene_names = adata_pt.var_names.to_numpy()[eligible_genes]
Y_genes = Y_all[:, eligible_genes].tocsr().astype(np.float64)
if gene_names.size < 20:
    raise ValueError(f'Only {gene_names.size} genes are eligible; the PT object looks wrong.')
Y_GENES_FINGERPRINT = digest(Y_genes)

ortholog_map = pd.read_csv(ORTHOLOG_MAP_PATH)
gene_universe_audit = pd.DataFrame([
    {'stage': 'accepted ortholog universe (assayed genes in the PT object)', 'n_genes': int(adata_pt.n_vars)},
    {'stage': 'accepted ortholog pairs in ortholog_map_used.csv', 'n_genes': int(len(ortholog_map))},
    {'stage': 'measured in both inputs', 'n_genes': int(measured_in_both.sum())},
    {'stage': f'detected in >= {UNIVERSE_CONFIG["min_detected_fraction"]:.0%} of PT structures '
              f'(>= {min_detected} structures)', 'n_genes': int((detected_structures >= min_detected).sum())},
    {'stage': 'eligible: measured in both AND detected', 'n_genes': int(gene_names.size)},
    {'stage': 'excluded: not measured in both inputs', 'n_genes': int((~measured_in_both).sum())},
    {'stage': 'excluded: below the detection floor', 'n_genes': int((detected_structures < min_detected).sum())},
])
gene_universe_audit['fraction_of_assayed'] = gene_universe_audit['n_genes'] / float(adata_pt.n_vars)
gene_universe_audit.to_csv(OUTPUT_DIR / 'gene_universe_audit.csv', index=False)
display(gene_universe_audit.round(4))
print(f'Eligible genes: {gene_names.size:,} of {adata_pt.n_vars:,} assayed '
      f'({int(ortholog_map["mouse_symbol"].isin(adata_pt.var_names).sum()):,} of '
      f'{len(ortholog_map):,} accepted ortholog pairs land in this assay)')
print('No mitochondrial, ribosomal or PT-axis gene was excluded, and no such flag was created.')


# %%
# Purpose: figure 1 - the cohort and the gene-universe audit.
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2),
                         gridspec_kw={'width_ratios': [1, 1.25]})
cohort = pt_cohort_summary.copy()
cohort['species_color'] = cohort['comparison_species'].map(SPECIES_COLORS)
axes[0].bar(cohort['sample'], cohort['n_pt_structures'], color=cohort['species_color'])
for position, row in enumerate(cohort.itertuples()):
    axes[0].text(position, row.n_pt_structures, f'{row.n_pt_structures:,}', ha='center', va='bottom',
                 fontsize=8)
axes[0].set_ylabel('PT structures')
axes[0].set_title(f'PT cohort ({adata_pt.n_obs:,} structures)')
axes[0].tick_params(axis='x', rotation=20)
handles = [plt.Rectangle((0, 0), 1, 1, color=colour) for colour in SPECIES_COLORS.values()]
axes[0].legend(handles, [f'{name} specimen' if name == 'mouse' else f'{name} sections'
                         for name in SPECIES_COLORS], frameon=False, fontsize=8)
axes[0].text(0.02, 0.02, 'human = 2 sections, 1 donor', transform=axes[0].transAxes, fontsize=8,
             style='italic')

funnel = gene_universe_audit.head(5)
axes[1].barh(funnel['stage'][::-1], funnel['n_genes'][::-1], color='#6B7A8F')
for position, (stage, count) in enumerate(zip(funnel['stage'][::-1], funnel['n_genes'][::-1])):
    axes[1].text(count, position, f' {count:,}', va='center', fontsize=8)
axes[1].set_xlabel('genes')
axes[1].set_title('Gene universe')
axes[1].tick_params(axis='y', labelsize=8)
fig.suptitle('Figure 1 - cohort and gene-universe audit', fontsize=12)
_save_figure(fig, 'fig01_cohort_and_gene_universe.png')


# %% [markdown]
# ## 3 - Common PT support and GAM fitting
#
# The common mouse-human support is the overlap of the two species' 1st-99th percentile DPT ranges,
# carried on a 101-point grid, with the same spline basis, knot count and smoothing convention
# notebook 03 validated.
#
# For every eligible gene this notebook fits a single-condition curve **per specimen** and then builds
# the specimen-balanced species curves:
#
# * `M_g(s)` - equally weighted mean of the mouse specimen curves,
# * `H_g(s)` - equally weighted mean of the human section curves.
#
# These balanced curves are the PRIMARY curves here. The per-feature smoothing parameter is the GCV
# choice from the pooled nested fit, reused so every specimen curve sits on the same smoothing scale
# as the pooled curves it is drawn against. Pooled-specimen curves are also kept, for the sensitivity
# analysis in section 10 only.

# %%
# Purpose: common support, grid, knots, and the pooled nested fit whose GCV lambda we reuse.
GAM_LAMBDA_GRID = GAM_CONFIG['lambda_grid']
SUPPORT_PERCENTILES = GAM_CONFIG['common_support_pct']
p_lo, p_hi = SUPPORT_PERCENTILES
lo = max(np.percentile(s[species == name], p_lo) for name in SPECIES_COLORS)
hi = min(np.percentile(s[species == name], p_hi) for name in SPECIES_COLORS)
if not lo < hi:
    raise ValueError(f'Mouse and human PT have no common DPT support: {lo:.4f} >= {hi:.4f}.')
grid = np.linspace(lo, hi, GAM_CONFIG['grid_points'])
knots = gam_internal_knots(s, basis_df=3 + GAM_CONFIG['n_internal_knots'])

# The pooled fit is only needed for its per-gene GCV lambda. The stage cache key assumes 03's exact
# conventions, so with them unchanged this is a cache hit and not a re-fit.
pooled_fit = cached_run_level_shape(
    Y_genes, s, comparison_indicator, knots, grid, GAM_LAMBDA_GRID,
    stage='section4_gene_fit', root=STAGE_CACHE_DIR, y_fingerprint=Y_GENES_FINGERPRINT,
    enabled=STAGE_CACHE_ENABLED,
)
gene_lambda_index = np.asarray(pooled_fit['lam_idx'], dtype=int)
print(f'Common mouse/human DPT support: [{lo:.4f}, {hi:.4f}]; grid {grid.size} points; '
      f'{knots.size} internal knots')
print(f'GCV lambda indices reused from the pooled fit: {np.unique(gene_lambda_index).size} distinct')


# %%
# Purpose: per-specimen GAM fits and the specimen-balanced species curves (cached, not refit).
def _fit_specimen_curves():
    """Fit one single-condition GAM per specimen for every eligible gene."""
    fitted = {}
    for name in sorted(species_by_sample):
        mask = samples == name
        fitted[name], _ = fit_single_condition_curves(
            Y_genes[mask], s[mask], knots, grid, GAM_LAMBDA_GRID, gene_lambda_index,
            support_pct=SUPPORT_PERCENTILES,
        )
    return fitted


specimen_curves = cached_payload(
    'pt_specimen_curves', _fit_specimen_curves, root=STAGE_CACHE_DIR,
    params={'grid': np.asarray(grid), 'lambda_grid': np.asarray(GAM_LAMBDA_GRID),
            'knots': np.asarray(knots), 'support_pct': list(SUPPORT_PERCENTILES)},
    inputs={'s': np.asarray(s), 'samples': np.asarray(samples), 'y': Y_GENES_FINGERPRINT,
            'lambda_index': gene_lambda_index},
    code=code_digest(fit_single_condition_curves), enabled=STAGE_CACHE_ENABLED,
)

mouse_curves = specimen_balanced_curves({name: specimen_curves[name] for name in mouse_specimens})
human_curves = specimen_balanced_curves({name: specimen_curves[name] for name in human_specimens})
# Pooled-specimen curves, kept for section 10A only; never the primary curves.
pooled_mouse_curves = np.asarray(pooled_fit['curve_healthy'], dtype=float)
pooled_human_curves = np.asarray(pooled_fit['curve_aki'], dtype=float)

if not np.isfinite(mouse_curves).all() or not np.isfinite(human_curves).all():
    raise ValueError('The balanced species curves are not finite over the common grid.')

print(f'Per-specimen fits: {len(specimen_curves)} specimens '
      f'({len(mouse_specimens)} mouse, {len(human_specimens)} human sections)')
print(f'Balanced curve matrices: mouse {mouse_curves.shape}, human {human_curves.shape}')
print('Balanced = equal weight per specimen. It stops one section dominating on structure count; it '
      'does NOT create independent human biological replication.')


# %% [markdown]
# ## 4 - Gene-level curve metrics
#
# For every eligible gene, with `D_g(s) = H_g(s) - M_g(s)`:
#
# ```
# level_shift_g = mean_s(D_g(s))          a global human-vs-mouse expression difference
# R_g(s)        = D_g(s) - level_shift_g  pseudospace-dependent remodelling
# ```
#
# This separation is essential and is the reason the primary cross-species analysis in section 6
# clusters `R_g(s)` and **never** the raw `D_g(s)`. A gene whose human and mouse curves differ by a
# constant offset has a large `level_shift` and no remodelling; clustering the raw difference would
# sort genes by how big their offset is and call that a program.
#
# These metrics describe the modules found later. They must not predefine biological response classes,
# and they are not used to select the genes that get clustered.

# %%
# Purpose: the gene-level metrics table, the primary curves, and the centred response matrix.
gene_curve_metrics = _curve_metric_table(mouse_curves, human_curves, grid, gene_names)
gene_curve_metrics.to_csv(OUTPUT_DIR / 'gene_curve_metrics.csv', index=False)

response_difference = difference_curves(mouse_curves, human_curves)               # D
response_centred = difference_curves(mouse_curves, human_curves, center=True)     # R = D - mean(D)

# `_curve_metric_table` returns one row per gene in `gene_names` order, and its effect merge is a left
# join that preserves that order, so these arrays are positionally aligned with `gene_names` and with
# the curve matrices. Deliberately no label-indexed reindex: a duplicated var name makes one raise.
gene_level_shift = gene_curve_metrics['level_shift'].to_numpy(dtype=float)
gene_response_rms = gene_curve_metrics['response_rms'].to_numpy(dtype=float)
gene_response_peak_to_peak = gene_curve_metrics['response_peak_to_peak'].to_numpy(dtype=float)
gene_mouse_amplitude = gene_curve_metrics['mouse_peak_to_peak'].to_numpy(dtype=float)
gene_human_amplitude = gene_curve_metrics['human_peak_to_peak'].to_numpy(dtype=float)
gene_curve_spearman = gene_curve_metrics['mouse_human_curve_spearman'].to_numpy(dtype=float)
if not np.allclose(gene_level_shift, np.nanmean(response_difference, axis=1)):
    raise ValueError('level_shift does not match the mean of the fitted difference curve.')

display(gene_curve_metrics.head(8)[[
    'gene', 'mouse_peak_to_peak', 'human_peak_to_peak', 'level_shift', 'response_rms',
    'response_peak_to_peak', 'mouse_human_curve_spearman', 'difference_type', 'pattern_status',
]].round(4))
print(f'Metrics computed for {len(gene_curve_metrics):,} genes.')
print('difference_type counts (a reading aid, not a gate): '
      f'{gene_curve_metrics["difference_type"].value_counts().to_dict()}')
print('pattern_status counts: '
      f'{gene_curve_metrics["pattern_status"].value_counts().to_dict()}')


# %% [markdown]
# ## 5 - Analysis A: normal mouse PT positional programs
#
# > Which recurring expression shapes organise the normal PT transcriptome?
#
# Only `M_g(s)` enters this section. Human expression and pathway information are not consulted while
# these modules are being discovered.
#
# The absolute amplitude gate comes first. A nearly flat curve must not be standardised and clustered:
# dividing it by its own tiny spread promotes numerical noise into an apparently strong shape. Genes
# that fail the gate are excluded from this clustering and reported - they are not silently absorbed.

# %%
# Purpose: the positional amplitude gate, with its distribution and exclusion count.
POSITIONAL_AMPLITUDE_FLOOR = UNIVERSE_CONFIG['positional_min_peak_to_peak']
positional_eligible, positional_columns = _curve_eligibility(
    mouse_curves, gene_mouse_amplitude, POSITIONAL_AMPLITUDE_FLOOR)
if positional_eligible.sum() < 20:
    raise ValueError('Too few genes pass the positional amplitude gate to cluster.')

fig, axes = plt.subplots(1, 2, figsize=(11, 3.6))
axes[0].hist(gene_mouse_amplitude, bins=80, color='#0072B2')
axes[0].axvline(POSITIONAL_AMPLITUDE_FLOOR, color='#D55E00', ls='--')
axes[0].set_yscale('log')
axes[0].set_xlabel('mouse curve peak-to-peak (log-normalised units)')
axes[0].set_ylabel('genes')
axes[0].set_title(f'Mouse amplitude; gate at {POSITIONAL_AMPLITUDE_FLOOR}')
axes[1].hist(np.log10(np.maximum(gene_human_amplitude - gene_mouse_amplitude, 1e-6)), bins=80,
             color='#6B7A8F')
axes[1].set_xlabel('log10(human amplitude - mouse amplitude)')
axes[1].set_ylabel('genes')
axes[1].set_title('Amplitude change (descriptive)')
fig.suptitle('Figure 2a - positional amplitude gate', fontsize=12)
_save_figure(fig, 'fig02a_positional_amplitude_gate.png')

print(f'Positional-clustering eligible: {int(positional_eligible.sum()):,} genes; '
      f'{int((~positional_eligible).sum()):,} excluded as essentially flat '
      f'(peak-to-peak < {POSITIONAL_AMPLITUDE_FLOOR}).')
print(f'Support this clustering was run on: {int(positional_columns.sum())} of {grid.size} '
      'grid points.')
print('The threshold is a named parameter, never tuned on pathway results.')


# %%
# Purpose: standardise, run PCA, and cluster the eligible mouse curves.
positional = _cluster_curves(
    mouse_curves[positional_eligible][:, positional_columns], grid[positional_columns],
    gene_names[positional_eligible], CLUSTER_CONFIG, 'P',
    representativeness=gene_mouse_amplitude[positional_eligible],
)
positional_labels = np.full(gene_names.size, 'unassigned', dtype=object)
positional_labels[positional_eligible] = positional['labels'].to_numpy()
positional_labels = pd.Series(positional_labels, index=gene_names, name='positional_module')

_pca_scores_frame(positional, gene_names[positional_eligible]).to_csv(
    OUTPUT_DIR / 'mouse_positional_pca_scores.csv', index=False)
positional['loadings'].to_csv(OUTPUT_DIR / 'mouse_positional_pca_loadings.csv', index=False)
positional['assignments'].to_csv(
    OUTPUT_DIR / 'mouse_positional_module_assignments.csv', index=False)
positional['centroids'].to_csv(OUTPUT_DIR / 'mouse_positional_module_centroids.csv', index=False)
_format_k_diagnostics(positional['diagnostics'], positional['k'], positional['best_silhouette'])[
    ['k', 'silhouette', 'smallest_module', 'largest_module_fraction', 'selected']
].to_csv(DIAGNOSTIC_DIR / 'mouse_positional_k_diagnostics.csv', index=False)

print(f'Retained PCs: {positional["n_pc"]} '
      f'(cumulative variance {positional["cumulative_variance"][positional["n_pc"] - 1]:.3f})')
print(f'Positional modules: {positional["k"]} over {int(positional_eligible.sum()):,} genes; '
      f'sizes {positional["catalog"]["n_genes"].sort_values(ascending=False).tolist()}')


# %%
# Purpose: figure 2 - positional PCA: variance explained and the leading loading curves.
cumulative = positional['cumulative_variance']
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4))
components = np.arange(1, min(len(cumulative), CLUSTER_CONFIG['pc_max_components']) + 1)
axes[0].bar(components, cumulative[:len(components)] - np.concatenate([[0], cumulative[:len(components) - 1]]),
            color='#6B7A8F', label='per component')
axes[0].plot(components, cumulative[:len(components)], color='#D55E00', marker='o', ms=4,
             label='cumulative')
axes[0].axvline(positional['n_pc'], color='k', ls=':', lw=1)
axes[0].annotate(f'retained: {positional["n_pc"]} PCs', (positional['n_pc'], 0.5),
                 xytext=(6, 0), textcoords='offset points', fontsize=8)
axes[0].set_xlabel('principal component')
axes[0].set_ylabel('variance explained')
axes[0].set_title('Mouse positional PCA variance')
axes[0].legend(frameon=False, fontsize=8)

loadings = positional['loadings']
for component in range(positional['n_pc']):
    axes[1].plot(loadings['pseudospace'], loadings[f'pc{component + 1}_loading'],
                 label=f'PC{component + 1}')
axes[1].axhline(0, color='k', lw=0.6)
axes[1].set_xlabel('shared PT DPT')
axes[1].set_ylabel('loading (unit-norm)')
axes[1].set_title('Leading loading curves = major modes of normal PT expression')
axes[1].legend(frameon=False, fontsize=8)
fig.suptitle('Figure 2 - normal mouse positional PCA', fontsize=12)
_save_figure(fig, 'fig02_positional_pca.png')

loading_pc = [f'pc{component + 1}_score' for component in range(positional['n_pc'])]
strongest = pd.DataFrame({
    'pc': loading_pc,
    'top_positive_genes': ['; '.join(np.asarray(gene_names)[positional_eligible][
        np.argsort(-positional['scores'][:, component])[:10]]) for component in range(positional['n_pc'])],
    'top_negative_genes': ['; '.join(np.asarray(gene_names)[positional_eligible][
        np.argsort(positional['scores'][:, component])[:10]]) for component in range(positional['n_pc'])],
})
strongest.to_csv(DIAGNOSTIC_DIR / 'mouse_positional_pc_extremes.csv', index=False)
display(strongest)


# %%
# Purpose: figure 3 - all positional-module centroids on one compact panel.
_raw_centroids = positional['centroids'][positional['centroids']['centroid'].eq('unscaled_centroid')]
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4), sharex=True)
for module, block in _raw_centroids.groupby('module'):
    colour = MODULE_CMAP(int(module[1:]) - 1)
    row = positional['catalog'].set_index('module').loc[module]
    axes[0].plot(block['pseudospace'], block['value'], color=colour, lw=2,
                 label=f"{module} ({row['descriptor']}, n={int(row['n_genes'])})")
    standard = positional['centroids'][
        positional['centroids']['module'].eq(module)
        & positional['centroids']['centroid'].eq('standardized_centroid')]
    axes[1].plot(standard['pseudospace'], standard['value'], color=colour, lw=2, label=module)
axes[0].set_xlabel('shared PT DPT')
axes[0].set_ylabel('fitted log-normalised expression')
axes[0].set_title('Unscaled mouse centroid per module')
axes[0].legend(frameon=False, fontsize=7)
axes[1].axhline(0, color='k', lw=0.6)
axes[1].set_xlabel('shared PT DPT')
axes[1].set_ylabel('z-scored over pseudospace')
axes[1].set_title('Standardised centroid (the shape that was clustered)')
fig.suptitle('Figure 3 - positional module centroids', fontsize=12)
_save_figure(fig, 'fig03_positional_module_centroids.png')

display(positional['catalog'][[
    'module', 'n_genes', 'descriptor', 'peak_position', 'monotonicity', 'median_amplitude',
    'median_member_to_centroid_correlation',
]].round(3))


# %% [markdown]
# ## 6 - Analysis B: human-vs-mouse remodelling programs
#
# This is the PRIMARY cross-species analysis.
#
# > Which recurring spatial patterns describe how gene expression along PT is remodelled in human
# > relative to healthy mouse?
#
# The clustering runs on the centred response `R_g(s) = D_g(s) - mean(D_g)`, so it captures positional
# remodelling independently of a constant global level shift. The raw human and mouse curves are
# **not** clustered jointly, and separate mouse and human clusterings are not used as the main species
# comparison.
#
# A gene whose cross-species difference is essentially constant or negligible must not be standardised
# and clustered: scaling a flat residual curve would manufacture a shape out of noise. The
# absolute-magnitude gate on `peak_to_peak(R_g)` is applied BEFORE scaling, and it is not a p-value.

# %%
# Purpose: the response gate, its distributions, and the exclusion accounting.
RESPONSE_PEAK_TO_PEAK_FLOOR = UNIVERSE_CONFIG['response_min_peak_to_peak']
response_eligible, response_columns = _curve_eligibility(
    response_centred, gene_response_peak_to_peak, RESPONSE_PEAK_TO_PEAK_FLOOR)
if response_eligible.sum() < 20:
    raise ValueError('Too few genes pass the response gate to cluster.')

fig, axes = plt.subplots(1, 3, figsize=(13, 3.5))
axes[0].hist(gene_response_peak_to_peak, bins=80, color='#D55E00')
axes[0].axvline(RESPONSE_PEAK_TO_PEAK_FLOOR, color='k', ls='--')
axes[0].set_xlabel('centred response peak-to-peak')
axes[0].set_ylabel('genes')
axes[0].set_title(f'Response amplitude; gate at {RESPONSE_PEAK_TO_PEAK_FLOOR}')
axes[1].hist(gene_response_rms, bins=80, color='#E15759')
axes[1].set_xlabel('centred response RMS')
axes[1].set_ylabel('genes')
axes[1].set_title('Centred response RMS')
axes[2].hist(gene_level_shift, bins=80, color='#6B7A8F')
axes[2].axvline(0, color='k', lw=0.8)
axes[2].set_xlabel('global human-minus-mouse level shift')
axes[2].set_ylabel('genes')
axes[2].set_title('Global level shift (a separate effect)')
fig.suptitle('Figure 4a - response gate distributions', fontsize=12)
_save_figure(fig, 'fig04a_response_gate_distributions.png')

print(f'Total eligible genes: {gene_names.size:,}')
print(f'Response-clustering eligible: {int(response_eligible.sum()):,}; support '
      f'{int(response_columns.sum())} of {grid.size} grid points')
print(f'Excluded because the cross-species difference is essentially constant or negligible: '
      f'{int((~response_eligible).sum()):,} '
      f'(peak-to-peak(R) < {RESPONSE_PEAK_TO_PEAK_FLOOR})')
print('A gene can carry a large global level shift and still fail the spatial-remodelling criterion. '
      'That is expected: it stays in gene_curve_metrics and is kept out of the response modules.')


# %%
# Purpose: standardise, run PCA, and cluster the eligible centred response curves.
response = _cluster_curves(
    response_centred[response_eligible][:, response_columns], grid[response_columns],
    gene_names[response_eligible], CLUSTER_CONFIG, 'R',
    representativeness=gene_response_rms[response_eligible],
)
response_labels = np.full(gene_names.size, 'unassigned', dtype=object)
response_labels[response_eligible] = response['labels'].to_numpy()
response_labels = pd.Series(response_labels, index=gene_names, name='response_module')

_pca_scores_frame(response, gene_names[response_eligible]).to_csv(
    OUTPUT_DIR / 'response_pca_scores.csv', index=False)
response['loadings'].to_csv(OUTPUT_DIR / 'response_pca_loadings.csv', index=False)
response['assignments'].to_csv(OUTPUT_DIR / 'response_module_assignments.csv', index=False)
response['centroids'].to_csv(OUTPUT_DIR / 'response_module_centroids.csv', index=False)
_format_k_diagnostics(response['diagnostics'], response['k'], response['best_silhouette'])[
    ['k', 'silhouette', 'smallest_module', 'largest_module_fraction', 'selected']
].to_csv(DIAGNOSTIC_DIR / 'response_k_diagnostics.csv', index=False)

print(f'Retained PCs: {response["n_pc"]} '
      f'(cumulative variance {response["cumulative_variance"][response["n_pc"] - 1]:.3f})')
print(f'Response modules: {response["k"]} over {int(response_eligible.sum()):,} genes; '
      f'sizes {response["catalog"]["n_genes"].sort_values(ascending=False).tolist()}')


# %%
# Purpose: figure 4 - response PCA: variance explained and the leading loading curves.
cumulative_r = response['cumulative_variance']
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4))
components_r = np.arange(1, min(len(cumulative_r), CLUSTER_CONFIG['pc_max_components']) + 1)
axes[0].bar(components_r,
            cumulative_r[:len(components_r)] - np.concatenate([[0], cumulative_r[:len(components_r) - 1]]),
            color='#6B7A8F', label='per component')
axes[0].plot(components_r, cumulative_r[:len(components_r)], color='#D55E00', marker='o', ms=4,
             label='cumulative')
axes[0].axvline(response['n_pc'], color='k', ls=':', lw=1)
axes[0].annotate(f'retained: {response["n_pc"]} PCs', (response['n_pc'], 0.5), xytext=(6, 0),
                 textcoords='offset points', fontsize=8)
axes[0].set_xlabel('principal component')
axes[0].set_ylabel('variance explained')
axes[0].set_title('Response PCA variance')
axes[0].legend(frameon=False, fontsize=8)

response_loadings = response['loadings']
for component in range(response['n_pc']):
    axes[1].plot(response_loadings['pseudospace'], response_loadings[f'pc{component + 1}_loading'],
                 label=f'PC{component + 1}')
axes[1].axhline(0, color='k', lw=0.6)
axes[1].set_xlabel('shared PT DPT')
axes[1].set_ylabel('loading (unit-norm)')
axes[1].set_title('Loading curves = continuous modes of cross-species remodelling')
axes[1].legend(frameon=False, fontsize=8)
fig.suptitle('Figure 4 - human-vs-mouse response PCA', fontsize=12)
_save_figure(fig, 'fig04_response_pca.png')

response_pc_extremes = pd.DataFrame({
    'pc': [f'pc{component + 1}_score' for component in range(response['n_pc'])],
    'top_positive_genes': ['; '.join(np.asarray(gene_names)[response_eligible][
        np.argsort(-response['scores'][:, component])[:10]]) for component in range(response['n_pc'])],
    'top_negative_genes': ['; '.join(np.asarray(gene_names)[response_eligible][
        np.argsort(response['scores'][:, component])[:10]]) for component in range(response['n_pc'])],
})
response_pc_extremes.to_csv(DIAGNOSTIC_DIR / 'response_pc_extremes.csv', index=False)
display(response_pc_extremes)
print('These components are scientific output, not preprocessing: they are the continuous axes of '
      'remodelling. No biological label is assigned before its loading curve has been inspected.')


# %%
# Purpose: the response-module statistics table, its centroids, and representative genes.
response_module_index = gene_names[response_eligible]
response_module_labels = response['labels']
response_catalog_stats = []
for module in sorted(set(response_module_labels), key=lambda name: int(name[1:])):
    members = (response_module_labels == module).to_numpy()
    response_catalog_stats.append({
        'module': module,
        'n_genes': int(members.sum()),
        'median_level_shift': float(np.nanmedian(gene_level_shift[response_eligible][members])),
        'median_response_rms': float(np.nanmedian(gene_response_rms[response_eligible][members])),
        'median_response_peak_to_peak': float(
            np.nanmedian(gene_response_peak_to_peak[response_eligible][members])),
        'median_early_difference': float(np.nanmedian(
            response_difference[response_eligible][members][:, _grid_thirds(grid)[0]])),
        'median_middle_difference': float(np.nanmedian(
            response_difference[response_eligible][members][:, _grid_thirds(grid)[1]])),
        'median_late_difference': float(np.nanmedian(
            response_difference[response_eligible][members][:, _grid_thirds(grid)[2]])),
        'median_mouse_amplitude': float(np.nanmedian(gene_mouse_amplitude[response_eligible][members])),
        'median_human_amplitude': float(np.nanmedian(gene_human_amplitude[response_eligible][members])),
        'median_mouse_human_curve_spearman': float(
            np.nanmedian(gene_curve_spearman[response_eligible][members])),
        # Read the assignments frame by MODULE, not positionally: it is written module block by module
        # block, so it is not aligned with the feature-ordered curve matrices.
        'median_member_to_centroid_correlation': float(
            response['assignments'].loc[response['assignments']['module'].eq(module),
                                        'member_to_centroid_correlation'].median()),
        'centroid_spearman_human_vs_mouse': safe_spearman(
            np.nanmean(mouse_curves[response_eligible][members], axis=0),
            np.nanmean(human_curves[response_eligible][members], axis=0)),
    })
response_catalog_stats = pd.DataFrame(response_catalog_stats).merge(
    response['catalog'][['module', 'descriptor', 'direction', 'peak_position', 'half_max_width',
                          'monotonicity', 'representative_genes']],
    on='module', how='left',
)

# The observable mouse and human centroids are also written to `response_module_centroids.csv`, so a
# centred difference is never interpreted on its own: it cannot say whether human rose, mouse fell, or
# both moved differently.

representatives = response['assignments'].copy()
# Gene-keyed dictionaries, not positional arrays: `assignments` is module-block ordered.
representatives['response_rms'] = representatives['gene'].map(
    dict(zip(response_module_index, gene_response_rms[response_eligible])))
representatives['level_shift'] = representatives['gene'].map(
    dict(zip(response_module_index, gene_level_shift[response_eligible])))
representatives = representatives.sort_values(
    ['module', 'representativeness'], ascending=[True, False])
representatives.to_csv(OUTPUT_DIR / 'response_module_representative_genes.csv', index=False)
representatives = representatives.groupby('module', sort=False).head(5)

response_centroid_rows = []
for module in response_catalog_stats['module']:
    members = (response_module_labels == module).to_numpy()
    for kind, values in (
        ('uncentred_human_minus_mouse_centroid', np.nanmean(response_difference[response_eligible][members], axis=0)),
        ('centred_response_centroid', np.nanmean(response_centred[response_eligible][members], axis=0)),
        ('mouse_centroid', np.nanmean(mouse_curves[response_eligible][members], axis=0)),
        ('human_centroid', np.nanmean(human_curves[response_eligible][members], axis=0)),
    ):
        response_centroid_rows.append(pd.DataFrame({
            'module': module, 'centroid': kind,
            'pseudospace': np.asarray(grid, dtype=float), 'value': values,
        }))
response_module_centroids = pd.concat(response_centroid_rows, ignore_index=True)
response_module_centroids.to_csv(OUTPUT_DIR / 'response_module_centroids.csv', index=False)

display(response_catalog_stats[[
    'module', 'n_genes', 'descriptor', 'median_level_shift', 'median_response_rms',
    'median_early_difference', 'median_middle_difference', 'median_late_difference',
    'median_mouse_human_curve_spearman', 'centroid_spearman_human_vs_mouse',
]].round(3))
print(f'Response-module statistics computed for {len(response_catalog_stats)} modules.')
print('Representative genes rank by correlation to the module centroid TIMES the response RMS, so a '
      'gene must both follow the pattern and carry a real effect. Global level shift alone never '
      'selects a representative.')


# %%
# Purpose: figure 5 - response module centroids, centred and uncentred.
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4), sharex=True)
for module in response_catalog_stats['module']:
    colour = MODULE_CMAP(int(module[1:]) - 1)
    row = response_catalog_stats.set_index('module').loc[module]
    centred = response_module_centroids[
        response_module_centroids['module'].eq(module)
        & response_module_centroids['centroid'].eq('centred_response_centroid')]
    uncentred = response_module_centroids[
        response_module_centroids['module'].eq(module)
        & response_module_centroids['centroid'].eq('uncentred_human_minus_mouse_centroid')]
    axes[0].plot(centred['pseudospace'], centred['value'], color=colour, lw=2,
                 label=f"{module} ({row['descriptor']}, n={int(row['n_genes'])})")
    axes[1].plot(uncentred['pseudospace'], uncentred['value'], color=colour, lw=2, label=module)
axes[0].axhline(0, color='k', lw=0.6)
axes[0].set_xlabel('shared PT DPT')
axes[0].set_ylabel('centred response (z-scored input)')
axes[0].set_title('Centred response centroid (the shape that was clustered)')
axes[0].legend(frameon=False, fontsize=7)
axes[1].axhline(0, color='k', lw=0.6)
axes[1].set_xlabel('shared PT DPT')
axes[1].set_ylabel('human - mouse')
axes[1].set_title('Uncentred human-minus-mouse centroid')
fig.suptitle('Figure 5 - response module centroids', fontsize=12)
_save_figure(fig, 'fig05_response_module_centroids.png')


# %%
# Purpose: figure 6 - for every response module: centroid, member spread, representatives, and the
# original human and mouse fitted curves for those representatives.
fig, axes = plt.subplots(2, len(response_catalog_stats),
                         figsize=(3.1 * len(response_catalog_stats), 6.4), sharex=True,
                         squeeze=False)
# A gene has TWO row numbers: one inside the module's own curve block (for the centred response)
# and one in the eligible-gene array (for the original human and mouse curves). Indexing the block
# with the global number is an out-of-bounds read that lands on a different gene.
module_positions = {
    module: {gene: index
             for index, gene in enumerate(response_module_index[
                 (response_module_labels == module).to_numpy()])}
    for module in response_catalog_stats['module']
}

for column, module in enumerate(response_catalog_stats['module']):
    members = (response_module_labels == module).to_numpy()
    colour = MODULE_CMAP(int(module[1:]) - 1)
    member_curves = response_centred[response_eligible][members]
    centroid = np.nanmean(member_curves, axis=0)
    step = max(1, member_curves.shape[0] // 200)
    for curve in member_curves[::step]:
        axes[0][column].plot(grid, curve, color='#BBBBBB', lw=0.4, alpha=0.5)
    axes[0][column].plot(grid, centroid, color=colour, lw=2.2, label='centroid')
    for _, row in representatives[representatives['module'].eq(module)].iterrows():
        in_module = module_positions[module][row['gene']]
        overall = int(np.flatnonzero(response_module_index == row['gene'])[0])
        axes[0][column].plot(grid, member_curves[in_module], lw=1.1, alpha=0.9)
        axes[1][column].plot(grid, mouse_curves[response_eligible][overall], lw=1.4,
                             color=SPECIES_COLORS['mouse'], ls='-')
        axes[1][column].plot(grid, human_curves[response_eligible][overall], lw=1.4,
                             color=SPECIES_COLORS['human'], ls='--')
    row = response_catalog_stats.set_index('module').loc[module]
    axes[0][column].set_title(f"{module}: {row['descriptor']}\nn={int(row['n_genes'])}", fontsize=8)
    axes[0][column].axhline(0, color='k', lw=0.5)
    axes[0][column].set_ylabel('centred response' if column == 0 else '')
    axes[1][column].set_xlabel('shared PT DPT')
    axes[1][column].set_ylabel('fitted lognorm' if column == 0 else '')
axes[1][0].plot([], [], color=SPECIES_COLORS['mouse'], label='mouse')
axes[1][0].plot([], [], color=SPECIES_COLORS['human'], ls='--', label='human')
axes[1][0].legend(frameon=False, fontsize=7)
fig.suptitle('Figure 6 - response modules: centroid, member spread, representatives, and the original '
             'mouse/human curves', fontsize=12)
_save_figure(fig, 'fig06_response_module_examples.png')

print('A centred difference curve alone cannot say whether human rose, mouse fell, or both moved '
      'differently. The bottom row always returns to the original fitted curves.')
display(representatives[['module', 'gene', 'member_to_centroid_correlation', 'representativeness',
                         'response_rms', 'level_shift']].round(3))


# %% [markdown]
# ## 7 - Relating the baseline programs to the remodelling programs
#
# Every gene assigned in **both** analyses now carries a `positional_module` and a `response_module`.
# The contingency table asks whether particular normal PT programs are disproportionately affected by
# particular kinds of human-vs-mouse remodelling. Observed/expected is preferred over raw counts
# because module sizes differ, and row/column margins alone would otherwise dominate.
#
# No interpretation is imposed in advance: "late-PT genes undergo late-selective loss" or "early
# programs are preserved" are hypotheses this table can support or contradict, not assumptions.

# %%
# Purpose: the joint gene table and the positional x response contingency table.
positional_response_gene_table = gene_curve_metrics.copy()
# Both label vectors are positionally aligned with `gene_names`, which is also the row order of
# `gene_curve_metrics`, so they are inserted as arrays rather than label-joined: a label join would
# raise on a duplicated var name.
positional_response_gene_table.insert(1, 'positional_module', positional_labels.to_numpy())
positional_response_gene_table.insert(2, 'response_module', response_labels.to_numpy())
positional_response_gene_table.to_csv(OUTPUT_DIR / 'positional_response_gene_table.csv', index=False)

assigned_both = (positional_labels.to_numpy() != 'unassigned') & (response_labels.to_numpy() != 'unassigned')
observed = pd.crosstab(
    pd.Series(positional_labels.to_numpy()[assigned_both], name='positional_module'),
    pd.Series(response_labels.to_numpy()[assigned_both], name='response_module'),
)
expected = np.outer(observed.sum(axis=1), observed.sum(axis=0)) / observed.to_numpy().sum()
positional_response_crosstab = pd.DataFrame({
    'positional_module': np.repeat(observed.index.to_numpy(), observed.shape[1]),
    'response_module': np.tile(observed.columns.to_numpy(), observed.shape[0]),
    'n_genes': observed.to_numpy().ravel(),
    'expected': expected.ravel(),
})
positional_response_crosstab['observed_over_expected'] = (
    positional_response_crosstab['n_genes'] / positional_response_crosstab['expected'])
positional_response_crosstab.to_csv(OUTPUT_DIR / 'positional_response_crosstab.csv', index=False)

print(f'Genes assigned in both analyses: {int(assigned_both.sum()):,} of {gene_names.size:,} eligible '
      f'({int((positional_labels.to_numpy() == "unassigned").sum()):,} unassigned positionally, '
      f'{int((response_labels.to_numpy() == "unassigned").sum()):,} unassigned in response).')
display(positional_response_crosstab.pivot(
    index='positional_module', columns='response_module', values='observed_over_expected').round(2))


# %%
# Purpose: figure 7 - positional x response heatmap on the observed/expected scale.
ratio = positional_response_crosstab.pivot(
    index='positional_module', columns='response_module', values='observed_over_expected')
counts = positional_response_crosstab.pivot(
    index='positional_module', columns='response_module', values='n_genes')
fig, axis = plt.subplots(figsize=(1.15 * ratio.shape[1] + 3.4, 0.7 * ratio.shape[0] + 3))
image = axis.imshow(ratio.to_numpy(), cmap='RdBu_r', vmin=0, vmax=2)
for row in range(ratio.shape[0]):
    for column in range(ratio.shape[1]):
        axis.text(column, row, f'{counts.to_numpy()[row, column]}\n{ratio.to_numpy()[row, column]:.2f}',
                  ha='center', va='center', fontsize=8,
                  color='white' if abs(ratio.to_numpy()[row, column] - 1) > 0.45 else 'black')
axis.set_xticks(range(ratio.shape[1]))
axis.set_xticklabels(ratio.columns)
axis.set_yticks(range(ratio.shape[0]))
axis.set_yticklabels(ratio.index)
axis.set_xlabel('response module (human remodelling)')
axis.set_ylabel('positional module (normal mouse program)')
fig.colorbar(image, ax=axis, label='observed / expected', shrink=0.7)
axis.set_title('Figure 7 - are particular normal PT programs remodelled in particular ways?\n'
               'cell text: gene count / observed-over-expected', fontsize=11)
_save_figure(fig, 'fig07_positional_response_heatmap.png')


# %% [markdown]
# ## 8 - Pathway enrichment, AFTER the clusters are frozen
#
# Pathway annotations are loaded only now, after every cluster assignment and every cluster-number
# decision is final. No pathway expression curve is built here and no pathway member's expression is
# aggregated: this notebook is deliberately not a pathway-first trajectory analysis.
#
# The direction of travel is `unsupervised gene GAM module -> pathway annotation`.
#
# Over-representation is Fisher/hypergeometric, and the universe is **the genes that could actually
# have entered that clustering** - positional-eligible genes for positional modules, response-eligible
# genes for response modules. Testing a module against genes that could never have been assigned to
# one inflates the expected overlap and understates every enrichment. Correction is BH across the full
# relevant module x pathway family, and near-identical Reactome child terms are grouped rather than
# reported as separate discoveries.

# %%
# Purpose: pathway membership through the accepted ortholog map, with explicit coverage stages.
pathway_coverage_parts = []
for library in PATHWAY_CONFIG['libraries']:
    library_path = PATHWAY_LIBRARY_DIR / f'{library}.json'
    if not library_path.exists():
        print(f'WARNING: missing pathway library {library_path.name} - skipped')
        continue
    pathway_coverage_parts.append(build_pathway_membership(
        json.loads(library_path.read_text()),
        adata_pt.var_names,
        ortholog_map=ortholog_map,
        library_name=library,
        min_genes=PATHWAY_CONFIG['min_genes'],
        max_genes=None,
        tested=gene_names,
    ))
if not pathway_coverage_parts:
    raise RuntimeError(
        f'No pathway library could be read from {_rel(PATHWAY_LIBRARY_DIR)}; the enrichment outputs '
        'in this notebook cannot be produced without them.'
    )
pathway_coverage = pd.concat(pathway_coverage_parts, ignore_index=True)
pathway_coverage.to_csv(DIAGNOSTIC_DIR / 'pathway_membership_coverage.csv', index=False)
retained_pathways = pathway_coverage[pathway_coverage['retained']].reset_index(drop=True)
if retained_pathways.empty:
    raise RuntimeError('No pathway retains enough assayed members to be tested.')

pathway_gene_sets = {f'{row.library}: {row.pathway}': list(row.genes_present)
                     for row in retained_pathways.itertuples()}
display(pathway_coverage.groupby('library', observed=True).agg(
    n_pathways=('pathway', 'size'), n_retained=('retained', 'sum'),
    median_members=('n_assayed', 'median')).reset_index())
print(f'Pathways retained: {len(retained_pathways):,} of {len(pathway_coverage):,} across '
      f'{len(pathway_coverage_parts)} libraries.')
print('The universe is the clustering-eligible gene set, per analysis - not every expressed gene.')


# %%
# Purpose: hypergeometric enrichment of positional modules and response modules.
positional_module_members = {
    module: list(members.index)
    for module, members in positional_labels.groupby(positional_labels, observed=True)
    if module != 'unassigned'
}
response_module_members = {
    module: list(members.index)
    for module, members in response_labels.groupby(response_labels, observed=True)
    if module != 'unassigned'
}
positional_background = gene_names[positional_eligible]
response_background = gene_names[response_eligible]


def _annotate_enrichment(frame, members):
    """Split the library: path label, add fold enrichment, and attach each pathway's redundancy group."""
    if frame.empty:
        return frame
    parts = frame['gene_set'].str.split(': ', n=1, expand=True)
    frame = frame.assign(library=parts[0], pathway=parts[1],
                         fold_enrichment=frame['n_overlap'] / frame['expected_overlap'])
    if len(members):
        _, redundancy = summarize_pathway_redundancy(
            members, gene_column='genes_present', label_columns=('library', 'pathway'),
            overlap_threshold=PATHWAY_CONFIG['overlap_threshold'])
        frame = frame.merge(redundancy, on=['library', 'pathway'], how='left')
    else:
        frame = frame.assign(redundancy_group=np.nan, group_size=np.nan)
    return frame.sort_values(['module', 'p_value_adjusted', 'fold_enrichment'],
                             ascending=[True, True, False]).reset_index(drop=True)


positional_enrichment = cached_frame(
    'pt_positional_pathway_enrichment',
    lambda: enrich_modules(positional_module_members, pathway_gene_sets,
                           background=positional_background),
    root=STAGE_CACHE_DIR,
    params={'logic': 1, 'min_genes': PATHWAY_CONFIG['min_genes']},
    inputs={'labels': digest(positional_labels.to_numpy()),
            'background': digest(np.asarray(positional_background)),
            'sets': digest(sorted((name, tuple(members)) for name, members in pathway_gene_sets.items()))},
    code=code_digest(enrich_modules), enabled=STAGE_CACHE_ENABLED,
)
positional_enrichment = _annotate_enrichment(positional_enrichment, retained_pathways)
positional_enrichment.to_csv(OUTPUT_DIR / 'positional_module_pathway_enrichment.csv', index=False)

response_enrichment = cached_frame(
    'pt_response_pathway_enrichment',
    lambda: enrich_modules(response_module_members, pathway_gene_sets, background=response_background),
    root=STAGE_CACHE_DIR,
    params={'logic': 1, 'min_genes': PATHWAY_CONFIG['min_genes']},
    inputs={'labels': digest(response_labels.to_numpy()),
            'background': digest(np.asarray(response_background)),
            'sets': digest(sorted((name, tuple(members)) for name, members in pathway_gene_sets.items()))},
    code=code_digest(enrich_modules), enabled=STAGE_CACHE_ENABLED,
)
response_enrichment = _annotate_enrichment(response_enrichment, retained_pathways)
response_enrichment.to_csv(OUTPUT_DIR / 'response_module_pathway_enrichment.csv', index=False)

for _label, _frame in (('positional', positional_enrichment), ('response', response_enrichment)):
    print(f'{_label}: {len(_frame):,} module x pathway pairs tested; '
          f'{int((_frame["p_value_adjusted"] < 0.05).sum()):,} at BH q < 0.05')
display(response_enrichment[response_enrichment['p_value_adjusted'] < 0.05].head(15)[[
    'module', 'library', 'pathway', 'n_module_genes', 'n_set_genes_in_background', 'n_overlap',
    'expected_overlap', 'fold_enrichment', 'p_value', 'p_value_adjusted', 'redundancy_group',
]].round(4))
print('Enrichment never revises the clustering: the modules were final before any pathway was read.')


# %%
# Purpose: figure 8 - the strongest non-redundant response-module pathway enrichments.
def _nonredundant_top(frame, limit, q_cut=0.05):
    """One row per (module, redundancy group): the strongest significant term in each group.

    Reactome is a hierarchy of nested sets, so without this a single biology appears as a dozen
    nearly identical child terms and the figure reports the library's structure, not the result.
    """
    significant = frame[frame['p_value_adjusted'].lt(q_cut)].copy()
    if significant.empty:
        return significant
    significant['redundancy_group'] = significant['redundancy_group'].fillna(
        pd.Series(np.arange(len(significant)), index=significant.index).astype(str))
    return significant.groupby(['module', 'redundancy_group'], observed=True).head(1) \
        .groupby('module', observed=True).head(limit).reset_index(drop=True)


response_top_pathways = _nonredundant_top(response_enrichment, PATHWAY_CONFIG['top_per_module'])
if response_top_pathways.empty:
    print('No response-module pathway passes BH q < 0.05; no enrichment figure was written.')
else:
    modules = sorted(response_top_pathways['module'], key=lambda name: int(name[1:]))
    modules = list(dict.fromkeys(modules))
    labels = response_top_pathways['gene_set'].to_numpy()
    rows = [modules.index(module) for module in response_top_pathways['module']]
    fig, axis = plt.subplots(figsize=(11, 0.34 * len(response_top_pathways) + 2.2))
    scores = -np.log10(np.maximum(response_top_pathways['p_value_adjusted'].to_numpy(), 1e-300))
    scatter = axis.scatter(
        rows, np.arange(len(response_top_pathways)),
        s=18 + 9 * response_top_pathways['n_overlap'].to_numpy(),
        c=scores, cmap='viridis', edgecolor='k', linewidth=0.3)
    axis.set_yticks(np.arange(len(response_top_pathways)))
    axis.set_yticklabels(labels, fontsize=7)
    axis.set_xticks(range(len(modules)))
    axis.set_xticklabels(modules)
    axis.set_xlabel('response module')
    axis.invert_yaxis()
    axis.grid(axis='x', ls=':', lw=0.5, alpha=0.6)
    fig.colorbar(scatter, ax=axis, label='-log10 BH q', shrink=0.6)
    axis.set_title('Figure 8 - non-redundant response-module pathway enrichment\n'
                   'dot size = overlap genes; one term per redundancy group', fontsize=11)
    _save_figure(fig, 'fig08_response_pathway_enrichment.png')
    display(response_top_pathways[[
        'module', 'library', 'pathway', 'n_overlap', 'fold_enrichment', 'p_value_adjusted',
    ]].round(4))


# %% [markdown]
# ## 9 - Pathway heterogeneity across response modules
#
# A pathway is not a program. Its members can concentrate in one coherent response module, spread
# across several, or occupy **opposing** modules - and opposing members are exactly the case an
# aggregate pathway-expression curve would cancel to nothing. This table is descriptive: it reports
# where each pathway's eligible genes land, without turning into a second significance ranking.

# %%
# Purpose: how each pathway's eligible genes distribute across the response modules.
response_module_by_gene = dict(zip(gene_names, response_labels.to_numpy()))
pathway_module_distribution = []
for name, members in pathway_gene_sets.items():
    eligible_members = [gene for gene in members if response_module_by_gene.get(gene) not in (None, 'unassigned')]
    if not eligible_members:
        continue
    counts = pd.Series([response_module_by_gene[gene] for gene in eligible_members]).value_counts()
    pathway_module_distribution.append({
        'pathway': name,
        'eligible_pathway_genes': len(eligible_members),
        'dominant_response_module': counts.index[0],
        'fraction_in_dominant_module': float(counts.iloc[0] / len(eligible_members)),
        'number_of_response_modules': int(len(counts)),
        'module_distribution': '; '.join(f'{module}:{count}' for module, count in counts.items()),
    })
pathway_module_distribution = pd.DataFrame(pathway_module_distribution).sort_values(
    ['eligible_pathway_genes', 'fraction_in_dominant_module'], ascending=[False, False]
).reset_index(drop=True)
pathway_module_distribution.to_csv(OUTPUT_DIR / 'pathway_response_module_distribution.csv', index=False)

_display_rows = pathway_module_distribution[pathway_module_distribution['eligible_pathway_genes'] >= 10]
display(_display_rows.head(15).round(3))
print(f'Pathways with >= 1 eligible member in a response module: {len(pathway_module_distribution):,}')
print('Concentrated (fraction 1.0): '
      f'{int((pathway_module_distribution["fraction_in_dominant_module"] == 1.0).sum()):,}; '
      'split across >= 3 modules: '
      f'{int((pathway_module_distribution["number_of_response_modules"] >= 3).sum()):,}')
print('Split and opposing pathways are the reason an aggregate pathway curve can cancel a real '
      'heterogeneous remodelling.')


# %% [markdown]
# ## 10 - Sensitivity analyses
#
# Kept compact, and none of them displaces the primary analysis. Throughout, the two human sections are
# treated as ONE donor: omitting a section is a **section**-sensitivity check, never independent
# biological replication.

# %%
# Purpose: 10A - pooled-specimen curves versus the primary specimen-balanced curves.
pooled_response_centred = difference_curves(pooled_mouse_curves, pooled_human_curves, center=True)
pooled_response_rms = np.sqrt(np.nanmean(pooled_response_centred ** 2, axis=1))
pooled_response_peak_to_peak = _peak_to_peak(pooled_response_centred)
pooled_eligible, pooled_columns = _curve_eligibility(
    pooled_response_centred, pooled_response_peak_to_peak, RESPONSE_PEAK_TO_PEAK_FLOOR)
pooled_response = _cluster_curves(
    pooled_response_centred[pooled_eligible][:, pooled_columns], grid[pooled_columns],
    gene_names[pooled_eligible], CLUSTER_CONFIG, 'R',
    representativeness=pooled_response_rms[pooled_eligible],
)

shared_genes = (response_eligible & pooled_eligible)
balanced_shared = response_labels.to_numpy()[shared_genes]
pooled_lookup = dict(zip(gene_names[pooled_eligible], pooled_response['labels'].to_numpy()))
pooled_shared = np.array([pooled_lookup.get(gene, 'unassigned') for gene in gene_names[shared_genes]])
pooled_vs_balanced_ari = float(adjusted_rand_score(balanced_shared, pooled_shared))

jaccard_rows = []
for module in sorted(set(balanced_shared), key=lambda name: int(name[1:])):
    members = set(gene_names[shared_genes][balanced_shared == module])
    best_module, best_jaccard = '', 0.0
    for candidate in sorted(set(pooled_shared)):
        other = set(gene_names[shared_genes][pooled_shared == candidate])
        union = len(members | other)
        score = len(members & other) / union if union else 0.0
        if score > best_jaccard:
            best_module, best_jaccard = candidate, score
    jaccard_rows.append({'module': module, 'balanced_n': len(members),
                         'best_pooled_module': best_module,
                         'best_jaccard_overlap': float(best_jaccard)})
pooled_centroid_matching = pd.DataFrame(jaccard_rows)

# Each run has its own eligible set, so membership is resolved by gene name through each run's
# own row positions and both centroids are built from the SAME genes. Indexing one run's boolean
# mask into another run's array is a silent mis-alignment when the two sets differ.
response_position = {gene: index for index, gene in enumerate(gene_names[response_eligible])}
pooled_position = {gene: index for index, gene in enumerate(gene_names[pooled_eligible])}
# `response_module_labels` is indexed by `gene_names[response_eligible]`, so the module gene lists
# are read from that same space - masking the full `gene_names` would be a length mismatch.
module_gene_lists = {
    module: gene_names[response_eligible][(response_module_labels == module).to_numpy()]
    for module in response_catalog_stats['module']
}

pooled_centroid_rows = []
for module, module_genes in module_gene_lists.items():
    # Named `matched_genes`, not `shared_genes`: that name already holds the eligible-set boolean
    # array this cell's summary reads, and a list here would shadow it into an AttributeError.
    matched_genes = [gene for gene in module_genes if gene in pooled_position]
    if not matched_genes:
        continue
    primary = np.nanmean(response_centred[response_eligible][
        [response_position[gene] for gene in matched_genes]], axis=0)
    pooled = np.nanmean(pooled_response_centred[pooled_eligible][
        [pooled_position[gene] for gene in matched_genes]], axis=0)
    pooled_centroid_rows.append({'module': module, 'centroid_spearman_pooled_vs_balanced':
                                 safe_spearman(primary, pooled)})
pooled_centroid_matching = pooled_centroid_matching.merge(
    pd.DataFrame(pooled_centroid_rows), on='module', how='left')

print(f'10A pooled vs specimen-balanced response clustering: ARI = {pooled_vs_balanced_ari:.3f} over '
      f'{int(shared_genes.sum()):,} shared genes; '
      f'median centroid correlation = '
      f'{pooled_centroid_matching["centroid_spearman_pooled_vs_balanced"].median():.3f}')
display(pooled_centroid_matching.round(3))
print('The specimen-balanced clustering stays primary: the pooled fit weights tubules, so a specimen '
      'with more structures shapes it more.')


# %%
# Purpose: 10B - specimen omission, recomputing the response programs one specimen at a time.
stability_runs = {'specimen_balanced_primary': response_labels}
module_centroid_stability = {module: [] for module in response_catalog_stats['module']}
omission_rows = []
for dropped in sorted(species_by_sample):
    kept = [name for name in sorted(species_by_sample) if name != dropped]
    kept_mouse = [name for name in kept if species_by_sample[name] == 'mouse']
    kept_human = [name for name in kept if species_by_sample[name] == 'human']
    if not kept_mouse or not kept_human:
        continue                      # dropping this section would collapse a side of the comparison
    dropped_mouse = specimen_balanced_curves({name: specimen_curves[name] for name in kept_mouse})
    dropped_human = specimen_balanced_curves({name: specimen_curves[name] for name in kept_human})
    dropped_response = difference_curves(dropped_mouse, dropped_human, center=True)
    dropped_ptp = _peak_to_peak(dropped_response)
    dropped_rms = np.sqrt(np.nanmean(dropped_response ** 2, axis=1))
    dropped_eligible, dropped_columns = _curve_eligibility(
        dropped_response, dropped_ptp, RESPONSE_PEAK_TO_PEAK_FLOOR)
    dropped_result = _cluster_curves(
        dropped_response[dropped_eligible][:, dropped_columns], grid[dropped_columns],
        gene_names[dropped_eligible], CLUSTER_CONFIG, 'R',
        representativeness=dropped_rms[dropped_eligible],
    )
    labels = np.full(gene_names.size, 'unassigned', dtype=object)
    labels[dropped_eligible] = dropped_result['labels'].to_numpy()
    stability_runs[f'without_{dropped}'] = pd.Series(labels, index=gene_names)
    # Does the major response centroid recur? Recompute each reference module's centroid from THIS
    # run's own curves, using the reference membership, so a permuted module name cannot hide it.
    dropped_position = {gene: index for index, gene in enumerate(gene_names[dropped_eligible])}
    correlations = []
    for module, module_genes in module_gene_lists.items():
        # The reference membership restricted to the genes THIS run could cluster, and both centroids
        # built from that same set, so a smaller eligible set cannot masquerade as a lost module.
        matched_genes = [gene for gene in module_genes if gene in dropped_position]
        if not matched_genes:
            correlations.append(np.nan)
            module_centroid_stability[module].append(np.nan)
            continue
        primary = np.nanmean(response_centred[response_eligible][
            [response_position[gene] for gene in matched_genes]], axis=0)
        recomputed = np.nanmean(dropped_response[dropped_eligible][
            [dropped_position[gene] for gene in matched_genes]], axis=0)
        value = safe_spearman(primary, recomputed)
        correlations.append(value)
        module_centroid_stability[module].append(value)
    omission_rows.append({
        'dropped_specimen': dropped,
        'dropped_species': species_by_sample[dropped],
        'n_genes_clustered': int(dropped_eligible.sum()),
        'selected_k': dropped_result['k'],
        'median_centroid_spearman_vs_primary': float(np.nanmedian(correlations)),
        'min_centroid_spearman_vs_primary': float(np.nanmin(correlations)),
    })
omission_sensitivity = pd.DataFrame(omission_rows)
display(omission_sensitivity.round(3))

# `module_stability` supplies the label-invariant run-level numbers (ARI, retention). The
# centroid-recursion column is merged in from the omission table; the primary run is its own
# reference, so its value is 1 by construction rather than a measurement.
run_stability = module_stability(stability_runs)
_centroid_by_run = {'specimen_balanced_primary': 1.0}
_centroid_by_run.update(dict(zip(
    'without_' + omission_sensitivity['dropped_specimen'].astype(str),
    omission_sensitivity['median_centroid_spearman_vs_primary'],
)))
run_stability['median_centroid_spearman_vs_primary'] = run_stability['run'].map(_centroid_by_run)
print('Specimen omission is a SECTION sensitivity for the human side, not biological replication.')
display(run_stability.round(3))


# %%
# Purpose: 10C - independent human clustering, as a secondary sensitivity analysis only.
human_amplitude = _peak_to_peak(human_curves)
human_eligible, human_columns = _curve_eligibility(
    human_curves, human_amplitude, POSITIONAL_AMPLITUDE_FLOOR)
human_result = _cluster_curves(
    human_curves[human_eligible][:, human_columns], grid[human_columns],
    gene_names[human_eligible], CLUSTER_CONFIG, 'H',
    representativeness=human_amplitude[human_eligible],
)
human_labels = np.full(gene_names.size, 'unassigned', dtype=object)
human_labels[human_eligible] = human_result['labels'].to_numpy()
human_labels = pd.Series(human_labels, index=gene_names, name='human_module')

# Both standardised matrices are compared on the grid columns BOTH runs cover, so a narrower
# support on one side cannot be read as a shape difference.
shared_columns = positional_columns & human_columns
positional_standardized = positional['standardized'][:, shared_columns[positional_columns]]
human_standardized = human_result['standardized'][:, shared_columns[human_columns]]

human_module_matching = []
for index, module in enumerate(positional['catalog']['module']):
    members = (positional_labels.to_numpy() == module)
    mouse_centroid = np.nanmean(positional_standardized[members[positional_eligible]], axis=0)
    for human_module in sorted(set(human_result['labels']), key=lambda name: int(name[1:])):
        human_members = (human_labels.to_numpy() == human_module)
        human_centroid = np.nanmean(human_standardized[human_members[human_eligible]], axis=0)
        overlap = set(gene_names[members]) & set(gene_names[human_members])
        union = set(gene_names[members]) | set(gene_names[human_members])
        human_module_matching.append({
            'positional_module': module,
            'human_module': human_module,
            'centroid_spearman': safe_spearman(mouse_centroid, human_centroid),
            'member_jaccard': (len(overlap) / len(union)) if union else 0.0,
        })
human_module_matching = pd.DataFrame(human_module_matching)
best_human_match = (human_module_matching.sort_values('centroid_spearman', ascending=False)
                    .groupby('positional_module', observed=True).head(1).reset_index(drop=True))

_mouse_shared = positional_labels.to_numpy() != 'unassigned'
_human_shared = human_labels.to_numpy() != 'unassigned'
_human_only_shared = _human_shared & np.isin(gene_names, gene_names[human_eligible])
_ari_mask = _mouse_shared & _human_shared
human_vs_mouse_ari = float(adjusted_rand_score(
    positional_labels.to_numpy()[_ari_mask], human_labels.to_numpy()[_ari_mask]))
human_clustering_summary = pd.DataFrame([{
    'human_modules': human_result['k'],
    'mouse_positional_modules': positional['k'],
    'human_shared_module_genes': int(_human_only_shared.sum()),
    'adjusted_rand_index_human_vs_mouse': human_vs_mouse_ari,
}])
human_clustering_summary.to_csv(DIAGNOSTIC_DIR / 'human_clustering_summary.csv', index=False)
human_module_matching.to_csv(DIAGNOSTIC_DIR / 'human_module_matching.csv', index=False)
display(human_clustering_summary.round(3))
display(best_human_match.round(3))
print('This asks only whether the broad positional architecture is recognisable independently in '
      'human. Cluster-number differences alone are not biology, and this is NOT the species '
      'comparison - that remains the centred human-minus-mouse response clustering.')


# %%
# Purpose: figure 9 - the stability summary, and the response-module stability table.
module_stability_lookup = pd.DataFrame([
    {'module': module,
     'mean_centroid_spearman_across_runs': float(np.nanmean(module_centroid_stability[module])),
     'min_centroid_spearman_across_runs': float(np.nanmin(module_centroid_stability[module]))}
    for module in response_catalog_stats['module']
])
run_stability.to_csv(OUTPUT_DIR / 'response_module_stability.csv', index=False)

fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.9))
runs = run_stability['run'].tolist()
axes[0].bar(runs, run_stability['adjusted_rand_index_vs_reference'], color='#4C9BD3')
axes[0].set_ylabel('ARI vs primary')
axes[0].set_title('Module agreement under perturbation')
axes[0].tick_params(axis='x', rotation=30)
axes[1].bar(runs, run_stability['mean_reference_module_retention'], color='#59A14F')
axes[1].set_ylabel('mean reference-module retention')
axes[1].set_title('Do reference modules stay together?')
axes[1].tick_params(axis='x', rotation=30)
axes[2].bar(module_stability_lookup['module'],
            module_stability_lookup['mean_centroid_spearman_across_runs'], color='#E15759')
axes[2].set_ylim(0, 1.05)
axes[2].set_ylabel('mean centroid Spearman')
axes[2].set_title('Do the major centroids recur?')
fig.suptitle('Figure 9 - stability summary (specimen omission and pooled-vs-balanced mixture)',
             fontsize=12)
_save_figure(fig, 'fig09_stability_summary.png')

display(run_stability.round(3))
display(module_stability_lookup.round(3))
print('Small or unstable modules must not be over-interpreted; the retention and centroid columns are '
      'the guard against that.')


# %% [markdown]
# ## 11 - Final module catalogs
#
# The catalogs are written now, after the stability analysis, so each response module can carry its
# stability beside its description. Descriptors are read off the actual centroids; small or unstable
# modules are labelled as such rather than over-interpreted.

# %%
# Purpose: the two final catalogs, each with its top non-redundant pathway enrichments.
def _top_enrichment_summary(frame, module, limit, q_cut=0.05):
    """The strongest non-redundant pathway terms of one module, as one readable string."""
    significant = _nonredundant_top(frame[frame['module'].eq(module)], limit, q_cut=q_cut)
    if significant.empty:
        return ''
    return '; '.join(
        f'{row.pathway} (q={row.p_value_adjusted:.1e}, O/E={row.fold_enrichment:.1f})'
        for row in significant.itertuples())


mouse_positional_module_catalog = positional['catalog'].copy()
mouse_positional_module_catalog['top_nonredundant_pathways'] = [
    _top_enrichment_summary(positional_enrichment, module, PATHWAY_CONFIG['top_per_module'])
    for module in mouse_positional_module_catalog['module']
]
mouse_positional_module_catalog = mouse_positional_module_catalog[[
    'module', 'n_genes', 'descriptor', 'direction', 'peak_position', 'median_amplitude',
    'median_member_to_centroid_correlation', 'representative_genes', 'top_nonredundant_pathways',
]]
mouse_positional_module_catalog.to_csv(OUTPUT_DIR / 'mouse_positional_module_catalog.csv', index=False)

response_module_catalog = response_catalog_stats.merge(module_stability_lookup, on='module', how='left')
response_module_catalog['top_nonredundant_pathways'] = [
    _top_enrichment_summary(response_enrichment, module, PATHWAY_CONFIG['top_per_module'])
    for module in response_module_catalog['module']
]
response_module_catalog['module_status'] = np.where(
    response_module_catalog['mean_centroid_spearman_across_runs'] >= 0.8, 'stable', 'interpret with care')
response_module_catalog = response_module_catalog[[
    'module', 'n_genes', 'descriptor', 'direction', 'peak_position', 'module_status',
    'mean_centroid_spearman_across_runs', 'min_centroid_spearman_across_runs',
    'representative_genes', 'median_response_rms', 'median_level_shift',
    'median_early_difference', 'median_middle_difference', 'median_late_difference',
    'median_mouse_amplitude', 'median_human_amplitude', 'median_mouse_human_curve_spearman',
    'median_member_to_centroid_correlation', 'top_nonredundant_pathways',
]]
response_module_catalog.to_csv(OUTPUT_DIR / 'response_module_catalog.csv', index=False)

pd.DataFrame([{'figure': name} for name in figure_index]).to_csv(
    DIAGNOSTIC_DIR / 'figure_index.csv', index=False)

print(f'Positional catalog: {len(mouse_positional_module_catalog)} modules.')
display(mouse_positional_module_catalog.round(3))
print(f'Response catalog: {len(response_module_catalog)} modules.')
display(response_module_catalog.drop(columns=['top_nonredundant_pathways']).round(3))
print('A small or unstable module is reported as such instead of being over-interpreted.')


# %% [markdown]
# ## 12 - Shared positional archetypes and cross-species program switching
#
# > Mouse and human PT may share a small repertoire of conserved positional-expression shapes, while
# > individual genes switch which positional program they occupy between species.
#
# Why this section exists, and why it is a different question from sections 5, 6 or 10C:
#
# * absolute normalised human and mouse expression levels are **not** assumed to be directly
#   comparable molecular abundances, so nothing below compares a human level with a mouse level;
# * the strongest current result is that mouse and human contain similar broad curve shapes while the
#   **independently** clustered gene memberships agree poorly - that is what section 10C's low ARI and
#   its centroid matching show, and it is exactly the label-matching problem this section attacks;
# * the current centred difference-curve response clusters (section 6) are strongly associated with
#   the original mouse positional programs (section 7's positional x response table), so their
#   interpretation may really be about which positional program a gene occupies in each species;
# * therefore test the hypothesis **directly**: do the same curve archetypes exist in both species,
#   and do individual genes switch between them?
#
# The unit of analysis is the **within-species standardised fitted curve**. Each species' curve is
# centred and scaled on its own, so a gene's mean level and its amplitude are removed and only its
# positional shape survives. Nothing here subtracts human from mouse: the only question is whether a
# gene's shape is conserved or changes.
#
# Three conventions carry through the whole section:
#
# * the archetypes are learned from **both species at once**, never as two separately named
#   clusterings that then have to be matched;
# * no pathway information enters until the archetypes, the transition matrix and the assignments are
#   frozen (12.10), exactly as sections 8 and 9 are kept downstream of the clustering;
# * every "switched" or "conserved" statement is **descriptive**. It is not a species-level test: the
#   human side is two sections of ONE donor, so the inference unit is the specimen and the comparison
#   is 2-vs-1-donor.

# %%
# Purpose: every knob section 12 uses, declared before any curve is standardised.
ARCHETYPE_CONFIG = {
    # Section 5's positional floor is the absolute amplitude a curve must carry in EACH species
    # before it may be standardised. The ladder is 12.12's sensitivity analysis, in multiples of it.
    'amplitude_floor': POSITIONAL_AMPLITUDE_FLOOR,
    'amplitude_floor_ladder': (1.0, 2.0, 4.0),
    'min_genes': 20,
    # Reading aids for the descriptive distribution of direct shape correlation. No biological
    # significance threshold is defined from these bins.
    'conservation_bins': (-np.inf, -0.5, 0.0, 0.5, 0.8, np.inf),
    'conservation_labels': ('< -0.5', '-0.5 to 0', '0 to 0.5', '0.5 to 0.8', '> 0.8'),
    # A transition, or a conserved set, enters 12.10's pathway test only above this size: a handful of
    # genes cannot support an over-representation claim, and testing every tiny cell would inflate the
    # family that BH corrects across.
    'pathway_min_transition_genes': 30,
    'coherence_min_pathway_genes': 15,
    # 12.11's empirical null: random gene sets matched to each pathway for size AND for approximate
    # mouse amplitude, so a low entropy is not read as concentration when size alone explains it.
    'null_permutations': 200,
    'seed': 0,
    'representative_transitions': 3,
    'representative_genes_per_transition': 4,
    'numerical_epsilon': 1e-12,
}
ARCHETYPE_FLOOR = ARCHETYPE_CONFIG['amplitude_floor']
print(f'Section 12 thresholds: amplitude floor per species {ARCHETYPE_FLOOR}, '
      f'ladder {ARCHETYPE_CONFIG["amplitude_floor_ladder"]}, '
      f'pathway-tested transitions >= {ARCHETYPE_CONFIG["pathway_min_transition_genes"]} genes.')


# %%
# Purpose: 12.1 - the shared shape space in one reusable unit: standardise both species, learn ONE
# joint PCA, cluster the stacked matrix once.
def _standardize_species(curves, columns):
    """Z_g(s) for one species: centre and scale THAT species' curve across the shared pseudospace.

    Centring removes the gene's mean level, scaling removes its amplitude, so positional shape is all
    that is left. The two species are standardised independently and are never subtracted from one
    another: their normalised levels are not assumed to be comparable molecular abundances.
    """
    return zscore_rows(np.asarray(curves, dtype=float)[:, columns])


def _genuine_shape_eligibility(mouse_block, human_block, grid, floor, min_finite_fraction=0.5):
    """Genes carrying real positional shape in BOTH species, plus the grid columns both are seen on.

    The SAME absolute amplitude floor section 5 applies before standardising is applied here twice,
    once per species. A nearly flat curve divided by its own tiny spread turns numerical noise into an
    apparently strong shape, so a gene that is essentially flat in either species is excluded rather
    than standardised. Returns ``(eligible, columns, mouse_ok, human_ok)``.
    """
    mouse_block = np.asarray(mouse_block, dtype=float)
    human_block = np.asarray(human_block, dtype=float)
    mouse_ok, mouse_columns = _curve_eligibility(mouse_block, _peak_to_peak(mouse_block), floor)
    human_ok, human_columns = _curve_eligibility(human_block, _peak_to_peak(human_block), floor)
    columns = mouse_columns & human_columns
    if columns.sum() < 3:                  # too little shared support left to describe any shape
        columns = np.ones(mouse_block.shape[1], dtype=bool)
    finite = (np.isfinite(mouse_block[:, columns]).all(axis=1)
              & np.isfinite(human_block[:, columns]).all(axis=1))
    return mouse_ok & human_ok & finite, columns, mouse_ok, human_ok


def _row_block_pearson(left, right, epsilon):
    """Per-row Pearson correlation of two aligned blocks, over the columns a row actually carries."""
    a = np.asarray(left, dtype=float)
    b = np.asarray(right, dtype=float)
    finite = np.isfinite(a) & np.isfinite(b)
    count = np.maximum(finite.sum(axis=1, keepdims=True), 1)
    mean_a = np.where(finite, a, 0.0).sum(axis=1, keepdims=True) / count
    mean_b = np.where(finite, b, 0.0).sum(axis=1, keepdims=True) / count
    a = np.where(finite, a - mean_a, 0.0)
    b = np.where(finite, b - mean_b, 0.0)
    numerator = (a * b).sum(axis=1)
    denominator = np.sqrt((a ** 2).sum(axis=1) * (b ** 2).sum(axis=1))
    return numerator / np.where(denominator > epsilon, denominator, np.nan)


def _median_or_nan(values):
    """Median of an array that may legitimately be empty (an archetype no gene switched out of)."""
    values = np.asarray(values, dtype=float)
    return float(np.nanmedian(values)) if values.size else np.nan


def _jaccard(left, right):
    left, right = set(left), set(right)
    union = left | right
    return len(left & right) / len(union) if union else 0.0


def _shared_archetype_solution(mouse_block, human_block, grid, feature_names, floor, min_genes):
    """One self-contained run of 12.1 -> 12.4, so 12.12 can repeat it at a stricter floor.

    The stacked matrix has two rows per gene - ``mouse:gene`` and ``human:gene`` - over the shared
    grid, and it is standardised, PCA'd and clustered exactly once. That joint step is what makes the
    archetypes COMMON to both species instead of two clusterings whose labels then need matching.
    Everything the section reports is returned rather than recomputed.
    """
    mouse_block = np.asarray(mouse_block, dtype=float)
    human_block = np.asarray(human_block, dtype=float)
    feature_names = np.asarray(feature_names)
    grid = np.asarray(grid, dtype=float)
    eligible, columns, mouse_ok, human_ok = _genuine_shape_eligibility(
        mouse_block, human_block, grid, floor)
    if int(eligible.sum()) < int(min_genes):
        raise ValueError(
            f'only {int(eligible.sum())} genes carry shape above amplitude {floor} in both species; '
            f'{int(min_genes)} are needed before a shared archetype space can be learned')

    mouse_standardized = _standardize_species(mouse_block, columns)[eligible]
    human_standardized = _standardize_species(human_block, columns)[eligible]
    standardized = np.vstack([mouse_standardized, human_standardized])
    genes = feature_names[eligible]
    row_ids = np.concatenate([np.asarray([f'mouse:{gene}' for gene in genes], dtype=object),
                              np.asarray([f'human:{gene}' for gene in genes], dtype=object)])
    row_species = np.array(['mouse'] * genes.size + ['human'] * genes.size)
    row_amplitude = np.concatenate([_peak_to_peak(mouse_block[eligible][:, columns]),
                                    _peak_to_peak(human_block[eligible][:, columns])])

    result = _cluster_curves(
        standardized, grid[columns], row_ids, CLUSTER_CONFIG, 'A',
        representativeness=row_amplitude,
    )
    labels = result['labels'].to_numpy()
    archetypes = sorted(set(labels), key=lambda name: int(name[1:]))
    is_mouse = row_species == 'mouse'
    centroid_combined, centroid_mouse, centroid_human = {}, {}, {}
    for archetype in archetypes:
        members = labels == archetype
        centroid_combined[archetype] = np.nanmean(standardized[members], axis=0)
        centroid_mouse[archetype] = np.nanmean(standardized[members & is_mouse], axis=0)
        centroid_human[archetype] = np.nanmean(standardized[members & ~is_mouse], axis=0)

    return {
        'eligible': eligible, 'columns': columns, 'mouse_ok': mouse_ok, 'human_ok': human_ok,
        'genes': genes, 'grid': grid[columns], 'standardized': standardized,
        'mouse_standardized': mouse_standardized, 'human_standardized': human_standardized,
        'row_ids': row_ids, 'row_species': row_species, 'labels': labels, 'archetypes': archetypes,
        'mouse_archetype': labels[is_mouse], 'human_archetype': labels[~is_mouse],
        'centroid_combined': centroid_combined, 'centroid_mouse': centroid_mouse,
        'centroid_human': centroid_human, 'pca': result['pca'], 'scores': result['scores'],
        'n_pc': result['n_pc'], 'cumulative_variance': result['cumulative_variance'],
        'diagnostics': result['diagnostics'], 'k': result['k'],
        'best_silhouette': result['best_silhouette'], 'loadings': result['loadings'],
        'floor': float(floor),
    }


shared = _shared_archetype_solution(
    mouse_curves, human_curves, grid, gene_names, ARCHETYPE_FLOOR, ARCHETYPE_CONFIG['min_genes'])
archetype_eligible = shared['eligible']
archetype_columns = shared['columns']
archetype_grid = shared['grid']
archetype_genes = shared['genes']
Z_mouse = shared['mouse_standardized']        # Z_mouse_g(s), eligible genes only
Z_human = shared['human_standardized']        # Z_human_g(s), eligible genes only
archetype_mouse_label = shared['mouse_archetype']
archetype_human_label = shared['human_archetype']
shared_archetypes = shared['archetypes']

archetype_eligibility_table = pd.DataFrame({
    'gene': np.asarray(gene_names),
    'mouse_peak_to_peak': gene_mouse_amplitude,
    'human_peak_to_peak': gene_human_amplitude,
    'eligible_mouse': shared['mouse_ok'],
    'eligible_human': shared['human_ok'],
    'eligible_both_species': archetype_eligible,
})
archetype_eligibility_table['exclusion_reason'] = np.where(
    archetype_eligible, 'eligible: real shape in both species',
    np.where(
        ~np.asarray(shared['mouse_ok']) & ~np.asarray(shared['human_ok']),
        'excluded: essentially flat in both species',
        np.where(np.asarray(shared['mouse_ok']),
                 'excluded: essentially flat in human only',
                 'excluded: essentially flat in mouse only')))

print(f'Genes in the eligible universe: {gene_names.size:,}')
print(f'  carrying shape above the floor in mouse: {int(shared["mouse_ok"].sum()):,}')
print(f'  carrying shape above the floor in human: {int(shared["human_ok"].sum()):,}')
print(f'  carrying shape above the floor in BOTH (enters section 12): '
      f'{int(archetype_eligible.sum()):,}')
print(f'  excluded as essentially flat in exactly ONE species (the other was fine): '
      f'{int((np.asarray(shared["mouse_ok"]) ^ np.asarray(shared["human_ok"])).sum()):,}')
print(f'  excluded as essentially flat in BOTH species: '
      f'{int((~np.asarray(shared["mouse_ok"]) & ~np.asarray(shared["human_ok"])).sum()):,}')
print(f'Shared support the shape space is defined on: {int(archetype_columns.sum())} of {grid.size} '
      'grid points.')
print(f'The floor is {ARCHETYPE_FLOOR} peak-to-peak per species - section 5\'s positional floor, so a '
      'gene excluded here is excluded for the same amplitude reason it would be there.')


# %%
# Purpose: 12.2 - direct per-gene mouse-versus-human shape similarity, computed before any clustering.
def _paired_shape_metrics(mouse_block, human_block, epsilon):
    """Shape-only metrics for each gene's mouse and human curve.

    Both blocks arrive z-scored WITHIN their own species, so the difference, the correlation and the
    RMS below describe positional shape alone: a level shift or an amplitude change cannot appear in
    any of them. Cosine similarity is deliberately not duplicated here - with both curves centred on
    zero, cosine and Pearson are the same number, so a separate column would only look like a second
    piece of evidence.
    """
    mouse = np.asarray(mouse_block, dtype=float)
    human = np.asarray(human_block, dtype=float)
    residual = human - mouse
    return pd.DataFrame({
        'shape_correlation': _row_block_pearson(mouse, human, epsilon),
        'shape_spearman': [safe_spearman(left, right) for left, right in zip(mouse, human)],
        'shape_rms_difference': np.sqrt((residual ** 2).mean(axis=1)),
        'shape_max_abs_difference': np.abs(residual).max(axis=1),
    })


archetype_shape_metrics = _paired_shape_metrics(
    Z_mouse, Z_human, ARCHETYPE_CONFIG['numerical_epsilon'])
archetype_shape_metrics.insert(0, 'gene', archetype_genes)
archetype_shape_correlation = archetype_shape_metrics['shape_correlation'].to_numpy(dtype=float)
archetype_shape_rms = archetype_shape_metrics['shape_rms_difference'].to_numpy(dtype=float)

# One row per gene in the universe: the eligibility status of 12.1 and the direct metrics of 12.2 in
# the same gene-level table. The merge is a left join on a unique key, so it preserves `gene_names`
# order - which is what lets the arrays above stay positionally aligned with the curve matrices.
shared_archetype_gene_metrics = archetype_eligibility_table.merge(
    archetype_shape_metrics, on='gene', how='left')
shared_archetype_gene_metrics.to_csv(OUTPUT_DIR / 'shared_archetype_gene_metrics.csv', index=False)

_conservation_bins = pd.cut(
    archetype_shape_correlation, bins=list(ARCHETYPE_CONFIG['conservation_bins']),
    labels=list(ARCHETYPE_CONFIG['conservation_labels']))
print(f'Genes with a comparable standardised shape in both species: '
      f'{archetype_shape_correlation.size:,}')
print(f'  median Pearson shape correlation:  {np.nanmedian(archetype_shape_correlation):.3f}')
print(f'  median Spearman shape correlation: '
      f'{np.nanmedian(archetype_shape_metrics["shape_spearman"]):.3f}')
print(f'  median RMS shape difference:       {np.nanmedian(archetype_shape_rms):.3f} z units')
print(f'  correlation > 0.8:  {np.nanmean(archetype_shape_correlation > 0.8):.1%}')
print(f'  correlation > 0.5:  {np.nanmean(archetype_shape_correlation > 0.5):.1%}')
print(f'  correlation in (-0.5, 0.5): {np.nanmean((archetype_shape_correlation > -0.5) & (archetype_shape_correlation < 0.5)):.1%}')
print(f'  correlation < 0:    {np.nanmean(archetype_shape_correlation < 0):.1%}')
print(f'  correlation < -0.5: {np.nanmean(archetype_shape_correlation < -0.5):.1%}')
print('This is descriptive. No biological significance threshold is defined from these values, and '
      'the metrics describe the modules found later rather than selecting the genes that get clustered.')


# %%
# Purpose: figure 10 - the distribution of direct gene-level human-mouse shape similarity.
fig, axes = plt.subplots(1, 4, figsize=(15.5, 3.5))
_pearson_finite = archetype_shape_correlation[np.isfinite(archetype_shape_correlation)]
_spearman_finite = archetype_shape_metrics['shape_spearman'].to_numpy(dtype=float)
_spearman_finite = _spearman_finite[np.isfinite(_spearman_finite)]
_rms_finite = archetype_shape_rms[np.isfinite(archetype_shape_rms)]

axes[0].hist(_pearson_finite, bins=80, color='#6B7A8F')
axes[0].axvline(0, color='k', lw=0.8)
axes[0].axvline(np.nanmedian(_pearson_finite), color='#D55E00', ls='--', lw=1.2,
                label=f'median {np.nanmedian(_pearson_finite):.2f}')
axes[0].set_xlabel('Pearson shape correlation')
axes[0].set_ylabel('genes')
axes[0].set_title('Direct shape conservation')
axes[0].legend(frameon=False, fontsize=8)
axes[1].hist(_spearman_finite, bins=80, color='#4C9BD3')
axes[1].axvline(0, color='k', lw=0.8)
axes[1].set_xlabel('Spearman shape correlation')
axes[1].set_ylabel('genes')
axes[1].set_title('Rank-based shape conservation')
axes[2].hist(_rms_finite, bins=80, color='#E15759')
axes[2].set_xlabel('RMS shape difference (z units)')
axes[2].set_ylabel('genes')
axes[2].set_title('Shape difference magnitude')
_bin_counts = _conservation_bins.value_counts().reindex(
    list(ARCHETYPE_CONFIG['conservation_labels'])).fillna(0)
axes[3].bar(range(len(_bin_counts)), _bin_counts.to_numpy(),
            color=plt.get_cmap('coolwarm')(np.linspace(0.05, 0.95, len(_bin_counts))))
axes[3].set_xticks(range(len(_bin_counts)))
axes[3].set_xticklabels(_bin_counts.index, rotation=30, ha='right', fontsize=8)
axes[3].set_ylabel('genes')
axes[3].set_title('Correlation bins')
for position, count in enumerate(_bin_counts.to_numpy()):
    axes[3].text(position, count, f'{int(count)}', ha='center', va='bottom', fontsize=7)
fig.suptitle('Figure 10 - direct gene-level mouse-vs-human shape similarity, before any clustering',
             fontsize=12)
_save_figure(fig, 'fig10_shared_shape_similarity.png')


# %% [markdown]
# ### 12.3 - one common positional shape space

# %%
# Purpose: 12.3 - learn ONE shape space jointly from both species' standardised curves.
shared_standardized = shared['standardized']
shared_pca = shared['pca']
shared_scores = shared['scores']
shared_n_pc = shared['n_pc']
shared_cumulative = shared['cumulative_variance']
shared_labels = shared['labels']
archetype_row_ids = shared['row_ids']
archetype_row_species = shared['row_species']
archetype_is_mouse = archetype_row_species == 'mouse'
archetype_transition_labels = np.asarray(
    [f'{mouse}->{human}' for mouse, human in zip(archetype_mouse_label, archetype_human_label)],
    dtype=object)

shared_pca_scores = pd.DataFrame({
    'row_id': archetype_row_ids,
    'gene': np.concatenate([archetype_genes, archetype_genes]),
    'species': archetype_row_species,
})
for component in range(shared_n_pc):
    shared_pca_scores[f'pc{component + 1}_score'] = shared_scores[:, component]

print(f'Combined standardised curve matrix: {shared_standardized.shape[0]:,} rows '
      f'({archetype_genes.size:,} genes x 2 species) x {shared_standardized.shape[1]} shared grid '
      'points.')
print(f'Retained shared PCs: {shared_n_pc} '
      f'(cumulative variance {shared_cumulative[shared_n_pc - 1]:.3f}; '
      f'PC1 alone {shared_cumulative[0]:.1%}).')
print('The PCA is learned across both species at once, so its loading curves are the shared repertoire '
      'of positional shapes - not a mouse space that human is projected into afterwards.')


# %%
# Purpose: figure 11 - shared positional PCA: variance explained and the leading loading curves.
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4))
_components = np.arange(1, min(len(shared_cumulative), CLUSTER_CONFIG['pc_max_components']) + 1)
axes[0].bar(_components,
            shared_cumulative[:_components.size]
            - np.concatenate([[0], shared_cumulative[:_components.size - 1]]),
            color='#6B7A8F', label='per component')
axes[0].plot(_components, shared_cumulative[:_components.size], color='#D55E00', marker='o', ms=4,
             label='cumulative')
axes[0].axvline(shared_n_pc, color='k', ls=':', lw=1)
axes[0].annotate(f'retained: {shared_n_pc} PCs', (shared_n_pc, 0.5), xytext=(6, 0),
                 textcoords='offset points', fontsize=8)
axes[0].set_xlabel('principal component')
axes[0].set_ylabel('variance explained')
axes[0].set_title('Shared PCA variance (both species)')
axes[0].legend(frameon=False, fontsize=8)

shared_loadings = shared['loadings']
for component in range(shared_n_pc):
    axes[1].plot(shared_loadings['pseudospace'], shared_loadings[f'pc{component + 1}_loading'],
                 label=f'PC{component + 1}')
axes[1].axhline(0, color='k', lw=0.6)
axes[1].set_xlabel('shared PT DPT')
axes[1].set_ylabel('loading (unit-norm)')
axes[1].set_title('Shared loading curves = the repertoire of positional shapes')
axes[1].legend(frameon=False, fontsize=8)
fig.suptitle('Figure 11 - one common positional shape space, learned jointly across species',
             fontsize=12)
_save_figure(fig, 'fig11_shared_archetype_pca.png')


# %% [markdown]
# ### 12.4 - cluster the shared space once

# %%
# Purpose: 12.4 - cluster the combined mouse+human curve observations ONCE in the retained shared PCs.
_k_sweep_figure(shared['diagnostics'],
                'Shared archetype clustering: k diagnostics '
                '(no species label and no pathway used to choose k)',
                'fig12a_shared_archetype_k_diagnostics.png')
display(_format_k_diagnostics(shared['diagnostics'], shared['k'], shared['best_silhouette'])[
    ['k', 'silhouette', 'smallest_module', 'largest_module_fraction', 'selected']
].round(4))
shared['diagnostics'][['k', 'silhouette', 'smallest_module', 'largest_module_fraction']].to_csv(
    DIAGNOSTIC_DIR / 'shared_archetype_k_diagnostics.csv', index=False)
print('k is chosen by the same rule as sections 5, 6 and 10C: best silhouette among the ks that avoid '
      'a tiny module, then the smallest k within tolerance. Species labels and pathway enrichment are '
      'not inputs to that choice.')


# %%
# Purpose: 12.4 - describe each shared archetype: combined, mouse-only and human-only centroids.
archetype_catalog_rows = []
for archetype in shared_archetypes:
    members = shared_labels == archetype
    mouse_members = members & archetype_is_mouse
    human_members = members & ~archetype_is_mouse
    combined_centroid = shared['centroid_combined'][archetype]
    descriptors = curve_descriptors(combined_centroid[None, :], archetype_grid).iloc[0].to_dict()
    direction = descriptors['direction']
    archetype_catalog_rows.append({
        'archetype': archetype,
        'n_curves': int(members.sum()),
        'n_mouse_curves': int(mouse_members.sum()),
        'n_human_curves': int(human_members.sum()),
        'fraction_mouse': float(mouse_members.sum() / members.sum()),
        'fraction_human': float(human_members.sum() / members.sum()),
        'peak_position': descriptors['peak_position'],
        'half_max_width': descriptors['half_max_width'],
        'monotonicity': descriptors['monotonicity'],
        'direction': direction,
        'descriptor': DIRECTION_DESCRIPTOR.get(direction, direction),
        # Within-archetype shape similarity: how tightly each species' own members follow the shared
        # centroid. A high value in both columns means the archetype is one shape, not an average of
        # two different ones.
        'median_mouse_member_to_centroid_spearman': _median_or_nan(
            _row_spearman(shared_standardized[mouse_members], combined_centroid)),
        'median_human_member_to_centroid_spearman': _median_or_nan(
            _row_spearman(shared_standardized[human_members], combined_centroid)),
        'centroid_spearman_mouse_vs_human': safe_spearman(
            shared['centroid_mouse'][archetype], shared['centroid_human'][archetype]),
    })
archetype_catalog = pd.DataFrame(archetype_catalog_rows)

shared_pca_scores.insert(3, 'archetype', shared_labels)
shared_pca_scores.to_csv(OUTPUT_DIR / 'shared_archetype_pca_scores.csv', index=False)
display(archetype_catalog[[
    'archetype', 'n_curves', 'n_mouse_curves', 'n_human_curves', 'fraction_mouse', 'descriptor',
    'median_mouse_member_to_centroid_spearman', 'median_human_member_to_centroid_spearman',
    'centroid_spearman_mouse_vs_human',
]].round(3))
print(f'Shared archetypes: {len(shared_archetypes)} ({", ".join(shared_archetypes)}) over '
      f'{shared_standardized.shape[0]:,} standardised curves.')
_represented_both = ((archetype_catalog['n_mouse_curves'] > 0)
                     & (archetype_catalog['n_human_curves'] > 0))
print(f'  represented by both species: {int(_represented_both.sum())} of {len(archetype_catalog)}')
print(f'  median mouse-vs-human centroid agreement within an archetype: '
      f'{archetype_catalog["centroid_spearman_mouse_vs_human"].median():.3f}')
print('An archetype that both species occupy AND whose mouse and human centroids agree is a shape the '
      'two species genuinely share; one they occupy with disagreeing centroids is not.')


# %%
# Purpose: figure 12 - all shared archetype centroids, with the mouse-only and human-only centroids
# of each archetype overlaid.
fig, axes = plt.subplots(1, 2, figsize=(13, 4.4), gridspec_kw={'width_ratios': [1, 1.1]})
for position, archetype in enumerate(shared_archetypes):
    colour = MODULE_CMAP(position % 10)
    row = archetype_catalog.set_index('archetype').loc[archetype]
    axes[0].plot(archetype_grid, shared['centroid_combined'][archetype], color=colour, lw=2,
                 label=f"{archetype} ({row['descriptor']}; m={int(row['n_mouse_curves'])}, "
                       f"h={int(row['n_human_curves'])})")
axes[0].axhline(0, color='k', lw=0.6)
axes[0].set_xlabel('shared PT DPT')
axes[0].set_ylabel('mean z-score over pseudospace')
axes[0].set_title('Shared archetype centroids (both species pooled)')
axes[0].legend(frameon=False, fontsize=6)
for position, archetype in enumerate(shared_archetypes):
    colour = MODULE_CMAP(position % 10)
    axes[1].plot(archetype_grid, shared['centroid_mouse'][archetype], color=colour, lw=1.8, ls='-',
                 label=f'{archetype} mouse')
    axes[1].plot(archetype_grid, shared['centroid_human'][archetype], color=colour, lw=1.8, ls='--')
axes[1].axhline(0, color='k', lw=0.6)
axes[1].set_xlabel('shared PT DPT')
axes[1].set_ylabel('mean z-score over pseudospace')
axes[1].set_title('Same archetype: mouse members (solid) vs human members (dashed)')
if len(shared_archetypes) <= 6:
    axes[1].legend(frameon=False, fontsize=7, ncol=2)
else:
    axes[1].text(0.02, 0.02, 'colour = archetype; solid = mouse, dashed = human',
                 transform=axes[1].transAxes, fontsize=7, style='italic')
fig.suptitle('Figure 12 - shared archetypes: are these shapes both species actually occupy?',
             fontsize=12)
_save_figure(fig, 'fig12b_shared_archetype_centroids.png')


# %%
# Purpose: 12.5 - every eligible gene gets a mouse archetype AND a human archetype.
positional_by_gene = dict(zip(gene_names, positional_labels.to_numpy()))
response_by_gene = dict(zip(gene_names, response_labels.to_numpy()))
human_independent_by_gene = dict(zip(gene_names, human_labels.to_numpy()))

archetype_assignments = pd.DataFrame({
    'gene': archetype_genes,
    'mouse_archetype': archetype_mouse_label,
    'human_archetype': archetype_human_label,
})
archetype_assignments['archetype_status'] = np.where(
    archetype_assignments['mouse_archetype'].to_numpy()
    == archetype_assignments['human_archetype'].to_numpy(),
    'conserved', 'switched')
archetype_assignments = archetype_assignments.merge(
    archetype_shape_metrics, on='gene', how='left')
archetype_assignments['mouse_amplitude'] = gene_mouse_amplitude[archetype_eligible]
archetype_assignments['human_amplitude'] = gene_human_amplitude[archetype_eligible]
archetype_assignments['response_peak_to_peak'] = gene_response_peak_to_peak[archetype_eligible]
archetype_assignments['response_rms'] = gene_response_rms[archetype_eligible]
archetype_assignments['level_shift'] = gene_level_shift[archetype_eligible]
archetype_assignments['transition'] = archetype_transition_labels
# The existing module labels are joined by gene name from the sections that own them; they are carried
# beside the archetypes, never re-derived, and they never influence an archetype.
archetype_assignments['positional_module'] = [positional_by_gene[gene] for gene in archetype_genes]
archetype_assignments['response_module'] = [response_by_gene[gene] for gene in archetype_genes]
archetype_assignments['human_independent_module'] = [
    human_independent_by_gene[gene] for gene in archetype_genes]
archetype_assignments = archetype_assignments[[
    'gene', 'mouse_archetype', 'human_archetype', 'shape_correlation', 'shape_rms_difference',
    'positional_module', 'response_module', 'archetype_status', 'transition', 'shape_spearman',
    'shape_max_abs_difference', 'human_independent_module', 'mouse_amplitude', 'human_amplitude',
    'response_peak_to_peak', 'response_rms', 'level_shift',
]]
if not (archetype_assignments['gene'].to_numpy() == archetype_genes).all():
    raise ValueError('The archetype assignment table is not aligned with the eligible gene order.')
archetype_assignments.to_csv(OUTPUT_DIR / 'shared_archetype_assignments.csv', index=False)

n_conserved = int((archetype_assignments['archetype_status'] == 'conserved').sum())
n_switched = int((archetype_assignments['archetype_status'] == 'switched').sum())
print(f'Genes carrying both a mouse and a human archetype: {len(archetype_assignments):,}')
print(f'  conserved archetype (mouse_archetype == human_archetype): {n_conserved:,} '
      f'({n_conserved / len(archetype_assignments):.1%})')
print(f'  switched archetype (mouse_archetype != human_archetype): {n_switched:,} '
      f'({n_switched / len(archetype_assignments):.1%})')
print('Conserved and switched are DESCRIPTIVE labels on two cluster assignments. They carry no '
      'p-value, no species-level inference, and the human side is two sections of one donor.')


# %% [markdown]
# ### 12.6 - the program-transition matrix

# %%
# Purpose: 12.6 - mouse archetype -> human archetype transitions: counts, row fractions, O/E and
# standardised residuals.
transition_counts = pd.crosstab(
    pd.Series(archetype_mouse_label, name='mouse_archetype'),
    pd.Series(archetype_human_label, name='human_archetype'),
).reindex(index=shared_archetypes, columns=shared_archetypes, fill_value=0)
transition_total = float(transition_counts.to_numpy().sum())
transition_expected = np.outer(
    transition_counts.sum(axis=1).to_numpy(dtype=float),
    transition_counts.sum(axis=0).to_numpy(dtype=float),
) / transition_total
transition_row_fractions = (
    transition_counts.to_numpy(dtype=float)
    / np.where(transition_counts.sum(axis=1).to_numpy(dtype=float)[:, None] > 0,
               transition_counts.sum(axis=1).to_numpy(dtype=float)[:, None], np.nan))

transition_rows = pd.DataFrame({
    'mouse_archetype': np.repeat(transition_counts.index.to_numpy(), transition_counts.shape[1]),
    'human_archetype': np.tile(transition_counts.columns.to_numpy(), transition_counts.shape[0]),
    'n_genes': transition_counts.to_numpy().ravel(),
    'row_fraction': transition_row_fractions.ravel(),
    'expected': transition_expected.ravel(),
})
transition_rows['observed_over_expected'] = (
    transition_rows['n_genes'] / transition_rows['expected'].replace(0.0, np.nan))
transition_rows['standardized_residual'] = (
    (transition_rows['n_genes'] - transition_rows['expected'])
    / np.sqrt(transition_rows['expected'].replace(0.0, np.nan)))
transition_rows['conserved'] = transition_rows['mouse_archetype'].eq(
    transition_rows['human_archetype'])
transition_rows.to_csv(OUTPUT_DIR / 'mouse_human_archetype_transitions.csv', index=False)

transition_diagonal = transition_rows[transition_rows['conserved']]
_off_diagonal = transition_rows[~transition_rows['conserved']].sort_values(
    'n_genes', ascending=False)
print(f'Transition matrix: {transition_counts.shape[0]} mouse archetypes x '
      f'{transition_counts.shape[1]} human archetypes over {int(transition_total):,} genes.')
print(f'  conserved (diagonal) genes: {int(transition_diagonal["n_genes"].sum()):,} '
      f'({int(transition_diagonal["n_genes"].sum()) / transition_total:.1%}); '
      f'diagonal O/E median {transition_diagonal["observed_over_expected"].median():.2f} '
      f'(range {transition_diagonal["observed_over_expected"].min():.2f}-'
      f'{transition_diagonal["observed_over_expected"].max():.2f})')
print('Most common off-diagonal transitions (the switching programs):')
display(_off_diagonal.head(8)[[
    'mouse_archetype', 'human_archetype', 'n_genes', 'row_fraction', 'expected',
    'observed_over_expected', 'standardized_residual',
]].round(3))
print('O/E compares each cell with the count its row and column margins alone would predict, so a '
      'large archetype does not dominate the reading. Standardised residuals are the cell-level '
      'counterpart and are reported with the counts rather than instead of them.')


# %%
# Purpose: figure 13 - the mouse -> human transition heatmap, on both the raw-count and the O/E scale.
transition_ratio = transition_rows.pivot(
    index='mouse_archetype', columns='human_archetype', values='observed_over_expected').reindex(
    index=shared_archetypes, columns=shared_archetypes)
transition_fraction = transition_rows.pivot(
    index='mouse_archetype', columns='human_archetype', values='row_fraction').reindex(
    index=shared_archetypes, columns=shared_archetypes)
_transition_counts_values = transition_counts.to_numpy(dtype=float)
_transition_ratio_values = transition_ratio.to_numpy(dtype=float)
_transition_fraction_values = transition_fraction.to_numpy(dtype=float)

fig, axes = plt.subplots(
    1, 2, figsize=(2.1 * transition_counts.shape[1] + 5.4, 0.74 * transition_counts.shape[0] + 3.2))
_counts_image = axes[0].imshow(_transition_counts_values, cmap='Blues', vmin=0)
for row in range(transition_counts.shape[0]):
    for column in range(transition_counts.shape[1]):
        axes[0].text(column, row, f'{int(_transition_counts_values[row, column])}',
                     ha='center', va='center', fontsize=8,
                     color='white' if _transition_counts_values[row, column]
                     > 0.55 * np.nanmax(_transition_counts_values) else 'black')
axes[0].set_title('Raw counts')
fig.colorbar(_counts_image, ax=axes[0], label='genes', shrink=0.7)
_ratio_image = axes[1].imshow(np.clip(np.nan_to_num(_transition_ratio_values, nan=0.0), 0, 2),
                              cmap='RdBu_r', vmin=0, vmax=2)
for row in range(transition_counts.shape[0]):
    for column in range(transition_counts.shape[1]):
        ratio = _transition_ratio_values[row, column]
        ratio_text = 'n/a' if not np.isfinite(ratio) else f'{ratio:.2f}'
        axes[1].text(column, row,
                     f'{_transition_fraction_values[row, column]:.2f}\nO/E {ratio_text}',
                     ha='center', va='center', fontsize=7,
                     color='white' if abs(np.nan_to_num(ratio, nan=1.0) - 1) > 0.45 else 'black')
axes[1].set_title('Row-normalised fraction / observed-over-expected')
fig.colorbar(_ratio_image, ax=axes[1], label='observed / expected', shrink=0.7)
for axis in axes:
    axis.set_xticks(range(transition_counts.shape[1]))
    axis.set_xticklabels(transition_counts.columns, rotation=45, ha='right')
    axis.set_yticks(range(transition_counts.shape[0]))
    axis.set_yticklabels(transition_counts.index)
    axis.set_xlabel('human archetype')
    axis.set_ylabel('mouse archetype')
fig.suptitle('Figure 13 - which positional programs are conserved, and which transitions are enriched',
             fontsize=12)
_save_figure(fig, 'fig13_archetype_transition_heatmap.png')


# %% [markdown]
# ### 12.7 - conservation at the archetype level

# %%
# Purpose: 12.7 - what each archetype retains, and what leaves it.
response_module_eligible = archetype_assignments['response_module'].to_numpy()
response_peak_to_peak_eligible = gene_response_peak_to_peak[archetype_eligible]
conservation_rows = []
for archetype in shared_archetypes:
    mouse_members = archetype_mouse_label == archetype
    human_members = archetype_human_label == archetype
    retained = mouse_members & human_members
    switched_out = mouse_members & ~human_members
    module_counts = pd.Series(response_module_eligible[mouse_members]).value_counts()
    conservation_rows.append({
        'archetype': archetype,
        'n_mouse_members': int(mouse_members.sum()),
        'n_human_members': int(human_members.sum()),
        'n_conserved_genes': int(retained.sum()),
        'fraction_mouse_members_retained': (float(retained.sum() / mouse_members.sum())
                                            if mouse_members.sum() else np.nan),
        'fraction_human_members_from_same_mouse_archetype': (
            float(retained.sum() / human_members.sum()) if human_members.sum() else np.nan),
        'median_shape_correlation_conserved': _median_or_nan(
            archetype_shape_correlation[retained]),
        'median_shape_correlation_switched_out': _median_or_nan(
            archetype_shape_correlation[switched_out]),
        'median_response_peak_to_peak': _median_or_nan(
            response_peak_to_peak_eligible[mouse_members]),
        'dominant_response_module': module_counts.index[0] if len(module_counts) else '',
        'fraction_in_dominant_response_module': (float(module_counts.iloc[0] / module_counts.sum())
                                                 if len(module_counts) else np.nan),
        'response_module_distribution': '; '.join(
            f'{module}:{count}' for module, count in module_counts.items()),
    })
archetype_conservation = pd.DataFrame(conservation_rows)
archetype_catalog = archetype_catalog.merge(archetype_conservation, on='archetype', how='left')
archetype_catalog.to_csv(OUTPUT_DIR / 'shared_archetype_catalog.csv', index=False)

display(archetype_catalog[[
    'archetype', 'n_mouse_members', 'n_human_members', 'n_conserved_genes',
    'fraction_mouse_members_retained', 'fraction_human_members_from_same_mouse_archetype',
    'median_shape_correlation_conserved', 'median_shape_correlation_switched_out',
    'median_response_peak_to_peak', 'dominant_response_module',
    'fraction_in_dominant_response_module',
]].round(3))
_structurally_shared = archetype_catalog['centroid_spearman_mouse_vs_human'] >= 0.8
_membership_rewired = archetype_catalog['fraction_mouse_members_retained'] < 0.5
print(f'Archetypes whose mouse and human centroids agree (Spearman >= 0.8): '
      f'{int(_structurally_shared.sum())} of {len(archetype_catalog)}')
print(f'Archetypes shared structurally but retaining < 50% of their mouse membership: '
      f'{int((_structurally_shared & _membership_rewired).sum())}')
print('A median conserved correlation far above a median switched-out correlation means the label '
      '"conserved" tracks a genuinely conserved shape rather than a convenience of the clustering. '
      'Where the two are similar, membership is changing without the shape changing much - or the '
      'assignment itself is close to a boundary, and the sensitivity analyses in 12.12/12.13 test it.')


# %% [markdown]
# ### 12.8 - the shared solution against the existing independent clustering
#
# Section 10C clusters human standardised curves independently of mouse and reports a low adjusted
# Rand index together with centroid matching. This section does not delete that: it asks the same
# question against the shared solution, which is the one that could resolve the label-matching problem.

# %%
# Purpose: 12.8 - ARI against the independent clusterings, and centroid matching against P and H.
positional_module_eligible = positional_labels.to_numpy()[archetype_eligible]
human_independent_eligible = human_labels.to_numpy()[archetype_eligible]
_mouse_assigned = positional_module_eligible != 'unassigned'
_human_assigned = human_independent_eligible != 'unassigned'
archetype_vs_independent = pd.DataFrame([
    {
        'comparison': 'mouse: independent positional modules vs shared mouse archetypes',
        'n_genes': int(_mouse_assigned.sum()),
        'adjusted_rand_index': float(adjusted_rand_score(
            positional_module_eligible[_mouse_assigned],
            archetype_mouse_label[_mouse_assigned])),
    },
    {
        'comparison': 'human: independent positional modules vs shared human archetypes',
        'n_genes': int(_human_assigned.sum()),
        'adjusted_rand_index': float(adjusted_rand_score(
            human_independent_eligible[_human_assigned],
            archetype_human_label[_human_assigned])),
    },
])
archetype_vs_independent.to_csv(
    DIAGNOSTIC_DIR / 'shared_archetype_vs_independent_clustering.csv', index=False)

# Each centroid is the mean of its members' standardised curves, and the three clusterings were
# standardised over their own support, so every comparison is made on the grid columns all three
# cover. Spearman is used because it reads the SHAPE of the two centroids, not their offset.
comparison_columns = positional_columns & human_columns & archetype_columns
if int(comparison_columns.sum()) < 3:
    comparison_columns = archetype_columns
archetype_subset = comparison_columns[archetype_columns]
positional_subset = comparison_columns[positional_columns]
human_subset = comparison_columns[human_columns]

independent_centroid_matching = []
for module in positional['catalog']['module']:
    members = positional_labels.to_numpy() == module
    if not members[positional_eligible].any():
        continue
    centroid = np.nanmean(
        positional['standardized'][members[positional_eligible]][:, positional_subset], axis=0)
    module_genes = set(gene_names[members])
    for archetype in shared_archetypes:
        archetype_members = archetype_mouse_label == archetype
        independent_centroid_matching.append({
            'independent_species': 'mouse',
            'independent_module': module,
            'archetype': archetype,
            'centroid_spearman': safe_spearman(
                centroid, shared['centroid_mouse'][archetype][archetype_subset]),
            'member_jaccard': _jaccard(module_genes, archetype_genes[archetype_members]),
        })
for module in sorted(set(human_result['labels']), key=lambda name: int(name[1:])):
    members = human_labels.to_numpy() == module
    if not members[human_eligible].any():
        continue
    centroid = np.nanmean(
        human_result['standardized'][members[human_eligible]][:, human_subset], axis=0)
    module_genes = set(gene_names[members])
    for archetype in shared_archetypes:
        archetype_members = archetype_human_label == archetype
        independent_centroid_matching.append({
            'independent_species': 'human',
            'independent_module': module,
            'archetype': archetype,
            'centroid_spearman': safe_spearman(
                centroid, shared['centroid_human'][archetype][archetype_subset]),
            'member_jaccard': _jaccard(module_genes, archetype_genes[archetype_members]),
        })
independent_centroid_matching = pd.DataFrame(independent_centroid_matching)
independent_centroid_matching.to_csv(
    DIAGNOSTIC_DIR / 'shared_archetype_vs_independent_centroids.csv', index=False)
best_independent_match = (
    independent_centroid_matching.sort_values('centroid_spearman', ascending=False)
    .groupby(['independent_species', 'independent_module'], observed=True).head(1)
    .reset_index(drop=True))
display(archetype_vs_independent.round(3))
display(best_independent_match.round(3))
for species_name in ('mouse', 'human'):
    _one = best_independent_match[best_independent_match['independent_species'].eq(species_name)]
    if len(_one):
        print(f'{species_name}: best centroid match per independent module, median Spearman '
              f'{_one["centroid_spearman"].median():.3f} '
              f'(min {_one["centroid_spearman"].min():.3f}); '
              f'median member Jaccard {_one["member_jaccard"].median():.3f}')
print('A low ARI alone is NOT the finding - two clusterings of the same curves can agree on shapes '
      'while disagreeing on labels, which is precisely section 10C\'s warning. Read the centroid '
      'correlations beside it, and remember that a high centroid correlation with a low Jaccard is '
      'the shared-archetype result, not a failure to reproduce it.')


# %% [markdown]
# ### 12.9 - what the existing response modules represent

# %%
# Purpose: 12.9 - re-read every response module in archetype-transition terms.
response_module_archetype_rows = []
for module in response_catalog_stats['module']:
    members = response_module_eligible == module
    if not members.any():
        continue
    mouse_distribution = pd.Series(archetype_mouse_label[members]).value_counts()
    human_distribution = pd.Series(archetype_human_label[members]).value_counts()
    transition_distribution = pd.Series(archetype_transition_labels[members]).value_counts()
    response_module_archetype_rows.append({
        'response_module': module,
        'n_genes': int(members.sum()),
        'dominant_mouse_archetype': mouse_distribution.index[0],
        'fraction_dominant_mouse_archetype': float(
            mouse_distribution.iloc[0] / mouse_distribution.sum()),
        'dominant_human_archetype': human_distribution.index[0],
        'fraction_dominant_human_archetype': float(
            human_distribution.iloc[0] / human_distribution.sum()),
        'most_common_transition': transition_distribution.index[0],
        'fraction_most_common_transition': float(
            transition_distribution.iloc[0] / transition_distribution.sum()),
        'fraction_switched': float(np.mean(
            archetype_mouse_label[members] != archetype_human_label[members])),
        'median_shape_correlation': _median_or_nan(archetype_shape_correlation[members]),
        'mouse_archetype_distribution': '; '.join(
            f'{archetype}:{count}' for archetype, count in mouse_distribution.items()),
        'human_archetype_distribution': '; '.join(
            f'{archetype}:{count}' for archetype, count in human_distribution.items()),
        'transition_distribution': '; '.join(
            f'{label}:{count}' for label, count in transition_distribution.items()),
    })
response_module_archetype_summary = pd.DataFrame(response_module_archetype_rows)
response_module_archetype_summary.to_csv(
    OUTPUT_DIR / 'response_module_archetype_summary.csv', index=False)
display(response_module_archetype_summary[[
    'response_module', 'n_genes', 'dominant_mouse_archetype', 'dominant_human_archetype',
    'most_common_transition', 'fraction_most_common_transition', 'fraction_switched',
    'median_shape_correlation',
]].round(3))
_dominant_transition_modules = response_module_archetype_summary[
    response_module_archetype_summary['fraction_most_common_transition'] >= 0.5]
print(f'Response modules whose members mostly share ONE positional transition (>= 50%): '
      f'{len(_dominant_transition_modules)} of {len(response_module_archetype_summary)}')
print('Where that fraction is high, the module is better described as a positional-program switch than '
      'as a generic difference-curve phenotype. Where it is low, the module mixes several programs and '
      'should not be given a single positional interpretation.')


# %%
# Purpose: figure 14 - response module x archetype transition, as a within-module fraction.
_transition_categories = [f'{mouse}->{human}' for mouse in shared_archetypes
                          for human in shared_archetypes]
response_transition_counts = pd.crosstab(
    pd.Series(response_module_eligible, name='response_module'),
    pd.Series(archetype_transition_labels, name='transition'),
).reindex(columns=_transition_categories, fill_value=0)
# Row order read off the module NUMBER rather than as text, so R10 follows R9 instead of R1; the
# genes with no response module are kept, and last, because they are part of the picture.
response_transition_counts = response_transition_counts.reindex(sorted(
    response_transition_counts.index,
    key=lambda name: (name == 'unassigned', int(name[1:]) if name[1:].isdigit() else 0)))
_response_transition_values = response_transition_counts.to_numpy(dtype=float)
_response_row_totals = _response_transition_values.sum(axis=1, keepdims=True)
_response_transition_fraction = _response_transition_values / np.where(
    _response_row_totals > 0, _response_row_totals, np.nan)

fig, axis = plt.subplots(
    figsize=(0.42 * len(_transition_categories) + 4.2, 0.62 * len(response_transition_counts) + 3))
_columns_to_show = np.flatnonzero(np.nan_to_num(_response_transition_fraction, nan=0.0).max(axis=0) >= 0.05)
if _columns_to_show.size == 0:
    _columns_to_show = np.arange(len(_transition_categories))
_shown = _response_transition_fraction[:, _columns_to_show]
_image = axis.imshow(_shown, cmap='Blues', vmin=0, vmax=max(0.05, float(np.nanmax(_shown))))
for row in range(_shown.shape[0]):
    for column in range(_shown.shape[1]):
        if np.isfinite(_shown[row, column]) and _shown[row, column] >= 0.05:
            axis.text(column, row, f'{_shown[row, column]:.2f}', ha='center', va='center',
                      fontsize=7, color='white' if _shown[row, column] > 0.55 * np.nanmax(_shown)
                      else 'black')
axis.set_xticks(range(_columns_to_show.size))
axis.set_xticklabels([f'{_transition_categories[position]}*' if
                      _transition_categories[position].split('->')[0]
                      == _transition_categories[position].split('->')[1]
                      else _transition_categories[position] for position in _columns_to_show],
                     rotation=60, ha='right', fontsize=7)
axis.set_yticks(range(_response_transition_fraction.shape[0]))
axis.set_yticklabels(response_transition_counts.index)
axis.set_xlabel('mouse archetype -> human archetype (* = conserved transition)')
axis.set_ylabel('existing response module (section 6)')
fig.colorbar(_image, ax=axis, label='fraction of the module', shrink=0.7)
axis.set_title('Figure 14 - do the existing response modules correspond to positional-program '
               'transitions?\ncolumns below 5% in every module are omitted; cell text = fraction',
               fontsize=11)
_save_figure(fig, 'fig14_response_module_transitions.png')


# %% [markdown]
# ### 12.10 - pathway enrichment of the transitions
#
# The archetypes and the transition matrix are frozen by now, so this is the first place pathway
# information is allowed in. No pathway expression is aggregated: every test is an over-representation
# test of a gene set defined by a transition or by conservation, against the genes that could actually
# have been assigned an archetype.

# %%
# Purpose: 12.10 - the transition and conservation gene sets, with their sizes and skip reasons.
transition_sets = {}
transition_set_meta = []
for archetype in shared_archetypes:
    mouse_members = archetype_mouse_label == archetype
    for set_type, members, label in (
        ('conserved within archetype', mouse_members & (archetype_human_label == archetype),
         f'conserved in {archetype}'),
        ('switched out of archetype', mouse_members & (archetype_human_label != archetype),
         f'switched out of {archetype}'),
    ):
        testable = int(members.sum()) >= ARCHETYPE_CONFIG['pathway_min_transition_genes']
        transition_set_meta.append({
            'transition_set': label, 'set_type': set_type, 'mouse_archetype': archetype,
            'human_archetype': archetype if set_type == 'conserved within archetype' else '',
            'n_genes': int(members.sum()), 'tested': bool(testable),
            'skipped_reason': '' if testable else
            f'below the {ARCHETYPE_CONFIG["pathway_min_transition_genes"]} gene minimum',
        })
        if testable:
            transition_sets[label] = list(archetype_genes[members])
for mouse_archetype in shared_archetypes:
    for human_archetype in shared_archetypes:
        members = (archetype_mouse_label == mouse_archetype) & (
            archetype_human_label == human_archetype)
        label = f'{mouse_archetype}->{human_archetype}'
        testable = int(members.sum()) >= ARCHETYPE_CONFIG['pathway_min_transition_genes']
        transition_set_meta.append({
            'transition_set': label, 'set_type': 'individual transition',
            'mouse_archetype': mouse_archetype, 'human_archetype': human_archetype,
            'n_genes': int(members.sum()), 'tested': bool(testable),
            'skipped_reason': '' if testable else
            f'below the {ARCHETYPE_CONFIG["pathway_min_transition_genes"]} gene minimum',
        })
        if testable:
            transition_sets[label] = list(archetype_genes[members])
transition_set_meta = pd.DataFrame(transition_set_meta)
transition_set_meta.to_csv(DIAGNOSTIC_DIR / 'shared_archetype_transition_set_meta.csv', index=False)
archetype_background = archetype_genes
print(f'Transition and conservation sets defined: {len(transition_set_meta)}; tested: '
      f'{len(transition_sets)}; skipped below the size minimum: '
      f'{int((~transition_set_meta["tested"]).sum())}')
print(f'Enrichment universe: the {archetype_background.size:,} genes eligible for the shared-archetype '
      'analysis - not every expressed gene.')

# Both branches bind `transition_enrichment`, so a run in which nothing is testable still produces the
# table with the schema below rather than a later NameError.
_transition_enrichment_columns = [
    'module', 'gene_set', 'n_module_genes', 'n_set_genes_in_background', 'n_overlap', 'overlap_genes',
    'expected_overlap', 'p_value', 'p_value_adjusted', 'library', 'pathway', 'fold_enrichment',
    'redundancy_group', 'group_size',
]
if transition_sets:
    transition_enrichment = cached_frame(
        'pt_archetype_transition_pathway_enrichment',
        lambda: enrich_modules(transition_sets, pathway_gene_sets, background=archetype_background),
        root=STAGE_CACHE_DIR,
        params={'logic': 1, 'min_genes': PATHWAY_CONFIG['min_genes'],
                'min_transition_genes': ARCHETYPE_CONFIG['pathway_min_transition_genes']},
        inputs={'sets': digest(sorted((name, tuple(members))
                                      for name, members in transition_sets.items())),
                'background': digest(np.asarray(archetype_background)),
                'pathway_sets': digest(sorted((name, tuple(members))
                                              for name, members in pathway_gene_sets.items()))},
        code=code_digest(enrich_modules), enabled=STAGE_CACHE_ENABLED,
    )
    transition_enrichment = _annotate_enrichment(transition_enrichment, retained_pathways)
else:
    transition_enrichment = pd.DataFrame(columns=_transition_enrichment_columns)
    print('No transition or conservation set reaches the minimum size, so no pathway test was run.')

transition_enrichment_out = transition_enrichment.rename(columns={
    'module': 'transition_set', 'n_module_genes': 'set_size',
    'n_set_genes_in_background': 'n_pathway_genes_in_universe',
    'fold_enrichment': 'observed_over_expected'}).merge(
    transition_set_meta[['transition_set', 'set_type', 'mouse_archetype', 'human_archetype']],
    on='transition_set', how='left')
transition_enrichment_out = transition_enrichment_out[[
    'transition_set', 'set_type', 'mouse_archetype', 'human_archetype', 'library', 'pathway',
    'n_overlap', 'set_size', 'n_pathway_genes_in_universe', 'expected_overlap',
    'observed_over_expected', 'p_value', 'p_value_adjusted', 'redundancy_group', 'group_size',
    'overlap_genes',
]]
transition_enrichment_out.to_csv(OUTPUT_DIR / 'archetype_transition_enrichment.csv', index=False)
if not transition_enrichment_out.empty:
    transition_enrichment_significant = transition_enrichment_out[
        transition_enrichment_out['p_value_adjusted'] < 0.05]
    print(f'Tested {len(transition_enrichment_out):,} transition x pathway pairs; '
          f'{len(transition_enrichment_significant):,} at BH q < 0.05.')
    display(transition_enrichment_significant.sort_values(
        'p_value_adjusted').head(15)[[
        'transition_set', 'set_type', 'library', 'pathway', 'n_overlap', 'set_size',
        'n_pathway_genes_in_universe', 'observed_over_expected', 'p_value', 'p_value_adjusted',
        'redundancy_group',
    ]].round(4))
else:
    transition_enrichment_significant = transition_enrichment_out
    print('No transition set could be tested, so there is no pathway result in this section.')
print('Correction is BH across this whole family of transition x pathway pairs, pooled over every '
      'tested set. These are exploratory set-level results from a 2-vs-1-donor comparison, not '
      'confirmatory claims, and they never revise the archetypes or the transitions.')


# %% [markdown]
# ### 12.11 - pathway coherence versus switching
#
# A pathway is not a program. Its members can sit in one positional archetype, spread across several,
# or make one dominant transition. Raw entropy cannot be read on its own - a 20-gene pathway looks
# concentrated whether or not it is - so every pathway statistic is compared with an empirical null of
# random gene sets matched for size AND for approximate mouse amplitude.

# %%
# Purpose: 12.11 - archetype entropy, transition entropy, conserved fraction and direct curve
# coherence per pathway, each against a size- and amplitude-matched null.
def _normalized_entropy(counts, n_categories):
    counts = np.asarray(counts, dtype=float)
    total = counts.sum()
    if total <= 0 or n_categories <= 1:
        return 0.0
    proportions = counts[counts > 0] / total
    return float(-(proportions * np.log(proportions)).sum() / np.log(n_categories))


_archetype_index = {archetype: position for position, archetype in enumerate(shared_archetypes)}
_n_archetypes = len(shared_archetypes)
_mouse_code = np.array([_archetype_index[value] for value in archetype_mouse_label], dtype=int)
_human_code = np.array([_archetype_index[value] for value in archetype_human_label], dtype=int)
_transition_code = _mouse_code * _n_archetypes + _human_code
_archetype_gene_position = {gene: position for position, gene in enumerate(archetype_genes)}


def _pathway_archetype_statistics(positions):
    """Where one gene set sits: archetype entropy per species, transition entropy, and coherence."""
    positions = np.asarray(positions, dtype=int)
    if positions.size == 0:
        return {'mouse_archetype_entropy': np.nan, 'human_archetype_entropy': np.nan,
                'transition_entropy': np.nan, 'fraction_conserved': np.nan,
                'dominant_transition': '', 'fraction_in_dominant_transition': np.nan,
                'median_shape_correlation': np.nan}
    mouse_counts = np.bincount(_mouse_code[positions], minlength=_n_archetypes)
    human_counts = np.bincount(_human_code[positions], minlength=_n_archetypes)
    transition_count = np.bincount(_transition_code[positions], minlength=_n_archetypes ** 2)
    dominant = int(np.argmax(transition_count))
    return {
        'mouse_archetype_entropy': _normalized_entropy(mouse_counts, _n_archetypes),
        'human_archetype_entropy': _normalized_entropy(human_counts, _n_archetypes),
        'transition_entropy': _normalized_entropy(transition_count, _n_archetypes ** 2),
        'fraction_conserved': float(np.mean(_mouse_code[positions] == _human_code[positions])),
        'dominant_transition': f'{shared_archetypes[dominant // _n_archetypes]}'
                               f'->{shared_archetypes[dominant % _n_archetypes]}',
        'fraction_in_dominant_transition': float(transition_count[dominant] / positions.size),
        'median_shape_correlation': _median_or_nan(archetype_shape_correlation[positions]),
    }


# Amplitude deciles over the eligible genes, so a null set can be drawn from the same amplitude range
# as the pathway it stands in for.
# `qcut` returns an ndarray rather than a Series when it is handed an ndarray, and a value landing on
# a dropped duplicate edge comes back as NaN; it is folded into the lowest decile, not discarded.
_amplitude_deciles = np.nan_to_num(np.asarray(
    pd.qcut(gene_mouse_amplitude[archetype_eligible], 10, labels=False, duplicates='drop'),
    dtype=float), nan=0.0).astype(int)
_decile_pool = {int(decile): np.flatnonzero(_amplitude_deciles == decile)
                for decile in np.unique(_amplitude_deciles)}
_archetype_rng = np.random.default_rng(ARCHETYPE_CONFIG['seed'])

coherence_rows = []
for pathway_name, pathway_members in pathway_gene_sets.items():
    pathway_positions = np.array(
        [_archetype_gene_position[gene] for gene in pathway_members
         if gene in _archetype_gene_position], dtype=int)
    if pathway_positions.size < ARCHETYPE_CONFIG['coherence_min_pathway_genes']:
        continue
    observed = _pathway_archetype_statistics(pathway_positions)
    pathway_deciles = _amplitude_deciles[pathway_positions]
    null_transition_entropy, null_correlation, null_conserved = [], [], []
    for _ in range(ARCHETYPE_CONFIG['null_permutations']):
        sampled = []
        for decile in np.unique(pathway_deciles):
            pool = _decile_pool[int(decile)]
            count = int((pathway_deciles == decile).sum())
            sampled.append(_archetype_rng.choice(pool, size=count, replace=count > pool.size))
        null_statistics = _pathway_archetype_statistics(
            np.concatenate(sampled) if sampled else np.array([], dtype=int))
        null_transition_entropy.append(null_statistics['transition_entropy'])
        null_correlation.append(null_statistics['median_shape_correlation'])
        null_conserved.append(null_statistics['fraction_conserved'])
    null_transition_entropy = np.asarray(null_transition_entropy, dtype=float)
    null_correlation = np.asarray(null_correlation, dtype=float)
    row = {
        'pathway': pathway_name,
        'eligible_gene_count': int(pathway_positions.size),
        **observed,
        'null_mean_transition_entropy': float(np.nanmean(null_transition_entropy)),
        'null_sd_transition_entropy': float(np.nanstd(null_transition_entropy, ddof=1)),
        'null_mean_shape_correlation': float(np.nanmean(null_correlation)),
        'null_sd_shape_correlation': float(np.nanstd(null_correlation, ddof=1)),
        'null_mean_fraction_conserved': float(np.nanmean(null_conserved)),
    }
    # Each z-score pairs an observed statistic with the null of the SAME statistic, so the observed
    # column and the null column are named explicitly rather than derived from one key.
    for observed_key, null_key in (('transition_entropy', 'transition_entropy'),
                                   ('median_shape_correlation', 'shape_correlation')):
        spread = row[f'null_sd_{null_key}']
        row[f'z_{null_key}'] = (
            float((row[observed_key] - row[f'null_mean_{null_key}']) / spread)
            if np.isfinite(spread) and spread > 0 else np.nan)
    row['p_low_transition_entropy'] = float(
        (1 + int((null_transition_entropy <= observed['transition_entropy']).sum()))
        / (null_transition_entropy.size + 1))
    row['p_high_shape_correlation'] = float(
        (1 + int((null_correlation >= observed['median_shape_correlation']).sum()))
        / (null_correlation.size + 1))
    coherence_rows.append(row)

pathway_archetype_coherence = pd.DataFrame(
    coherence_rows,
    columns=[
        'pathway', 'eligible_gene_count', 'fraction_conserved', 'dominant_transition',
        'fraction_in_dominant_transition', 'mouse_archetype_entropy', 'human_archetype_entropy',
        'transition_entropy', 'median_shape_correlation', 'null_mean_transition_entropy',
        'null_sd_transition_entropy', 'z_transition_entropy', 'null_mean_shape_correlation',
        'null_sd_shape_correlation', 'z_shape_correlation', 'null_mean_fraction_conserved',
        'p_low_transition_entropy', 'p_high_shape_correlation',
    ]).sort_values(['eligible_gene_count'], ascending=False).reset_index(drop=True)
pathway_archetype_coherence.to_csv(OUTPUT_DIR / 'pathway_archetype_coherence.csv', index=False)

if pathway_archetype_coherence.empty:
    print(f'No pathway had >= {ARCHETYPE_CONFIG["coherence_min_pathway_genes"]} eligible genes, so the '
          'coherence analysis is empty.')
else:
    display(pathway_archetype_coherence.head(12)[[
        'pathway', 'eligible_gene_count', 'fraction_conserved', 'dominant_transition',
        'fraction_in_dominant_transition', 'mouse_archetype_entropy', 'human_archetype_entropy',
        'transition_entropy', 'median_shape_correlation', 'z_transition_entropy',
        'z_shape_correlation', 'p_low_transition_entropy',
    ]].round(3))
    print(f'Pathways assessed: {len(pathway_archetype_coherence):,} '
          f'(>= {ARCHETYPE_CONFIG["coherence_min_pathway_genes"]} eligible genes each).')
    print('A high transition entropy against its matched null means the pathway\'s genes do NOT make '
          'one positional transition: they split across programs, which is a case an aggregate '
          'pathway curve would average away. A low entropy means one dominant transition.')


# %%
# Purpose: figure 15 - pathway coherence: concentration against the matched null, and the extremes.
if pathway_archetype_coherence.empty:
    print('Figure 15 skipped: no pathway reached the eligible-gene minimum.')
else:
    _ranked_by_heterogeneity = pathway_archetype_coherence.dropna(
        subset=['z_transition_entropy']).sort_values('z_transition_entropy', ascending=False)
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 4.6),
                             gridspec_kw={'width_ratios': [1, 1.15]})
    _scatter = axes[0].scatter(
        pathway_archetype_coherence['transition_entropy'],
        pathway_archetype_coherence['fraction_conserved'],
        s=18 + 2.4 * pathway_archetype_coherence['eligible_gene_count'],
        c=pathway_archetype_coherence['median_shape_correlation'], cmap='coolwarm',
        vmin=-1, vmax=1, edgecolor='k', linewidth=0.3)
    axes[0].set_xlabel('normalised transition entropy')
    axes[0].set_ylabel('fraction of the pathway conserved')
    axes[0].set_title('One program, or several?')
    fig.colorbar(_scatter, ax=axes[0], label='median shape correlation', shrink=0.7)
    for _, row in pd.concat([_ranked_by_heterogeneity.head(5),
                             _ranked_by_heterogeneity.tail(5)]).iterrows():
        axes[0].annotate(row['pathway'].split(': ')[-1][:26],
                         (row['transition_entropy'], row['fraction_conserved']), fontsize=6)
    _extremes = pd.concat([_ranked_by_heterogeneity.head(8),
                           _ranked_by_heterogeneity.tail(8)]).drop_duplicates('pathway')
    _extremes = _extremes.sort_values('z_transition_entropy')
    axes[1].barh(
        [name.split(': ')[-1][:34] for name in _extremes['pathway']],
        _extremes['z_transition_entropy'],
        color=np.where(_extremes['z_transition_entropy'] > 0, '#E15759', '#4C9BD3'))
    axes[1].axvline(0, color='k', lw=0.8)
    axes[1].set_xlabel('transition entropy z-score vs matched null')
    axes[1].set_title('Pathways most split (>0) and most coherent (<0)')
    axes[1].tick_params(axis='y', labelsize=7)
    fig.suptitle('Figure 15 - pathway coherence versus positional-program switching '
                 '(null matched for size and mouse amplitude)', fontsize=12)
    _save_figure(fig, 'fig15_pathway_coherence.png')


# %% [markdown]
# ### 12.12 - threshold sensitivity
#
# The current response analysis contains many genes with small response amplitudes, so the question is
# whether the archetypes and the transitions survive a stricter amplitude floor. Only the cheap part is
# repeated - eligibility, the joint PCA and the single clustering. No downstream enrichment is rerun at
# every rung, and every rung keeps one fixed definition of what a transition is by matching solutions on
# centroid SHAPE rather than on the arbitrary cluster names.

# %%
# Purpose: 12.12 - repeat the shared-archetype solution at 1x, 2x and 4x the amplitude floor.
sensitivity_rows = []
sensitivity_persistence_rows = []
_primary_tracked_transitions = {}
for row in transition_rows.sort_values('n_genes', ascending=False).head(
        ARCHETYPE_CONFIG['representative_transitions']).itertuples():
    _primary_tracked_transitions[f'{row.mouse_archetype}->{row.human_archetype}'] = (
        (archetype_mouse_label == row.mouse_archetype)
        & (archetype_human_label == row.human_archetype))

for multiplier in ARCHETYPE_CONFIG['amplitude_floor_ladder']:
    floor = ARCHETYPE_CONFIG['amplitude_floor'] * float(multiplier)
    row = {
        'amplitude_floor_multiplier': float(multiplier),
        'amplitude_floor': float(floor),
        'status': 'evaluated',
        'n_eligible_genes': np.nan,
        'n_archetypes': np.nan,
        'selected_k_best_silhouette': np.nan,
        'smallest_archetype': np.nan,
        'largest_archetype_fraction': np.nan,
        'median_shape_correlation': np.nan,
        'fraction_conserved': np.nan,
        'dominant_transition': '',
        'fraction_in_dominant_transition': np.nan,
        'median_best_centroid_spearman_vs_primary': np.nan,
        'median_fraction_primary_transitions_reproduced': np.nan,
    }
    try:
        threshold_solution = _shared_archetype_solution(
            mouse_curves, human_curves, grid, gene_names, floor, ARCHETYPE_CONFIG['min_genes'])
    except ValueError as error:
        row['status'] = f'not evaluable: {error}'
        sensitivity_rows.append(row)
        continue

    # Match this rung's archetypes to the primary ones on centroid SHAPE. Comparing labels would be
    # meaningless: each run names its own clusters A1..Ak by centroid peak position.
    _rung_common = archetype_columns & threshold_solution['columns']
    if int(_rung_common.sum()) < 3:
        _rung_common = archetype_columns
    _primary_subset = _rung_common[archetype_columns]
    _rung_subset = _rung_common[threshold_solution['columns']]
    _matching = pd.DataFrame([
        {'primary_archetype': primary_archetype, 'rung_archetype': rung_archetype,
         'centroid_spearman': safe_spearman(
             shared['centroid_combined'][primary_archetype][_primary_subset],
             threshold_solution['centroid_combined'][rung_archetype][_rung_subset])}
        for primary_archetype in shared_archetypes
        for rung_archetype in threshold_solution['archetypes']
    ])
    _best_match = (_matching.sort_values('centroid_spearman', ascending=False)
                   .groupby('primary_archetype', observed=True).head(1))
    _rung_to_primary = dict(zip(_best_match['rung_archetype'], _best_match['primary_archetype']))
    row['median_best_centroid_spearman_vs_primary'] = float(_best_match['centroid_spearman'].median())

    _rung_position = {gene: position
                      for position, gene in enumerate(threshold_solution['genes'])}
    _mapped_mouse = np.asarray(
        [_rung_to_primary.get(value, '?') for value in threshold_solution['mouse_archetype']],
        dtype=object)
    _mapped_human = np.asarray(
        [_rung_to_primary.get(value, '?') for value in threshold_solution['human_archetype']],
        dtype=object)
    _reproduced_fractions = []
    for label, members in _primary_tracked_transitions.items():
        mouse_name, human_name = label.split('->')
        positions = np.asarray([_rung_position[gene] for gene in archetype_genes[members]
                                if gene in _rung_position], dtype=int)
        if positions.size == 0:
            sensitivity_persistence_rows.append({
                'amplitude_floor': float(floor), 'transition': label, 'n_genes_available': 0,
                'fraction_reproducing': np.nan})
            continue
        reproduced = ((_mapped_mouse[positions] == mouse_name)
                      & (_mapped_human[positions] == human_name))
        fraction_reproduced = float(reproduced.mean())
        _reproduced_fractions.append(fraction_reproduced)
        sensitivity_persistence_rows.append({
            'amplitude_floor': float(floor), 'transition': label,
            'n_genes_available': int(positions.size),
            'fraction_reproducing': fraction_reproduced})

    _rung_sizes = pd.Series(threshold_solution['mouse_archetype']).value_counts()
    _rung_transitions = pd.Series(np.asarray(
        [f'{mouse}->{human}' for mouse, human in zip(
            threshold_solution['mouse_archetype'], threshold_solution['human_archetype'])],
        dtype=object)).value_counts()
    row.update({
        'n_eligible_genes': int(threshold_solution['eligible'].sum()),
        'n_archetypes': len(threshold_solution['archetypes']),
        'selected_k_best_silhouette': float(threshold_solution['best_silhouette']),
        'smallest_archetype': int(_rung_sizes.min()),
        'largest_archetype_fraction': float(_rung_sizes.max() / _rung_sizes.sum()),
        'median_shape_correlation': float(np.nanmedian(_row_block_pearson(
            threshold_solution['mouse_standardized'], threshold_solution['human_standardized'],
            ARCHETYPE_CONFIG['numerical_epsilon']))),
        'fraction_conserved': float(np.mean(
            np.asarray(threshold_solution['mouse_archetype'])
            == np.asarray(threshold_solution['human_archetype']))),
        # Expressed in the PRIMARY archetype names wherever this rung's archetype matched one, so the
        # column is comparable with the primary transition table instead of silently mixing two
        # clusterings' private labels. An unmatched rung archetype keeps its own name, marked.
        'dominant_transition': '->'.join(
            _rung_to_primary.get(name, f'{name}(unmatched)')
            for name in _rung_transitions.index[0].split('->')),
        'fraction_in_dominant_transition': float(
            _rung_transitions.iloc[0] / _rung_transitions.sum()),
        'median_fraction_primary_transitions_reproduced': (
            float(np.nanmedian(_reproduced_fractions)) if _reproduced_fractions else np.nan),
    })
    sensitivity_rows.append(row)

shared_archetype_threshold_sensitivity = pd.DataFrame(sensitivity_rows)
shared_archetype_threshold_sensitivity.to_csv(
    OUTPUT_DIR / 'shared_archetype_threshold_sensitivity.csv', index=False)
shared_archetype_threshold_persistence = pd.DataFrame(sensitivity_persistence_rows)
shared_archetype_threshold_persistence.to_csv(
    DIAGNOSTIC_DIR / 'shared_archetype_threshold_transition_persistence.csv', index=False)

display(shared_archetype_threshold_sensitivity[[
    'amplitude_floor', 'status', 'n_eligible_genes', 'n_archetypes', 'smallest_archetype',
    'largest_archetype_fraction', 'median_shape_correlation', 'fraction_conserved',
    'dominant_transition', 'fraction_in_dominant_transition',
    'median_best_centroid_spearman_vs_primary',
    'median_fraction_primary_transitions_reproduced',
]].round(3))
display(shared_archetype_threshold_persistence.round(3))
print('The 1x rung is the primary analysis re-derived, so it must reproduce itself exactly (centroid '
      'agreement 1.0 and every tracked transition at 1.0) - it is here as a self-consistency check. '
      'A rung that cannot be evaluated is reported as such rather than dropped: too few genes above the '
      'floor is itself the answer to the sensitivity question.')


# %% [markdown]
# ### 12.13 - specimen-pair consistency of the switching
#
# This is technical and specimen-level robustness only. The two human sections are sections of ONE
# donor: pairs are not biological replicates, and nothing here is a species-level test. Each mouse
# specimen x human section pair is standardised on its own and **projected into the primary shared
# PCA**, then assigned to the nearest primary archetype centroid, so there is one fixed archetype
# definition throughout and no reclustering per pair.

# %%
# Purpose: 12.13 - project every specimen pair into the primary archetype model and score consistency.
primary_raw_scores = shared_pca.transform(shared['standardized'])[:, :shared_n_pc]
primary_whitening_spread = primary_raw_scores.std(axis=0, ddof=1)
if not np.allclose(primary_raw_scores / primary_whitening_spread, shared_scores, atol=1e-8):
    raise ValueError('The whitening recovered from the primary PCA does not reproduce the primary '
                     'scores, so projecting into them would silently use a different space.')
archetype_centroid_scores = pd.DataFrame(
    shared_scores,
    columns=[f'pc{index + 1}_score' for index in range(shared_n_pc)],
).groupby(shared_labels).mean().reindex(shared_archetypes).to_numpy()


def _project_to_primary(curves, columns, pca, whitening, n_pc):
    """One specimen's curves -> standardised, then projected into the primary shared PC space.

    A single specimen's curve is undefined outside its own pseudospace support. Those columns are set
    to the row's own centred mean (0) rather than dropped, because dropping them would move the row
    into a different feature space than the PCA the archetypes were learned in.
    """
    standardized = zscore_rows(np.asarray(curves, dtype=float)[:, columns])
    standardized = np.where(np.isfinite(standardized), standardized, 0.0)
    return standardized, pca.transform(standardized)[:, :n_pc] / whitening


def _nearest_archetype(scores, centroid_scores):
    distance = ((scores[:, None, :] - centroid_scores[None, :, :]) ** 2).sum(axis=2)
    return np.asarray(shared_archetypes, dtype=object)[np.argmin(distance, axis=1)]


specimen_pair_reproduction = {}
specimen_pair_correlations = {}
for mouse_name in mouse_specimens:
    mouse_pair_standardized, mouse_pair_scores = _project_to_primary(
        specimen_curves[mouse_name], archetype_columns, shared_pca, primary_whitening_spread,
        shared_n_pc)
    mouse_pair_labels = _nearest_archetype(mouse_pair_scores, archetype_centroid_scores)
    for human_name in human_specimens:
        human_pair_standardized, human_pair_scores = _project_to_primary(
            specimen_curves[human_name], archetype_columns, shared_pca, primary_whitening_spread,
            shared_n_pc)
        human_pair_labels = _nearest_archetype(human_pair_scores, archetype_centroid_scores)
        reproduced_mouse = mouse_pair_labels[archetype_eligible] == archetype_mouse_label
        reproduced_human = human_pair_labels[archetype_eligible] == archetype_human_label
        pair_name = f'{mouse_name}|{human_name}'
        specimen_pair_reproduction[pair_name] = {
            'mouse': reproduced_mouse,
            'human': reproduced_human,
            'transition': reproduced_mouse & reproduced_human,
        }
        specimen_pair_correlations[pair_name] = _row_block_pearson(
            mouse_pair_standardized[archetype_eligible], human_pair_standardized[archetype_eligible],
            ARCHETYPE_CONFIG['numerical_epsilon'])

specimen_pair_names = sorted(specimen_pair_reproduction)
_pair_mouse = np.asarray([specimen_pair_reproduction[name]['mouse'] for name in specimen_pair_names],
                         dtype=float)
_pair_human = np.asarray([specimen_pair_reproduction[name]['human'] for name in specimen_pair_names],
                         dtype=float)
_pair_transition = np.asarray(
    [specimen_pair_reproduction[name]['transition'] for name in specimen_pair_names], dtype=float)
_pair_correlation = np.asarray([specimen_pair_correlations[name] for name in specimen_pair_names],
                               dtype=float)

specimen_consistency = pd.DataFrame({
    'gene': archetype_genes,
    'mouse_archetype': archetype_mouse_label,
    'human_archetype': archetype_human_label,
    'transition': archetype_transition_labels,
    'archetype_status': archetype_assignments['archetype_status'].to_numpy(),
    'n_specimen_pairs': len(specimen_pair_names),
    'fraction_specimen_pairs_reproducing_mouse_archetype': np.nanmean(_pair_mouse, axis=0),
    'fraction_specimen_pairs_reproducing_human_archetype': np.nanmean(_pair_human, axis=0),
    'fraction_specimen_pairs_reproducing_transition': np.nanmean(_pair_transition, axis=0),
    'median_specimen_pair_shape_correlation': np.nanmedian(_pair_correlation, axis=0),
    'primary_shape_correlation': archetype_shape_correlation,
})
specimen_consistency.to_csv(
    OUTPUT_DIR / 'shared_archetype_specimen_consistency.csv', index=False)
transition_consistency = specimen_consistency.groupby('transition', observed=True).agg(
    n_genes=('gene', 'size'),
    median_fraction_reproducing_mouse_archetype=(
        'fraction_specimen_pairs_reproducing_mouse_archetype', 'median'),
    median_fraction_reproducing_human_archetype=(
        'fraction_specimen_pairs_reproducing_human_archetype', 'median'),
    median_fraction_reproducing_transition=(
        'fraction_specimen_pairs_reproducing_transition', 'median'),
    median_specimen_pair_shape_correlation=('median_specimen_pair_shape_correlation', 'median'),
    median_primary_shape_correlation=('primary_shape_correlation', 'median'),
).reset_index().sort_values('n_genes', ascending=False)
transition_consistency.to_csv(
    OUTPUT_DIR / 'shared_archetype_specimen_consistency_by_transition.csv', index=False)

print(f'Specimen pairs evaluated: {len(specimen_pair_names)} '
      f'({len(mouse_specimens)} mouse specimens x {len(human_specimens)} human sections): '
      f'{", ".join(specimen_pair_names)}')
print(f'  median fraction of pairs reproducing the primary MOUSE archetype: '
      f'{np.nanmedian(specimen_consistency["fraction_specimen_pairs_reproducing_mouse_archetype"]):.3f}')
print(f'  median fraction of pairs reproducing the primary HUMAN archetype: '
      f'{np.nanmedian(specimen_consistency["fraction_specimen_pairs_reproducing_human_archetype"]):.3f}')
print(f'  median fraction of pairs reproducing the primary TRANSITION: '
      f'{np.nanmedian(specimen_consistency["fraction_specimen_pairs_reproducing_transition"]):.3f}')
print(f'  median pair-level mouse-vs-human shape correlation: '
      f'{np.nanmedian(specimen_consistency["median_specimen_pair_shape_correlation"]):.3f} '
      f'(the specimen-balanced value is '
      f'{np.nanmedian(archetype_shape_correlation):.3f})')
display(transition_consistency.head(10).round(3))
print('A transition whose consistency is near chance is a boundary effect of the specimen-balanced '
      'curves, not a biological finding, and 12.17 downgrades it accordingly. The two human sections '
      'remain ONE donor: these pairs are technical robustness, never biological replication.')


# %%
# Purpose: figure 16 - compact threshold and specimen-pair consistency summary.
_sensitivity_evaluated = shared_archetype_threshold_sensitivity[
    shared_archetype_threshold_sensitivity['status'].eq('evaluated')]
fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))

axes[0].bar(_sensitivity_evaluated['amplitude_floor'].astype(str),
            _sensitivity_evaluated['n_eligible_genes'], color='#6B7A8F', label='eligible genes')
_twin = axes[0].twinx()
_twin.plot(_sensitivity_evaluated['amplitude_floor'].astype(str),
           _sensitivity_evaluated['n_archetypes'], color='#D55E00', marker='o', ms=5,
           label='archetypes')
_twin.set_ylabel('shared archetypes', color='#D55E00')
_twin.tick_params(axis='y', colors='#D55E00')
axes[0].set_xlabel('amplitude floor per species')
axes[0].set_ylabel('genes eligible in both species')
axes[0].set_title('Do the archetypes persist?')

axes[1].plot(_sensitivity_evaluated['amplitude_floor'].astype(str),
             _sensitivity_evaluated['median_best_centroid_spearman_vs_primary'], marker='o', ms=5,
             color='#4C9BD3', label='centroid agreement with primary')
axes[1].plot(_sensitivity_evaluated['amplitude_floor'].astype(str),
             _sensitivity_evaluated['median_fraction_primary_transitions_reproduced'], marker='s',
             ms=5, color='#E15759', label='primary transitions reproduced')
axes[1].plot(_sensitivity_evaluated['amplitude_floor'].astype(str),
             _sensitivity_evaluated['fraction_conserved'], marker='^', ms=5, color='#59A14F',
             label='fraction conserved')
axes[1].set_ylim(-0.02, 1.05)
axes[1].set_xlabel('amplitude floor per species')
axes[1].set_ylabel('agreement with the primary solution')
axes[1].set_title('Do the transitions persist?')
axes[1].legend(frameon=False, fontsize=7)

_transition_consistency_shown = transition_consistency[
    transition_consistency['n_genes'] >= ARCHETYPE_CONFIG['pathway_min_transition_genes']].copy()
if _transition_consistency_shown.empty:
    _transition_consistency_shown = transition_consistency.head(10).copy()
_transition_consistency_shown = _transition_consistency_shown.sort_values(
    'median_fraction_reproducing_transition', ascending=False).head(14)
axes[2].barh(_transition_consistency_shown['transition'][::-1],
             _transition_consistency_shown['median_fraction_reproducing_transition'][::-1],
             color='#59A14F')
axes[2].set_xlim(0, 1.05)
axes[2].set_xlabel('median fraction of specimen pairs reproducing the transition')
axes[2].set_title('Specimen-pair consistency')
axes[2].tick_params(axis='y', labelsize=7)
fig.suptitle('Figure 16 - sensitivity to the amplitude floor and to the specimen pairing', fontsize=12)
_save_figure(fig, 'fig16_archetype_sensitivity_and_specimen_consistency.png')


# %% [markdown]
# ### 12.14 - representative genes for the transitions the data actually produced

# %%
# Purpose: 12.14 - a few high-interest transitions, with the standardised curves and the original fits.
top_transitions = transition_rows[
    transition_rows['n_genes'] >= ARCHETYPE_CONFIG['pathway_min_transition_genes']].sort_values(
    'n_genes', ascending=False).head(ARCHETYPE_CONFIG['representative_transitions'])
if top_transitions.empty:
    print('No transition reaches the minimum size, so no representative-gene figure was written.')
else:
    fig, axes = plt.subplots(len(top_transitions), 2,
                             figsize=(11.5, 2.7 * len(top_transitions)), squeeze=False,
                             sharex=True)
    for row_index, transition_row in enumerate(top_transitions.itertuples()):
        label = f'{transition_row.mouse_archetype}->{transition_row.human_archetype}'
        members = ((archetype_mouse_label == transition_row.mouse_archetype)
                   & (archetype_human_label == transition_row.human_archetype))
        positions = np.flatnonzero(members)
        # Representative of a TRANSITION means it carries real amplitude in both species, not merely
        # that it sits near a centroid boundary.
        rank = np.minimum(gene_mouse_amplitude[archetype_eligible][positions],
                          gene_human_amplitude[archetype_eligible][positions])
        chosen = positions[np.argsort(-rank)][
            :ARCHETYPE_CONFIG['representative_genes_per_transition']]
        for position in chosen:
            axes[row_index][0].plot(archetype_grid, Z_mouse[position], lw=1.5,
                                    color=SPECIES_COLORS['mouse'], alpha=0.85)
            axes[row_index][0].plot(archetype_grid, Z_human[position], lw=1.5, ls='--',
                                    color=SPECIES_COLORS['human'], alpha=0.85)
            # The original curves span the whole grid while `archetype_grid` is the shared-column
            # subset, so the same column mask has to be applied before they share an x axis.
            axes[row_index][1].plot(
                archetype_grid, mouse_curves[archetype_eligible][position][archetype_columns],
                lw=1.4, color=SPECIES_COLORS['mouse'], alpha=0.85)
            axes[row_index][1].plot(
                archetype_grid, human_curves[archetype_eligible][position][archetype_columns],
                lw=1.4, ls='--', color=SPECIES_COLORS['human'], alpha=0.85)
        axes[row_index][0].annotate(
            '; '.join(archetype_genes[chosen]), xy=(0.01, 0.02), xycoords='axes fraction',
            fontsize=6, style='italic')
        axes[row_index][0].axhline(0, color='k', lw=0.5)
        axes[row_index][1].axhline(0, color='k', lw=0.5)
        axes[row_index][0].set_ylabel(f'{label}\nz-score', fontsize=8)
        axes[row_index][1].set_ylabel('fitted lognorm', fontsize=8)
        axes[row_index][0].set_title(
            f'{label}: {transition_row.mouse_archetype} mouse -> {transition_row.human_archetype} '
            f'human, n={int(transition_row.n_genes)} (standardised)', fontsize=8)
        axes[row_index][1].set_title('Original unstandardised fits for the same genes', fontsize=8)
    for axis in axes[-1]:
        axis.set_xlabel('shared PT DPT')
    axes[0][0].plot([], [], color=SPECIES_COLORS['mouse'], label='mouse')
    axes[0][0].plot([], [], color=SPECIES_COLORS['human'], ls='--', label='human')
    axes[0][0].legend(frameon=False, fontsize=7)
    fig.suptitle('Figure 17 - representative genes per transition: standardised shapes, and the '
                 'original fits behind them', fontsize=12)
    _save_figure(fig, 'fig17_representative_transition_genes.png')


# %% [markdown]
# ### 12.17 - what section 12 supports
#
# Written from the objects this notebook computed. Nothing here restates a number it did not print.

# %%
# Purpose: 12.17 - the new section's conclusion, derived from the executed outputs.
print('Section 12 - shared positional archetypes and cross-species program switching')
print()

print(f'1. Shared positional archetypes. {len(shared_archetypes)} archetypes '
      f'({", ".join(shared_archetypes)}) were learned jointly from '
      f'{shared_standardized.shape[0]:,} within-species standardised curves '
      f'({archetype_genes.size:,} genes x 2 species) on {shared_standardized.shape[1]} shared grid '
      f'points. {int(_represented_both.sum())} of them are occupied by both species. Median '
      f'agreement between the mouse-only and human-only centroid of an archetype: '
      f'{archetype_catalog["centroid_spearman_mouse_vs_human"].median():.3f}.')

print(f'2. Same archetype in both species. {n_conserved:,} of {len(archetype_assignments):,} genes '
      f'({n_conserved / len(archetype_assignments):.1%}) keep the same archetype; {n_switched:,} '
      f'({n_switched / len(archetype_assignments):.1%}) switch. The median direct mouse-vs-human '
      f'shape correlation is {np.nanmedian(archetype_shape_correlation):.3f} overall, '
      f'{archetype_catalog["median_shape_correlation_conserved"].median():.3f} for genes that stay '
      f'and {archetype_catalog["median_shape_correlation_switched_out"].median():.3f} for genes that '
      'leave their mouse archetype.')

_top_transition_rows = transition_rows.sort_values('n_genes', ascending=False).head(3)
print('3. Most common mouse -> human transitions. ' + '; '.join(
    f'{row.mouse_archetype}->{row.human_archetype} '
    f'({int(row.n_genes)} genes, {row.row_fraction:.1%} of the row, O/E '
    f'{row.observed_over_expected:.2f})' for row in _top_transition_rows.itertuples()) + '.')
print('   By enrichment rather than by size: ' + '; '.join(
    f'{row.mouse_archetype}->{row.human_archetype} (O/E {row.observed_over_expected:.2f}, '
    f'{int(row.n_genes)} genes)' for row in _off_diagonal.sort_values(
        'observed_over_expected', ascending=False).head(3).itertuples()) + '.')

if response_module_archetype_summary.empty:
    print('4. Response modules: no response module carried an eligible gene.')
else:
    print(f'4. Response modules as positional-program transitions. '
          f'{len(_dominant_transition_modules)} of {len(response_module_archetype_summary)} response '
          f'modules have >= 50% of their members in ONE transition, and the median module puts '
          f'{response_module_archetype_summary["fraction_most_common_transition"].median():.1%} of '
          f'its genes in its most common transition. The most concentrated: ' + '; '.join(
              f'{row.response_module} -> {row.most_common_transition} '
              f'({row.fraction_most_common_transition:.1%})'
              for row in response_module_archetype_summary.sort_values(
                  'fraction_most_common_transition', ascending=False).head(3).itertuples()) + '.')

if transition_enrichment_significant.empty:
    print('5. Pathways by conserved versus switching programs: no transition x pathway pair reached '
          'BH q < 0.05, so this section identifies no pathway-specific program transition.')
else:
    _conserved_only = transition_enrichment_significant[
        transition_enrichment_significant['set_type'].eq('conserved within archetype')]
    _switching_only = transition_enrichment_significant[
        ~transition_enrichment_significant['set_type'].eq('conserved within archetype')]
    print(f'5. Pathways by conserved versus switching programs. At BH q < 0.05 there are '
          f'{len(_conserved_only)} conserved-set and {len(_switching_only)} switching/transition '
          f'results. Strongest conserved: ' + ('; '.join(
              f'{row.pathway} in {row.transition_set} (q={row.p_value_adjusted:.1e}, '
              f'O/E={row.observed_over_expected:.1f})' for row in _conserved_only.sort_values(
                  'p_value_adjusted').head(3).itertuples()) or 'none') + '. Strongest switching: ' +
          ('; '.join(f'{row.pathway} in {row.transition_set} (q={row.p_value_adjusted:.1e}, '
                     f'O/E={row.observed_over_expected:.1f})' for row in _switching_only.sort_values(
                         'p_value_adjusted').head(3).itertuples()) or 'none') + '.')

if pathway_archetype_coherence.empty:
    print('6. Pathway heterogeneity: no pathway had enough eligible genes to assess.')
else:
    _heterogeneous = pathway_archetype_coherence.sort_values(
        'z_transition_entropy', ascending=False).head(3)
    print(f'6. Pathways whose genes are internally heterogeneous. Highest transition entropy against '
          f'the size- and amplitude-matched null: ' + '; '.join(
              f'{row.pathway.split(": ")[-1]} (z={row.z_transition_entropy:.2f}, '
              f'fraction conserved {row.fraction_conserved:.2f}, n={row.eligible_gene_count})'
              for row in _heterogeneous.itertuples()) +
          '. These pathways span several positional programs at once, so an aggregate pathway curve '
          'would average opposing shapes and can look modest while the members are not.')

_strictest_evaluated = _sensitivity_evaluated.sort_values('amplitude_floor').tail(1)
print(f'7. Robustness. Median specimen-pair reproduction of a primary transition is '
      f'{np.nanmedian(specimen_consistency["fraction_specimen_pairs_reproducing_transition"]):.3f} '
      f'(mouse archetype '
      f'{np.nanmedian(specimen_consistency["fraction_specimen_pairs_reproducing_mouse_archetype"]):.3f}, '
      f'human archetype '
      f'{np.nanmedian(specimen_consistency["fraction_specimen_pairs_reproducing_human_archetype"]):.3f}). '
      + (f'At the strictest evaluable floor ({_strictest_evaluated["amplitude_floor"].iloc[0]}) '
         f'{int(_strictest_evaluated["n_archetypes"].iloc[0])} archetypes remain over '
         f'{int(_strictest_evaluated["n_eligible_genes"].iloc[0]):,} genes, centroid agreement with '
         f'the primary solution is '
         f'{_strictest_evaluated["median_best_centroid_spearman_vs_primary"].iloc[0]:.3f} and the '
         f'fraction conserved is '
         f'{_strictest_evaluated["fraction_conserved"].iloc[0]:.1%}.'
         if len(_sensitivity_evaluated) and np.isfinite(
             _strictest_evaluated["median_best_centroid_spearman_vs_primary"].iloc[0])
         else 'No stricter floor could be evaluated.'))

print()
print('Reading rules specific to section 12:')
print('  - the analysis is within-species POSITIONAL SHAPE. It never compares a human level with a')
print('    mouse level and never subtracts one species from the other, so a level shift or an')
print('    amplitude change cannot masquerade as a shape change here.')
print('  - the archetypes were learned jointly from both species and were frozen before any pathway')
print('    was read; the number of archetypes was never chosen on enrichment.')
print('  - "conserved" and "switched" are descriptive labels on two cluster assignments, not a test of')
print('    a species effect. The human side is two sections of ONE donor, so the inference unit is the')
print('    specimen and the comparison is 2-vs-1-donor.')
print('  - a shared archetype with poorly conserved membership is a distinct finding, not a failed')
print('    analysis: it separates "the shape is available in both species" from "the same genes use it".')
print('  - the two human sections are never described as biological replicates, in this section or')
print('    anywhere else.')

pd.DataFrame([{'figure': name} for name in figure_index]).to_csv(
    DIAGNOSTIC_DIR / 'figure_index.csv', index=False)
_new_files = sorted(name for name in (
    'shared_archetype_gene_metrics.csv', 'shared_archetype_pca_scores.csv',
    'shared_archetype_assignments.csv', 'shared_archetype_catalog.csv',
    'mouse_human_archetype_transitions.csv', 'archetype_transition_enrichment.csv',
    'response_module_archetype_summary.csv', 'pathway_archetype_coherence.csv',
    'shared_archetype_threshold_sensitivity.csv', 'shared_archetype_specimen_consistency.csv',
    'shared_archetype_specimen_consistency_by_transition.csv'))
print()
print(f'Section 12 wrote {len(_new_files)} tables into {_rel(OUTPUT_DIR)} and '
      f'{len([name for name in figure_index if name.startswith(("fig1", "fig12"))])} figures; '
      'figure_index.csv was rewritten so the index covers every figure the notebook produced.')


# %%
# Purpose: close with what was written, and what must not be read into it.
written = sorted(path.name for path in OUTPUT_DIR.iterdir() if path.is_file())
print(f'Tables and figures written to {_rel(OUTPUT_DIR)}:')
for name in written:
    print('  ', name)
print()
print(f'PT structures: {adata_pt.n_obs:,}; eligible genes: {gene_names.size:,}; '
      f'positional modules: {positional["k"]} (k diagnostics in {_rel(DIAGNOSTIC_DIR)}); '
      f'response modules: {response["k"]}.')
print()
print('Reading rules for these numbers:')
print('  - the coordinate and the PT subset are upstream (notebook 03). Nothing here re-derives')
print('    Harmony, the diffusion map, DPT, the PT trajectory or the ortholog map.')
print('  - genes absent from either input feature list never enter a fit or an enrichment background')
print('    (measured_in_both_inputs): a structural zero is not a measurement.')
print('  - the primary cross-species analysis clusters the CENTRED human-minus-mouse response. A large')
print('    global level shift is reported in the metrics table, not turned into a response module.')
print('  - the inference unit is the specimen: two mouse specimens versus two sections of ONE human')
print('    donor, so these are descriptive summaries, not species-level tests.')
print('  - pathway annotations were loaded only after the clustering was final, and never revise it.')
print('  - the PCA axes are reported as biological output (loading curves and per-PC gene extremes),')
print('    not merely as preprocessing for the clustering.')
