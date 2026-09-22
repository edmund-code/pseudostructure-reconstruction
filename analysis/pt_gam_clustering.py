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
