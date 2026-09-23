# %% [markdown]
# # Genes2Genes alignment of healthy mouse and human PT trajectories
#
# **Question.** Which one-to-one ortholog trajectories are conserved or divergent, where do they diverge, and what biology characterizes their alignment patterns?
#
# Mouse is always the **reference** and human the **query**. Inputs are PT structures from notebook 03; both human sections are healthy cortex from **one donor**, including the section named `HUK1_MED1`. All cross-species and pathway results are descriptive. [Genes2Genes v0.2.0](https://github.com/Teichlab/Genes2Genes) performs distributional interpolation and five-state gene alignment on individual structures. See the [method paper](https://doi.org/10.1038/s41592-024-02378-4). Run notebook 03 first; notebook 07 is needed only for the optional mouse P-module bridge.
#
# The notebook is deliberately output-free in Git. Run it top to bottom in the analysis environment with Genes2Genes installed. The full gene universe and four pairwise runs can take hours; completed gene chunks are cached in the results directory.

# %% [markdown]
# ## 1. Inputs and gene universe

# %%
import argparse
import hashlib
import json
from importlib.metadata import version as package_version
import os
import sys
from pathlib import Path

CACHE_ROOT = Path(os.environ.get('PSEUDOSPACE_CACHE_ROOT', '/tmp/pseudospace_pt_genes2genes_cache'))
for subdir in ('matplotlib', 'numba', 'xdg'):
    (CACHE_ROOT / subdir).mkdir(parents=True, exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR', str(CACHE_ROOT / 'matplotlib'))
os.environ.setdefault('NUMBA_CACHE_DIR', str(CACHE_ROOT / 'numba'))
os.environ.setdefault('XDG_CACHE_HOME', str(CACHE_ROOT / 'xdg'))

def _project_dir():
    starts = [Path.cwd().resolve()]
    if '__file__' in globals():
        starts.insert(0, Path(__file__).resolve().parent)
    for start in starts:
        for candidate in (start, *start.parents):
            if (candidate / 'pseudospace').is_dir() and (candidate / 'data').is_dir():
                return candidate
    raise RuntimeError('Cannot locate the pseudospace project root.')

PROJECT_DIR = _project_dir()
if str(PROJECT_DIR) not in sys.path:
    sys.path.insert(0, str(PROJECT_DIR))
parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--data-root', type=Path)
parser.add_argument('--results-root', type=Path)
args, _ = parser.parse_known_args()
DATA_ROOT = Path(args.data_root or os.environ.get('PSEUDOSPACE_DATA_ROOT') or PROJECT_DIR / 'data').expanduser().resolve()
RESULTS_ROOT = Path(args.results_root or os.environ.get('PSEUDOSPACE_RESULTS_ROOT') or PROJECT_DIR / 'results').expanduser().resolve()
UPSTREAM = RESULTS_ROOT / 'human_vs_healthy_mouse'
OUTPUT_DIR = RESULTS_ROOT / 'pt_genes2genes'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
INPUT_PATH = UPSTREAM / 'cross_species_pt_dpt.h5ad'
MAP_PATH = UPSTREAM / 'ortholog_map_used.csv'
ATLAS_PATH = RESULTS_ROOT / 'pt_gam_clustering' / 'mouse_positional_module_assignments.csv'
for required in (INPUT_PATH, MAP_PATH):
    if not required.exists():
        raise FileNotFoundError(f'Missing {required}; run notebook 03 first.')
NOTEBOOK_LOGIC_VERSION = '09-g2g-v1'
HUMAN_MARKER_TRIAL_VERSION = '09-human-marker-trial-v1'
G2G_BINS = 14  # explicit API argument; 14 is the tutorial's demonstrated interpolation resolution
G2G_CHUNK_GENES = 128  # memory ceiling: G2G densifies each input matrix; raise only after profiling RAM
RNG = 0
print(f'Input: {INPUT_PATH}; output: {OUTPUT_DIR}')

# %%
import numpy as np
import warnings
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import entropy, hypergeom, mannwhitneyu, spearmanr
from sklearn.metrics import silhouette_score
from statsmodels.stats.multitest import multipletests
from IPython.display import display
from pseudospace.pathways import build_pathway_membership

try:
    from genes2genes.Main import RefQueryAligner
    import genes2genes
    import Levenshtein
except ImportError as exc:
    raise ImportError('Install the pinned Genes2Genes analysis dependency from environment.yml.') from exc

warnings.filterwarnings('ignore', message='divide by zero encountered in log',
                        category=RuntimeWarning, module=r'genes2genes\.OrgAlign')
assert genes2genes.__version__ == '0.2.0', 'This notebook was checked against Genes2Genes 0.2.0.'
SPECIMENS = {'mouse': ('Ctrl1A2', 'Ctrl1A4'), 'human': ('HUK1_COR1', 'HUK1_MED1')}
COLORS = {'mouse': '#0072B2', 'human': '#D55E00'}
MARKERS = {'S1': ('Slc5a2', 'Slc5a12', 'Gatm', 'Lrp2', 'Cubn', 'Slc34a1'),
           'S3': ('Slc22a7', 'Cyp7b1', 'Slc7a13', 'Slc6a18', 'Acsm3')}
HUMAN_TRIAL_MARKERS = {'S1': ('SLC5A2', 'SLC5A12'),
                       'S2': ('SLC22A6', 'SLC13A3', 'ACSM3'),
                       'S3': ('SLC22A7', 'SLC7A13', 'AGXT', 'DCXR')}

def savefig(fig, name):
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / name, dpi=180, bbox_inches='tight')
    plt.show()
    plt.close(fig)

def dense_column(matrix, index):
    value = matrix[:, index]
    return np.asarray(value.toarray() if sparse.issparse(value) else value).ravel()


# %%
adata = sc.read_h5ad(INPUT_PATH)
required_obs = {'comparison_species', 'sample', 'broad_tubule_marker_call'}
if not required_obs.issubset(adata.obs) or 'measured_in_both_inputs' not in adata.var or 'lognorm' not in adata.layers:
    raise ValueError('Notebook 03 PT artifact lacks species/sample/PT labels, availability, or lognorm.')
mask = adata.obs['broad_tubule_marker_call'].astype(str).eq('PT')
adata = adata[mask].copy()
actual = {species: tuple(sorted(adata.obs.loc[adata.obs.comparison_species.astype(str).eq(species), 'sample'].astype(str).unique())) for species in SPECIMENS}
if actual != SPECIMENS:
    raise ValueError(f'Expected the specified two mouse specimens and human sections; found {actual}.')
if not adata.var_names.is_unique or not adata.obs_names.is_unique:
    raise ValueError('G2G needs unique gene and structure names.')
orthologs = pd.read_csv(MAP_PATH)
if not {'mouse_symbol', 'human_symbol'}.issubset(orthologs):
    raise ValueError('Accepted ortholog map has unexpected columns.')
measured = adata.var['measured_in_both_inputs'].to_numpy(dtype=bool)
lognorm = adata.layers['lognorm'].tocsr() if sparse.issparse(adata.layers['lognorm']) else sparse.csr_matrix(adata.layers['lognorm'])
if lognorm.data.size and (not np.isfinite(lognorm.data).all() or lognorm.data.min() < 0):
    raise ValueError('lognorm must be finite and nonnegative.')
detection = np.asarray((lognorm > 0).sum(axis=0)).ravel() / adata.n_obs
eligible = measured & (detection >= 0.02)
universe = pd.DataFrame({'gene': adata.var_names, 'assayed': True, 'measured_in_both': measured,
                         'pt_detected_fraction': detection, 'passes_detection': detection >= 0.02,
                         'g2g_eligible': eligible})
universe.to_csv(OUTPUT_DIR / 'g2g_gene_universe.csv', index=False)
genes = adata.var_names[eligible].tolist()
if len(genes) < 20:
    raise ValueError('Too few G2G genes after the measured-in-both and 2% detection rules.')
print({'accepted_orthologs': len(orthologs),
       'assayed_orthologs': int(orthologs.mouse_symbol.isin(adata.var_names).sum()),
       'measured_in_both': int(measured.sum()),
       'pass_detection': int((detection >= 0.02).sum()), 'final_g2g_genes': len(genes)})
display(adata.obs.groupby(['comparison_species', 'sample'], observed=True).size().rename('pt_structures').reset_index())

# %% [markdown]
# ## 2. Species-specific PT pseudotimes and minimal validation
#
# Notebook 02 uses a different mouse structure set, so a direct row-for-row transfer is unsafe. Compute PT-only DPT within each species from lognorm expression, without cross-species Harmony. The species-specific PT coordinates are cached by input artifact and parameters so repeated notebook runs reuse exactly the same times. The S1/S3 marker axis selects an early root and orients DPT. The gate below checks each of the four sections separately.
#
# The existing `orient_and_normalize` helper uses the oriented DPT p5 and p95 as the 0 and 1 anchors and clips the tails. This robust within-species scaling is recorded because endpoint occupancy affects G2G interpolation.

# %%
from pseudospace.trajectory import orient_and_normalize
from pseudospace.stage_cache import cached_payload

def marker_score(a, labels):
    present = [g for g in labels if g in a.var_names]
    if len(present) < 2:
        raise ValueError(f'Only {len(present)} marker genes are assayed: {labels}')
    x = a[:, present].X
    x = x.toarray() if sparse.issparse(x) else np.asarray(x)
    sd = x.std(axis=0)
    use = sd > 0
    if use.sum() < 2:
        raise ValueError(f'Too few variable markers in {present}')
    return ((x[:, use] - x[:, use].mean(axis=0)) / sd[use]).mean(axis=1)

def species_pt(species):
    a = adata[adata.obs.comparison_species.astype(str).eq(species), genes].copy()
    a.X = a.layers['lognorm'].copy()
    a.obs['S1_score'] = marker_score(a, MARKERS['S1'])
    a.obs['S3_score'] = marker_score(a, MARKERS['S3'])
    axis = (a.obs.S3_score - a.obs.S1_score).to_numpy()
    a.obs['marker_axis'] = axis
    def fit_time():
        # PCA is fitted within species; no joint embedding or cross-species coordinate enters DPT.
        sc.pp.pca(a, n_comps=30, random_state=RNG)
        sc.pp.neighbors(a, n_neighbors=30, use_rep='X_pca', random_state=RNG)
        sc.tl.diffmap(a)
        candidates = np.flatnonzero(axis <= np.quantile(axis, 0.05))
        position = a.obsm['X_pca'][candidates]
        root = candidates[np.linalg.norm(position - np.median(position, axis=0), axis=1).argmin()]
        a.uns['iroot'] = int(root)
        sc.tl.dpt(a)
        return {'time': orient_and_normalize(a.obs.dpt_pseudotime.to_numpy(), axis, min_valid=10),
                'root': np.asarray(int(root))}

    artifact = INPUT_PATH.stat()
    payload = cached_payload(
        f'pt_time_{species}', fit_time, root=OUTPUT_DIR / 'stage_cache',
        params={'logic': NOTEBOOK_LOGIC_VERSION, 'n_pcs': 30, 'n_neighbors': 30,
                'random_state': RNG, 'markers': MARKERS, 'scanpy': package_version('scanpy')},
        inputs={'input_size': artifact.st_size, 'input_mtime_ns': artifact.st_mtime_ns,
                'obs_names': a.obs_names.to_numpy(), 'genes': genes}, code=NOTEBOOK_LOGIC_VERSION)
    time = np.asarray(payload['time'], dtype=np.float64)
    root = int(payload['root'])
    if len(time) != a.n_obs or not 0 <= root < a.n_obs or not np.isfinite(time).all() or (time < 0).any() or (time > 1).any() or len(np.unique(time)) < 10:
        raise ValueError(f'{species}: invalid or collapsed DPT pseudotime cache.')
    a.obs['time'] = time
    a.uns['iroot'] = root
    print(f'{species}: root {a.obs_names[root]}, DPT-marker Spearman {spearmanr(a.obs.time, axis).statistic:.3f}')
    return a

adata_mouse = species_pt('mouse')
adata_human = species_pt('human')

# %%
qc = []
fig, axes = plt.subplots(2, 3, figsize=(12, 6), sharex='col')
for row, (species, a) in enumerate((('mouse', adata_mouse), ('human', adata_human))):
    for sample, group in a.obs.groupby('sample', observed=True):
        group = group.sort_values('time')
        axes[row, 0].hist(group.time, bins=20, alpha=.45, label=str(sample))
        for col, marker in ((1, 'S1_score'), (2, 'S3_score')):
            axes[row, col].scatter(group.time, group[marker], s=2, alpha=.06, label=str(sample))
        values = {'species': species, 'sample': str(sample), 'n_structures': len(group),
                  'time_min': group.time.min(), 'time_median': group.time.median(), 'time_max': group.time.max(),
                  'rho_marker_axis': spearmanr(group.time, group.marker_axis).statistic,
                  'rho_S1': spearmanr(group.time, group.S1_score).statistic,
                  'rho_S3': spearmanr(group.time, group.S3_score).statistic}
        qc.append(values)
    for col, title in enumerate(('PT time occupancy', 'S1 score', 'S3 score')):
        axes[row, col].set_title(f'{species}: {title}')
        axes[row, col].set_xlabel('species-specific PT time')
    axes[row, 0].legend(fontsize=7)
qc = pd.DataFrame(qc)
qc.to_csv(OUTPUT_DIR / 'g2g_pseudotime_orientation.csv', index=False)
savefig(fig, 'fig01_species_pt_orientation.png')
display(qc.round(3))
if ((qc.rho_marker_axis <= 0) | (qc.rho_S1 >= 0) | (qc.rho_S3 <= 0)).any():
    raise RuntimeError('S1 -> S3 orientation failed in at least one specimen/section; G2G is not biologically interpretable.')

# %% [markdown]
# ### Human-only marker trial (separate DPT diagnostic)
#
# The proposed human S1/S2/S3 symbols are mapped through the accepted ortholog table to the shared gene names. S1 and S3 choose the trial root and orientation exactly as in the existing DPT; S2 is checked for an intermediate peak. A missing marker is recorded rather than silently substituted. This trial writes new QC files and a figure, while the original human coordinate and all downstream G2G alignments remain unchanged.

# %%
human_trial_coverage = []
human_trial_genes = {}
eligible_genes = set(genes)
for stage, symbols in HUMAN_TRIAL_MARKERS.items():
    present = []
    for symbol in symbols:
        matches = orthologs.loc[orthologs.human_symbol.eq(symbol), 'mouse_symbol']
        if len(matches) != 1:
            raise ValueError(f'{symbol}: expected one accepted human-to-mouse ortholog, found {len(matches)}')
        shared_gene = str(matches.iloc[0])
        available = shared_gene in eligible_genes
        human_trial_coverage.append({'stage': stage, 'human_symbol': symbol,
                                     'shared_gene': shared_gene, 'in_pt_gene_set': available})
        if available:
            present.append(shared_gene)
    if len(present) < 2:
        raise ValueError(f'Human {stage} trial has fewer than two available markers: {present}')
    human_trial_genes[stage] = present
human_trial_coverage = pd.DataFrame(human_trial_coverage)
human_trial_coverage.to_csv(OUTPUT_DIR / 'g2g_human_marker_trial_coverage.csv', index=False)
display(human_trial_coverage)

human_trial = adata_human.copy()
if not human_trial.obs_names.equals(adata_human.obs_names):
    raise ValueError('Human trial and baseline PT structures differ or are out of order.')
for stage, names in human_trial_genes.items():
    human_trial.obs[f'{stage}_trial_score'] = marker_score(human_trial, names)
trial_axis = (human_trial.obs.S3_trial_score - human_trial.obs.S1_trial_score).to_numpy()
human_trial.obs['trial_marker_axis'] = trial_axis

def fit_human_marker_trial():
    sc.pp.pca(human_trial, n_comps=30, random_state=RNG)
    sc.pp.neighbors(human_trial, n_neighbors=30, use_rep='X_pca', random_state=RNG)
    sc.tl.diffmap(human_trial)
    candidates = np.flatnonzero(trial_axis <= np.quantile(trial_axis, 0.05))
    position = human_trial.obsm['X_pca'][candidates]
    root = candidates[np.linalg.norm(position - np.median(position, axis=0), axis=1).argmin()]
    human_trial.uns['iroot'] = int(root)
    sc.tl.dpt(human_trial)
    return {'time': orient_and_normalize(human_trial.obs.dpt_pseudotime.to_numpy(), trial_axis, min_valid=10),
            'root': np.asarray(int(root))}

artifact = INPUT_PATH.stat()
trial_payload = cached_payload(
    'pt_time_human_marker_trial', fit_human_marker_trial, root=OUTPUT_DIR / 'stage_cache',
    params={'logic': HUMAN_MARKER_TRIAL_VERSION, 'n_pcs': 30, 'n_neighbors': 30,
            'random_state': RNG, 'markers': human_trial_genes, 'scanpy': package_version('scanpy')},
    inputs={'input_size': artifact.st_size, 'input_mtime_ns': artifact.st_mtime_ns,
            'obs_names': human_trial.obs_names.to_numpy(), 'genes': genes}, code=HUMAN_MARKER_TRIAL_VERSION)
trial_time = np.asarray(trial_payload['time'], dtype=np.float64)
trial_root = int(trial_payload['root'])
if (len(trial_time) != human_trial.n_obs or not 0 <= trial_root < human_trial.n_obs
        or not np.isfinite(trial_time).all() or (trial_time < 0).any() or (trial_time > 1).any()
        or len(np.unique(trial_time)) < 10):
    raise ValueError('Human marker trial produced invalid or collapsed DPT pseudotime.')
human_trial.obs['time'] = trial_time
human_trial.uns['iroot'] = trial_root

trial_qc = []
for sample, group in human_trial.obs.groupby('sample', observed=True):
    old = adata_human.obs.loc[group.index]
    bins = pd.cut(group.time, bins=np.linspace(0, 1, 11), labels=False, include_lowest=True)
    s2_by_bin = group.groupby(bins, observed=True).S2_trial_score.mean()
    trial_qc.append({'sample': str(sample), 'n_structures': len(group),
                     'baseline_root': str(adata_human.obs_names[adata_human.uns['iroot']]),
                     'trial_root': str(human_trial.obs_names[trial_root]),
                     'rho_baseline_axis_baseline_time': spearmanr(old.time, old.marker_axis).statistic,
                     'rho_trial_axis_baseline_time': spearmanr(old.time, group.trial_marker_axis).statistic,
                     'rho_trial_axis_trial_time': spearmanr(group.time, group.trial_marker_axis).statistic,
                     'rho_trial_S1': spearmanr(group.time, group.S1_trial_score).statistic,
                     'rho_trial_S2': spearmanr(group.time, group.S2_trial_score).statistic,
                     'rho_trial_S3': spearmanr(group.time, group.S3_trial_score).statistic,
                     'S2_peak_decile': int(s2_by_bin.idxmax()) + 1,
                     'rho_baseline_vs_trial_time': spearmanr(old.time, group.time).statistic})
trial_qc = pd.DataFrame(trial_qc)
trial_qc.to_csv(OUTPUT_DIR / 'g2g_human_marker_trial_qc.csv', index=False)
display(trial_qc.round(3))

fig, axes = plt.subplots(1, 2, figsize=(11, 4))
for sample, group in human_trial.obs.groupby('sample', observed=True):
    axes[0].scatter(adata_human.obs.loc[group.index, 'time'], group.time, s=2, alpha=.1, label=str(sample))
    for stage, color in (('S1', '#0072B2'), ('S2', '#E69F00'), ('S3', '#D55E00')):
        bins = pd.cut(group.time, bins=np.linspace(0, 1, 11), labels=False, include_lowest=True)
        means = group.groupby(bins, observed=True)[f'{stage}_trial_score'].mean()
        axes[1].plot((means.index.to_numpy(dtype=float) + .5) / 10, means.to_numpy(),
                     color=color, linestyle='-' if str(sample) == SPECIMENS['human'][0] else '--',
                     label=f'{sample} {stage}')
axes[0].plot([0, 1], [0, 1], color='black', linestyle='--', linewidth=.8)
axes[0].set(xlabel='existing human PT time', ylabel='trial human PT time', title='Coordinate agreement')
axes[1].set(xlabel='trial human PT time', ylabel='mean z-scored marker signal', title='S1 / S2 / S3 progression')
for ax in axes:
    ax.legend(fontsize=7)
savefig(fig, 'fig01b_human_marker_trial.png')
print('Human marker trial is diagnostic only; downstream G2G still uses the existing human PT time.')

# %% [markdown]
# ## 3. Prepare G2G objects

# %%
adata_mouse_g2g = adata_mouse[:, genes].copy()
adata_human_g2g = adata_human[:, genes].copy()
for a in (adata_mouse_g2g, adata_human_g2g):
    a.X = a.layers['lognorm'].tocsr().astype(np.float32) if sparse.issparse(a.layers['lognorm']) else sparse.csr_matrix(a.layers['lognorm'], dtype=np.float32)
    assert a.var_names.tolist() == genes and a.obs.time.between(0, 1).all()
    assert np.isfinite(a.X.data).all() and (a.X.data >= 0).all()
print(f'Mouse reference {adata_mouse_g2g.shape}; human query {adata_human_g2g.shape}; same ordered genes: {len(genes)}')


# %% [markdown]
# ## 4. Gene-level G2G alignment
#
# G2G v0.2.0 has no default interpolation-point count; this notebook passes 14, as in its tutorial. Kernel mode, state parameters, and sequential alignment stay at package defaults. Chunk size affects memory and resume behavior only; every gene is interpolated and aligned from its individual PT structures. Cached chunks are keyed by the input artifact, ordered genes, pseudotimes, version, and sample pair. `G2G` reports a log2 ratio of **mean lognorm values**, reference/query, so a positive value means higher mean lognorm in mouse; it is neither an absolute abundance fold change nor a species inference statistic.
#
# `human_minus_mouse_mean_lognorm` is the difference between pooled structure means within each species; its scale is normalized expression and it is secondary to alignment.

# %%
def align_g2g(ref, query, label):
    if ref.var_names.tolist() != query.var_names.tolist():
        raise ValueError('Reference and query genes differ or are out of order.')
    ordered = ref.var_names.tolist()
    signature = hashlib.sha256()
    signature.update(json.dumps({'logic': NOTEBOOK_LOGIC_VERSION, 'label': label,
                                 'input_size': INPUT_PATH.stat().st_size,
                                 'input_mtime_ns': INPUT_PATH.stat().st_mtime_ns,
                                 'g2g_version': genes2genes.__version__, 'bins': G2G_BINS,
                                 'genes': ordered}, sort_keys=True).encode())
    signature.update(ref.obs.time.to_numpy(dtype=np.float64).tobytes())
    signature.update(query.obs.time.to_numpy(dtype=np.float64).tobytes())
    cache = OUTPUT_DIR / 'stage_cache' / signature.hexdigest()[:20]
    cache.mkdir(parents=True, exist_ok=True)
    chunks = []
    for start in range(0, len(ordered), G2G_CHUNK_GENES):
        names = ordered[start:start + G2G_CHUNK_GENES]
        target = cache / f'{start:05d}.csv'
        if target.exists():
            part = pd.read_csv(target)
            if part.gene.tolist() != names:
                raise ValueError(f'Corrupt or stale G2G chunk: {target}')
            print(f'[stage cache] hit {label} {start}:{start + len(names)}')
        else:
            # G2G 0.2.0 densifies both AnnData inputs internally; chunking bounds peak RAM.
            fit = RefQueryAligner(ref[:, names].copy(), query[:, names].copy(), names, G2G_BINS)
            fit.align_all_pairs()  # G2G sequential default
            rows = []
            for gene in names:
                result = fit.results_map[gene]
                mouse_mean = float(np.asarray(fit.ref_mat[gene]).mean())
                human_mean = float(np.asarray(fit.query_mat[gene]).mean())
                rows.append({'gene': gene, 'alignment_string': result.alignment_str,
                             'alignment_similarity': result.match_percentage / 100,
                             'opt_alignment_cost': result.fwd_DP.opt_cost,
                             'g2g_log2fc_mouse_over_human': np.log2(mouse_mean / human_mean) if human_mean > 0 and mouse_mean > 0 else np.nan,
                             'human_minus_mouse_mean_lognorm': human_mean - mouse_mean})
            part = pd.DataFrame(rows)
            part.to_csv(target, index=False)
            del fit
            print(f'[stage cache] wrote {label} {start}:{start + len(names)}')
        chunks.append(part)
    frame = pd.concat(chunks, ignore_index=True)
    if frame.gene.tolist() != ordered or frame.alignment_string.isna().any():
        raise ValueError(f'Incomplete G2G results for {label}.')
    return frame

primary = align_g2g(adata_mouse_g2g, adata_human_g2g, 'combined')
# Reproduce get_stat_df's reference/query formula after loading the cached alignments;
# retain +/-inf when one species has zero mean instead of losing direction.
with np.errstate(divide='ignore', invalid='ignore'):
    primary['g2g_log2fc_mouse_over_human'] = np.log2(
        np.asarray(adata_mouse_g2g.X.mean(axis=0)).ravel() /
        np.asarray(adata_human_g2g.X.mean(axis=0)).ravel())
primary = primary.merge(universe[['gene', 'pt_detected_fraction']], on='gene', validate='one_to_one')
primary.to_csv(OUTPUT_DIR / 'g2g_gene_alignment_statistics.csv', index=False)
primary[['gene', 'alignment_string']].to_csv(OUTPUT_DIR / 'g2g_alignment_strings.csv', index=False)
print(f'G2G aligned {len(primary):,} genes; mouse=reference, human=query.')

# %% [markdown]
# ## 5. Alignment similarity and level-vs-dynamics comparison
#
# Alignment similarity counts `M`, `V`, and `W` as matching states. The 50% and 80% marks are **descriptive** display bands, not significance calls. The horizontal level band below uses the observed median absolute level difference as a display guide; it does not define a new eligibility filter.

# %%
sim = primary.alignment_similarity.to_numpy(dtype=float)
print({'median_similarity': np.median(sim), 'q1': np.quantile(sim, .25), 'q3': np.quantile(sim, .75),
       'fraction_at_least_80pct': np.mean(sim >= .8), 'fraction_below_50pct': np.mean(sim < .5)})
ranked_columns = ['gene', 'alignment_similarity', 'alignment_string', 'human_minus_mouse_mean_lognorm']
display(primary.nsmallest(12, 'alignment_similarity')[ranked_columns])
display(primary.nlargest(12, 'alignment_similarity')[ranked_columns])
fig, ax = plt.subplots(figsize=(7, 4))
ax.hist(sim, bins=np.linspace(0, 1, 41), color='#527A88')
for cut in (.5, .8): ax.axvline(cut, color='black', linestyle='--', alpha=.55)
ax.set(xlabel='G2G alignment similarity (M+V+W fraction)', ylabel='genes', title='PT ortholog trajectory alignment')
savefig(fig, 'fig02_alignment_similarity_distribution.png')
level = primary.human_minus_mouse_mean_lognorm.to_numpy(dtype=float)
level_guide = np.median(np.abs(level))
fig, ax = plt.subplots(figsize=(7, 5))
ax.scatter(level, sim, s=7, alpha=.22, color='#315F78')
ax.axhline(.5, color='black', linestyle='--', alpha=.5)
ax.axvline(-level_guide, color='black', linestyle=':', alpha=.5)
ax.axvline(level_guide, color='black', linestyle=':', alpha=.5)
quadrants = {'similar level + conserved': (np.abs(level) <= level_guide) & (sim >= .8),
             'different level + conserved': (np.abs(level) > level_guide) & (sim >= .8),
             'similar level + divergent': (np.abs(level) <= level_guide) & (sim < .5),
             'different level + divergent': (np.abs(level) > level_guide) & (sim < .5)}
examples = []
for label, eligible_rows in quadrants.items():
    candidates = primary.loc[eligible_rows]
    if candidates.empty:
        continue
    row = candidates.nlargest(1, 'pt_detected_fraction').iloc[0]
    examples.append({'category': label, 'gene': row.gene, 'alignment_similarity': row.alignment_similarity,
                     'human_minus_mouse_mean_lognorm': row.human_minus_mouse_mean_lognorm})
    ax.scatter(row.human_minus_mouse_mean_lognorm, row.alignment_similarity, s=55,
               color=plt.get_cmap('tab10')(len(examples) - 1), edgecolor='black',
               label=f'{row.gene}: {label}')
ax.legend(loc='lower left', fontsize=7, frameon=True)
pd.DataFrame(examples).to_csv(OUTPUT_DIR / 'g2g_level_dynamics_examples.csv', index=False)
ax.set(xlabel='Human minus mouse mean lognorm (normalized level)', ylabel='G2G alignment similarity',
       title='Level difference and trajectory alignment are separate')
savefig(fig, 'fig03_level_vs_alignment.png')

# %% [markdown]
# ## 6. G2G alignment-pattern clusters
#
# Cluster only G2G five-state strings with the package's normalized Levenshtein distance, using average-linkage hierarchical clustering. Identical strings share a node; the full gene count is used for cluster sizes and sampled silhouette diagnostics. Resolution is chosen before pathway files are opened. At most eight clusters are considered to keep the catalog interpretable.

# %%
strings = primary.alignment_string.astype(str).to_numpy()
unique, inverse, counts = np.unique(strings, return_inverse=True, return_counts=True)
if len(unique) < 3:
    raise ValueError('Too few distinct G2G strings to cluster.')
# O(U^2) exact distances among U unique strings; if this exceeds RAM, use a disk-backed
# condensed distance vector rather than subsampling genes from the biological universe.
condensed = np.empty(len(unique) * (len(unique) - 1) // 2, dtype=np.float64)
p = 0
for i in range(len(unique) - 1):
    for j in range(i + 1, len(unique)):
        condensed[p] = Levenshtein.distance(unique[i], unique[j]) / max(len(unique[i]), len(unique[j]))
        p += 1
tree = linkage(condensed, method='average')
probe = np.random.default_rng(RNG).choice(len(primary), size=min(900, len(primary)), replace=False)
probe_strings = strings[probe]
probe_distance = np.zeros((len(probe), len(probe)), dtype=np.float32)
for i in range(len(probe) - 1):
    for j in range(i + 1, len(probe)):
        probe_distance[i, j] = probe_distance[j, i] = Levenshtein.distance(probe_strings[i], probe_strings[j]) / max(len(probe_strings[i]), len(probe_strings[j]))
resolution = []
for k in range(2, min(8, len(unique) - 1) + 1):
    labels = fcluster(tree, k, criterion='maxclust')[inverse]
    sizes = pd.Series(labels).value_counts()
    sampled = labels[probe]
    score = silhouette_score(probe_distance, sampled, metric='precomputed') if 1 < len(np.unique(sampled)) < len(probe) else np.nan
    resolution.append({'k': k, 'silhouette': score, 'smallest_cluster': sizes.min(),
                       'largest_cluster_fraction': sizes.max() / len(labels)})
resolution = pd.DataFrame(resolution)
valid = resolution[resolution.smallest_cluster >= max(25, .01 * len(primary))]
if valid.empty:
    print('No cut meets the 1% minimum-cluster guide; inspect the rare alignment cluster in the catalog.')
    valid = resolution
best = valid.silhouette.max()
chosen_k = int(valid[valid.silhouette >= best - .02].k.min())
labels = fcluster(tree, chosen_k, criterion='maxclust')[inverse]
# Name clusters by decreasing size, before any biological description.
name_map = {old: f'G{i + 1}' for i, old in enumerate(pd.Series(labels).value_counts().index)}
primary['g2g_cluster'] = [name_map[x] for x in labels]
resolution['selected'] = resolution.k.eq(chosen_k)
resolution.to_csv(OUTPUT_DIR / 'g2g_cluster_resolution.csv', index=False)
assignments = primary[['gene', 'g2g_cluster', 'alignment_string', 'alignment_similarity']].copy()
assignments.to_csv(OUTPUT_DIR / 'g2g_cluster_assignments.csv', index=False)
display(resolution.round(3))


# %%
def state_summary(values):
    strings_here = list(values)
    chars = ''.join(strings_here)
    lengths = np.array([len(s) for s in strings_here])
    match = sum(ch in 'MVW' for ch in chars) / len(chars)
    warp = sum(ch in 'VW' for ch in chars) / len(chars)
    mismatch = sum(ch in 'ID' for ch in chars) / len(chars)
    positions = [i / max(len(s) - 1, 1) for s in strings_here for i, ch in enumerate(s) if ch in 'ID']
    location = 'none' if not positions else ('early' if np.median(positions) < 1 / 3 else 'late' if np.median(positions) > 2 / 3 else 'mid/mixed')
    consensus = pd.Series(strings_here).value_counts().index[0]  # observed modal string
    representative_match = sum(ch in 'MVW' for ch in consensus) / len(consensus)
    representative_warp = sum(ch in 'VW' for ch in consensus) / len(consensus)
    representative_mismatch = sum(ch in 'ID' for ch in consensus) / len(consensus)
    # Name the modal observed string; pooled state fractions remain separate catalog columns.
    descriptor = ('predominantly warped but matched' if representative_match >= .8 and representative_warp >= .3 else
                  'broadly conserved' if representative_match >= .8 and representative_mismatch < .1 else
                  'broadly divergent' if representative_mismatch >= .5 else
                  f'{location} mismatch' if location != 'none' else 'mixed alignment')
    return consensus, match, warp, mismatch, location, descriptor
assert state_summary(['MMMM', 'MMID'])[3] == 0.25
catalog_rows = []
for cluster, group in primary.groupby('g2g_cluster', sort=True):
    consensus, match, warp, mismatch, location, descriptor = state_summary(group.alignment_string)
    catalog_rows.append({'g2g_cluster': cluster, 'n_genes': len(group),
                         'median_alignment_similarity': group.alignment_similarity.median(),
                         'representative_alignment': consensus, 'fraction_match_states': match,
                         'fraction_warp_states': warp, 'fraction_mismatch_states': mismatch,
                         'mismatch_location': location, 'alignment_label': descriptor,
                         'median_human_minus_mouse_mean_lognorm': group.human_minus_mouse_mean_lognorm.median()})
catalog = pd.DataFrame(catalog_rows)
catalog.to_csv(OUTPUT_DIR / 'g2g_cluster_catalog.csv', index=False)
display(catalog)
# Fraction of each state across relative string position is a compact alignment clustergram.
state_order = 'MVWID'
profile = np.zeros((len(catalog), 20, len(state_order)))
for row, cluster in enumerate(catalog.g2g_cluster):
    for s in primary.loc[primary.g2g_cluster.eq(cluster), 'alignment_string']:
        for i, ch in enumerate(s):
            if ch in state_order:
                profile[row, min(19, int(20 * i / len(s))), state_order.index(ch)] += 1
profile /= np.maximum(profile.sum(axis=2, keepdims=True), 1)
fig, axes = plt.subplots(1, len(state_order), figsize=(14, max(3, .5 * len(catalog) + 1)), sharey=True)
for j, state in enumerate(state_order):
    image = axes[j].imshow(profile[:, :, j], aspect='auto', vmin=0, vmax=1, cmap='viridis')
    axes[j].set_title(state)
    axes[j].set_xlabel('relative alignment position')
    axes[j].set_yticks(range(len(catalog)), catalog.g2g_cluster)
fig.colorbar(image, ax=axes, label='state fraction', shrink=.7)
savefig(fig, 'fig04_alignment_pattern_clusters.png')

# %%
# Pick observed strings closest to each modal string; expression and pairwise consistency
# will be checked in the candidate tables below. Plot raw structures and G2G's interpolation.
representatives = {}
for row in catalog.itertuples():
    candidates = primary[primary.g2g_cluster.eq(row.g2g_cluster)].copy()
    candidates['distance_to_consensus'] = [Levenshtein.distance(s, row.representative_alignment) / max(len(s), len(row.representative_alignment)) for s in candidates.alignment_string]
    candidates = candidates.sort_values(['distance_to_consensus', 'pt_detected_fraction'], ascending=[True, False])
    representatives[row.g2g_cluster] = candidates.head(3).gene.tolist()
for cluster, selected in representatives.items():
    fit = RefQueryAligner(adata_mouse_g2g[:, selected].copy(), adata_human_g2g[:, selected].copy(), selected, G2G_BINS)
    fit.align_all_pairs()
    fig, axes = plt.subplots(len(selected), 3, figsize=(11, 2.8 * len(selected)), squeeze=False)
    for row, gene in enumerate(selected):
        result = fit.results_map[gene]
        for col, (a, color, title) in enumerate(((adata_mouse_g2g, COLORS['mouse'], 'mouse raw'),
                                                  (adata_human_g2g, COLORS['human'], 'human raw'))):
            values = dense_column(a.X, a.var_names.get_loc(gene))
            axes[row, col].scatter(a.obs.time, values, s=2, alpha=.08, color=color)
            axes[row, col].set(title=f'{gene}: {title}', xlabel='species PT time', ylabel='lognorm')
        axes[row, 2].plot(result.S.time_points, result.S.mean_trend, color=COLORS['mouse'], label='mouse G2G mean')
        axes[row, 2].plot(result.T.time_points, result.T.mean_trend, color=COLORS['human'], label='human G2G mean')
        axes[row, 2].fill_between(result.S.time_points, result.S.mean_trend - result.S.std_trend, result.S.mean_trend + result.S.std_trend, color=COLORS['mouse'], alpha=.12)
        axes[row, 2].fill_between(result.T.time_points, result.T.mean_trend - result.T.std_trend, result.T.mean_trend + result.T.std_trend, color=COLORS['human'], alpha=.12)
        axes[row, 2].set(title=f'G2G {result.alignment_str}', xlabel='PT time', ylabel='interpolated lognorm')
        axes[row, 2].title.set_fontsize(8)
        axes[row, 2].legend(fontsize=7)
    savefig(fig, f'fig05_representative_alignments_{cluster}.png')
    del fit

# %% [markdown]
# ## 7. Normal mouse PT program × G2G phenotype
#
# Use notebook 07's saved P assignments as annotation only. Its unavailable genes remain `unassigned`; P modules are not recomputed here.

# %%
if ATLAS_PATH.exists():
    atlas = pd.read_csv(ATLAS_PATH)
    if not {'gene', 'module'}.issubset(atlas) or atlas.gene.duplicated().any():
        raise ValueError('Unexpected or duplicated mouse positional assignment table.')
    primary = primary.merge(atlas[['gene', 'module']].rename(columns={'module': 'mouse_positional_module'}), on='gene', how='left', validate='one_to_one')
    primary['mouse_positional_module'] = primary.mouse_positional_module.fillna('unassigned')
    p_assigned = primary[primary.mouse_positional_module.isin([f'P{i}' for i in range(1, 6)])]
    table = pd.crosstab(p_assigned.mouse_positional_module, p_assigned.g2g_cluster)
    expected = np.outer(table.sum(axis=1), table.sum(axis=0)) / table.values.sum()
    oe = table.values / np.maximum(expected, 1e-12)
    cross = pd.DataFrame([{'mouse_positional_module': p, 'g2g_cluster': g, 'count': int(table.loc[p, g]),
                           'row_fraction': table.loc[p, g] / table.loc[p].sum(),
                           'observed_expected_ratio': oe[i, j]}
                          for i, p in enumerate(table.index) for j, g in enumerate(table.columns)])
    fig, ax = plt.subplots(figsize=(max(5, .7 * len(table.columns) + 2), max(3, .5 * len(table.index) + 1)))
    image = ax.imshow(np.log2(np.maximum(oe, 1e-6)), cmap='coolwarm', vmin=-2, vmax=2, aspect='auto')
    ax.set_xticks(range(len(table.columns)), table.columns)
    ax.set_yticks(range(len(table.index)), table.index)
    ax.set(xlabel='G2G alignment cluster', ylabel='mouse positional module', title='Mouse P program × G2G cluster (log2 O/E)')
    fig.colorbar(image, ax=ax, label='log2 observed / expected')
    savefig(fig, 'fig06_mouse_p_by_g2g_oe.png')
else:
    print(f'Optional notebook 07 atlas unavailable: {ATLAS_PATH}')
    primary['mouse_positional_module'] = 'unavailable'
    cross = pd.DataFrame(columns=['mouse_positional_module', 'g2g_cluster', 'count', 'row_fraction', 'observed_expected_ratio'])
cross.to_csv(OUTPUT_DIR / 'g2g_mouse_positional_crosstab.csv', index=False)
display(cross.head(15))

# %% [markdown]
# ## 8. Pathway interpretation
#
# Ranked pathway analysis is primary: a two-sided Mann–Whitney rank test asks whether members lie toward the divergent or conserved end of the **continuous** similarity distribution. Its p-values describe competitive gene-rank enrichment, not biological species replication. ORA is secondary and uses all G2G-eligible genes as background. Benjamini–Hochberg correction is applied over all tests in each analysis family. Pathway coherence reports dispersion and cluster entropy against size-matched random gene sets so a mixed pathway is not mislabeled as uniformly remodeled.

# %%
libraries = ('Reactome_2022', 'MSigDB_Hallmark_2020', 'KEGG_2019_Mouse')
pathway_dir = DATA_ROOT / 'mouse_vs_human' / 'pathway_gene_sets'
coverage = []
for library in libraries:
    path = pathway_dir / f'{library}.json'
    if not path.exists():
        raise FileNotFoundError(f'Required pathway library is missing: {path}')
    coverage.append(build_pathway_membership(json.loads(path.read_text()), adata.var_names,
                                             ortholog_map=orthologs, library_name=library,
                                             min_genes=10, tested=genes))
pathways = pd.concat(coverage, ignore_index=True)
pathways.to_csv(OUTPUT_DIR / 'g2g_pathway_coverage.csv', index=False)
pathways = pathways[pathways.retained].copy()
if pathways.empty:
    raise ValueError('No pathway has ten eligible genes.')
statistics = primary.set_index('gene')
rank_rows, ora_rows, coherence_rows = [], [], []
for row in pathways.itertuples():
    members = list(row.genes_present)
    inside = statistics.loc[members, 'alignment_similarity'].to_numpy(dtype=float)
    outside = statistics.loc[~statistics.index.isin(members), 'alignment_similarity'].to_numpy(dtype=float)
    u, p = mannwhitneyu(inside, outside, alternative='two-sided')
    auc = u / (len(inside) * len(outside))
    direction = 'conserved' if auc > .5 else 'divergent'
    rank_rows.append({'library': row.library, 'pathway': row.pathway, 'n_eligible_genes': len(members),
                      'median_similarity': np.median(inside), 'rank_auc': auc,
                      'rank_direction': direction, 'p_value': p})
    cluster_counts = statistics.loc[members, 'g2g_cluster'].value_counts().reindex(catalog.g2g_cluster, fill_value=0)
    probs = cluster_counts.to_numpy() / len(members)
    coherence_rows.append({'library': row.library, 'pathway': row.pathway, 'n_eligible_genes': len(members),
                           'median_similarity': np.median(inside), 'similarity_iqr': np.subtract(*np.percentile(inside, [75, 25])),
                           'cluster_entropy': entropy(probs),
                           'normalized_cluster_entropy': entropy(probs) / np.log(len(catalog)) if len(catalog) > 1 else 0,
                           'cluster_counts': json.dumps(cluster_counts.to_dict())})
    for cluster, group in primary.groupby('g2g_cluster', sort=True):
        n_overlap = len(set(members) & set(group.gene))
        ora_rows.append({'library': row.library, 'pathway': row.pathway, 'g2g_cluster': cluster,
                         'n_eligible_genes': len(members), 'n_cluster_genes': len(group),
                         'n_overlap': n_overlap,
                         'observed_expected_ratio': n_overlap / (len(members) * len(group) / len(genes)),
                         'p_value': hypergeom.sf(n_overlap - 1, len(genes), len(members), len(group))})
ranked_pathways = pd.DataFrame(rank_rows)
ranked_pathways['q_value'] = multipletests(ranked_pathways.p_value, method='fdr_bh')[1]
ranked_pathways = ranked_pathways.sort_values(['q_value', 'p_value'])
cluster_pathways = pd.DataFrame(ora_rows)
cluster_pathways['q_value'] = multipletests(cluster_pathways.p_value, method='fdr_bh')[1]
cluster_pathways = cluster_pathways.sort_values(['q_value', 'p_value'])
coherence = pd.DataFrame(coherence_rows)
# A size-matched null accounts for the dominant-cluster frequency and small set sizes.
cluster_codes = pd.Categorical(primary.g2g_cluster, categories=catalog.g2g_cluster).codes
null_rng = np.random.default_rng(RNG)
null_entropy = {}
for n in coherence.n_eligible_genes.unique():
    null_entropy[n] = np.array([entropy(np.bincount(cluster_codes[null_rng.choice(len(genes), size=n, replace=False)], minlength=len(catalog)) / n)
                                for _ in range(200)])
coherence['size_matched_null_median_entropy'] = [np.median(null_entropy[n]) for n in coherence.n_eligible_genes]
coherence['entropy_below_null_fraction'] = [np.mean(null_entropy[n] >= observed) for n, observed in zip(coherence.n_eligible_genes, coherence.cluster_entropy)]
ranked_pathways.to_csv(OUTPUT_DIR / 'g2g_ranked_pathway_enrichment.csv', index=False)
cluster_pathways.to_csv(OUTPUT_DIR / 'g2g_cluster_pathway_enrichment.csv', index=False)
coherence.to_csv(OUTPUT_DIR / 'g2g_pathway_coherence.csv', index=False)
display(ranked_pathways.head(15))
display(cluster_pathways.head(15))

# %%
top = (ranked_pathways.groupby('rank_direction', group_keys=False)
       .head(8).sort_values('rank_auc'))
fig, ax = plt.subplots(figsize=(9, max(4, .38 * len(top) + 1)))
ax.barh([f'{r.library}: {r.pathway}' for r in top.itertuples()],
        top.rank_auc - .5, color=np.where(top.rank_auc < .5, '#D55E00', '#0072B2'))
ax.axvline(0, color='black', linewidth=.8)
ax.set(xlabel='Pathway gene rank AUC − 0.5 (divergent ← → conserved)',
       title='Top continuous G2G similarity pathway shifts (BH across libraries)')
savefig(fig, 'fig07_ranked_pathway_enrichment.png')

# %% [markdown]
# ## 9. Specimen-pair robustness
#
# Four technical comparisons use the **same species-specific pseudotime** but one mouse specimen and one human section at a time. Each gene keeps four similarities, their median, range, and SD. The two human sections are one donor and never count as independent biological replication.

# %%
pairs = {}
for mouse in SPECIMENS['mouse']:
    for human in SPECIMENS['human']:
        label = f'{mouse}__{human}'
        ref = adata_mouse_g2g[adata_mouse_g2g.obs['sample'].astype(str).eq(mouse)].copy()
        query = adata_human_g2g[adata_human_g2g.obs['sample'].astype(str).eq(human)].copy()
        pairs[label] = align_g2g(ref, query, label).set_index('gene')['alignment_similarity']
consistency = primary[['gene', 'g2g_cluster', 'alignment_similarity']].rename(columns={'alignment_similarity': 'primary_alignment_similarity'}).set_index('gene')
for label, values in pairs.items():
    consistency[label] = values
pair_cols = list(pairs)
consistency['pairwise_median_similarity'] = consistency[pair_cols].median(axis=1)
consistency['pairwise_range_similarity'] = consistency[pair_cols].max(axis=1) - consistency[pair_cols].min(axis=1)
consistency['pairwise_sd_similarity'] = consistency[pair_cols].std(axis=1, ddof=0)
consistency = consistency.reset_index()
consistency.to_csv(OUTPUT_DIR / 'g2g_specimen_pair_consistency.csv', index=False)
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
axes[0].scatter(consistency.primary_alignment_similarity, consistency.pairwise_median_similarity, s=6, alpha=.2)
axes[0].plot([0, 1], [0, 1], color='black', linestyle='--', linewidth=.8)
axes[0].set(xlabel='combined G2G similarity', ylabel='median pairwise similarity')
axes[1].boxplot([consistency.loc[consistency.g2g_cluster.eq(g), 'pairwise_range_similarity'] for g in catalog.g2g_cluster], tick_labels=catalog.g2g_cluster)
axes[1].set(ylabel='range across four pairs', xlabel='G2G cluster')
savefig(fig, 'fig08_specimen_pair_consistency.png')

# %% [markdown]
# ## 10. Biological summary
#
# Candidate lists require detection in **each of the four specimens/sections** (at least 5% of PT structures) and agreement across all four technical pairings. These are conservative follow-up priorities, not statistical discoveries. Genes with little normalized mean-level difference remain eligible for the divergent list.
#
# For the follow-up tables only, raw ten-bin mean amplitude must reach 0.1 lognorm in at least one species for divergent candidates and in both species for conserved-dynamics candidates. This check is downstream of G2G and changes no gene alignment or cluster.

# %%
mouse_detect = np.asarray((adata_mouse_g2g.X > 0).sum(axis=0)).ravel() / adata_mouse_g2g.n_obs
human_detect = np.asarray((adata_human_g2g.X > 0).sum(axis=0)).ravel() / adata_human_g2g.n_obs
candidate = primary.merge(consistency.drop(columns=['g2g_cluster']), on='gene', validate='one_to_one')
assert candidate.gene.tolist() == genes, 'Detection vectors require the original gene order.'
candidate['mouse_detected_fraction'] = mouse_detect
candidate['human_detected_fraction'] = human_detect
per_section_detection = []
for a in (adata_mouse_g2g, adata_human_g2g):
    for sample in a.obs['sample'].astype(str).unique():
        subset = a.X[a.obs['sample'].astype(str).eq(sample).to_numpy()]
        per_section_detection.append(np.asarray((subset > 0).sum(axis=0)).ravel() / subset.shape[0])
candidate['min_specimen_detected_fraction'] = np.min(per_section_detection, axis=0)
adequate = candidate.min_specimen_detected_fraction.ge(.05)
for species, a in (('mouse', adata_mouse_g2g), ('human', adata_human_g2g)):
    bins = np.minimum((a.obs.time.to_numpy(dtype=float) * 10).astype(int), 9)
    bin_means = np.vstack([np.asarray(a.X[bins == b].mean(axis=0)).ravel() for b in range(10) if np.any(bins == b)])
    candidate[f'{species}_raw_decile_amplitude_lognorm'] = np.ptp(bin_means, axis=0)
candidate['representative_alignment_type'] = candidate.g2g_cluster.map(catalog.set_index('g2g_cluster').alignment_label)
candidate['pairwise_consistency'] = candidate.pairwise_range_similarity
# Presentation-only confidence bands: primary and every pair must agree on the same side.
divergent = candidate[adequate & candidate[['mouse_raw_decile_amplitude_lognorm', 'human_raw_decile_amplitude_lognorm']].max(axis=1).ge(.1) & candidate.alignment_similarity.lt(.5) & candidate[pair_cols].lt(.5).all(axis=1)].sort_values(['alignment_similarity', 'pairwise_range_similarity'])
conserved = candidate[adequate & candidate[['mouse_raw_decile_amplitude_lognorm', 'human_raw_decile_amplitude_lognorm']].min(axis=1).ge(.1) & candidate.alignment_similarity.ge(.8) & candidate[pair_cols].ge(.8).all(axis=1)].sort_values(['alignment_similarity', 'pairwise_range_similarity'], ascending=[False, True])
fields = ['gene', 'g2g_cluster', 'alignment_similarity', 'representative_alignment_type',
          'human_minus_mouse_mean_lognorm', 'mouse_positional_module', 'mouse_detected_fraction',
          'human_detected_fraction', 'min_specimen_detected_fraction', 'mouse_raw_decile_amplitude_lognorm', 'human_raw_decile_amplitude_lognorm', 'pairwise_median_similarity', 'pairwise_range_similarity', *pair_cols]
divergent[fields].to_csv(OUTPUT_DIR / 'g2g_high_confidence_divergent_genes.csv', index=False)
conserved[fields].to_csv(OUTPUT_DIR / 'g2g_high_confidence_conserved_genes.csv', index=False)
print(f'High-confidence divergent: {len(divergent)}; conserved: {len(conserved)}')
display(divergent[fields].head(15))
display(conserved[fields].head(15))

# %%
print('A. Overall trajectory conservation: median similarity', round(float(np.median(sim)), 3),
      '; IQR', tuple(np.round(np.quantile(sim, [.25, .75]), 3)))
print('B. Major mismatch patterns:')
display(catalog[['g2g_cluster', 'n_genes', 'alignment_label', 'median_alignment_similarity', 'mismatch_location']])
print('C. Mouse positional programs:')
if cross.empty:
    print('P atlas unavailable; run notebook 07 to add this bridge.')
else:
    display(cross[cross.mouse_positional_module.ne('unassigned')].sort_values('observed_expected_ratio', ascending=False).head(10))
print('D. Divergent pathways with coherent alignment-cluster membership:')
coherent_divergent = ranked_pathways[ranked_pathways.rank_direction.eq('divergent') & ranked_pathways.q_value.lt(.05)].merge(
    coherence[['library', 'pathway', 'similarity_iqr', 'normalized_cluster_entropy', 'entropy_below_null_fraction']],
    on=['library', 'pathway'])
coherent_divergent = coherent_divergent[coherent_divergent.entropy_below_null_fraction.ge(.95)]
display(coherent_divergent.sort_values(['q_value', 'normalized_cluster_entropy']).head(10))
if coherent_divergent.empty:
    print('No pathway meets both the exploratory ranked-enrichment and size-matched coherence criteria.')
print('E. Divergent candidates with comparatively small normalized mean-level differences:')
display(divergent.iloc[np.argsort(np.abs(divergent.human_minus_mouse_mean_lognorm.to_numpy()))[:10]][fields])
print('F. Four-pair technical consistency: median range', round(consistency.pairwise_range_similarity.median(), 3))
display(consistency.groupby('g2g_cluster').agg(n_genes=('gene', 'size'),
    median_pairwise_range=('pairwise_range_similarity', 'median'),
    median_pairwise_similarity=('pairwise_median_similarity', 'median')).reset_index())
print('The two human sections are healthy cortex from one donor. These are descriptive patterns, not population-level species tests.')
print('Representative genes across the four technical pairings:')
rep_genes = [gene for selected in representatives.values() for gene in selected]
display(consistency[consistency.gene.isin(rep_genes)][['gene', 'g2g_cluster', 'primary_alignment_similarity', *pair_cols, 'pairwise_range_similarity']])

# %% [markdown]
# ## 11. Fixed-pseudospace comparison on the shared Harmony PT coordinate
#
# **Sensitivity assumption.** Notebook 03's already oriented, shared PT DPT represents corresponding mouse and human PT positions. Here mouse $s_j$ is compared only with human $s_j$ on one common grid. No coordinate is refitted and no gene-specific alignment is performed. The same 9,904 G2G-eligible genes, log1p-normalized expression, specimens, and pathway libraries are retained. The two human cortex sections come from one donor; four pairings are spatial/technical sensitivity checks.
#
# G2G v0.2.0's `DP5.compute_cell` computes an MML compression from *generated* distributions inside its DP object. Isolating it would initialize the alignment machinery and make the local statistic depend on synthetic draws. This section therefore uses a **G2G-inspired fixed-position distribution comparison**: the same Gaussian kernel form as G2G's nonadaptive interpolator (width 0.1 on the original [0,1] coordinate), weighted local means and variances, and symmetric Gaussian KL (Jeffreys divergence). This is **not** the exact G2G MML score. Higher divergence means greater mismatch at the same position. A variance floor of 0.05 lognorm units prevents almost constant genes from producing singular Gaussian costs; it is a numerical regularizer, not a significance threshold. All 9,904 genes stay in the analysis.

# %%
# Exact pairing to the existing primary result; fail before any new calculation if it drifts.
assert len(genes) == 9904 and primary.gene.tolist() == genes
assert adata_mouse_g2g.var_names.tolist() == genes == adata_human_g2g.var_names.tolist()
coord_name = 'shared_pseudospace'
if coord_name not in adata.obs:
    raise ValueError('Notebook 03 shared PT coordinate is absent.')
shared = adata.obs[coord_name].astype(float)
if not np.isfinite(shared).all() or shared.nunique() < 15:
    raise ValueError('Invalid or collapsed shared PT coordinate.')
if shared.min() < 0 or shared.max() > 1:
    shared = (shared - shared.min()) / (shared.max() - shared.min())
for a in (adata_mouse_g2g, adata_human_g2g):
    a.obs['fixed_s'] = shared.loc[a.obs_names].to_numpy(dtype=float)
assert np.isfinite(adata_mouse_g2g.obs.fixed_s).all() and np.isfinite(adata_human_g2g.obs.fixed_s).all()

# Common support: intersection of the pooled species p5–p95 intervals. The
# 14 equally spaced positions use G2G's interpolation resolution. No tail with
# sparse representation in either species enters the comparison.
mouse_s = adata_mouse_g2g.obs.fixed_s.to_numpy(dtype=float)
human_s = adata_human_g2g.obs.fixed_s.to_numpy(dtype=float)
support = (max(np.quantile(mouse_s, .05), np.quantile(human_s, .05)),
           min(np.quantile(mouse_s, .95), np.quantile(human_s, .95)))
if support[1] <= support[0] or support[1] - support[0] < .1:
    raise ValueError(f'Insufficient shared PT support: {support}')
fixed_grid = np.linspace(*support, G2G_BINS)

def kernel_weights(a):
    s = a.obs.fixed_s.to_numpy(dtype=float)
    w = np.exp(-((fixed_grid[:, None] - s[None, :]) / .1) ** 2)
    w /= w.sum(axis=1, keepdims=True)
    neff = 1 / np.square(w).sum(axis=1)
    return w, neff

for species, a in (('mouse', adata_mouse_g2g), ('human', adata_human_g2g)):
    _, neff = kernel_weights(a)
    if neff.min() < 20:
        raise ValueError(f'{species} has too few effective structures on common support: {neff.min():.1f}')
    print(f'{species}: shared s p5–p95 [{np.quantile(a.obs.fixed_s,.05):.3f}, {np.quantile(a.obs.fixed_s,.95):.3f}], minimum effective N {neff.min():.1f}')
pd.DataFrame({'grid_index': np.arange(len(fixed_grid)), 'fixed_grid': fixed_grid}).to_csv(OUTPUT_DIR / 'fixed_grid.csv', index=False)
print(f'Common p5–p95 intersection: {support}; K={len(fixed_grid)}')

# %% [markdown]
# ### 11.1. Check shared-coordinate orientation in each specimen
#
# Mouse markers are the panel already used upstream. Human symbols are translated through the accepted ortholog map, and unavailable genes are reported. Conserved anchors are a small prespecified early/late candidate list; an anchor is retained only if its observed Spearman direction agrees in **both** species. This validates plausibility, not exact anatomical homology. A weak or reversed human trend warns but does not suppress the sensitivity analysis.

# %%
from scipy.stats import pearsonr
mouse_stage = {'early': MARKERS['S1'],
               'mid': ('Slc22a6', 'Slc13a3', 'Acsm3'),
               'late': MARKERS['S3']}
human_stage_symbols = {'early': ('SLC5A2','SLC5A12','SLC6A19','SLC7A7','SLC4A4'),
                       'mid': HUMAN_TRIAL_MARKERS['S2'], 'late': HUMAN_TRIAL_MARKERS['S3']}
map_human = orthologs.drop_duplicates('human_symbol').set_index('human_symbol').mouse_symbol
human_stage = {stage: tuple(map_human.loc[g] for g in symbols if g in map_human.index and map_human.loc[g] in genes)
               for stage, symbols in human_stage_symbols.items()}
anchor_candidates = {'early': ('Slc5a2','Slc5a12'), 'late': ('Slc22a7','Slc7a13')}
anchor_stage = {}
for stage, symbols in anchor_candidates.items():
    direction = -1 if stage == 'early' else 1
    anchor_stage[stage] = tuple(g for g in symbols if g in genes and all(
        np.isfinite(spearmanr(a.obs.fixed_s, dense_column(a.X, a.var_names.get_loc(g))).statistic)
        and direction * spearmanr(a.obs.fixed_s, dense_column(a.X, a.var_names.get_loc(g))).statistic > .1
        for a in (adata_mouse_g2g, adata_human_g2g)))
print('Human marker mapping:', human_stage, '; direction-supported conserved anchors:', anchor_stage)

validation_rows = []
fig, axes = plt.subplots(4, 3, figsize=(12, 12), sharex=True)
for r, (species, a, stages) in enumerate((('mouse', adata_mouse_g2g, mouse_stage),
                                         ('human', adata_human_g2g, human_stage))):
    for sample in SPECIMENS[species]:
        sub = a[a.obs['sample'].astype(str).eq(sample)]
        row = r * 2 + SPECIMENS[species].index(sample)
        for col, stage in enumerate(('early','mid','late')):
            markers = [g for g in stages[stage] if g in genes]
            if not markers:
                raise ValueError(f'No assayed {species} {stage} PT markers.')
            x = sub[:, markers].X
            x = x.toarray() if sparse.issparse(x) else np.asarray(x)
            std = x.std(axis=0)
            score = ((x[:, std > 0] - x[:, std > 0].mean(axis=0)) / std[std > 0]).mean(axis=1)
            s = sub.obs.fixed_s.to_numpy(dtype=float)
            rho = spearmanr(s, score).statistic
            validation_rows.append({'species': species, 'sample': sample, 'stage': stage,
                                    'markers': ','.join(markers), 'n_markers': len(markers), 'spearman_rho': rho})
            axes[row, col].scatter(s, score, s=2, alpha=.08, color=COLORS[species])
            edges = np.linspace(s.min(), s.max(), 16)
            centers = (edges[:-1] + edges[1:]) / 2
            trend = [np.median(score[(s >= lo) & (s < hi)]) if np.any((s >= lo) & (s < hi)) else np.nan
                     for lo, hi in zip(edges[:-1], edges[1:])]
            axes[row, col].plot(centers, trend, color='black', lw=1.2)
            axes[row, col].set(title=f'{sample}: {stage} (ρ={rho:.2f})', xlabel='shared PT s', ylabel='marker z-score')
validation = pd.DataFrame(validation_rows)
validation.to_csv(OUTPUT_DIR / 'fixed_shared_coordinate_validation.csv', index=False)
savefig(fig, 'fig09_fixed_coordinate_validation.png')
display(validation)
human_qc = validation[validation.species.eq('human')].pivot(index='sample', columns='stage', values='spearman_rho')
if (human_qc.early.ge(0) | human_qc.late.le(0)).any():
    print('WARNING: Human early→late marker ordering is inconsistent. Fixed-coordinate results are sensitivity analysis only; anatomical correspondence is uncertain.')
if len(anchor_stage['early']) < 2 or len(anchor_stage['late']) < 2:
    print('WARNING: Fewer than two direction-supported conserved anchors in one stage; anatomical correspondence remains provisional.')
if human_qc.mid.lt(0).all():
    print('WARNING: Human S2 marker scores lack a clear intermediate trend in the validation plot. Early/late orientation is supported, but full S1→S2→S3 correspondence is uncertain.')

# %% [markdown]
# ### 11.2. Local distributions and diagonal mismatch
#
# For each species, the 14 × structures weight matrix uses G2G's nonadaptive Gaussian kernel form. Sparse matrix products yield the weighted first and second expression moments for **every** eligible gene. Effective sample size is $1/\sum_i w_i^2$. The mismatch is half the two directional KL divergences between local Gaussian approximations, evaluated only at matching grid indices. These are descriptive distribution distances; the kernel treats individual structures as observations, and the pooled statistic does not represent donor-level replication.

# %%
FIXED_SD_FLOOR = .05

def local_fixed(a):
    w, neff = kernel_weights(a)
    x = a.X.tocsr() if sparse.issparse(a.X) else sparse.csr_matrix(a.X)
    mu = np.asarray(w @ x)
    second = np.asarray(w @ x.power(2))
    variance = np.maximum(second - mu**2, 0)
    return mu, np.sqrt(variance), neff

def fixed_distance(mouse, human):
    mu_m, sd_m, n_m = local_fixed(mouse)
    mu_h, sd_h, n_h = local_fixed(human)
    vm = np.maximum(sd_m**2, FIXED_SD_FLOOR**2)
    vh = np.maximum(sd_h**2, FIXED_SD_FLOOR**2)
    # Symmetric Gaussian KL / Jeffreys divergence: equal distributions give 0.
    d = .25 * (vm/vh + vh/vm - 2 + (mu_m-mu_h)**2 * (1/vm + 1/vh))
    if not np.isfinite(d).all() or (d < -1e-10).any():
        raise ValueError('Fixed-position distance is invalid.')
    return np.maximum(d, 0), (mu_m, sd_m, n_m), (mu_h, sd_h, n_h)

# Runnable invariants: symmetry, zero self-distance, and strictly diagonal grid comparison.
_demo = np.array([[0., 1.]])
_demo_dist, _, _ = fixed_distance(adata_mouse_g2g[:, genes[:2]], adata_mouse_g2g[:, genes[:2]])
assert np.allclose(_demo_dist, 0, atol=1e-10)
fixed_profile, mouse_local, human_local = fixed_distance(adata_mouse_g2g, adata_human_g2g)
assert fixed_profile.shape == (len(fixed_grid), len(genes))
np.savez_compressed(OUTPUT_DIR / 'fixed_local_distributions.npz', fixed_grid=fixed_grid,
                    genes=np.asarray(genes), mu_mouse=mouse_local[0], sd_mouse=mouse_local[1], neff_mouse=mouse_local[2],
                    mu_human=human_local[0], sd_human=human_local[1], neff_human=human_local[2])
profile_table = pd.DataFrame(fixed_profile.T, columns=[f's{j:02d}' for j in range(len(fixed_grid))])
profile_table.insert(0, 'gene', genes)
profile_table.to_csv(OUTPUT_DIR / 'fixed_pseudospace_mismatch_profiles.csv', index=False)
regions = np.array_split(np.arange(len(fixed_grid)), 3)
fixed = pd.DataFrame({'gene': genes, 'fixed_overall_mismatch': fixed_profile.mean(axis=0),
                      'fixed_similarity': -fixed_profile.mean(axis=0),
                      'fixed_max_mismatch': fixed_profile.max(axis=0),
                      'fixed_max_mismatch_position': fixed_grid[fixed_profile.argmax(axis=0)],
                      'fixed_early_mismatch': fixed_profile[regions[0]].mean(axis=0),
                      'fixed_mid_mismatch': fixed_profile[regions[1]].mean(axis=0),
                      'fixed_late_mismatch': fixed_profile[regions[2]].mean(axis=0)})
fixed = fixed.merge(candidate[['gene','human_minus_mouse_mean_lognorm','mouse_detected_fraction',
    'human_detected_fraction','min_specimen_detected_fraction','mouse_raw_decile_amplitude_lognorm',
    'human_raw_decile_amplitude_lognorm']], on='gene', validate='one_to_one')
assert fixed.gene.tolist() == primary.gene.tolist()
fixed['level_matched'] = fixed.human_minus_mouse_mean_lognorm.abs().le(.15)
fixed.to_csv(OUTPUT_DIR / 'fixed_pseudospace_gene_statistics.csv', index=False)
print(f'Level-matched (|human–mouse mean lognorm| ≤ 0.15): {fixed.level_matched.sum():,}/{len(fixed):,}')
print('Fixed mismatch vs absolute mean level difference, Spearman ρ:',
      spearmanr(fixed.fixed_overall_mismatch, fixed.human_minus_mouse_mean_lognorm.abs()).statistic)
fig, ax = plt.subplots(figsize=(6,4))
ax.scatter(fixed.human_minus_mouse_mean_lognorm.abs(), fixed.fixed_overall_mismatch, s=4, alpha=.12)
ax.axvline(.15, color='black', ls='--')
ax.set(xlabel='|human − mouse mean lognorm|', ylabel='fixed overall mismatch (Jeffreys)', yscale='symlog', title='Level difference remains visible')
savefig(fig, 'fig10_fixed_mismatch_vs_level.png')

# %% [markdown]
# ### 11.3. Positional profiles and technical pairings
#
# Continuous early/mid/late burdens stay primary. Profiles are ordered by dominant third, then overall mismatch for display; no hard fixed clusters are imposed. Four specimen pairings use this **same** shared coordinate and grid. Pairwise ranks quantify stability within the one-donor design.

# %%
dominant = np.argmax(fixed[['fixed_early_mismatch','fixed_mid_mismatch','fixed_late_mismatch']].to_numpy(), axis=1)
fixed['dominant_fixed_region'] = np.array(['early','mid','late'])[dominant]
order = np.lexsort((-fixed.fixed_overall_mismatch.to_numpy(), dominant))
fig, ax = plt.subplots(figsize=(8,6))
image = ax.imshow(np.log1p(fixed_profile[:, order].T), aspect='auto', interpolation='nearest', cmap='magma',
                  extent=[fixed_grid[0],fixed_grid[-1],len(genes),0])
ax.set(xlabel='shared PT s', ylabel='genes ordered by dominant region and burden', title='Fixed-position mismatch profiles')
fig.colorbar(image, ax=ax, label='log(1 + Jeffreys)')
savefig(fig, 'fig11_fixed_mismatch_heatmap.png')

pair_fixed = {}
for mouse in SPECIMENS['mouse']:
    m = adata_mouse_g2g[adata_mouse_g2g.obs['sample'].astype(str).eq(mouse)]
    for human in SPECIMENS['human']:
        h = adata_human_g2g[adata_human_g2g.obs['sample'].astype(str).eq(human)]
        label = f'{mouse}__{human}'
        distance, _, _ = fixed_distance(m, h)
        pair_fixed[label] = distance.mean(axis=0)
fixed_pairs = fixed[['gene','fixed_overall_mismatch']].copy()
for label, values in pair_fixed.items():
    fixed_pairs[label] = values
fixed_pair_cols = list(pair_fixed)
fixed_pairs['pairwise_median_mismatch'] = fixed_pairs[fixed_pair_cols].median(axis=1)
fixed_pairs['pairwise_range_mismatch'] = fixed_pairs[fixed_pair_cols].max(axis=1) - fixed_pairs[fixed_pair_cols].min(axis=1)
fixed_pairs['pairwise_sd_mismatch'] = fixed_pairs[fixed_pair_cols].std(axis=1, ddof=0)
rank_matrix = fixed_pairs[fixed_pair_cols].rank(pct=True, ascending=True)
fixed_pairs['pairwise_min_divergence_percentile'] = rank_matrix.min(axis=1)
fixed_pairs['pairwise_rank_sd'] = rank_matrix.std(axis=1, ddof=0)
fixed_pairs.to_csv(OUTPUT_DIR / 'fixed_specimen_pair_consistency.csv', index=False)

# %% [markdown]
# ### 11.4. Same pathway libraries and rank-test framework
#
# For each pathway, a two-sided Mann–Whitney test compares member ranks with the other eligible genes. AUC above 0.5 means stronger fixed mismatch. BH correction covers all libraries and regions within each analysis family. The p and q values describe competitive **gene-rank** enrichment, not biological replication. Coherence uses the dominant positional third or a broad label assigned from continuous burdens; its matched-null entropy is descriptive.

# %%
fixed_stats = fixed.set_index('gene')
path_rows = []
for region, column in [('overall','fixed_overall_mismatch'), ('early','fixed_early_mismatch'),
                       ('mid','fixed_mid_mismatch'), ('late','fixed_late_mismatch')]:
    for row in pathways.itertuples():
        members = list(row.genes_present)
        inside = fixed_stats.loc[members, column].to_numpy(float)
        outside = fixed_stats.loc[~fixed_stats.index.isin(members), column].to_numpy(float)
        u, p = mannwhitneyu(inside, outside, alternative='two-sided')
        path_rows.append({'region':region,'library':row.library,'pathway':row.pathway,
            'n_eligible_genes':len(members),'rank_auc':u/(len(inside)*len(outside)),
            'rank_effect':u/(len(inside)*len(outside))-.5,
            'rank_direction':'divergent' if u/(len(inside)*len(outside))>.5 else 'conserved',
            'p_value':p,'median_fixed_mismatch':np.median(inside),
            'median_g2g_similarity':statistics.loc[members,'alignment_similarity'].median()})
fixed_path_all = pd.DataFrame(path_rows)
fixed_path_all['q_value'] = np.nan
for family, idx in fixed_path_all.groupby(fixed_path_all.region.eq('overall')).groups.items():
    fixed_path_all.loc[idx,'q_value'] = multipletests(fixed_path_all.loc[idx,'p_value'], method='fdr_bh')[1]
fixed_ranked = fixed_path_all[fixed_path_all.region.eq('overall')].drop(columns='region').sort_values('q_value')
fixed_positional = fixed_path_all[~fixed_path_all.region.eq('overall')].sort_values(['region','q_value'])
fixed_ranked.to_csv(OUTPUT_DIR / 'fixed_ranked_pathway_enrichment.csv', index=False)
fixed_positional.to_csv(OUTPUT_DIR / 'fixed_positional_pathway_enrichment.csv', index=False)

burden = fixed[['fixed_early_mismatch','fixed_mid_mismatch','fixed_late_mismatch']].to_numpy()
# Broad means comparable burden throughout PT. It is a descriptive phenotype.
shape = np.where(burden.max(axis=1) <= 1.5*np.maximum(burden.min(axis=1), 1e-12), 'broad',
                 fixed.dominant_fixed_region)
fixed['fixed_profile_shape'] = shape
shape_classes = ['early','mid','late','broad']
shape_codes = pd.Categorical(shape, categories=shape_classes).codes
null_rng = np.random.default_rng(RNG)
null_by_size = {}
coherence_fixed_rows = []
for row in pathways.itertuples():
    members = list(row.genes_present)
    loc = fixed_stats.index.get_indexer(members)
    values = fixed.fixed_overall_mismatch.to_numpy()[loc]
    counts = np.bincount(shape_codes[loc], minlength=4)
    observed_entropy = entropy(counts/counts.sum())
    n = len(members)
    if n not in null_by_size:
        null_by_size[n] = np.array([entropy(np.bincount(shape_codes[null_rng.choice(len(genes),n,replace=False)],minlength=4)/n)
                                    for _ in range(200)])
    regional = burden[loc].mean(axis=0)
    coherence_fixed_rows.append({'library':row.library,'pathway':row.pathway,'n_eligible_genes':n,
        'median_fixed_mismatch':np.median(values),'fixed_mismatch_iqr':np.subtract(*np.percentile(values,[75,25])),
        'dominant_fixed_region':shape_classes[np.argmax(regional)],
        **{f'fraction_{name}':counts[i]/n for i,name in enumerate(shape_classes)},
        'fixed_profile_entropy':observed_entropy,
        'fixed_profile_entropy_null_z':(observed_entropy-null_by_size[n].mean())/max(null_by_size[n].std(),1e-12)})
fixed_coherence = pd.DataFrame(coherence_fixed_rows)
fixed_coherence = fixed_coherence.merge(coherence[['library','pathway','normalized_cluster_entropy','similarity_iqr']],
                                      on=['library','pathway'],validate='one_to_one')
fixed_coherence.to_csv(OUTPUT_DIR / 'fixed_pathway_coherence.csv',index=False)
display(fixed_ranked.head(12))

# %% [markdown]
# ### 11.5. Gene phenotypes and ranking agreement
#
# The primary table keeps continuous statistics. Display categories use the fixed top/bottom quartiles and the existing G2G similarity guides (high ≥0.8, low <0.5); intermediate genes remain uncertain. We repeat key counts across fixed top 20–30% and G2G high 0.7–0.9 to expose threshold dependence. Candidate sets also require all four pairings to be strong, ≥5% detection in every specimen, and ≥0.1 lognorm decile amplitude in at least one species. None of these guides is a significance cutoff.
#
# **Ranking caveat:** G2G similarity has many exact ties at zero. Deterministic top-N overlap is reported as requested, with a tie-expanded overlap that includes every gene sharing the cutoff score. Interpret the tie-aware full-universe Spearman correlation first.

# %%
comparison = candidate.merge(fixed.drop(columns=['human_minus_mouse_mean_lognorm','mouse_detected_fraction',
    'human_detected_fraction','min_specimen_detected_fraction','mouse_raw_decile_amplitude_lognorm',
    'human_raw_decile_amplitude_lognorm']),on='gene',validate='one_to_one')
comparison = comparison.merge(fixed_pairs.drop(columns='fixed_overall_mismatch').rename(
    columns={name: f'fixed_{name}' for name in fixed_pair_cols}),on='gene',validate='one_to_one')
assert comparison.gene.tolist() == genes
strings = comparison.alignment_string.astype(str)
comparison['g2g_match_fraction'] = strings.map(lambda s: sum(c=='M' for c in s)/len(s))
comparison['g2g_warp_fraction'] = strings.map(lambda s: sum(c in 'VW' for c in s)/len(s))
comparison['g2g_mismatch_fraction'] = strings.map(lambda s: sum(c in 'ID' for c in s)/len(s))
q25,q75 = fixed.fixed_overall_mismatch.quantile([.25,.75])
comparison['fixed_strong'] = comparison.fixed_overall_mismatch.ge(q75)
comparison['fixed_low'] = comparison.fixed_overall_mismatch.le(q25)
comparison['g2g_high'] = comparison.alignment_similarity.ge(.8)
comparison['g2g_low'] = comparison.alignment_similarity.lt(.5)
comparison['fixed_pairwise_strong'] = comparison.pairwise_min_divergence_percentile.ge(.75)
comparison['technical_quality'] = comparison.min_specimen_detected_fraction.ge(.05) & comparison[[
    'mouse_raw_decile_amplitude_lognorm','human_raw_decile_amplitude_lognorm']].max(axis=1).ge(.1)
comparison['high_confidence_fixed_divergent'] = comparison.fixed_strong & comparison.fixed_pairwise_strong & comparison.technical_quality
conditions = [comparison.fixed_low & comparison.g2g_high,
              comparison.fixed_strong & comparison.g2g_high,
              comparison.fixed_strong & comparison.g2g_low,
              comparison.fixed_low & comparison.g2g_low]
labels = ['conserved by both','warp-rescuable positional difference','divergent by both','method-discordant / uncertain']
comparison['comparison_phenotype'] = np.select(conditions, labels, default='intermediate / uncertain')
comparison.to_csv(OUTPUT_DIR / 'g2g_vs_fixed_gene_comparison.csv',index=False)
comparison[['gene','comparison_phenotype','alignment_similarity','g2g_warp_fraction',
            'fixed_overall_mismatch','fixed_pairwise_strong','technical_quality','level_matched']].to_csv(
                OUTPUT_DIR / 'g2g_vs_fixed_phenotype_assignments.csv',index=False)
print('G2G similarity vs fixed mismatch, Pearson/Spearman:',
      pearsonr(comparison.alignment_similarity,comparison.fixed_overall_mismatch).statistic,
      spearmanr(comparison.alignment_similarity,comparison.fixed_overall_mismatch).statistic)
display(comparison.comparison_phenotype.value_counts().rename_axis('phenotype').reset_index(name='genes'))
threshold_rows=[]
for fixed_quantile in (.7,.75,.8):
    strong=comparison.fixed_overall_mismatch.ge(comparison.fixed_overall_mismatch.quantile(fixed_quantile))
    for high in (.7,.8,.9):
        conserved=comparison.alignment_similarity.ge(high)
        threshold_rows.append({'fixed_high_quantile':fixed_quantile,'g2g_high_similarity':high,
            'n_g2g_high':int(conserved.sum()),'fraction_g2g_high_fixed_strong':float(strong[conserved].mean()),
            'n_warp_rich_high':int((conserved & comparison.g2g_warp_fraction.ge(.3)).sum()),
            'fraction_warp_rich_high_fixed_strong':float(strong[conserved & comparison.g2g_warp_fraction.ge(.3)].mean())})
threshold_sensitivity=pd.DataFrame(threshold_rows)
threshold_sensitivity.to_csv(OUTPUT_DIR/'fixed_threshold_sensitivity.csv',index=False)
display(threshold_sensitivity)
fig,ax=plt.subplots(figsize=(7,5))
for name,group in comparison.groupby('comparison_phenotype'):
    ax.scatter(group.alignment_similarity,group.fixed_overall_mismatch,s=5,alpha=.22,label=name)
ax.axhline(q75,color='gray',ls='--');ax.axhline(q25,color='gray',ls=':')
ax.set(xlabel='G2G similarity',ylabel='fixed mismatch (Jeffreys)',yscale='symlog',title='Alignment freedom versus same-position difference')
ax.legend(fontsize=7,markerscale=2)
savefig(fig,'fig12_g2g_vs_fixed_genes.png')
cluster_fixed = comparison.groupby('g2g_cluster').agg(n_genes=('gene','size'),
    median_fixed=('fixed_overall_mismatch','median'),median_early=('fixed_early_mismatch','median'),
    median_mid=('fixed_mid_mismatch','median'),median_late=('fixed_late_mismatch','median'),
    median_warp_fraction=('g2g_warp_fraction','median')).reset_index()
display(cluster_fixed)
fig,ax=plt.subplots(figsize=(7,3.5))
ax.imshow(np.log1p(cluster_fixed[['median_early','median_mid','median_late']]),aspect='auto',cmap='magma')
ax.set(xticks=range(3),xticklabels=['early','mid','late'],yticks=range(len(cluster_fixed)),
       yticklabels=cluster_fixed.g2g_cluster,title='G2G cluster × fixed mismatch')
savefig(fig,'fig13_g2g_cluster_fixed_mismatch.png')

rank_rows=[]
for label,sub in [('all',comparison),('level_matched',comparison[comparison.level_matched])]:
    g2g_order=sub.sort_values(['alignment_similarity','gene'],kind='stable').gene.tolist()
    fixed_order=sub.sort_values(['fixed_overall_mismatch','gene'],ascending=[False,True],kind='stable').gene.tolist()
    for n in (100,250,500,1000):
        size=min(n,len(sub)); overlap=set(g2g_order[:size]) & set(fixed_order[:size])
        gr={g:i for i,g in enumerate(g2g_order)};fr={g:i for i,g in enumerate(fixed_order)}
        cutoff=sub.set_index('gene').loc[g2g_order[size-1],'alignment_similarity']
        tie_expanded=set(sub.loc[sub.alignment_similarity.le(cutoff),'gene'])
        tie_overlap=len(tie_expanded & set(fixed_order[:size]))
        rank_rows.append({'subset':label,'top_n_requested':n,'top_n_used':size,'overlap':len(overlap),
            'jaccard':len(overlap)/(2*size-len(overlap)),
            'g2g_cutoff_similarity':cutoff,'g2g_tie_expanded_size':len(tie_expanded),
            'tie_expanded_overlap':tie_overlap,
            'tie_expanded_jaccard':tie_overlap/(len(tie_expanded)+size-tie_overlap),
            'shared_gene_rank_spearman':spearmanr([gr[g] for g in overlap],[fr[g] for g in overlap]).statistic if len(overlap)>2 else np.nan})
rank_overlap=pd.DataFrame(rank_rows)
rank_overlap.to_csv(OUTPUT_DIR/'g2g_vs_fixed_top_gene_overlap.csv',index=False)
display(rank_overlap)
fig,ax=plt.subplots(figsize=(6,4))
for name,group in rank_overlap.groupby('subset'):
    ax.plot(group.top_n_requested,group.jaccard,marker='o',label=name)
ax.set(xlabel='top divergent genes',ylabel='Jaccard overlap',ylim=(0,1),title='G2G vs fixed divergent gene overlap')
ax.legend();savefig(fig,'fig14_g2g_fixed_gene_overlap.png')

# %% [markdown]
# ### 11.6. Pathway concordance and candidate programs
#
# G2G divergence effect is `0.5 − rank_auc`; fixed divergence effect is `rank_auc − 0.5`, so positive values have the same interpretation. ORA uses the unchanged 9,904-gene background. Candidate positional relocation requires fixed divergence, G2G matching, warp-rich alignment, technical consistency, detection, and amplitude. Non-alignable remodeling requires divergence in both methods and four-pair agreement. Both are descriptive follow-up sets.

# %%
path_compare=ranked_pathways.merge(fixed_ranked,on=['library','pathway','n_eligible_genes'],suffixes=('_g2g','_fixed'),validate='one_to_one')
path_compare['g2g_divergence_effect']=.5-path_compare.rank_auc_g2g
path_compare['fixed_divergence_effect']=path_compare.rank_auc_fixed-.5
path_compare['g2g_divergence_rank']=path_compare.g2g_divergence_effect.rank(ascending=False,method='min')
path_compare['fixed_divergence_rank']=path_compare.fixed_divergence_effect.rank(ascending=False,method='min')
path_compare['rank_difference']=path_compare.fixed_divergence_rank-path_compare.g2g_divergence_rank
path_compare.to_csv(OUTPUT_DIR/'g2g_vs_fixed_pathway_comparison.csv',index=False)
print('Pathway divergence-effect Spearman:',spearmanr(path_compare.g2g_divergence_effect,path_compare.fixed_divergence_effect).statistic)
fig,ax=plt.subplots(figsize=(7,6))
ax.scatter(path_compare.g2g_divergence_effect,path_compare.fixed_divergence_effect,s=8,alpha=.25)
for _,row in pd.concat([path_compare.nlargest(3,'g2g_divergence_effect'),
                        path_compare.nlargest(3,'fixed_divergence_effect'),
                        path_compare.nlargest(3,'rank_difference'),
                        path_compare.nsmallest(3,'rank_difference')]).drop_duplicates(['library','pathway']).iterrows():
    ax.annotate(row.pathway[:34],(row.g2g_divergence_effect,row.fixed_divergence_effect),fontsize=6)
ax.axhline(0,color='gray',lw=.7);ax.axvline(0,color='gray',lw=.7)
ax.set(xlabel='G2G pathway divergence effect',ylabel='fixed-coordinate pathway divergence effect',title='Pathway effect agreement')
savefig(fig,'fig15_g2g_fixed_pathway_effect.png')

# Prespecified themes, including possible discordance; exact pathway names are
# selected by substring, and every matching tested pathway is shown.
themes={'biological oxidations':'biological oxidation','xenobiotic metabolism':'xenobiotic',
'phase-II conjugation':'conjugation','glutathione metabolism':'glutathione',
'oxidative phosphorylation':'oxidative phosphorylation','amino-acid metabolism':'amino acid',
'fatty-acid metabolism':'fatty acid','small-molecule transport':'transport',
'rRNA processing':'rrna processing','RNA modification':'rna modification',
'chromatin/transcription':'chromatin|transcription'}
major_rows=[]
for theme,pattern in themes.items():
    hits=path_compare[path_compare.pathway.str.contains(pattern,case=False,regex=True)]
    for row in hits.itertuples():
        positional=fixed_positional[(fixed_positional.library==row.library)&(fixed_positional.pathway==row.pathway)]
        coh=fixed_coherence[(fixed_coherence.library==row.library)&(fixed_coherence.pathway==row.pathway)].iloc[0]
        major_rows.append({'theme':theme,'library':row.library,'pathway':row.pathway,
            'g2g_effect':row.g2g_divergence_effect,'g2g_q':row.q_value_g2g,
            'fixed_effect':row.fixed_divergence_effect,'fixed_q':row.q_value_fixed,
            **{f'{r.region}_effect':r.rank_effect for r in positional.itertuples()},
            'fixed_coherence_iqr':coh.fixed_mismatch_iqr,
            'interpretation':('concordant divergent' if row.g2g_divergence_effect>0 and row.fixed_divergence_effect>0 else
                              'concordant conserved' if row.g2g_divergence_effect<0 and row.fixed_divergence_effect<0 else 'method-discordant')})
major=pd.DataFrame(major_rows)
major.to_csv(OUTPUT_DIR/'g2g_vs_fixed_major_pathway_themes.csv',index=False)
display(major.sort_values(['theme','g2g_q']).groupby('theme').head(1))

quality=comparison.technical_quality & comparison.fixed_pairwise_strong
warp_set=comparison.loc[quality & comparison.fixed_strong & comparison.g2g_high & comparison.g2g_warp_fraction.ge(.3),'gene']
remodel_set=comparison.loc[quality & comparison.fixed_strong & comparison.g2g_low &
    comparison[pair_cols].lt(.5).all(axis=1),'gene']
def phenotype_ora(members, output, background=None):
    background=set(genes) if background is None else set(background)
    chosen=set(members) & background; rows=[]
    for row in pathways.itertuples():
        pathway_genes=set(row.genes_present) & background; overlap=len(chosen & pathway_genes)
        rows.append({'library':row.library,'pathway':row.pathway,'n_eligible_genes':len(pathway_genes),
            'n_selected':len(chosen),'n_overlap':overlap,
            'fold_enrichment':overlap/(len(chosen)*len(pathway_genes)/len(background)) if chosen and pathway_genes else np.nan,
            'p_value':hypergeom.sf(overlap-1,len(background),len(pathway_genes),len(chosen)) if chosen and pathway_genes else 1.})
    out=pd.DataFrame(rows);out['q_value']=multipletests(out.p_value,method='fdr_bh')[1]
    out.sort_values(['q_value','p_value']).to_csv(OUTPUT_DIR/output,index=False)
    return out.sort_values(['q_value','p_value'])
warp_path=phenotype_ora(warp_set,'warp_rescuable_pathway_enrichment.csv')
remodel_path=phenotype_ora(remodel_set,'nonalignable_remodeling_pathway_enrichment.csv')
print(f'Candidate warp-rescuable genes: {len(warp_set)}; candidate non-alignable genes: {len(remodel_set)}')
display(warp_path.head(8));display(remodel_path.head(8))

# %% [markdown]
# ### 11.7. Representative profiles and quantitative synthesis
#
# Examples are selected by high detection within each observed category, so they illustrate the distinction rather than serve as independent validation. The left panels use the **shared** fixed coordinate and local distribution moments. The right panels show the existing G2G species-specific trajectories and alignment string; only this display invokes G2G alignment. Conclusions below are calculated from observed output rather than specified in advance.

# %%
example_genes=[]
for label in labels:
    subset=comparison[comparison.comparison_phenotype.eq(label) & comparison.technical_quality]
    if not subset.empty:
        example_genes.append((label,subset.sort_values('min_specimen_detected_fraction',ascending=False).iloc[0].gene))
if example_genes:
    selected=[g for _,g in example_genes]
    fit=RefQueryAligner(adata_mouse_g2g[:,selected].copy(),adata_human_g2g[:,selected].copy(),selected,G2G_BINS)
    fit.align_all_pairs()  # visualization of the pre-existing G2G method only
    fig,axes=plt.subplots(len(selected),2,figsize=(11,3*len(selected)),squeeze=False)
    for r,(label,gene) in enumerate(example_genes):
        i=genes.index(gene)
        for local,color,name in ((mouse_local,COLORS['mouse'],'mouse'),(human_local,COLORS['human'],'human')):
            axes[r,0].plot(fixed_grid,local[0][:,i],color=color,label=name)
            axes[r,0].fill_between(fixed_grid,local[0][:,i]-local[1][:,i],local[0][:,i]+local[1][:,i],color=color,alpha=.12)
        result=fit.results_map[gene]
        axes[r,1].plot(result.S.time_points,result.S.mean_trend,color=COLORS['mouse'],label='mouse G2G')
        axes[r,1].plot(result.T.time_points,result.T.mean_trend,color=COLORS['human'],label='human G2G')
        axes[r,0].set(title=f'{gene}: {label}',xlabel='shared fixed PT s',ylabel='local lognorm')
        axes[r,1].set(title=f'G2G alignment: {result.alignment_str}',xlabel='species-specific PT time',ylabel='G2G lognorm')
        axes[r,0].legend(fontsize=7);axes[r,1].legend(fontsize=7)
    savefig(fig,'fig16_fixed_vs_g2g_representative_genes.png')
    del fit

level_sub=comparison[comparison.level_matched]
print('A. Gene ranking agreement: Spearman ρ =',round(spearmanr(comparison.alignment_similarity,
      comparison.fixed_overall_mismatch).statistic,3),'(negative means concordant divergence).')
display(rank_overlap)
print('B. Pathway divergence-effect Spearman ρ =',round(spearmanr(path_compare.g2g_divergence_effect,
      path_compare.fixed_divergence_effect).statistic,3))
print('Prespecified major theme directions:')
display(major.groupby(['theme','interpretation']).size().rename('tested_pathways').reset_index())
print('C. Fraction of G2G-high genes in fixed top quartile:',
      round(comparison.loc[comparison.g2g_high,'fixed_strong'].mean(),3),
      '; among warp-rich G2G-high genes:',
      round(comparison.loc[comparison.g2g_high & comparison.g2g_warp_fraction.ge(.3),'fixed_strong'].mean(),3))
print('D. Fraction of G2G-low genes strongly fixed-divergent and technically consistent:',
      round(comparison.loc[comparison.g2g_low,'high_confidence_fixed_divergent'].mean(),3))
print('Among the existing high-confidence G2G-divergent genes, fraction strongly fixed-divergent:',
      round(comparison.set_index('gene').loc[divergent.gene,'high_confidence_fixed_divergent'].mean(),3))
print('Among warp-rich G2G genes (V/W ≥0.3), fraction strongly fixed-divergent:',
      round(comparison.loc[comparison.g2g_warp_fraction.ge(.3),'fixed_strong'].mean(),3))
print('Level-matched G2G/fixed Spearman ρ:',
      round(spearmanr(level_sub.alignment_similarity,level_sub.fixed_overall_mismatch).statistic,3))
print('E. Level-matched genes:',len(level_sub),'; warp-rescuable fraction:',
      round((level_sub.comparison_phenotype=='warp-rescuable positional difference').mean(),3),
      '; divergent-by-both fraction:',round((level_sub.comparison_phenotype=='divergent by both').mean(),3))
display(level_sub.sort_values('fixed_overall_mismatch',ascending=False)[[
    'gene','alignment_similarity','fixed_overall_mismatch','g2g_warp_fraction','comparison_phenotype']].head(15))
level_top=set(level_sub.sort_values('fixed_overall_mismatch',ascending=False).head(max(100,int(.1*len(level_sub)))).gene)
display(phenotype_ora(level_top,'fixed_level_matched_top_pathway_enrichment.csv',background=level_sub.gene).head(12))
print('Interpretation: pathway effects and full-universe gene ranks are concordant, while top-N G2G lists are tie-sensitive. The four-pair candidate warp-rescuable set may be empty. Shared-coordinate validation limits anatomical claims; two human sections are one donor. All pathway tests remain exploratory.')

# %% [markdown]
# ### Observed synthesis from this run
#
# - **Gene rankings:** G2G similarity and fixed mismatch have Spearman ρ = −0.839 across 9,904 genes (−0.656 among 5,607 level-matched genes). Exact top-100 overlap is 10 genes, but all top-1,000 G2G ranks fall inside a 1,257-gene tie at similarity zero. The tie-expanded top-100 set contains 98 of the fixed top 100, so the raw top-N count cannot be read as strong biological disagreement.
# - **Pathways:** divergence-effect Spearman ρ = 0.855 across 1,508 tested pathways. Biological oxidations, xenobiotic/phase-II and glutathione metabolism, oxidative phosphorylation, amino-acid and fatty-acid metabolism remain divergent in both analyses. rRNA processing and RNA modification remain relatively conserved. Transport and chromatin/transcription themes contain mixed pathway-level directions; the complete, unfiltered comparison table retains these discordances.
# - **Warping sensitivity:** 17 of 4,163 G2G-high genes lie in the fixed top mismatch quartile (0.4%); none pass the additional warp-rich, four-pair, detection, and amplitude candidate rules. Under these descriptive guides, there is little evidence that broad G2G conservation is mainly rescued by warping. This does not rule out individual positionally shifted genes.
# - **Robust divergence:** 2,191 genes meet both broad divergence guides; 865 pass the stricter four-pair non-alignable candidate rules. Among the existing high-confidence G2G-divergent list, 46.4% are also strongly fixed-divergent with technical consistency. Mean level contributes to fixed mismatch (Spearman ρ = 0.690 with absolute mean lognorm difference), so level-matched results are essential context.
# - **Biological interpretation:** prohibiting warping does **not materially reverse the main pathway directions in this run**, though gene-level candidate membership and intensity vary. This is a sensitivity analysis under a shared-coordinate assumption. Human early and late PT markers have the expected directions, but human S2 markers show no clear intermediate trend and only one proposed early conserved anchor passes the directional check. Anatomical relocation remains provisional. Both human sections are healthy cortex from one donor; pathway q-values describe gene-rank tests, not species-level replication.

# %% [markdown]
# ## 12. Prioritizing trajectory-specific biology beyond bulk PT differences
#
# Run notebook 10 after notebook 09, then run this section. Notebook 10 uses the same eligible genes and pathway membership; this section reads its saved results **only after** all notebook-09 G2G and fixed-coordinate analyses. It does not refit, filter, or re-rank genes before discovery.
#
# "Strong" here means the top 10% of the tested pathway **effect ranks** within each method. Bulk uses the absolute signed pathway rank effect; trajectory uses the positive fixed-coordinate mismatch rank effect. The four labels (`bulk_supported`, `trajectory_specific`, `both`, `weak_in_both`) are mutually exclusive descriptive rank bins, not claims of absence or presence of a population-level species effect. The table also retains both methods' competitive q-values. Ranking outside bulk's top decile is a downstream priority view, not a new test family.
#
# The follow-up shortlist uses notebook 10's existing continuous-only candidate screen: top-decile fixed-coordinate signal, weak bulk and three-bin ranks/effects, adequate member detection, four-pair technical support, and size-matched low profile entropy. Notebook 09's existing G2G-cluster and fixed-profile coherence calculations are reported alongside member genes. A low profile-entropy null z indicates more shared mismatch location than size-matched random sets. G2G alignment and fixed-coordinate mismatch remain distinct measurements. The two human sections are healthy cortex from **one donor**; pathway tests are exploratory competitive gene-rank tests, not donor-level inference. These are candidates for biological follow-up, not validated mechanisms.

# %%
BULK_OUTPUT = RESULTS_ROOT / 'pt_continuous_vs_conventional'
needed = ('pathway_method_comparison.csv', 'continuous_only_pathway_candidates.csv',
          'gene_method_rank_comparison.csv')
missing = [name for name in needed if not (BULK_OUTPUT / name).exists()]
if missing:
    print(f'Run notebook 10, then rerun section 12; missing: {missing}')
else:
    keys = ['library', 'pathway', 'n_eligible_genes']
    bulk_path = pd.read_csv(BULK_OUTPUT / needed[0])
    bulk_candidates = pd.read_csv(BULK_OUTPUT / needed[1])
    bulk_genes = pd.read_csv(BULK_OUTPUT / needed[2])
    if bulk_path.duplicated(keys).any() or bulk_candidates.duplicated(keys).any():
        raise ValueError('Notebook 10 pathway keys must be unique.')
    if set(map(tuple, bulk_path[keys].to_numpy())) != set(map(tuple, path_compare[keys].to_numpy())):
        raise ValueError('Notebook 09/10 pathway universes differ; rerun notebook 10.')
    if bulk_genes.gene.isna().any() or bulk_genes.gene.duplicated().any() or set(bulk_genes.gene) != set(genes):
        raise ValueError('Notebook 09/10 eligible gene universes differ; rerun notebook 10.')
    bulk_fields = keys + ['bulk_rank_effect', 'bulk_abs_rank_effect', 'bulk_q_value',
                          'bulk_rank', 'bulk_ordinal', 'continuous_rank_effect',
                          'continuous_q_value', 'continuous_rank', 'continuous_ordinal',
                          'segment_rank', 'fixed_profile_entropy_null_z',
                          'fraction_technically_robust', 'fraction_adequately_detected']
    trajectory_bulk = path_compare.merge(bulk_path[bulk_fields], on=keys, validate='one_to_one')
    if not (np.allclose(trajectory_bulk.fixed_divergence_effect, trajectory_bulk.continuous_rank_effect) and
            np.allclose(trajectory_bulk.q_value_fixed, trajectory_bulk.continuous_q_value)):
        raise ValueError('Notebook 10 contains stale fixed-coordinate pathway results; rerun it.')
    top_n = max(1, round(len(trajectory_bulk) * .10))
    bulk_strong = trajectory_bulk.bulk_ordinal.le(top_n)
    trajectory_strong = trajectory_bulk.continuous_ordinal.le(top_n)
    trajectory_bulk['bulk_context'] = np.select(
        [bulk_strong & trajectory_strong, bulk_strong, trajectory_strong],
        ['both', 'bulk_supported', 'trajectory_specific'], default='weak_in_both')
    assert bulk_strong.sum() == trajectory_strong.sum() == top_n
    trajectory_bulk.to_csv(OUTPUT_DIR / 'g2g_bulk_pathway_context.csv', index=False)
    outside_bulk = trajectory_bulk.loc[~bulk_strong].sort_values(
        ['fixed_divergence_effect', 'library', 'pathway'], ascending=[False, True, True]).copy()
    outside_bulk.insert(0, 'trajectory_rank_without_bulk', np.arange(1, len(outside_bulk) + 1))
    outside_bulk.to_csv(OUTPUT_DIR / 'g2g_pathways_ranked_outside_bulk_top_decile.csv', index=False)
    display(trajectory_bulk.bulk_context.value_counts().rename_axis('bulk_context').reset_index(name='pathways'))
    display(outside_bulk[['trajectory_rank_without_bulk', 'library', 'pathway',
                          'fixed_divergence_effect', 'q_value_fixed', 'bulk_rank',
                          'bulk_abs_rank_effect', 'bulk_q_value', 'bulk_context']].head(15))

    follow_up = outside_bulk.merge(bulk_candidates[keys], on=keys, how='inner', validate='one_to_one')
    if not follow_up.bulk_context.eq('trajectory_specific').all():
        raise ValueError('Notebook 10 candidate screen disagrees with the rank categories; rerun notebook 10.')
    follow_up = follow_up.head(3).merge(
        fixed_coherence[['library', 'pathway', 'fixed_profile_entropy_null_z',
                         'fixed_mismatch_iqr', 'dominant_fixed_region']],
        on=['library', 'pathway'], suffixes=('', '_09'), validate='one_to_one').merge(
        coherence[['library', 'pathway', 'normalized_cluster_entropy',
                   'entropy_below_null_fraction']], on=['library', 'pathway'], validate='one_to_one')
    if not np.allclose(follow_up.fixed_profile_entropy_null_z,
                       follow_up.fixed_profile_entropy_null_z_09):
        raise ValueError('Notebook 10 contains stale pathway coherence; rerun it.')
    follow_up = follow_up.drop(columns='fixed_profile_entropy_null_z_09')
    follow_up.to_csv(OUTPUT_DIR / 'g2g_trajectory_specific_follow_up.csv', index=False)
    display(follow_up[['library', 'pathway', 'continuous_rank', 'bulk_rank', 'segment_rank',
                       'q_value_fixed', 'bulk_q_value', 'fixed_profile_entropy_null_z',
                       'entropy_below_null_fraction', 'fraction_technically_robust',
                       'fraction_adequately_detected']])

    fig, ax = plt.subplots(figsize=(7, 6))
    palette = {'both': '#6A3D9A', 'bulk_supported': '#D55E00',
               'trajectory_specific': '#0072B2', 'weak_in_both': '#999999'}
    for label, group in trajectory_bulk.groupby('bulk_context'):
        ax.scatter(group.bulk_rank, group.continuous_rank, s=12, alpha=.55,
                   color=palette[label], label=f'{label} ({len(group)})')
    ax.axvline(top_n, color='black', lw=.7, ls='--')
    ax.axhline(top_n, color='black', lw=.7, ls='--')
    for row in follow_up.itertuples():
        ax.annotate(row.pathway[:30], (row.bulk_rank, row.continuous_rank), fontsize=7)
    ax.set(xscale='log', yscale='log', xlabel='bulk pathway rank (absolute signed effect)',
           ylabel='fixed-coordinate divergence rank',
           title='PT pathway ranks: bulk versus trajectory')
    ax.invert_xaxis(); ax.invert_yaxis()
    ax.legend(fontsize=7)
    savefig(fig, 'fig17_bulk_vs_trajectory_pathway_ranks.png')

    if follow_up.empty:
        print('No pathway passed notebook 10 candidate screen; no biological shortlist assigned.')
    else:
        member_lookup = pathways.set_index(['library', 'pathway']).genes_present
        evidence = comparison.merge(bulk_genes[['gene', 'bulk_effect', 'bulk_abs_effect']],
                                    on='gene', validate='one_to_one').set_index('gene')
        member_rows = []
        for row in follow_up.itertuples():
            members = list(member_lookup.loc[(row.library, row.pathway)])
            if set(members) - set(evidence.index):
                raise ValueError(f'Missing member gene evidence for {row.pathway}.')
            detail = evidence.loc[members, ['bulk_effect', 'bulk_abs_effect',
                'alignment_similarity', 'g2g_warp_fraction', 'fixed_overall_mismatch',
                'fixed_early_mismatch', 'fixed_mid_mismatch', 'fixed_late_mismatch',
                'pairwise_min_divergence_percentile', 'min_specimen_detected_fraction']].reset_index()
            detail.insert(0, 'pathway', row.pathway)
            detail.insert(0, 'library', row.library)
            member_rows.append(detail)
        member_detail = pd.concat(member_rows, ignore_index=True)
        member_detail.to_csv(OUTPUT_DIR / 'g2g_trajectory_specific_member_evidence.csv', index=False)
        display(member_detail.sort_values('fixed_overall_mismatch', ascending=False).groupby(
            ['library', 'pathway'], sort=False).head(5))
        gene_to_index = {gene: i for i, gene in enumerate(genes)}
        fig, axes = plt.subplots(len(follow_up), 3, figsize=(12, 3 * len(follow_up)), squeeze=False)
        for r, row in enumerate(follow_up.itertuples()):
            detail = member_detail.loc[(member_detail.library.eq(row.library) &
                                        member_detail.pathway.eq(row.pathway) &
                                        member_detail.min_specimen_detected_fraction.ge(.05))]
            detail = detail.sort_values('fixed_overall_mismatch', ascending=False).head(3)
            for c, gene_row in enumerate(detail.itertuples()):
                i = gene_to_index[gene_row.gene]
                axes[r, c].plot(fixed_grid, mouse_local[0][:, i], color=COLORS['mouse'], label='mouse')
                axes[r, c].plot(fixed_grid, human_local[0][:, i], color=COLORS['human'], label='human')
                axes[r, c].set(title=(f'{gene_row.gene}: bulk {gene_row.bulk_effect:+.2f}, '
                                      f'fixed {gene_row.fixed_overall_mismatch:.2f}, '
                                      f'G2G {gene_row.alignment_similarity:.2f}'),
                               xlabel='shared PT s', ylabel='local lognorm')
            for c in range(len(detail), 3):
                axes[r, c].set_visible(False)
            axes[r, 0].annotate(row.pathway[:55], (0, 1), xycoords='axes fraction',
                                xytext=(0, 13), textcoords='offset points', fontsize=8)
        axes[0, 0].legend(fontsize=7)
        savefig(fig, 'fig18_trajectory_specific_pathway_member_curves.png')
        print('Shortlist is descriptive; inspect full member evidence and curves before biological interpretation.')
