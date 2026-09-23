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
# Notebook 02 uses a different mouse structure set, so a direct row-for-row transfer is unsafe. Recompute PT-only DPT within each species from lognorm expression, without cross-species Harmony. The S1/S3 marker axis selects an early root and orients DPT. The gate below checks each of the four sections separately.
#
# The existing `orient_and_normalize` helper uses the oriented DPT p5 and p95 as the 0 and 1 anchors and clips the tails. This robust within-species scaling is recorded because endpoint occupancy affects G2G interpolation.

# %%
from pseudospace.trajectory import orient_and_normalize

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
    # PCA is fitted within species; no joint embedding or cross-species coordinate enters DPT.
    sc.pp.pca(a, n_comps=30, random_state=RNG)
    sc.pp.neighbors(a, n_neighbors=30, use_rep='X_pca', random_state=RNG)
    sc.tl.diffmap(a)
    candidates = np.flatnonzero(axis <= np.quantile(axis, 0.05))
    position = a.obsm['X_pca'][candidates]
    root = candidates[np.linalg.norm(position - np.median(position, axis=0), axis=1).argmin()]
    a.uns['iroot'] = int(root)
    sc.tl.dpt(a)
    a.obs['time'] = orient_and_normalize(a.obs.dpt_pseudotime.to_numpy(), axis, min_valid=10)
    if not np.isfinite(a.obs.time).all() or a.obs.time.nunique() < 10:
        raise ValueError(f'{species}: invalid or collapsed DPT pseudotime.')
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
