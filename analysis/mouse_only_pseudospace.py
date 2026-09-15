# %% [markdown]
# # Mouse-only nephron pseudospace reconstruction
#
# Samples: `Ctrl1A2`, `Ctrl1A4` (Control) and `IR2A2`, `IR2A4` (AKI / ischemia-reperfusion).
#
# Notebook port of `6_mouse_only_pseudospace.py`; the code is identical and this notebook is
# run top-to-bottom from the repository root in the `pseudospace` kernel.
#
# ```
# Pipeline
#     Section 0  config: paths, parameters, marker panels
#     Section 1  pass-1 Harmony integration over the 4 mouse samples (theta = HARMONY_THETA)
#     Section 2  cell typing on the pass-1 embedding (ONE manual checkpoint), non-tubule removal,
#                pass-2 Harmony on the tubule subset (labels are NOT recomputed), Scanpy DPT
#     Section 3  marker heatmaps on the global DPT and on per-family subset DPTs
#     Section 4  healthy vs AKI level/shape decomposition over the PT cohort
#     Section 5  three-axis concordance (DPT vs marker axis vs distance-to-glomerulus)
#     Section 6  QuPath export + spatial validation hooks (opt-in)
#
# Differences from the v4 notebook
#     * NEPHRON_AXIS_MARKERS is the cell-typing panel and the source of the fine `segment_class`
#       vocabulary, the coarse rollup, and the early->late trajectory axis.
#     * Cell typing happens ONCE, on the pass-1 embedding. The pass-2 Harmony re-integrates the
#       tubule subset for a cleaner DPT graph but does not re-cluster or re-label anything.
#     * The injury / failed-repair diagnostics section is removed.
#     * The curated chart-derived heatmap panels are carried over verbatim, and all six of them
#       (PT, thin limb, TAL, DCT, collecting region, collecting-duct cell type) are plotted.
#
# Labels written to obs
#     segment_class               fine call, one of NEPHRON_SEGMENT_ORDER or a non-tubule class
#     coarse_class                rollup, one of COARSE_ORDER (or a non-tubule class in pass-1)
#     broad_tubule_marker_call    alias of coarse_class; the name the pseudospace package,
#                                 the QuPath export and the spatial validation script expect
#     celltype_primary/_secondary per-cell marker-score argmax (coarse vocabulary) - an
#                                 independent second opinion on the cluster-level label
#     celltype_*_segment          the same argmax at fine resolution
#     celltype_margin             score gap between the best and runner-up panel (confidence)
# ```

# %% [markdown]
# # Section 0 - configuration
#
# Every path, parameter and marker panel used below. Later sections reference these names only.

# %%
import os
import argparse
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings('ignore', category=FutureWarning)

# ---------------------------------------------------------------------------
# Cache redirects. These MUST be set before matplotlib (and, transitively, numba via umap) is
# imported -- both read their cache-dir environment variable once, at import time.
# ---------------------------------------------------------------------------
CACHE_ROOT = Path('/tmp/pseudospace_mouse_only_cache')
for sub in ['matplotlib', 'numba', 'xdg']:
    (CACHE_ROOT / sub).mkdir(parents=True, exist_ok=True)
os.environ.setdefault('NUMBA_CACHE_DIR', str(CACHE_ROOT / 'numba'))
os.environ.setdefault('MPLCONFIGDIR', str(CACHE_ROOT / 'matplotlib'))
os.environ.setdefault('XDG_CACHE_HOME', str(CACHE_ROOT / 'xdg'))

try:                                    # notebook / interactive-window rendering
    from IPython.display import display, Image
except ImportError:                     # plain `python script.py`
    def display(*objects):
        for obj in objects:
            print(obj)

    def Image(filename=None, **_kwargs):     # noqa: N802 - mirrors the IPython name
        return f'<image: {filename}>'

_INTERACTIVE = hasattr(sys, 'ps1') or 'ipykernel' in sys.modules
if not _INTERACTIVE:
    import matplotlib
    matplotlib.use('Agg')               # no blocking windows during an end-to-end run


# ---------------------------------------------------------------------------
# Paths -- inputs stay under data/, ALL outputs go under RESULTS_DIR.
# PROJECT_DIR is found from this file's location (or the cwd in an interactive
# window), so the script does not care what directory it is launched from.
# ---------------------------------------------------------------------------
def _find_project_dir() -> Path:
    starts = []
    if '__file__' in globals():
        starts.append(Path(__file__).resolve().parent)
    starts.append(Path.cwd().resolve())
    for start in starts:
        for candidate in (start, *start.parents):
            if (candidate / 'pseudospace').is_dir() and (candidate / 'data').is_dir():
                return candidate
    raise RuntimeError(
        'Could not locate the repository root (a directory containing both pseudospace/ and '
        f'data/). Searched upward from {[str(s) for s in starts]}.'
    )


PROJECT_DIR = _find_project_dir()
if str(PROJECT_DIR) not in sys.path:
    sys.path.insert(0, str(PROJECT_DIR))            # makes `import pseudospace` work anywhere

def _workflow_roots(project_dir: Path) -> tuple[Path, Path]:
    """Resolve private inputs and generated outputs without baking in a machine path.

    Command-line values take precedence over environment variables. ``parse_known_args`` keeps
    the cell-marked script usable in Jupyter, whose kernel adds arguments of its own.
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--data-root', type=Path)
    parser.add_argument('--results-root', type=Path)
    args, _ = parser.parse_known_args()
    data_root = args.data_root or os.environ.get('PSEUDOSPACE_DATA_ROOT') or project_dir / 'data'
    results_root = (args.results_root or os.environ.get('PSEUDOSPACE_RESULTS_ROOT')
                    or project_dir / 'results')
    return Path(data_root).expanduser().resolve(), Path(results_root).expanduser().resolve()


DATA_DIR, RESULTS_ROOT = _workflow_roots(PROJECT_DIR)

TUBULE_BY_GENE_DIR = DATA_DIR / 'tubule_by_gene'                          # input
PATHWAY_LIBRARY_DIR = DATA_DIR / 'mouse_vs_human' / 'pathway_gene_sets'   # input (read-only)

RESULTS_DIR = RESULTS_ROOT / 'mouse_only_v5'
HARMONY_OUTPUT_PATH = RESULTS_DIR / 'all_mouse_tubules_harmony_pass1.h5ad'   # pass-1, all cells
PASS2_HARMONY_OUTPUT_PATH = RESULTS_DIR / 'all_mouse_tubules_harmony.h5ad'   # pass-2, tubule-only
DPT_OUTPUT_PATH = RESULTS_DIR / 'all_mouse_tubules_scanpy_dpt.h5ad'
CELLTYPING_DIR = RESULTS_DIR / 'celltyping'
HEATMAP_OUTPUT_DIR = RESULTS_DIR / 'heatmaps'
HEALTHY_VS_AKI_OUTPUT_DIR = RESULTS_DIR / 'healthy_vs_aki'
CONCORDANCE_DIR = RESULTS_DIR / 'concordance'

for directory in [RESULTS_DIR, CELLTYPING_DIR, HEATMAP_OUTPUT_DIR,
                  HEALTHY_VS_AKI_OUTPUT_DIR, CONCORDANCE_DIR]:
    directory.mkdir(parents=True, exist_ok=True)

# Do not inherit a personal R package library; use the active conda environment's pinned Harmony.
os.environ.pop('R_LIBS_USER', None)
# NOTE: the matplotlib/numba/XDG cache redirects live at the top of this file, above the first
# matplotlib import -- they are no-ops if set here.

# ---------------------------------------------------------------------------
# Sample scope
# ---------------------------------------------------------------------------
MOUSE_SAMPLES = ['Ctrl1A2', 'Ctrl1A4', 'IR2A2', 'IR2A4']
BATCH_KEY = 'sample'
RANDOM_STATE = 0

# ---------------------------------------------------------------------------
# Harmony / HVG params. BOTH Harmony passes use INTERSECTION HVGs
# (condition-invariant features; avoids feeding condition-associated genes to Harmony).
# ---------------------------------------------------------------------------
HARMONY_HVG_MIN_MEAN = 0.0125
HARMONY_HVG_MAX_MEAN = 3
HARMONY_HVG_MIN_DISP = 0.5

HARMONY_N_PCS = 50
HARMONY_THETA = 6                    # used for BOTH passes
HARMONY_LAMBDA = 1
HARMONY_MAX_ITER = 30
HARMONY_TAU = 0
HARMONY_PCA_N_COMPS = 50
HARMONY_NEIGHBORS_N = 25
HARMONY_UMAP_MIN_DIST = 0.3
HARMONY_UMAP_SPREAD = 1.0

# ---------------------------------------------------------------------------
# Gene filter / clustering / DPT params
# ---------------------------------------------------------------------------
MIN_GENES_PER_TUBULE = 100            # tubule filter: minimum detected genes per structure
MIN_GENE_TUBULE_FRACTION = 0.05
MIN_GENE_TOTAL_COUNTS = 20
NORMALIZE_TARGET_SUM = 1e4
N_HVGS = 2000

N_NEIGHBORS = 30
N_PCS_FALLBACK = 30
N_DIFFMAP_COMPONENTS_TO_TEST = 8
MIN_VALID_FOR_SPEARMAN = 100
EIGENVALUE_FLOOR = 0.5

COARSE_RESOLUTION = 0.7              # the ONLY clustering run; labelled at the manual checkpoint

# ---------------------------------------------------------------------------
# Topology check (Section 2f.1). DPT assumes the population is a genuine continuum; nothing in
# DPT tests that. PAGA connectivity between anatomically adjacent segments is the test: a
# consecutive pair below this threshold means DPT is INTERPOLATING ACROSS A GAP there, not
# measuring a transition. The nephron is a real risk case -- PT -> thin limb is a sharp
# epithelial boundary, and the collecting duct is a different developmental lineage
# (ureteric bud) from the nephron proper (metanephric mesenchyme).
# ---------------------------------------------------------------------------
PAGA_CONNECTIVITY_THRESHOLD = 0.05

# ---------------------------------------------------------------------------
# Integration QC (Section 2d.1). Batch-mixing metrics are trivially maximised by over-correction,
# so batch ASW is only ever read NEXT TO bio ASW. theta=6 is well above Harmony's default of 2;
# these two numbers, uncorrected vs corrected, are what justify that choice.
# ---------------------------------------------------------------------------
INTEGRATION_QC_MAX_CELLS = 15000     # silhouette is O(n^2); subsample above this


# ---------------------------------------------------------------------------
# Heatmap params
# ---------------------------------------------------------------------------
TRIM_FRACTION = 0.05
N_BINS = 120
SMOOTH_SIGMA = 2.5
Z_CLIP = 2


# %%
TUBULE_BY_GENE_DIR

# %% [markdown]
# ## Section 0.1 - nephron axis markers and the label vocabulary
#
# `NEPHRON_AXIS_MARKERS` is the cell-typing panel: it defines the dotplot, the fine
# `segment_class` vocabulary, the coarse rollup, and the early->late trajectory axis.

# %%
NEPHRON_AXIS_MARKERS = {
    'Podocyte':     ['Nphs1', 'Nphs2', 'Podxl'],
    'PT-S1':        ['Slc5a2', 'Slc5a12', 'Gatm'],
    'PT-S2':        ['Slc22a6', 'Slc13a3', 'Cyp2e1'],
    'PT-S3':        ['Slc7a13', 'Slc22a7', 'Cyp7b1'],   # + Agt/Slc5a1 if desired
    'DTL1':         ['Corin', 'Uncx', 'Slc14a2'],
    'DTL2':         ['Fst', 'Cdh13', 'Stk32a'],
    'DTL3':         ['Nr2e3', 'Dmkn'],
    'ATL':          ['Clcnka', 'Sptssb', 'Akr1b3'],
    'mTAL':         ['Slc12a1', 'Umod', 'Cldn10', 'Ptger3'],   # medullary-enriched
    'cTAL':         ['Slc12a1', 'Umod', 'Cldn16', 'Kcnj10'],   # cortical-enriched
    'Macula-densa': ['Nos1', 'Ptgs2', 'Pappa2'],               # drop Slc12a1 as an ID gene
    'DCT1':         ['Slc12a3', 'Pvalb', 'Trpm6', 'Egf'],      # Mg cassette
    'DCT2':         ['Slc12a3', 'Trpv5', 'Slc8a1', 'S100g'],   # Ca cassette, no Aqp2
    'CNT':          ['Calb1', 'Slc8a1', 'Aqp2', 'Rhcg'],       # Aqp2+ but Slc12a3-
    'CCD':          ['Aqp2', 'Aqp3', 'Kit'],                   # cortical CD
    'OMCD':         ['Aqp2', 'Atp6v0d2', 'Rhcg'],              # (see caveat)
    'IMCD':         ['Aqp2', 'Aqp4', 'Slc14a2', 'Wnt7b'],      # Aqp4/UT-A1 terminal
}
# Anatomical order, as written above: glomerulus -> PT -> thin limb -> TAL -> DCT -> CNT/CD.
NEPHRON_SEGMENT_ORDER = list(NEPHRON_AXIS_MARKERS)

# Genes DELIBERATELY shared by more than one segment above. Unlike the six curated heatmap
# panels (validated further down by _validate_no_repeated_assignments, which forbids repeats),
# this panel is allowed to repeat a gene -- some markers really do span adjacent segments. The
# declaration below is what makes that intentional rather than silent: an undeclared repeat is
# a typo and raises, and the early/late anchor check right after TOTAL_POSITION_MARKERS reads
# this to report which anchors are not segment-exclusive.
NEPHRON_AXIS_MARKER_OVERLAPS = {
    'Slc12a1': ({'mTAL', 'cTAL'}, 'NKCC2 spans the whole TAL; mTAL/cTAL split is Cldn10 vs Cldn16'),
    'Umod':    ({'mTAL', 'cTAL'}, 'uromodulin is pan-TAL'),
    'Slc12a3': ({'DCT1', 'DCT2'}, 'NCC is pan-DCT; DCT1/DCT2 split is the Mg vs Ca cassette'),
    'Slc8a1':  ({'DCT2', 'CNT'}, 'NCX1 spans the Ca cassette shared by DCT2 and CNT'),
    'Aqp2':    ({'CNT', 'CCD', 'OMCD', 'IMCD'}, 'defines the principal-cell lineage, not a region'),
    'Rhcg':    ({'CNT', 'OMCD'}, 'ammonium transport in both the connecting segment and the OMCD'),
    # The one overlap that spans DISTANT segments, and the reason this declaration exists.
    'Slc14a2': ({'DTL1', 'IMCD'}, 'UT-A2 in the descending thin limb vs UT-A1 in the terminal '
                                  'IMCD -- same gene, opposite ends of the nephron'),
}

_segment_gene_map = {}
for _seg, _genes in NEPHRON_AXIS_MARKERS.items():
    for _g in _genes:
        _segment_gene_map.setdefault(_g, set()).add(_seg)
_undeclared = {g: segs for g, segs in _segment_gene_map.items()
               if len(segs) > 1 and segs != NEPHRON_AXIS_MARKER_OVERLAPS.get(g, (None,))[0]}
if _undeclared:
    raise ValueError(
        'NEPHRON_AXIS_MARKERS contains cross-segment gene assignments that are not declared in '
        f'NEPHRON_AXIS_MARKER_OVERLAPS: { {g: sorted(s) for g, s in _undeclared.items()} }. '
        'Either fix the typo, or add the gene to NEPHRON_AXIS_MARKER_OVERLAPS with a reason.'
    )

# Non-epithelial / non-nephron panels. There is no 'Glomerulus' entry here -- the Podocyte row
# of NEPHRON_AXIS_MARKERS covers it and rolls up to the 'Glomerulus' class below.
NON_TUBULE_MARKERS = {
    'Vessel':        ['Pecam1', 'Cdh5', 'Kdr', 'Emcn', 'Flt1', 'Egfl7'],
    'Stroma':        ['Col1a1', 'Col3a1', 'Col1a2', 'Dcn', 'Pdgfrb', 'Mgp'],
    'SmoothMuscle':  ['Acta2', 'Myh11', 'Tagln', 'Rgs5', 'Pdgfra'],
    'Immune':        ['Ptprc', 'C1qa', 'C1qb', 'Cd52', 'Lyz2', 'Cd74'],
}

# Fine segment -> coarse family, the coarse family list, the class vocabulary and the display
# order all come from the shared `pseudospace.vocabulary` module, so notebooks 02, 03 and 04
# cannot drift apart on any of them. Only the marker panels above stay local: they are
# cell-typing evidence, not vocabulary.
from pseudospace import vocabulary as segment_vocabulary

if NEPHRON_SEGMENT_ORDER != list(segment_vocabulary.FINE_SEGMENTS):
    raise ValueError(
        'NEPHRON_AXIS_MARKERS keys no longer match pseudospace.vocabulary.FINE_SEGMENTS: '
        f'{NEPHRON_SEGMENT_ORDER} vs {list(segment_vocabulary.FINE_SEGMENTS)}'
    )

# 'AL' = ascending limb: ATL, mTAL, cTAL and the macula densa are not separable at this
# resolution, so they share one class.
SEGMENT_TO_COARSE = dict(segment_vocabulary.SEGMENT_TO_COARSE)
COARSE_ORDER = list(segment_vocabulary.COARSE_FAMILIES)   # the nephron continuum, in order
KEEP_TUBULE_CLASSES = list(segment_vocabulary.KEEP_TUBULE_CLASSES)
REMOVE_CLASSES = list(segment_vocabulary.REMOVE_CLASSES)

# A coarse family is a LABEL IN ITS OWN RIGHT, not only a rollup target. At the resolutions this
# pipeline actually clusters at most clusters carry coarse identity only: cluster markers like
# Pck1/Slc34a1/Slc4a4 say "proximal tubule" without distinguishing S1 from S2 from S3. Forcing a
# fine label onto such a cluster invents a distinction the data does not support, and it distorts
# the pseudospace axis those labels annotate. So each coarse family maps to ITSELF and
# 'Glomerulus' is directly typable alongside 'Podocyte': mixing 'PT' for one cluster and 'DCT2'
# for another is expected, not a compromise. Both mappings live in the module.
#
# LABEL_VOCABULARY is every label you may type at the manual checkpoint. SEGMENT_DISPLAY_ORDER
# places each coarse family immediately before its own first fine member so the anatomical
# left-to-right reading survives: ... PT, PT-S1, PT-S2, PT-S3, DTL, DTL1, ...
LABEL_VOCABULARY = list(segment_vocabulary.LABEL_VOCABULARY)
SEGMENT_DISPLAY_ORDER = list(segment_vocabulary.SEGMENT_DISPLAY_ORDER)


# ---------------------------------------------------------------------------
# Injury panels (Section 4.0b). NOT used for cell typing or for building the pseudospace axis --
# these exist purely to test whether a Section 4 "shape" hit is confounded by injury state.
#
# Why this matters: in mouse IRI, proximal tubule cells DE-DIFFERENTIATE and lose the very
# segment-identity markers (Slc5a2, Slc34a1, Slc22a6, ...) that define the pseudospace axis.
# So an injured S2 tubule can move along the axis because it is injured, not because it sits
# somewhere else in the nephron. Worse, failed-repair PT forms a BRANCH off the repair
# trajectory (Kirita 2021), and DPT runs with n_branchings=0, so those tubules get projected
# onto the linear axis at an essentially arbitrary position.
#   Kirita et al. 2020 PNAS 117(27):15874  -- cell profiling of mouse AKI
#   Kirita et al. 2021 PNAS 118(3):e2026684118 -- PT cell states, FR-PTC branch
#
# Two panels, deliberately kept separate: they peak at different times post-IRI, so collapsing
# them into one score would blur an acute (day ~2) response against a maladaptive (week+) one.
FR_PTC_MARKERS = ['Vcam1', 'Dcdc2a', 'Krt20', 'Sox9']            # failed-repair / maladaptive
ACUTE_INJURY_MARKERS = ['Havcr1', 'Lcn2', 'Spp1']                # Kim1 / NGAL / osteopontin

# Early->late anchors for the marker axis that orients DPT (start of PT -> end of the CD).
AXIS_EARLY_SEGMENTS = ('PT-S1', 'PT-S2')
AXIS_LATE_SEGMENTS = ('OMCD', 'IMCD')
TOTAL_POSITION_MARKERS = {
    'early': sorted({g for seg in AXIS_EARLY_SEGMENTS for g in NEPHRON_AXIS_MARKERS[seg]}),
    'late': sorted({g for seg in AXIS_LATE_SEGMENTS for g in NEPHRON_AXIS_MARKERS[seg]}),
}

# Anchor purity. An anchor gene that also marks a segment OUTSIDE its own anchor block pulls
# tubules from that segment toward the wrong end of the axis. Sharing with a NEIGHBOURING
# segment is harmless (Aqp2 in CNT/CCD sits right next to the OMCD/IMCD late block); sharing
# across the nephron is not (Slc14a2 marks DTL1, ~12 segments upstream of the IMCD, so every
# DTL1 tubule scores on the late panel). The distance is reported so the difference is visible
# rather than buried in the panel definition.
#
# Set STRICT_AXIS_ANCHORS = True to DROP long-range contaminants from the anchor panels. This
# changes total_marker_axis, and therefore the DPT root, the orientation and every downstream
# number -- it is off by default so the committed results stay reproducible.
STRICT_AXIS_ANCHORS = False
AXIS_ANCHOR_MAX_SEGMENT_DISTANCE = 3      # segments; beyond this an overlap counts as long-range

_anchor_blocks = {'early': AXIS_EARLY_SEGMENTS, 'late': AXIS_LATE_SEGMENTS}
axis_anchor_contamination = []
for _role, _block in _anchor_blocks.items():
    _block_pos = [NEPHRON_SEGMENT_ORDER.index(s) for s in _block]
    for _gene in TOTAL_POSITION_MARKERS[_role]:
        for _other in sorted(_segment_gene_map[_gene] - set(_block)):
            _dist = min(abs(NEPHRON_SEGMENT_ORDER.index(_other) - p) for p in _block_pos)
            axis_anchor_contamination.append(
                {'role': _role, 'gene': _gene, 'also_marks': _other, 'segment_distance': _dist,
                 'long_range': _dist > AXIS_ANCHOR_MAX_SEGMENT_DISTANCE})

_long_range = {r['gene'] for r in axis_anchor_contamination if r['long_range']}
if _long_range:
    print(f'AXIS ANCHOR WARNING: {sorted(_long_range)} appear in the early/late anchor panels but '
          f'also mark segments more than {AXIS_ANCHOR_MAX_SEGMENT_DISTANCE} positions away:')
    for _r in axis_anchor_contamination:
        if _r['long_range']:
            print(f"    {_r['gene']:>9s} ({_r['role']:>5s} anchor) also marks {_r['also_marks']} "
                  f"-- {_r['segment_distance']} segments away"
                  f"  [{NEPHRON_AXIS_MARKER_OVERLAPS.get(_r['gene'], (None, 'undeclared'))[1]}]")
    if STRICT_AXIS_ANCHORS:
        for _role in TOTAL_POSITION_MARKERS:
            TOTAL_POSITION_MARKERS[_role] = [g for g in TOTAL_POSITION_MARKERS[_role]
                                             if g not in _long_range]
        print('  STRICT_AXIS_ANCHORS=True -> dropped from the anchor panels.')
    else:
        print('  STRICT_AXIS_ANCHORS=False -> KEPT. Tubules of those segments are pulled toward '
              'the wrong end of the marker axis; that axis roots and orients DPT, so read the '
              'PAGA check (2f.1) and the segment-ordering plot (2f) with this in mind.')

# Strip / annotation colors for both the fine and the coarse vocabulary.
SEGMENT_STRIP_COLORS = {
    'Podocyte': '#9E9E9E',
    'PT-S1': '#4C9BD3', 'PT-S2': '#5B8E55', 'PT-S3': '#D6B48A',
    'DTL1': '#B39DDB', 'DTL2': '#9575CD', 'DTL3': '#7E57C2',
    'ATL': '#26A69A', 'mTAL': '#F58518', 'cTAL': '#E8A33D', 'Macula-densa': '#C7522A',
    'DCT1': '#D81B60', 'DCT2': '#EC6EA5',
    'CNT': '#76B7B2', 'CCD': '#4E79A7', 'OMCD': '#9C755F', 'IMCD': '#593C8F',
    'PT': '#4C9BD3', 'DTL': '#7E57C2', 'AL': '#F58518', 'DCT': '#D81B60', 'CNT_CD': '#8D6E63',
    'Glomerulus': '#9E9E9E',
}

print(f'Nephron segments: {len(NEPHRON_SEGMENT_ORDER)}  ->  coarse classes: {COARSE_ORDER}')
print(f"Marker axis early: {TOTAL_POSITION_MARKERS['early']}")
print(f"Marker axis late:  {TOTAL_POSITION_MARKERS['late']}")

# %% [markdown]
# ## Section 0.2 - curated heatmap panels
#
# Carried over verbatim from the notebook: chart- and paper-derived, with no gene assigned to
# more than one module within a panel. These drive the Section 3 heatmap rows only; cell typing
# uses NEPHRON_AXIS_MARKERS above.

# %%
# ==============================================================================
# PROXIMAL TUBULE POSITIONAL MARKERS
# ==============================================================================
# Assignment:
#   S1  => Slc5a2/Slc5a12 (S1-dominant glucose reabsorption)
#   S2  => Slc22a6/Slc13a3 (S2-enriched organic-anion + dicarboxylate transport)
#   S3  => Slc7a13/Slc22a7 (canonical S3 apical markers)
# NOTE: Slc7a13, Slc22a7, Cyp7b1, Slc7a12 are SEX-DIMORPHIC in S3
#       (male- vs female-restricted). If reps mix sexes these behave
#       inconsistently and should be down-weighted as identity genes.
#       (Ransick et al., Dev Cell 2019)
PT_MARKER_GROUPS = {
    'PT-S1': [
        {'label': 'SLC5A2 / Slc5a2',  'candidates': ['SLC5A2', 'Slc5a2', 'slc5a2']},
        {'label': 'SLC5A12 / Slc5a12', 'candidates': ['SLC5A12', 'Slc5a12', 'slc5a12']},
        {'label': 'GATM / Gatm',       'candidates': ['GATM', 'Gatm', 'gatm']},
    ],
    'PT-S2': [
        {'label': 'SLC22A6 / Slc22a6', 'candidates': ['SLC22A6', 'Slc22a6', 'slc22a6']},
        {'label': 'SLC13A3 / Slc13a3', 'candidates': ['SLC13A3', 'Slc13a3', 'slc13a3']},
        {'label': 'SLC1A1 / Slc1a1',   'candidates': ['SLC1A1', 'Slc1a1', 'slc1a1']},
        {'label': 'CYP2E1 / Cyp2e1',   'candidates': ['CYP2E1', 'Cyp2e1', 'cyp2e1']},
    ],
    'PT-S3': [
        {'label': 'SLC7A13 / Slc7a13', 'candidates': ['SLC7A13', 'Slc7a13', 'slc7a13']},  # moved from S2
        {'label': 'SLC22A7 / Slc22a7', 'candidates': ['SLC22A7', 'Slc22a7', 'slc22a7']},
        {'label': 'SLC5A1 / Slc5a1',   'candidates': ['SLC5A1', 'Slc5a1', 'slc5a1']},
        {'label': 'CYP7B1 / Cyp7b1',   'candidates': ['CYP7B1', 'Cyp7b1', 'cyp7b1']},     # male-restricted
    ],
}

PT_GROUP_ORDER = ['PT-S1', 'PT-S2', 'PT-S3']

PT_GROUP_COLORS = {
    'PT-S1': '#4C9BD3',
    'PT-S2': '#5B8E55',
    'PT-S3': '#D6B48A',
}

# PT identity markers that do not resolve S1/S2/S3.
# S1/S2 is a continuum; expect adjacent blending there regardless of panel.
PT_PAN_MARKERS = [
    {'label': 'SLC34A1 / Slc34a1', 'candidates': ['SLC34A1', 'Slc34a1', 'slc34a1']},  # all PT except medullary S3
    {'label': 'LRP2 / Lrp2',       'candidates': ['LRP2', 'Lrp2', 'lrp2']},
    {'label': 'SLC4A4 / Slc4a4',   'candidates': ['SLC4A4', 'Slc4a4', 'slc4a4']},
    {'label': 'SLC3A1 / Slc3a1',   'candidates': ['SLC3A1', 'Slc3a1', 'slc3a1']},
    {'label': 'SLC22A5 / Slc22a5', 'candidates': ['SLC22A5', 'Slc22a5', 'slc22a5']},
]


# ==============================================================================
# THIN-LIMB SUBSEGMENTS — CHEN ET AL., JASN 2021
# ==============================================================================
# DTL3 (Nr2e3/Dmkn) are lower-confidence, data-driven markers — validate in
# your own data before trusting them as identity genes.
THIN_LIMB_MARKER_GROUPS = {
    'DTL1': [
        {'label': 'CORIN / Corin',   'candidates': ['CORIN', 'Corin', 'corin']},
        {'label': 'UNCX / Uncx',     'candidates': ['UNCX', 'Uncx', 'uncx']},
        {'label': 'SLC14A2 / Slc14a2', 'candidates': ['SLC14A2', 'Slc14a2', 'slc14a2']},  # urea channel, DTL1
    ],
    'DTL2': [
        {'label': 'FST / Fst',       'candidates': ['FST', 'Fst', 'fst']},
        {'label': 'CDH13 / Cdh13',   'candidates': ['CDH13', 'Cdh13', 'cdh13']},
        {'label': 'STK32A / Stk32a', 'candidates': ['STK32A', 'Stk32a', 'stk32a']},
    ],
    'DTL3': [
        {'label': 'NR2E3 / Nr2e3',   'candidates': ['NR2E3', 'Nr2e3', 'nr2e3']},
        {'label': 'DMKN / Dmkn',     'candidates': ['DMKN', 'Dmkn', 'dmkn']},
    ],
    'ATL': [
        {'label': 'CLCNKA / Clcnka', 'candidates': ['CLCNKA', 'Clcnka', 'clcnka']},
        {'label': 'SPTSSB / Sptssb', 'candidates': ['SPTSSB', 'Sptssb', 'sptssb']},
        {'label': 'IGFBP2 / Igfbp2', 'candidates': ['IGFBP2', 'Igfbp2', 'igfbp2']},
        {'label': 'CCDC3 / Ccdc3',   'candidates': ['CCDC3', 'Ccdc3', 'ccdc3']},
    ],
}

THIN_LIMB_GROUP_ORDER = ['DTL1', 'DTL2', 'DTL3', 'ATL']

THIN_LIMB_GROUP_COLORS = {
    'DTL1': '#8AC926',
    'DTL2': '#55A630',
    'DTL3': '#2B9348',
    'ATL':  '#007F5F',
}


# ==============================================================================
# MEDULLARY VERSUS CORTICAL TAL  +  MACULA DENSA
# ==============================================================================
# CAVEAT: mTAL vs cTAL is a MOSAIC/GRADIENT, not on/off. No gene cleanly
#   separates them. Cldn10 spans the whole TAL but is medullary-enriched;
#   Cldn16/Kcnj10/Casr/Pth1r are cortex + outer-stripe enriched.
#   Sycn/Shd retained as experimentally-validated regional markers.
#   (Chevalier/Lee JCI Insight 2024; Milatz PMC11011785)
# Macula densa is a cTAL sub-population; give it MD-restricted identity genes
#   (Nos1/Ptgs2/Pappa2) — do NOT use Slc12a1 as an MD ID gene.
TAL_MARKER_GROUPS = {
    'mTAL': [
        {'label': 'CLDN10 / Cldn10',   'candidates': ['CLDN10', 'Cldn10', 'cldn10']},   # medullary-enriched
        {'label': 'SYCN / Sycn',       'candidates': ['SYCN', 'Sycn', 'sycn']},
        {'label': 'TMEM86A / Tmem86a', 'candidates': ['TMEM86A', 'Tmem86a', 'tmem86a']},
    ],
    'cTAL': [
        {'label': 'CLDN16 / Cldn16',   'candidates': ['CLDN16', 'Cldn16', 'cldn16']},   # cortical-enriched
        {'label': 'KCNJ10 / Kcnj10',   'candidates': ['KCNJ10', 'Kcnj10', 'kcnj10']},
        {'label': 'SHD / Shd',         'candidates': ['SHD', 'Shd', 'shd']},
    ],
    'Macula-densa': [
        {'label': 'NOS1 / Nos1',       'candidates': ['NOS1', 'Nos1', 'nos1']},
        {'label': 'PTGS2 / Ptgs2',     'candidates': ['PTGS2', 'Ptgs2', 'ptgs2', 'COX2', 'Cox2']},
        {'label': 'PAPPA2 / Pappa2',   'candidates': ['PAPPA2', 'Pappa2', 'pappa2']},
    ],
}

TAL_GROUP_ORDER = ['mTAL', 'cTAL', 'Macula-densa']

TAL_GROUP_COLORS = {
    'mTAL': '#4C78A8',
    'cTAL': '#F58518',
    'Macula-densa': '#E45756',
}

# TAL identity markers that do not separate mTAL/cTAL/MD.
TAL_PAN_MARKERS = [
    {'label': 'SLC12A1 / Slc12a1', 'candidates': ['SLC12A1', 'Slc12a1', 'slc12a1']},
    {'label': 'UMOD / Umod',       'candidates': ['UMOD', 'Umod', 'umod']},
]


# ==============================================================================
# DCT1 VERSUS DCT2  (Mg cassette vs Ca cassette)
# ==============================================================================
# Both first satisfy pan-DCT Slc12a3.
#   DCT1: Mg cassette (Pvalb, Trpm6, Egf, Erbb4)
#   DCT2: Ca cassette (Trpv5, Slc8a1, Calb1, S100g), Pvalb low/absent
# The Ca cassette is SHARED with CNT — the only clean DCT2-vs-CNT separators
# are Slc12a3 (DCT2+, CNT-) and Aqp2 (CNT+, DCT2-).
# (Grimm/McDonough; Poulsen PMC11000743)
DCT_MARKER_GROUPS = {
    'DCT1': [
        {'label': 'PVALB / Pvalb',   'candidates': ['PVALB', 'Pvalb', 'pvalb']},
        {'label': 'TRPM6 / Trpm6',   'candidates': ['TRPM6', 'Trpm6', 'trpm6']},
        {'label': 'EGF / Egf',       'candidates': ['EGF', 'Egf', 'egf']},
        {'label': 'ERBB4 / Erbb4',   'candidates': ['ERBB4', 'Erbb4', 'erbb4']},
    ],
    'DCT2': [
        {'label': 'TRPV5 / Trpv5',   'candidates': ['TRPV5', 'Trpv5', 'trpv5']},
        {'label': 'SLC8A1 / Slc8a1', 'candidates': ['SLC8A1', 'Slc8a1', 'slc8a1']},
    ],
}

DCT_GROUP_ORDER = ['DCT1', 'DCT2']

DCT_GROUP_COLORS = {
    'DCT1': '#4C78A8',
    'DCT2': '#F58518',
}

DCT_PAN_MARKERS = [
    {'label': 'SLC12A3 / Slc12a3', 'candidates': ['SLC12A3', 'Slc12a3', 'slc12a3']},
]

# Ca-cassette genes that also occur in CNT — supporting evidence only,
# never scored as uniquely DCT2-positive.
DCT2_SUPPORT_MARKERS = [
    {'label': 'CALB1 / Calb1',   'candidates': ['CALB1', 'Calb1', 'calb1']},
    {'label': 'S100G / S100g',   'candidates': ['S100G', 'S100g', 's100g']},
    {'label': 'VDR / Vdr',       'candidates': ['VDR', 'Vdr', 'vdr']},
]

DCT_CLASSIFICATION_RULES = {
    'DCT1': {
        'required_positive': ['SLC12A3', 'PVALB'],
        'support_positive': ['TRPM6', 'EGF', 'ERBB4'],
        'expected_negative': ['AQP2'],
    },
    'DCT2': {
        'required_positive': ['SLC12A3', 'TRPV5'],
        'support_positive': ['SLC8A1', 'CALB1', 'S100G'],
        'expected_negative': ['PVALB', 'AQP2'],  # AQP2- separates DCT2 from CNT
    },
}


# ==============================================================================
# CNT AND COLLECTING-DUCT POSITIONAL REGIONS
# ==============================================================================
# CAVEAT: CCD vs OMCD has essentially NO unique transcript — same principal
#   cell at different cortico-medullary depth. Do not expect a clean split;
#   treat CCD/OMCD confusion as expected. IMCD is the only CD region with a
#   real discriminator (Aqp4 basolateral + Slc14a2/UT-A1 terminal).
# CNT: Aqp2+ but Slc12a3- ; Rhcg adds CNT specificity. Calb1/S100g are shared
#   with DCT2 (support only).
# (Chen JASN 2021; Nielsen AJP-Renal 2005 for Aqp3/Aqp4 gradient)
COLLECTING_REGION_MARKER_GROUPS = {
    'CNT': [
        {'label': 'RHCG / Rhcg',   'candidates': ['RHCG', 'Rhcg', 'rhcg']},
        {'label': 'SLC8A1 / Slc8a1', 'candidates': ['SLC8A1', 'Slc8a1', 'slc8a1']},  # NCX1, connecting segment
    ],
    'CCD': [
        {'label': 'AQP3 / Aqp3',   'candidates': ['AQP3', 'Aqp3', 'aqp3']},          # CD-wide, higher cortical
    ],
    'OMCD': [
        {'label': 'ATP6V0D2 / Atp6v0d2', 'candidates': ['ATP6V0D2', 'Atp6v0d2', 'atp6v0d2']},  # IC-heavy region; weak PC-region ID
    ],
    'IMCD': [
        {'label': 'AQP4 / Aqp4',     'candidates': ['AQP4', 'Aqp4', 'aqp4']},        # IMCD-restricted basolateral
        {'label': 'SLC14A2 / Slc14a2', 'candidates': ['SLC14A2', 'Slc14a2', 'slc14a2']},  # UT-A1, terminal IMCD
        {'label': 'WNT7B / Wnt7b',   'candidates': ['WNT7B', 'Wnt7b', 'wnt7b']},
    ],
}

COLLECTING_REGION_GROUP_ORDER = ['CNT', 'CCD', 'OMCD', 'IMCD']

COLLECTING_REGION_GROUP_COLORS = {
    'CNT':  '#76B7B2',
    'CCD':  '#4E79A7',
    'OMCD': '#9C755F',
    'IMCD': '#593C8F',
}

# Aqp2 spans CNT->IMCD and defines the PC lineage, not a single region.
COLLECTING_REGION_PAN_MARKERS = [
    {'label': 'AQP2 / Aqp2', 'candidates': ['AQP2', 'Aqp2', 'aqp2']},
]


# ==============================================================================
# COLLECTING-DUCT CELL IDENTITIES  (orthogonal to position)
# ==============================================================================
# Interspersed A-IC / B-IC will otherwise be misassigned into whichever CD
# region they sit in. Gate ICs out (or classify them) BEFORE the positional
# CD call above. (Chen PNAS 2017; Nanami Nat Commun 2021)
COLLECTING_CELLTYPE_MARKER_GROUPS = {
    'CD-PC': [
        {'label': 'AQP2 / Aqp2',     'candidates': ['AQP2', 'Aqp2', 'aqp2']},
        {'label': 'AQP3 / Aqp3',     'candidates': ['AQP3', 'Aqp3', 'aqp3']},
        {'label': 'SCNN1G / Scnn1g', 'candidates': ['SCNN1G', 'Scnn1g', 'scnn1g']},
        {'label': 'FXYD4 / Fxyd4',   'candidates': ['FXYD4', 'Fxyd4', 'fxyd4']},
    ],
    'CD-A-IC': [
        {'label': 'SLC4A1 / Slc4a1',   'candidates': ['SLC4A1', 'Slc4a1', 'slc4a1']},
        {'label': 'AQP6 / Aqp6',       'candidates': ['AQP6', 'Aqp6', 'aqp6']},
        {'label': 'SLC26A7 / Slc26a7', 'candidates': ['SLC26A7', 'Slc26a7', 'slc26a7']},
        {'label': 'KIT / Kit',         'candidates': ['KIT', 'Kit', 'kit']},
    ],
    'CD-B-IC': [
        {'label': 'SLC26A4 / Slc26a4', 'candidates': ['SLC26A4', 'Slc26a4', 'slc26a4']},
        {'label': 'SLC4A9 / Slc4a9',   'candidates': ['SLC4A9', 'Slc4a9', 'slc4a9']},
        {'label': 'INSRR / Insrr',     'candidates': ['INSRR', 'Insrr', 'insrr']},
    ],
}

COLLECTING_CELLTYPE_GROUP_ORDER = ['CD-PC', 'CD-A-IC', 'CD-B-IC']

COLLECTING_CELLTYPE_GROUP_COLORS = {
    'CD-PC':   '#4E79A7',
    'CD-A-IC': '#E15759',
    'CD-B-IC': '#B07AA1',
}

# Pan-IC (identify IC lineage, do not split A vs B).
IC_PAN_MARKERS = [
    {'label': 'FOXI1 / Foxi1',       'candidates': ['FOXI1', 'Foxi1', 'foxi1']},
    {'label': 'ATP6V1B1 / Atp6v1b1', 'candidates': ['ATP6V1B1', 'Atp6v1b1', 'atp6v1b1']},
    {'label': 'ATP6V1G3 / Atp6v1g3', 'candidates': ['ATP6V1G3', 'Atp6v1g3', 'atp6v1g3']},
    {'label': 'ATP6V0D2 / Atp6v0d2', 'candidates': ['ATP6V0D2', 'Atp6v0d2', 'atp6v0d2']},
]


# ==============================================================================
# CONTINUOUS AXIS MARKERS
# ==============================================================================
# One non-repeated representative identity marker per positional segment.
# Contrast-based calls flagged inline.
TOTAL_NEPHRON_MARKERS = [
    {'label': 'Nphs1 (Glom)',           'candidates': ['NPHS1', 'Nphs1', 'nephrin']},
    {'label': 'Nphs2 (Glom)',           'candidates': ['NPHS2', 'Nphs2', 'podocin']},

    {'label': 'Slc5a2 (PT-S1)',         'candidates': ['SLC5A2', 'Slc5a2', 'slc5a2']},
    {'label': 'Slc22a6 (PT-S2)',        'candidates': ['SLC22A6', 'Slc22a6', 'slc22a6']},
    {'label': 'Slc7a13 (PT-S3; sex-dim)', 'candidates': ['SLC7A13', 'Slc7a13', 'slc7a13']},

    {'label': 'Corin (DTL1)',           'candidates': ['CORIN', 'Corin', 'corin']},
    {'label': 'Cdh13 (DTL2)',           'candidates': ['CDH13', 'Cdh13', 'cdh13']},
    {'label': 'Nr2e3 (DTL3; low-conf)', 'candidates': ['NR2E3', 'Nr2e3', 'nr2e3']},
    {'label': 'Sptssb (ATL)',           'candidates': ['SPTSSB', 'Sptssb', 'sptssb']},

    {'label': 'Cldn10 (mTAL-enriched)', 'candidates': ['CLDN10', 'Cldn10', 'cldn10']},
    {'label': 'Cldn16 (cTAL-enriched)', 'candidates': ['CLDN16', 'Cldn16', 'cldn16']},
    {'label': 'Nos1 (Macula densa)',    'candidates': ['NOS1', 'Nos1', 'nos1']},

    {'label': 'Pvalb (DCT1)',           'candidates': ['PVALB', 'Pvalb', 'pvalb']},
    {'label': 'Trpv5 (DCT2; Slc12a3+/Pvalb-)', 'candidates': ['TRPV5', 'Trpv5', 'trpv5']},

    {'label': 'Rhcg (CNT; Aqp2+/Slc12a3-)', 'candidates': ['RHCG', 'Rhcg', 'rhcg']},
    {'label': 'Aqp3 (CCD; ~= OMCD)',    'candidates': ['AQP3', 'Aqp3', 'aqp3']},
    {'label': 'Aqp2 (OMCD; not unique vs CCD)', 'candidates': ['AQP2', 'Aqp2', 'aqp2']},
    {'label': 'Slc14a2 (IMCD; UT-A1)',  'candidates': ['SLC14A2', 'Slc14a2', 'slc14a2']},

    {'label': 'Slc4a1 (A-IC)',          'candidates': ['SLC4A1', 'Slc4a1', 'slc4a1']},
    {'label': 'Slc26a4 (B-IC)',         'candidates': ['SLC26A4', 'Slc26a4', 'slc26a4']},
]


# ==============================================================================
# BOUNDARIES WITH NO UNIQUE MARKER (expected confusion in the matrix)
# ==============================================================================
# These are true cortico-medullary / axial gradients. A confusion matrix that
# blends ONLY these adjacent pairs is the correct result, not a bug.
GRADIENT_BOUNDARIES = [
    ('PT-S1', 'PT-S2'),   # convoluted continuum
    ('PT-S2', 'PT-S3'),
    ('mTAL',  'cTAL'),    # mosaic claudin expression, not on/off
    ('DCT2',  'CNT'),     # shared Ca cassette; split only via Slc12a3/Aqp2
    ('CCD',   'OMCD'),    # same PC, different depth — no unique transcript
]


# ==============================================================================
# PROVENANCE
# ==============================================================================
MARKER_PROVENANCE = {
    'PT': {
        'paper': 'Ransick et al., Single-Cell Profiling Reveals Sex, Lineage, '
                 'and Regional Diversity in the Mouse Kidney',
        'note': 'Slc7a13 assigned to S3 (canonical S3 apical marker); '
                'S3 markers are sex-dimorphic.',
        'url': 'https://www.cell.com/developmental-cell/fulltext/S1534-5807(19)30814-7',
    },
    'DTL1_DTL2_DTL3_ATL': {
        'paper': 'Chen et al., A Comprehensive Map of mRNAs and Their Isoforms '
                 'across All 14 Renal Tubule Segments of Mouse',
        'url': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC8017530/',
    },
    'mTAL_cTAL_MD': {
        'paper': 'Lee et al., Distinct cell types along thick ascending limb '
                 '(JCI Insight 2024); Milatz, Claudin-10 TAL gene expression '
                 '(2024); Sycn/Shd from Atlas of Gene Expression Mouse Kidney',
        'note': 'mTAL/cTAL is a mosaic gradient (Cldn10 vs Cldn16/Kcnj10). '
                'MD uses Nos1/Ptgs2/Pappa2, NOT Slc12a1.',
        'url': 'https://insight.jci.org/articles/view/190992',
    },
    'DCT1_DCT2': {
        'paper': 'Heterogeneity of Distal Convoluted Tubule Cells '
                 '(Mg vs Ca cassette)',
        'note': 'DCT1 = Mg cassette (Pvalb/Trpm6/Egf/Erbb4); '
                'DCT2 = Ca cassette (Trpv5/Slc8a1/Calb1/S100g).',
        'url': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC11000743/',
    },
    'CNT_CCD_OMCD_IMCD': {
        'paper': 'Chen et al., 14 Renal Tubule Segments; '
                 'Nielsen, Aquaporin CD gradient (AJP-Renal 2005)',
        'note': 'Atp4a (gastric H/K-ATPase) removed from CCD; Hsd11b2 removed '
                'from OMCD (spans DCT2->CNT->CD PCs); Aqp4/Slc14a2 restricted '
                'to IMCD; CCD vs OMCD has no unique transcript.',
        'url': 'https://journals.physiology.org/doi/10.1152/ajprenal.00273.2005',
    },
    'CD_PC_AIC_BIC': {
        'paper': 'Chen et al., Transcriptomes of Major Renal Collecting Duct '
                 'Cell Types (PNAS 2017); Nanami, IC subtypes (Nat Commun 2021)',
        'url': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC5699061/',
    },
}


# ==============================================================================
# VALIDATION
# ==============================================================================
def _canonical_gene(marker):
    return marker['candidates'][0].upper()


def _validate_no_repeated_assignments(grouped_markers, panel_name):
    seen = {}

    for group, markers in grouped_markers.items():
        for marker in markers:
            gene = _canonical_gene(marker)

            if gene in seen:
                raise ValueError(
                    f'{gene} is assigned to both {seen[gene]!r} and '
                    f'{group!r} in {panel_name}.'
                )

            seen[gene] = group

    return seen


_validate_no_repeated_assignments(
    PT_MARKER_GROUPS,
    'PT_MARKER_GROUPS',
)
_validate_no_repeated_assignments(
    THIN_LIMB_MARKER_GROUPS,
    'THIN_LIMB_MARKER_GROUPS',
)
_validate_no_repeated_assignments(
    TAL_MARKER_GROUPS,
    'TAL_MARKER_GROUPS',
)
_validate_no_repeated_assignments(
    DCT_MARKER_GROUPS,
    'DCT_MARKER_GROUPS',
)
_validate_no_repeated_assignments(
    COLLECTING_REGION_MARKER_GROUPS,
    'COLLECTING_REGION_MARKER_GROUPS',
)
_validate_no_repeated_assignments(
    COLLECTING_CELLTYPE_MARKER_GROUPS,
    'COLLECTING_CELLTYPE_MARKER_GROUPS',
)

print(
    'Marker panels defined with non-repeated assignments; '
    'pan-segment markers and contrast-based calls are stored separately.'
)

# %%
# Panel -> the coarse class(es) whose cells it describes. THIN_LIMB spans two classes because
# DTL1-3 roll up to 'DTL' while ATL rolls up to 'AL'.
HEATMAP_PANELS = {
    'PT': {'markers': PT_MARKER_GROUPS, 'order': PT_GROUP_ORDER,
           'colors': PT_GROUP_COLORS, 'classes': ['PT'], 'z_clip': Z_CLIP},
    'THIN_LIMB': {'markers': THIN_LIMB_MARKER_GROUPS, 'order': THIN_LIMB_GROUP_ORDER,
                  'colors': THIN_LIMB_GROUP_COLORS, 'classes': ['DTL', 'AL'], 'z_clip': Z_CLIP},
    'TAL': {'markers': TAL_MARKER_GROUPS, 'order': TAL_GROUP_ORDER,
            'colors': TAL_GROUP_COLORS, 'classes': ['AL'], 'z_clip': Z_CLIP},
    'DCT': {'markers': DCT_MARKER_GROUPS, 'order': DCT_GROUP_ORDER,
            'colors': DCT_GROUP_COLORS, 'classes': ['DCT'], 'z_clip': Z_CLIP},
    'CD_REGION': {'markers': COLLECTING_REGION_MARKER_GROUPS, 'order': COLLECTING_REGION_GROUP_ORDER,
                  'colors': COLLECTING_REGION_GROUP_COLORS, 'classes': ['CNT_CD'], 'z_clip': Z_CLIP},
    'CD_CELLTYPE': {'markers': COLLECTING_CELLTYPE_MARKER_GROUPS, 'order': COLLECTING_CELLTYPE_GROUP_ORDER,
                    'colors': COLLECTING_CELLTYPE_GROUP_COLORS, 'classes': ['CNT_CD'], 'z_clip': Z_CLIP},
}

# One recomputed within-family DPT per coarse class, anchored by that family's panel.
SUBSET_DPT_SPECS = {
    'PT': {'panel': 'PT', 'pseudotime_col': 'pt_subset_scanpy_dpt'},
    'DTL': {'panel': 'THIN_LIMB', 'pseudotime_col': 'dtl_subset_scanpy_dpt'},
    'AL': {'panel': 'TAL', 'pseudotime_col': 'al_subset_scanpy_dpt'},
    'DCT': {'panel': 'DCT', 'pseudotime_col': 'dct_subset_scanpy_dpt'},
    'CNT_CD': {'panel': 'CD_REGION', 'pseudotime_col': 'cnt_cd_subset_scanpy_dpt'},
}


def _panel_genes(*panels):
    """Every candidate symbol mentioned by a marker panel (dict-of-lists or list-of-dicts)."""
    genes = set()
    for panel in panels:
        entries = panel.values() if isinstance(panel, dict) else [panel]
        for entry in entries:
            for item in entry:
                if isinstance(item, dict):
                    genes.update(item['candidates'])
                else:
                    genes.add(item)
    return genes


# Curated markers are never gated out by the expression filter in Section 2a.
CURATED_MARKER_WHITELIST = sorted(_panel_genes(
    NEPHRON_AXIS_MARKERS, NON_TUBULE_MARKERS,
    PT_MARKER_GROUPS, THIN_LIMB_MARKER_GROUPS, TAL_MARKER_GROUPS, DCT_MARKER_GROUPS,
    COLLECTING_REGION_MARKER_GROUPS, COLLECTING_CELLTYPE_MARKER_GROUPS,
    PT_PAN_MARKERS, TAL_PAN_MARKERS, DCT_PAN_MARKERS, COLLECTING_REGION_PAN_MARKERS,
    IC_PAN_MARKERS, TOTAL_NEPHRON_MARKERS,
))

# Genes protected from the Section 2a expression filter but NOT added to the clustering feature
# space -- kept separate from CURATED_MARKER_WHITELIST precisely because that one feeds
# `selected_for_clustering`, and adding to it would move the Leiden result and invalidate
# COARSE_LABELS. The injury panels are diagnostic only (Sections 4.1b / 4.6b-d) and must never
# influence cell typing or the pseudospace axis, but they still have to survive the gene filter:
# Vcam1 and Dcdc2a sit around 8% of tubules, not far above MIN_GENE_TUBULE_FRACTION, and losing
# them would silently gut the injury-confounding analyses AND the 4.6d circularity guard.
DIAGNOSTIC_GENE_WHITELIST = sorted(_panel_genes(FR_PTC_MARKERS, ACUTE_INJURY_MARKERS))
GENE_FILTER_WHITELIST = sorted(set(CURATED_MARKER_WHITELIST) | set(DIAGNOSTIC_GENE_WHITELIST))

print(f'Results directory:           {RESULTS_DIR.relative_to(PROJECT_DIR)}')
print(f'Mouse samples:               {MOUSE_SAMPLES}')
print(f'Pass-1 Harmony output:       {HARMONY_OUTPUT_PATH.relative_to(PROJECT_DIR)}')
print(f'Pass-2 Harmony output:       {PASS2_HARMONY_OUTPUT_PATH.relative_to(PROJECT_DIR)}')
print(f'DPT output:                  {DPT_OUTPUT_PATH.relative_to(PROJECT_DIR)}')
print(f'Curated markers whitelisted: {len(CURATED_MARKER_WHITELIST)} '
      f'(clustering features + gene filter)')
print(f'Diagnostic genes whitelisted: {len(DIAGNOSTIC_GENE_WHITELIST)} '
      f'(gene filter only): {DIAGNOSTIC_GENE_WHITELIST}')

# %% [markdown]
# # Section 1 - pass-1 Harmony integration
#
# Native mouse gene symbols throughout: no HCOP table, no ortholog mapping, no human samples.

# %%
# Third-party imports for the whole script live here -- Sections 3-5 rely on these names and do
# not re-import them.
import scanpy as sc
import matplotlib.pyplot as plt
import rpy2.robjects as ro
from scipy import sparse
from scipy.ndimage import gaussian_filter1d
from scipy.spatial import cKDTree
from matplotlib.colors import ListedColormap
from scipy.stats import ks_2samp, spearmanr, wasserstein_distance

from pseudospace.io_qc import (
    annotate_mito_ribo_mouse_symbols,
    combine_adatas,
    load_and_process_mouse_samples,
    sample_name_from_path,
)
from pseudospace.harmony import run_harmony_rpy2, select_harmony_hvgs_by_condition
from pseudospace.markers import (
    assign_cluster_labels,
    available_markers_by_group,
    build_gene_lookup,
    resolve_available_marker_groups,
)
from pseudospace.heatmaps import plot_marker_heatmap
from pseudospace.trajectory import (
    choose_diffusion_component,
    choose_root_global,
    choose_root_pt_cluster,
    orient_and_normalize,
    recompute_subset_dpt,
    save_total_pseudotime_anndata,
)

# Keep R on the active conda environment's library path; do not prepend ~/R/library.
plt.rcParams['figure.dpi'] = 120

# The cache redirects at the top of this file only work if they run BEFORE matplotlib and numba
# are first imported -- both read their cache-dir env var once, at import time. Verify rather
# than trust: moving that block below an import turns it into a silent no-op.
import matplotlib as _mpl

# Compare RESOLVED paths: on macOS /tmp is a symlink to /private/tmp, so matplotlib reports the
# resolved form and a plain string prefix test would warn on every run.
_mpl_cfg = Path(_mpl.get_configdir()).resolve()
if not _mpl_cfg.is_relative_to(CACHE_ROOT.resolve()):
    warnings.warn(
        f'matplotlib is caching to {_mpl_cfg}, not under {CACHE_ROOT.resolve()}. The MPLCONFIGDIR '
        '/ NUMBA_CACHE_DIR / XDG_CACHE_HOME block at the top of this file must execute before the '
        'first matplotlib import -- check that nothing imports matplotlib above it, and that '
        'MPLCONFIGDIR is not already set in the shell environment.', stacklevel=2)
else:
    print(f'Cache redirects active ({CACHE_ROOT}).')

try:
    import igraph  # noqa: F401
except ImportError as exc:
    raise ImportError(
        'python-igraph is required for sc.tl.leiden(flavor="igraph"). Install it in this '
        'environment before running Section 2.'
    ) from exc

print('rpy2 / R library path configured.')

# %%
files = sorted(path.name for path in TUBULE_BY_GENE_DIR.glob('*_tubule_by_gene_caleb.h5ad'))
files = [f for f in files if 'HUK' not in f]
if not files:
    raise FileNotFoundError(f'No mouse tubule-by-gene files found in {TUBULE_BY_GENE_DIR}')

found_samples = sorted(sample_name_from_path(Path(f)) for f in files)
if found_samples != sorted(MOUSE_SAMPLES):
    raise FileNotFoundError(
        f'Expected exactly {sorted(MOUSE_SAMPLES)} mouse samples, found {found_samples}. '
        'Check data/tubule_by_gene/ for missing *_tubule_by_gene_caleb.h5ad files.'
    )
print(f'Discovered {len(files)} mouse sample files: {found_samples}')

mouse_adatas = load_and_process_mouse_samples(TUBULE_BY_GENE_DIR, files)
adata_combined = combine_adatas(mouse_adatas)
adata_combined = annotate_mito_ribo_mouse_symbols(adata_combined)
print(f'\nCombined: {adata_combined.n_obs:,} tubules x {adata_combined.n_vars:,} genes')
print(adata_combined.obs.groupby('condition')['sample'].value_counts())

# %%
# Capture the per-tubule counts BEFORE the tubule filter runs: once it does, the removed
# tubules are gone and the distribution it acted on can no longer be plotted. `total_counts`
# is the raw UMI sum of the structure -- the quantity the gene-level support filter
# (MIN_GENE_TUBULE_FRACTION, MIN_GENE_TOTAL_COUNTS) is built from.
_tubule_qc_before = pd.DataFrame({
    'n_genes_by_counts': np.asarray((adata_combined.X > 0).sum(axis=1)).ravel(),
    'total_counts': np.asarray(adata_combined.X.sum(axis=1)).ravel(),
    'sample': adata_combined.obs['sample'].astype(str).to_numpy(),
    'condition': adata_combined.obs['condition'].astype(str).to_numpy(),
}, index=adata_combined.obs_names)
print(f'Before cell filtering: {adata_combined.n_obs} tubules x {adata_combined.n_vars} genes')
sc.pp.filter_cells(adata_combined, min_genes=MIN_GENES_PER_TUBULE)
print(f'After cell filtering:  {adata_combined.n_obs} tubules x {adata_combined.n_vars} genes')

# The same measurement once the filter has run, for the before/after plot below. Taken here,
# before the gene filter drops columns, so it reflects the tubule filter alone.
_tubule_qc_after = pd.DataFrame({
    'n_genes_by_counts': np.asarray((adata_combined.X > 0).sum(axis=1)).ravel(),
    'total_counts': np.asarray(adata_combined.X.sum(axis=1)).ravel(),
    'sample': adata_combined.obs['sample'].astype(str).to_numpy(),
    'condition': adata_combined.obs['condition'].astype(str).to_numpy(),
}, index=adata_combined.obs_names)

# --- Low-support gene filter, now BEFORE Harmony ---------------------------------------------
# This ran after Harmony, on the object read back from disk, so X_harmony was built with genes the
# later gene-level analysis then discarded. Notebook 03 filters ahead of Harmony; this matches it.
# The whitelist is unchanged mouse-specific content: curated and diagnostic markers stay in the
# matrix even when they fail the support thresholds.
_gene_counts = adata_combined.X if sparse.issparse(adata_combined.X) else np.asarray(adata_combined.X)
_gene_expressed = np.asarray((_gene_counts > 0).sum(axis=0)).ravel()
_gene_total_counts = np.asarray(_gene_counts.sum(axis=0)).ravel()
_min_cells = max(1, int(np.ceil(MIN_GENE_TUBULE_FRACTION * adata_combined.n_obs)))
gene_keep = (_gene_expressed >= _min_cells) & (_gene_total_counts >= MIN_GENE_TOTAL_COUNTS)
_whitelist_hits = adata_combined.var_names.str.upper().isin({g.upper() for g in GENE_FILTER_WHITELIST})
gene_keep = gene_keep | np.asarray(_whitelist_hits)
adata_combined.var['n_tubules_expressed'] = _gene_expressed
adata_combined.var['total_counts'] = _gene_total_counts
adata_combined.var['passes_expression_count_filter'] = gene_keep
print(f'Gene filter before Harmony keeps {int(gene_keep.sum()):,}/{adata_combined.n_vars:,} genes '
      f'(>= {MIN_GENE_TUBULE_FRACTION:.0%} of tubules and >= {MIN_GENE_TOTAL_COUNTS} counts, '
      f'plus {int(np.asarray(_whitelist_hits).sum())} whitelisted markers)')
adata_combined = adata_combined[:, gene_keep].copy()

adata_combined.layers['counts'] = adata_combined.X.copy()
sc.pp.normalize_total(adata_combined, target_sum=NORMALIZE_TARGET_SUM)
sc.pp.log1p(adata_combined)
adata_combined.layers['lognorm'] = adata_combined.X.copy()
print("Normalized; counts in layers['counts'], log-normalized values in layers['lognorm'].")

# %%
OUTPUT_DIR = CELLTYPING_DIR

# --- Per-tubule size around the tubule filter: before and after, per sample ---------------------
# Several downstream steps are sensitive to tubule size (HVG selection, the subset DPTs, the
# gene-level shape analysis), so cutting one sample harder than another is a real confound rather
# than a cosmetic QC number. `_tubule_qc_before` is the pre-filter state captured in the cell
# above; `_tubule_qc_after` is the same measurement once the filter has run.
_before = _tubule_qc_before
_before['passes_min_genes'] = _before['n_genes_by_counts'] >= MIN_GENES_PER_TUBULE
_after = _tubule_qc_after
_before.to_csv(OUTPUT_DIR / 'tubule_counts_before_filter.csv')
_after.to_csv(OUTPUT_DIR / 'tubule_counts_after_filter.csv')
_retention = _before.groupby('sample', observed=True)['passes_min_genes'].agg(['size', 'sum', 'mean'])
_retention['fraction_removed'] = 1 - _retention['mean']
_retention.to_csv(OUTPUT_DIR / 'tubule_retention_by_sample.csv')
print(f"Tubule filter (>= {MIN_GENES_PER_TUBULE} detected genes): "
      f"{int(_before['passes_min_genes'].sum()):,}/{len(_before):,} tubules pass, "
      f"{1 - _before['passes_min_genes'].mean():.2%} removed")
print(_retention.to_string())

_samples = sorted(_before['sample'].unique())
_colors = dict(zip(_samples, plt.get_cmap('tab10').colors))
_GENE_COLUMN = 'shared_n_genes' if 'shared_n_genes' in _before.columns else 'n_genes_by_counts'
_NOUN = 'structure' if _GENE_COLUMN == 'shared_n_genes' else 'tubule'
_GENE_AXIS = ('detected genes per structure (shared ortholog space)' if _GENE_COLUMN == 'shared_n_genes'
              else 'detected genes per tubule')

fig, axes = plt.subplots(2, 3, figsize=(21, 10))

for column, (frame, state) in enumerate(((_before, 'before'), (_after, 'after'))):
    for sample in _samples:
        sub = frame[frame['sample'] == sample]
        axes[0, column].hist(sub[_GENE_COLUMN], bins=100, histtype='step', lw=1.2,
                             color=_colors[sample], label=f'{sample} (n={len(sub):,})')
    if state == 'before':
        axes[0, column].axvline(MIN_GENES_PER_TUBULE, color='crimson', ls='--', lw=1.6,
                                label=f'threshold = {MIN_GENES_PER_TUBULE} genes')
    axes[0, column].set_yscale('log')
    axes[0, column].set_xlabel(_GENE_AXIS)
    axes[0, column].set_ylabel(f'{_NOUN}s (log scale)')
    axes[0, column].set_title(f'Detected genes per {_NOUN}, {state} filtering')
    axes[0, column].legend(fontsize=7)

    for sample in _samples:
        sub = frame[frame['sample'] == sample]
        axes[1, column].hist(np.log10(sub['total_counts'].to_numpy() + 1), bins=100,
                             histtype='step', lw=1.2, color=_colors[sample], label=sample)
    axes[1, column].set_xlabel(f'log10(total counts + 1) per {_NOUN}')
    axes[1, column].set_ylabel(f'{_NOUN}s')
    axes[1, column].set_title(f'Total counts per {_NOUN}, {state} filtering')
    axes[1, column].legend(fontsize=7)

# How much each sample lost, and what that did to its size distribution.
_removed = (1 - _retention['mean'].reindex(_samples).to_numpy()) * 100
axes[0, 2].bar(_samples, _removed, color=[_colors[sample] for sample in _samples])
for index, value in enumerate(_removed):
    axes[0, 2].text(index, value, f'{value:.0f}%', ha='center', va='bottom', fontsize=9)
axes[0, 2].set_ylim(0, max(_removed) * 1.25)
axes[0, 2].set_ylabel(f'% of {_NOUN}s removed')
axes[0, 2].set_title('Removed by the tubule filter, per sample')
axes[0, 2].tick_params(axis='x', rotation=30)

_positions = np.arange(len(_samples)) * 3.0
axes[1, 2].boxplot([np.log10(_before.loc[_before['sample'] == sample, 'total_counts'].to_numpy() + 1)
                    for sample in _samples], positions=_positions - 0.6, widths=0.9,
                   showfliers=False, patch_artist=True,
                   boxprops=dict(facecolor='#CCCCCC', edgecolor='#888888'),
                   medianprops=dict(color='black'))
axes[1, 2].boxplot([np.log10(_after.loc[_after['sample'] == sample, 'total_counts'].to_numpy() + 1)
                    for sample in _samples], positions=_positions + 0.6, widths=0.9,
                   showfliers=False, patch_artist=True,
                   boxprops=dict(facecolor='#7FB3D5', edgecolor='#2E6DA4'),
                   medianprops=dict(color='black'))
axes[1, 2].set_xticks(_positions)
axes[1, 2].set_xticklabels(_samples, rotation=30, ha='right')
axes[1, 2].set_ylabel('log10(total counts + 1)')
axes[1, 2].set_title('Size per sample, before (grey) vs after (blue)')

fig.suptitle(f'{_NOUN.capitalize()} size around the tubule filter '
             f'(>= {MIN_GENES_PER_TUBULE} detected genes)', fontsize=14)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(OUTPUT_DIR / 'tubule_count_distribution.png', dpi=200, bbox_inches='tight')
plt.show()


# %%
adata_combined = select_harmony_hvgs_by_condition(
    adata_combined,
    group_key='condition',
    groups=('Control', 'IR'),
    mode='intersection',
    min_mean=HARMONY_HVG_MIN_MEAN,
    max_mean=HARMONY_HVG_MAX_MEAN,
    min_disp=HARMONY_HVG_MIN_DISP,
)
hvg_col = 'highly_variable_for_harmony'
print(f'\nFinal Harmony features: {int(adata_combined.var[hvg_col].sum())}')
print('HVG condition-count distribution:')
print(adata_combined.var['hvg_condition_count'].value_counts().sort_index())

# %%
# Uncorrected PCA + pre-Harmony UMAP (the "before" panel for the integration check).
adata_hvg = adata_combined[:, adata_combined.var[hvg_col]].copy()
adata_hvg.X = adata_hvg.layers['lognorm'].copy()
sc.tl.pca(adata_hvg, n_comps=HARMONY_PCA_N_COMPS, random_state=RANDOM_STATE)

adata_combined.obsm['X_pca'] = adata_hvg.obsm['X_pca'].copy()
sc.pp.neighbors(adata_combined, use_rep='X_pca', n_neighbors=HARMONY_NEIGHBORS_N,
                key_added='pre_harmony', random_state=RANDOM_STATE)
sc.tl.umap(adata_combined, min_dist=HARMONY_UMAP_MIN_DIST, spread=HARMONY_UMAP_SPREAD,
           neighbors_key='pre_harmony', random_state=RANDOM_STATE)
adata_combined.obsm['X_umap_pre_harmony'] = adata_combined.obsm['X_umap'].copy()

sc.pl.umap(adata_combined, color=['sample', 'condition'], ncols=2, frameon=False,
           title=['Pre-Harmony: sample', 'Pre-Harmony: condition'])

# %%
print(f'\nRunning Harmony on {adata_hvg.n_obs} tubules x {adata_hvg.n_vars} selected features '
      f'(theta={HARMONY_THETA})')
print(adata_hvg.obs['sample'].value_counts())

# Seed R's global RNG before every harmony call. harmony initialises its soft-kmeans centroids
# randomly, so an unseeded run returns a DIFFERENT X_harmony each time -- which moves the
# neighbour graph, which moves Leiden, which invalidates the hand-written COARSE_LABELS below.
# rpy2 shares one R session, so setting the seed here governs the RunHarmony call that follows.
ro.r(f'set.seed({RANDOM_STATE})')
adata_hvg = run_harmony_rpy2(
    adata_hvg, batch_key=BATCH_KEY, n_pcs=HARMONY_N_PCS, theta=HARMONY_THETA,
    lambda_val=HARMONY_LAMBDA, max_iter=HARMONY_MAX_ITER, tau=HARMONY_TAU,
)
adata_combined.obsm['X_harmony'] = adata_hvg.obsm['X_harmony'].copy()
adata_combined.obsm['X_pca'] = adata_hvg.obsm['X_pca'].copy()

sc.pp.neighbors(adata_combined, use_rep='X_harmony', n_neighbors=HARMONY_NEIGHBORS_N,
                random_state=RANDOM_STATE)
sc.tl.umap(adata_combined, min_dist=HARMONY_UMAP_MIN_DIST, spread=HARMONY_UMAP_SPREAD,
           random_state=RANDOM_STATE)
sc.pl.umap(adata_combined, color=['sample', 'condition'], ncols=2, frameon=False)

# %%
adata_combined.write(HARMONY_OUTPUT_PATH)
print(f'\nSaved pass-1 harmonized AnnData to: {HARMONY_OUTPUT_PATH.relative_to(PROJECT_DIR)}')
print('\n--- Section 1 summary ---')
print(f'Tubules: {adata_combined.n_obs:,}')
print(f'Genes (native mouse space): {adata_combined.n_vars:,}')
print(f'Harmony features: {int(adata_combined.var[hvg_col].sum()):,}')
print(f'X_harmony shape: {adata_combined.obsm["X_harmony"].shape}')

# %% [markdown]
# # Section 2 - cell typing, pass-2 Harmony, Scanpy DPT
#
# 2a whitelist-protected gene filter -> 2b coarse Leiden + NEPHRON_AXIS_MARKERS dotplot ->
# **manual checkpoint** -> 2c rollup + non-tubule removal -> 2d pass-2 Harmony on the tubule
# subset (no re-clustering, no re-labelling) -> 2e marker axis -> 2f diffmap + root + DPT.

# %%
# --- Section 2a: whitelist-protected gene filter ---
adata_all = sc.read_h5ad(HARMONY_OUTPUT_PATH)
print(f'Loaded pass-1 object: {adata_all.n_obs:,} tubules x {adata_all.n_vars:,} genes')

if 'counts' in adata_all.layers:
    counts = adata_all.layers['counts'].copy()
else:
    counts = adata_all.X.copy()
    print('No counts layer found; using current X as the raw-count source')
adata_all.layers['counts'] = counts.copy()

# The low-support gene filter now runs BEFORE Harmony (Section 1), matching notebook 03, so this
# object arrives already filtered. Assert that rather than re-applying the filter: a silent
# re-application would hide a Section 1 mistake instead of failing on it.
if 'passes_expression_count_filter' not in adata_all.var.columns:
    raise RuntimeError(
        'Pass-1 object carries no gene filter: Section 1 must apply it before Harmony. '
        f'Rebuild {HARMONY_OUTPUT_PATH.name} before continuing.'
    )
_still_filtered = int((~adata_all.var['passes_expression_count_filter'].astype(bool)).sum())
if _still_filtered:
    raise RuntimeError(
        f'Pass-1 object still holds {_still_filtered:,} filtered genes; Section 1 did not subset them.'
    )
print(f'Gene filter already applied before Harmony: {adata_all.n_vars:,} genes retained '
      f'(>= {MIN_GENE_TUBULE_FRACTION:.0%} of tubules and >= {MIN_GENE_TOTAL_COUNTS} counts, '
      f'plus whitelisted markers)')

adata_all.X = adata_all.layers['counts'].copy()
adata_all.uns.pop('log1p', None)
sc.pp.normalize_total(adata_all, target_sum=NORMALIZE_TARGET_SUM)
sc.pp.log1p(adata_all)
adata_all.layers['lognorm'] = adata_all.X.copy()

# Snapshot before HVG subsetting -- this is the expression matrix written to the DPT file.
adata_all_expression_filtered = adata_all.copy()
print(f'After gene filter + normalization: {adata_all.n_obs:,} tubules x {adata_all.n_vars:,} genes')

# %%
COARSE_RESOLUTION

# %%
# --- Section 2b: the ONE clustering run, on the pass-1 Harmony embedding ---
gene_lookup = build_gene_lookup(adata_all)

sc.pp.highly_variable_genes(adata_all, n_top_genes=min(N_HVGS, adata_all.n_vars), flavor='seurat')
whitelist_present = {gene_lookup[g.upper()] for g in CURATED_MARKER_WHITELIST if g.upper() in gene_lookup}
feature_mask = (adata_all.var['highly_variable'].to_numpy(dtype=bool)
                | adata_all.var_names.isin(whitelist_present))
adata_all.var['selected_for_clustering'] = feature_mask
adata_cluster = adata_all[:, feature_mask].copy()
gene_lookup = build_gene_lookup(adata_cluster)
print(f'Clustering feature space: {int(feature_mask.sum()):,} genes '
      f'({int(adata_all.var["highly_variable"].sum()):,} HVGs + curated markers)')

sc.pp.neighbors(adata_cluster, n_neighbors=N_NEIGHBORS, use_rep='X_harmony', random_state=RANDOM_STATE)
sc.tl.leiden(adata_cluster, resolution=COARSE_RESOLUTION, key_added='leiden_coarse',
             flavor='igraph', n_iterations=2, directed=False, random_state=RANDOM_STATE)
n_clusters = adata_cluster.obs['leiden_coarse'].nunique()
print(f'Leiden (res={COARSE_RESOLUTION}): {n_clusters} clusters')

# Fingerprint of everything COARSE_LABELS is keyed against. The hand-written label dict below is
# valid only for THIS clustering: change the resolution, the feature space, the gene filter or
# the random state and cluster id '3' stops meaning what the dict says it means. The membership
# hash catches the dangerous case the id-coverage guards cannot -- same cluster COUNT, different
# cluster CONTENTS, which passes every other check while silently relabelling the whole dataset.
import hashlib

CLUSTERING_FINGERPRINT = {
    'n_cells': int(adata_cluster.n_obs),
    'n_features': int(adata_cluster.n_vars),
    'resolution': COARSE_RESOLUTION,
    'n_neighbors': N_NEIGHBORS,
    'random_state': RANDOM_STATE,
    'n_clusters': int(n_clusters),
    'membership_sha1': hashlib.sha1(
        ','.join(adata_cluster.obs['leiden_coarse'].astype(str)).encode()).hexdigest()[:12],
}
print('Clustering fingerprint:')
for _k, _v in CLUSTERING_FINGERPRINT.items():
    print(f'    {_k:>16s}: {_v}')
sc.tl.rank_genes_groups(adata_cluster, 'leiden_coarse', method='wilcoxon', pts=True,
                        key_added='rank_coarse')

# %%
# Dotplot over the full nephron axis + the non-tubule panels: the visual basis for the labels.
coarse_panel = {**NEPHRON_AXIS_MARKERS, **NON_TUBULE_MARKERS}
resolved_panel, missing_panel = available_markers_by_group(adata_cluster, coarse_panel, gene_lookup)
if missing_panel:
    print('Markers absent from the gene space (skipped):', missing_panel)

sc.pl.dotplot(adata_cluster, resolved_panel, groupby='leiden_coarse', standard_scale='var', show=False)
plt.savefig(CELLTYPING_DIR / 'coarse_marker_dotplot.png', dpi=150, bbox_inches='tight')
plt.show()

sc.pl.embedding(adata_cluster, basis='umap', color='leiden_coarse', frameon=False,
                legend_loc='on data', show=False)
plt.savefig(CELLTYPING_DIR / 'coarse_leiden_umap.png', dpi=150, bbox_inches='tight')
plt.show()

top_de = sc.get.rank_genes_groups_df(adata_cluster, group=None, key='rank_coarse')
top5 = (top_de.sort_values(['group', 'scores'], ascending=[True, False])
        .groupby('group').head(5).groupby('group')['names'].apply(lambda s: '; '.join(s)))
top5.to_csv(CELLTYPING_DIR / 'coarse_cluster_top_markers.csv')
print(top5.to_string())

# %%
# Per-cluster mean score for every panel in the dotplot -- a SUGGESTION to seed COARSE_LABELS.
# It is not authoritative: read the dotplot, then decide.
score_columns = {}
for panel_name, genes in resolved_panel.items():
    score_name = f'_score_{panel_name}'
    sc.tl.score_genes(adata_cluster, genes, score_name=score_name, use_raw=False,
                      random_state=RANDOM_STATE)
    score_columns[panel_name] = score_name

panel_scores = (adata_cluster.obs.groupby('leiden_coarse', observed=True)[list(score_columns.values())]
                .mean().rename(columns={v: k for k, v in score_columns.items()}))
panel_scores.to_csv(CELLTYPING_DIR / 'coarse_cluster_panel_scores.csv')
suggested_labels = {str(cluster): str(panel_scores.columns[int(np.argmax(row.to_numpy()))])
                    for cluster, row in panel_scores.iterrows()}

# Per-cell argmax from the same scores: a label-independent second opinion, and the confidence
# margin between the best and runner-up panel. `celltype_primary`/`celltype_secondary` are rolled
# up to the coarse vocabulary (the spatial validation script expects 'Glomerulus'/'Vessel' there);
# the fine per-cell call is kept alongside in the *_segment columns.
score_matrix = adata_cluster.obs[list(score_columns.values())].to_numpy(dtype=float)
panel_names = np.array(list(score_columns))
ranked = np.argsort(score_matrix, axis=1)
primary_segment = panel_names[ranked[:, -1]]
secondary_segment = panel_names[ranked[:, -2]]
adata_cluster.obs['celltype_primary_segment'] = primary_segment
adata_cluster.obs['celltype_secondary_segment'] = secondary_segment
adata_cluster.obs['celltype_primary'] = [SEGMENT_TO_COARSE.get(p, p) for p in primary_segment]
adata_cluster.obs['celltype_secondary'] = [SEGMENT_TO_COARSE.get(p, p) for p in secondary_segment]
adata_cluster.obs['celltype_margin'] = (
    np.take_along_axis(score_matrix, ranked[:, -1:], 1).ravel()
    - np.take_along_axis(score_matrix, ranked[:, -2:-1], 1).ravel())
PER_CELL_CALL_COLUMNS = ['celltype_primary', 'celltype_secondary',
                         'celltype_primary_segment', 'celltype_secondary_segment',
                         'celltype_margin']
adata_cluster.obs.drop(columns=list(score_columns.values()), inplace=True)

display(panel_scores.round(3))
print('\nPer-cell argmax vs cluster id (rows = leiden_coarse):')
display(pd.crosstab(adata_cluster.obs['leiden_coarse'], adata_cluster.obs['celltype_primary_segment']))
print('\nSuggested COARSE_LABELS (argmax panel score -- VERIFY against the dotplot):')
print('COARSE_LABELS = {')
for cluster, label in suggested_labels.items():
    print(f"    {cluster!r}: {label!r},")
print('}')

# %% [markdown]
# ## Manual checkpoint - label every cluster
#
# Inspect `celltyping/coarse_marker_dotplot.png`, `coarse_cluster_top_markers.csv` and the
# suggestion printed above, then fill in `COARSE_LABELS` below. Valid values are the 17 fine
# segments (`Podocyte, PT-S1 ... IMCD`), the COARSE families (`PT, DTL, AL, DCT, CNT_CD`),
# `Glomerulus`, plus `Vessel, Stroma, SmoothMuscle, Immune, Unassigned`. Label each cluster at
# the granularity its markers actually support -- a coarse label is the CORRECT answer when a
# cluster does not resolve a sub-segment, not a fallback.
# This is the only labelling step in the pipeline: the pass-2 Harmony below re-integrates the
# tubule subset but does NOT re-cluster or re-label.

# %%
# Manual labels for the reproducible R-harmony 2.0.5 run (fingerprint below). Re-reviewed against
# celltyping/coarse_marker_dotplot.png and coarse_cluster_top_markers.csv after the low-support gene
# filter moved ahead of Harmony (Section 1), which changed the clustering while leaving the cell set
# and the cluster count unchanged. The Leiden IDs shifted, so every entry was re-derived from the
# current run's panel scores and DE; the non-nephron programs are unchanged, only renumbered:
#   9 = Cryab/Vim/Tnc/Thbs1 (Stroma), 8 = Krt19/Sprr1a/Krt18 (Unassigned).
# Two calls are the least certain and worth your eye: 6/7 split the thick ascending limb (cTAL +1.63
# vs mTAL +1.48; mTAL +1.59 vs cTAL +1.40) and 10 is the weakest collecting-duct call (IMCD +1.15
# over ATL +0.95, with Aqp2 in the DE).
# They are valid ONLY for the fingerprint below; a changed clustering must be reviewed.
COARSE_LABELS = {
    '0': 'PT-S1',        # panel +1.66; DE Alpl/Slc5a2/Gatm/Slc34a1
    '1': 'PT-S2',        # panel +1.32; DE Slc22a6/Slc13a3
    '2': 'Podocyte',     # panel +2.33; DE Podxl/Synpo/Nphs2
    '3': 'DCT2',         # panel DCT2 +1.44 over CNT +1.10; DE Slc8a1/Clcnkb
    '4': 'PT-S3',        # panel +0.91; DE Slc6a18/Napsa/Mep1a
    '5': 'CCD',          # panel +1.26; DE Aqp2/Aqp3/Hsd11b2
    '6': 'cTAL',         # panel cTAL +1.63 vs mTAL +1.48; DE Umod/Slc12a1
    '7': 'mTAL',         # panel mTAL +1.59 vs cTAL +1.40; DE Umod/Ppp1r1a
    '8': 'Unassigned',   # DE Krt19/Sprr1a/Krt18; best panel only ATL +0.80
    '9': 'Stroma',       # DE Cryab/Vim/Tnc/Thbs1; all panels <= +0.40
    '10': 'IMCD',        # panel IMCD +1.15 over ATL +0.95; DE Aqp2 -- weakest call
    '11': 'SmoothMuscle',  # panel +2.44; DE Acta2/Myh11/Tagln
}

# The reference run used R harmony 2.0.5 and Scanpy/Leiden resolution 0.7. Do not disable this
# guard: the labels above are keyed by Leiden ID, not by marker identity.
COARSE_LABELS_FINGERPRINT = {
    'n_cells': 43848,
    'n_clusters': 12,
    'resolution': 0.7,
    'n_neighbors': 30,
    'random_state': 0,
    'membership_sha1': '86657cf8af29',
}

if not COARSE_LABELS:
    raise ValueError('Fill in COARSE_LABELS from the marker dotplot before proceeding.')
_unknown = sorted(set(COARSE_LABELS.values()) - set(LABEL_VOCABULARY))
if _unknown:
    raise ValueError(f'COARSE_LABELS contains labels outside LABEL_VOCABULARY: {_unknown}')
_unlabelled = sorted(set(adata_cluster.obs['leiden_coarse'].astype(str)) - set(COARSE_LABELS))
if _unlabelled:
    raise ValueError(f'Clusters {_unlabelled} have no entry in COARSE_LABELS; review the dotplot.')
_stale = sorted(set(COARSE_LABELS) - set(adata_cluster.obs['leiden_coarse'].astype(str)))
if _stale:
    raise ValueError(f'COARSE_LABELS has stale cluster IDs {_stale}; review the dotplot.')
_diff = {k: (v, CLUSTERING_FINGERPRINT[k]) for k, v in COARSE_LABELS_FINGERPRINT.items()
         if v != CLUSTERING_FINGERPRINT[k]}
if _diff:
    raise ValueError(
        'This clustering does not match the one COARSE_LABELS was written against. '
        + '; '.join(f'{k}: {a} -> {b}' for k, (a, b) in sorted(_diff.items()))
        + '. Re-read celltyping/coarse_marker_dotplot.png and update both labels and fingerprint.'
    )
print(f'COARSE_LABELS pinned to clustering {CLUSTERING_FINGERPRINT["membership_sha1"]}.')

# %%
# --- Section 2c: fine label -> coarse rollup -> drop non-tubule cells ---
assign_cluster_labels(adata_cluster, 'leiden_coarse', COARSE_LABELS, 'segment_class')

adata_all.obs['segment_class'] = adata_cluster.obs['segment_class'].reindex(adata_all.obs_names).astype(str).to_numpy()
adata_all.obs['leiden_coarse'] = adata_cluster.obs['leiden_coarse'].reindex(adata_all.obs_names).astype(str).to_numpy()
for col in PER_CELL_CALL_COLUMNS:
    values = adata_cluster.obs[col].reindex(adata_all.obs_names)
    adata_all.obs[col] = values.to_numpy() if col == 'celltype_margin' else values.astype(str).to_numpy()
# Non-tubule labels pass through the rollup unchanged; anything unmapped becomes 'Unassigned'.
adata_all.obs['coarse_class'] = (
    adata_all.obs['segment_class']
    .map(segment_vocabulary.coarse_for)
    .astype(str)
)
# Compatibility alias: the pseudospace package, the QuPath export and the spatial validation
# script all key off this name.
adata_all.obs['broad_tubule_marker_call'] = adata_all.obs['coarse_class']

print(pd.crosstab(adata_all.obs['segment_class'], adata_all.obs['coarse_class']).to_string())

# Which parts of the nephron did the labels actually cover? An absent family is not an error --
# downstream code skips it -- but it is skipped SILENTLY, three sections later, so surface it now.
_present_coarse = set(adata_all.obs['coarse_class'].astype(str))
_missing_families = [f for f in COARSE_ORDER if f not in _present_coarse]
if _missing_families:
    print(f'\nNOTE: no cluster was labelled as {_missing_families}. Consequences: the subset DPT '
          f'for {_missing_families} is skipped (Section 3), and any heatmap panel scoped to those '
          'classes plots only its remaining families. If a family should be there, it is either '
          'merged into a neighbour at this resolution or genuinely absent from the sections.')
if 'Glomerulus' not in _present_coarse:
    print('\nNOTE: no cluster was labelled Podocyte, so there is no Glomerulus class. The '
          'physical-axis sensitivity analysis (Section 4.9) and the three-axis concordance '
          '(Section 5) both derive their depth proxy from glomerular positions and will be '
          'SKIPPED. Label a podocyte cluster if you need those checks.')
sc.pl.embedding(adata_cluster, basis='umap', color=['leiden_coarse'], frameon=False, show=False)
plt.savefig(CELLTYPING_DIR / 'labelled_umap.png', dpi=150, bbox_inches='tight')
plt.show()

# %%
# Persist the labels onto the pass-1 file (ALL cells, including the removed classes) so the
# QuPath export and the spatial validation keep their glomerulus + vessel locations.
adata_pass1 = sc.read_h5ad(HARMONY_OUTPUT_PATH)
for col in ['segment_class', 'coarse_class', 'broad_tubule_marker_call', 'leiden_coarse',
            *PER_CELL_CALL_COLUMNS]:
    if col == 'celltype_margin':
        adata_pass1.obs[col] = adata_all.obs.loc[adata_pass1.obs_names, col].to_numpy(dtype=float)
    else:
        adata_pass1.obs[col] = adata_all.obs.loc[adata_pass1.obs_names, col].astype(str).to_numpy()
adata_pass1.write(HARMONY_OUTPUT_PATH)
print(f'Wrote labels to {HARMONY_OUTPUT_PATH.name} (all {adata_pass1.n_obs:,} cells).')
del adata_pass1

keep_mask = adata_all.obs['coarse_class'].isin(KEEP_TUBULE_CLASSES).to_numpy()
for cls in REMOVE_CLASSES:
    n_removed = int((adata_all.obs['coarse_class'] == cls).sum())
    if n_removed:
        print(f'  removing {n_removed:,} {cls} cells')
adata_tubule = adata_all[keep_mask].copy()
adata_tubule.obs['coarse_class'] = pd.Categorical(
    adata_tubule.obs['coarse_class'], categories=COARSE_ORDER, ordered=True)
adata_tubule.obs['segment_class'] = pd.Categorical(
    adata_tubule.obs['segment_class'],
    categories=[s for s in SEGMENT_DISPLAY_ORDER if s in set(adata_tubule.obs['segment_class'])],
    ordered=True)
print(f'Tubule-only object: {adata_tubule.n_obs:,} tubules x {adata_tubule.n_vars:,} genes')
print(adata_tubule.obs['coarse_class'].value_counts().reindex(COARSE_ORDER).to_string())

# %%
# --- Section 2d: pass-2 Harmony on the tubule subset (labels are carried, never recomputed) ---
adata_tubule = select_harmony_hvgs_by_condition(
    adata_tubule, group_key='condition', groups=('Control', 'IR'), mode='intersection',
    min_mean=HARMONY_HVG_MIN_MEAN, max_mean=HARMONY_HVG_MAX_MEAN, min_disp=HARMONY_HVG_MIN_DISP)
tub_hvg = 'highly_variable_for_harmony'
print(f'Pass-2 HVGs (intersection, tubule-only): {int(adata_tubule.var[tub_hvg].sum())}')

adata_tubule_hvg = adata_tubule[:, adata_tubule.var[tub_hvg]].copy()
adata_tubule_hvg.X = adata_tubule_hvg.layers['lognorm'].copy()
sc.tl.pca(adata_tubule_hvg, n_comps=HARMONY_PCA_N_COMPS, random_state=RANDOM_STATE)
adata_tubule.obsm['X_pca'] = adata_tubule_hvg.obsm['X_pca'].copy()

ro.r(f'set.seed({RANDOM_STATE})')          # see the note at the Section 1 harmony call
adata_tubule_hvg = run_harmony_rpy2(
    adata_tubule_hvg, batch_key=BATCH_KEY, n_pcs=HARMONY_N_PCS, theta=HARMONY_THETA,
    lambda_val=HARMONY_LAMBDA, max_iter=HARMONY_MAX_ITER, tau=HARMONY_TAU)
adata_tubule.obsm['X_harmony'] = adata_tubule_hvg.obsm['X_harmony'].copy()

sc.pp.neighbors(adata_tubule, use_rep='X_harmony', n_neighbors=HARMONY_NEIGHBORS_N,
                random_state=RANDOM_STATE)
sc.tl.umap(adata_tubule, min_dist=HARMONY_UMAP_MIN_DIST, spread=HARMONY_UMAP_SPREAD,
           random_state=RANDOM_STATE)
adata_tubule.write(PASS2_HARMONY_OUTPUT_PATH)

sc.pl.embedding(adata_tubule, basis='umap', color=['sample', 'condition', 'coarse_class', 'segment_class'],
                ncols=2, frameon=False, show=False)
plt.savefig(CELLTYPING_DIR / 'pass2_umap.png', dpi=150, bbox_inches='tight')
plt.close()
print(f'Saved pass-2 tubule-only Harmony object to {PASS2_HARMONY_OUTPUT_PATH.name}')

# %%
# --- Section 2d.1: integration QC -- is theta=6 over-correcting? ---
# Harmony's default theta is 2; this pipeline uses 6, on the embedding that DPT then consumes.
# Over-correction is the silent failure of integration: it snaps continuous gradients toward
# shared anchors, which is exactly the structure pseudospace is trying to measure. Batch-mixing
# metrics are trivially maximised by over-correction (a method that destroys all structure mixes
# batches perfectly), so batch ASW is only meaningful READ NEXT TO bio ASW.
#
# Read it as: batch ASW should FALL (batches mixed) while bio ASW should HOLD (families kept).
# batch ASW falling and bio ASW also falling = over-correction; revisit HARMONY_THETA.
from sklearn.metrics import silhouette_score

_rng = np.random.default_rng(RANDOM_STATE)
_qc_idx = (_rng.choice(adata_tubule.n_obs, INTEGRATION_QC_MAX_CELLS, replace=False)
           if adata_tubule.n_obs > INTEGRATION_QC_MAX_CELLS else np.arange(adata_tubule.n_obs))
_batch_labels = adata_tubule.obs[BATCH_KEY].astype(str).to_numpy()[_qc_idx]
_bio_labels = adata_tubule.obs['coarse_class'].astype(str).to_numpy()[_qc_idx]

integration_qc = pd.DataFrame([
    {'embedding': name,
     'batch_asw_lower_is_better': float(silhouette_score(adata_tubule.obsm[name][_qc_idx], _batch_labels)),
     'bio_asw_higher_is_better': float(silhouette_score(adata_tubule.obsm[name][_qc_idx], _bio_labels))}
    for name in ['X_pca', 'X_harmony'] if name in adata_tubule.obsm
])
integration_qc['theta'] = HARMONY_THETA
integration_qc['n_cells_scored'] = len(_qc_idx)
integration_qc.to_csv(CELLTYPING_DIR / 'integration_qc_silhouette.csv', index=False)
display(integration_qc.round(4))

if {'X_pca', 'X_harmony'} <= set(adata_tubule.obsm):
    _pca_row = integration_qc.set_index('embedding').loc['X_pca']
    _har_row = integration_qc.set_index('embedding').loc['X_harmony']
    _batch_gain = _pca_row['batch_asw_lower_is_better'] - _har_row['batch_asw_lower_is_better']
    _bio_cost = _pca_row['bio_asw_higher_is_better'] - _har_row['bio_asw_higher_is_better']
    print(f'\nHarmony (theta={HARMONY_THETA}) removed {_batch_gain:+.4f} batch ASW '
          f'at a cost of {_bio_cost:+.4f} bio ASW.')
    if _bio_cost > 0 and _batch_gain > 0 and _bio_cost > 0.5 * _batch_gain:
        print('WARNING: bio ASW fell by more than half the batch ASW gain -- this is the '
              f'over-correction signature. Consider lowering HARMONY_THETA (currently '
              f'{HARMONY_THETA}; Harmony default is 2) and re-checking the DPT/segment ordering.')
    else:
        print('Bio conservation held relative to the batch-mixing gain; theta looks defensible.')
print('NOTE: this is ASW only. scib_metrics (iLISI/cLISI/isolated-label F1) is not installed in '
      'this environment; `pip install scib-metrics` for the full scIB panel.')

# %%
# --- Section 2e: the continuum + the early->late marker axis ---
CONTINUUM_PREFIX = 'total'
tubule_lookup = build_gene_lookup(adata_tubule)

exclude_mask = adata_tubule.obs['coarse_class'].astype(str).isin(['Unassigned', 'nan', 'None']).to_numpy()
adata_total = adata_tubule[~exclude_mask].copy()
print(f'TOTAL continuum: {adata_total.n_obs:,} tubules (excluded {int(exclude_mask.sum()):,})')

for role in ['early', 'late']:
    genes = [tubule_lookup[g.upper()] for g in TOTAL_POSITION_MARKERS[role] if g.upper() in tubule_lookup]
    if not genes:
        raise ValueError(f'No {role} axis markers present in the gene space: {TOTAL_POSITION_MARKERS[role]}')
    print(f'  {role} axis genes ({len(genes)}): {genes}')
    sc.tl.score_genes(adata_total, genes, score_name=f'{CONTINUUM_PREFIX}_{role}_trajectory_score',
                      use_raw=False, random_state=RANDOM_STATE)

adata_total.obs[f'{CONTINUUM_PREFIX}_marker_axis'] = (
    adata_total.obs[f'{CONTINUUM_PREFIX}_late_trajectory_score'].astype(float)
    - adata_total.obs[f'{CONTINUUM_PREFIX}_early_trajectory_score'].astype(float))
# Strip label for the TOTAL heatmap: the fine segment call.
adata_total.obs['total_segment_marker_call'] = adata_total.obs['segment_class'].astype(str).values

# %%
# --- Section 2f: diffmap + PT-anchored root + Scanpy DPT on the PASS-2 embedding ---
marker_axis = adata_total.obs[f'{CONTINUUM_PREFIX}_marker_axis'].to_numpy(dtype=float)
use_rep = 'X_harmony'
n_neighbors = max(N_NEIGHBORS, int(np.sqrt(adata_total.n_obs)))
neighbors_key = 'trajectory_neighbors'

sc.pp.neighbors(adata_total, n_neighbors=n_neighbors, use_rep=use_rep,
                key_added=neighbors_key, random_state=RANDOM_STATE)

# %%
# --- Section 2f.1: PAGA topology check -- IS there a continuum to order? ---
# DPT will happily return a smooth axis over discrete cell types or over noise; no part of DPT
# tests whether a continuum exists. That judgement has to be made first, and PAGA connectivity
# between anatomically adjacent segments is how. Where a consecutive pair is disconnected, DPT
# is interpolating across a void, and any gene trend reported across that stretch is an artifact
# of the interpolation rather than a measured transition.
sc.tl.paga(adata_total, groups='segment_class', neighbors_key=neighbors_key)
sc.pl.paga(adata_total, threshold=PAGA_CONNECTIVITY_THRESHOLD, color='segment_class',
           frameon=False, show=False)
plt.savefig(CELLTYPING_DIR / 'paga_segment_topology.png', dpi=150, bbox_inches='tight')
plt.close()

_paga_conn = pd.DataFrame(
    np.asarray(adata_total.uns['paga']['connectivities'].todense()),
    index=adata_total.obs['segment_class'].cat.categories,
    columns=adata_total.obs['segment_class'].cat.categories)
_paga_conn.to_csv(CELLTYPING_DIR / 'paga_connectivity_matrix.csv')

# Connectivity between segments that are ADJACENT in the anatomy -- the pairs the pseudospace
# axis claims to traverse.
_present_order = [s for s in NEPHRON_SEGMENT_ORDER if s in _paga_conn.index]
consecutive_paga = pd.DataFrame([
    {'from': a, 'to': b, 'paga_connectivity': float(_paga_conn.loc[a, b]),
     'connected': bool(_paga_conn.loc[a, b] >= PAGA_CONNECTIVITY_THRESHOLD)}
    for a, b in zip(_present_order[:-1], _present_order[1:])])
consecutive_paga.to_csv(CELLTYPING_DIR / 'paga_consecutive_segment_connectivity.csv', index=False)
display(consecutive_paga.round(4))

_gaps = consecutive_paga[~consecutive_paga['connected']]
if len(_gaps):
    print(f'\nWARNING: {len(_gaps)} anatomically-consecutive segment pair(s) fall below the PAGA '
          f'connectivity threshold ({PAGA_CONNECTIVITY_THRESHOLD}):')
    for _, r in _gaps.iterrows():
        print(f"    {r['from']:>14s} -> {r['to']:<14s}  connectivity={r['paga_connectivity']:.4f}")
    print('  DPT INTERPOLATES across these boundaries rather than measuring a transition. Treat '
          '  gene trends spanning them as interpolation artifacts, and say so in the Methods.\n'
          '  Expected offenders: PT->DTL (sharp corticomedullary epithelial boundary) and the\n'
          '  nephron->collecting-duct junction (distinct developmental lineages: metanephric\n'
          '  mesenchyme vs ureteric bud).')
else:
    print(f'\nAll {len(consecutive_paga)} consecutive segment pairs exceed the connectivity '
          'threshold: the anatomical ordering is supported by a connected graph.')

# %%
sc.tl.diffmap(adata_total, neighbors_key=neighbors_key, random_state=RANDOM_STATE)
# DIAGNOSTIC ONLY. This reports which diffusion component tracks the marker axis best; it does
# NOT feed sc.tl.dpt below, which uses scanpy's default n_dcs across all components. To let the
# selected component drive the pseudotime instead, pass n_dcs=diagnostic_component + 1 to
# sc.tl.dpt -- that changes total_scanpy_dpt, so it is not done implicitly.
diagnostic_component, component_df = choose_diffusion_component(
    adata_total, marker_axis, n_components_to_test=N_DIFFMAP_COMPONENTS_TO_TEST,
    eigenvalue_floor=EIGENVALUE_FLOOR, min_valid=MIN_VALID_FOR_SPEARMAN)
print(f'Diffusion component best correlated with the marker axis: {diagnostic_component} '
      '(diagnostic; sc.tl.dpt below uses all components)')
display(component_df)

# Root at the medoid of the lowest-marker-axis PT Leiden cluster.
iroot, root_cluster = choose_root_pt_cluster(
    adata_total, marker_axis, leiden_col='leiden_coarse', rep_key=use_rep,
    family_col='broad_tubule_marker_call', pt_label='PT')
if iroot is None:
    iroot = choose_root_global(marker_axis, adata_total, use_rep, bottom_quantile=0.01)
    root_cluster = 'fallback (global bottom-1%)'
adata_total.uns['iroot'] = iroot
print(f'Root: {adata_total.obs_names[iroot]} (leiden_coarse={root_cluster}, '
      f'segment={adata_total.obs["segment_class"].iloc[iroot]}), marker_axis={marker_axis[iroot]:.3f}')

sc.tl.dpt(adata_total, neighbors_key=neighbors_key)
total_scanpy_dpt = orient_and_normalize(
    adata_total.obs['dpt_pseudotime'].to_numpy(dtype=float), marker_axis,
    min_valid=MIN_VALID_FOR_SPEARMAN)
adata_total.obs['total_scanpy_dpt'] = total_scanpy_dpt

# sc.tl.dpt assigns inf to cells in graph components disconnected from the root, and
# orient_and_normalize maps every non-finite value to NaN. Drop those cells from the continuum
# here: total_scanpy_dpt is contractually a finite [0, 1] coordinate for Sections 3-5, and a NaN
# reaching Section 4's common-support percentiles would silently poison the whole grid.
_nonfinite = ~np.isfinite(total_scanpy_dpt)
if _nonfinite.any():
    print(f'WARNING: {int(_nonfinite.sum()):,} tubules ({_nonfinite.mean():.2%}) have a '
          'non-finite DPT (disconnected from the root component) and are dropped from the '
          'continuum. A large fraction means the trajectory graph is fragmented -- raise '
          f'n_neighbors (currently {n_neighbors}) and re-run.')
    adata_total = adata_total[~_nonfinite].copy()
    adata_total.uns.pop('iroot', None)          # positional index is stale after subsetting
    marker_axis = marker_axis[~_nonfinite]
    total_scanpy_dpt = total_scanpy_dpt[~_nonfinite]

print('Spearman(total_scanpy_dpt, marker_axis) =',
      round(spearmanr(total_scanpy_dpt, marker_axis, nan_policy='omit').correlation, 3))

sc.pl.embedding(adata_total, basis='umap', color=['total_scanpy_dpt', 'segment_class'],
                cmap='viridis', frameon=False, ncols=2, show=False)
plt.savefig(CELLTYPING_DIR / 'total_dpt_umap.png', dpi=150, bbox_inches='tight')
plt.show()

# %%
# Does the DPT order the segments the way the anatomy says it should?
segment_dpt = (adata_total.obs.groupby('segment_class', observed=True)['total_scanpy_dpt']
               .agg(['size', 'mean', 'median']).reindex(
                   [s for s in SEGMENT_DISPLAY_ORDER if s in set(adata_total.obs['segment_class'])]))
segment_dpt.to_csv(CELLTYPING_DIR / 'dpt_by_segment.csv')
display(segment_dpt.round(3))

fig, ax = plt.subplots(figsize=(9, 4))
order = list(segment_dpt.index)
ax.boxplot([adata_total.obs.loc[adata_total.obs['segment_class'].astype(str) == seg,
                                'total_scanpy_dpt'].to_numpy(dtype=float) for seg in order],
           showfliers=False)
# Set the tick labels separately: boxplot's `labels=` was renamed `tick_labels` in
# matplotlib 3.9 and removed in 3.11.
ax.set_xticks(range(1, len(order) + 1))
ax.set_xticklabels(order)
ax.set_ylabel('total_scanpy_dpt')
ax.set_title('Pseudospace by segment (anatomical order left to right)')
plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
fig.savefig(CELLTYPING_DIR / 'dpt_by_segment.png', dpi=150, bbox_inches='tight')
plt.show()

# %%
# --- Section 2g: save the canonical DPT file ---
DPT_EXTRA_OBS_COLS = ['segment_class', 'coarse_class', 'leiden_coarse',
                      'celltype_primary_segment', 'celltype_secondary_segment',
                      f'{CONTINUUM_PREFIX}_marker_axis',
                      # The two halves of the marker axis are carried separately so Section 4
                      # can test whether AKI shifts the axis's own BASIS (de-differentiation)
                      # rather than only shifting position along it.
                      f'{CONTINUUM_PREFIX}_early_trajectory_score',
                      f'{CONTINUUM_PREFIX}_late_trajectory_score']

adata_dpt_saved = save_total_pseudotime_anndata(
    adata_total, 'total_scanpy_dpt', DPT_OUTPUT_PATH,
    expression_filtered_adata=adata_all_expression_filtered,
    annotation_adata=adata_total, project_dir=PROJECT_DIR,
    forbidden_labels=REMOVE_CLASSES + ['Glomeruli', 'nan', 'none'],
    assert_no_nan=True,
    extra_obs_cols=DPT_EXTRA_OBS_COLS)

_positions = adata_total.obs_names.get_indexer(adata_dpt_saved.obs_names)
if (_positions < 0).any():
    raise RuntimeError('DPT file contains cells missing from adata_total; cannot map the embedding.')

# save_total_pseudotime_anndata's `extra_obs_cols` only KEEPS columns that are already on the
# expression snapshot -- and that snapshot was taken in Section 2a, before any label existed. So
# # copy the label columns across from adata_total explicitly.
for _col in DPT_EXTRA_OBS_COLS:
    if _col not in adata_total.obs.columns or _col in adata_dpt_saved.obs.columns:
        continue
    _values = adata_total.obs[_col].to_numpy()[_positions]
    adata_dpt_saved.obs[_col] = (_values.astype(float)
                                 if pd.api.types.is_numeric_dtype(adata_total.obs[_col])
                                 else _values.astype(str))

# Carry the PASS-2 embedding onto the DPT file. Without this the file inherits the pass-1
# coordinates from adata_all_expression_filtered, and the subset DPTs below would be built on a
# different embedding than the global DPT.
for _key in ['X_harmony', 'X_pca', 'X_umap']:
    if _key in adata_total.obsm:
        adata_dpt_saved.obsm[_key] = np.asarray(adata_total.obsm[_key])[_positions]
adata_dpt_saved.write(DPT_OUTPUT_PATH)
print(f'DPT obs columns: {list(adata_dpt_saved.obs.columns)}')

print('\n--- Section 2 summary ---')
print(f'Labels: manual, on the pass-1 embedding (no re-typing after pass-2). Root: {root_cluster}')
print(f'Tubule continuum: {adata_total.n_obs:,} tubules')
print(adata_total.obs['coarse_class'].value_counts().reindex(COARSE_ORDER).to_string())
print(f'total_scanpy_dpt range: [{np.nanmin(total_scanpy_dpt):.3f}, {np.nanmax(total_scanpy_dpt):.3f}]')

# %% [markdown]
# # Section 3 - marker heatmaps
#
# Curated panels on the global pseudospace, then on per-family recomputed pseudospace. Only
# `mouse_control` and `mouse_aki` subsets are produced (no human samples in this pipeline).

# %%
OUTPUT_DIR = HEATMAP_OUTPUT_DIR                     # plot_marker_heatmap saves here

adata_heatmap = sc.read_h5ad(DPT_OUTPUT_PATH)
if 'lognorm' in adata_heatmap.layers:
    adata_heatmap.X = adata_heatmap.layers['lognorm'].copy()

_leaked = sorted(set(adata_heatmap.obs['coarse_class'].astype(str)) & set(REMOVE_CLASSES))
if _leaked:
    raise AssertionError(
        f'{DPT_OUTPUT_PATH.name} contains non-tubule classes {_leaked}; re-run Section 2.')

BASE_OBS_COLUMNS = ['sample', 'species', 'condition', 'segment_class', 'coarse_class',
                    'broad_tubule_marker_call', 'total_segment_marker_call', 'total_scanpy_dpt']
missing_obs = [col for col in BASE_OBS_COLUMNS if col not in adata_heatmap.obs.columns]
if missing_obs:
    raise KeyError(f'{DPT_OUTPUT_PATH.name} is missing required obs columns: {missing_obs}')

mouse_control_mask = adata_heatmap.obs['condition'].astype(str).str.lower().eq('control').to_numpy()
mouse_aki_mask = adata_heatmap.obs['condition'].astype(str).str.lower().eq('ir').to_numpy()
CONDITIONS = [('Mouse_Control', 'Mouse Control', mouse_control_mask),
              ('Mouse_AKI', 'Mouse AKI', mouse_aki_mask)]

print(f'Loaded {DPT_OUTPUT_PATH.name}: {adata_heatmap.n_obs:,} tubules x {adata_heatmap.n_vars:,} genes')
print(f'mouse_control: {int(mouse_control_mask.sum()):,} tubules')
print(f'mouse_aki:     {int(mouse_aki_mask.sum()):,} tubules')
print(f'Heatmaps will be saved to: {HEATMAP_OUTPUT_DIR.relative_to(PROJECT_DIR)}')

# %%
# TOTAL panel: one row per continuous-axis marker, colored by the segment named in its label.
TOTAL_MARKER_GROUPS = {
    row['label']: [{'label': row['label'], 'candidates': row['candidates']}]
    for row in TOTAL_NEPHRON_MARKERS
}
TOTAL_GROUP_ORDER = [row['label'] for row in TOTAL_NEPHRON_MARKERS]


def _segment_from_label(label):
    """'Sycn (mTAL)' -> 'mTAL'; 'Trpv5 (DCT2-enriched)' -> 'DCT2'; 'Nphs1 (Glom)' -> 'Podocyte'."""
    match = re.search(r'\(([^);]+)', label)
    segment = match.group(1).strip() if match else label
    segment = re.sub(r'-enriched$', '', segment)
    return {'Glom': 'Podocyte'}.get(segment, segment)


TOTAL_GROUP_COLORS = {
    group: SEGMENT_STRIP_COLORS.get(_segment_from_label(group), '#4D4D4D')
    for group in TOTAL_GROUP_ORDER
}
STRIP_ORDER = [s for s in SEGMENT_DISPLAY_ORDER
               if s in set(adata_heatmap.obs['total_segment_marker_call'].astype(str))]

# %%
# Marker resolution check -- every panel must resolve before any plotting runs.
panels_to_check = {'TOTAL': (TOTAL_MARKER_GROUPS, TOTAL_GROUP_ORDER)}
panels_to_check.update({name: (spec['markers'], spec['order']) for name, spec in HEATMAP_PANELS.items()})

any_missing = False
for panel_name, (marker_groups, group_order) in panels_to_check.items():
    _, _, _, missing = resolve_available_marker_groups(adata_heatmap, marker_groups, group_order)
    n_total = sum(len(v) for v in marker_groups.values())
    status = f'{panel_name}: {n_total - len(missing)}/{n_total} markers resolved'
    if missing:
        status += f'; MISSING: {missing}'
        any_missing = True
    print(status)

if any_missing:
    warnings.warn('One or more marker panels are missing genes in the mouse-only gene space.')
else:
    print('\nAll marker panels resolved with no missing genes.')

# %%
# --- Heatmaps on the global nephron coordinate ---
total_heatmap_results = {}
for condition_key, condition_label, condition_mask in CONDITIONS:
    total_heatmap_results[f'{condition_key}_TOTAL'] = plot_marker_heatmap(
        adata_heatmap, TOTAL_MARKER_GROUPS, TOTAL_GROUP_ORDER, TOTAL_GROUP_COLORS,
        output_dir=OUTPUT_DIR, project_dir=PROJECT_DIR,
        pseudotime_col='total_scanpy_dpt', cluster_col='coarse_class', cluster_value=None,
        mask=condition_mask, target_label=condition_label,
        strip_col='total_segment_marker_call', strip_order=STRIP_ORDER,
        strip_colors=SEGMENT_STRIP_COLORS,
        title=f'{condition_label} nephron markers across total Scanpy DPT',
        output_name=f'{condition_key.lower()}_total_marker_heatmap.png',
        n_bins=N_BINS, smooth_sigma=SMOOTH_SIGMA, trim_fraction=TRIM_FRACTION, z_clip=3,
        show_module_legend=False, show_group_labels=False, figsize=(18, 11),
    )

# Unsmoothed companion -- guards against smoothing-manufactured monotonicity.
plot_marker_heatmap(
    adata_heatmap, TOTAL_MARKER_GROUPS, TOTAL_GROUP_ORDER, TOTAL_GROUP_COLORS,
    output_dir=OUTPUT_DIR, project_dir=PROJECT_DIR,
    pseudotime_col='total_scanpy_dpt', cluster_col='coarse_class', cluster_value=None,
    mask=mouse_control_mask, target_label='Mouse Control (unsmoothed)',
    strip_col='total_segment_marker_call', strip_order=STRIP_ORDER,
    strip_colors=SEGMENT_STRIP_COLORS,
    title='Mouse Control nephron markers across total DPT (UNSMOOTHED)',
    output_name='mouse_control_total_marker_heatmap_unsmoothed.png',
    n_bins=N_BINS, smooth_sigma=0, trim_fraction=TRIM_FRACTION, z_clip=3,
    show_module_legend=False, show_group_labels=False, figsize=(18, 11))

total_summary = pd.DataFrame(total_heatmap_results).T
total_summary.index.name = 'heatmap'
total_summary.to_csv(OUTPUT_DIR / 'total_nephron_heatmap_summary.csv')
display(total_summary[['n_cells_before_trim', 'n_cells_trimmed_per_end', 'n_cells', 'cluster_col',
                       'cluster_value', 'pseudotime_col', 'dpt_min', 'dpt_max', 'output_path']])

# %%
# --- Each curated panel, on its own family, along the GLOBAL pseudospace ---
coarse_labels_str = adata_heatmap.obs['coarse_class'].astype(str)
global_segment_results = {}

for condition_key, condition_label, condition_mask in CONDITIONS:
    for panel_name, spec in HEATMAP_PANELS.items():
        class_mask = coarse_labels_str.isin(spec['classes']).to_numpy()
        try:
            global_segment_results[f'{condition_key}_{panel_name}'] = plot_marker_heatmap(
                adata_heatmap, spec['markers'], spec['order'], spec['colors'],
                output_dir=OUTPUT_DIR, project_dir=PROJECT_DIR,
                pseudotime_col='total_scanpy_dpt', cluster_col='coarse_class', cluster_value=None,
                mask=condition_mask & class_mask, target_label=condition_label,
                strip_col='total_segment_marker_call', strip_order=STRIP_ORDER,
                strip_colors=SEGMENT_STRIP_COLORS,
                title=f'{condition_label} {panel_name} markers on global nephron Scanpy DPT '
                      f'({"+".join(spec["classes"])})',
                output_name=f'{condition_key.lower()}_{panel_name.lower()}_global_dpt_marker_heatmap.png',
                n_bins=N_BINS, smooth_sigma=SMOOTH_SIGMA, trim_fraction=TRIM_FRACTION,
                z_clip=spec['z_clip'],
            )
        except (ValueError, KeyError) as exc:
            print(f'SKIPPED {condition_key}/{panel_name}: {exc}')

global_segment_summary = pd.DataFrame(global_segment_results).T
global_segment_summary.index.name = 'heatmap'
global_segment_summary.to_csv(OUTPUT_DIR / 'global_segment_heatmap_summary.csv')
display(global_segment_summary[['n_cells_before_trim', 'n_cells_trimmed_per_end', 'n_cells',
                                'pseudotime_col', 'dpt_min', 'dpt_max', 'output_path']])

# %% [markdown]
# ## Family-specific diffusion pseudotime
#
# Diffusion pseudotime rebuilt from scratch inside each coarse family, anchored by that family's
# marker panel. These columns are additive: the global `total_scanpy_dpt` is never overwritten.

# %%
subset_dpt_diagnostics = []
subset_dpt_available = {}
for family, spec in SUBSET_DPT_SPECS.items():
    panel = HEATMAP_PANELS[spec['panel']]
    try:
        subset_dpt_diagnostics.append(
            recompute_subset_dpt(
                adata_heatmap, segment=family, marker_groups=panel['markers'],
                group_order=panel['order'], output_col=spec['pseudotime_col'],
                n_neighbors=N_NEIGHBORS, random_state=RANDOM_STATE))
        subset_dpt_available[family] = spec
    except (ValueError, KeyError) as exc:
        print(f'SKIPPED subset DPT for {family}: {exc}')

subset_dpt_diagnostics = pd.DataFrame(subset_dpt_diagnostics)
subset_dpt_diagnostics.to_csv(OUTPUT_DIR / 'subset_recomputed_dpt_diagnostics.csv', index=False)
display(subset_dpt_diagnostics)

# %%
# Subset DPT on the embedding, and against the global coordinate.
for basis, obsm_key in [('umap', 'X_umap'), ('pca', 'X_pca')]:
    if obsm_key not in adata_heatmap.obsm or not subset_dpt_available:
        continue
    coords = np.asarray(adata_heatmap.obsm[obsm_key])
    fig, axes = plt.subplots(len(subset_dpt_available), 2,
                             figsize=(12, 4 * len(subset_dpt_available)),
                             constrained_layout=True, squeeze=False)
    for row_i, (family, spec) in enumerate(subset_dpt_available.items()):
        subset_mask = coarse_labels_str.eq(family).to_numpy()
        subset_dpt = adata_heatmap.obs.loc[subset_mask, spec['pseudotime_col']].to_numpy(dtype=float)
        total_dpt = adata_heatmap.obs.loc[subset_mask, 'total_scanpy_dpt'].to_numpy(dtype=float)

        axes[row_i, 0].scatter(coords[:, 0], coords[:, 1], s=1, color='#E0E0E0', alpha=0.2,
                               linewidths=0, rasterized=True)
        scatter = axes[row_i, 0].scatter(coords[subset_mask, 0], coords[subset_mask, 1],
                                         c=subset_dpt, cmap='viridis', vmin=0, vmax=1, s=4,
                                         linewidths=0, rasterized=True)
        axes[row_i, 0].set_title(f'{family}: subset-recomputed DPT')
        axes[row_i, 0].set_xticks([])
        axes[row_i, 0].set_yticks([])
        fig.colorbar(scatter, ax=axes[row_i, 0], label=spec['pseudotime_col'])

        finite = np.isfinite(total_dpt) & np.isfinite(subset_dpt)
        axes[row_i, 1].hexbin(total_dpt[finite], subset_dpt[finite], gridsize=55, mincnt=1, cmap='magma')
        corr = spearmanr(total_dpt[finite], subset_dpt[finite]).correlation
        axes[row_i, 1].set_title(f'{family}: total vs subset DPT (Spearman={corr:.2f})')
        axes[row_i, 1].set_xlabel('total_scanpy_dpt')
        axes[row_i, 1].set_ylabel(spec['pseudotime_col'])

    fig.savefig(OUTPUT_DIR / f'subset_recomputed_dpt_diagnostics_{basis}.png', dpi=220,
                bbox_inches='tight')
    plt.show()

# %%
# --- Each curated panel on its family's OWN recomputed pseudospace ---
subset_dpt_results = {}
for condition_key, condition_label, condition_mask in CONDITIONS:
    for family, spec in subset_dpt_available.items():
        panel = HEATMAP_PANELS[spec['panel']]
        try:
            subset_dpt_results[f'{condition_key}_{family}'] = plot_marker_heatmap(
                adata_heatmap, panel['markers'], panel['order'], panel['colors'],
                output_dir=OUTPUT_DIR, project_dir=PROJECT_DIR,
                pseudotime_col=spec['pseudotime_col'], cluster_col='coarse_class', cluster_value=family,
                mask=condition_mask, target_label=condition_label,
                strip_col='total_segment_marker_call', strip_order=STRIP_ORDER,
                strip_colors=SEGMENT_STRIP_COLORS,
                title=f'{condition_label} {family} markers across subset-recomputed Scanpy DPT',
                output_name=f'{condition_key.lower()}_{family.lower()}_subset_recomputed_dpt_marker_heatmap.png',
                n_bins=N_BINS, smooth_sigma=SMOOTH_SIGMA, trim_fraction=TRIM_FRACTION,
                z_clip=panel['z_clip'])
        except (ValueError, KeyError) as exc:
            print(f'SKIPPED {condition_key}/{family} subset heatmap: {exc}')

subset_dpt_summary = pd.DataFrame(subset_dpt_results).T
subset_dpt_summary.index.name = 'heatmap'
subset_dpt_summary.to_csv(OUTPUT_DIR / 'subset_recomputed_dpt_heatmap_summary.csv')
display(subset_dpt_summary[['n_cells_before_trim', 'n_cells_trimmed_per_end', 'n_cells',
                            'cluster_value', 'pseudotime_col', 'dpt_min', 'dpt_max', 'output_path']])

# %%
# Persist the subset-recomputed DPT columns into the DPT file (in-memory only until now).
_subset_cols = [spec['pseudotime_col'] for spec in subset_dpt_available.values()
                if spec['pseudotime_col'] in adata_heatmap.obs]
if _subset_cols:
    _dpt = sc.read_h5ad(DPT_OUTPUT_PATH)
    for col in _subset_cols:
        _dpt.obs[col] = adata_heatmap.obs.loc[_dpt.obs_names, col].to_numpy(dtype=float)
    _dpt.write(DPT_OUTPUT_PATH)
    print(f'Wrote subset-DPT columns {_subset_cols} into {DPT_OUTPUT_PATH.name}')

print('\n--- Section 3 summary ---')
print(f'Mouse Control tubules: {int(mouse_control_mask.sum()):,}')
print(f'Mouse AKI tubules:     {int(mouse_aki_mask.sum()):,}')
print(f'total_scanpy_dpt range: [{adata_heatmap.obs["total_scanpy_dpt"].min():.3f}, '
      f'{adata_heatmap.obs["total_scanpy_dpt"].max():.3f}]')
for family, spec in subset_dpt_available.items():
    vals = adata_heatmap.obs[spec['pseudotime_col']].dropna()
    print(f'{spec["pseudotime_col"]} range: [{vals.min():.3f}, {vals.max():.3f}] (n={len(vals):,})')
print(f'Wrote {len(list(HEATMAP_OUTPUT_DIR.glob("*.png")))} PNGs and '
      f'{len(list(HEATMAP_OUTPUT_DIR.glob("*.csv")))} CSVs to '
      f'{HEATMAP_OUTPUT_DIR.relative_to(PROJECT_DIR)}')

# %% [markdown]
# # Section 4 - healthy vs AKI (mouse-only PT tubules)
#
# Within the PT cohort, each gene and each pathway module is decomposed along the shared
# pseudospace (`total_scanpy_dpt`, already in [0, 1]) into three nested penalized-B-spline GAMs
# on all PT tubules, with condition indicator c (Control -> healthy=0, IR -> aki=1):
#
# * M0 (shared): intercept + f(s)
# * M1 (level):  intercept + f(s) + beta*c
# * M2 (shape):  intercept + f(s) + beta*c + g(s)*c   (condition x basis interaction)
#
# Smoothing lambda is GCV-selected per gene/module on M2 and reused for M0/M1. Nested F-tests
# give level_F (M1 vs M0) and shape_F (M2 vs M1). Per-condition curves come from the joint M2
# model on the common-support grid (intersection of the two conditions' [p1, p99]), so shape is
# never compared where one condition has no cells. Primary ranking metric is shape_rms.
#
# Significance is honest for the 2-vs-2 specimen design: the healthy/aki label is permuted
# across the 4 samples (all C(4,2)=6 relabelings). Pseudospace is never shuffled within a
# sample.
#
# > **Cohort disclaimer.** 2 control + 2 IR specimens. Results are descriptive and
# > effect-size-ranked, not confirmatory; the exact permutation p-floor is 2/6 (a feature and its condition-swapped mirror tie), so `sample_perm_p` is
# > a calibration/sanity check only.

# %%
# NB: bh_adjust is deliberately NOT imported. With 2v2 specimens sample_perm_p is floored at
# 2/6 and takes ~3 distinct values, so BH over it is meaningless; results are ranked by the
# shape_rms effect size instead (see Section 4.7).
from pseudospace.stats_gam import (
    as_csr,
    zscore_rows,
    gam_internal_knots,
    safe_spearman,
    resolve_present,
)

# ============================================================================
# Section 4.0 -- imports + GAM/stats utilities (copied from notebook 5) + output dir
# ============================================================================
# numpy/pandas/scanpy/pyplot, spearmanr, ks_2samp, wasserstein_distance and gaussian_filter1d
# all come from the Section 0 / Section 1 import cells; display / Image come from the
# fallback-aware import in Section 0.
import json

HEALTHY_VS_AKI_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# --- utilities copied verbatim from 5_mouse_healthy_vs_aki.ipynb (CONFIG defaults inlined) ---

# %%
# ----------------------------------------------------------------------------
# Section 4 configuration
# ----------------------------------------------------------------------------
SECTION4_CONFIG = {
    'n_internal_knots': 9,                     # target ~8-10 (vs old default of 3) -> less oversmoothing
    'lambda_grid': np.logspace(-3, 3, 13),     # GCV smoothing search grid (alpha in make_gam_penalty)
    'common_support_pct': (1, 99),             # per-condition [p1, p99]; grid = intersection
    'grid_points': 101,
    'min_detected_fraction': 0.02,             # tested-gene background: detected in >= 2% of PT cells
    'min_mean_expression': 0.0,
    'pathway_min_genes': 10,
    'pathway_max_genes': 150,
    'pathway_libraries': ['Reactome_2022', 'MSigDB_Hallmark_2020', 'KEGG_2019_Mouse'],
    'ranking_metric': 'shape_rms',
    'random_state': RANDOM_STATE,
    # visualization layer
    'N_TOP_PLOT': 12,
    'heatmap_n_bins': N_BINS,                  # reuse NB3 (120)
    'heatmap_trim_fraction': TRIM_FRACTION,    # reuse NB3 (0.05)
    'heatmap_smooth_sigma': SMOOTH_SIGMA,      # reuse NB3 (2.5)
    'diverging_cmap': 'bwr',
    'condition_colors': {'healthy': 'tab:blue', 'aki': 'tab:orange'},
}
MARKER_GROUPS = {
    'S1': ['Slc5a2', 'Slc5a12', 'Gatm'],       # Gatm moved here (early PT, not S2)
    'S2': ['Slc22a6', 'Slc13a3', 'Cyp2e1'],    # dropped Slc34a1 (S1/S2-shared) + Gatm
    'S3': ['Slc22a7', 'Slc5a10', 'Slc7a13'],   # sex-neutral anchors + Slc7a13 (male-biased); dropped Acsm3
}
for key, val in SECTION4_CONFIG.items():
    print(f'{key:>22}: {val}')

# %%
# ----------------------------------------------------------------------------
# Section 4.0c -- shared plotting helpers
# Defined up here because Sections 4.6b, 4.8 and 4.9 all use them.
# ----------------------------------------------------------------------------
N_TOP = SECTION4_CONFIG['N_TOP_PLOT']
COL_H = SECTION4_CONFIG['condition_colors']['healthy']
COL_A = SECTION4_CONFIG['condition_colors']['aki']
figure_index = []


def _grid_axes(n_panels, ncols=3, panel_h=3.0, panel_w=5.0):
    nrows = max(1, int(np.ceil(n_panels / ncols)))
    fig, axes = plt.subplots(nrows, ncols, figsize=(panel_w * ncols, panel_h * nrows), squeeze=False)
    return fig, axes.ravel()


def _save(fig, name):
    path = HEALTHY_VS_AKI_OUTPUT_DIR / name
    fig.tight_layout()
    fig.savefig(path, bbox_inches='tight', dpi=130)
    plt.close(fig)
    figure_index.append(name)
    display(Image(filename=str(path)))


def binned_matrix_profile(mat, s_vals, lo_, hi_, n_bins, sigma):
    """mat: (cells, cols) dense; returns (cols, n_bins) interpolated + smoothed binned means."""
    edges = np.linspace(lo_, hi_, n_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    idx = np.clip(np.digitize(s_vals, edges) - 1, 0, n_bins - 1)
    inb = (s_vals >= lo_) & (s_vals <= hi_)
    out = np.full((mat.shape[1], n_bins), np.nan)
    for b in range(n_bins):
        m = inb & (idx == b)
        if m.any():
            out[:, b] = mat[m].mean(axis=0)
    for i in range(out.shape[0]):
        row = out[i]
        finite = np.isfinite(row)
        if not finite.any():
            continue
        if not finite.all():
            row = np.interp(centers, centers[finite], row[finite])
        out[i] = gaussian_filter1d(row, sigma)
    return out, centers


# %%
from pseudospace.pathways import (
    build_pathway_membership,
    member_gene_evidence,
    summarize_pathway_redundancy,
)
from pseudospace.levelshape import (
    run_level_shape,
    sample_perm_pvalues,
    summarize_curve_effects,
)

# %% [markdown]
# ## 4.1 - PT cohort, shared coordinate, common-support grid, tested genes

# %%
# ----------------------------------------------------------------------------
# Section 4.1 -- PT cohort, shared coordinate, common-support grid, tested genes
# ----------------------------------------------------------------------------
adata_full = sc.read_h5ad(DPT_OUTPUT_PATH)
if 'lognorm' in adata_full.layers:
    adata_full.X = adata_full.layers['lognorm'].copy()

call = adata_full.obs['broad_tubule_marker_call'].astype(str)
_unexpected = sorted(set(call) - set(KEEP_TUBULE_CLASSES))
if _unexpected:
    raise ValueError(
        f'{DPT_OUTPUT_PATH.name} contains non-tubule classes {_unexpected}; this file must hold '
        f'only {KEEP_TUBULE_CLASSES}. Re-run Section 2, which filters these out before saving '
        'the DPT AnnData.'
    )

# PT subset: cluster-derived broad_tubule_marker_call == 'PT'.
adata_pt = adata_full[call == 'PT'].copy()
if adata_pt.n_obs == 0:
    raise ValueError("No cells labelled 'PT' in the DPT file; check COARSE_LABELS in Section 2.")

# Condition: Control -> healthy (0), IR -> aki (1); matched case-insensitively.
cond_lower = adata_pt.obs['condition'].astype(str).str.lower()
c = (cond_lower == 'ir').astype(float).to_numpy()
samples = adata_pt.obs['sample'].astype(str).to_numpy()
control_samples = sorted(set(samples[c == 0]))
aki_samples = sorted(set(samples[c == 1]))

# Shared pseudospace = total_scanpy_dpt (already [0,1]); do NOT rerank per condition.
s = adata_pt.obs['total_scanpy_dpt'].to_numpy(dtype=float)
if not np.all(np.isfinite(s) & (s >= 0) & (s <= 1)):
    raise ValueError(
        f'total_scanpy_dpt must be a finite [0, 1] coordinate, but '
        f'{int((~np.isfinite(s)).sum()):,} of {s.size:,} PT values are non-finite and '
        f'{int((np.isfinite(s) & ((s < 0) | (s > 1))).sum()):,} fall outside [0, 1]. '
        'Section 2f drops non-finite DPT cells -- regenerate the DPT file.'
    )

# Common-support grid = intersection of per-condition [p1, p99].
p_lo, p_hi = SECTION4_CONFIG['common_support_pct']
lo = max(np.percentile(s[c == 0], p_lo), np.percentile(s[c == 1], p_lo))
hi = min(np.percentile(s[c == 0], p_hi), np.percentile(s[c == 1], p_hi))
if not lo < hi:
    raise ValueError(
        f'Empty common support between conditions: healthy and AKI [p{p_lo}, p{p_hi}] ranges '
        f'give lo={lo:.4f} >= hi={hi:.4f}. The two conditions occupy disjoint stretches of '
        'pseudospace, so no shape comparison is possible.'
    )
grid = np.linspace(lo, hi, SECTION4_CONFIG['grid_points'])
knots = gam_internal_knots(s, basis_df=3 + SECTION4_CONFIG['n_internal_knots'])

# Tested-gene background: detection filter over PT cells.
Y_all = as_csr(adata_pt.layers['lognorm'])
detected = np.asarray((Y_all > 0).sum(axis=0)).ravel()
gene_mean_all = np.asarray(Y_all.mean(axis=0)).ravel()
min_detected = int(np.ceil(SECTION4_CONFIG['min_detected_fraction'] * adata_pt.n_obs))
tested = (detected >= min_detected) & (gene_mean_all >= SECTION4_CONFIG['min_mean_expression'])
gene_names = adata_pt.var_names.to_numpy()[tested]
gene_local = {g_: i for i, g_ in enumerate(gene_names)}
Y_genes = Y_all[:, tested].tocsr()
Y_csc = Y_genes.tocsc()

# Per-gene mean/std over PT cells (for module z-scores + z-score plotting).
gene_mean = np.asarray(Y_genes.mean(axis=0)).ravel()
gene_sq = np.asarray(Y_genes.multiply(Y_genes).mean(axis=0)).ravel()
gene_std = np.sqrt(np.maximum(gene_sq - gene_mean ** 2, 0.0))
gene_std[gene_std == 0] = 1.0

# %%
# ----------------------------------------------------------------------------
# Section 4.1b -- injury scores on the PT cohort (diagnostic; never feeds the axis)
# ----------------------------------------------------------------------------
INJURY_PANELS = {'fr_ptc': FR_PTC_MARKERS, 'acute_injury': ACUTE_INJURY_MARKERS}
injury_scores = {}
injury_panel_report = []
for _panel, _genes in INJURY_PANELS.items():
    _present = resolve_present(_genes, adata_pt.var_names.to_numpy())
    injury_panel_report.append({'panel': _panel, 'requested': '; '.join(_genes),
                                'present': '; '.join(_present), 'n_present': len(_present)})
    if not _present:
        print(f'WARNING: no {_panel} markers ({_genes}) present; skipping that score.')
        continue
    sc.tl.score_genes(adata_pt, _present, score_name=f'{_panel}_score', use_raw=False,
                      random_state=RANDOM_STATE)
    injury_scores[_panel] = adata_pt.obs[f'{_panel}_score'].to_numpy(dtype=float)
display(pd.DataFrame(injury_panel_report))

n_healthy = int((c == 0).sum())
n_aki = int((c == 1).sum())
print(f'PT tubules: {adata_pt.n_obs:,}  (healthy={n_healthy:,}, aki={n_aki:,})')
print(f'Healthy samples: {control_samples}    AKI samples: {aki_samples}')
print(f'Tested genes: {len(gene_names):,} / {adata_pt.n_vars:,}')
print(f'Common support: [{lo:.3f}, {hi:.3f}]   internal knots: {len(knots)}  ({knots.round(3)})')

# %% [markdown]
# ## 4.3 - gene-level nested level/shape decomposition + effect sizes

# %%
# ----------------------------------------------------------------------------
# Section 4.3 -- gene-level nested level/shape decomposition + effect sizes
# ----------------------------------------------------------------------------
gene_ls = run_level_shape(Y_genes, s, c, knots, grid, SECTION4_CONFIG['lambda_grid'])
print(f'Gene-level level/shape fit complete for {len(gene_names):,} genes.')
print(f"  shape_rms:    [{gene_ls['shape_rms'].min():.4f}, {gene_ls['shape_rms'].max():.4f}]")
print(f"  |level_eff|:  [{np.abs(gene_ls['level_effect']).min():.4f}, {np.abs(gene_ls['level_effect']).max():.4f}]")

# %%
# ----------------------------------------------------------------------------
# Section 4.4 -- sample-level condition-label permutation (honest for n=2 vs 2)
# ----------------------------------------------------------------------------
gene_shape_p, gene_level_p, splits, true_idx = sample_perm_pvalues(
    Y_genes, s, samples, knots, grid, SECTION4_CONFIG['lambda_grid'],
    gene_ls['lam_idx'], gene_ls['sl2'], gene_ls['p_b'], control_samples)

print('=' * 78)
print('SAMPLE-LEVEL PERMUTATION -- CALIBRATION/SANITY CHECK, NOT A POWERED TEST.')
print('With 2 healthy vs 2 AKI specimens there are only C(4,2)=6 sample relabelings;')
print('p = #{splits with permuted shape_rms >= observed}/6 (inclusive tail; ties count, so the')
print('observed labeling is included and p is never 0). The 6 splits form 3 mirror-image')
print('(healthy/aki swap) pairs; shape_rms is symmetric under swap, so a maximal feature stops at')
print('2/6 -- RANK BY EFFECT SIZE (shape_rms), not by p.')
print('=' * 78)
print('distinct gene sample_perm_p values:', np.unique(np.round(gene_shape_p, 4)))

# %% [markdown]
# ## 4.5 - pathway module scores + level/shape GAM

# %%
# ----------------------------------------------------------------------------
# Section 4.5 -- pathway module scores + same level/shape GAM (dense path)
# NB: Reactome/Hallmark/KEGG are human-symbol Enrichr libraries; genes are matched to mouse
# var_names by case-insensitive symbol collision, against the full tested mouse gene background.
# Module score = mean of per-gene z-scored lognorm (z-scored over the PT subset only).
# ----------------------------------------------------------------------------
# Membership is built by the same helper the cross-species workflow uses, with the coverage stages
# recorded. This notebook is mouse-only, so no ortholog table is needed; the upper bound is no longer
# applied to membership (it removed whole Hallmark sets before they were scored) and the count of
# pathways such a bound would remove is reported instead.
pathway_coverage_frames = []
symbol_mapping_report = []
for library in SECTION4_CONFIG['pathway_libraries']:
    gene_sets = json.loads((PATHWAY_LIBRARY_DIR / f'{library}.json').read_text())
    membership = build_pathway_membership(
        gene_sets, adata_pt.var_names,
        library_name=library,
        min_genes=SECTION4_CONFIG['pathway_min_genes'],
        max_genes=None,
        tested=gene_names,
    )
    pathway_coverage_frames.append(membership)
    symbol_mapping_report.append({
        'library': library,
        'n_pathways': int(len(membership)),
        'n_retained': int(membership['retained'].sum()),
        'n_dropped_too_few_genes': int((~membership['retained']).sum()),
        'n_over_previous_cap': int(
            (membership['n_assayed'] > SECTION4_CONFIG['pathway_max_genes']).sum()),
        'unique_symbols_requested': int(membership['n_requested'].sum()),
        'unique_symbols_matched': int(membership['n_assayed'].sum()),
        'symbol_match_rate': float(membership['n_assayed'].sum()
                                   / max(membership['n_requested'].sum(), 1))})
pathway_coverage = pd.concat(pathway_coverage_frames, ignore_index=True)
pathway_coverage.to_csv(HEALTHY_VS_AKI_OUTPUT_DIR / 'pathway_membership_coverage.csv', index=False)
retained_pathways = pathway_coverage[pathway_coverage['retained']].reset_index(drop=True)
P = len(retained_pathways)
print(f'Pathways retained: {P:,} of {len(pathway_coverage):,}; a '
      f'{SECTION4_CONFIG["pathway_max_genes"]}-member upper bound would remove '
      f'{int((pathway_coverage["n_assayed"] > SECTION4_CONFIG["pathway_max_genes"]).sum()):,}')

# Symbols are matched to the TESTED PT gene background by case-insensitive collision, and a set
# member can fail to match for two quite different reasons:
#   (1) species -- Reactome_2022 and MSigDB_Hallmark_2020 are HUMAN-symbol libraries, so any gene
#       whose mouse ortholog is not spelled identically is lost. KEGG_2019_Mouse is already mouse,
#       so this does not apply to it.
#   (2) detection -- the gene is simply not expressed in >= min_detected_fraction of PT cells and
#       never entered `gene_names` in the first place. This applies to EVERY library.
# The rate below is the combined effect and cannot separate the two; it is a floor on coverage,
# not an ortholog-mapping statistic. Reported because the loss used to be silent, and a low rate
# means the module score describes a biased subset of the set rather than the set.
symbol_mapping = pd.DataFrame(symbol_mapping_report)
symbol_mapping['human_symbol_library'] = ~symbol_mapping['library'].str.contains('Mouse', case=False)
symbol_mapping.to_csv(HEALTHY_VS_AKI_OUTPUT_DIR / 'pathway_symbol_mapping_report.csv', index=False)
display(symbol_mapping.round(3))
for _r in symbol_mapping_report:
    if _r['symbol_match_rate'] < 0.5:
        _cause = ('cross-species symbol mismatch and/or PT detection filtering'
                  if 'mouse' not in _r['library'].lower() else
                  'PT detection filtering (this library is already mouse-symbol, so species is '
                  'NOT the cause)')
        print(f"WARNING: only {_r['symbol_match_rate']:.0%} of {_r['library']} symbols matched the "
              f'tested PT gene background -- attributable to {_cause}. Module scores from this '
              'library cover a biased subset of each gene set.')

module_scores = np.empty((adata_pt.n_obs, P), dtype=np.float64)
for j, members in enumerate(retained_pathways['genes_present']):
    idxs = [gene_local[g_] for g_ in members]
    sub = Y_csc[:, idxs].toarray()
    module_scores[:, j] = ((sub - gene_mean[idxs]) / gene_std[idxs]).mean(axis=1)

path_ls = run_level_shape(module_scores, s, c, knots, grid, SECTION4_CONFIG['lambda_grid'])
path_shape_p, path_level_p, _, _ = sample_perm_pvalues(
    module_scores, s, samples, knots, grid, SECTION4_CONFIG['lambda_grid'],
    path_ls['lam_idx'], path_ls['sl2'], path_ls['p_b'], control_samples)
print(f'Pathway modules tested: {P}  '
      f"(shape_rms max={path_ls['shape_rms'].max():.4f})")

# %% [markdown]
# ## 4.6 - Analysis B: pseudospace-structure diagnostics (secondary)

# %%
# ----------------------------------------------------------------------------
# Section 4.6 -- Analysis B: pseudospace-structure diagnostics (secondary)
# Flags whether a shape hit could be driven by density redistribution rather than expression.
# ----------------------------------------------------------------------------
s_h = s[c == 0]
s_a = s[c == 1]
wass = float(wasserstein_distance(s_h, s_a))
ks_stat, ks_p = ks_2samp(s_h, s_a)


# NB: this is a PT-INTERNAL S1->S3 axis built from MARKER_GROUPS, and is NOT the global
# early->late axis of Section 2e (`total_marker_axis`, from TOTAL_POSITION_MARKERS). Different
# anchors, different scope -- hence the distinct name, so the two can never be confused or
# accidentally shadow one another.
def _pt_s1s3_score(group):
    present = resolve_present(MARKER_GROUPS[group], gene_names)
    if not present:
        return np.zeros(adata_pt.n_obs), present
    idxs = [gene_local[g_] for g_ in present]
    sub = Y_csc[:, idxs].toarray()
    return ((sub - gene_mean[idxs]) / gene_std[idxs]).mean(axis=1), present


s1_axis, s1_present = _pt_s1s3_score('S1')
s3_axis, s3_present = _pt_s1s3_score('S3')
pt_s1s3_axis = s3_axis - s1_axis                             # low (S1) -> high (S3)
mono_healthy = safe_spearman(pt_s1s3_axis[c == 0], s[c == 0])
mono_aki = safe_spearman(pt_s1s3_axis[c == 1], s[c == 1])

structure_diag = pd.DataFrame([
    {'metric': 'wasserstein_distance_healthy_aki', 'value': wass},
    {'metric': 'ks_statistic', 'value': float(ks_stat)},
    {'metric': 'ks_pvalue', 'value': float(ks_p)},
    {'metric': 'marker_axis_monotonicity_spearman_healthy', 'value': mono_healthy},
    {'metric': 'marker_axis_monotonicity_spearman_aki', 'value': mono_aki},
    {'metric': 's1_markers_present', 'value': ';'.join(s1_present)},
    {'metric': 's3_markers_present', 'value': ';'.join(s3_present)},
    {'metric': 'n_cells_healthy', 'value': n_healthy},
    {'metric': 'n_cells_aki', 'value': n_aki},
])
structure_diag.to_csv(HEALTHY_VS_AKI_OUTPUT_DIR / 'pseudospace_structure_diagnostics.csv', index=False)

# per-condition density on the common grid (also reused for the optional overlay)
dens_edges = np.linspace(lo, hi, SECTION4_CONFIG['grid_points'] + 1)
dens_centers = 0.5 * (dens_edges[:-1] + dens_edges[1:])
dens_healthy, _ = np.histogram(s_h, bins=dens_edges, density=True)
dens_aki, _ = np.histogram(s_a, bins=dens_edges, density=True)

fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
axes[0].fill_between(dens_centers, dens_healthy, alpha=0.45,
                     color=SECTION4_CONFIG['condition_colors']['healthy'], label='healthy')
axes[0].fill_between(dens_centers, dens_aki, alpha=0.45,
                     color=SECTION4_CONFIG['condition_colors']['aki'], label='aki')
axes[0].set_xlabel('Shared pseudospace')
axes[0].set_ylabel('PT cell density')
axes[0].set_title(f'PT density  (Wasserstein={wass:.3f}, KS p={ks_p:.2g})')
axes[0].legend()
rng = np.random.default_rng(SECTION4_CONFIG['random_state'])
for lab, mask, col in [('healthy', c == 0, SECTION4_CONFIG['condition_colors']['healthy']),
                       ('aki', c == 1, SECTION4_CONFIG['condition_colors']['aki'])]:
    idx = np.flatnonzero(mask)
    sel = rng.choice(idx, size=min(3000, idx.size), replace=False)
    axes[1].scatter(s[sel], pt_s1s3_axis[sel], s=4, alpha=0.25, color=col, label=lab)
axes[1].set_xlabel('Shared pseudospace')
axes[1].set_ylabel('S1->S3 marker axis  (z(S3) - z(S1))')
axes[1].set_title(f'Marker-axis monotonicity  (healthy rho={mono_healthy:.2f}, aki rho={mono_aki:.2f})')
axes[1].legend()
fig.tight_layout()
fig.savefig(HEALTHY_VS_AKI_OUTPUT_DIR / 'pseudospace_structure_diagnostics.png', bbox_inches='tight', dpi=130)
plt.close(fig)
display(Image(filename=str(HEALTHY_VS_AKI_OUTPUT_DIR / 'pseudospace_structure_diagnostics.png')))
print('Analysis B diagnostics written.')
print(structure_diag.to_string(index=False))

# %% [markdown]
# ## 4.6b - Analysis C: is the pseudospace axis itself condition-dependent?
#
# Analysis B asks whether AKI redistributes tubules *along* the axis. Analysis C asks the harder
# question: does AKI move tubules along the axis for reasons that have nothing to do with
# position? Two ways that happens, both documented in mouse IRI:
#
# 1. **De-differentiation.** Injured PT loses the segment-identity markers the axis is built
#    from, so the axis's own basis shifts under condition (Kirita 2020 PNAS 117:15874).
# 2. **Failed-repair PT is a branch, not a point.** FR-PTC leaves the repair trajectory
#    (Kirita 2021 PNAS 118:e2026684118). DPT ran with n_branchings=0, so those tubules are
#    projected onto the linear axis at an arbitrary position and then averaged into the AKI
#    condition curve.
#
# A shape hit whose signal sits where injury-high tubules pile up is not separable from injury.

# %%
# ----------------------------------------------------------------------------
# Section 4.6b -- Analysis C: injury confounding of the coordinate
# ----------------------------------------------------------------------------
analysis_c_rows = []

# (1) Does the axis's own basis shift with condition? The early-marker score is half of what
# defines the coordinate; if AKI depresses it globally, the coordinate is not condition-invariant.
for _basis_col in [f'{CONTINUUM_PREFIX}_early_trajectory_score',
                   f'{CONTINUUM_PREFIX}_late_trajectory_score']:
    if _basis_col not in adata_pt.obs.columns:
        print(f'NOTE: {_basis_col} absent from the DPT file; regenerate Section 2g to enable the '
              'axis-basis check.')
        continue
    _v = adata_pt.obs[_basis_col].to_numpy(dtype=float)
    _mh, _ma = float(np.nanmean(_v[c == 0])), float(np.nanmean(_v[c == 1]))
    _pooled_sd = float(np.nanstd(_v)) or 1.0
    analysis_c_rows.append({
        'metric': f'{_basis_col}__healthy_mean', 'value': _mh})
    analysis_c_rows.append({
        'metric': f'{_basis_col}__aki_mean', 'value': _ma})
    analysis_c_rows.append({
        'metric': f'{_basis_col}__aki_minus_healthy_cohens_d', 'value': (_ma - _mh) / _pooled_sd})

# (2) Where do injury-high tubules sit on the axis?
for _panel, _score in injury_scores.items():
    analysis_c_rows.append({'metric': f'{_panel}__spearman_with_pseudospace',
                            'value': safe_spearman(_score, s)})
    analysis_c_rows.append({'metric': f'{_panel}__aki_minus_healthy_cohens_d',
                            'value': (float(np.nanmean(_score[c == 1])) - float(np.nanmean(_score[c == 0])))
                                     / (float(np.nanstd(_score)) or 1.0)})
    _hi = _score >= np.nanquantile(_score, 0.90)          # top-decile injury tubules
    analysis_c_rows.append({'metric': f'{_panel}__top_decile_median_pseudospace',
                            'value': float(np.nanmedian(s[_hi]))})
    analysis_c_rows.append({'metric': f'{_panel}__all_median_pseudospace',
                            'value': float(np.nanmedian(s))})

analysis_c = pd.DataFrame(analysis_c_rows)
analysis_c.to_csv(HEALTHY_VS_AKI_OUTPUT_DIR / 'injury_confounding_diagnostics.csv', index=False)
display(analysis_c.round(4))

if injury_scores:
    _n_panels = len(injury_scores)
    fig, axes = plt.subplots(2, _n_panels, figsize=(6 * _n_panels, 8), squeeze=False)
    for _j, (_panel, _score) in enumerate(injury_scores.items()):
        # Row 0: injury score along pseudospace, per condition (binned mean +/- IQR band).
        ax = axes[0, _j]
        for _lab, _m, _col in [('healthy', c == 0, COL_H), ('aki', c == 1, COL_A)]:
            _prof, _ctr = binned_matrix_profile(_score[_m][:, None], s[_m], lo, hi,
                                                SECTION4_CONFIG['heatmap_n_bins'],
                                                SECTION4_CONFIG['heatmap_smooth_sigma'])
            ax.plot(_ctr, _prof[0], color=_col, lw=2, label=_lab)
        ax.set_title(f'{_panel} score along pseudospace')
        ax.set_xlabel('Shared pseudospace')
        ax.set_ylabel(f'{_panel} score')
        ax.legend(fontsize=8)

        # Row 1: where the injury-high decile sits, against the overall distribution.
        ax = axes[1, _j]
        _hi = _score >= np.nanquantile(_score, 0.90)
        ax.hist(s, bins=60, range=(lo, hi), color='#CCCCCC', label='all PT')
        ax.hist(s[_hi], bins=60, range=(lo, hi), color='firebrick', alpha=0.75,
                label=f'{_panel} top decile')
        ax.set_xlabel('Shared pseudospace')
        ax.set_ylabel('PT tubules')
        ax.legend(fontsize=8)
    fig.suptitle('Analysis C: injury confounding of the pseudospace coordinate', y=1.01, fontsize=12)
    _save(fig, 'injury_confounding_diagnostics.png')

    print('\nHOW TO READ THIS: if a panel\'s top decile piles up in one stretch of pseudospace,\n'
          'every Section 4 shape hit whose signal sits in that stretch is confounded with injury\n'
          'state and cannot be attributed to a positional shift. Cross-check the top shape genes\n'
          'in Section 4.7 against those pseudospace ranges before interpreting them.')

# %% [markdown]
# ## 4.6c - Pseudospace-matched injury test (de-confounds position and sex)
#
# A raw gene<->injury correlation is confounded: injured tubules pile up at high pseudospace, so a
# late-segment gene looks "injury-associated" for free. This test asks the honest question -- at
# MATCHED position, do injured tubules differ from uninjured ones? -- via
#   (a) partial Spearman of gene~injury controlling for pseudospace (and optionally sex), and
#   (b) within-pseudospace-bin association computed INSIDE AKI, where injury score actually varies.
#
# If the matched effect collapses toward 0 relative to the naive effect, the gene was position-
# confounded. If it collapses only after adding the sex covariate, it was sex-confounded. A gene
# survives only if the matched, sex-adjusted, within-AKI association stays strong AND consistent in
# sign across bins.
#
# CAVEAT: with 2 healthy + 2 AKI specimens the significance ceiling is set by specimens, not
# tubules. Read the effect sizes and sign-consistency; treat any tubule-level p as descriptive.

# %%
# ----------------------------------------------------------------------------
# Section 4.6c -- position/sex-matched injury association
# ----------------------------------------------------------------------------
from scipy.stats import rankdata

import scipy.sparse as sp

def _dense_col(mat, j):
    col = mat[:, j]
    if sp.issparse(col):
        col = col.toarray()
    return np.asarray(col, dtype=float).ravel()

# Known sex-biased mouse PT genes (from the sex-dimorphism literature). Used only to FLAG hits so
# you don't over-interpret them; the sex covariate below is the real correction.
SEX_BIASED_PT_GENES = {
    # male-biased
    'Cyp4b1', 'Cyp2j13', 'Cyp7b1', 'Cyp2e1', 'Slc22a22', 'Slc22a28', 'Slc22a30',
    'Acsm3', 'Slc7a13', 'Atp11a', 'Akr1c21', 'Slco1a1',
    # female-biased
    'Slc7a12', 'Kynu', 'Prlr', 'Cyp4a14',
}

MATCH_CONFIG = {
    'n_bins': SECTION4_CONFIG.get('heatmap_n_bins', 20),
    'min_per_bin': 20,          # skip bins too sparse to estimate a within-bin rho
    'collapse_ratio': 0.5,      # |partial| < collapse_ratio * |naive|  => "confounded"
    'strong_rho': 0.10,         # |matched rho| above this = worth keeping
    'sign_consistency': 0.70,   # fraction of used bins agreeing in sign
}

# Optional per-specimen sex, e.g. {'H1': 'M', 'H2': 'F', 'A1': 'M', 'A2': 'F'}.
# If left as None the sex covariate is skipped (but genes are still name-flagged).
SPECIMEN_SEX = globals().get('SPECIMEN_SEX', None)
if SPECIMEN_SEX is not None:
    sex_code = np.array([1.0 if SPECIMEN_SEX.get(sp_, 'U') == 'M' else 0.0 for sp_ in samples])
else:
    sex_code = None


def _rank_residual(x, covars):
    """Residual of rank(x) after least-squares regression on ranked covariates (+ intercept)."""
    xr = rankdata(x).astype(float)
    cols = [rankdata(cv).astype(float) for cv in covars]
    C = np.column_stack(cols + [np.ones_like(xr)])
    beta, *_ = np.linalg.lstsq(C, xr, rcond=None)
    return xr - C @ beta


def partial_spearman(x, y, covars, mask):
    """Spearman(x, y) controlling for covars, computed on `mask` rows only."""
    m = mask & np.isfinite(x) & np.isfinite(y)
    for cv in covars:
        m = m & np.isfinite(cv)
    if m.sum() < 10:
        return np.nan
    rx = _rank_residual(x[m], [cv[m] for cv in covars])
    ry = _rank_residual(y[m], [cv[m] for cv in covars])
    if np.std(rx) == 0 or np.std(ry) == 0:
        return np.nan
    return float(np.corrcoef(rx, ry)[0, 1])


def within_bin_association(gene, injury, s_all, mask, n_bins, min_per_bin):
    """Per-bin Spearman(gene, injury) inside `mask`; returns (n-weighted mean rho, sign-consistency, n_bins)."""
    edges = np.linspace(lo, hi, n_bins + 1)
    rhos, ns = [], []
    for b in range(n_bins):
        mb = mask & (s_all >= edges[b]) & (s_all < edges[b + 1])
        if mb.sum() < min_per_bin:
            continue
        r = safe_spearman(gene[mb], injury[mb])
        if np.isfinite(r):
            rhos.append(r)
            ns.append(mb.sum())
    if not rhos:
        return np.nan, np.nan, 0
    rhos, ns = np.array(rhos), np.array(ns, dtype=float)
    mean_rho = float(np.sum(rhos * ns) / ns.sum())
    dominant = np.sign(mean_rho) if mean_rho != 0 else 1.0
    consistency = float(np.mean(np.sign(rhos) == dominant))
    return mean_rho, consistency, len(rhos)


# genes: top-N by shape_rms
gene_top = np.argsort(gene_ls['shape_rms'])[::-1][:N_TOP]
gt_names = gene_names[gene_top]
gt_h_raw = gene_ls['curve_healthy'][gene_top]
gt_a_raw = gene_ls['curve_aki'][gene_top]
gt_mean = gene_mean[gene_top][:, None]
gt_std = gene_std[gene_top][:, None]
gt_h_z = (gt_h_raw - gt_mean) / gt_std
gt_a_z = (gt_a_raw - gt_mean) / gt_std

# Run on the top shape-ranked genes. gene_top / gt_names are DEFINED just above, in this cell --
# Section 4.8 reuses them, not the other way round. Nothing here depends on 4.8 having run.
mask_aki = c == 1
matched_rows = []
for gi, gname in zip(gene_top, gt_names):
    g = _dense_col(Y_genes, gi)
    for _panel, _score in injury_scores.items():
        naive = safe_spearman(g, _score)                                   # confounded baseline
        part_pos = partial_spearman(g, _score, [s], np.ones_like(g, bool)) # control position
        covars_ps = [s] if sex_code is None else [s, sex_code]
        part_pos_sex = partial_spearman(g, _score, covars_ps, np.ones_like(g, bool))
        aki_part = partial_spearman(g, _score, [s], mask_aki)              # within AKI, matched position
        wb_rho, wb_sign, wb_nbins = within_bin_association(
            g, _score, s, mask_aki, MATCH_CONFIG['n_bins'], MATCH_CONFIG['min_per_bin'])

        # verdicts
        pos_conf = (abs(part_pos) < MATCH_CONFIG['collapse_ratio'] * abs(naive)) if np.isfinite(naive) and naive != 0 else np.nan
        sex_flag = gname in SEX_BIASED_PT_GENES
        sex_conf = (sex_code is not None and np.isfinite(part_pos) and np.isfinite(part_pos_sex)
                    and abs(part_pos_sex) < MATCH_CONFIG['collapse_ratio'] * abs(part_pos))
        survives = (np.isfinite(aki_part) and abs(aki_part) >= MATCH_CONFIG['strong_rho']
                    and np.isfinite(wb_sign) and wb_sign >= MATCH_CONFIG['sign_consistency']
                    and not sex_conf)

        if survives:
            verdict = 'injury-associated (matched)'
        elif sex_conf:
            verdict = 'sex-confounded'
        elif pos_conf is True:
            verdict = 'position-confounded'
        else:
            verdict = 'ambiguous / weak'

        matched_rows.append({
            'gene': gname, 'panel': _panel,
            'naive_rho_all': naive,
            'partial_rho_pos': part_pos,
            'partial_rho_pos_sex': part_pos_sex,
            'aki_partial_rho_pos': aki_part,
            'within_bin_mean_rho_aki': wb_rho,
            'within_bin_sign_consistency': wb_sign,
            'within_bin_n_bins': wb_nbins,
            'sex_biased_gene': sex_flag,
            'verdict': verdict,
        })

matched = pd.DataFrame(matched_rows)
matched.to_csv(HEALTHY_VS_AKI_OUTPUT_DIR / 'injury_matched_diagnostics.csv', index=False)
display(matched.round(4))

# ---- Figure: naive vs matched effect, and within-bin profiles for a few example genes ----------
if len(matched):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # (left) naive vs position-matched partial rho -- points falling toward the x=0 line collapsed
    ax = axes[0]
    for _panel in injury_scores:
        sub = matched[matched['panel'] == _panel]
        ax.scatter(sub['naive_rho_all'], sub['partial_rho_pos'], s=28, alpha=0.8, label=_panel)
    _lim = np.nanmax(np.abs([matched['naive_rho_all'], matched['partial_rho_pos']])) * 1.1
    ax.plot([-_lim, _lim], [-_lim, _lim], color='grey', ls='--', lw=1, label='no confound (y=x)')
    ax.axhline(0, color='k', lw=0.8)
    ax.set_xlabel('naive Spearman (gene vs injury, ALL tubules)')
    ax.set_ylabel('partial Spearman (control pseudospace)')
    ax.set_title('Collapse toward y=0 => position-confounded')
    ax.legend(fontsize=8)

    # (right) within-AKI, within-bin rho vs pseudospace for the strongest surviving/failing genes
    ax = axes[1]
    _first_panel = next(iter(injury_scores))
    _score0 = injury_scores[_first_panel]
    edges = np.linspace(lo, hi, MATCH_CONFIG['n_bins'] + 1)
    ctrs = 0.5 * (edges[:-1] + edges[1:])
    _examples = matched[matched['panel'] == _first_panel].reindex(
        matched[matched['panel'] == _first_panel]['aki_partial_rho_pos'].abs().sort_values(ascending=False).index
    ).head(4)
    for _, row in _examples.iterrows():
        g = _dense_col(Y_genes, gene_top[list(gt_names).index(row['gene'])])
        prof = []
        for b in range(MATCH_CONFIG['n_bins']):
            mb = mask_aki & (s >= edges[b]) & (s < edges[b + 1])
            prof.append(safe_spearman(g[mb], _score0[mb]) if mb.sum() >= MATCH_CONFIG['min_per_bin'] else np.nan)
        ax.plot(ctrs, prof, marker='o', ms=3, lw=1.2,
                label=f"{row['gene']} [{row['verdict']}]")
    ax.axhline(0, color='k', lw=0.8)
    ax.set_xlabel('Shared pseudospace (AKI only)')
    ax.set_ylabel(f'within-bin Spearman (gene vs {_first_panel})')
    ax.set_title('Consistent nonzero across bins = real matched association')
    ax.legend(fontsize=7)

    fig.suptitle('Analysis C+: pseudospace/sex-matched injury association', y=1.02, fontsize=12)
    _save(fig, 'injury_matched_diagnostics.png')

    print('\nHOW TO READ THIS:')
    print('  * naive_rho_all      = confounded baseline (gene vs injury, all tubules).')
    print('  * partial_rho_pos    = same, controlling pseudospace. If it collapses toward 0 the')
    print('                         gene was a POSITION effect masquerading as injury.')
    print('  * partial_rho_pos_sex= adds a sex covariate. A drop here means SEX confounding')
    print('                         (expected for Cyp4b1, Cyp7b1, Cyp2e1, Slc22a28/30, ...).')
    print('  * aki_partial_rho_pos + within_bin_* = the real test, computed INSIDE AKI at matched')
    print('                         position. Only genes strong AND sign-consistent here survive.')
    print('  * Set SPECIMEN_SEX = {"H1":"M",...} to activate the sex correction (currently '
          f'{"ON" if sex_code is not None else "OFF"}).')

# %% [markdown]
# ## 4.6d - Analysis C+: coordinate-defining vs coordinate-independent injury signal
#
# 4.6c asks whether a gene<->injury association survives position matching. This asks the follow-up
# question: *even if it survives, is the gene allowed to be evidence?* Three ways it is not:
#
#   1. the gene is a member of the injury panel it is being tested against (circular);
#   2. the gene is part of the marker basis that defines/orients the pseudospace coordinate, so
#      de-differentiation under injury bends the ruler rather than moving the tubule;
#   3. the association collapses once position is controlled (position-confounded).
#
# What is left -- strong inside AKI at matched position, sign-consistent across bins, not panel-
# derived, not axis-basis -- is the shortlist worth chasing.
#
# Sex is deliberately NOT modeled here: this cohort is all-male, so a sex covariate has zero
# variance and the design matrix is rank-deficient. See 4.6c for the sex-aware variant.

# %%
# ----------------------------------------------------------------------------
# Section 4.6d -- Analysis C+: adjudicate the top shape hits
# Self-contained: re-derives its own gene ranking from gene_ls (Section 4.3) so it does not depend
# on Section 4.8 having been run first, and prefixes every helper with `_adj_` so it cannot shadow
# the 4.6c definitions of partial_spearman / within_bin_association / MATCH_CONFIG.
# ----------------------------------------------------------------------------
ADJUDICATE_CONFIG = {
    'strong_rho':       0.15,   # |partial rho| must clear this to count as "associated"
    'sign_consistency': 0.80,   # fraction of within-bin rho sharing the mean sign
    'collapse_ratio':   0.50,   # |partial| < ratio * |naive| => position-confounded
    'min_per_bin':      15,     # min AKI tubules per pseudospace bin
    'n_bins':           60,     # within-bin sweep resolution over the common support [lo, hi]
}


def _adj_norm(g):
    """Case-insensitive gene key."""
    return str(g).strip().lower()


# --- gene sets, derived from the globals that actually define the analysis -------------------
# Anything that defines or orients the coordinate: the S1/S2/S3 marker panel (Section 4.0),
# the early/late anchors used to orient DPT (Section 0.1), plus the S1/S2-shared pan-PT anchor
# Slc34a1, which was dropped from MARKER_GROUPS but still loads on the same axis.
ADJ_AXIS_BASIS_EXTRA = ['Slc34a1']
ADJ_AXIS_BASIS_GENES = (
    {_adj_norm(g) for grp in MARKER_GROUPS.values() for g in grp}
    | {_adj_norm(g) for g in TOTAL_POSITION_MARKERS['early']}
    | {_adj_norm(g) for g in TOTAL_POSITION_MARKERS['late']}
    | {_adj_norm(g) for g in ADJ_AXIS_BASIS_EXTRA}
)

# Circularity guard: read the panel membership off INJURY_PANELS (Section 4.1b), never a
# hand-copied list -- a stale copy silently disables the guard for whatever it omits.
_adj_unknown_panels = sorted(set(injury_scores) - set(INJURY_PANELS))
if _adj_unknown_panels:
    raise KeyError(
        f'injury_scores contains panels {_adj_unknown_panels} that are absent from INJURY_PANELS, '
        'so the circularity guard cannot know which genes define them. Add them to INJURY_PANELS '
        '(Section 4.1b) before running this cell.'
    )
ADJ_PANEL_GENES = {
    panel: {_adj_norm(g) for g in INJURY_PANELS[panel]} for panel in injury_scores
}

if SPECIMEN_SEX is not None and len(set(SPECIMEN_SEX.values())) > 1:
    print('WARNING: SPECIMEN_SEX reports more than one sex, but this section does not model sex. '
          'Use the 4.6c table (partial_rho_pos_sex) for the sex-adjusted verdicts.')


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
def _adj_rank_residual(v, C):
    """Residual of rank(v) after OLS on the ranked columns of C (+ intercept). C may be None."""
    rv = rankdata(v).astype(float)
    if C is None or C.size == 0:
        return rv - rv.mean()
    A = np.column_stack(
        [np.ones(rv.size)] + [rankdata(C[:, k]).astype(float) for k in range(C.shape[1])])
    beta, *_ = np.linalg.lstsq(A, rv, rcond=None)
    return rv - A @ beta


def _adj_partial_spearman(y, x, covars=None):
    """
    Partial Spearman(y, x) controlling for `covars` (1D array, (n, k) array, or None for the
    plain Spearman). Rows non-finite in any column are dropped first. Returns nan if degenerate.
    """
    y = np.asarray(y, float).ravel()
    x = np.asarray(x, float).ravel()
    stack = [y, x]
    if covars is not None:
        C = np.asarray(covars, float)
        if C.ndim == 1:
            C = C[:, None]
        elif C.shape[0] != y.size and C.shape[1] == y.size:
            C = C.T
        if C.shape[0] != y.size:
            raise ValueError(f'covars has {C.shape[0]} rows but y has {y.size}')
        stack += [C[:, k] for k in range(C.shape[1])]
    M = np.column_stack(stack)
    ok = np.all(np.isfinite(M), axis=1)
    if ok.sum() < 10:
        return np.nan
    M = M[ok]
    Cv = M[:, 2:] if M.shape[1] > 2 else None
    ey = _adj_rank_residual(M[:, 0], Cv)
    ex = _adj_rank_residual(M[:, 1], Cv)
    if ey.std() < 1e-12 or ex.std() < 1e-12:
        return np.nan
    return float(np.corrcoef(ey, ex)[0, 1])


def _adj_within_bin(y, x, pos, n_bins, min_per_bin, lo_=None, hi_=None):
    """
    Spearman(y, x) inside equal-width pseudospace bins spanning [lo_, hi_] (the common support by
    default, matching 4.6c). Inputs are already AKI-only. Cells outside the span are dropped, not
    clipped into the edge bins. Returns (n-weighted mean rho, sign consistency, n bins used).
    """
    y = np.asarray(y, float); x = np.asarray(x, float); pos = np.asarray(pos, float)
    lo_ = lo if lo_ is None else lo_
    hi_ = hi if hi_ is None else hi_
    edges = np.linspace(lo_, hi_, n_bins + 1)
    finite = np.isfinite(y) & np.isfinite(x) & np.isfinite(pos)
    rhos, ns = [], []
    for b in range(n_bins):
        upper = (pos <= edges[b + 1]) if b == n_bins - 1 else (pos < edges[b + 1])
        m = finite & (pos >= edges[b]) & upper
        if m.sum() < min_per_bin:
            continue
        r = safe_spearman(y[m], x[m])
        if np.isfinite(r):
            rhos.append(r)
            ns.append(int(m.sum()))
    if not rhos:
        return np.nan, np.nan, 0
    rhos = np.asarray(rhos, float)
    ns = np.asarray(ns, float)
    mean_rho = float(np.sum(rhos * ns) / ns.sum())
    dom_sign = np.sign(mean_rho) if mean_rho != 0 else 1.0
    sign_consistency = float(np.mean(np.sign(rhos) == dom_sign))
    return mean_rho, sign_consistency, len(rhos)


# -----------------------------------------------------------------------------
# Adjudication
# -----------------------------------------------------------------------------
def adjudicate_shape_hits(gene_idx, gene_labels, cfg=ADJUDICATE_CONFIG):
    """
    Classify each (gene, injury panel) pair using the module-level PT cohort
    (Y_genes / s / mask_aki / injury_scores).
    """
    pos_aki = s[mask_aki]
    rows = []

    for gi, gname in zip(gene_idx, gene_labels):
        g_all = _dense_col(Y_genes, gi)
        g_aki = g_all[mask_aki]
        gkey = _adj_norm(gname)
        is_axis_basis = gkey in ADJ_AXIS_BASIS_GENES

        for panel, score in injury_scores.items():
            score = np.asarray(score, float)
            s_aki_score = score[mask_aki]

            # --- circularity guard: the gene helps define this injury score ---------------
            if gkey in ADJ_PANEL_GENES[panel]:
                rows.append(dict(
                    gene=gname, panel=panel,
                    verdict='excluded (panel member / circular)',
                    naive_rho_all=np.nan, partial_rho_pos=np.nan,
                    aki_partial_rho_pos=np.nan, within_bin_mean_rho=np.nan,
                    within_bin_sign_consistency=np.nan, within_bin_n_bins=0,
                    is_axis_basis=is_axis_basis))
                continue

            # --- statistics ---------------------------------------------------------------
            naive    = safe_spearman(g_all, score)                             # confounded baseline
            part_pos = _adj_partial_spearman(g_all, score, covars=s)           # control position
            aki_part = _adj_partial_spearman(g_aki, s_aki_score, covars=pos_aki)  # AKI-only, matched
            wb_mean, wb_sign, wb_n = _adj_within_bin(
                g_aki, s_aki_score, pos_aki, cfg['n_bins'], cfg['min_per_bin'])

            # --- survival test ------------------------------------------------------------
            survives = bool(
                np.isfinite(aki_part) and abs(aki_part) >= cfg['strong_rho']
                and np.isfinite(wb_sign) and wb_sign >= cfg['sign_consistency']
                and np.isfinite(wb_mean) and abs(wb_mean) >= cfg['strong_rho']
            )
            collapsed = bool(
                np.isfinite(naive) and naive != 0
                and np.isfinite(part_pos)
                and abs(part_pos) < cfg['collapse_ratio'] * abs(naive)
            )

            # --- verdict (the coordinate split) -------------------------------------------
            if survives and is_axis_basis:
                verdict = 'injury-associated BUT axis-basis (de-diff confounds coordinate)'
            elif survives:
                verdict = 'coordinate-independent injury signal (SHORTLIST)'
            elif collapsed:
                verdict = 'position-confounded'
            else:
                verdict = 'weak / ambiguous'

            rows.append(dict(
                gene=gname, panel=panel, verdict=verdict,
                naive_rho_all=naive, partial_rho_pos=part_pos,
                aki_partial_rho_pos=aki_part,
                within_bin_mean_rho=wb_mean,
                within_bin_sign_consistency=wb_sign,
                within_bin_n_bins=wb_n,
                is_axis_basis=is_axis_basis))

    # Explicit columns so an empty `rows` (no injury panel resolved) still yields a well-formed
    # frame rather than a column-less one that blows up on sort_values/`['gene']` downstream.
    df = pd.DataFrame(rows, columns=[
        'gene', 'panel', 'verdict', 'naive_rho_all', 'partial_rho_pos', 'aki_partial_rho_pos',
        'within_bin_mean_rho', 'within_bin_sign_consistency', 'within_bin_n_bins',
        'is_axis_basis'])
    order = ['coordinate-independent injury signal (SHORTLIST)',
             'injury-associated BUT axis-basis (de-diff confounds coordinate)',
             'position-confounded', 'weak / ambiguous',
             'excluded (panel member / circular)']
    df['verdict'] = pd.Categorical(df['verdict'], categories=order, ordered=True)
    return df.sort_values(['verdict', 'panel', 'gene']).reset_index(drop=True)


# -----------------------------------------------------------------------------
# Run + report
# -----------------------------------------------------------------------------
# Same ranking as Section 4.8 (top-N by shape_rms), recomputed locally so cell order cannot matter.
adj_gene_top = np.argsort(gene_ls['shape_rms'])[::-1][:N_TOP]
adj_gt_names = gene_names[adj_gene_top]

verdict_df = adjudicate_shape_hits(adj_gene_top, adj_gt_names)
verdict_df.to_csv(HEALTHY_VS_AKI_OUTPUT_DIR / 'injury_shape_adjudication.csv', index=False)

print('\n=== Analysis C+ verdicts (all-male cohort; sex not modeled) ===')
print(verdict_df.to_string(index=False,
                           float_format=lambda x: f'{x: .3f}' if np.isfinite(x) else '   nan'))

print('\n--- counts by verdict ---')
print(verdict_df['verdict'].value_counts().to_string())

shortlist = verdict_df[
    verdict_df['verdict'] == 'coordinate-independent injury signal (SHORTLIST)']
print(f"\nCoordinate-independent shortlist ({shortlist['gene'].nunique()} genes): "
      f"{sorted(map(str, shortlist['gene'].unique()))}")
print('Axis-basis genes excluded from that shortlist by construction: '
      f'{sorted(ADJ_AXIS_BASIS_GENES)}')

# %%
# ----------------------------------------------------------------------------
# Section 4.7 -- result CSVs (ranked by shape_rms)
# ----------------------------------------------------------------------------
gene_results = pd.DataFrame({
    'gene': gene_names,
    # NB: these two F statistics are computed with n = number of TUBULES, but the effective n
    # for a condition contrast is 4 specimens. They are retained for diagnostics and are NOT
    # calibrated tests -- hence the explicit column names. Rank by shape_rms.
    'level_F_cellwise_uncalibrated': gene_ls['level_F'],
    'level_effect': gene_ls['level_effect'],
    'level_abs': np.abs(gene_ls['level_effect']),
    'm2_condition_coef': gene_ls['m2_condition_coef'],
    'shape_F_cellwise_uncalibrated': gene_ls['shape_F'],
    'shape_rms': gene_ls['shape_rms'],
    'condition_effect_rms': gene_ls['condition_effect_rms'],
    'condition_effect_max_abs': gene_ls['condition_effect_max_abs'],
    'curve_spearman': gene_ls['curve_spearman'],
    'amplitude_healthy': gene_ls['amplitude_healthy'],
    'amplitude_aki': gene_ls['amplitude_aki'],
    'sample_perm_p': gene_shape_p,  # FDR omitted: p floored at 2/6 (n=2v2) -> BH non-functional
    'level_perm_p': gene_level_p,
    'n_cells_healthy': n_healthy,
    'n_cells_aki': n_aki,
})
# `shape_rms` alone counts a weaker gradient as a shape change, and `condition_effect_rms` mixes a
# vertical offset with a redistribution. The split adds level, amplitude and pattern explicitly.
gene_curve_effects = summarize_curve_effects(
    gene_ls['curve_healthy'], gene_ls['curve_aki'], feature_names=gene_names
)
gene_results = gene_results.merge(
    gene_curve_effects.rename(columns={'feature': 'gene'}), on='gene', how='left'
).sort_values('shape_rms', ascending=False).reset_index(drop=True)
top_shape_genes = gene_results.head(20)
print('Top-20 shape-ranked genes: median level fraction '
      f"{top_shape_genes['level_fraction'].median():.3f}; "
      f"types {top_shape_genes['difference_type'].value_counts().to_dict()}")
gene_results.to_csv(HEALTHY_VS_AKI_OUTPUT_DIR / 'gene_level_shape_results.csv', index=False)

pathway_results = retained_pathways.drop(columns='genes_present').assign(
    level_F_cellwise_uncalibrated=path_ls['level_F'],
    level_effect=path_ls['level_effect'],
    level_abs=np.abs(path_ls['level_effect']),
    m2_condition_coef=path_ls['m2_condition_coef'],
    shape_F_cellwise_uncalibrated=path_ls['shape_F'],
    shape_rms=path_ls['shape_rms'],
    condition_effect_rms=path_ls['condition_effect_rms'],
    condition_effect_max_abs=path_ls['condition_effect_max_abs'],
    curve_spearman=path_ls['curve_spearman'],
    sample_perm_p=path_shape_p,
    level_perm_p=path_level_p,
    n_cells_healthy=n_healthy,
    n_cells_aki=n_aki,
).sort_values('shape_rms', ascending=False).reset_index(drop=True)
pathway_results.to_csv(HEALTHY_VS_AKI_OUTPUT_DIR / 'pathway_level_shape_results.csv', index=False)

display(gene_results.head(10))
display(pathway_results[['library', 'pathway', 'n_genes_present', 'level_effect',
                         'shape_rms', 'sample_perm_p']].head(10))

# %%
# Leave-one-specimen-out robustness of the top shape hits (honest ceiling at n=2 vs 2).
from pseudospace.levelshape import loso_shape_stability
top_idx = np.argsort(gene_ls['shape_rms'])[::-1][:SECTION4_CONFIG['N_TOP_PLOT']]
loso = loso_shape_stability(Y_genes, s, c, samples, knots, grid, SECTION4_CONFIG['lambda_grid'],
                            top_idx, feature_names=gene_names)
loso.to_csv(HEALTHY_VS_AKI_OUTPUT_DIR / 'loso_shape_stability_genes.csv', index=False)
feats = list(gene_names[top_idx])
fig, ax = plt.subplots(figsize=(8, 0.4 * len(feats) + 1.5))
for j, feat in enumerate(feats):
    sub = loso[loso['feature'] == feat]
    ax.plot(sub['loso_shape_rms'], [j] * len(sub), 'o', color='steelblue', alpha=0.7)
    ax.plot([sub['baseline_shape_rms'].iloc[0]], [j], 'D', color='firebrick')
ax.set_yticks(range(len(feats))); ax.set_yticklabels(feats, fontsize=8)
ax.set_xlabel('shape_rms  (diamond = all specimens; circles = leave-one-out)')
ax.set_title('Leave-one-specimen-out robustness of top shape genes')
ax.invert_yaxis()
fig.savefig(HEALTHY_VS_AKI_OUTPUT_DIR / 'loso_shape_stability_genes.png', dpi=130, bbox_inches='tight')
plt.show()
display(loso.head(12))


# %% [markdown]
# ## 4.8 - visualization layer (reuses fitted curves + result CSVs; no stat recompute)

# %%
# ----------------------------------------------------------------------------
# Section 4.8 -- VISUALIZATION LAYER (reuses fitted curves + result CSVs; no stat recompute)
# Group 1 -- familiar per-condition trajectory views
# ----------------------------------------------------------------------------
def plot_trajectory_grid(labels, curve_h, curve_a, ylabel, suptitle, name,
                         densities=None, specimen_curves=None):
    """Pooled condition curves, optionally overlaid with per-specimen curves.

    ``specimen_curves`` is ``{specimen: (curves, condition_code)}`` where ``curves`` is
    (n_labels, len(grid)) with NaN outside that specimen's own support. Drawn as thin lines under
    the thick pooled curves: with 2 healthy and 2 AKI specimens the pooled pair alone implies far
    more confidence than the design supports, and the condition gap is only credible where one
    condition's thin lines stay clear of the other's.
    """
    fig, axes = _grid_axes(len(labels))
    for ax, i in zip(axes, range(len(labels))):
        if specimen_curves is not None:
            for _spec, (_curves, _code) in specimen_curves.items():
                ax.plot(grid, _curves[i], color=COL_H if _code == 0 else COL_A,
                        lw=0.9, alpha=0.65, ls='--', zorder=1)
        ax.plot(grid, curve_h[i], color=COL_H, lw=2.4, label='healthy (pooled)', zorder=3)
        ax.plot(grid, curve_a[i], color=COL_A, lw=2.4, label='aki (pooled)', zorder=3)
        if densities is not None:
            twin = ax.twinx()
            twin.fill_between(dens_centers, densities[0], color=COL_H, alpha=0.08)
            twin.fill_between(dens_centers, densities[1], color=COL_A, alpha=0.08)
            twin.set_yticks([])
        ax.set_title(str(labels[i])[:48], fontsize=9)
        ax.set_xlabel('Shared pseudospace')
        ax.set_ylabel(ylabel)
    handles, lbls = axes[0].get_legend_handles_labels()
    if specimen_curves is not None:
        handles.append(plt.Line2D([], [], color='grey', lw=0.9, ls='--'))
        lbls.append('individual specimens')
    axes[0].legend(handles, lbls, fontsize=7)
    for ax in axes[len(labels):]:
        ax.axis('off')
    fig.suptitle(suptitle, y=1.005, fontsize=12)
    _save(fig, name)



# ---- Per-specimen curves for the top-N genes -------------------------------------------------
# The pooled 2-curve view is the single most over-confident figure this pipeline produces: it
# renders a 2-vs-2 specimen comparison as if it were two well-sampled populations. Refit the
# shared-smooth model within each specimen (reusing that gene's pooled GCV lambda, so spread is
# between-specimen variation and not a smoothing difference) and draw all four underneath.
from pseudospace.levelshape import fit_single_condition_curves

specimen_codes = {sp_: int(c[samples == sp_][0]) for sp_ in sorted(set(samples))}
gene_specimen_curves = {}
for _spec, _code in specimen_codes.items():
    _m = samples == _spec
    _curves, _ = fit_single_condition_curves(
        Y_genes[_m][:, gene_top], s[_m], knots, grid, SECTION4_CONFIG['lambda_grid'],
        gene_ls['lam_idx'][gene_top], support_pct=SECTION4_CONFIG['common_support_pct'])
    gene_specimen_curves[_spec] = (_curves, _code)
    print(f'  per-specimen fit: {_spec} ({"aki" if _code else "healthy"}), n={int(_m.sum()):,}')

_gs_z = {k: ((v - gt_mean) / gt_std, cd) for k, (v, cd) in gene_specimen_curves.items()}

plot_trajectory_grid(gt_names, gt_h_z, gt_a_z, 'Mean gene z-score',
                     'Top shape-ranked gene trajectories (z-score, per-specimen overlay)',
                     'top_shape_trajectories_zscore_genes.png',
                     densities=(dens_healthy, dens_aki), specimen_curves=_gs_z)
plot_trajectory_grid(gt_names, gt_h_raw, gt_a_raw, 'Log-normalized expression',
                     'Top shape-ranked gene trajectories (raw log-normalized, per-specimen overlay)',
                     'top_shape_trajectories_rawlognorm_genes.png',
                     specimen_curves=gene_specimen_curves)

# pathways: top-N by shape_rms; z-score curves reuse module fit
path_top = np.argsort(path_ls['shape_rms'])[::-1][:N_TOP]
pt_labels = [f"{retained_pathways.loc[i, 'library'].split('_')[0]}: {retained_pathways.loc[i, 'pathway']}"
             for i in path_top]
pt_h_z = path_ls['curve_healthy'][path_top]
pt_a_z = path_ls['curve_aki'][path_top]
path_specimen_curves = {}
for _spec, _code in specimen_codes.items():
    _m = samples == _spec
    _curves, _ = fit_single_condition_curves(
        module_scores[_m][:, path_top], s[_m], knots, grid, SECTION4_CONFIG['lambda_grid'],
        path_ls['lam_idx'][path_top], support_pct=SECTION4_CONFIG['common_support_pct'])
    path_specimen_curves[_spec] = (_curves, _code)
plot_trajectory_grid(pt_labels, pt_h_z, pt_a_z, 'Mean gene z-score',
                     'Top shape-ranked pathway trajectories (module z-score, per-specimen overlay)',
                     'top_shape_trajectories_zscore_pathways.png',
                     specimen_curves=path_specimen_curves)

# %%
# ----------------------------------------------------------------------------
# Between-specimen spread vs the condition gap -- the number the 4-curve figure shows visually.
#
# separation_t is a SPECIMEN-level two-sample t statistic, computed on the 4 specimen curves
# (not the ~20,000 tubules): the condition gap divided by the pooled between-specimen SD. This
# is the honest unit of replication for a condition contrast in this design.
#
# Calibration on simulated 2v2 data (12 seeds, 4000 cells, 4 specimens):
#     no true condition effect -> separation_t median 1.18, 90th pct 2.69, max 2.78
#     true shape effect        -> separation_t median 31.6, min 20.3
# The reference line is the two-sided t critical value at df = 2, i.e. 4.303: it flagged 0/12
# null features and 12/12 real ones. Descriptive, not a p-value -- averaging |t| over the grid
# is not itself t-distributed.
# NOTE: the pooled SD uses ddof=1. With n=2 per group, ddof=0 understates the spread by sqrt(2)
# and inflates every ratio.
# ----------------------------------------------------------------------------
from scipy.stats import t as _t_dist

SEPARATION_T_REFERENCE = float(_t_dist.ppf(0.975, df=2))          # 4.303


def specimen_separation(specimen_curves, labels):
    rows = []
    _h = [v for v, cd in specimen_curves.values() if cd == 0]
    _a = [v for v, cd in specimen_curves.values() if cd == 1]
    n_h, n_a = len(_h), len(_a)
    if n_h < 2 or n_a < 2:
        raise ValueError(f'specimen_separation needs >=2 specimens per condition, got '
                         f'{n_h} healthy and {n_a} aki.')
    for i, label in enumerate(labels):
        hs = np.vstack([x[i] for x in _h])
        as_ = np.vstack([x[i] for x in _a])
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', category=RuntimeWarning)
            gap = np.nanmean(as_, axis=0) - np.nanmean(hs, axis=0)
            var_h = np.nanvar(hs, axis=0, ddof=1)
            var_a = np.nanvar(as_, axis=0, ddof=1)
            pooled_sd = np.sqrt(((n_h - 1) * var_h + (n_a - 1) * var_a) / (n_h + n_a - 2))
            se = pooled_sd * np.sqrt(1.0 / n_h + 1.0 / n_a)
            t_grid = np.abs(gap) / np.where(se > 0, se, np.nan)
            mean_t = float(np.nanmean(t_grid))
        rows.append({'feature': str(label),
                     'mean_condition_gap': float(np.nanmean(np.abs(gap))),
                     'mean_between_specimen_sd': float(np.nanmean(pooled_sd)),
                     'separation_t': mean_t,
                     'clears_df2_reference': bool(mean_t >= SEPARATION_T_REFERENCE)})
    return pd.DataFrame(rows)


gene_separation = specimen_separation(gene_specimen_curves, gt_names)
gene_separation.to_csv(HEALTHY_VS_AKI_OUTPUT_DIR / 'specimen_separation_genes.csv', index=False)
display(gene_separation.round(3))
_weak = gene_separation[~gene_separation['clears_df2_reference']]
if len(_weak):
    print(f'\nWARNING: {len(_weak)}/{len(gene_separation)} top shape genes fall below the '
          f'specimen-level reference (separation_t < {SEPARATION_T_REFERENCE:.2f}):')
    print('   ' + ', '.join(_weak['feature'].tolist()))
    print('  Their condition gap is not large relative to the variation between animals of the '
          'SAME condition, so they are not interpretable as condition effects at n=2 vs 2 -- '
          'regardless of how high shape_rms is.')

# pathways raw log-normalized: refit M2 on raw mean-lognorm module for the top-N (for plotting only)
raw_top = np.empty((adata_pt.n_obs, len(path_top)))
for jj, pi in enumerate(path_top):
    idxs = [gene_local[g_] for g_ in retained_pathways.loc[pi, 'genes_present']]
    raw_top[:, jj] = Y_csc[:, idxs].toarray().mean(axis=1)
raw_fit = run_level_shape(raw_top, s, c, knots, grid, SECTION4_CONFIG['lambda_grid'])
plot_trajectory_grid(pt_labels, raw_fit['curve_healthy'], raw_fit['curve_aki'],
                     'Mean raw log-normalized expression',
                     'Top shape-ranked pathway trajectories (raw log-normalized)',
                     'top_shape_trajectories_rawlognorm_pathways.png')

# canonical S1/S2/S3 marker validation (per condition), reusing gene_ls fitted curves
fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharex=True)
for ax, group in zip(axes, ['S1', 'S2', 'S3']):
    for gene in resolve_present(MARKER_GROUPS[group], gene_names):
        gi = gene_local[gene]
        ax.plot(grid, gene_ls['curve_healthy'][gi], color=COL_H, lw=2)
        ax.plot(grid, gene_ls['curve_aki'][gi], color=COL_A, lw=2, ls='--')
        ax.text(grid[-1], gene_ls['curve_healthy'][gi][-1], f' {gene}', fontsize=8, va='center')
    ax.set_title(f'{group} markers')
    ax.set_xlabel('Shared pseudospace')
    ax.set_ylabel('GAM-fitted log-normalized expression')
axes[0].plot([], [], color=COL_H, lw=2, label='healthy')
axes[0].plot([], [], color=COL_A, lw=2, ls='--', label='aki')
axes[0].legend(fontsize=8)
fig.suptitle('Canonical PT S1/S2/S3 marker validation', y=1.02)
_save(fig, 'marker_validation_S1S2S3.png')
print('Group 1 figures written:', len(figure_index))


# %%
# ----------------------------------------------------------------------------
# Group 2 -- level/shape decomposition views (the point of the method)
# ----------------------------------------------------------------------------
def level_shape_scatter(df, name, title):
    fig, ax = plt.subplots(figsize=(7.5, 6))
    x = df['level_effect'].to_numpy()
    y = df['shape_rms'].to_numpy()
    color = -np.log10(np.clip(df['sample_perm_p'].to_numpy(), 1e-6, 1))
    sc_ = ax.scatter(x, y, c=color, cmap='viridis', s=18, alpha=0.7)
    fig.colorbar(sc_, ax=ax, label='-log10(sample_perm_p)')
    y_ref = np.quantile(y, 0.98)
    x_ref = np.quantile(np.abs(x), 0.98)
    ax.axhline(y_ref, color='grey', ls=':', lw=1)
    ax.axvline(x_ref, color='grey', ls=':', lw=1)
    ax.axvline(-x_ref, color='grey', ls=':', lw=1)
    ax.text(0.98, 0.02, 'high |x|, low y = pure level shift (DE-like)\nhigh y = shape change',
            transform=ax.transAxes, ha='right', va='bottom', fontsize=8,
            bbox=dict(boxstyle='round', fc='white', ec='grey', alpha=0.8))
    label_col = 'gene' if 'gene' in df.columns else 'pathway'
    for _, row in df.sort_values('shape_rms', ascending=False).head(N_TOP).iterrows():
        ax.annotate(str(row[label_col])[:26], (row['level_effect'], row['shape_rms']),
                    fontsize=7, xytext=(3, 3), textcoords='offset points')
    ax.set_xlabel('level_effect  (mean healthy->aki vertical gap)')
    ax.set_ylabel('shape_rms  (offset-removed curve difference)')
    ax.set_title(title + '\nDE would surface only the x-axis (level) hits', fontsize=10)
    _save(fig, name)


level_shape_scatter(gene_results, 'level_vs_shape_scatter_genes.png', 'Genes: level vs shape')
level_shape_scatter(pathway_results, 'level_vs_shape_scatter_pathways.png', 'Pathways: level vs shape')


def offset_removed_grid(labels, curve_h, curve_a, name, suptitle):
    fig, axes = _grid_axes(len(labels))
    for ax, i in zip(axes, range(len(labels))):
        diff = curve_a[i] - curve_h[i]
        ax.plot(grid, diff - diff.mean(), color='firebrick', lw=2)
        ax.axhline(0, color='grey', lw=1, ls='--')
        ax.set_title(str(labels[i])[:48], fontsize=9)
        ax.set_xlabel('Shared pseudospace')
        ax.set_ylabel('(aki - healthy) - mean gap')
    for ax in axes[len(labels):]:
        ax.axis('off')
    fig.suptitle(suptitle, y=1.005, fontsize=12)
    _save(fig, name)


offset_removed_grid(gt_names, gt_h_raw, gt_a_raw, 'offset_removed_curves_genes.png',
                    'Offset-removed difference curves (top shape genes)')
offset_removed_grid(pt_labels, pt_h_z, pt_a_z, 'offset_removed_curves_pathways.png',
                    'Offset-removed difference curves (top shape pathways)')


def delta_pseudospace_heatmap(labels, mat_cols, name, suptitle):
    n_bins = SECTION4_CONFIG['heatmap_n_bins']
    sigma = SECTION4_CONFIG['heatmap_smooth_sigma']
    prof_h, centers = binned_matrix_profile(mat_cols[c == 0], s[c == 0], lo, hi, n_bins, sigma)
    prof_a, _ = binned_matrix_profile(mat_cols[c == 1], s[c == 1], lo, hi, n_bins, sigma)
    delta = prof_a - prof_h
    order = np.argsort(np.argmax(np.abs(delta), axis=1))          # order rows by where |delta| peaks
    delta = delta[order]
    labels = [labels[i] for i in order]
    # NB: `x or 1.0` does NOT work here -- NaN is truthy, so an all-NaN delta would give
    # vmin=vmax=nan and a blank panel.
    vmax = float(np.nanpercentile(np.abs(delta), 99)) if np.isfinite(delta).any() else 1.0
    if not np.isfinite(vmax) or vmax == 0:
        vmax = 1.0
    strip, _ = binned_matrix_profile(pt_s1s3_axis[:, None], s, lo, hi, n_bins, sigma)
    fig, (ax_s, ax) = plt.subplots(2, 1, figsize=(11, 0.32 * len(labels) + 2.0),
                                   gridspec_kw={'height_ratios': [1, max(len(labels), 6)]},
                                   sharex=True)
    ax_s.imshow(strip, aspect='auto', cmap='PuOr', extent=[lo, hi, 0, 1])
    ax_s.set_yticks([])
    ax_s.set_title('S1 -> S3 marker axis (PuOr)', fontsize=8)
    im = ax.imshow(delta, aspect='auto', cmap=SECTION4_CONFIG['diverging_cmap'],
                   vmin=-vmax, vmax=vmax, extent=[lo, hi, len(labels), 0])
    ax.set_yticks(np.arange(len(labels)) + 0.5)
    ax.set_yticklabels([str(x)[:40] for x in labels], fontsize=7)
    ax.set_xlabel('Shared pseudospace')
    fig.colorbar(im, ax=ax, label='aki - healthy (binned)')
    fig.suptitle(suptitle, y=1.01, fontsize=12)
    _save(fig, name)


delta_pseudospace_heatmap(list(gt_names), Y_csc[:, gene_top].toarray(),
                          'delta_pseudospace_heatmap_genes.png',
                          'Delta (aki - healthy) pseudospace heatmap: top shape genes')
delta_pseudospace_heatmap(list(pt_labels), module_scores[:, path_top],
                          'delta_pseudospace_heatmap_pathways.png',
                          'Delta (aki - healthy) pseudospace heatmap: top shape pathways')


def paired_condition_heatmaps(labels, curve_h, curve_a, name, suptitle):
    """Paired condition heatmaps, twice: a common scale and a shape-only view.

    Standardising each condition on its own (the previous behaviour) removes that condition's level
    and amplitude, so a collapsed AKI gradient can look as strong as the control gradient. The two
    views answer different questions and must be read together:
      * ``<name>`` - both conditions on ONE colour scale in fitted lognorm units: level and amplitude
        survive, so a weaker AKI gradient reads as weaker;
      * ``<name>_shape`` - each condition standardised on its own: patterns are comparable in shape
        only, and colour is NOT an expression amount.
    """
    stem, dot, extension = str(name).rpartition('.')
    stem = stem if dot else str(name)
    extension = extension if dot else 'png'

    common_vmin = float(np.nanmin([np.nanmin(curve_h), np.nanmin(curve_a)]))
    common_vmax = float(np.nanmax([np.nanmax(curve_h), np.nanmax(curve_a)]))
    fig, axes = plt.subplots(2, 1, figsize=(11, 0.32 * len(labels) + 2.5), sharex=True)
    for ax, mat, title in [(axes[0], curve_h, 'healthy'), (axes[1], curve_a, 'aki')]:
        im = ax.imshow(mat, aspect='auto', cmap='viridis', vmin=common_vmin, vmax=common_vmax,
                       extent=[lo, hi, len(labels), 0])
        ax.set_yticks(np.arange(len(labels)) + 0.5)
        ax.set_yticklabels([str(x)[:40] for x in labels], fontsize=7)
        ax.set_ylabel(title)
        fig.colorbar(im, ax=ax, label='fitted log-normalised expression')
    axes[1].set_xlabel('Shared pseudospace')
    fig.suptitle(suptitle + '\ncommon expression scale: level and amplitude preserved',
                 y=1.01, fontsize=12)
    _save(fig, f'{stem}.{extension}')

    zh = zscore_rows(curve_h)
    za = zscore_rows(curve_a)
    fig, axes = plt.subplots(2, 1, figsize=(11, 0.32 * len(labels) + 2.5), sharex=True)
    for ax, mat, title in [(axes[0], zh, 'healthy'), (axes[1], za, 'aki')]:
        im = ax.imshow(mat, aspect='auto', cmap='bwr', vmin=-Z_CLIP, vmax=Z_CLIP,
                       extent=[lo, hi, len(labels), 0])
        ax.set_yticks(np.arange(len(labels)) + 0.5)
        ax.set_yticklabels([str(x)[:40] for x in labels], fontsize=7)
        ax.set_ylabel(title)
        fig.colorbar(im, ax=ax, label='per-row z-score')
    axes[1].set_xlabel('Shared pseudospace')
    fig.suptitle(suptitle + '\nshape only: each condition standardised on its own '
                            '(colour is not an expression amount)', y=1.01, fontsize=12)
    _save(fig, f'{stem}_shape.{extension}')


paired_condition_heatmaps(list(gt_names), gt_h_raw, gt_a_raw,
                          'paired_condition_heatmaps_genes.png',
                          'Paired condition heatmaps (top shape genes, shared order)')
paired_condition_heatmaps(list(pt_labels), pt_h_z, pt_a_z,
                          'paired_condition_heatmaps_pathways.png',
                          'Paired condition heatmaps (top shape pathways, shared order)')
print('Group 2 figures written; total so far:', len(figure_index))


# %%
# ----------------------------------------------------------------------------
# Group 3 -- summary / diagnostic views
# ----------------------------------------------------------------------------
def volcano(df, name, title):
    fig, ax = plt.subplots(figsize=(7.5, 6))
    x = df['shape_rms'].to_numpy()
    y = -np.log10(np.clip(df['sample_perm_p'].to_numpy(), 1e-6, 1))
    ax.scatter(x, y, s=16, alpha=0.6, color='slategrey')
    label_col = 'gene' if 'gene' in df.columns else 'pathway'
    for _, row in df.head(N_TOP).iterrows():
        ax.annotate(str(row[label_col])[:24], (row['shape_rms'], -np.log10(max(row['sample_perm_p'], 1e-6))),
                    fontsize=7, xytext=(3, 3), textcoords='offset points')
    ax.set_xlabel('shape_rms (effect size)')
    ax.set_ylabel('-log10(sample_perm_p)')
    ax.text(0.02, 0.98, 'CAVEAT: 2v2 design -> exact permutation p-floor is coarse (2/6);\n'
                        'this axis is a calibration check, not a powered test. Rank by shape_rms.',
            transform=ax.transAxes, ha='left', va='top', fontsize=8,
            bbox=dict(boxstyle='round', fc='lightyellow', ec='grey', alpha=0.9))
    ax.set_title(title)
    _save(fig, name)


volcano(gene_results, 'volcano_shape_genes.png', 'Genes: shape effect vs permutation p')
volcano(pathway_results, 'volcano_shape_pathways.png', 'Pathways: shape effect vs permutation p')


def rank_barplot(df, name, title):
    top = df.head(N_TOP).iloc[::-1]
    label_col = 'gene' if 'gene' in df.columns else 'pathway'
    big_level = np.abs(top['level_effect']) >= np.median(np.abs(df['level_effect']))
    colors = ['darkorange' if b else 'steelblue' for b in big_level]
    fig, ax = plt.subplots(figsize=(8, 0.42 * len(top) + 1.5))
    ax.barh([str(x)[:44] for x in top[label_col]], top['shape_rms'], color=colors)
    ax.set_xlabel('shape_rms')
    ax.set_title(title)
    ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color='darkorange'),
                       plt.Rectangle((0, 0), 1, 1, color='steelblue')],
              labels=['also large level shift', 'shape-dominant'], fontsize=8, loc='lower right')
    _save(fig, name)


rank_barplot(gene_results, 'rank_barplot_shape_genes.png', 'Top shape-ranked genes')
rank_barplot(pathway_results, 'rank_barplot_shape_pathways.png', 'Top shape-ranked pathways')

# pathway shape summary matrix: pathway x [level_effect, shape_rms, shape_F_cellwise, sample_perm_p]
summary_cols = ['level_effect', 'shape_rms', 'shape_F_cellwise_uncalibrated', 'sample_perm_p']
top_paths = pathway_results.head(N_TOP)
mat = top_paths[summary_cols].to_numpy(dtype=float)
norm = np.zeros_like(mat)
for j in range(mat.shape[1]):
    col = mat[:, j]
    rng_ = np.ptp(col)
    norm[:, j] = (col - col.min()) / rng_ if rng_ > 0 else 0.0
fig, ax = plt.subplots(figsize=(7.5, 0.45 * len(top_paths) + 1.8))
im = ax.imshow(norm, aspect='auto', cmap='magma')
ax.set_xticks(range(len(summary_cols)))
ax.set_xticklabels(summary_cols, rotation=30, ha='right')
ax.set_yticks(range(len(top_paths)))
ax.set_yticklabels([str(p)[:44] for p in top_paths['pathway']], fontsize=7)
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        ax.text(j, i, f'{mat[i, j]:.2g}', ha='center', va='center', fontsize=7,
                color='white' if norm[i, j] < 0.6 else 'black')
ax.set_title('Top pathway shape summary (color = column min-max normalized)')
_save(fig, 'pathway_shape_summary_matrix.png')
print('Group 3 figures written; total:', len(figure_index))

# %% [markdown]
# ## 4.9 - sensitivity: does the result survive on a PHYSICAL, injury-independent axis?
#
# `total_scanpy_dpt` is built from segment-identity markers, and AKI degrades those markers
# (Section 4.6b). The strongest available check is therefore to re-run the whole level/shape
# decomposition against a coordinate that injury cannot move: physical depth in the tissue.
#
# Depth proxy: glomeruli exist only in the cortex, so the mean distance to the k nearest
# glomeruli is a smooth cortex -> medulla coordinate. It is preferred over distance to the single
# nearest glomerulus, which is noisy in the convoluted cortical labyrinth (a DCT can sit closer
# to a glomerulus than its own PT). Rescaled per specimen, because section size, orientation and
# the cut plane differ between specimens and the GAM basis is defined on [0, 1].
#
# The headline number is the Spearman correlation between the two coordinates' `shape_rms`
# rankings. High = the finding is a property of the tissue, not of the marker-derived axis.

# %%
# ----------------------------------------------------------------------------
# Section 4.9 -- physical-axis sensitivity analysis
# ----------------------------------------------------------------------------
GLOM_DENSITY_K = 10
N_TOP_PATH_SENSITIVITY = 60      # how many top-DPT pathways to retest and report on
PATH_SENSITIVITY_FAIL_RANK = 200 # physical-axis rank past which a top-DPT pathway is a failure
CURVE_SPEARMAN_AMP_RATIO = 0.5   # shape_rms/min(amplitude) above this -> curve_spearman is noise

_pass1 = sc.read_h5ad(HARMONY_OUTPUT_PATH)          # all cells, incl. the Glomerulus class
_pt_names = adata_pt.obs_names
_depth = np.full(adata_pt.n_obs, np.nan)
_depth_report = []

for _sample in sorted(set(samples)):
    _m = samples == _sample
    _p1 = _pass1.obs[_pass1.obs['sample'].astype(str) == _sample]
    _glom = _p1[_p1['coarse_class'].astype(str) == 'Glomerulus']
    if len(_glom) < 3:
        _depth_report.append({'sample': _sample, 'n_glomeruli': len(_glom), 'k_used': 0,
                              'status': 'SKIPPED (<3 glomeruli)'})
        continue
    _k = min(GLOM_DENSITY_K, len(_glom))
    _tree = cKDTree(_glom[['x_centroid', 'y_centroid']].to_numpy())
    _xy = _pass1.obs.loc[_pt_names[_m], ['x_centroid', 'y_centroid']].to_numpy()
    _d, _ = _tree.query(_xy, k=_k)
    _raw = _d.mean(axis=1) if _k > 1 else _d.ravel()     # larger = further from cortex = deeper
    _p1v, _p99v = np.percentile(_raw, [1, 99])
    _depth[_m] = np.clip((_raw - _p1v) / max(_p99v - _p1v, 1e-8), 0, 1)
    _depth_report.append({'sample': _sample, 'n_glomeruli': len(_glom), 'k_used': _k,
                          'status': 'ok', 'raw_p1': float(_p1v), 'raw_p99': float(_p99v)})

display(pd.DataFrame(_depth_report))
del _pass1

_ok = np.isfinite(_depth)
if _ok.sum() < 0.5 * adata_pt.n_obs:
    print(f'SKIPPING the physical-axis sensitivity analysis: depth is defined for only '
          f'{int(_ok.sum()):,}/{adata_pt.n_obs:,} PT tubules. Too few specimens carry enough '
          'glomeruli for the proxy to be meaningful.')
    physical_sensitivity = None
else:
    print(f'Depth proxy defined for {int(_ok.sum()):,}/{adata_pt.n_obs:,} PT tubules.')
    print(f'Spearman(total_scanpy_dpt, physical depth) = {safe_spearman(s[_ok], _depth[_ok]):+.3f}'
          '   <- independent support for the DPT axis; low values are the finding, not a bug.')

    _dlo = max(np.percentile(_depth[_ok & (c == 0)], p_lo), np.percentile(_depth[_ok & (c == 1)], p_lo))
    _dhi = min(np.percentile(_depth[_ok & (c == 0)], p_hi), np.percentile(_depth[_ok & (c == 1)], p_hi))
    _dgrid = np.linspace(_dlo, _dhi, SECTION4_CONFIG['grid_points'])
    _dknots = gam_internal_knots(_depth[_ok], basis_df=3 + SECTION4_CONFIG['n_internal_knots'])

    phys_ls = run_level_shape(Y_genes[_ok], _depth[_ok], c[_ok], _dknots, _dgrid,
                              SECTION4_CONFIG['lambda_grid'])
    phys_shape_p, phys_level_p, _, _ = sample_perm_pvalues(
        Y_genes[_ok], _depth[_ok], samples[_ok], _dknots, _dgrid, SECTION4_CONFIG['lambda_grid'],
        phys_ls['lam_idx'], phys_ls['sl2'], phys_ls['p_b'], control_samples)

    physical_sensitivity = pd.DataFrame({
        'gene': gene_names,
        'shape_rms_dpt_axis': gene_ls['shape_rms'],
        'shape_rms_physical_axis': phys_ls['shape_rms'],
        'level_effect_dpt_axis': gene_ls['level_effect'],
        'level_effect_physical_axis': phys_ls['level_effect'],
        'sample_perm_p_physical_axis': phys_shape_p,
        'level_perm_p_physical_axis': phys_level_p,
    })
    physical_sensitivity['rank_dpt'] = (-physical_sensitivity['shape_rms_dpt_axis']).rank(method='min')
    physical_sensitivity['rank_physical'] = (-physical_sensitivity['shape_rms_physical_axis']).rank(method='min')
    physical_sensitivity = physical_sensitivity.sort_values('rank_dpt').reset_index(drop=True)
    physical_sensitivity.to_csv(
        HEALTHY_VS_AKI_OUTPUT_DIR / 'physical_axis_sensitivity_genes.csv', index=False)

    # --- same check for the pathway module scores (Section 4.5) -------------------
    # Pathway shape transfers MUCH worse than gene shape (rho ~0.34 vs ~0.70): module
    # scores average away the per-gene gradients, so what survives on the DPT axis is
    # often just where injured cells land in DPT space. Amplitudes are exported here
    # (unlike in pathway_level_shape_results.csv) because curve_spearman is
    # uninterpretable when shape_rms is a large fraction of the curve amplitude.
    phys_path_ls = run_level_shape(module_scores[_ok], _depth[_ok], c[_ok], _dknots, _dgrid,
                                   SECTION4_CONFIG['lambda_grid'])
    phys_path_shape_p, phys_path_level_p, _, _ = sample_perm_pvalues(
        module_scores[_ok], _depth[_ok], samples[_ok], _dknots, _dgrid,
        SECTION4_CONFIG['lambda_grid'], phys_path_ls['lam_idx'], phys_path_ls['sl2'],
        phys_path_ls['p_b'], control_samples)

    physical_sensitivity_paths = pd.DataFrame({
        'library': retained_pathways['library'],
        'pathway': retained_pathways['pathway'],
        'n_genes_present': retained_pathways['n_genes_present'],
        'shape_rms_dpt_axis': path_ls['shape_rms'],
        'shape_rms_physical_axis': phys_path_ls['shape_rms'],
        'level_effect_dpt_axis': path_ls['level_effect'],
        'level_effect_physical_axis': phys_path_ls['level_effect'],
        'curve_spearman_dpt_axis': path_ls['curve_spearman'],
        'curve_spearman_physical_axis': phys_path_ls['curve_spearman'],
        'amplitude_healthy_dpt': path_ls['amplitude_healthy'],
        'amplitude_aki_dpt': path_ls['amplitude_aki'],
        'amplitude_healthy_physical': phys_path_ls['amplitude_healthy'],
        'amplitude_aki_physical': phys_path_ls['amplitude_aki'],
        'sample_perm_p_physical_axis': phys_path_shape_p,
        'level_perm_p_physical_axis': phys_path_level_p,
    })
    physical_sensitivity_paths['rank_dpt'] = (
        -physical_sensitivity_paths['shape_rms_dpt_axis']).rank(method='min')
    physical_sensitivity_paths['rank_physical'] = (
        -physical_sensitivity_paths['shape_rms_physical_axis']).rank(method='min')
    physical_sensitivity_paths = physical_sensitivity_paths.sort_values(
        'rank_dpt').reset_index(drop=True)
    physical_sensitivity_paths.to_csv(
        HEALTHY_VS_AKI_OUTPUT_DIR / 'physical_axis_sensitivity_pathways.csv', index=False)

    _rho_path = safe_spearman(physical_sensitivity_paths['shape_rms_dpt_axis'].to_numpy(),
                              physical_sensitivity_paths['shape_rms_physical_axis'].to_numpy())

    _rho = safe_spearman(physical_sensitivity['shape_rms_dpt_axis'].to_numpy(),
                         physical_sensitivity['shape_rms_physical_axis'].to_numpy())
    _top_dpt = set(gene_names[np.argsort(gene_ls['shape_rms'])[::-1][:N_TOP]])
    _top_phys = set(gene_names[np.argsort(phys_ls['shape_rms'])[::-1][:N_TOP]])
    _overlap = _top_dpt & _top_phys

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
    axes[0].scatter(s[_ok], _depth[_ok], s=3, alpha=0.15, linewidths=0, rasterized=True,
                    c=np.where(c[_ok] == 0, 0, 1), cmap=ListedColormap([COL_H, COL_A]))
    axes[0].set_xlabel('total_scanpy_dpt (marker-derived)')
    axes[0].set_ylabel('physical depth proxy (glomerular density)')
    axes[0].set_title(f'Coordinate agreement (Spearman={safe_spearman(s[_ok], _depth[_ok]):+.2f})')
    axes[1].scatter(physical_sensitivity['shape_rms_dpt_axis'],
                    physical_sensitivity['shape_rms_physical_axis'], s=10, alpha=0.5,
                    color='slategrey')
    for _, r in physical_sensitivity.head(N_TOP).iterrows():
        axes[1].annotate(str(r['gene'])[:18], (r['shape_rms_dpt_axis'], r['shape_rms_physical_axis']),
                         fontsize=7, xytext=(3, 3), textcoords='offset points')
        axes[1].scatter([r['shape_rms_dpt_axis']], [r['shape_rms_physical_axis']], s=26,
                        color='firebrick', zorder=3)
    axes[1].set_xlabel('shape_rms on the DPT axis')
    axes[1].set_ylabel('shape_rms on the physical axis')
    axes[1].set_title(f'Shape effect agreement (Spearman={_rho:+.2f})')
    fig.suptitle('Section 4.10: does the shape result survive an injury-independent coordinate?',
                 y=1.02, fontsize=12)
    _save(fig, 'physical_axis_sensitivity.png')

    print(f'\nRank concordance of shape_rms across coordinates: Spearman = {_rho:+.3f}')
    print(f'Top-{N_TOP} overlap: {len(_overlap)}/{N_TOP}'
          + (f'  -> {sorted(_overlap)}' if _overlap else ''))
    if _rho < 0.3 or len(_overlap) < N_TOP // 3:
        print('WARNING: the shape ranking does NOT transfer to a physical coordinate. The top '
              'hits are properties of the marker-derived axis, not of position in the tissue. '
              'Report them as such, or prefer the physical-axis result.')
    else:
        print('The shape ranking largely transfers to a physical coordinate -- the finding is not '
              'an artifact of building the axis out of injury-sensitive identity markers.')

    # --- pathway arm: concordance + which top-DPT pathways fail the coordinate swap ---
    _p_top = physical_sensitivity_paths.head(N_TOP_PATH_SENSITIVITY)
    _p_failed = _p_top[_p_top['rank_physical'] > PATH_SENSITIVITY_FAIL_RANK]
    print(f'\nPathway modules -- shape_rms concordance across coordinates: '
          f'Spearman = {_rho_path:+.3f}   (genes: {_rho:+.3f})')
    print(f'level_effect concordance: Spearman = '
          f"{safe_spearman(physical_sensitivity_paths['level_effect_dpt_axis'].to_numpy(), physical_sensitivity_paths['level_effect_physical_axis'].to_numpy()):+.3f}"
          '   <- level is near coordinate-invariant; only SHAPE is at risk.')
    print(f'{len(_p_failed)}/{len(_p_top)} of the top-{N_TOP_PATH_SENSITIVITY} DPT pathways fall '
          f'below rank {PATH_SENSITIVITY_FAIL_RANK} on the physical axis:')
    for _, _r in _p_failed.iterrows():
        print(f"   rank {int(_r['rank_dpt']):>3} -> {int(_r['rank_physical']):>4}   "
              f"shape {_r['shape_rms_dpt_axis']:.3f} -> {_r['shape_rms_physical_axis']:.3f}   "
              f"{str(_r['pathway'])[:62]}")
    if _rho_path < _rho:
        print('\nNOTE: pathway shape transfers WORSE than gene shape. Module scores average away '
              'the per-gene gradients, so a pathway can rank highly on the DPT axis largely '
              'because injury repositions cells along that axis. Prefer the physical-axis '
              'ranking for pathway-level claims.')

    # curve_spearman is a rank correlation between two fitted curves; when shape_rms is a
    # large fraction of the smaller curve amplitude, the curve is too flat for the rank
    # correlation to mean anything. Flag those rows rather than letting them be read as
    # profile "inversions".
    _amp_min = physical_sensitivity_paths[['amplitude_healthy_dpt', 'amplitude_aki_dpt']].min(axis=1)
    _shape_over_amp = physical_sensitivity_paths['shape_rms_dpt_axis'] / np.maximum(_amp_min, 1e-9)
    _unreliable = (_shape_over_amp > CURVE_SPEARMAN_AMP_RATIO).sum()
    print(f'\ncurve_spearman reliability: {_unreliable}/{len(physical_sensitivity_paths)} pathway '
          f'modules have shape_rms > {CURVE_SPEARMAN_AMP_RATIO:g}x the smaller curve amplitude. '
          'curve_spearman is NOT interpretable for those rows.')


# %% [markdown]
# ## 4.10 - summary

# %%
# ----------------------------------------------------------------------------
# Section 4.10 -- summary + figure index
# ----------------------------------------------------------------------------
print('=' * 78)
print('SECTION 4 SUMMARY -- mouse-only healthy vs AKI (PT tubules)')
print('=' * 78)
print(f'PT tubules: healthy={n_healthy:,}  aki={n_aki:,}  (samples {control_samples} vs {aki_samples})')
print(f'Genes tested: {len(gene_names):,}   Pathway modules tested: {P:,}')
print(f'Common support: [{lo:.3f}, {hi:.3f}]   internal knots: {len(knots)}')
print('\nTop 10 genes by shape_rms:')
print(gene_results.head(10)[['gene', 'level_effect', 'shape_rms', 'sample_perm_p']].to_string(index=False))
print('\nTop 10 pathways by shape_rms:')
print(pathway_results.head(10)[['library', 'pathway', 'level_effect', 'shape_rms', 'sample_perm_p']].to_string(index=False))
print('\n--- Robustness checks (read these BEFORE the gene list above) ---')
_weak_n = int((~gene_separation['clears_df2_reference']).sum())
print(f'Between-specimen separation: {len(gene_separation) - _weak_n}/{len(gene_separation)} top '
      f'shape genes clear the specimen-level reference (separation_t >= '
      f'{SEPARATION_T_REFERENCE:.2f}), i.e. have a condition gap exceeding the spread between '
      'same-condition specimens.')
if _weak_n:
    print(f'  {_weak_n} do NOT clear that bar and are not interpretable as condition effects.')

if injury_scores:
    for _panel, _score in injury_scores.items():
        _hi = _score >= np.nanquantile(_score, 0.90)
        print(f'Injury confounding ({_panel}): top-decile tubules sit at median pseudospace '
              f'{np.nanmedian(s[_hi]):.3f} vs {np.nanmedian(s):.3f} overall; '
              f'rho(score, pseudospace) = {safe_spearman(_score, s):+.3f}.')
    print('  Shape hits whose signal coincides with the injury-high stretch are confounded.')

if physical_sensitivity is not None:
    _rho_final = safe_spearman(physical_sensitivity['shape_rms_dpt_axis'].to_numpy(),
                               physical_sensitivity['shape_rms_physical_axis'].to_numpy())
    _rho_final_path = safe_spearman(
        physical_sensitivity_paths['shape_rms_dpt_axis'].to_numpy(),
        physical_sensitivity_paths['shape_rms_physical_axis'].to_numpy())
    print(f'Physical-axis sensitivity: shape_rms rank concordance Spearman = {_rho_final:+.3f} '
          f'(genes) and {_rho_final_path:+.3f} (pathway modules) between the marker-derived DPT '
          'axis and the glomerular-density depth proxy. Pathway-level shape claims rest on the '
          'weaker of the two -- prefer the physical-axis ranking for them.')
else:
    print('Physical-axis sensitivity: NOT RUN (too few glomeruli for the depth proxy).')

print('\nCOHORT DISCLAIMER: 2 healthy + 2 AKI specimens. Results are DESCRIPTIVE and')
print('EFFECT-SIZE-RANKED (shape_rms). A formal 2v2 significance test is not achievable')
print('(exact permutation p-floor 2/6); sample_perm_p is a calibration/sanity check only.')
print('Note also that *_F_cellwise_uncalibrated in the CSVs are CELL-level F statistics: their')
print('effective n is 4 specimens, not the tubule count. Do not read them as tests.')
print('\nFigures written to', str(HEALTHY_VS_AKI_OUTPUT_DIR.relative_to(PROJECT_DIR)) + ':')
for name in figure_index:
    print('  -', name)
csvs = ['gene_level_shape_results.csv', 'pathway_level_shape_results.csv', 'pseudospace_structure_diagnostics.csv']
if physical_sensitivity is not None:
    csvs += ['physical_axis_sensitivity_genes.csv', 'physical_axis_sensitivity_pathways.csv']
print('CSVs written:', ', '.join(csvs))

# %% [markdown]
# ## 4.11 - collaborator shortlist: robust S1/S2/S3-peaking PT gene programs
#
# **Question.** Which reproducible PT expression programs peak in S1, S2 or S3, and which show
# the largest descriptive Control-versus-AKI trajectory remodeling?
#
# This is an exploratory prioritization layer for kidney collaborators, not a new significance
# analysis. It balances sample-consistent segment specificity, `shape_rms`, and detection; it
# excludes axis-defining markers, known sex-biased PT genes, and mitochondrial/ribosomal genes
# from the primary shortlist. Similar sample-balanced trajectories are clustered *within* each
# peak segment. The exports include the full candidate table and plotted source data so every
# inclusion can be audited. All selection, clustering, plotting, and export code is defined in
# the notebook cells below; this section has no external analysis-helper dependency.

# %%
"""Build a collaborator-facing summary of spatially patterned PT genes.

The primary shortlist is exploratory.  It prioritizes genes with a reproducible S1/S2/S3 peak
across specimens and a large healthy-versus-AKI trajectory-shape effect, then clusters similar
sample-balanced expression profiles within each peak segment.  It does not treat tubules as
independent biological replicates and does not convert the 2-vs-2 design into significance claims.
"""

import os
import re
import sys
from pathlib import Path
from typing import Iterable

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.cluster.hierarchy import fcluster, leaves_list, linkage
from scipy.ndimage import gaussian_filter1d
from scipy.spatial.distance import pdist


SEGMENTS = ("PT-S1", "PT-S2", "PT-S3")
SEGMENT_COLORS = {"PT-S1": "#3B6FB6", "PT-S2": "#D88932", "PT-S3": "#3D9970"}
CONDITION_COLORS = {"Control": "#0072B2", "AKI": "#D55E00"}

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 7,
    "axes.titlesize": 8,
    "axes.labelsize": 7,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.7,
    "legend.frameon": False,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})


def _norm_gene(gene: str) -> str:
    return str(gene).strip().upper()


def _as_csr(matrix):
    return matrix.tocsr() if sparse.issparse(matrix) else sparse.csr_matrix(np.asarray(matrix))


def _dense_mean(matrix, mask: np.ndarray) -> np.ndarray:
    if not np.any(mask):
        return np.full(matrix.shape[1], np.nan)
    return np.asarray(matrix[mask].mean(axis=0)).ravel()


def _rank01(values: pd.Series) -> pd.Series:
    return values.rank(method="average", pct=True).fillna(0.0)


def _interpolate_and_smooth(curve: np.ndarray, sigma: float = 1.0) -> np.ndarray:
    x = np.arange(curve.size)
    finite = np.isfinite(curve)
    if finite.sum() == 0:
        return np.zeros_like(curve)
    if finite.sum() == 1:
        curve = np.full_like(curve, curve[finite][0])
    elif not finite.all():
        if not np.all(np.diff(x[finite]) > 0):
            raise ValueError("Interpolation coordinates must be strictly increasing.")
        curve = np.interp(x, x[finite], curve[finite])
    return gaussian_filter1d(curve, sigma=sigma, mode="nearest")


def _load_alignment_gate():
    """Load the render-time figure QA gate from its standard installation or an override."""
    override = os.environ.get("NATURE_FIGURE_QA_DIR")
    qa_dir = Path(override).expanduser() if override else (
        Path.home() / ".codex" / "skills" / "nature-figure" / "scripts"
    )
    if not (qa_dir / "audit_panel_alignment.py").exists():
        raise FileNotFoundError(
            "Figure alignment QA is unavailable. Set NATURE_FIGURE_QA_DIR to the directory "
            "containing audit_panel_alignment.py."
        )
    if str(qa_dir) not in sys.path:
        sys.path.insert(0, str(qa_dir))
    from audit_panel_alignment import require_matplotlib_panel_alignment

    return require_matplotlib_panel_alignment


def _export_figure(fig, stem: Path, *, axes, panel_ids, row_groups, exclude_axes=()):
    """Run alignment QA, then export editable and collaborator-preview formats."""
    require_matplotlib_panel_alignment = _load_alignment_gate()
    fig.canvas.draw()
    require_matplotlib_panel_alignment(
        fig,
        axes=axes,
        panel_ids=panel_ids,
        row_groups=row_groups,
        exclude_axes=exclude_axes,
        require_panel_labels=True,
        strict=True,
        tolerance_pt=1.5,
        gutter_tolerance_pt=1.5,
        json_out=stem.with_suffix(".alignment.json"),
        overlay_svg=stem.with_suffix(".alignment.svg"),
    )
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"})


def _cluster_within_segments(shortlist: pd.DataFrame, profiles: np.ndarray,
                             modules_per_segment: int) -> tuple[pd.DataFrame, list[int]]:
    """Cluster genes within their robust peak segment and return a stable display order."""
    shortlist = shortlist.copy()
    shortlist["module"] = ""
    display_order: list[int] = []
    for segment in SEGMENTS:
        positions = np.flatnonzero(shortlist["peak_segment"].to_numpy() == segment)
        if positions.size == 0:
            continue
        local = profiles[positions]
        if positions.size < 3:
            raw_clusters = np.ones(positions.size, dtype=int)
            leaves = np.arange(positions.size)
        else:
            distances = np.nan_to_num(pdist(local, metric="correlation"), nan=0.0)
            tree = linkage(distances, method="average", optimal_ordering=True)
            raw_clusters = fcluster(tree, t=min(modules_per_segment, positions.size),
                                    criterion="maxclust")
            leaves = leaves_list(tree)
            # A one- or two-gene branch is an outlier, not a collaborator-facing program. Keep
            # the hierarchical gene order but collapse such unstable splits to the segment module.
            cluster_sizes = np.bincount(raw_clusters)[1:]
            if modules_per_segment == 1 or np.any(cluster_sizes < 3):
                raw_clusters = np.ones(positions.size, dtype=int)

        cluster_ids = sorted(
            np.unique(raw_clusters),
            key=lambda cluster: np.argmax(np.mean(local[raw_clusters == cluster], axis=0)),
        )
        cluster_to_letter = {cluster: chr(ord("A") + i) for i, cluster in enumerate(cluster_ids)}
        for local_i, global_i in enumerate(positions):
            shortlist.loc[global_i, "module"] = f"{segment.replace('PT-', '')}-{cluster_to_letter[raw_clusters[local_i]]}"

        for cluster in cluster_ids:
            members = [i for i in leaves if raw_clusters[i] == cluster]
            display_order.extend(positions[members].tolist())
    return shortlist, display_order



# %%
def build_pt_collaborator_summary(
    dpt_path: str | Path,
    shape_results_path: str | Path,
    output_dir: str | Path,
    *,
    axis_basis_genes: Iterable[str] = (),
    known_sex_biased_genes: Iterable[str] = (),
    n_per_segment: int = 12,
    modules_per_segment: int = 2,
    n_bins: int = 48,
    minimum_detection_fraction: float = 0.10,
    minimum_peak_agreement: float = 0.75,
) -> dict[str, object]:
    """Create the PT program shortlist, source-data tables, and two presentation figures."""
    dpt_path = Path(dpt_path)
    shape_results_path = Path(shape_results_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(dpt_path)
    pt_mask = adata.obs["broad_tubule_marker_call"].astype(str).eq("PT").to_numpy()
    pt = adata[pt_mask].copy()
    if pt.n_obs == 0:
        raise ValueError("The DPT file contains no PT observations.")
    if "lognorm" not in pt.layers:
        raise KeyError("The DPT file must contain layers['lognorm'] for expression summaries.")

    segment = pt.obs["segment_class"].astype(str).to_numpy()
    unexpected = sorted(set(segment) - set(SEGMENTS))
    if unexpected:
        raise ValueError(f"PT subset contains unexpected segment labels: {unexpected}")
    sample = pt.obs["sample"].astype(str).to_numpy()
    condition_raw = pt.obs["condition"].astype(str).str.lower().to_numpy()
    condition = np.where(condition_raw == "ir", "AKI", "Control")
    position = pt.obs["total_scanpy_dpt"].to_numpy(dtype=float)
    if not np.all(np.isfinite(position)):
        raise ValueError("PT DPT positions must all be finite.")

    shape = pd.read_csv(shape_results_path).drop_duplicates("gene").set_index("gene")
    var_lookup = {str(gene): i for i, gene in enumerate(pt.var_names)}
    genes = np.array([gene for gene in shape.index.astype(str) if gene in var_lookup], dtype=object)
    if genes.size == 0:
        raise ValueError("No genes overlap between the DPT file and shape-results table.")
    indices = np.array([var_lookup[gene] for gene in genes], dtype=int)
    expression = _as_csr(pt.layers["lognorm"][:, indices])
    detection_fraction = np.asarray((expression > 0).sum(axis=0)).ravel() / pt.n_obs

    sample_names = sorted(set(sample))
    balanced = np.full((len(sample_names), len(SEGMENTS), genes.size), np.nan)
    for sample_i, sample_name in enumerate(sample_names):
        for segment_i, segment_name in enumerate(SEGMENTS):
            mask = (sample == sample_name) & (segment == segment_name)
            balanced[sample_i, segment_i] = _dense_mean(expression, mask)
    segment_means = np.nanmean(balanced, axis=0)
    if np.any(~np.isfinite(segment_means)):
        raise ValueError("At least one sample/segment combination is empty; balanced peak calls fail.")

    segment_sd = segment_means.std(axis=0)
    segment_sd[segment_sd == 0] = 1.0
    segment_z = (segment_means - segment_means.mean(axis=0)) / segment_sd
    peak_index = np.argmax(segment_means, axis=0)
    peak_segment = np.array(SEGMENTS, dtype=object)[peak_index]
    sorted_z = np.sort(segment_z, axis=0)
    peak_margin = sorted_z[-1] - sorted_z[-2]
    per_sample_peak = np.argmax(balanced, axis=1)
    peak_agreement = np.mean(per_sample_peak == peak_index[None, :], axis=0)

    metrics = pd.DataFrame({
        "gene": genes,
        "peak_segment": peak_segment,
        "peak_agreement": peak_agreement,
        "segment_specificity_z_margin": peak_margin,
        "detection_fraction": detection_fraction,
    })
    for segment_i, segment_name in enumerate(SEGMENTS):
        metrics[f"balanced_mean_{segment_name}"] = segment_means[segment_i]
    metrics = metrics.join(shape.reset_index().set_index("gene"), on="gene", rsuffix="_shape")

    axis_set = {_norm_gene(gene) for gene in axis_basis_genes}
    sex_set = {_norm_gene(gene) for gene in known_sex_biased_genes}
    gene_upper = metrics["gene"].map(_norm_gene)
    metrics["axis_basis_gene"] = gene_upper.isin(axis_set)
    metrics["known_sex_biased_gene"] = gene_upper.isin(sex_set)
    metrics["technical_gene"] = metrics["gene"].str.match(r"^(mt-|Rpl|Rps)", case=False)
    metrics["specificity_percentile"] = _rank01(metrics["segment_specificity_z_margin"])
    metrics["shape_percentile"] = _rank01(metrics["shape_rms"])
    metrics["detection_percentile"] = _rank01(metrics["detection_fraction"])
    metrics["priority_score"] = (
        0.55 * metrics["specificity_percentile"]
        + 0.35 * metrics["shape_percentile"]
        + 0.10 * metrics["detection_percentile"]
    )
    metrics["primary_eligible"] = (
        (metrics["detection_fraction"] >= minimum_detection_fraction)
        & (metrics["peak_agreement"] >= minimum_peak_agreement)
        & ~metrics["axis_basis_gene"]
        & ~metrics["known_sex_biased_gene"]
        & ~metrics["technical_gene"]
    )

    chosen = []
    for segment_name in SEGMENTS:
        group = metrics[metrics["primary_eligible"] & metrics["peak_segment"].eq(segment_name)]
        chosen.append(group.nlargest(n_per_segment, "priority_score"))
    shortlist = pd.concat(chosen, ignore_index=True)
    if shortlist.empty:
        raise ValueError("No genes passed the collaborator-shortlist filters.")

    selected_lookup = {gene: i for i, gene in enumerate(genes)}
    selected_columns = np.array([selected_lookup[gene] for gene in shortlist["gene"]], dtype=int)
    selected_expression = expression[:, selected_columns].tocsr()
    gene_mean = np.asarray(selected_expression.mean(axis=0)).ravel()
    gene_sq = np.asarray(selected_expression.multiply(selected_expression).mean(axis=0)).ravel()
    gene_std = np.sqrt(np.maximum(gene_sq - gene_mean ** 2, 1e-12))

    # The figure compares specimens only over their shared p1-p99 PT support. Selection metrics
    # above still use every PT observation; this trim only prevents extrapolated visual curves.
    low = max(np.percentile(position[sample == name], 1) for name in sample_names)
    high = min(np.percentile(position[sample == name], 99) for name in sample_names)
    if not low < high:
        raise ValueError("Specimens have no shared p1-p99 PT DPT support.")
    edges = np.linspace(low, high, n_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    sample_curves = np.full((len(sample_names), len(shortlist), n_bins), np.nan)
    source_rows = []
    for sample_i, sample_name in enumerate(sample_names):
        for bin_i in range(n_bins):
            right = position <= edges[bin_i + 1] if bin_i == n_bins - 1 else position < edges[bin_i + 1]
            mask = (sample == sample_name) & (position >= edges[bin_i]) & right
            raw = _dense_mean(selected_expression, mask)
            sample_curves[sample_i, :, bin_i] = raw
        for gene_i in range(len(shortlist)):
            sample_curves[sample_i, gene_i] = _interpolate_and_smooth(
                sample_curves[sample_i, gene_i], sigma=1.0
            )

    z_curves = (sample_curves - gene_mean[None, :, None]) / gene_std[None, :, None]
    pooled_profile = np.nanmean(z_curves, axis=0)
    shortlist, display_order = _cluster_within_segments(
        shortlist, pooled_profile, modules_per_segment
    )
    shortlist["display_order"] = np.arange(len(shortlist))
    shortlist.loc[display_order, "display_order"] = np.arange(len(display_order))
    shortlist = shortlist.sort_values("display_order").reset_index(drop=True)

    # Reindex all trajectory arrays to the final clustered order.
    inverse = np.array(display_order, dtype=int)
    sample_curves = sample_curves[:, inverse]
    z_curves = z_curves[:, inverse]
    pooled_profile = pooled_profile[inverse]

    for sample_i, sample_name in enumerate(sample_names):
        sample_condition = condition[sample == sample_name][0]
        for gene_i, row in shortlist.iterrows():
            for bin_i, center in enumerate(centers):
                source_rows.append({
                    "sample": sample_name,
                    "condition": sample_condition,
                    "gene": row["gene"],
                    "peak_segment": row["peak_segment"],
                    "module": row["module"],
                    "dpt_bin_center": center,
                    "smoothed_log_normalized_mean": sample_curves[sample_i, gene_i, bin_i],
                    "smoothed_gene_z_score": z_curves[sample_i, gene_i, bin_i],
                })
    source_data = pd.DataFrame(source_rows)

    shortlist.to_csv(output_dir / "pt_gene_program_shortlist.csv", index=False)
    metrics.sort_values("priority_score", ascending=False).to_csv(
        output_dir / "pt_gene_program_all_candidates.csv", index=False
    )
    source_data.to_csv(output_dir / "pt_gene_program_binned_source_data.csv", index=False)

    # Figure 1: hero heatmap plus aligned evidence bars.
    fig = plt.figure(figsize=(7.2, 8.6), constrained_layout=True)
    grid_spec = fig.add_gridspec(1, 7, width_ratios=[1, 1, 1, 1, 1, 0.9, 0.9])
    ax_heat = fig.add_subplot(grid_spec[0, :5])
    ax_shape = fig.add_subplot(grid_spec[0, 5], sharey=ax_heat)
    ax_specificity = fig.add_subplot(grid_spec[0, 6], sharey=ax_heat)
    row_positions = np.arange(len(shortlist))
    image = ax_heat.imshow(
        pooled_profile, aspect="auto", interpolation="nearest", cmap="RdBu_r", vmin=-2.5, vmax=2.5,
        extent=[centers[0], centers[-1], len(shortlist) - 0.5, -0.5], rasterized=True,
    )
    module_start = ~shortlist["module"].duplicated()
    display_labels = [
        f"{module}  {gene}" if is_start else f"      {gene}"
        for gene, module, is_start in zip(shortlist["gene"], shortlist["module"], module_start)
    ]
    ax_heat.set_yticks(row_positions)
    ax_heat.set_yticklabels(display_labels)
    for label, peak in zip(ax_heat.get_yticklabels(), shortlist["peak_segment"]):
        label.set_color(SEGMENT_COLORS[peak])
    ax_heat.set_xlabel("Total DPT within shared PT support")
    ax_heat.set_ylabel("Exploratory gene shortlist")

    segment_medians = {name: float(np.median(position[segment == name])) for name in SEGMENTS}
    boundaries = [(segment_medians[SEGMENTS[i]] + segment_medians[SEGMENTS[i + 1]]) / 2
                  for i in range(2)]
    for boundary in boundaries:
        ax_heat.axvline(boundary, color="black", lw=0.6, ls="--", alpha=0.7)
    zone_edges = [centers[0], *boundaries, centers[-1]]
    for i, segment_name in enumerate(SEGMENTS):
        middle = 0.5 * (zone_edges[i] + zone_edges[i + 1])
        ax_heat.text(middle, -1.35, segment_name.replace("PT-", ""), ha="center", va="bottom",
                     color=SEGMENT_COLORS[segment_name], fontweight="bold", clip_on=False)

    changes = np.flatnonzero(shortlist["module"].to_numpy()[1:] != shortlist["module"].to_numpy()[:-1]) + 0.5
    for boundary in changes:
        for axis in (ax_heat, ax_shape, ax_specificity):
            axis.axhline(boundary, color="white" if axis is ax_heat else "0.75", lw=1.0)
    bar_colors = [SEGMENT_COLORS[value] for value in shortlist["peak_segment"]]
    ax_shape.barh(row_positions, shortlist["shape_rms"], color=bar_colors, height=0.72)
    ax_specificity.barh(row_positions, shortlist["segment_specificity_z_margin"],
                        color=bar_colors, height=0.72)
    ax_shape.set_title("AKI shape\neffect", pad=9)
    ax_specificity.set_title("Segment\nspecificity", pad=9)
    ax_shape.set_xlabel("RMS")
    ax_specificity.set_xlabel("z margin")
    for axis in (ax_shape, ax_specificity):
        axis.tick_params(axis="y", left=False, labelleft=False)
        axis.set_ylim(len(shortlist) - 0.5, -0.5)
        axis.grid(axis="x", color="0.9", lw=0.5)
        axis.set_axisbelow(True)

    colorbar_axis = ax_heat.inset_axes([0.02, -0.11, 0.36, 0.018])
    colorbar = fig.colorbar(image, cax=colorbar_axis, orientation="horizontal")
    colorbar.set_label("Gene-wise z-score", labelpad=1)
    colorbar.ax.tick_params(labelsize=6, pad=1)
    for panel, axis in zip("abc", (ax_heat, ax_shape, ax_specificity)):
        axis.text(-0.08, 1.025, panel, transform=axis.transAxes, fontsize=8,
                  fontweight="bold", va="bottom", ha="left")
    fig.suptitle(
        "Robust PT gene programs for collaborator review",
        fontsize=10, fontweight="bold",
    )
    heatmap_stem = output_dir / "pt_gene_program_heatmap"
    _export_figure(
        fig, heatmap_stem,
        axes=[ax_heat, ax_shape, ax_specificity],
        panel_ids=["a", "b", "c"],
        row_groups=[["a", "b", "c"]],
        exclude_axes=[colorbar_axis],
    )
    plt.show()
    plt.close(fig)

    # Figure 2: one curve panel per module, retaining all four specimen trajectories.
    modules = list(dict.fromkeys(shortlist["module"]))
    n_columns = len(modules) if len(modules) <= 3 else 2
    n_rows = int(np.ceil(len(modules) / n_columns))
    fig, axes = plt.subplots(
        n_rows, n_columns, figsize=(7.2, max(2.8, 2.25 * n_rows)),
        sharex=True, sharey=True, constrained_layout=True, squeeze=False,
    )
    flat_axes = axes.ravel()
    panel_ids = []
    module_rows = []
    for panel_i, (axis, module) in enumerate(zip(flat_axes, modules)):
        members = np.flatnonzero(shortlist["module"].to_numpy() == module)
        per_sample_module = np.nanmean(z_curves[:, members], axis=1)
        for sample_i, sample_name in enumerate(sample_names):
            sample_condition = condition[sample == sample_name][0]
            axis.plot(centers, per_sample_module[sample_i], ls="--", lw=0.8, alpha=0.62,
                      color=CONDITION_COLORS[sample_condition])
            for bin_i, center in enumerate(centers):
                module_rows.append({
                    "module": module,
                    "sample": sample_name,
                    "condition": sample_condition,
                    "dpt_bin_center": center,
                    "module_mean_gene_z_score": per_sample_module[sample_i, bin_i],
                })
        for condition_name in ("Control", "AKI"):
            sample_indices = [i for i, name in enumerate(sample_names)
                              if condition[sample == name][0] == condition_name]
            axis.plot(centers, np.nanmean(per_sample_module[sample_indices], axis=0),
                      lw=2.0, color=CONDITION_COLORS[condition_name], label=condition_name)
        for boundary in boundaries:
            axis.axvline(boundary, color="0.75", lw=0.55, ls=":")
        peak_segment_name = shortlist.loc[members[0], "peak_segment"]
        axis.set_title(f"{module} · {len(members)} genes", color=SEGMENT_COLORS[peak_segment_name])
        axis.axhline(0, color="0.8", lw=0.5)
        axis.set_xlabel("Total DPT")
        axis.set_ylabel("Module mean z-score")
        panel_id = chr(ord("a") + panel_i)
        panel_ids.append(panel_id)
        axis.text(-0.08, 1.03, panel_id, transform=axis.transAxes, fontsize=8,
                  fontweight="bold", va="bottom", ha="left")
    for axis in flat_axes[len(modules):]:
        axis.set_visible(False)
    flat_axes[0].text(0.98, 0.92, "Control", transform=flat_axes[0].transAxes,
                      color=CONDITION_COLORS["Control"], ha="right", fontweight="bold")
    flat_axes[0].text(0.98, 0.85, "AKI", transform=flat_axes[0].transAxes,
                      color=CONDITION_COLORS["AKI"], ha="right", fontweight="bold")
    fig.suptitle("PT gene programs across Control and AKI specimens\n"
                 "Dashed lines show individual specimens",
                 fontsize=9, fontweight="bold")
    curves_stem = output_dir / "pt_gene_program_module_curves"
    visible_axes = list(flat_axes[:len(modules)])
    row_groups = []
    for row_i in range(n_rows):
        group = panel_ids[row_i * n_columns:(row_i + 1) * n_columns]
        if len(group) > 1:
            row_groups.append(group)
    _export_figure(
        fig, curves_stem,
        axes=visible_axes,
        panel_ids=panel_ids,
        row_groups=row_groups,
    )
    plt.show()
    plt.close(fig)
    module_source = pd.DataFrame(module_rows)
    module_source.to_csv(output_dir / "pt_gene_program_module_source_data.csv", index=False)

    notes = (
        "PT collaborator gene-program summary\n"
        f"PT observations used for selection: {pt.n_obs:,}\n"
        f"Specimens: {len(sample_names)} ({', '.join(sample_names)})\n"
        f"Genes selected: {len(shortlist)} ({n_per_segment} requested per segment)\n"
        f"Curve display support: p1-p99 intersection across specimens [{low:.4f}, {high:.4f}]\n"
        f"Primary filters: detection >= {minimum_detection_fraction:.0%}; peak agreement >= "
        f"{minimum_peak_agreement:.0%}; exclude axis-basis, known sex-biased, mitochondrial, "
        "and ribosomal genes.\n"
        "Priority score: 55% segment-specificity percentile + 35% AKI shape-effect percentile "
        "+ 10% detection percentile.\n"
        "Inference boundary: exploratory effect-size ranking from 2 Control and 2 AKI specimens; "
        "not a significance screen. Selection and visualization use the same cohort.\n"
    )
    (output_dir / "README.txt").write_text(notes)
    print(notes)
    print("Modules:")
    print(shortlist.groupby(["peak_segment", "module"], sort=False)["gene"]
          .apply(lambda values: ", ".join(values)).to_string())

    return {
        "shortlist": shortlist,
        "all_candidates": metrics,
        "binned_source_data": source_data,
        "module_source_data": module_source,
        "output_dir": output_dir,
    }



# %%
PT_COLLABORATOR_DIR = RESULTS_DIR / 'pt_collaborator_summary'
pt_collaborator_results = build_pt_collaborator_summary(
    DPT_OUTPUT_PATH,
    HEALTHY_VS_AKI_OUTPUT_DIR / 'gene_level_shape_results.csv',
    PT_COLLABORATOR_DIR,
    axis_basis_genes=ADJ_AXIS_BASIS_GENES,
    known_sex_biased_genes=SEX_BIASED_PT_GENES,
    n_per_segment=12,
    modules_per_segment=2,
)
display(pt_collaborator_results['shortlist'][[
    'gene', 'peak_segment', 'module', 'peak_agreement',
    'segment_specificity_z_margin', 'shape_rms', 'level_effect', 'priority_score',
]])

# %% [markdown]
# ## 4.12 - Curve modules, signed rankings and pathway redundancy (mouse-only)
#
# The same shape-first machinery the cross-species workflow uses, applied to the control-vs-AKI
# contrast:
#
# * **positional** modules cluster the **control** curves: which genes share a normal PT pattern,
#   discovered before the injury effect is consulted, so a large injury effect is not a requirement for
#   a positional program.
# * **response** modules cluster the mean-centred **difference** curves, which separates broad
#   suppression, selective loss of a late program, early induction and peak relocation.
# * stability under specimen omission (recomputed without each specimen) and under a different data
#   mixture (pooled fit vs specimen-balanced curves).
# * module-pathway overlap tested against the genes eligible for module discovery, BH-corrected across
#   all pairs. Five signed rankings (level, amplitude, redistribution, early, late) are tested
#   competitively as well, with a correlation-aware p so a pathway is not credited for one correlated
#   block.
# * the prioritized pathways are checked for redundancy and for member-gene support, so a claim that
#   rests on a single gene, or on a pathway that is another pathway's member set, is visible.
#

# %%
# Purpose: curve modules, signed rankings and pathway redundancy for the mouse contrast.
from pseudospace.enrichment import camera_like_enrichment, signed_gene_rankings
from pseudospace.modules import (
    difference_curves,
    discover_curve_modules,
    enrich_modules,
    module_stability,
)
from pseudospace.specimen import specimen_balanced_curves

# Per-specimen curves for every tested gene: a pooled fit cannot be split per specimen because the
# condition indicator is constant inside one specimen.
specimen_full_curves = {}
for specimen in sorted(set(samples)):
    mask = samples == specimen
    curves, _ = fit_single_condition_curves(
        Y_genes[mask], s[mask], knots, grid, SECTION4_CONFIG['lambda_grid'],
        gene_ls['lam_idx'], support_pct=SECTION4_CONFIG['common_support_pct'],
    )
    specimen_full_curves[specimen] = curves
balanced_healthy = specimen_balanced_curves(
    {name: curves for name, curves in specimen_full_curves.items() if name in control_samples})
balanced_aki = specimen_balanced_curves(
    {name: curves for name, curves in specimen_full_curves.items() if name in aki_samples})

MODULE_CONFIG = {'max_distance': 0.4, 'min_amplitude': 0.05, 'min_features': 10}
positional_labels, positional_modules = discover_curve_modules(
    gene_ls['curve_healthy'], grid, max_distance=MODULE_CONFIG['max_distance'],
    min_amplitude=MODULE_CONFIG['min_amplitude'], min_features=MODULE_CONFIG['min_features'],
    feature_names=gene_names,
)
response_delta = difference_curves(gene_ls['curve_healthy'], gene_ls['curve_aki'])
response_labels, response_modules = discover_curve_modules(
    difference_curves(gene_ls['curve_healthy'], gene_ls['curve_aki'], center=True), grid,
    max_distance=MODULE_CONFIG['max_distance'], min_amplitude=MODULE_CONFIG['min_amplitude'],
    min_features=MODULE_CONFIG['min_features'], feature_names=gene_names,
)
mouse_gene_modules = pd.DataFrame({
    'gene': gene_names,
    'positional_module': positional_labels.to_numpy(),
    'response_module': response_labels.to_numpy(),
    'delta_mean_aki_minus_healthy': np.nanmean(response_delta, axis=1),
})
mouse_gene_modules.to_csv(HEALTHY_VS_AKI_OUTPUT_DIR / 'module_gene_assignment.csv', index=False)
positional_modules.to_csv(HEALTHY_VS_AKI_OUTPUT_DIR / 'module_positional_catalog.csv', index=False)
response_modules.to_csv(HEALTHY_VS_AKI_OUTPUT_DIR / 'module_response_catalog.csv', index=False)
display(positional_modules)
display(response_modules)


def _response_modules(curve_reference, curve_comparison):
    labels, _ = discover_curve_modules(
        difference_curves(curve_reference, curve_comparison, center=True), grid,
        max_distance=MODULE_CONFIG['max_distance'], min_amplitude=MODULE_CONFIG['min_amplitude'],
        min_features=MODULE_CONFIG['min_features'], feature_names=gene_names,
    )
    return labels


stability_runs = {'pooled_fit': response_labels}
for dropped in sorted(specimen_full_curves):
    kept = [name for name in specimen_full_curves if name != dropped]
    kept_control = [name for name in kept if name in control_samples]
    kept_aki = [name for name in kept if name in aki_samples]
    if not kept_control or not kept_aki:
        continue
    stability_runs[f'without_{dropped}'] = _response_modules(
        specimen_balanced_curves({name: specimen_full_curves[name] for name in kept_control}),
        specimen_balanced_curves({name: specimen_full_curves[name] for name in kept_aki}),
    )
stability_runs['specimen_balanced_mixture'] = _response_modules(balanced_healthy, balanced_aki)
response_stability = module_stability(stability_runs)
response_stability.to_csv(HEALTHY_VS_AKI_OUTPUT_DIR / 'module_response_stability.csv', index=False)
display(response_stability.round(3))

pathway_gene_sets = {
    f'{row.library}: {row.pathway}': row.genes_present
    for row in retained_pathways.itertuples()
}
module_enrichment = enrich_modules(
    {module: members.index
     for module, members in response_labels.groupby(response_labels, observed=True)
     if module != 'unassigned'},
    pathway_gene_sets,
    background=gene_names,
)
module_enrichment.to_csv(HEALTHY_VS_AKI_OUTPUT_DIR / 'module_pathway_enrichment.csv', index=False)
print('Response modules:', int(response_labels[response_labels != 'unassigned'].nunique()),
      'covering', int((response_labels != 'unassigned').sum()), 'genes; module-pathway pairs:',
      len(module_enrichment), '; significant after BH:',
      int((module_enrichment['p_value_adjusted'] < 0.05).sum()))
display(module_enrichment.head(10)[[
    'module', 'gene_set', 'n_module_genes', 'n_overlap', 'expected_overlap',
    'p_value', 'p_value_adjusted',
]].round(4))

signed_rankings = signed_gene_rankings(
    gene_ls['curve_healthy'], gene_ls['curve_aki'], grid, gene_names=gene_names
).merge(gene_results[['gene', 'level_fraction', 'shape_fraction', 'pattern_rms_z',
                      'difference_type']], on='gene', how='left')
signed_rankings.to_csv(HEALTHY_VS_AKI_OUTPUT_DIR / 'gene_signed_rankings.csv', index=False)
SIGNED_QUESTIONS = {
    'level_shift': 'level_effect',
    'amplitude_change': 'amplitude_log2_ratio',
    'redistribution': 'redistribution',
    'early_contrast': 'early_delta',
    'late_contrast': 'late_delta',
}
signed_enrichment_summary = []
for question, column in SIGNED_QUESTIONS.items():
    statistics = signed_rankings[['gene', column]].dropna().set_index('gene')[column]
    table = camera_like_enrichment(
        statistics, pathway_gene_sets, expression=Y_genes, gene_names=gene_names,
        background=gene_names, n_permutations=500,
        min_set_size=SECTION4_CONFIG['pathway_min_genes'],
    )
    table.insert(0, 'question', question)
    table.to_csv(HEALTHY_VS_AKI_OUTPUT_DIR / f'signed_enrichment_{question}.csv', index=False)
    signed_enrichment_summary.append({
        'question': question,
        'n_genes_ranked': int(len(statistics)),
        'n_sets_tested': int(table['p_value_permutation'].notna().sum()),
        'n_significant_after_BH': int((table['p_value_permutation_adjusted'] < 0.05).sum()),
        'median_correlation_inflation': float(table['correlation_inflation'].median()),
    })
signed_enrichment_summary = pd.DataFrame(signed_enrichment_summary)
signed_enrichment_summary.to_csv(
    HEALTHY_VS_AKI_OUTPUT_DIR / 'signed_enrichment_summary.csv', index=False)
display(signed_enrichment_summary)

# Redundancy and member-gene support for the prioritized pathways the visualization layer ranks.
prioritized_pathways = pathway_results.head(40).copy()
prioritized_pathways['member_list'] = prioritized_pathways['genes_present'].str.split('; ').map(
    lambda members: [member for member in members if member])
redundancy_pairs, redundancy_groups = summarize_pathway_redundancy(
    prioritized_pathways, gene_column='member_list', label_columns=('library', 'pathway'),
    overlap_threshold=0.6,
)
redundancy_pairs.to_csv(HEALTHY_VS_AKI_OUTPUT_DIR / 'pathway_redundancy_pairs.csv', index=False)
redundancy_groups.to_csv(HEALTHY_VS_AKI_OUTPUT_DIR / 'pathway_redundancy_groups.csv', index=False)
print(f'Overlapping prioritized pathway pairs: {len(redundancy_pairs)}; largest group: '
      f'{int(redundancy_groups["group_size"].max()) if len(redundancy_groups) else 0}')
if len(redundancy_pairs):
    display(redundancy_pairs.sort_values('n_shared', ascending=False).head(10))

member_evidence = []
for row in prioritized_pathways.head(20).itertuples():
    evidence = member_gene_evidence(row.member_list, gene_results, effect_column='level_effect')
    member_evidence.append({'library': row.library, 'pathway': row.pathway,
                            'shape_rms': row.shape_rms, **evidence.to_dict()})
member_evidence = pd.DataFrame(member_evidence)
member_evidence.to_csv(HEALTHY_VS_AKI_OUTPUT_DIR / 'pathway_top_member_evidence.csv', index=False)
display(member_evidence[[
    'pathway', 'n_members_present', 'fraction_members_agreeing', 'strongest_member',
    'member_effect_mean', 'member_effect_mean_without_strongest', 'sign_flips_without_strongest',
]].round(3))


# %% [markdown]
# ## 4.13 - Magnitude and matched-position checks (mouse-only)
#
# "AKI changed this gene" is also a magnitude claim, so the same checks the cross-species workflow
# needs apply here: detection and abundance per condition, expression ratio against abundance,
# within-condition early-to-late gradients (which cancel a constant gene-specific offset), a
# pseudospace-matched contrast computed per specimen pair, and a capture summary.
#

# %%
# Purpose: magnitude and matched-position checks for the mouse contrast.
from pseudospace.specimen import pseudobulk_profiles

condition_is_aki = (c == 1)
detection_abundance = pd.DataFrame({
    'gene': gene_names,
    'detected_fraction_healthy': np.asarray((Y_genes[~condition_is_aki] > 0).mean(axis=0)).ravel(),
    'detected_fraction_aki': np.asarray((Y_genes[condition_is_aki] > 0).mean(axis=0)).ravel(),
    'mean_lognorm_healthy': np.asarray(Y_genes[~condition_is_aki].mean(axis=0)).ravel(),
    'mean_lognorm_aki': np.asarray(Y_genes[condition_is_aki].mean(axis=0)).ravel(),
})
detection_abundance['log2_ratio_aki_over_healthy'] = np.log2(
    (detection_abundance['mean_lognorm_aki'] + 1e-3)
    / (detection_abundance['mean_lognorm_healthy'] + 1e-3)
)
detection_abundance['mean_abundance'] = 0.5 * (
    detection_abundance['mean_lognorm_healthy'] + detection_abundance['mean_lognorm_aki']
)
detection_abundance.to_csv(
    HEALTHY_VS_AKI_OUTPUT_DIR / 'condition_detection_abundance.csv', index=False)
abundance_bins = pd.qcut(detection_abundance['mean_abundance'], 10, duplicates='drop')
ratio_by_abundance = detection_abundance.groupby(abundance_bins, observed=True)[
    'log2_ratio_aki_over_healthy'
].agg(['size', 'median', 'mean', 'std']).reset_index()
ratio_by_abundance['abundance_midpoint'] = ratio_by_abundance['mean_abundance'].apply(
    lambda interval: float(interval.mid))
ratio_by_abundance.to_csv(HEALTHY_VS_AKI_OUTPUT_DIR / 'ratio_vs_abundance_bins.csv', index=False)
fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
axes[0].scatter(detection_abundance['mean_abundance'],
                detection_abundance['log2_ratio_aki_over_healthy'], s=3, alpha=0.25,
                color='#4C72B0', linewidths=0)
axes[0].axhline(0.0, color='black', lw=0.8)
axes[0].plot(ratio_by_abundance['abundance_midpoint'], ratio_by_abundance['median'],
             color='crimson', lw=2, marker='o', ms=4, label='median per abundance decile')
axes[0].legend(frameon=False, fontsize=8)
axes[0].set_xlabel('Mean fitted lognorm expression (both conditions)')
axes[0].set_ylabel('log2(AKI / healthy)')
axes[0].set_title('Ratio versus abundance')
axes[1].hist(detection_abundance['log2_ratio_aki_over_healthy'].dropna(), bins=80,
             color='#4C72B0')
axes[1].axvline(0.0, color='black', lw=0.8)
axes[1].set_xlabel('log2(AKI / healthy)')
axes[1].set_title('Distribution of mean-expression ratios')
fig.suptitle('Magnitude checks: is the ratio abundance-dependent?')
_save(fig, 'ratio_vs_abundance.png')

early_third = grid <= np.quantile(grid, 1 / 3)
late_third = grid >= np.quantile(grid, 2 / 3)
within_condition_gradients = pd.DataFrame({
    'gene': gene_names,
    'healthy_early_to_late': np.nanmean(balanced_healthy[:, late_third], axis=1)
        - np.nanmean(balanced_healthy[:, early_third], axis=1),
    'aki_early_to_late': np.nanmean(balanced_aki[:, late_third], axis=1)
        - np.nanmean(balanced_aki[:, early_third], axis=1),
})
within_condition_gradients['gradient_difference_aki_minus_healthy'] = (
    within_condition_gradients['aki_early_to_late']
    - within_condition_gradients['healthy_early_to_late']
)
within_condition_gradients = within_condition_gradients.merge(
    gene_results[['gene', 'level_effect', 'difference_type']], on='gene', how='left')
within_condition_gradients.to_csv(
    HEALTHY_VS_AKI_OUTPUT_DIR / 'within_condition_gradients.csv', index=False)
finite_gradients = within_condition_gradients.dropna(
    subset=['healthy_early_to_late', 'aki_early_to_late'])
print('Within-condition early-to-late gradients: Spearman(healthy, aki) =',
      round(spearmanr(finite_gradients['healthy_early_to_late'],
                      finite_gradients['aki_early_to_late']).correlation, 3),
      '; opposite-sign gradients:',
      int((np.sign(finite_gradients['healthy_early_to_late'])
           != np.sign(finite_gradients['aki_early_to_late'])).sum()),
      'of', len(finite_gradients))

# Pseudospace-matched, per-specimen contrast: each AKI specimen against each control specimen inside
# the same pseudospace bin, on within-row fractions, so the units are specimens.
pseudobulk = pseudobulk_profiles(
    adata_pt.layers['counts'], samples, s, n_bins=20, gene_names=list(adata_pt.var_names)
)
pseudobulk.to_csv(HEALTHY_VS_AKI_OUTPUT_DIR / 'pseudobulk_profiles.csv', index=False)
gene_columns = [column for column in pseudobulk if column.startswith('count_')]
counts_block = pseudobulk[gene_columns].to_numpy(dtype=float)
fractions = pd.DataFrame(
    counts_block / np.maximum(counts_block.sum(axis=1), 1.0)[:, None],
    columns=[column.replace('count_', '') for column in gene_columns],
)
fractions['specimen'] = pseudobulk['specimen'].to_numpy()
fractions['pseudospace_bin'] = pseudobulk['pseudospace_bin'].to_numpy()
pair_rows = []
for bin_index, block in fractions.groupby('pseudospace_bin', observed=True):
    aki_rows = block[block['specimen'].isin(aki_samples)]
    control_rows = block[block['specimen'].isin(control_samples)]
    if aki_rows.empty or control_rows.empty:
        continue
    genes = [column for column in block.columns if column not in ('specimen', 'pseudospace_bin')]
    aki_values = aki_rows[genes].to_numpy(dtype=float)
    control_values = control_rows[genes].to_numpy(dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        ratios = np.log2((aki_values[:, None, :] + 1e-6) / (control_values[None, :, :] + 1e-6))
    for aki_index, aki_name in enumerate(aki_rows['specimen']):
        for control_index, control_name in enumerate(control_rows['specimen']):
            pair_rows.append(pd.DataFrame({
                'bin': int(bin_index),
                'aki_specimen': aki_name,
                'control_specimen': control_name,
                'gene': genes,
                'log2_ratio': ratios[aki_index, control_index],
            }))
if pair_rows:
    matched_ratios = pd.concat(pair_rows, ignore_index=True)
    matched_ratios.to_csv(
        HEALTHY_VS_AKI_OUTPUT_DIR / 'matched_region_specimen_ratios.csv', index=False)
    per_gene_matched = matched_ratios.groupby('gene')['log2_ratio'].agg(
        ['median', 'std', 'size']).reset_index()
    per_gene_matched['n_specimen_pairs'] = matched_ratios[
        ['aki_specimen', 'control_specimen']].drop_duplicates().shape[0]
    per_gene_matched.to_csv(HEALTHY_VS_AKI_OUTPUT_DIR / 'matched_region_per_gene.csv', index=False)
    display(per_gene_matched.sort_values('median', ascending=False).head(10).round(3))
    print('Matched-position contrast rows:', len(matched_ratios), '; specimen pairs:',
          int(per_gene_matched['n_specimen_pairs'].iloc[0]))
else:
    print('Matched-position contrast skipped: no bin holds both an AKI and a control specimen.')

capture_columns = [column for column in ('n_genes_by_counts', 'total_counts', 'n_spots',
                                         'pct_counts_mt')
                   if column in adata_pt.obs.columns]
capture_summary = adata_pt.obs.groupby(['condition', 'sample'], observed=True)[
    capture_columns].median().reset_index()
capture_summary.to_csv(HEALTHY_VS_AKI_OUTPUT_DIR / 'capture_summary.csv', index=False)
display(capture_summary)


# %% [markdown]
# # Section 5 - three-axis concordance (non-circular fidelity check)
#
# Per specimen, agreement among the DPT pseudospace, the early->late marker axis, and physical
# distance-to-nearest-glomerulus (from the pass-1 centroids + `coarse_class`). Agreement across
# three independent axes is the strongest fidelity evidence achievable at n=2 vs 2.

# %%
dpt_ad = sc.read_h5ad(DPT_OUTPUT_PATH)
if 'lognorm' in dpt_ad.layers:
    dpt_ad.X = dpt_ad.layers['lognorm'].copy()
pass1 = sc.read_h5ad(HARMONY_OUTPUT_PATH)   # all cells: coarse_class + x/y centroids


def _score(adata_obj, genes, role):
    """Mean marker score for one end of the axis. random_state is required: score_genes draws a
    random control gene set, so without it this axis (and the whole concordance table) changes
    between runs."""
    present = [g for g in genes if g in adata_obj.var_names]
    if not present:
        raise ValueError(f'None of the {role} axis markers {genes} are in {DPT_OUTPUT_PATH.name}')
    sc.tl.score_genes(adata_obj, present, score_name='_tmp', use_raw=False,
                      random_state=RANDOM_STATE)
    values = adata_obj.obs['_tmp'].to_numpy(dtype=float)
    del adata_obj.obs['_tmp']
    return values


dpt_ad.obs['total_marker_axis'] = (_score(dpt_ad, TOTAL_POSITION_MARKERS['late'], 'late')
                                   - _score(dpt_ad, TOTAL_POSITION_MARKERS['early'], 'early'))

def _sp(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    return float(spearmanr(a[ok], b[ok]).correlation) if ok.sum() > 10 else np.nan

rows = []
for sample in sorted(dpt_ad.obs['sample'].astype(str).unique()):
    m = (dpt_ad.obs['sample'].astype(str) == sample).to_numpy()
    obs = dpt_ad.obs[m]
    p1 = pass1.obs[pass1.obs['sample'].astype(str) == sample]
    glom = p1[p1['coarse_class'].astype(str) == 'Glomerulus']
    if len(glom) < 3:
        print(f'{sample}: <3 glomeruli, skipping'); continue
    tree = cKDTree(glom[['x_centroid', 'y_centroid']].to_numpy())
    xy = pass1.obs.loc[obs.index, ['x_centroid', 'y_centroid']].to_numpy()
    dist, _ = tree.query(xy, k=1)
    dpt = obs['total_scanpy_dpt'].to_numpy(dtype=float)
    ma = obs['total_marker_axis'].to_numpy(dtype=float)
    rows.append({'sample': sample, 'n': int(m.sum()),
                 'rho_dpt_vs_dist_glom': _sp(dpt, dist),
                 'rho_marker_vs_dist_glom': _sp(ma, dist),
                 'rho_dpt_vs_marker': _sp(dpt, ma)})
conc = pd.DataFrame(rows)
conc.to_csv(CONCORDANCE_DIR / 'three_axis_concordance.csv', index=False)
display(conc)
fig, ax = plt.subplots(figsize=(7, 4))
conc.set_index('sample')[['rho_dpt_vs_dist_glom', 'rho_marker_vs_dist_glom', 'rho_dpt_vs_marker']].plot(kind='bar', ax=ax)
ax.axhline(0, color='k', lw=0.8); ax.set_ylabel('Spearman rho')
ax.set_title('Three-axis concordance per specimen')
fig.savefig(CONCORDANCE_DIR / 'three_axis_concordance.png', dpi=130, bbox_inches='tight'); plt.show()

# %% [markdown]
# # Section 6 - QuPath export and spatial validation (opt-in)
#
# Both steps live in `analysis/scripts/` and accept explicit data/results roots, so the
# module constants are pointed at this run's outputs before `main()` is called. Nothing here
# runs unless you flip the flags: the export writes GeoJSON, and the spatial validation needs a
# round of manual curation in QuPath between its `--dry-run` and its real run.

# %%
RUN_QUPATH_EXPORT = True
RUN_SPATIAL_VALIDATION = True
SPATIAL_VALIDATION_DRY_RUN = True      # seed the manual-review GeoJSON, then curate in QuPath

import importlib.util


def _load_script_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


if RUN_QUPATH_EXPORT:
    sys.path.insert(0, str(PROJECT_DIR / 'analysis' / 'scripts'))
    qupath_export = _load_script_module(
        'export_kept_tubules_to_qupath_geojson',
        PROJECT_DIR / 'analysis' / 'scripts' / 'export_kept_tubules_to_qupath_geojson.py')
    # The label source is the pass-1 file: it still carries Glomerulus and Vessel.
    qupath_export.HARMONY_RELATIVE_PATH = HARMONY_OUTPUT_PATH.relative_to(PROJECT_DIR)
    qupath_export.OUTPUT_RELATIVE_DIR = (RESULTS_DIR / 'qupath_geojson').relative_to(PROJECT_DIR)
    qupath_export.main([
        '--project-root', str(PROJECT_DIR),
        '--label-column', 'coarse_class',      # or 'segment_class' for the 17-way overlay
    ])

if RUN_SPATIAL_VALIDATION:
    sys.path.insert(0, str(PROJECT_DIR / 'analysis' / 'scripts'))
    spatial_validation = _load_script_module(
        'spatial_validation_of_pseudospace',
        PROJECT_DIR / 'analysis' / 'scripts' / 'spatial_validation_of_pseudospace.py')
    spatial_validation.HARMONY_RELATIVE_PATH = HARMONY_OUTPUT_PATH.relative_to(PROJECT_DIR)
    spatial_validation.DPT_RELATIVE_PATH = DPT_OUTPUT_PATH.relative_to(PROJECT_DIR)
    spatial_validation.OUTPUT_RELATIVE_DIR = (RESULTS_DIR / 'spatial_validation').relative_to(PROJECT_DIR)
    argv = ['--project-root', str(PROJECT_DIR)]
    if SPATIAL_VALIDATION_DRY_RUN:
        argv.append('--dry-run')
    spatial_validation.main(argv)

if not (RUN_QUPATH_EXPORT or RUN_SPATIAL_VALIDATION):
    print('Section 6 skipped. Set RUN_QUPATH_EXPORT / RUN_SPATIAL_VALIDATION to True to run it.')

# %%
PROJECT_DIR
