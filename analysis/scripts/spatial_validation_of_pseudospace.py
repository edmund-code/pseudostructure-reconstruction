"""Validate pseudospace against independent spatial axes: distance to nearest confidently-
labeled Glomerulus, and distance from tissue edge. Headline methods-validation figure for the
paper's claim that computational pseudospace recovers nephron anatomy from spatial
transcriptomics. Also cross-checks pseudospace-derived PT/TAL/DCT sub-segment boundaries against
an independent marker-panel argmax method (Step 8).

STEP 1 GLOMERULUS/VESSEL LABELING IS MANUALLY CURATED IN QUPATH, NOT AUTOMATED.
An automated statistical filter (score percentile + score-gap + celltype_primary thresholds) was
tried first and rejected: visual inspection in QuPath repeatedly showed it excluding real
glomeruli as "likely vessels," even after two rounds of loosening. The workflow is now:
  1. `--dry-run` seeds one GeoJSON per sample (results/mouse_only_v5/spatial_validation/
     {sample}_glomerulus_vessel_manual_review.geojson), restricted to cells in the Glomerulus
     DE-cluster (broad_tubule_marker_call == 'Glomerulus') whose ORIGINAL per-cell classification
     (celltype_primary, from notebook 6) is 'Glomerulus' or 'Vessel' (cells whose celltype_primary
     is neither are dropped). Labeled/colored via qupath_export.LABEL_COLORS -- the same
     Glomerulus/Vessel classes and colors already used in your existing full-annotation overlay --
     so they render as ordinary QuPath classifications, not a special QC scheme. The recomputed
     Glomerulus_broad_score-minus-Vessel_broad_score gap is attached as a
     `glom_vs_vessel_gap` property for reference in QuPath tooltips only; it no longer drives any
     automated decision. A sample already having a review file is never overwritten (manual edits
     are never clobbered by re-running --dry-run).
  2. You open each file in QuPath, correct any misclassified features by hand -- including
     renaming the class entirely (e.g. 'Glomerulus' -> 'Glomeruli' to match your segmentation
     vocabulary; 'Glomerulus'/'Glomeruli' are accepted interchangeably, case-insensitively,
     everywhere, and normalized internally to 'Glomeruli' as the canonical form for all output
     going forward) or reclassifying a feature entirely outside Glomerulus/Vessel scope (e.g. to
     'Tubules', 'Blood Cells' -- such features are simply dropped from the high-confidence set,
     not treated as an error) -- then save back to the exact same path (same filename).
  3. Running WITHOUT --dry-run reads the manually-curated GeoJSON back in. QuPath strips ALL
     custom properties on export (confirmed empirically: feature_index/source_sample/
     glom_vs_vessel_gap do not survive -- only objectType/name/classification/measurements do),
     so features are re-identified by nearest-centroid spatial match against the raw segmentation
     GeoJSON (verified exact, distance 0.0, and unique for real edited files -- geometry itself is
     preserved by QuPath even though properties are not). A manual annotation that doesn't match
     any raw polygon within tolerance (e.g. a brand-new hand-drawn shape) has no corresponding
     gene-expression data anywhere in this pipeline and fails loudly rather than being silently
     dropped or included. Before anything downstream runs, a plain-text diff comparing this
     manually-curated set against a fresh in-memory-only reconstruction of the automated seed
     (never written to disk) is saved to manual_annotation_diff.txt and printed to console. If a
     sample's review file doesn't exist yet at read-back time, it is auto-generated from
     celltype_primary on the spot and a loud warning is printed that those labels were never
     manually reviewed. Pass --use-manual-annotations on a future run to skip straight to reading
     pre-existing review files (no Harmony reload, no diff report) once you're done curating.

Read-only on all inputs (h5ads, raw GeoJSON, sibling QuPath export script). Only writes under
results/mouse_only_v5/spatial_validation/ -- and never overwrites an existing
*_glomerulus_vessel_manual_review.geojson. Run with the `agproject` conda env
with compatible shapely, scanpy, anndata, and scikit-learn versions.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import shapely
from scipy import stats
from scipy.spatial import cKDTree
from shapely.geometry import MultiPoint, shape
from sklearn.metrics import adjusted_rand_score

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
import export_kept_tubules_to_qupath_geojson as qupath_export  # noqa: E402

# ---------------------------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------------------------

DPT_FILENAME = 'all_mouse_tubules_scanpy_dpt.h5ad'
HARMONY_FILENAME = 'all_mouse_tubules_harmony_pass1.h5ad'
OUTPUT_DIRNAME = 'spatial_validation'

GENE_ALIAS_COLUMNS = ['gene_symbol', 'gene_symbols', 'symbol', 'gene_name', 'feature_name', 'features']

# Notebook 6 cell 17's BROAD_SEGMENT_MARKERS, Glomerulus/Vessel entries only (verbatim).
DISAMBIGUATION_MARKERS = {
    'Glomerulus': ['Nphs1', 'Nphs2', 'Podxl'],
    'Vessel': ['Pecam1', 'Cdh5', 'Kdr', 'Emcn'],
}

# Verified byte-for-byte against notebook 6's current (updated) PT_MARKER_GROUPS /
# TAL_MARKER_GROUPS / DCT_MARKER_GROUPS via direct grep of the live notebook source
# (2026-07-16) -- not the older commented-out version earlier in the same cell. Notebook 6 is
# not importable (it's a notebook, not a module), hence hardcoded here. Short codes used as dict
# keys in place of the notebook's longer descriptive TAL group names ('mTAL / early TAL-like'
# etc.) -- gene sets are unchanged.
PT_MARKER_GROUPS = {
    'PT-S1': ['Slc5a2', 'Slc5a12', 'Slc22a12', 'Slc17a3'],
    'PT-S2': ['Slc22a8', 'Slc22a6', 'Slc13a3', 'Slc22a2', 'Slc34a1'],
    'PT-S3': ['Slc7a13', 'Slc5a8', 'Slc38a3', 'Hao2', 'Pck1', 'Acsm3'],
}
TAL_MARKER_GROUPS = {
    'mTAL': ['Slc12a1', 'Umod', 'Kcnj1', 'Bsnd', 'Clcnkb', 'Cldn10'],
    'cTAL': ['Cldn16', 'Cldn19', 'Casr', 'Pth1r', 'Vdr'],
    'MD': ['Nos1', 'Ptgs2'],
}
DCT_MARKER_GROUPS = {
    'DCT1': ['Slc12a3', 'Pvalb', 'Kcnj10'],
    'DCT2': ['Calb1', 'Hsd11b2', 'Trpv5'],
}
SUBSEGMENT_ORDER = {'PT': ['PT-S1', 'PT-S2', 'PT-S3'], 'TAL': ['mTAL', 'cTAL', 'MD'], 'DCT': ['DCT1', 'DCT2']}
SUBSEGMENT_MARKER_GROUPS = {'PT': PT_MARKER_GROUPS, 'TAL': TAL_MARKER_GROUPS, 'DCT': DCT_MARKER_GROUPS}

SEGMENT_ORDER = ['PT', 'TL', 'TAL', 'DCT', 'CNT_CD']

N_BOOTSTRAP = 1000  # fixed, not CLI-parameterized (per spec)
CI_LEVEL = 0.95

MANUAL_REVIEW_FILENAME_TEMPLATE = '{sample}_glomerulus_vessel_manual_review.geojson'
DIFF_REPORT_FILENAME = 'manual_annotation_diff.txt'

# Both the h5ad's own vocabulary ('Glomerulus') and the user's segmentation vocabulary
# ('Glomeruli', adopted here in QuPath) refer to the same class. 'Glomeruli' is canonical for
# all output going forward (see module docstring); this map accepts either, case-insensitively,
# on read.
GLOM_LABEL_ALIASES = {'glomerulus': 'Glomeruli', 'glomeruli': 'Glomeruli', 'vessel': 'Vessel'}
# qupath_export.LABEL_COLORS is keyed by 'Glomerulus' (singular); reuse its color under the new
# canonical 'Glomeruli' key so newly-written seed files keep the same color the user already
# associates with this class from their other QuPath overlays.
GLOM_VESSEL_SEED_COLORS = {
    'Glomeruli': qupath_export.LABEL_COLORS['Glomerulus'],
    'Vessel': qupath_export.LABEL_COLORS['Vessel'],
}
# Pixel tolerance for centroid-based re-identification of manually-edited QuPath features against
# the raw segmentation GeoJSON (QuPath strips all custom properties on export/re-import, so
# feature identity can't be read back directly -- see read_manual_glomerulus_labels). Verified
# empirically to be an exact, distance-0.0 match for real edited files; this margin only guards
# against floating-point rounding in QuPath's own export.
MANUAL_MATCH_TOLERANCE_PX = 1.0

SUBSEGMENT_COLORS = {
    'PT-S1': [150, 200, 255], 'PT-S2': [70, 130, 220], 'PT-S3': [20, 60, 150],   # light->dark blue
    'mTAL': [255, 200, 130], 'cTAL': [230, 130, 40], 'MD': [160, 70, 10],        # light->dark orange
    'DCT1': [160, 220, 150], 'DCT2': [50, 140, 60],                             # light->dark green
}  # fallback for anything else (TL, CNT_CD, Transitional): qupath_export.LABEL_COLORS, then
   # qupath_export.hash_fallback_color for anything genuinely unrecognized.

TEST_SPECS = [  # 3 axes x 2 subsets, full factorial
    dict(axis='dist_glom_nearest', tubule_subset='all'),
    dict(axis='dist_glom_k3', tubule_subset='all'),
    dict(axis='dist_edge', tubule_subset='all'),
    dict(axis='dist_glom_nearest', tubule_subset='PT'),
    dict(axis='dist_glom_k3', tubule_subset='PT'),
    dict(axis='dist_edge', tubule_subset='PT'),
]


# ---------------------------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------------------------

def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root', type=Path, default=qupath_export.DEFAULT_PROJECT_ROOT)
    parser.add_argument('--data-root', type=Path,
                        help='Directory containing private v4 GeoJSON files (default: '
                             '<project-root>/data).')
    parser.add_argument('--results-root', type=Path,
                        help='Directory containing mouse_only_v5 outputs (default: '
                             '<project-root>/results).')
    parser.add_argument('--samples', nargs='+', default=list(qupath_export.SAMPLE_TO_RAW_GEOJSON.keys()))
    parser.add_argument('--n-permutations', type=int, default=10_000)
    parser.add_argument('--k-glom', type=int, default=3)
    parser.add_argument('--dry-run', action='store_true',
                         help='Seed the per-sample Glomerulus/Vessel manual-review GeoJSON '
                              '(skipping samples that already have one) and exit, for manual '
                              'correction in QuPath before the rest of the pipeline runs.')
    parser.add_argument('--random-seed', type=int, default=0,
                         help='Seeds permutation/bootstrap RNG and sc.tl.score_genes for '
                              'reproducibility. Not in the original spec; added deliberately.')
    parser.add_argument('--use-manual-annotations', action='store_true',
                         help='Skip the Step 1 automated-seed reconstruction and diff report '
                              'entirely (no Harmony reload) and go straight to reading the '
                              'existing *_glomerulus_vessel_manual_review.geojson files as ground '
                              'truth. For fast re-runs once manual curation is finalized. Mutually '
                              'exclusive with --dry-run.')
    return parser.parse_args(argv)


# ---------------------------------------------------------------------------------------------
# Gene resolution (ported from notebook 6 cell 17, reused by Step 1 and Step 8)
# ---------------------------------------------------------------------------------------------

def build_gene_lookup(adata: ad.AnnData) -> dict:
    lookup = {}
    for var_name in adata.var_names:
        lookup.setdefault(str(var_name).upper(), str(var_name))
    for col in GENE_ALIAS_COLUMNS:
        if col in adata.var.columns:
            for var_name, alias in zip(adata.var_names, adata.var[col]):
                if pd.isna(alias):
                    continue
                lookup.setdefault(str(alias).upper(), str(var_name))
    return lookup


def resolve_marker_genes(adata: ad.AnnData, marker_dict: dict, label_prefix: str = '',
                          fatal_on_empty: bool = True):
    """Resolves each family's gene panel via build_gene_lookup. If fatal_on_empty, sys.exit(1)
    on any zero-gene family (Step 1's Glomerulus/Vessel disambiguation). If not fatal_on_empty,
    returns None if ANY family in marker_dict resolved to zero genes (Step 8's per-outer-family
    sub-panel resolution, where the caller skips Method B for that whole family instead)."""
    lookup = build_gene_lookup(adata)
    resolved = {}
    any_empty = False
    for family, genes in marker_dict.items():
        found = [lookup[g.upper()] for g in genes if g.upper() in lookup]
        missing = [g for g in genes if g.upper() not in lookup]
        msg = f'  {label_prefix}{family}: resolved {len(found)}/{len(genes)} markers'
        if missing:
            msg += f'; missing: {missing}'
        print(msg)
        if not found:
            any_empty = True
            if fatal_on_empty:
                print(f'FATAL: {label_prefix}{family} resolved to zero marker genes.')
                sys.exit(1)
        resolved[family] = found
    if any_empty and not fatal_on_empty:
        return None
    return resolved


def compute_broad_scores(adata: ad.AnnData, resolved_markers: dict, random_state: int,
                          score_name_template: str, pool_label: str) -> list:
    """Runs sc.tl.score_genes per family in resolved_markers, in place on adata.obs. Reused by
    Step 1 (Glomerulus/Vessel, score_name_template='{family}_broad_score') and Step 8
    (per sub-segment, score_name_template='{family}_score')."""
    if 'lognorm' not in adata.layers:
        print(f"FATAL: adata.layers['lognorm'] missing -- cannot verify .X basis before "
              f"computing {pool_label} scores.")
        sys.exit(1)
    x_slice = adata.X[:5]
    ln_slice = adata.layers['lognorm'][:5]
    if sp.issparse(x_slice):
        x_slice = x_slice.toarray()
    if sp.issparse(ln_slice):
        ln_slice = ln_slice.toarray()
    if not np.allclose(np.asarray(x_slice), np.asarray(ln_slice), equal_nan=True):
        print(f"FATAL: adata.X does not match adata.layers['lognorm'] on the first 5 rows -- "
              f"score_genes would run on an unexpected basis. Aborting before computing "
              f"{pool_label} scores.")
        sys.exit(1)

    print(f"Note: {pool_label} scores computed over full {adata.n_vars}-gene control pool "
          f"(persisted h5ad has no HVG mask on .X), not notebook 6 Section 2's HVG-subsetted "
          f"pool. Score magnitudes will differ slightly from notebook 6 obs values, but relative "
          f"ranking within each labeled cluster -- which drives the downstream filter -- is "
          f"preserved.")

    score_cols = []
    for family, genes in resolved_markers.items():
        score_col = score_name_template.format(family=family)
        sc.tl.score_genes(adata, genes, score_name=score_col, use_raw=False, random_state=random_state)
        score_cols.append(score_col)
    return score_cols


# ---------------------------------------------------------------------------------------------
# Shared GeoJSON writer (Step 1's QC file and Step 8's two sub-segment files)
# ---------------------------------------------------------------------------------------------

def write_labeled_geojson(data_root: Path, output_path: Path, labeled_obs: pd.DataFrame,
                           samples: list, label_column: str, color_map: dict,
                           extra_property_columns: dict | None = None,
                           extra_constant_properties: dict | None = None) -> int:
    """Combined FeatureCollection across samples, joined via qupath_export.parse_feature_index
    (the confirmed f'{sample}_unit_{feature_index}' key). label_column drives
    classification.name/color (falling back to qupath_export.LABEL_COLORS, then
    qupath_export.hash_fallback_color, for labels not in color_map). extra_property_columns:
    {property_name: obs_column_name} carried through per-row. extra_constant_properties: fixed
    key/value pairs added to every feature. Returns total feature count written."""
    all_features = []
    top_level_keys = None
    for sample in samples:
        raw_payload = qupath_export.load_raw_geojson(data_root, sample)
        if top_level_keys is None:
            top_level_keys = {k: v for k, v in raw_payload.items() if k != 'features'}
        raw_features = raw_payload.get('features', [])
        sample_obs = labeled_obs.loc[labeled_obs['sample'].astype(str) == sample]
        for obs_name, row in sample_obs.iterrows():
            feature_index = qupath_export.parse_feature_index(obs_name, sample)
            if not (0 <= feature_index < len(raw_features)):
                print(f'FATAL: {obs_name}: feature_index {feature_index} out of range for '
                      f'sample {sample!r} (n_features={len(raw_features)}).')
                sys.exit(1)
            feature = raw_features[feature_index]
            label = str(row[label_column])
            color = color_map.get(label)
            if color is None:
                color = qupath_export.LABEL_COLORS.get(label)
            if color is None:
                color = qupath_export.hash_fallback_color(label)
            properties = dict(feature.get('properties', {}))
            properties['classification'] = {'name': label, 'color': color}
            properties['name'] = label
            properties['source_sample'] = sample
            # Explicit, so a file can be read back (e.g. after manual editing in QuPath) via
            # (source_sample, feature_index) rather than relying on array position, which QuPath
            # may not preserve on re-export.
            properties['feature_index'] = feature_index
            if extra_property_columns:
                for prop_name, obs_col in extra_property_columns.items():
                    val = row[obs_col]
                    if pd.isna(val):
                        properties[prop_name] = None
                    elif isinstance(val, (int, float, np.floating, np.integer)):
                        properties[prop_name] = float(val)
                    else:
                        properties[prop_name] = str(val)
            if extra_constant_properties:
                properties.update(extra_constant_properties)
            new_feature = dict(feature)
            new_feature['properties'] = properties
            all_features.append(new_feature)

    output_payload = dict(top_level_keys or {})
    output_payload['features'] = all_features
    with output_path.open('w', encoding='utf-8') as f:
        json.dump(output_payload, f, ensure_ascii=False)
    return len(all_features)


# ---------------------------------------------------------------------------------------------
# Step 1 -- Glomerulus/Vessel labeling (manually curated in QuPath; see module docstring).
# Three automated filter designs (score-percentile + gap + celltype_primary + celltype_secondary;
# then dropping celltype_secondary; then dropping the score percentile too) were each tried and
# each still visually excluded real glomeruli as "likely vessels" on inspection in QuPath. The
# seed below no longer makes an automated accept/reject decision -- it only narrows the
# Glomerulus DE-cluster down to cells with an unambiguous original per-cell label (Glomerulus or
# Vessel) for a human to review and correct.
# ---------------------------------------------------------------------------------------------

def normalize_glom_vessel_label(label) -> str | None:
    """Case-insensitive: accepts both 'Glomerulus' (the h5ad's/notebook 6's own vocabulary) and
    'Glomeruli' (the user's segmentation-label vocabulary, adopted in QuPath) as the same class,
    plus 'Vessel'. Returns the canonical form ('Glomeruli' or 'Vessel') or None if label doesn't
    match either (e.g. a feature manually reclassified entirely outside Glomerulus/Vessel scope,
    such as 'Tubules' or 'Blood Cells' -- not an error, just out of scope for this file)."""
    if label is None:
        return None
    return GLOM_LABEL_ALIASES.get(str(label).strip().lower())


def prepare_glomerulus_seed_labels(glom_obs: pd.DataFrame, samples: list) -> pd.DataFrame:
    """Restricts the Glomerulus DE-cluster to cells whose original per-cell classification
    (celltype_primary, from notebook 6) is 'Glomerulus' or 'Vessel' -- dropping cells whose
    per-cell call is neither. glom_qc_class is the normalized canonical label ('Glomeruli'/
    'Vessel'), used for both newly-written seed files and everywhere downstream.
    glom_qc_class_raw preserves the literal celltype_primary string ('Glomerulus'/'Vessel', the
    h5ad's own vocabulary) for later comparison in the manual-annotation diff report. Also
    carries glom_vs_vessel_gap (recomputed score gap) for reference only."""
    glom_obs = glom_obs.copy()

    print('\nPer-sample Glomerulus_broad_score / Vessel_broad_score distributions '
          '(Glomerulus-labeled cluster; reference only, not used to filter):')
    for sample in samples:
        sub = glom_obs.loc[glom_obs['sample'].astype(str) == sample]
        if sub.empty:
            print(f'  {sample}: no Glomerulus-cluster cells')
            continue
        g_q = np.percentile(sub['Glomerulus_broad_score'], [0, 25, 50, 75, 100])
        v_q = np.percentile(sub['Vessel_broad_score'], [0, 25, 50, 75, 100])
        print(f'  {sample}: n={len(sub)}  Glom[min/25/50/75/max]={np.round(g_q, 4).tolist()}  '
              f'Vessel[min/25/50/75/max]={np.round(v_q, 4).tolist()}')

    glom_obs['glom_vs_vessel_gap'] = glom_obs['Glomerulus_broad_score'] - glom_obs['Vessel_broad_score']

    keep_mask = glom_obs['celltype_primary'].astype(str).isin(['Glomerulus', 'Vessel'])
    glom_obs = glom_obs.loc[keep_mask].copy()
    glom_obs['glom_qc_class_raw'] = glom_obs['celltype_primary'].astype(str)
    glom_obs['glom_qc_class'] = glom_obs['glom_qc_class_raw'].map(normalize_glom_vessel_label)

    print('\nPer-sample seed label counts (celltype_primary == Glomerulus or Vessel; other '
          'per-cell calls dropped):')
    any_empty = False
    for sample in samples:
        sub = glom_obs.loc[glom_obs['sample'].astype(str) == sample]
        n_glom = int((sub['glom_qc_class'] == 'Glomeruli').sum())
        n_vessel = int((sub['glom_qc_class'] == 'Vessel').sum())
        print(f'  {sample}: n_labeled_glomeruli={n_glom}  n_labeled_vessel={n_vessel}')
        if n_glom + n_vessel == 0:
            any_empty = True
    if any_empty:
        print('FATAL: at least one sample has zero Glomerulus/Vessel seed cells.')
        sys.exit(1)
    return glom_obs


def match_manual_features_to_raw_index(features: list, raw_centroids: dict, sample: str,
                                        tolerance: float = MANUAL_MATCH_TOLERANCE_PX) -> list:
    """Re-identifies each manually-edited GeoJSON feature against the raw segmentation GeoJSON's
    precomputed centroids (raw_centroids: {feature_index: (x, y)}), by nearest-centroid spatial
    match. Necessary because QuPath strips ALL custom properties (feature_index, source_sample,
    etc.) on export/re-import -- confirmed empirically -- while it preserves geometry exactly
    (also confirmed empirically: distance 0.0, unique, for every feature in real edited files).
    Returns a list of (feature, matched_feature_index) in file order. FATAL if any feature has no
    raw match within `tolerance` pixels (e.g. a brand-new hand-drawn annotation -- it has no
    corresponding entry anywhere else in this pipeline's gene-expression data and cannot be used)
    or if two features match the same raw index (ambiguous/duplicate annotation)."""
    raw_indices = list(raw_centroids.keys())
    raw_xy = np.array([raw_centroids[i] for i in raw_indices])
    tree = cKDTree(raw_xy)
    results = []
    matched_position_by_index = {}
    for position, feature in enumerate(features):
        centroid = shape(feature['geometry']).centroid
        dist, nn = tree.query([centroid.x, centroid.y], k=1)
        matched_feature_index = raw_indices[nn]
        if dist > tolerance:
            print(f'FATAL: {sample}: manually-edited feature at position {position} has no '
                  f'matching raw polygon within {tolerance}px (nearest is {dist:.2f}px away). '
                  f'This looks like a newly hand-drawn or reshaped annotation, which has no '
                  f'corresponding tubule anywhere else in this pipeline (no gene expression, no '
                  f'pseudospace value) and cannot be used. Delete it in QuPath if unintended, or '
                  f'investigate further if you meant to add a new segmentation region.')
            sys.exit(1)
        if matched_feature_index in matched_position_by_index:
            print(f'FATAL: {sample}: manually-edited features at positions '
                  f'{matched_position_by_index[matched_feature_index]} and {position} both '
                  f'matched raw feature_index {matched_feature_index} -- ambiguous/duplicate '
                  f'annotation, cannot disambiguate.')
            sys.exit(1)
        matched_position_by_index[matched_feature_index] = position
        results.append((feature, matched_feature_index))
    return results


def read_manual_glomerulus_labels(output_dir: Path, samples: list, centroids_by_sample: dict):
    """Reads back the manually-curated per-sample Glomerulus/Vessel GeoJSON files. Features are
    re-identified via match_manual_features_to_raw_index (nearest-centroid spatial match, since
    QuPath strips custom properties on export). Returns (classified_df, dropped_df):
      classified_df -- obs_name-indexed, columns sample/glom_qc_class (normalized canonical
        'Glomeruli'/'Vessel')/glom_qc_class_raw (literal on-disk classification.name) -- the
        ground truth used by the rest of the pipeline.
      dropped_df -- obs_name-indexed, columns sample/glom_qc_class_raw -- features the user
        manually reclassified to something outside Glomerulus/Vessel scope entirely (e.g.
        'Tubules', 'Blood Cells'); not an error, just excluded from the high-confidence set."""
    kept_rows = []
    dropped_rows = []
    for sample in samples:
        path = output_dir / MANUAL_REVIEW_FILENAME_TEMPLATE.format(sample=sample)
        if not path.exists():
            print(f'FATAL: manual review file not found for sample {sample!r}: {path}\n'
                  f'Run with --dry-run first to generate it, correct it in QuPath, save it back '
                  f'to this exact path, then re-run without --dry-run.')
            sys.exit(1)
        try:
            payload = json.loads(path.read_text())
        except json.JSONDecodeError as exc:
            print(f'FATAL: {sample}: {path} is not valid JSON ({exc}).')
            sys.exit(1)
        features = payload.get('features', [])
        if not features:
            print(f'FATAL: {sample}: {path} has zero features -- malformed or emptied by mistake.')
            sys.exit(1)
        if sample not in centroids_by_sample:
            print(f'FATAL: no raw centroids loaded for sample {sample!r}; cannot re-identify '
                  f'manually-edited features.')
            sys.exit(1)

        matches = match_manual_features_to_raw_index(features, centroids_by_sample[sample], sample)
        for feature, feature_index in matches:
            props = feature.get('properties', {}) or {}
            classification = props.get('classification')
            label = classification.get('name') if isinstance(classification, dict) else props.get('name')
            if label is None:
                print(f'FATAL: {sample}: feature matched to feature_index {feature_index} has no '
                      f'classification.name -- malformed.')
                sys.exit(1)
            obs_name = f'{sample}_unit_{feature_index}'
            normalized = normalize_glom_vessel_label(label)
            if normalized is None:
                dropped_rows.append(dict(obs_name=obs_name, sample=sample, glom_qc_class_raw=str(label)))
            else:
                kept_rows.append(dict(obs_name=obs_name, sample=sample, glom_qc_class=normalized,
                                       glom_qc_class_raw=str(label)))

        n_glom = sum(1 for r in kept_rows if r['sample'] == sample and r['glom_qc_class'] == 'Glomeruli')
        n_vessel = sum(1 for r in kept_rows if r['sample'] == sample and r['glom_qc_class'] == 'Vessel')
        n_dropped = sum(1 for r in dropped_rows if r['sample'] == sample)
        print(f'  {sample}: read back n_glomeruli={n_glom}  n_vessel={n_vessel}  '
              f'n_dropped_other_class={n_dropped}  from {path}')

    classified_df = (pd.DataFrame(kept_rows).set_index('obs_name') if kept_rows
                      else pd.DataFrame(columns=['sample', 'glom_qc_class', 'glom_qc_class_raw']))
    dropped_df = (pd.DataFrame(dropped_rows).set_index('obs_name') if dropped_rows
                  else pd.DataFrame(columns=['sample', 'glom_qc_class_raw']))
    if classified_df.empty or classified_df['sample'].nunique() < len(samples):
        print('FATAL: at least one sample produced zero Glomerulus/Vessel rows after reading '
              'back manual labels.')
        sys.exit(1)
    return classified_df, dropped_df


def write_manual_annotation_diff_report(seed_df: pd.DataFrame, classified_df: pd.DataFrame,
                                         dropped_df: pd.DataFrame, samples: list,
                                         manual_paths: dict, output_path: Path) -> None:
    """Compares the in-memory 'auto' reconstruction (seed_df, from prepare_glomerulus_seed_labels
    -- rebuilt fresh for this comparison only, NEVER written to disk) against the actual on-disk
    manually-curated files (classified_df + dropped_df, from read_manual_glomerulus_labels).
    Writes a plain-text report to output_path and prints the identical text to console."""
    lines = []
    emit = lines.append

    auto_small = pd.DataFrame({'sample': seed_df['sample'], 'label_norm': seed_df['glom_qc_class']})
    manual_small = pd.DataFrame({'sample': classified_df['sample'], 'label_norm': classified_df['glom_qc_class']})

    emit('=' * 78)
    emit('MANUAL ANNOTATION DIFF REPORT')
    emit(f'Generated: {datetime.now().isoformat(timespec="seconds")}')
    emit("Auto = in-memory-only reconstruction of Step 1's current seed logic "
         '(celltype_primary-derived), rebuilt fresh for this comparison; NOT written to disk, so '
         'your manual edits are never overwritten.')
    emit('Manual = the on-disk, user-curated files:')
    for sample in samples:
        emit(f'  {sample}: {manual_paths[sample]}')
    emit(f'Total features: auto={len(auto_small)}  manual (Glomeruli/Vessel)={len(manual_small)}  '
         f'manual (reclassified outside scope)={len(dropped_df)}')
    emit('=' * 78)

    # --- Label vocabulary changes ---
    emit('\n--- LABEL VOCABULARY CHANGES ---')
    auto_vocab = seed_df['glom_qc_class_raw'].astype(str).value_counts()
    manual_raw_all = pd.concat([classified_df['glom_qc_class_raw'], dropped_df['glom_qc_class_raw']])
    manual_vocab = manual_raw_all.astype(str).value_counts()
    emit('Auto (raw celltype_primary values -- the h5ad/notebook 6 vocabulary):')
    for label, count in auto_vocab.items():
        emit(f'  {label}: {count}')
    emit('Manual (raw on-disk classification.name values, all classes including anything outside '
         'Glomerulus/Vessel scope):')
    for label, count in manual_vocab.items():
        emit(f'  {label}: {count}')
    if 'Glomerulus' in auto_vocab.index and 'Glomeruli' in manual_vocab.index:
        emit("\nNOTE: the manual file renamed the 'Glomerulus' class to 'Glomeruli' (matching the "
             "user's original segmentation label vocabulary). This is a label-text change only -- "
             "both normalize to the same canonical class ('Glomeruli', now used throughout the "
             "pipeline and for all future auto-generated seed files) and are NOT counted as "
             "reclassifications in the section below.")
    other_labels = [lbl for lbl in manual_vocab.index if normalize_glom_vessel_label(lbl) is None]
    if other_labels:
        emit(f"\nNOTE: {int(manual_vocab.loc[other_labels].sum())} feature(s) were manually "
             f"reclassified entirely outside Glomerulus/Vessel scope ({', '.join(other_labels)}) "
             f"-- these are excluded from the high-confidence set used downstream, not an error.")

    # --- Per-feature classification changes (normalized comparison only) ---
    joined = auto_small.join(manual_small, how='inner', lsuffix='_auto', rsuffix='_manual')
    changed = joined.loc[joined['label_norm_auto'] != joined['label_norm_manual']]
    emit('\n--- PER-FEATURE CLASSIFICATION CHANGES (normalized: Glomeruli vs Vessel only; the '
         'Glomerulus/Glomeruli vocabulary rename above is not counted here) ---')
    if changed.empty:
        emit('  (none -- no feature changed between Glomeruli and Vessel)')
    else:
        for (auto_lbl, manual_lbl), group in changed.groupby(['label_norm_auto', 'label_norm_manual']):
            emit(f'\nauto={auto_lbl} -> manual={manual_lbl}: {len(group)} feature(s)')
            for obs_name, row in group.iterrows():
                sample = row['sample_auto']
                feature_index = qupath_export.parse_feature_index(obs_name, sample)
                emit(f'  sample={sample}, feature_index={feature_index}, auto={auto_lbl}, manual={manual_lbl}')

    # --- Features added/removed ---
    emit('\n--- FEATURES ADDED BY MANUAL EDIT (present in manual, not matched to any auto seed '
         'feature) ---')
    emit('  none possible by construction: an unmatched manual annotation fails the pipeline '
         '(see match_manual_features_to_raw_index) before this report is reached.')

    only_in_auto = auto_small.index.difference(manual_small.index.union(dropped_df.index))
    emit('\n--- FEATURES REMOVED BY MANUAL EDIT ---')
    emit('(a) present in auto reconstruction, absent entirely from the manual file (deleted in QuPath):')
    if len(only_in_auto) == 0:
        emit('  none')
    else:
        for obs_name in only_in_auto:
            sample = auto_small.loc[obs_name, 'sample']
            feature_index = qupath_export.parse_feature_index(obs_name, sample)
            emit(f'  sample={sample}, feature_index={feature_index}, auto={auto_small.loc[obs_name, "label_norm"]}')
    emit('(b) present in auto reconstruction, reclassified in the manual file to something '
         'outside Glomerulus/Vessel scope:')
    reclassified_out = dropped_df.loc[dropped_df.index.intersection(auto_small.index)]
    if reclassified_out.empty:
        emit('  none')
    else:
        for obs_name, row in reclassified_out.iterrows():
            sample = row['sample']
            feature_index = qupath_export.parse_feature_index(obs_name, sample)
            auto_lbl = auto_small.loc[obs_name, 'label_norm']
            emit(f'  sample={sample}, feature_index={feature_index}, auto={auto_lbl}, '
                 f'manual_reclassified_to={row["glom_qc_class_raw"]}')

    # --- Per-sample summary table ---
    emit('\n--- PER-SAMPLE SUMMARY ---')
    emit(f'{"sample":<10} {"n_features_auto":>16} {"n_features_manual":>18} {"n_class_changes":>16} '
         f'{"n_added":>8} {"n_removed":>10}   manual class counts')
    per_sample_changed = {}
    for sample in samples:
        n_features_auto = int((auto_small['sample'] == sample).sum())
        n_features_manual = int((manual_small['sample'] == sample).sum())
        n_class_changes = int((changed['sample_auto'] == sample).sum()) if not changed.empty else 0
        per_sample_changed[sample] = n_class_changes
        n_removed = int((only_in_auto.to_series().map(lambda n: auto_small.loc[n, 'sample']) == sample).sum()) \
            if len(only_in_auto) else 0
        n_removed += int((reclassified_out['sample'] == sample).sum()) if not reclassified_out.empty else 0
        n_added = 0  # always 0 by construction -- see "FEATURES ADDED" section above
        manual_counts = classified_df.loc[classified_df['sample'] == sample, 'glom_qc_class'].value_counts()
        counts_str = ', '.join(f'{k}={v}' for k, v in manual_counts.items())
        emit(f'{sample:<10} {n_features_auto:>16} {n_features_manual:>18} {n_class_changes:>16} '
             f'{n_added:>8} {n_removed:>10}   {counts_str}')

    # --- Interpretation notes ---
    emit('\n--- INTERPRETATION NOTES ---')
    total_auto = len(auto_small)
    total_changed = len(changed)
    pct_changed = 100.0 * total_changed / total_auto if total_auto else 0.0
    to_glomeruli = int((changed['label_norm_manual'] == 'Glomeruli').sum()) if not changed.empty else 0
    to_vessel = int((changed['label_norm_manual'] == 'Vessel').sum()) if not changed.empty else 0
    if total_changed == 0:
        symmetry_note = 'No reclassifications occurred.'
    elif abs(to_glomeruli - to_vessel) < 0.2 * max(to_glomeruli, to_vessel, 1):
        symmetry_note = 'Disagreements were roughly symmetric between the two directions.'
    elif to_glomeruli > to_vessel:
        symmetry_note = ('Disagreements were predominantly one-directional: the automated seed '
                          'under-called real glomeruli (mislabeled them Vessel).')
    else:
        symmetry_note = ('Disagreements were predominantly one-directional: the automated seed '
                          'over-called vessels as glomeruli.')
    max_sample = max(per_sample_changed, key=per_sample_changed.get) if per_sample_changed else None
    max_sample_note = ''
    if max_sample is not None and per_sample_changed[max_sample] > 0:
        max_sample_note = (f' {max_sample} had the most manual corrections '
                            f'({per_sample_changed[max_sample]} features)')
        if max_sample == 'IR2A2':
            max_sample_note += (', consistent with IR2A2 being called out elsewhere in this '
                                 'pipeline as a potential outlier -- worth checking whether this '
                                 'reflects batch/segmentation-quality issues specific to it.')
        else:
            max_sample_note += ' -- worth a quick look at whether this reflects a batch/segmentation-quality issue.'
    emit(f'Overall, {total_changed}/{total_auto} ({pct_changed:.1f}%) of the automated seed labels '
         f'were manually corrected between Glomeruli and Vessel ({to_glomeruli} corrected from '
         f'auto=Vessel to manual=Glomeruli, {to_vessel} corrected from auto=Glomeruli to '
         f'manual=Vessel). {symmetry_note}{max_sample_note}')

    report_text = '\n'.join(lines)
    output_path.write_text(report_text + '\n')
    print(report_text)


# ---------------------------------------------------------------------------------------------
# Step 2 -- Centroids + tissue edge
# ---------------------------------------------------------------------------------------------

def load_all_centroids(data_root: Path, samples: list) -> dict:
    centroids = {}
    for sample in samples:
        payload = qupath_export.load_raw_geojson(data_root, sample)
        features = payload.get('features', [])
        sample_centroids = {}
        for i, feature in enumerate(features):
            c = shape(feature['geometry']).centroid
            sample_centroids[i] = (c.x, c.y)
        centroids[sample] = sample_centroids
        print(f'  {sample}: {len(sample_centroids)} centroids computed from raw GeoJSON')
    return centroids


def compute_tissue_edge_polygon(centroid_xy: np.ndarray, ratio: float = 0.05):
    mp = MultiPoint(centroid_xy)
    if hasattr(shapely, 'concave_hull'):
        poly = shapely.concave_hull(mp, ratio=ratio, allow_holes=False)
        print(f'  tissue edge: shapely.concave_hull(ratio={ratio}, allow_holes=False), '
              f'{len(centroid_xy)} points -> polygon area={poly.area:.1f}')
        return poly
    try:
        import alphashape
        poly = alphashape.alphashape(centroid_xy, alpha=ratio)
        print(f'  tissue edge: alphashape.alphashape(alpha={ratio}) [fallback path, untested '
              f'in the agproject environment since shapely.concave_hull is available there], '
              f'{len(centroid_xy)} points -> polygon area={poly.area:.1f}')
        return poly
    except ImportError:
        print('FATAL: neither shapely.concave_hull nor alphashape is available in this environment.')
        sys.exit(1)


def attach_centroids(obs: pd.DataFrame, samples: list, centroids_by_sample: dict) -> pd.DataFrame:
    """Adds centroid_x/centroid_y computed fresh from raw GeoJSON via the join key -- does not
    read existing x_centroid/y_centroid columns or obsm['spatial']."""
    obs = obs.copy()
    cx = np.full(len(obs), np.nan)
    cy = np.full(len(obs), np.nan)
    obs_names = obs.index.to_numpy()
    sample_values = obs['sample'].astype(str).to_numpy()
    missing = []
    for i, (obs_name, sample) in enumerate(zip(obs_names, sample_values)):
        if sample not in centroids_by_sample:
            missing.append(f'{obs_name}: unknown sample {sample!r}')
            continue
        try:
            feature_index = qupath_export.parse_feature_index(obs_name, sample)
        except ValueError as exc:
            missing.append(f'{obs_name}: {exc}')
            continue
        centroid = centroids_by_sample[sample].get(feature_index)
        if centroid is None:
            missing.append(f'{obs_name}: feature_index {feature_index} not in centroid dict for {sample}')
            continue
        cx[i], cy[i] = centroid
    if missing:
        print(f'FATAL: {len(missing)} centroid join mismatch(es). First 20:')
        for m in missing[:20]:
            print('  -', m)
        sys.exit(1)
    obs['centroid_x'] = cx
    obs['centroid_y'] = cy
    return obs


# ---------------------------------------------------------------------------------------------
# Step 3 -- Distances
# ---------------------------------------------------------------------------------------------

def compute_glom_distances(centroid_xy: np.ndarray, hc_glom_xy: np.ndarray, k: int, sample: str):
    n_hc = len(hc_glom_xy)
    if n_hc == 0:
        print(f'FATAL: sample {sample!r} has zero high-confidence glomeruli -- cannot compute distances.')
        sys.exit(1)
    k_eff = min(k, n_hc)
    if k_eff < k:
        print(f'  WARNING: sample {sample!r} has only {n_hc} high-confidence glomeruli (< k={k}); '
              f'falling back to k={k_eff} for dist_glom_k3.')
    tree = cKDTree(hc_glom_xy)
    dists, _ = tree.query(centroid_xy, k=k_eff)
    if k_eff == 1:
        dist_nearest = np.asarray(dists).reshape(-1)
        dist_k = dist_nearest.copy()
    else:
        dist_nearest = dists[:, 0]
        dist_k = dists.mean(axis=1)
    return dist_nearest, dist_k


def compute_edge_distances(centroid_xy: np.ndarray, tissue_polygon) -> np.ndarray:
    """Unsigned distance to the tissue-boundary polygon (no inside/outside sign logic --
    removed as gratuitous for a correlation test and a source of edge-case ambiguity)."""
    points = shapely.points(centroid_xy[:, 0], centroid_xy[:, 1])
    boundary = tissue_polygon.boundary
    return shapely.distance(points, boundary)


# ---------------------------------------------------------------------------------------------
# Step 4 -- Per-sample Spearman + permutation + bootstrap (vectorized)
# ---------------------------------------------------------------------------------------------

def _rank_zscore(a: np.ndarray) -> np.ndarray:
    r = stats.rankdata(a)
    return (r - r.mean()) / r.std(ddof=0)


def spearman_permutation_test(x: np.ndarray, y: np.ndarray, n_permutations: int,
                               rng: np.random.Generator, chunk_size: int = 2000):
    """Ranks computed once (shuffling never changes rank values); each permutation reduces to a
    single BLAS-backed matrix-vector product instead of a full spearmanr recomputation."""
    n = len(x)
    zx = _rank_zscore(x)
    zy = _rank_zscore(y)
    rho_obs = float(zx @ zy / n)
    count_ge = 0
    remaining = n_permutations
    while remaining > 0:
        chunk = min(chunk_size, remaining)
        idx = np.argsort(rng.random((chunk, n)), axis=1)
        zy_perm = zy[idx]
        rho_perm = zy_perm @ zx / n
        count_ge += int(np.sum(np.abs(rho_perm) >= abs(rho_obs)))
        remaining -= chunk
    pvalue = (count_ge + 1) / (n_permutations + 1)
    return rho_obs, pvalue


def spearman_bootstrap_ci(x: np.ndarray, y: np.ndarray, n_bootstrap: int,
                           rng: np.random.Generator, chunk_size: int = 200,
                           ci_level: float = CI_LEVEL):
    """Resampling with replacement changes tie structure, so ranks are recomputed per batch
    (vectorized ordinal ranking via double-argsort -- immaterial deviation from scipy's
    average-tie ranking here since all inputs are continuous floats with no real ties)."""
    n = len(x)
    rho_boot = np.empty(n_bootstrap)
    pos = 0
    remaining = n_bootstrap
    while remaining > 0:
        chunk = min(chunk_size, remaining)
        idx = rng.integers(0, n, size=(chunk, n))
        xb, yb = x[idx], y[idx]
        rxb = np.argsort(np.argsort(xb, axis=1), axis=1).astype(float)
        ryb = np.argsort(np.argsort(yb, axis=1), axis=1).astype(float)
        zxb = (rxb - rxb.mean(axis=1, keepdims=True)) / rxb.std(axis=1, ddof=0, keepdims=True)
        zyb = (ryb - ryb.mean(axis=1, keepdims=True)) / ryb.std(axis=1, ddof=0, keepdims=True)
        rho_boot[pos:pos + chunk] = (zxb * zyb).mean(axis=1)
        pos += chunk
        remaining -= chunk
    lo_pct, hi_pct = 100 * (1 - ci_level) / 2, 100 * (1 + ci_level) / 2
    ci_lo, ci_hi = np.percentile(rho_boot, [lo_pct, hi_pct])
    return float(ci_lo), float(ci_hi)


def run_correlation_test(x: np.ndarray, y: np.ndarray, n_permutations: int, n_bootstrap: int,
                          rng: np.random.Generator) -> dict:
    if np.any(~np.isfinite(x)) or np.any(~np.isfinite(y)):
        print('FATAL: NaN/inf values passed to run_correlation_test -- upstream join/distance '
              'computation must be complete for all rows.')
        sys.exit(1)
    spearman_rho = float(stats.spearmanr(x, y).statistic)
    _, perm_pvalue = spearman_permutation_test(x, y, n_permutations, rng)
    ci_lo, ci_hi = spearman_bootstrap_ci(x, y, n_bootstrap, rng)
    return dict(spearman_rho=spearman_rho, n_tubules=len(x), perm_pvalue=perm_pvalue,
                ci95_lo=ci_lo, ci95_hi=ci_hi)


def run_all_sample_tests(dpt_obs: pd.DataFrame, samples: list, n_permutations: int,
                          rng: np.random.Generator) -> pd.DataFrame:
    has_pt_subset_col = 'pt_subset_scanpy_dpt' in dpt_obs.columns
    rows = []
    for sample in samples:
        sample_obs = dpt_obs.loc[dpt_obs['sample'].astype(str) == sample]
        for spec in TEST_SPECS:
            axis, tubule_subset = spec['axis'], spec['tubule_subset']
            if tubule_subset == 'PT':
                subset_obs = sample_obs.loc[sample_obs['broad_tubule_marker_call'].astype(str) == 'PT']
                if has_pt_subset_col and subset_obs['pt_subset_scanpy_dpt'].notna().any():
                    dpt_col, dpt_source = 'pt_subset_scanpy_dpt', 'pt_subset_scanpy_dpt'
                else:
                    dpt_col, dpt_source = 'total_scanpy_dpt', 'total_scanpy_dpt_restricted'
            else:
                subset_obs = sample_obs
                dpt_col, dpt_source = 'total_scanpy_dpt', 'total_scanpy_dpt'

            x = subset_obs[axis].to_numpy(dtype=float)
            y = subset_obs[dpt_col].to_numpy(dtype=float)
            result = run_correlation_test(x, y, n_permutations, N_BOOTSTRAP, rng)
            rows.append(dict(sample=sample, axis=axis, tubule_subset=tubule_subset,
                              dpt_source=dpt_source, **result))
    results_df = pd.DataFrame(rows)

    warn_rows = results_df[(results_df['tubule_subset'] == 'PT') &
                            (results_df['axis'] == 'dist_glom_nearest') &
                            (results_df['spearman_rho'] <= 0)]
    if not warn_rows.empty:
        print('\n' + '!' * 78)
        print('!!! WARNING: PT vs dist_glom_nearest Spearman rho is <= 0 for the following '
              'sample(s) -- this is the key validation claim and should be positive: '
              + ', '.join(warn_rows['sample'].tolist()))
        print('!' * 78)
    else:
        print('\nPT vs dist_glom_nearest Spearman rho is positive for all samples (key '
              'validation claim holds).')
    return results_df


# ---------------------------------------------------------------------------------------------
# Step 5 -- Pooled Fisher-z
# ---------------------------------------------------------------------------------------------

def pooled_fisher_z(rhos: np.ndarray, ns: np.ndarray, samples: list) -> dict:
    z = np.arctanh(np.clip(rhos, -0.999999, 0.999999))
    z_bar = z.mean()  # unweighted mean across samples, per original formula
    n_eff = float(np.sum(ns - 3))
    se = 1.0 / np.sqrt(n_eff)
    pooled_rho = float(np.tanh(z_bar))
    ci_lo, ci_hi = float(np.tanh(z_bar - 1.96 * se)), float(np.tanh(z_bar + 1.96 * se))
    min_i, max_i = int(np.argmin(rhos)), int(np.argmax(rhos))
    return dict(pooled_rho=pooled_rho, ci95_lo=ci_lo, ci95_hi=ci_hi, n_effective=n_eff,
                min_rho=float(rhos[min_i]), min_sample=samples[min_i],
                max_rho=float(rhos[max_i]), max_sample=samples[max_i])


def compute_pooled_correlations(results_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (axis, tubule_subset), group in results_df.groupby(['axis', 'tubule_subset']):
        group = group.sort_values('sample')
        pooled = pooled_fisher_z(group['spearman_rho'].to_numpy(), group['n_tubules'].to_numpy(),
                                  group['sample'].tolist())
        rows.append(dict(axis=axis, tubule_subset=tubule_subset, **pooled))
    pooled_df = pd.DataFrame(rows)

    key_row = pooled_df[(pooled_df['axis'] == 'dist_glom_nearest') & (pooled_df['tubule_subset'] == 'PT')]
    if not key_row.empty:
        ci_lo, ci_hi = key_row.iloc[0]['ci95_lo'], key_row.iloc[0]['ci95_hi']
        if ci_lo > 0:
            print(f'\nPASS: pooled within-PT vs dist_glom_nearest Spearman rho 95% CI '
                  f'[{ci_lo:.3f}, {ci_hi:.3f}] excludes 0.')
        else:
            print('\n' + '!' * 78)
            print(f'!!! WARNING: pooled within-PT vs dist_glom_nearest Spearman rho 95% CI '
                  f'[{ci_lo:.3f}, {ci_hi:.3f}] does NOT exclude 0.')
            print('!' * 78)
    return pooled_df


# ---------------------------------------------------------------------------------------------
# Step 6 -- Segment-order test
# ---------------------------------------------------------------------------------------------

def segment_order_test(dpt_obs: pd.DataFrame, samples: list) -> pd.DataFrame:
    rows = []
    for sample in samples:
        sample_obs = dpt_obs.loc[dpt_obs['sample'].astype(str) == sample]
        groups = {seg: sample_obs.loc[sample_obs['broad_tubule_marker_call'].astype(str) == seg,
                                       'dist_glom_nearest'].to_numpy() for seg in SEGMENT_ORDER}
        row = {'sample': sample}
        for seg in SEGMENT_ORDER:
            row[f'{seg}_median_dist_glom_nearest'] = float(np.median(groups[seg])) if len(groups[seg]) else np.nan

        nonempty_groups = [groups[seg] for seg in SEGMENT_ORDER if len(groups[seg]) > 0]
        if len(nonempty_groups) >= 2:
            kw_stat, kw_p = stats.kruskal(*nonempty_groups)
        else:
            print(f'  WARNING: {sample}: fewer than 2 non-empty segments, skipping Kruskal-Wallis test.')
            kw_stat, kw_p = np.nan, np.nan
        row['kruskal_stat'] = float(kw_stat) if np.isfinite(kw_stat) else np.nan
        row['kruskal_pvalue'] = float(kw_p) if np.isfinite(kw_p) else np.nan

        pt_vals = groups['PT']
        for seg in SEGMENT_ORDER:
            if seg == 'PT':
                continue
            seg_vals = groups[seg]
            if len(pt_vals) > 0 and len(seg_vals) > 0:
                mwu_stat, mwu_p = stats.mannwhitneyu(pt_vals, seg_vals, alternative='two-sided')
            else:
                mwu_stat, mwu_p = np.nan, np.nan
            row[f'mwu_PT_vs_{seg}_stat'] = float(mwu_stat) if np.isfinite(mwu_stat) else np.nan
            row[f'mwu_PT_vs_{seg}_pvalue'] = float(mwu_p) if np.isfinite(mwu_p) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------------------------
# Step 7 -- Figures
# ---------------------------------------------------------------------------------------------

def _sample_grid(samples: list, figsize=(11, 9), ncols: int = 2):
    nrows = int(np.ceil(len(samples) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    return fig, np.asarray(axes).flatten()


def _segment_color(seg: str):
    rgb = qupath_export.LABEL_COLORS.get(seg, qupath_export.hash_fallback_color(seg))
    return tuple(c / 255 for c in rgb)


def fig_pseudospace_vs_glom_scatter(dpt_obs: pd.DataFrame, samples: list, results_df: pd.DataFrame,
                                     out_path: Path) -> None:
    fig, axes = _sample_grid(samples)
    for ax, sample in zip(axes, samples):
        sub = dpt_obs.loc[dpt_obs['sample'].astype(str) == sample]
        pt_mask = sub['broad_tubule_marker_call'].astype(str) == 'PT'
        ax.scatter(sub.loc[~pt_mask, 'dist_glom_nearest'], sub.loc[~pt_mask, 'total_scanpy_dpt'],
                   s=4, alpha=0.35, color='#999999', label='other tubules', linewidths=0)
        ax.scatter(sub.loc[pt_mask, 'dist_glom_nearest'], sub.loc[pt_mask, 'total_scanpy_dpt'],
                   s=4, alpha=0.6, color='#3366CC', label='PT', linewidths=0)
        row = results_df[(results_df['sample'] == sample) & (results_df['axis'] == 'dist_glom_nearest')
                          & (results_df['tubule_subset'] == 'PT')]
        rho = row.iloc[0]['spearman_rho'] if not row.empty else np.nan
        pval = row.iloc[0]['perm_pvalue'] if not row.empty else np.nan
        ax.set_title(f'{sample}  (PT rho={rho:.2f}, p={pval:.1e})')
        ax.set_xlabel('distance to nearest high-confidence glomerulus')
        ax.set_ylabel('total_scanpy_dpt')
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right')
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close(fig)


def fig_pseudospace_vs_edge_scatter(dpt_obs: pd.DataFrame, samples: list, results_df: pd.DataFrame,
                                     out_path: Path) -> None:
    fig, axes = _sample_grid(samples)
    for ax, sample in zip(axes, samples):
        sub = dpt_obs.loc[dpt_obs['sample'].astype(str) == sample]
        pt_mask = sub['broad_tubule_marker_call'].astype(str) == 'PT'
        ax.scatter(sub.loc[~pt_mask, 'dist_edge'], sub.loc[~pt_mask, 'total_scanpy_dpt'],
                   s=4, alpha=0.35, color='#999999', label='other tubules', linewidths=0)
        ax.scatter(sub.loc[pt_mask, 'dist_edge'], sub.loc[pt_mask, 'total_scanpy_dpt'],
                   s=4, alpha=0.6, color='#3366CC', label='PT', linewidths=0)
        row = results_df[(results_df['sample'] == sample) & (results_df['axis'] == 'dist_edge')
                          & (results_df['tubule_subset'] == 'PT')]
        rho = row.iloc[0]['spearman_rho'] if not row.empty else np.nan
        pval = row.iloc[0]['perm_pvalue'] if not row.empty else np.nan
        ax.set_title(f'{sample}  (PT rho={rho:.2f}, p={pval:.1e})')
        ax.set_xlabel('distance from tissue edge')
        ax.set_ylabel('total_scanpy_dpt')
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right')
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close(fig)


def fig_segment_order_boxplot(dpt_obs: pd.DataFrame, samples: list, segment_df: pd.DataFrame,
                               out_path: Path) -> None:
    fig, axes = _sample_grid(samples)
    colors = [_segment_color(seg) for seg in SEGMENT_ORDER]
    for ax, sample in zip(axes, samples):
        sub = dpt_obs.loc[dpt_obs['sample'].astype(str) == sample]
        data = [sub.loc[sub['broad_tubule_marker_call'].astype(str) == seg, 'dist_glom_nearest'].to_numpy()
                for seg in SEGMENT_ORDER]
        bp = ax.boxplot(data, labels=SEGMENT_ORDER, patch_artist=True, showfliers=False)
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
        kw_row = segment_df.loc[segment_df['sample'] == sample]
        kw_p = kw_row.iloc[0]['kruskal_pvalue'] if not kw_row.empty else np.nan
        ax.set_title(f'{sample}  (Kruskal-Wallis p={kw_p:.1e})')
        ax.set_ylabel('distance to nearest high-confidence glomerulus')
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close(fig)


def fig_forest_plot(results_df: pd.DataFrame, pooled_df: pd.DataFrame, out_path: Path) -> None:
    order = [f"{s['axis']} | {s['tubule_subset']}" for s in TEST_SPECS]
    pooled_df = pooled_df.copy()
    pooled_df['combo'] = pd.Categorical(pooled_df['axis'] + ' | ' + pooled_df['tubule_subset'],
                                         categories=order, ordered=True)
    pooled_df = pooled_df.sort_values('combo')

    fig, ax = plt.subplots(figsize=(7, 4))
    y_positions = np.arange(len(pooled_df))
    xerr = np.vstack([pooled_df['pooled_rho'] - pooled_df['ci95_lo'],
                       pooled_df['ci95_hi'] - pooled_df['pooled_rho']])
    ax.errorbar(pooled_df['pooled_rho'], y_positions, xerr=xerr, fmt='o', color='#3366CC',
                ecolor='#3366CC', capsize=3, zorder=3)
    for y, (_, row) in zip(y_positions, pooled_df.iterrows()):
        per_sample = results_df.loc[(results_df['axis'] == row['axis']) &
                                     (results_df['tubule_subset'] == row['tubule_subset']),
                                     'spearman_rho']
        ax.scatter(per_sample, np.full(len(per_sample), y), s=15, alpha=0.5, color='#999999', zorder=1)
    ax.axvline(0, color='gray', linestyle='--', linewidth=1, zorder=0)
    ax.set_yticks(y_positions)
    ax.set_yticklabels([str(c) for c in pooled_df['combo']])
    ax.set_xlabel('pooled Spearman rho (95% CI, Fisher-z)')
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close(fig)


def fig_qc_map(dpt_obs: pd.DataFrame, classified: pd.DataFrame, samples: list, out_path: Path) -> None:
    fig, axes = plt.subplots(len(samples), 2, figsize=(11, 4.2 * len(samples)))
    axes = np.atleast_2d(axes)
    vmin, vmax = dpt_obs['total_scanpy_dpt'].min(), dpt_obs['total_scanpy_dpt'].max()
    sc_handle = None
    for row_i, sample in enumerate(samples):
        sub = dpt_obs.loc[dpt_obs['sample'].astype(str) == sample]
        hc = classified.loc[(classified['sample'].astype(str) == sample) &
                             (classified['glom_qc_class'] == 'Glomeruli')]

        ax_left = axes[row_i, 0]
        sc_handle = ax_left.scatter(sub['centroid_x'], sub['centroid_y'], c=sub['total_scanpy_dpt'],
                                     cmap='viridis', vmin=vmin, vmax=vmax, s=4, linewidths=0)
        ax_left.scatter(hc['centroid_x'], hc['centroid_y'], marker='*', s=60, facecolor='white',
                         edgecolor='black', linewidths=0.8, zorder=5, label='high-confidence glomerulus')
        ax_left.set_title(f'{sample}: pseudospace')
        ax_left.set_aspect('equal')
        # display-only convention (image-pixel y-down); does NOT affect the underlying
        # dist_edge/dist_glom_* math, which is computed on raw, uninverted coordinates.
        ax_left.invert_yaxis()

        ax_right = axes[row_i, 1]
        for seg in SEGMENT_ORDER:
            seg_sub = sub.loc[sub['broad_tubule_marker_call'].astype(str) == seg]
            ax_right.scatter(seg_sub['centroid_x'], seg_sub['centroid_y'], s=4, linewidths=0,
                              color=_segment_color(seg), label=seg)
        ax_right.scatter(hc['centroid_x'], hc['centroid_y'], marker='*', s=60, facecolor='white',
                          edgecolor='black', linewidths=0.8, zorder=5)
        ax_right.set_title(f'{sample}: broad_tubule_marker_call')
        ax_right.set_aspect('equal')
        ax_right.invert_yaxis()
        if row_i == 0:
            ax_right.legend(loc='upper right', fontsize=7, markerscale=2)

    fig.colorbar(sc_handle, ax=list(axes[:, 0]), label='total_scanpy_dpt', shrink=0.6)
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close(fig)


# ---------------------------------------------------------------------------------------------
# Step 8 -- Segment-level sub-segment cross-check (dual method)
# ---------------------------------------------------------------------------------------------

def resolve_subsegment_markers(adata_dpt: ad.AnnData) -> dict:
    resolved = {}
    for family, marker_groups in SUBSEGMENT_MARKER_GROUPS.items():
        resolved[family] = resolve_marker_genes(adata_dpt, marker_groups, label_prefix=f'{family} ',
                                                  fatal_on_empty=False)
        if resolved[family] is None:
            print(f'  WARNING: {family}: at least one sub-panel resolved to zero genes -- '
                  f'skipping Method B for the whole {family} family (Method A falls back to an '
                  f'equal-count split for it).')
    return resolved


def compute_subsegment_scores(adata_dpt: ad.AnnData, resolved_subsegment_markers: dict,
                               random_state: int) -> dict:
    score_cols_by_family = {}
    for family, resolved in resolved_subsegment_markers.items():
        if resolved is None:
            continue
        compute_broad_scores(adata_dpt, resolved, random_state, '{family}_score',
                              f'{family} sub-segment')
        score_cols_by_family[family] = [f'{sp}_score' for sp in resolved]
    return score_cols_by_family


def assign_marker_argmax_labels(dpt_obs: pd.DataFrame, score_cols_by_family: dict) -> pd.Series:
    labels = dpt_obs['broad_tubule_marker_call'].astype(str).copy()
    for family, subpanels in SUBSEGMENT_ORDER.items():
        if family not in score_cols_by_family:
            continue  # Method B skipped for this family; broad label stands in
        mask = dpt_obs['broad_tubule_marker_call'].astype(str) == family
        score_matrix = dpt_obs.loc[mask, [f'{sp}_score' for sp in subpanels]].to_numpy()
        argmax_idx = np.argmax(score_matrix, axis=1)
        labels.loc[mask] = [subpanels[i] for i in argmax_idx]
    return labels


def compute_within_family_orientation(dpt_obs: pd.DataFrame, family: str, marker_labels: pd.Series) -> float:
    order = SUBSEGMENT_ORDER[family]
    mask = dpt_obs['broad_tubule_marker_call'].astype(str) == family
    ordinal = marker_labels.loc[mask].map({seg: i + 1 for i, seg in enumerate(order)})
    dpt_vals = dpt_obs.loc[mask, 'total_scanpy_dpt']
    valid = ordinal.notna()
    if valid.sum() < 10:
        return np.nan
    return float(stats.spearmanr(dpt_vals[valid], ordinal[valid]).statistic)


def assign_pseudospace_labels(dpt_obs: pd.DataFrame, score_cols_by_family: dict,
                               marker_labels: pd.Series):
    """Method A: pseudospace bins with cutpoints dynamically derived from Method B's within-
    family sub-segment proportions (not equal-count tertiles), pooled across all 4 samples.
    PT is assumed correctly oriented (Section 3 orients the global pseudospace axis PT-early by
    construction). TAL/DCT get an orientation check first (C1): if the within-family
    pseudospace-vs-marker-order correlation is negative, the cutpoint order is reversed; if too
    weak (|rho|<0.1), Method A falls back to an equal-count split for that family."""
    labels = dpt_obs['broad_tubule_marker_call'].astype(str).copy()
    dynamic_source = pd.Series('n/a', index=dpt_obs.index)

    for family, order in SUBSEGMENT_ORDER.items():
        mask = dpt_obs['broad_tubule_marker_call'].astype(str) == family
        family_idx = dpt_obs.index[mask]
        dpt_vals = dpt_obs.loc[mask, 'total_scanpy_dpt']

        use_order = order
        if family not in score_cols_by_family:
            reason = 'equal_tertile_fallback'
        elif family == 'PT':
            reason = 'dynamic_from_method_b'
        else:
            rho = compute_within_family_orientation(dpt_obs, family, marker_labels)
            print(f'  {family}: within-family pseudospace-vs-marker-order Spearman rho = {rho:.3f}')
            if not np.isfinite(rho) or abs(rho) < 0.1:
                print(f'  WARNING: {family}: within-family pseudospace ordering too weak for '
                      f'reliable sub-segment definition (|rho|<0.1) -- falling back to '
                      f'equal-count split.')
                reason = 'weak_within_family_orientation_fallback'
            else:
                reason = 'dynamic_from_method_b'
                if rho < 0:
                    use_order = list(reversed(order))

        if reason == 'dynamic_from_method_b':
            proportions = marker_labels.loc[family_idx].value_counts(normalize=True)
            cum = 0.0
            cutpoints = []
            for seg in use_order[:-1]:
                cum += proportions.get(seg, 0.0)
                cutpoints.append(np.percentile(dpt_vals, cum * 100))
            bin_labels = use_order
        else:
            n_groups = len(order)
            edges_pct = np.linspace(0, 100, n_groups + 1)
            cutpoints = [np.percentile(dpt_vals, p) for p in edges_pct[1:-1]]
            bin_labels = order

        bin_edges = [-np.inf] + list(cutpoints) + [np.inf]
        seg_assignment = pd.cut(dpt_vals, bins=bin_edges, labels=bin_labels, include_lowest=True)
        labels.loc[family_idx] = seg_assignment.astype(str)
        dynamic_source.loc[family_idx] = reason

    return labels, dynamic_source


def compute_segment_agreement(dpt_obs: pd.DataFrame, marker_labels: pd.Series,
                               pseudospace_labels: pd.Series, dynamic_source: pd.Series) -> pd.DataFrame:
    rows = []
    for family, order in SUBSEGMENT_ORDER.items():
        mask = dpt_obs['broad_tubule_marker_call'].astype(str) == family
        m_labels, p_labels = marker_labels.loc[mask], pseudospace_labels.loc[mask]
        confusion = pd.crosstab(p_labels, m_labels).reindex(index=order, columns=order, fill_value=0)
        print(f'\nConfusion matrix for {family} (rows=pseudospace, cols=marker):')
        print(confusion.to_string())

        n_cells = int(mask.sum())
        agreement_pct = float(100.0 * np.trace(confusion.to_numpy()) / n_cells) if n_cells else np.nan
        # ARI here is between two labelings with matched class marginals by construction (Method
        # A's cutpoints use Method B's proportions). Interpretation is cell-level agreement, not
        # independent-classifier agreement. A high ARI reflects agreement on which cells belong
        # to each sub-segment, not on the sub-segment proportions themselves.
        ari = float(adjusted_rand_score(m_labels, p_labels)) if n_cells else np.nan
        source_mode = dynamic_source.loc[mask].mode()
        method_a_source = source_mode.iloc[0] if not source_mode.empty else 'n/a'
        rows.append(dict(family=family, agreement_pct=agreement_pct, ari=ari, n_cells=n_cells,
                          method_a_source=method_a_source))
        if agreement_pct < 60 or ari < 0.3:
            print(f'  WARNING: {family} sub-segment definition unstable -- inspect visually '
                  f'(agreement={agreement_pct:.1f}%, ARI={ari:.3f}).')
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------------------------

def main(argv=None):
    args = parse_args(argv)
    project_root = args.project_root.resolve()
    data_root = (args.data_root or project_root / 'data').resolve()
    results_root = (args.results_root or project_root / 'results').resolve()
    run_root = results_root / 'mouse_only_v5'
    output_dir = run_root / OUTPUT_DIRNAME
    output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.random_seed)

    harmony_path = run_root / HARMONY_FILENAME
    dpt_path = run_root / DPT_FILENAME
    if not harmony_path.exists():
        print(f'FATAL: Harmony h5ad not found at {harmony_path}')
        sys.exit(1)
    if not dpt_path.exists():
        print(f'FATAL: DPT h5ad not found at {dpt_path}')
        sys.exit(1)
    for sample in args.samples:
        if sample not in qupath_export.SAMPLE_TO_RAW_GEOJSON:
            print(f'FATAL: no known raw GeoJSON mapping for sample {sample!r}. '
                  f'Known samples: {list(qupath_export.SAMPLE_TO_RAW_GEOJSON.keys())}')
            sys.exit(1)
        if '_' in sample:
            print(f'FATAL: sample name {sample!r} contains an underscore; the join parser '
                  f'splits on the literal "{{sample}}_unit_" prefix, so this would be ambiguous.')
            sys.exit(1)
    if args.dry_run and args.use_manual_annotations:
        print('FATAL: --dry-run and --use-manual-annotations are mutually exclusive (--dry-run '
              'generates seed files for review; --use-manual-annotations skips straight to '
              'reading already-reviewed files).')
        sys.exit(1)

    # Centroids are needed both to re-identify manually-edited QuPath features (Step 1, since
    # QuPath strips all custom properties on export) and for the tissue-edge hull + DPT/classified
    # attachment (Step 2) -- computed once, up front, and reused for both.
    print('=' * 78)
    print('STEP 0 -- RAW CENTROIDS (needed by Step 1 for re-identifying manually-edited features)')
    print('=' * 78)
    centroids = load_all_centroids(data_root, args.samples)

    print('\n' + '=' * 78)
    print('STEP 1 -- GLOMERULUS/VESSEL LABELING (manually curated in QuPath -- see module docstring)')
    print('=' * 78)

    if args.use_manual_annotations:
        print('\n' + '!' * 78)
        print('!!! --use-manual-annotations is set: skipping the automated seed reconstruction '
              'and manual_annotation_diff.txt generation entirely. Reading pre-existing '
              '*_glomerulus_vessel_manual_review.geojson files directly as ground truth. !!!')
        print('!' * 78)
        print()
        classified, _dropped_df = read_manual_glomerulus_labels(output_dir, args.samples, centroids)
    else:
        adata_h = ad.read_h5ad(harmony_path)
        for col in ('broad_tubule_marker_call', 'celltype_primary', 'sample'):
            if col not in adata_h.obs.columns:
                print(f'FATAL: Harmony obs missing required column {col!r}.')
                sys.exit(1)

        resolved_disambig = resolve_marker_genes(adata_h, DISAMBIGUATION_MARKERS)
        compute_broad_scores(adata_h, resolved_disambig, args.random_seed,
                              '{family}_broad_score', 'Glomerulus/Vessel')

        glom_mask = adata_h.obs['broad_tubule_marker_call'].astype(str) == 'Glomerulus'
        glom_obs = adata_h.obs.loc[glom_mask].copy()
        if glom_obs.empty:
            print('FATAL: no cells with broad_tubule_marker_call == "Glomerulus" found in Harmony obs.')
            sys.exit(1)

        seed_df = prepare_glomerulus_seed_labels(glom_obs, args.samples)

        # One file per sample: coordinates are sample-specific pixel space, so a single file
        # combining all 4 samples' features cannot be meaningfully imported into any one QuPath
        # image (each image would show three other samples' tubules overlaid on it too). A sample
        # whose review file already exists is never overwritten, so manual corrections are never
        # clobbered by re-running this step.
        print()
        manual_paths = {sample: output_dir / MANUAL_REVIEW_FILENAME_TEMPLATE.format(sample=sample)
                         for sample in args.samples}
        freshly_seeded = set()
        for sample in args.samples:
            review_path = manual_paths[sample]
            if review_path.exists():
                print(f'  {sample}: {review_path.name} already exists -- leaving it untouched '
                      f'(may contain your manual corrections).')
                continue
            n_written = write_labeled_geojson(
                data_root, review_path, seed_df, [sample],
                label_column='glom_qc_class', color_map=GLOM_VESSEL_SEED_COLORS,
                extra_property_columns={'glom_vs_vessel_gap': 'glom_vs_vessel_gap'},
            )
            print(f'  {sample}: wrote {n_written} seed features to {review_path}')
            freshly_seeded.add(sample)

        print('\n' + '=' * 78)
        print('QuPath usage:')
        print('For each sample, open that sample\'s histology image in QuPath, then File -> '
              'Object data -> Import objects, choose the matching '
              '{sample}_glomerulus_vessel_manual_review.geojson. Correct any misclassified '
              'Glomeruli/Vessel features by hand (renaming the class or moving a feature to a '
              'different class is fine), then save back to the exact same path/filename.')

        if args.dry_run:
            print('\nCorrect each sample\'s *_glomerulus_vessel_manual_review.geojson in QuPath '
                  'and save it back to the same path, then re-run without --dry-run to use your '
                  'corrected labels for the rest of the pipeline.')
            return

        if freshly_seeded:
            print('\n' + '!' * 78)
            print('!!! WARNING: the following sample(s) had no existing manual-review file, so '
                  'one was just generated from celltype_primary and is being used AS-IS, '
                  'unreviewed: ' + ', '.join(sorted(freshly_seeded)))
            print('!' * 78)

        print()
        classified, dropped_df = read_manual_glomerulus_labels(output_dir, args.samples, centroids)

        diff_path = output_dir / DIFF_REPORT_FILENAME
        print('\n' + '=' * 78)
        print('MANUAL ANNOTATION DIFF REPORT')
        print('=' * 78)
        write_manual_annotation_diff_report(seed_df, classified, dropped_df, args.samples,
                                             manual_paths, diff_path)
        print(f'\nWrote {diff_path}')
        print("\nUsing user's manually-edited glomerulus labels from "
              "*_glomerulus_vessel_manual_review.geojson (not the script's auto-classification). "
              "See manual_annotation_diff.txt for the differences.")

        # Free the large (~2.6GB resident) Harmony expression matrix now that Step 1 is done --
        # `classified` (built from the manual-review GeoJSON, not from adata_h) is all that's
        # needed downstream.
        del adata_h

    print('\n' + '=' * 78)
    print('STEP 2 -- TISSUE EDGE + CENTROID ATTACHMENT')
    print('=' * 78)
    tissue_polys = {}
    for sample in args.samples:
        xy = np.array(list(centroids[sample].values()))
        tissue_polys[sample] = compute_tissue_edge_polygon(xy)

    adata_dpt = ad.read_h5ad(dpt_path)
    for col in ('sample', 'broad_tubule_marker_call', 'total_scanpy_dpt'):
        if col not in adata_dpt.obs.columns:
            print(f'FATAL: DPT obs missing required column {col!r}.')
            sys.exit(1)

    adata_dpt.obs = attach_centroids(adata_dpt.obs, args.samples, centroids)
    classified = attach_centroids(classified, args.samples, centroids)

    print('\n' + '=' * 78)
    print('STEP 3 -- DISTANCE COMPUTATIONS')
    print('=' * 78)
    for sample in args.samples:
        dpt_mask = adata_dpt.obs['sample'].astype(str) == sample
        dpt_xy = adata_dpt.obs.loc[dpt_mask, ['centroid_x', 'centroid_y']].to_numpy()

        hc_mask = ((classified['sample'].astype(str) == sample)
                   & (classified['glom_qc_class'] == 'Glomeruli'))
        hc_xy = classified.loc[hc_mask, ['centroid_x', 'centroid_y']].to_numpy()

        dist_nearest, dist_k3 = compute_glom_distances(dpt_xy, hc_xy, args.k_glom, sample)
        dist_edge = compute_edge_distances(dpt_xy, tissue_polys[sample])

        idx = adata_dpt.obs.index[dpt_mask]
        adata_dpt.obs.loc[idx, 'dist_glom_nearest'] = dist_nearest
        adata_dpt.obs.loc[idx, 'dist_glom_k3'] = dist_k3
        adata_dpt.obs.loc[idx, 'dist_edge'] = dist_edge
        print(f'  {sample}: n_tubules={int(dpt_mask.sum())}  n_high_confidence_glom={len(hc_xy)}  '
              f'dist_glom_nearest median={np.median(dist_nearest):.1f}  '
              f'dist_edge median={np.median(dist_edge):.1f}')

    print('\n' + '=' * 78)
    print('STEP 4 -- PER-SAMPLE SPEARMAN + PERMUTATION + BOOTSTRAP')
    print('=' * 78)
    results_df = run_all_sample_tests(adata_dpt.obs, args.samples, args.n_permutations, rng)
    results_path = output_dir / 'spatial_validation_results.csv'
    results_df.to_csv(results_path, index=False)
    print(f'Wrote {results_path}')

    print('\n' + '=' * 78)
    print('STEP 5 -- POOLED FISHER-Z CORRELATIONS')
    print('=' * 78)
    pooled_df = compute_pooled_correlations(results_df)
    pooled_path = output_dir / 'pooled_correlations.csv'
    pooled_df.to_csv(pooled_path, index=False)
    print(f'Wrote {pooled_path}')

    print('\n' + '=' * 78)
    print('STEP 6 -- SEGMENT-ORDER TEST')
    print('=' * 78)
    segment_df = segment_order_test(adata_dpt.obs, args.samples)
    segment_path = output_dir / 'segment_order_test.csv'
    segment_df.to_csv(segment_path, index=False)
    print(f'Wrote {segment_path}')

    print('\n' + '=' * 78)
    print('STEP 7 -- FIGURES')
    print('=' * 78)
    fig_paths = []
    p = output_dir / 'fig1_pseudospace_vs_glom_scatter.png'
    fig_pseudospace_vs_glom_scatter(adata_dpt.obs, args.samples, results_df, p)
    fig_paths.append(p)
    p = output_dir / 'fig2_pseudospace_vs_edge_scatter.png'
    fig_pseudospace_vs_edge_scatter(adata_dpt.obs, args.samples, results_df, p)
    fig_paths.append(p)
    p = output_dir / 'fig3_segment_order_boxplot.png'
    fig_segment_order_boxplot(adata_dpt.obs, args.samples, segment_df, p)
    fig_paths.append(p)
    p = output_dir / 'fig4_forest_plot.png'
    fig_forest_plot(results_df, pooled_df, p)
    fig_paths.append(p)
    p = output_dir / 'fig5_qc_map.png'
    fig_qc_map(adata_dpt.obs, classified, args.samples, p)
    fig_paths.append(p)
    for p in fig_paths:
        print(f'Wrote {p}')

    print('\n' + '=' * 78)
    print('STEP 8 -- SEGMENT-LEVEL SUB-SEGMENT CROSS-CHECK')
    print('=' * 78)
    resolved_subsegment = resolve_subsegment_markers(adata_dpt)
    score_cols_by_family = compute_subsegment_scores(adata_dpt, resolved_subsegment, args.random_seed)
    marker_labels = assign_marker_argmax_labels(adata_dpt.obs, score_cols_by_family)
    pseudospace_labels, dynamic_source = assign_pseudospace_labels(adata_dpt.obs, score_cols_by_family,
                                                                     marker_labels)
    adata_dpt.obs['segment_label_marker'] = marker_labels
    adata_dpt.obs['segment_label_pseudospace'] = pseudospace_labels
    adata_dpt.obs['dynamic_source'] = dynamic_source

    agreement_df = compute_segment_agreement(adata_dpt.obs, marker_labels, pseudospace_labels, dynamic_source)
    agreement_path = output_dir / 'segment_agreement_summary.csv'
    agreement_df.to_csv(agreement_path, index=False)
    print(f'Wrote {agreement_path}')

    pt_row = agreement_df.loc[agreement_df['family'] == 'PT']
    if not pt_row.empty and pt_row.iloc[0]['agreement_pct'] < 60:
        print('\n' + '!' * 78)
        print(f"!!! WARNING: PT sub-segment agreement is {pt_row.iloc[0]['agreement_pct']:.1f}% "
              f"(< 60%) -- primary validation target unstable.")
        print('!' * 78)

    # One file per sample, per method -- same reasoning as Step 1's QC GeoJSON: QuPath needs
    # per-image files, coordinates are sample-specific pixel space.
    print()
    for sample in args.samples:
        marker_path = output_dir / f'{sample}_segment_labels_marker.geojson'
        n_marker = write_labeled_geojson(
            data_root, marker_path, adata_dpt.obs, [sample],
            label_column='segment_label_marker', color_map=SUBSEGMENT_COLORS,
            extra_property_columns={'total_scanpy_dpt': 'total_scanpy_dpt'},
            extra_constant_properties={'method': 'marker_argmax'},
        )
        print(f'Wrote {n_marker} features to {marker_path}')

        pseudospace_path = output_dir / f'{sample}_segment_labels_pseudospace.geojson'
        n_pseudospace = write_labeled_geojson(
            data_root, pseudospace_path, adata_dpt.obs, [sample],
            label_column='segment_label_pseudospace', color_map=SUBSEGMENT_COLORS,
            extra_property_columns={'total_scanpy_dpt': 'total_scanpy_dpt', 'dynamic_source': 'dynamic_source'},
            extra_constant_properties={'method': 'pseudospace_dynamic'},
        )
        print(f'Wrote {n_pseudospace} features to {pseudospace_path}')

    print('\nFor each sample, load that sample\'s {sample}_segment_labels_pseudospace.geojson '
          'and {sample}_segment_labels_marker.geojson side by side in QuPath (or as two overlays '
          'on the same image). Regions where the two disagree highlight boundary tubules or '
          'misplacements on the pseudospace axis -- these are diagnostic, not errors.')

    print('\n' + '=' * 78)
    print('DONE. Files written under', output_dir)
    print('=' * 78)


if __name__ == '__main__':
    main()
