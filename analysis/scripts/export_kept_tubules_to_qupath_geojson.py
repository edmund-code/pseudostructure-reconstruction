"""Export QC-kept, cluster-labeled mouse tubules back to per-sample GeoJSON for QuPath overlay.

Reads the labeled AnnData from the configured `results-root` and the original
raw per-sample tubule segmentation GeoJSON files, and writes one GeoJSON per sample containing
only the QC-kept tubules, annotated with a QuPath-readable classification (name + color) derived
from obs['broad_tubule_marker_call']. Geometry is passed through unmodified.

Output format: QuPath GeoJSON, i.e. a FeatureCollection of Features whose properties carry
objectType / isLocked / measurements (a list of {name, value}) and a classification of the form
{"name": <label>, "colorRGB": <signed int32 ARGB>}. The raw segmentation files already use this
shape, so only the classification is retargeted; assert_qupath_conformant() re-checks the whole
payload before anything is written. Note colorRGB, NOT an [r, g, b] array under "color" -- QuPath
does not read that, and a file using it imports with the classification stripped of its color.

The label source is the Harmony object, not the DPT object: Glomerulus and Vessel tubules are
excluded from the DPT continuum (the DPT artifact only contains
the five tubule families + Transitional), but they are legitimate labeling-control classes that
should still appear on the QuPath overlay -- the Harmony object (updated in place by notebook 6
Section 2 with the finalized cluster-annotation columns, without re-running the Harmony
integration itself) is the only saved artifact that still carries labels for all seven classes.

Read-only on all inputs. Only writes under the configured results root.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import anndata as ad
import pandas as pd

# The requested raw-segmentation path pattern (data/segmentations/raw_segmentation/...) does not
# exist on disk. These are the actual raw files, matching notebook 1's NEW_SEGMENTATION_SAMPLE_MAP.
#
# These MUST be the same files notebook 1 Part 1 enumerated to assign obs['feature_index'] -- the
# v4 segmentations. The older *_processed.geojson files are a different segmentation entirely
# (e.g. Ctrl_1A2: 21,951 features vs v4's 13,630, and classes 'Tubules'/'Vessel' rather than
# 'tubule_proximal'). Because they hold MORE features than v4, every v4 feature_index still
# lands in range, so a wrong-file join passes the range/uniqueness checks below and silently
# attaches correct labels to the wrong polygons. verify_centroids_match() is what actually
# catches it.
SAMPLE_TO_RAW_GEOJSON = {
    'Ctrl1A2': 'Ctrl_1A2_v4.geojson',
    'Ctrl1A4': 'Ctrl_1A4_v4.geojson',
    'IR2A2': 'IR_2A2_v4.geojson',
    'IR2A4': 'IR_2A4_v4.geojson',
}

DEFAULT_SAMPLES = list(SAMPLE_TO_RAW_GEOJSON.keys())
# Defaults track the CURRENT run (results/mouse_only_v5/), which notebook 6 also assigns onto
# these names before calling main(). They previously pointed into a retired result location -- a
# superseded location still holding Jul-2026 exports built from the old *_processed.geojson
# segmentation. Leaving the defaults there meant a standalone run either failed on a missing
# input or, worse, left stale wrong-polygon files sitting next to the real ones to be picked up
# by hand. Pass --project-root and reassign these if you need a different run.
HARMONY_RELATIVE_PATH = Path('results/mouse_only_v5/all_mouse_tubules_harmony_pass1.h5ad')
OUTPUT_RELATIVE_DIR = Path('results/mouse_only_v5/qupath_geojson')

# Keyed on the notebook's coarse vocabulary: COARSE_ORDER (PT, DTL, AL, DCT, CNT_CD) plus the
# REMOVE_CLASSES that the pass-1 object still carries. 'AL' is the ascending limb (ATL/mTAL/cTAL/
# macula densa are not separable at this resolution) and 'DTL' the descending thin limb -- these
# replaced the older 'TAL' / 'TL' spellings, which are kept below as spelling traps.
LABEL_COLORS = {
    # nephron continuum, in anatomical order
    'PT': [230, 25, 75],
    'DTL': [245, 130, 48],
    'AL': [60, 180, 75],
    'DCT': [0, 130, 200],
    'CNT_CD': [145, 30, 180],
    # non-tubule / control classes
    'Glomerulus': [255, 225, 25],
    'Vessel': [0, 255, 255],
    'Stroma': [170, 110, 40],
    'SmoothMuscle': [240, 50, 230],
    'Immune': [128, 128, 0],
    'Unassigned': [128, 128, 128],
    'Transitional': [160, 160, 160],
}
# Spelling variants that would silently miss LABEL_COLORS and fall back to the hash color even
# though a canonical entry already exists for the intended category -- worth a loud warning.
LABEL_SPELLING_TRAPS = {
    'TAL': 'AL',
    'TL': 'DTL',
    'Thin_Limb': 'DTL',
    'Ambiguous': 'Unassigned',
}


def hash_fallback_color(label: str) -> list[int]:
    digest = hashlib.md5(label.encode('utf-8')).digest()
    # Keep channels away from the extremes (0 / 255) so fallback colors stay visually distinct
    # from pure black/white and from each other.
    return [32 + (b % 192) for b in digest[:3]]


def color_for_label(label: str) -> list[int]:
    if label in LABEL_COLORS:
        return LABEL_COLORS[label]
    return hash_fallback_color(label)


def pack_color_rgb(rgb: list[int]) -> int:
    """Pack [r, g, b] into the signed 32-bit ARGB int QuPath stores as classification.colorRGB.

    QuPath's PathClass GeoJSON field is 'colorRGB', an integer produced by ColorTools.packARGB
    with alpha forced to 0xFF -- so it is ALWAYS negative when read back as a signed int (the
    alpha byte sets the sign bit). A [r, g, b] array under a 'color' key is not a format QuPath
    reads inside 'classification': the class imports with no color and QuPath silently assigns
    its own, which is what made the previous output look wrong on import. The raw segmentation
    GeoJSONs this script copies from already use colorRGB (e.g. -65281 = opaque magenta for the
    'Tubules' class), so this restores the convention rather than inventing one.
    """
    r, g, b = (int(c) & 0xFF for c in rgb)
    packed = (0xFF << 24) | (r << 16) | (g << 8) | b
    return packed - (1 << 32) if packed >= (1 << 31) else packed


DEFAULT_PROJECT_ROOT = Path(__file__).resolve().parent.parent


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root', type=Path, default=DEFAULT_PROJECT_ROOT,
                         help=f'Repo root (default: {DEFAULT_PROJECT_ROOT}, derived from this '
                              f'script\'s own location so it works regardless of cwd).')
    parser.add_argument('--samples', nargs='+', default=DEFAULT_SAMPLES)
    parser.add_argument('--label-column', default='broad_tubule_marker_call')
    parser.add_argument('--data-root', type=Path,
                        help='Directory containing the private v4 GeoJSON files (default: '
                             '<project-root>/data).')
    parser.add_argument('--results-root', type=Path,
                        help='Directory containing mouse_only_v5 outputs (default: '
                             '<project-root>/results).')
    parser.add_argument('--dry-run', action='store_true',
                         help='Run identity-key inspection + join validation only; write no files.')
    return parser.parse_args(argv)


def load_raw_geojson(data_root: Path, sample: str) -> dict:
    raw_path = data_root / SAMPLE_TO_RAW_GEOJSON[sample]
    if not raw_path.exists():
        raise FileNotFoundError(f'Raw GeoJSON for sample {sample!r} not found at {raw_path}')
    return json.loads(raw_path.read_text())


def parse_feature_index(obs_name: str, sample: str) -> int:
    prefix = f'{sample}_unit_'
    if not obs_name.startswith(prefix):
        raise ValueError(f'obs_name {obs_name!r} does not start with expected prefix {prefix!r}')
    suffix = obs_name[len(prefix):]
    if not suffix.isdigit():
        raise ValueError(f'obs_name {obs_name!r} suffix {suffix!r} is not a plain integer')
    return int(suffix)


def verify_centroids_match(obs: pd.DataFrame, sample: str, data_root: Path,
                            n_check: int = 400, tol: float = 1.0) -> tuple[int, int, float]:
    """Prove the GeoJSON being read is the one obs['feature_index'] was assigned against.

    Notebook 1 Part 1 stores each tubule's polygon centroid in obs as x_centroid / y_centroid at
    the same time it assigns feature_index by enumerating the segmentation file. Re-deriving the
    centroid of feature[feature_index] and comparing is therefore a direct test of the join: it
    agrees to sub-pixel precision on the right file and diverges wildly on any other.

    This exists because the range/uniqueness checks cannot catch a wrong file. Every older
    *_processed.geojson holds MORE features than its v4 counterpart, so v4 indices all resolve,
    stay unique, and produce a clean-looking report while pointing at unrelated polygons.
    """
    from shapely.geometry import shape

    sample_obs = obs[obs['sample'].astype(str) == sample]
    if not {'x_centroid', 'y_centroid'}.issubset(sample_obs.columns):
        print(f'  {sample}: obs has no x_centroid/y_centroid -- centroid check SKIPPED')
        return 0, 0, float('nan')

    features = load_raw_geojson(data_root, sample).get('features', [])
    stride = max(1, len(sample_obs) // n_check)
    checked = 0
    mismatches = []
    max_dev = 0.0
    for obs_name, row in list(sample_obs.iloc[::stride].iterrows())[:n_check]:
        feature_index = parse_feature_index(obs_name, sample)
        if not (0 <= feature_index < len(features)):
            mismatches.append((obs_name, feature_index, float('inf')))
            continue
        geometry = shape(features[feature_index]['geometry'])
        if not geometry.is_valid:
            geometry = geometry.buffer(0)
        if geometry.is_empty:
            continue
        centroid = geometry.centroid
        deviation = ((centroid.x - float(row['x_centroid'])) ** 2
                     + (centroid.y - float(row['y_centroid'])) ** 2) ** 0.5
        checked += 1
        max_dev = max(max_dev, deviation)
        if deviation > tol:
            mismatches.append((obs_name, feature_index, deviation))

    return checked, len(mismatches), max_dev


def inspect_and_validate(obs: pd.DataFrame, samples: list[str], data_root: Path, label_column: str):
    print('=' * 78)
    print('STEP 1 -- IDENTITY-KEY INSPECTION')
    print('=' * 78)
    print('obs.columns:', obs.columns.tolist())
    print('obs_names[:3]:', obs.index[:3].tolist())
    print("unique obs['sample']:", sorted(obs['sample'].astype(str).unique().tolist()))
    print("\nobs['sample'].value_counts():")
    print(obs['sample'].value_counts().to_string())
    print(f"\nobs[{label_column!r}].value_counts(dropna=False):")
    print(obs[label_column].value_counts(dropna=False).to_string())

    for sample in samples:
        if '_' in sample:
            print(f'\nFATAL: sample name {sample!r} contains an underscore; the join parser splits '
                  f'on the literal "{{sample}}_unit_" prefix, so this would be ambiguous.')
            sys.exit(1)
        sample_obs_names = obs.index[obs['sample'].astype(str) == sample]
        print(f"\nFirst 3 obs_names for sample {sample!r}: {sample_obs_names[:3].tolist()}")

    first_sample = samples[0]
    raw_payload = load_raw_geojson(data_root, first_sample)
    features = raw_payload.get('features', [])
    print(f'\nRaw GeoJSON for {first_sample!r}: top-level keys = {list(raw_payload.keys())}, '
          f'n_features = {len(features)}')
    for i, feature in enumerate(features[:2]):
        print(f'  feature {i}: keys={list(feature.keys())}, properties={feature.get("properties")}')

    print('\nChosen join strategy: obs_name -> strip "{sample}_unit_" prefix -> int(feature_index) '
          '-> positional (enumerate) index into that sample\'s raw GeoJSON "features" list.')
    print('(No raw GeoJSON in this repo has properties["id"]; notebook 1\'s load_raw_segments always '
          'falls back to unit_id = f"unit_{feature_index}", so this is the proven join key, not a guess.)')

    print('\n' + '=' * 78)
    print('STEP 1 -- JOIN VALIDATION (must be 100% resolved before any file is written)')
    print('=' * 78)
    any_failure = False
    n_features_by_sample = {}
    for sample in samples:
        sample_mask = obs['sample'].astype(str) == sample
        kept_obs_names = obs.index[sample_mask].tolist()
        raw_payload = load_raw_geojson(data_root, sample)
        n_features = len(raw_payload.get('features', []))
        n_features_by_sample[sample] = n_features

        mismatches = []
        seen_indices = set()
        for obs_name in kept_obs_names:
            try:
                feature_index = parse_feature_index(obs_name, sample)
            except ValueError as exc:
                mismatches.append(f'{obs_name}: unparseable ({exc})')
                continue
            if not (0 <= feature_index < n_features):
                mismatches.append(f'{obs_name}: feature_index {feature_index} out of range '
                                   f'[0, {n_features})')
                continue
            if feature_index in seen_indices:
                mismatches.append(f'{obs_name}: duplicate feature_index {feature_index}')
                continue
            seen_indices.add(feature_index)

        matched = len(seen_indices)
        kept = len(kept_obs_names)
        status = 'OK' if matched == kept and not mismatches else 'MISMATCH'
        print(f'  {sample}: kept_obs={kept}  matched={matched}  raw_features={n_features}  [{status}]')
        if mismatches:
            any_failure = True
            print(f'    {len(mismatches)} mismatch(es) for {sample}:')
            for m in mismatches[:20]:
                print('      -', m)
            if len(mismatches) > 20:
                print(f'      ... and {len(mismatches) - 20} more')

    if any_failure:
        print('\nFATAL: join validation failed for at least one sample. Aborting before writing '
              'any output -- a wrong join must never produce partial/incorrect overlays.')
        sys.exit(1)

    print('\nJoin validation passed for all samples: every kept obs tubule resolves to a unique, '
          'in-range GeoJSON feature.')

    print('\n' + '=' * 78)
    print('STEP 1b -- CENTROID VERIFICATION (proves it is the RIGHT segmentation file)')
    print('=' * 78)
    centroid_failure = False
    for sample in samples:
        checked, n_bad, max_dev = verify_centroids_match(obs, sample, data_root)
        if not checked:
            continue
        source = SAMPLE_TO_RAW_GEOJSON[sample]
        if n_bad:
            centroid_failure = True
            print(f'  {sample}: {n_bad}/{checked} sampled tubules DISAGREE with {source} '
                  f'(max deviation {max_dev:.1f} px)  [MISMATCH]')
        else:
            print(f'  {sample}: {checked}/{checked} sampled centroids match {source} '
                  f'(max deviation {max_dev:.3g} px)  [OK]')

    if centroid_failure:
        print('\nFATAL: obs centroids do not match the segmentation being read, so '
              "obs['feature_index'] indexes a DIFFERENT file than the one loaded here. Point "
              'SAMPLE_TO_RAW_GEOJSON at the segmentation notebook 1 Part 1 actually used (the '
              '*_v4.geojson files) and re-run. Aborting before writing any output.')
        sys.exit(1)

    return n_features_by_sample


def assert_qupath_conformant(payload: dict, sample: str) -> None:
    """Fail before writing if the payload deviates from what QuPath's importer reads.

    Checked against a known-good QuPath export: a FeatureCollection whose every feature carries a
    string id, a Polygon/MultiPolygon geometry with closed rings, and properties holding
    objectType, isLocked, a {name, colorRGB:int} classification, and measurements as a list of
    {name, value} objects. The colorRGB check is the one that matters most -- an [r,g,b] array
    under a 'color' key is accepted by json.dump and rejected by QuPath, so it can only be caught
    here, not by a schema-free write.
    """
    problems = []
    if payload.get('type') != 'FeatureCollection':
        problems.append(f"top-level type is {payload.get('type')!r}, expected 'FeatureCollection'")

    seen_ids = set()
    for i, feature in enumerate(payload.get('features', [])):
        where = f'feature {i}'
        if not isinstance(feature.get('id'), str):
            problems.append(f'{where}: id is {type(feature.get("id")).__name__}, expected str')
        elif feature['id'] in seen_ids:
            problems.append(f'{where}: duplicate id {feature["id"]!r}')
        else:
            seen_ids.add(feature['id'])

        geometry = feature.get('geometry') or {}
        if geometry.get('type') not in ('Polygon', 'MultiPolygon'):
            problems.append(f'{where}: geometry type {geometry.get("type")!r} not Polygon/MultiPolygon')
        else:
            polygons = ([geometry['coordinates']] if geometry['type'] == 'Polygon'
                        else geometry['coordinates'])
            for polygon in polygons:
                for ring in polygon:
                    if len(ring) < 4 or ring[0] != ring[-1]:
                        problems.append(f'{where}: unclosed or degenerate ring')
                        break

        props = feature.get('properties') or {}
        if props.get('objectType') not in ('annotation', 'detection', 'cell', 'tile'):
            problems.append(f'{where}: objectType {props.get("objectType")!r} not a QuPath type')
        if not isinstance(props.get('isLocked'), bool):
            problems.append(f'{where}: isLocked is not a bool')

        classification = props.get('classification')
        if not isinstance(classification, dict):
            problems.append(f'{where}: classification missing')
        else:
            if not isinstance(classification.get('name'), str) or not classification['name']:
                problems.append(f'{where}: classification.name missing or not a str')
            color = classification.get('colorRGB')
            if not isinstance(color, int) or isinstance(color, bool):
                problems.append(f'{where}: classification.colorRGB is '
                                f'{type(color).__name__}, expected int '
                                f'(a [r,g,b] list under "color" is NOT read by QuPath)')
            elif not (-(1 << 31) <= color < (1 << 31)):
                problems.append(f'{where}: classification.colorRGB {color} outside signed int32')

        measurements = props.get('measurements')
        if measurements is not None:
            if not isinstance(measurements, list):
                problems.append(f'{where}: measurements is {type(measurements).__name__}, '
                                f'expected a list of {{name, value}} objects')
            else:
                for m in measurements:
                    if not (isinstance(m, dict) and isinstance(m.get('name'), str)
                            and isinstance(m.get('value'), (int, float))):
                        problems.append(f'{where}: malformed measurement {m!r}')
                        break

        if len(problems) > 20:
            break

    if problems:
        print(f'\nFATAL: [{sample}] output is not QuPath-conformant; refusing to write.')
        for p in problems[:20]:
            print('  -', p)
        sys.exit(1)


def export_sample(data_root: Path, output_dir: Path, sample: str, obs: pd.DataFrame,
                   label_column: str, n_features_input: int) -> dict:
    sample_mask = obs['sample'].astype(str) == sample
    sample_obs = obs.loc[sample_mask, [label_column]]

    n_missing_label = int(sample_obs[label_column].isna().sum())
    if n_missing_label:
        print(f'  [{sample}] skipping {n_missing_label} tubule(s) with missing {label_column!r}')

    kept_dict: dict[int, str] = {}
    for obs_name, label in sample_obs[label_column].items():
        if pd.isna(label):
            continue
        feature_index = parse_feature_index(obs_name, sample)
        kept_dict[feature_index] = str(label)

    present_labels = sample_obs[label_column].dropna().astype(str)
    label_counts = present_labels.value_counts()
    unrecognized = [lbl for lbl in label_counts.index if lbl not in LABEL_COLORS]
    if unrecognized:
        print(f'  [{sample}] labels not in LABEL_COLORS (will use MD5 fallback color):')
        for lbl in unrecognized:
            print(f'      {lbl!r}: {int(label_counts[lbl])}')
            if lbl in LABEL_SPELLING_TRAPS:
                print(f'      *** WARNING: {lbl!r} looks like a spelling variant of the canonical '
                      f'label {LABEL_SPELLING_TRAPS[lbl]!r} -- check for a category-name mismatch '
                      f'against LABEL_COLORS before trusting the fallback color. ***')

    raw_payload = load_raw_geojson(data_root, sample)
    kept_features = []
    for feature_index, feature in enumerate(raw_payload.get('features', [])):
        if feature_index not in kept_dict:
            continue
        label = kept_dict[feature_index]
        properties = dict(feature.get('properties', {}))
        # objectType / isLocked / measurements come through from the raw file already in QuPath's
        # shape; only the classification is retargeted from 'Tubules' to the cluster label.
        properties['objectType'] = properties.get('objectType', 'annotation')
        properties.setdefault('isLocked', False)
        properties['classification'] = {'name': label,
                                        'colorRGB': pack_color_rgb(color_for_label(label))}
        properties['name'] = label
        properties['kept_tubule'] = True
        properties['source_sample'] = sample
        # feature['id'] (top-level, sibling of 'properties') and 'geometry' are copied through
        # unchanged below -- only 'properties' is replaced.
        new_feature = dict(feature)
        new_feature['properties'] = properties
        kept_features.append(new_feature)

    output_payload = {k: v for k, v in raw_payload.items() if k != 'features'}
    output_payload['type'] = 'FeatureCollection'
    output_payload['features'] = kept_features

    assert_qupath_conformant(output_payload, sample)

    output_path = output_dir / f'{sample}_kept_tubules_labeled.geojson'
    with output_path.open('w', encoding='utf-8') as f:
        json.dump(output_payload, f, ensure_ascii=False)

    written_label_counts = pd.Series([f['properties']['name'] for f in kept_features]).value_counts()
    return {
        'sample': sample,
        'n_features_input': n_features_input,
        'n_kept_obs': len(kept_dict),
        'n_features_written': len(kept_features),
        'label_counts': written_label_counts.to_dict(),
        'output_path': output_path,
    }


def main(argv=None):
    args = parse_args(argv)
    project_root: Path = args.project_root.resolve()
    data_root = (args.data_root or project_root / 'data').resolve()
    results_root = (args.results_root or project_root / 'results').resolve()
    harmony_path = results_root / 'mouse_only_v5' / HARMONY_RELATIVE_PATH.name
    if not harmony_path.exists():
        print(f'FATAL: input AnnData not found at {harmony_path}')
        sys.exit(1)

    adata = ad.read_h5ad(harmony_path, backed='r')
    obs = adata.obs
    if 'broad_tubule_marker_call' not in obs.columns:
        print(f'FATAL: {harmony_path.name} has no broad_tubule_marker_call column. '
              f'Re-run notebook 6 Section 2, which updates this file in place with the finalized '
              f'cluster-annotation columns after annotation.')
        sys.exit(1)

    for sample in args.samples:
        if sample not in SAMPLE_TO_RAW_GEOJSON:
            print(f'FATAL: no known raw GeoJSON mapping for sample {sample!r}. '
                  f'Known samples: {list(SAMPLE_TO_RAW_GEOJSON.keys())}')
            sys.exit(1)

    n_features_by_sample = inspect_and_validate(obs, args.samples, data_root, args.label_column)

    if args.dry_run:
        print('\n--dry-run: no files written.')
        return

    output_dir = results_root / 'mouse_only_v5' / OUTPUT_RELATIVE_DIR.name
    output_dir.mkdir(parents=True, exist_ok=True)

    print('\n' + '=' * 78)
    print('STEP 2/3 -- PER-SAMPLE EXPORT')
    print('=' * 78)
    reports = []
    for sample in args.samples:
        print(f'\n[{sample}]')
        report = export_sample(data_root, output_dir, sample, obs, args.label_column,
                               n_features_by_sample[sample])
        reports.append(report)
        print(f'  n_features_input={report["n_features_input"]}  '
              f'n_kept_obs={report["n_kept_obs"]}  '
              f'n_features_written={report["n_features_written"]}')
        print(f'  label counts written: {report["label_counts"]}')
        print(f'  output: {report["output_path"]}')

    print('\n' + '=' * 78)
    print('QuPath usage:')
    for report in reports:
        sample = report['sample']
        print(f'In QuPath, open the sample\'s histology image, then File → Object data → '
              f'Import objects, choose the matching {sample}_kept_tubules_labeled.geojson. '
              f'Class colors and names render on the tubule outlines.')


if __name__ == '__main__':
    main()
