"""Notebook 70: run the frozen ``pseudospace_reconstruction`` package (stages 01-04) and describe its PT coordinate.

The package assumes one gene naming system. For the human vs mouse dataset, PIPELINE_GUIDE.md says to map both
species to shared one-to-one ortholog names before stage 02. This module does that mapping and the bookkeeping
around the package run. It never patches the package.

- ``shared_ortholog_genes``: the accepted HCOP one-to-one pairs (``ortholog_map_used.csv``) whose mouse symbol is
  a stage-01 gene of every mouse sample and whose human symbol is a stage-01 gene of every human sample.
- ``to_shared_genes``: one sample's stage-01 units restricted to those genes, under mouse symbols.
- ``write_shared_bin_dir``: a Space Ranger bin directory holding the same bins on the shared genes. Stage 02 needs
  it because its split test (doublet tier 2) re-reads each sample's bin matrix, and ``build_bin_table`` requires
  that the bin totals sum to the unit totals.
- ``match_polygons``: pairs old and new segmentation polygons by IoU, to carry old structure ids over.
- ``adjacent_aucs``, ``pt_boundary``, ``geojson_schema``, ``uns_safe``: descriptive checks and export helpers.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

MATRIX_FILE = 'filtered_feature_bc_matrix.h5'


# ---------------------------------------------------------------------------------------------------
# shared one-to-one orthologs
# ---------------------------------------------------------------------------------------------------
def validate_ortholog_map(orthologs: pd.DataFrame) -> pd.DataFrame:
    """The accepted map as ``pathway_inputs.rebuild_pt_expression`` checks it: nonempty, one-to-one in both columns
    (case-insensitive), no blank symbol. Returns the two symbol columns as strings."""
    columns = ['mouse_symbol', 'human_symbol']
    if (not set(columns + ['mapping_status']).issubset(orthologs) or orthologs.empty
            or orthologs[columns].isna().any().any()
            or any(orthologs[c].astype(str).str.upper().duplicated().any() for c in columns)
            or any(orthologs[c].astype(str).str.strip().eq('').any() for c in columns)):
        raise ValueError('Expected a nonempty accepted one-to-one ortholog map.')
    return orthologs[columns].astype(str).reset_index(drop=True)


def shared_ortholog_genes(orthologs: pd.DataFrame, mouse_genes, human_genes) -> pd.DataFrame:
    """Shared one-to-one orthologs of the stage-01 gene lists.

    ``mouse_genes`` and ``human_genes`` are lists of var-name sequences, one per sample of that species. A mouse
    symbol matches exactly (``rebuild_pt_expression``). A human symbol matches case-insensitively
    (``cross_species.map_human_to_mouse_space``) and must resolve to the same spelling in every human sample. Rows
    follow the first mouse sample's gene order. Columns: ``mouse_symbol``, ``human_symbol`` (as in the map),
    ``human_source`` (the spelling in the human data).
    """
    pairs = validate_ortholog_map(orthologs)
    if not mouse_genes or not human_genes:
        raise ValueError('Need at least one mouse and one human gene list.')
    mouse_sets = [set(map(str, genes)) for genes in mouse_genes]
    lookups = []
    for genes in human_genes:
        upper = pd.Series([str(g) for g in genes])
        clash = upper[upper.str.upper().duplicated(keep=False)]
        lookups.append({g.upper(): g for g in upper if g not in set(clash)})
        mapped_clash = set(clash.str.upper()) & set(pairs.human_symbol.str.upper())
        if mapped_clash:
            raise ValueError(f'Human genes differing only in case match an ortholog: {sorted(mapped_clash)[:5]}')
    keep = pairs.mouse_symbol.map(lambda m: all(m in s for s in mouse_sets)).to_numpy(bool)
    source = pairs.human_symbol.str.upper().map(lambda h: {lk.get(h) for lk in lookups})
    keep &= source.map(lambda v: None not in v).to_numpy(bool)
    if (source[keep].map(len) != 1).any():
        raise ValueError('A human gene is spelled differently across human samples.')
    table = pairs[keep].assign(human_source=source[keep].map(lambda v: next(iter(v))).to_numpy())
    order = pd.Index(list(map(str, mouse_genes[0])))
    table = table.iloc[np.argsort(order.get_indexer(table.mouse_symbol), kind='stable')].reset_index(drop=True)
    if table.empty:
        raise ValueError('No shared one-to-one ortholog.')
    return table


def to_shared_genes(adata, table: pd.DataFrame, species: str):
    """One sample's units restricted to the shared genes, renamed to mouse symbols (``var['source_symbol']`` keeps
    the native name). X, obs, obsm and uns are carried over unchanged."""
    if species not in ('mouse', 'human'):
        raise ValueError("species must be 'mouse' or 'human'")
    source = table['mouse_symbol' if species == 'mouse' else 'human_source'].astype(str)
    columns = pd.Index(adata.var_names.astype(str)).get_indexer(source)
    if (columns < 0).any():
        raise ValueError(f'{int((columns < 0).sum())} shared genes are missing from this sample.')
    out = adata[:, columns].copy()
    out.var['source_symbol'] = source.to_numpy()
    out.var_names = pd.Index(table.mouse_symbol.astype(str).to_numpy())
    return out


def h5_gene_names(h5_path) -> pd.Index:
    """Gene names of a 10x HDF5 matrix after ``var_names_make_unique``, as ``io.visium_hd.load_bin_matrix`` names
    them (the stage-01 var names)."""
    import h5py
    from anndata.utils import make_index_unique

    with h5py.File(h5_path, 'r') as handle:
        names = pd.Index([n.decode() for n in handle['matrix/features/name'][:]])
    return make_index_unique(names)


def write_shared_bin_dir(source_dir, target_dir, source_names, target_names) -> Path:
    """Write a bin directory whose matrix holds only ``source_names`` (stage-01 spellings), renamed ``target_names``.

    The barcodes, their order and every count of the kept genes are unchanged; feature ids, types and genomes are
    copied. ``spatial/scalefactors_json.json`` and the tissue positions are symbolic links to the source files.
    """
    import h5py
    from scipy import sparse

    source_dir, target_dir = Path(source_dir), Path(target_dir)
    names = h5_gene_names(source_dir / MATRIX_FILE)
    rows = names.get_indexer(pd.Index(list(map(str, source_names))))
    if (rows < 0).any() or len(rows) != len(target_names) or pd.Index(target_names).duplicated().any():
        raise ValueError('source_names must be genes of the matrix and target_names unique, one per source name.')
    (target_dir / 'spatial').mkdir(parents=True, exist_ok=True)
    with h5py.File(source_dir / MATRIX_FILE, 'r') as src:
        m = src['matrix']
        n_genes, n_barcodes = (int(v) for v in m['shape'][:])
        # 10x stores genes x barcodes as CSC; the same arrays read as CSR are barcodes x genes.
        counts = sparse.csr_matrix((m['data'][:], m['indices'][:], m['indptr'][:]), shape=(n_barcodes, n_genes))
        kept = counts[:, rows].tocsr()
        kept.sort_indices()
        features = {key: m['features'][key][:][rows] for key in ('id', 'feature_type', 'genome')}
        tmp = target_dir / (MATRIX_FILE + '.tmp')
        with h5py.File(tmp, 'w') as dst:
            for key, value in src.attrs.items():
                dst.attrs[key] = value
            g = dst.create_group('matrix')
            g.create_dataset('barcodes', data=m['barcodes'][:], compression='gzip')
            g.create_dataset('data', data=kept.data.astype(m['data'].dtype), compression='gzip')   # CSC of genes x barcodes
            g.create_dataset('indices', data=kept.indices.astype(m['indices'].dtype), compression='gzip')
            g.create_dataset('indptr', data=kept.indptr.astype(m['indptr'].dtype), compression='gzip')
            g.create_dataset('shape', data=np.array([len(rows), n_barcodes], dtype=m['shape'].dtype))
            f = g.create_group('features')
            f.create_dataset('_all_tag_keys', data=np.array([b'genome']))
            f.create_dataset('name', data=np.array([str(n).encode() for n in target_names]))
            for key, value in features.items():
                f.create_dataset(key, data=value)
    tmp.replace(target_dir / MATRIX_FILE)
    for name in ('scalefactors_json.json', 'tissue_positions.parquet', 'tissue_positions.csv'):
        source = (source_dir / 'spatial' / name).resolve()
        link = target_dir / 'spatial' / name
        if source.exists() and not link.exists():
            link.symlink_to(source)
    return target_dir


# ---------------------------------------------------------------------------------------------------
# old structure ids
# ---------------------------------------------------------------------------------------------------
def match_polygons(old, new, *, min_iou: float = .5) -> pd.DataFrame:
    """Pair old and new polygons by intersection over union.

    A pair is kept when IoU >= ``min_iou`` and each polygon is the other's best match. With IoU above 0.5 a
    polygon has at most one such partner among polygons that do not overlap each other, so the mutual-best rule
    only matters for overlapping segmentations. Returns ``old``, ``new`` (positions in the inputs), ``iou`` and
    ``centroid_distance``.
    """
    import shapely
    from shapely.strtree import STRtree

    old, new = np.asarray(old, dtype=object), np.asarray(new, dtype=object)
    empty = pd.DataFrame({'old': pd.Series(dtype=int), 'new': pd.Series(dtype=int),
                          'iou': pd.Series(dtype=float), 'centroid_distance': pd.Series(dtype=float)})
    if not len(old) or not len(new):
        return empty
    left, right = STRtree(new).query(old, predicate='intersects')
    if not left.size:
        return empty
    a, b = old[left], new[right]
    inter = shapely.area(shapely.intersection(a, b))
    iou = inter / (shapely.area(a) + shapely.area(b) - inter)
    pairs = pd.DataFrame({'old': left, 'new': right, 'iou': iou})
    best_old = pairs.loc[pairs.groupby('old').iou.idxmax()]
    best_new = pairs.loc[pairs.groupby('new').iou.idxmax()]
    kept = best_old.merge(best_new, on=['old', 'new', 'iou'])
    kept = kept[kept.iou >= min_iou].sort_values('old').reset_index(drop=True)
    kept['centroid_distance'] = shapely.distance(shapely.centroid(old[kept.old.to_numpy()]),
                                                 shapely.centroid(new[kept.new.to_numpy()]))
    return kept


def geojson_schema(path) -> dict:
    """What the package reads from a segmentation GeoJSON: feature count, classes, and how many features carry a
    QuPath ``classification.name``, a ``Score`` measurement and ``properties.id``."""
    payload = json.loads(Path(path).read_text())
    features = payload.get('features', []) if isinstance(payload, dict) else list(payload)
    props = [f.get('properties') or {} for f in features]
    classes = [(p.get('classification') or {}).get('name') if isinstance(p.get('classification'), dict)
               else p.get('classification') for p in props]
    score = [any(m.get('name') == 'Score' for m in (p.get('measurements') or [])) for p in props]
    return {'n_features': len(features), 'n_with_class': int(sum(c is not None for c in classes)),
            'classes': dict(pd.Series([c for c in classes if c is not None]).value_counts()),
            'n_with_score': int(sum(score)), 'n_with_id': int(sum('id' in p for p in props))}


# ---------------------------------------------------------------------------------------------------
# descriptive checks
# ---------------------------------------------------------------------------------------------------
def adjacent_aucs(position, segment, order=('PT-S1', 'PT-S2', 'PT-S3')) -> dict:
    """P(a unit of the later segment lies further along than one of the earlier), ties counted half, per adjacent
    pair of ``order``. NaN when a segment is absent."""
    from scipy.stats import mannwhitneyu

    position, segment = np.asarray(position, float), np.asarray(segment).astype(str)
    out = {}
    for early, late in zip(order[:-1], order[1:]):
        x, y = position[segment == late], position[segment == early]
        out[f'{early}->{late}'] = (float(mannwhitneyu(x, y).statistic / (len(x) * len(y)))
                                   if len(x) and len(y) else np.nan)
    return out


def pt_boundary(position, is_pt) -> dict:
    """The cut c maximising the balanced accuracy of "PT iff position <= c" on a curve that starts in the PT.

    Returns ``cut``, ``balanced_accuracy``, ``pt_beyond`` (PT units with position > c) and ``other_before``
    (non-PT units with position <= c).
    """
    position, is_pt = np.asarray(position, float), np.asarray(is_pt, bool)
    if not is_pt.any() or is_pt.all():
        raise ValueError('Need PT and non-PT units.')
    order = np.argsort(position, kind='stable')
    pos, pt = position[order], is_pt[order]
    # candidate cuts after each distinct position
    tp = np.cumsum(pt)
    fp = np.cumsum(~pt)
    last = np.r_[pos[1:] != pos[:-1], True]
    sens = tp[last] / pt.sum()
    spec = 1 - fp[last] / (~pt).sum()
    best = int(np.argmax(sens + spec))
    cut = float(pos[last][best])
    return {'cut': cut, 'balanced_accuracy': float((sens[best] + spec[best]) / 2),
            'pt_beyond': int((is_pt & (position > cut)).sum()), 'other_before': int((~is_pt & (position <= cut)).sum())}


def uns_safe(value):
    """A copy of nested config values that AnnData can write to HDF5: None -> 'null', tuples -> lists, paths and
    other scalars -> str when they are not numbers or booleans."""
    if isinstance(value, dict):
        return {str(k): uns_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        items = [uns_safe(v) for v in value]
        return items if all(isinstance(v, str) for v in items) or all(
            isinstance(v, (int, float)) and not isinstance(v, bool) for v in items) else json.dumps(items)
    if value is None:
        return 'null'
    if isinstance(value, (bool, int, float, str, np.integer, np.floating, np.bool_)):
        return value
    return str(value)
