"""Inputs of the label-based species comparison (notebook 73), read from the package pseudospace.

Notebooks 42, 43, 46, 47, 48 B, 49 and 64 compared human and mouse PT by S1/S2/S3 label on notebook 13's structures.
Notebook 73 reruns them on the package's units, counts and automatic labels. ``package_inputs.read_pt_units`` checks
the export; this module adds what the label-based analyses read beyond it:

- ``label_units``: the PT units with notebook 43's objects: every unit (no coordinate trim), the gene space, the
  library over it, and log1p CP10k. ``exclude_genes`` (the mt-encoded genes in the primary arm) leave both the gene
  space and the library, which is notebook 69's rule;
- ``all_units``: every QC unit, with its label, family, polygon index, centroid and, for given genes, raw counts. These
  give 42's glomeruli and spillover and 43's depth anchors. They come from ``all_units_positions.csv`` joined to the
  run's ``02_qc/units_qc.h5ad``;
- ``cortical_depth``: 43's depth, the mean centroid distance to the k nearest anchors of the same specimen;
- ``native_sources``: the stage-01 unit matrices (native symbols, every gene) with an obs frame in the format
  ``pathway_inputs.rebuild_pt_expression`` expects, for 64's alternative ortholog pairings;
- ``harmony_features``, ``label_panel_genes`` and ``check_fingerprint`` (the coordinator's label sign-off).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from .coordinate_evaluation import log_normalize
from .package_inputs import read_pt_units

SPECIES = {'Mouse': 'mouse', 'Human': 'human', 'Control': 'mouse', 'IR': 'mouse'}


def _export(results_root, dataset):
    return Path(results_root) / 'package_pseudospace' / dataset


def _run(results_root, dataset):
    return Path(results_root) / 'package_pseudospace' / 'runs' / dataset


def label_units(results_root, dataset='hvm', *, exclude_genes=(), gene_order=None):
    """The PT units of one export for label-based analyses: notebook 43's objects, with no coordinate filter.

    Returns a dict:
    - obs: sample, group, species, segment (S1/S2/S3), feature_index, centroids, n_spots, total_counts,
      old_structure_id;
    - gene_space: the measured genes (all of them unless a var ``measured_in_both`` says otherwise) minus
      ``exclude_genes``;
    - counts (CSR), library (the total over the gene space) and lognorm (log1p CP10k on that library);
    - excluded (the excluded genes present in the export), provenance and fingerprint (pass-1 clustering).

    ``gene_order`` (e.g. the ortholog map's mouse symbols, the order ``rebuild_pt_expression`` gives) reorders the gene
    space; genes it does not list keep the export's order after it. Seeded matched draws (notebook 48 B) depend on the
    gene order, so the old notebooks' order is kept.
    """
    units = read_pt_units(_export(results_root, dataset) / 'pt_units.h5ad', dataset)
    obs = units.obs
    frame = pd.DataFrame({'sample': obs['sample'].astype(str), 'group': obs.group.astype(str)}, index=obs.index.astype(str))
    frame['species'] = frame.group.map(SPECIES)
    frame['segment'] = obs.segment_class.astype(str).str.replace('PT-', '', regex=False).to_numpy()
    for column in ('feature_index', 'x_centroid', 'y_centroid', 'n_spots', 'total_counts', 'old_structure_id'):
        if column in obs:
            frame[column] = obs[column].to_numpy()
    names = pd.Index(units.var_names.astype(str))
    measured = (units.var.measured_in_both.astype(bool).to_numpy() if 'measured_in_both' in units.var
                else np.ones(len(names), bool))
    drop = names.isin(set(map(str, exclude_genes)))
    keep = np.flatnonzero(measured & ~drop)
    if gene_order is not None:
        rank = pd.Series(np.arange(len(gene_order)), index=pd.Index(map(str, gene_order))).groupby(level=0).first()
        position = rank.reindex(names[keep]).to_numpy(float)
        keep = keep[np.lexsort((np.arange(len(keep)), np.where(np.isnan(position), np.inf, position)))]
    counts = sparse.csr_matrix(units.X[:, keep], dtype=float)
    library = np.asarray(counts.sum(axis=1)).ravel()
    provenance = dict(units.uns.get('provenance', {}))
    return {'obs': frame, 'gene_space': names[keep], 'counts': counts, 'library': library,
            'lognorm': log_normalize(counts, library), 'excluded': list(names[measured & drop]),
            'provenance': provenance, 'fingerprint': dict(provenance.get('fingerprints', {}).get('pass1', {}))}


def all_units(results_root, dataset='hvm', *, genes=(), exclude_genes=()):
    """Every QC-retained unit: label and family (``all_units_positions.csv``) joined to ``02_qc/units_qc.h5ad``.

    Returns (obs, counts, library): obs with sample, label, family, feature_index, centroids, n_spots, total_counts
    and seg_class; raw counts (CSR) for ``genes`` that the run measured, as a (matrix, gene index) pair; and each
    unit's library, the total over every QC gene minus ``exclude_genes``.
    """
    import anndata as ad

    positions = pd.read_csv(_export(results_root, dataset) / 'all_units_positions.csv', index_col=0)
    positions.index = positions.index.astype(str)
    qc = ad.read_h5ad(_run(results_root, dataset) / '02_qc' / 'units_qc.h5ad')
    names = pd.Index(qc.obs_names.astype(str))
    missing = positions.index.difference(names)
    if len(missing):
        raise ValueError(f'{len(missing)} units of all_units_positions.csv are missing from units_qc.h5ad.')
    rows = names.get_indexer(positions.index)
    columns = ['feature_index', 'x_centroid', 'y_centroid', 'n_spots', 'total_counts', 'seg_class']
    obs = positions[['sample', 'label', 'family']].astype(str).copy()
    for column in columns:
        obs[column] = qc.obs[column].to_numpy()[rows]
    if not (obs['sample'].to_numpy() == qc.obs['sample'].astype(str).to_numpy()[rows]).all():
        raise ValueError('all_units_positions.csv and units_qc.h5ad disagree on the samples.')
    raw = sparse.csr_matrix(qc.layers['counts'], dtype=float)[rows]
    var = pd.Index(qc.var_names.astype(str))
    kept = ~var.isin(set(map(str, exclude_genes)))
    library = np.asarray(raw[:, kept].sum(axis=1)).ravel()
    present = pd.Index([g for g in map(str, genes) if g in set(var)])
    return obs, (raw[:, var.get_indexer(present)], present), library


def cortical_depth(units, anchors, samples, *, k=3):
    """Notebook 43's depth: mean distance from each unit centroid to the k nearest anchor centroids of its sample.

    ``units`` and ``anchors`` have sample, x_centroid and y_centroid. Units of other samples get NaN.
    """
    from scipy.spatial import cKDTree

    values = np.full(len(units), np.nan)
    unit_sample = units['sample'].astype(str).to_numpy()
    for name in samples:
        own = anchors[anchors['sample'].astype(str).eq(name)]
        if len(own) < k:
            raise ValueError(f'{name}: fewer than {k} anchors.')
        mask = unit_sample == name
        distance, _ = cKDTree(own[['x_centroid', 'y_centroid']].to_numpy(float)).query(
            units.loc[mask, ['x_centroid', 'y_centroid']].to_numpy(float), k=k)
        values[mask] = np.asarray(distance).reshape(mask.sum(), -1).mean(axis=1)
    return values


def native_sources(results_root, dataset, obs):
    """(pt_obs, sources) for ``rebuild_pt_expression``: the stage-01 matrix of each sample (native symbols, all genes)
    and the units of ``obs`` (ids ``<sample>_<stage-01 id>``) with a placeholder coordinate."""
    sources = {s: _run(results_root, dataset) / '01_units' / f'{s}_units.h5ad' for s in obs['sample'].unique()}
    pt_obs = pd.DataFrame({'sample': obs['sample'].astype(str), 'comparison_species': obs.species.astype(str),
                           'shared_pseudospace': .5, 'segment_class': 'PT-' + obs.segment.astype(str),
                           'broad_tubule_marker_call': 'PT'}, index=obs.index)
    return pt_obs, sources


def native_panel(path):
    """Gene ids (symbol -> Ensembl id) of one stage-01 matrix."""
    import anndata as ad

    native = ad.read_h5ad(path, backed='r')
    ids = native.var['gene_ids'].astype(str).copy()
    native.file.close()
    return ids


def harmony_features(results_root, dataset, which='pass1'):
    """The genes of one pass's integration: ``pass1_harmony_features.csv``, or for pass 2 (which has no such table)
    the var flag ``highly_variable_for_harmony`` of ``pass2.h5ad``."""
    import anndata as ad

    folder = _run(results_root, dataset) / '03_integration'
    if which == 'pass1':
        return set(pd.read_csv(folder / 'pass1_harmony_features.csv', index_col=0).index.astype(str))
    if which != 'pass2':
        raise ValueError("which must be 'pass1' or 'pass2'.")
    stored = ad.read_h5ad(folder / 'pass2.h5ad', backed='r')
    genes = set(stored.var_names[stored.var['highly_variable_for_harmony'].to_numpy(bool)].astype(str))
    stored.file.close()
    return genes


def label_panel_genes(provenance, labels=('PT-S1', 'PT-S2', 'PT-S3')):
    """The marker-panel genes that named the given labels, from the export's provenance."""
    panels = provenance['config']['marker_panels_yaml']['label_panels']
    return {str(g) for label in labels for g in panels[label]}


def check_fingerprint(results_root, dataset, provenance, confirmed):
    """The coordinator's label sign-off: n_cells, n_clusters, membership_sha1 and the pass-1 label map must equal
    ``confirmed``, both in the run's ``clustering_fingerprints.json`` and (for the first three) in the export."""
    if not confirmed:
        raise ValueError('No confirmed label fingerprint: the coordinator has not signed off the labels.')
    record = json.loads((_run(results_root, dataset) / '03_integration' / 'clustering_fingerprints.json').read_text())
    observed = {'n_cells': record['pass1']['n_cells'], 'n_clusters': record['pass1']['n_clusters'],
                'membership_sha1': record['pass1']['membership_sha1'],
                'labels_pass1': {str(k): str(v) for k, v in record['labels_pass1'].items()}}
    exported = provenance.get('fingerprints', {}).get('pass1', {})
    for key, value in confirmed.items():
        value = {str(k): str(v) for k, v in value.items()} if isinstance(value, dict) else value
        if observed.get(key) != value:
            raise ValueError(f'Label fingerprint {key} differs from the confirmed value: {observed.get(key)} != {value}')
        if key != 'labels_pass1' and exported.get(key) != value:
            raise ValueError(f'The export was built from other labels ({key}: {exported.get(key)}).')
    return observed
