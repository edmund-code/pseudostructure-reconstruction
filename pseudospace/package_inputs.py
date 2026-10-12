"""Notebook 69's PT pathway inputs, built from the package pseudospace (notebook 71).

The ``pseudospace_reconstruction`` package (stages 01-04) is the frozen pseudospace for both datasets. Its export,
one ``pt_units.h5ad`` per dataset, holds the PT units after the package's QC with raw integer counts in X (mouse
symbols; human vs mouse keeps the shared one-to-one orthologs and a ``human_symbol`` column) and obs ``sample``,
``group``, ``segment_class`` (the package's PT-S1/S2/S3) and ``pseudospace`` (stage 04, increasing S1 -> S3).

Two builders return the objects notebook 69 analyses, with notebook 69's rules unchanged:
- ``hvm_inputs``: human vs mouse, the dict of ``coordinate_inputs.pt_inputs`` (notebook 37's construction);
- ``aki_inputs``: AKI vs control, notebook 69's cell (notebook 45's construction): the ortholog panel, the library
  over the measured panel genes, and the universe restricted to the human vs mouse universe.

Shared rules: each specimen's [1st, 99th] percentile common support; the library and the normalization on the
measured ortholog panel; ``exclude_genes`` (notebook 69: the ``mt-`` genes) leave the library and the universe, and
the lognorm matrix keeps the full library; the 2% equal-specimen detection rule; 10-300 tested members per pathway
(Reactome 2022, Hallmark 2020, KEGG 2019 Mouse, mapped through the ortholog map); the covariate strata.

Position scale. The NB gene model, its denominator, the grid, the zone cuts and the drift windows are invariant to an
affine change of position, but the clustering method's ``CompareParams`` windows (0.02) and the node spacing (1/29)
are in axis units. ``rescale='support'`` (the default) min-max scales the units inside the common support to [0, 1]
(notebook 69's AKI rule), so those parameters keep their meaning as fractions of the trimmed PT whatever span of the
nephron curve the PT occupies. ``rescale='unit'`` is ``pt_inputs``' rule (min-max over all units, only when one lies
outside [0, 1]); it exists to check the human vs mouse builder against ``pt_inputs`` on the old coordinate.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from .coordinate_inputs import LIBRARIES, common_support_mask
from .pathway_remodeling import gene_covariates, matching_strata
from .pathways import build_pathway_membership

PT_SEGMENTS = ('PT-S1', 'PT-S2', 'PT-S3')
COHORTS = {   # dataset -> {group label: its specimens}, reference group first
    'hvm': {'Mouse': ['Ctrl1A2', 'Ctrl1A4'], 'Human': ['HUK1_COR1', 'HUK1_MED1']},
    'aki': {'Control': ['Ctrl1A2', 'Ctrl1A4'], 'IR': ['IR2A2', 'IR2A4']},
}
TARGET_SUM = 1e4


def read_pt_units(units, dataset):
    """Load (or take) one ``pt_units.h5ad`` and check its contract; returns the AnnData with CSR float counts in X."""
    import anndata as ad

    if dataset not in COHORTS:
        raise ValueError(f'dataset must be one of {sorted(COHORTS)}')
    adata = ad.read_h5ad(units) if isinstance(units, (str, Path)) else units.copy()
    missing = {'sample', 'group', 'segment_class', 'pseudospace'} - set(adata.obs)
    if missing:
        raise ValueError(f'pt_units obs lacks {sorted(missing)}')
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError('pt_units needs unique unit ids and gene symbols.')
    groups = adata.obs.groupby(adata.obs['sample'].astype(str), observed=True)['group'].agg(
        lambda v: sorted(v.astype(str).unique()))
    expected = {s: [g] for g, names in COHORTS[dataset].items() for s in names}
    if groups.to_dict() != expected:
        raise ValueError(f'Unexpected cohort for {dataset}: {groups.to_dict()}')
    if not adata.obs.segment_class.astype(str).isin(PT_SEGMENTS).all():
        raise ValueError(f'Every unit must be labelled one of {PT_SEGMENTS}.')
    if not np.isfinite(adata.obs.pseudospace.to_numpy(float)).all():
        raise ValueError('The pseudospace must be finite.')
    x = sparse.csr_matrix(adata.X, dtype=float)
    if (x.data < 0).any() or not np.isfinite(x.data).all() or not np.allclose(x.data, np.round(x.data), rtol=0, atol=1e-8):
        raise ValueError('X must hold nonnegative raw integer counts.')
    adata.X = x
    return adata


def support_positions(raw, specimen, segment, *, rescale='support'):
    """(mask, position): the units inside every specimen's [1st, 99th] percentile range and their analysed position.

    ``'support'`` min-max scales the kept units to [0, 1]; ``'unit'`` is ``pt_inputs``' rule. Both trim with
    ``common_support_mask``, which an affine map does not change, and require the S1 median below the S3 median.
    """
    raw = np.asarray(raw, float)
    segment = pd.Series(np.asarray(segment).astype(str))
    medians = pd.Series(raw).groupby(segment.to_numpy()).median()
    if not medians.get('PT-S1', np.inf) < medians.get('PT-S3', -np.inf):
        raise ValueError('The coordinate must run from S1 to S3.')
    if rescale == 'unit':
        values = (raw - raw.min()) / np.ptp(raw) if raw.min() < 0 or raw.max() > 1 else raw
        mask = common_support_mask(values, specimen)
        return mask, values[mask]
    if rescale != 'support':
        raise ValueError("rescale must be 'support' or 'unit'.")
    mask = common_support_mask(raw, specimen)
    kept = raw[mask]
    if not np.ptp(kept) > 0:
        raise ValueError('The common support is a single point.')
    return mask, (kept - kept.min()) / np.ptp(kept)


def _lognorm(counts, library):
    if (library <= 0).any():
        raise ValueError('A retained PT unit has zero counts in the measured panel.')
    y = sparse.csr_matrix(counts.multiply((TARGET_SUM / library)[:, None]))
    y.data = np.log1p(y.data)
    return y


def _pathway_sets(data_root, measured, tested, orthologs, set_sizes):
    """Notebook 37's pathways: library members mapped through the ortholog map onto ``measured``, kept when 10-300
    (``set_sizes``) of them are ``tested``; {library::pathway: tested members}."""
    folder = Path(data_root) / 'mouse_vs_human' / 'pathway_gene_sets'
    coverage = pd.concat([build_pathway_membership(
        json.loads((folder / f'{lib}.json').read_text()), measured, ortholog_map=orthologs, library_name=lib,
        tested=tested, min_genes=1) for lib in LIBRARIES], ignore_index=True)
    coverage['pathway_id'] = coverage.library + '::' + coverage.pathway
    coverage = coverage[coverage.n_tested.between(set_sizes[0], min(set_sizes[1], len(tested) - 1))]
    return {row.pathway_id: list(row.genes_present) for row in coverage.itertuples()}


def _grid(position, specimen, n_points=61):
    """pt_inputs' grid: the overlap of the specimens' full ranges (notebook 37; not the gene model's grid)."""
    names = np.unique(specimen)
    lo = max(position[specimen == n].min() for n in names)
    hi = min(position[specimen == n].max() for n in names)
    return np.linspace(lo, hi, n_points)


def _common(adata, dataset, rescale):
    obs = adata.obs
    specimen_all = obs['sample'].astype(str).to_numpy()
    raw = obs.pseudospace.to_numpy(float)
    mask, position = support_positions(raw, specimen_all, obs.segment_class, rescale=rescale)
    case_label = list(COHORTS[dataset])[1]
    kept = obs[mask].copy()
    kept['position'] = position
    return {'obs': kept, 'position': position, 'position_raw': raw[mask],
            'specimen': specimen_all[mask].astype(str), 'segment': obs.segment_class.astype(str).to_numpy()[mask],
            'case': obs.group.astype(str).eq(case_label).to_numpy()[mask],
            'n_units': {'before support': int(len(obs)), 'after support': int(mask.sum())},
            'support_raw': (float(raw[mask].min()), float(raw[mask].max()))}, mask


def hvm_inputs(units, data_root, orthologs, *, exclude_genes=(), detection=.02, set_sizes=(10, 300), strata_bins=3,
               rescale='support'):
    """Human vs mouse PT inputs from the package export: ``pt_inputs``' dict, with notebook 69's ``exclude_genes``.

    Every var gene of the export is a shared one-to-one ortholog (checked against ``orthologs``), measured in both
    inputs unless a var column ``measured_in_both`` says otherwise. The library is the total over the measured genes
    without ``exclude_genes``; ``Y`` is
    log1p(count / full library x 1e4). The universe is the measured genes, minus ``exclude_genes``, with an
    equal-specimen detection >= ``detection``. Extra keys: case (= human), library_with_excluded, excluded_library,
    excluded_universe, position_raw, support_raw, n_units.
    """
    adata = read_pt_units(units, 'hvm')
    pairs = orthologs.set_index(orthologs.mouse_symbol.astype(str)).human_symbol.astype(str)
    if 'human_symbol' not in adata.var or not adata.var_names.isin(pairs.index).all():
        raise ValueError('Every human vs mouse gene must be a pair of the ortholog map, with its human_symbol.')
    if not (adata.var.human_symbol.astype(str).to_numpy() == pairs.reindex(adata.var_names).to_numpy()).all():
        raise ValueError("The export's human symbols disagree with the ortholog map.")
    out, mask = _common(adata, 'hvm', rescale)
    names = pd.Index(adata.var_names.astype(str))
    # A var column measured_in_both, when the export carries one, marks structural zeros (pt_inputs' rule).
    measured = (adata.var.measured_in_both.astype(bool).to_numpy() if 'measured_in_both' in adata.var
                else np.ones(len(names), bool))
    all_counts = adata.X[mask]
    excluded = names.isin(set(map(str, exclude_genes)))
    library_full = np.asarray(all_counts[:, measured].sum(axis=1)).ravel()
    expression = _lognorm(all_counts, library_full)
    universe = gene_covariates(expression, out['position'], out['case'], out['specimen'], names)
    universe['measured_in_both'] = measured & ~excluded
    universe['eligible'] = universe.measured_in_both & universe.detection.ge(detection)
    eligible = universe.eligible.to_numpy()
    genes = universe.loc[eligible, 'gene'].tolist()
    covariates = universe.set_index('gene').loc[genes]
    out.update({
        'Y': expression[:, eligible], 'counts': all_counts[:, eligible], 'human': out['case'], 'genes': genes,
        'covariates': covariates, 'strata': matching_strata(covariates, bins=strata_bins),
        'gene_sets': _pathway_sets(data_root, list(names[measured & ~excluded]), genes, orthologs, set_sizes),
        'grid': _grid(out['position'], out['specimen']),
        'library': np.asarray(all_counts[:, measured & ~excluded].sum(axis=1)).ravel(),
        'library_with_excluded': library_full, 'excluded_library': list(names[measured & excluded]),
        'excluded_universe': list(names[measured & excluded & universe.detection.ge(detection).to_numpy()])})
    return out


def aki_inputs(units, data_root, orthologs, universe, *, exclude_genes=(), detection=.02, set_sizes=(10, 300),
               strata_bins=3, rescale='support'):
    """AKI vs control PT inputs from the package export: notebook 69's construction (notebook 45's).

    The ortholog panel (var genes that are mouse symbols of ``orthologs``) measured in the kept units (any count) gives
    the library, without ``exclude_genes``, and ``Y`` (full library). The candidate genes are ``universe`` (the human
    vs mouse universe; notebook 69 used notebook 12's) among the measured ones; equal-specimen detection >=
    ``detection`` and not excluded makes the universe. Pathways: notebook 37's on ``universe`` minus
    ``exclude_genes``, then restricted to the AKI universe with the same 10-300 rule. Covariates are notebook 45's
    (equal-specimen mean, detection, share of segments with a detection). Same keys as ``hvm_inputs`` except human.
    """
    adata = read_pt_units(units, 'aki')
    out, mask = _common(adata, 'aki', rescale)
    panel = adata.var_names.isin(set(orthologs.mouse_symbol.astype(str)))
    raw = adata.X[mask][:, panel]
    var = pd.Index(adata.var_names[panel].astype(str))
    measured = sorted(set(var[np.asarray(raw.sum(axis=0)).ravel() > 0]))
    counts_all = sparse.csr_matrix(raw[:, var.get_indexer(measured)], dtype=float)
    measured = pd.Index(measured)
    drop = set(map(str, exclude_genes))
    excluded = measured.isin(drop)
    library_full = np.asarray(counts_all.sum(axis=1)).ravel()
    expression = _lognorm(counts_all, library_full)
    pool = [g for g in map(str, universe) if g in set(measured)]
    x = expression[:, measured.get_indexer(pool)]
    specimen, segment = out['specimen'], out['segment']
    per_donor = {}
    for name in np.unique(specimen):
        rows, seg = x[specimen == name], segment[specimen == name]
        present = np.array([np.asarray((rows[seg == s] > 0).sum(axis=0)).ravel() > 0 for s in np.unique(seg)])
        per_donor[name] = (np.asarray(rows.mean(axis=0)).ravel(), np.asarray((rows > 0).mean(axis=0)).ravel(),
                           present.mean(axis=0))
    covariates = pd.DataFrame({key: np.mean([per_donor[d][k] for d in per_donor], axis=0)
                               for k, key in enumerate(('mean_expression', 'detection', 'coverage'))}, index=pool)
    detected = list(covariates.index[covariates.detection.ge(detection)])
    genes = [g for g in detected if g not in drop]
    sets = _pathway_sets(data_root, [g for g in map(str, universe) if g not in drop],
                         [g for g in map(str, universe) if g not in drop], orthologs, set_sizes)
    keep = set(genes)
    sets = {p: [g for g in v if g in keep] for p, v in sets.items()}
    sets = {p: v for p, v in sets.items() if set_sizes[0] <= len(v) <= min(set_sizes[1], len(keep) - 1)}
    columns = measured.get_indexer(genes)
    covariates = covariates.loc[genes]
    out.update({
        'Y': expression[:, columns], 'counts': counts_all[:, columns], 'genes': genes, 'covariates': covariates,
        'strata': matching_strata(covariates, bins=strata_bins), 'gene_sets': sets,
        'grid': _grid(out['position'], specimen),
        'library': library_full - np.asarray(counts_all[:, excluded].sum(axis=1)).ravel(),
        'library_with_excluded': library_full, 'excluded_library': list(measured[excluded]),
        'excluded_universe': [g for g in detected if g in drop]})
    return out
