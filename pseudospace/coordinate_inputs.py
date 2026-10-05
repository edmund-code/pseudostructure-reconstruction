"""Notebook-37 PT analysis inputs on any supplied coordinate (notebooks 56-60).

The coordinate is a parameter: a CSV with a structure id column and a position column. Everything
else (PT structures, common support, normalization, the 2% equal-specimen detection rule, the
10-300-gene pathway rule and the covariate strata) follows notebook 37 unchanged, so results on
SCF13, DPT13, SCF13-ED and later coordinates are directly comparable.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

EXPECTED = {'mouse': ['Ctrl1A2', 'Ctrl1A4'], 'human': ['HUK1_COR1', 'HUK1_MED1']}
LIBRARIES = ('Reactome_2022', 'MSigDB_Hallmark_2020', 'KEGG_2019_Mouse')


def read_coordinate(path, column='position'):
    """Structure id -> position, from a CSV with a `structure` or `structure_id` column."""
    table = pd.read_csv(path)
    id_column = next((c for c in ('structure', 'structure_id') if c in table), None)
    if id_column is None or column not in table:
        raise ValueError(f'{path} needs a structure or structure_id column and a {column!r} column.')
    if table[id_column].duplicated().any():
        raise ValueError(f'{path} has duplicated structure ids.')
    values = table.set_index(table[id_column].astype(str))[column].astype(float)
    if values.isna().all():
        raise ValueError(f'{path}: column {column!r} is empty.')
    return values.dropna()


def common_support_mask(position, specimen, quantiles=(.01, .99)):
    """Structures inside every specimen's [1st, 99th] percentile interval (notebook 37)."""
    s = pd.Series(np.asarray(position, float))
    groups = s.groupby(np.asarray(specimen).astype(str))
    lo = groups.quantile(quantiles[0]).max()
    hi = groups.quantile(quantiles[1]).min()
    return ((s >= lo) & (s <= hi)).to_numpy()


def pt_inputs(results_root, data_root, coordinate, *, detection=.02, set_sizes=(10, 300), strata_bins=3):
    """Rebuild notebook 37's expression, gene universe, pathways and strata on ``coordinate``.

    ``coordinate`` maps structure ids to positions; it is rescaled to [0, 1] when needed and
    must run from S1 to S3. Returns a dict with obs, Y (sparse lognorm, eligible genes), counts,
    position, human, specimen, segment, genes, gene_sets, strata, covariates, grid and library.
    """
    import anndata as ad
    from scipy import sparse

    from .pathway_inputs import rebuild_pt_expression
    from .pathway_remodeling import gene_covariates, matching_strata
    from .pathways import build_pathway_membership

    results_root, data_root = Path(results_root), Path(data_root)
    upstream13 = results_root / 'minimal_pt_scfates'
    saved = ad.read_h5ad(upstream13 / 'cross_species_pt_scfates.h5ad', backed='r')
    obs = saved.obs.copy()
    saved.file.close()
    obs = obs[obs.broad_tubule_marker_call.astype(str).eq('PT')].copy()
    actual = obs.groupby('comparison_species', observed=True)['sample'].apply(
        lambda v: sorted(v.astype(str).unique())).to_dict()
    if actual != EXPECTED:
        raise ValueError(f'Unexpected cohort: {actual}')
    mapped = pd.Series(coordinate, dtype=float).reindex(obs.index.astype(str))
    obs = obs[mapped.notna().to_numpy()].copy()
    values = mapped.dropna().to_numpy(float)
    if values.min() < 0 or values.max() > 1:
        values = (values - values.min()) / (values.max() - values.min())
    obs['shared_pseudospace'] = values
    medians = obs.groupby('segment_class', observed=True).shared_pseudospace.median()
    if not medians['PT-S1'] < medians['PT-S3']:
        raise ValueError('The coordinate must run from S1 to S3.')
    obs = obs[common_support_mask(obs.shared_pseudospace, obs['sample'])].copy()
    orthologs = pd.read_csv(upstream13 / 'ortholog_map_used.csv')
    source_paths = {name: data_root / 'tubule_by_gene' / f'{name}_tubule_by_gene_caleb.h5ad'
                    for names in EXPECTED.values() for name in names}
    adata = rebuild_pt_expression(obs, source_paths, orthologs, target_sum=1e4)
    position = adata.obs.shared_pseudospace.to_numpy(float)
    human = adata.obs.comparison_species.astype(str).eq('human').to_numpy()
    specimen = adata.obs['sample'].astype(str).to_numpy(dtype=str)
    segment = adata.obs.segment_class.astype(str).to_numpy(dtype=str)
    names = EXPECTED['mouse'] + EXPECTED['human']
    lo = max(position[specimen == n].min() for n in names)
    hi = min(position[specimen == n].max() for n in names)

    expression = sparse.csr_matrix(adata.layers['lognorm'], dtype=float)
    universe = gene_covariates(expression, position, human, specimen, adata.var_names)
    universe['measured_in_both'] = adata.var.measured_in_both_inputs.to_numpy()
    universe['eligible'] = universe.measured_in_both & universe.detection.ge(detection)
    genes = universe.loc[universe.eligible, 'gene'].tolist()
    eligible = universe.eligible.to_numpy()
    all_counts = sparse.csr_matrix(adata.layers['counts'], dtype=float)
    library = np.asarray(all_counts[:, universe.measured_in_both.to_numpy()].sum(axis=1)).ravel()
    covariates = universe.set_index('gene').loc[genes]
    coverage = pd.concat([build_pathway_membership(
        json.loads((data_root / 'mouse_vs_human' / 'pathway_gene_sets' / f'{lib}.json').read_text()),
        universe.loc[universe.measured_in_both, 'gene'], ortholog_map=orthologs, library_name=lib,
        tested=genes, min_genes=1) for lib in LIBRARIES], ignore_index=True)
    coverage['pathway_id'] = coverage.library + '::' + coverage.pathway
    coverage = coverage[coverage.n_tested.between(set_sizes[0], min(set_sizes[1], len(genes) - 1))]
    return {'obs': adata.obs, 'Y': expression[:, eligible], 'counts': all_counts[:, eligible],
            'position': position, 'human': human, 'specimen': specimen, 'segment': segment,
            'genes': genes, 'covariates': covariates, 'strata': matching_strata(covariates, bins=strata_bins),
            'gene_sets': {row.pathway_id: list(row.genes_present) for row in coverage.itertuples()},
            'grid': np.linspace(lo, hi, 61), 'library': library}
