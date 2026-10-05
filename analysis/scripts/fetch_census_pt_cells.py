"""Nucleus-level PT raw counts from CELLxGENE Census (public data) over the accepted ortholog panel.

Input to notebook 45 (many-draw specimen null). Needs the `cellxgene-census` package, which is not
part of the project environment; run it with a separate interpreter, as for
`fetch_census_pt_segments.py`:

    census_env/bin/python analysis/scripts/fetch_census_pt_cells.py --data-root DATA --results-root RESULTS

It reads notebook 13's ortholog map from the results root and writes one AnnData per species with raw
counts (nuclei x ortholog genes, mouse symbols as var names) and donor/sex/segment/stage metadata to
DATA/external/census_pt_cells/.
"""
import argparse
import os
from pathlib import Path

import anndata as ad
import cellxgene_census
import numpy as np
import pandas as pd
import scipy.sparse as sp

parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
parser.add_argument('--data-root', type=Path, default=os.environ.get('PSEUDOSPACE_DATA_ROOT'))
parser.add_argument('--results-root', type=Path, default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT'))
args = parser.parse_args()
if args.data_root is None or args.results_root is None:
    parser.error('Pass --data-root and --results-root (or set PSEUDOSPACE_DATA_ROOT / PSEUDOSPACE_RESULTS_ROOT).')
DATA, RESULTS = args.data_root.expanduser(), args.results_root.expanduser()
OUT = DATA / 'external' / 'census_pt_cells'
OUT.mkdir(parents=True, exist_ok=True)
VERSION = '2025-11-08'
orth = pd.read_csv(RESULTS / 'minimal_pt_scfates' / 'ortholog_map_used.csv')
CELL_TYPES = ['epithelial cell of proximal tubule segment 1', 'epithelial cell of proximal tubule segment 2',
              'epithelial cell of proximal tubule segment 3', 'kidney proximal convoluted tubule epithelial cell']
SETS = {
    'mouse': dict(org='mus_musculus', dataset='25818bf7-e2a7-41ec-8ff2-bc369c0ff4f5',
                  genes=sorted(orth.mouse_symbol), to_mouse=None),
    'human': dict(org='homo_sapiens', dataset='09b518f9-da64-44cc-aec8-70a89d55611f',
                  genes=sorted(orth.human_symbol), to_mouse=dict(zip(orth.human_symbol, orth.mouse_symbol))),
}
with cellxgene_census.open_soma(census_version=VERSION) as census:
    for label, cfg in SETS.items():
        obs_filter = (f"dataset_id == '{cfg['dataset']}' and is_primary_data == True and cell_type in {CELL_TYPES}")
        adata = cellxgene_census.get_anndata(census, organism=cfg['org'], obs_value_filter=obs_filter,
                                             var_value_filter=f"feature_name in {cfg['genes']}", X_name='raw',
                                             obs_column_names=['dataset_id', 'donor_id', 'sex', 'cell_type',
                                                               'development_stage', 'raw_sum'])
        symbols = adata.var.feature_name.astype(str).to_numpy()
        if cfg['to_mouse']:
            symbols = np.array([cfg['to_mouse'][s] for s in symbols])
        # Duplicated symbols (several feature ids for one name) collapse by summation.
        unique, codes = np.unique(symbols, return_inverse=True)
        collapse = sp.csr_matrix((np.ones(len(codes)), (np.arange(len(codes)), codes)), shape=(len(codes), len(unique)))
        counts = (sp.csr_matrix(adata.X, dtype=np.float32) @ collapse).tocsr().astype(np.float32)
        obs = adata.obs[['donor_id', 'sex', 'cell_type', 'development_stage', 'raw_sum']].astype(
            {'donor_id': str, 'sex': str, 'cell_type': str, 'development_stage': str})
        obs.index = [f'{label}_{i}' for i in adata.obs.index.astype(str)]
        out = ad.AnnData(X=counts, obs=obs, var=pd.DataFrame(index=pd.Index(unique, name='mouse_symbol')))
        out.uns['source'] = {'census_version': VERSION, 'dataset_id': cfg['dataset'], 'X': 'raw counts'}
        out.write_h5ad(OUT / f'{label}_pt_cells.h5ad', compression='gzip')
        print(label, out.shape, out.obs.groupby(['sex', 'cell_type']).size().to_string())
(OUT / 'SOURCE.txt').write_text(
    f'CELLxGENE Census {VERSION}; public data.\n'
    'mouse: dataset 25818bf7-e2a7-41ec-8ff2-bc369c0ff4f5 (Multi-omic and spatial analysis of mouse kidneys highlights sex-specific ...)\n'
    'human: dataset 09b518f9-da64-44cc-aec8-70a89d55611f (Multimodal benchmarking dataset for renal cortex characterization)\n'
    'Nucleus-level raw counts of PT cells (segment 1/2/3, convoluted PT) over the accepted ortholog panel, '
    'mouse symbols as var names; obs carries donor_id, sex, cell_type, development_stage and raw_sum.\n'
    'Written by analysis/scripts/fetch_census_pt_cells.py (notebook 45 input).\n')
