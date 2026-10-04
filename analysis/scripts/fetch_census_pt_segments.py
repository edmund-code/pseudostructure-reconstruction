"""Per-donor PT segment pseudobulk from CELLxGENE Census (public data) for the PT gene universe.

Input to notebook 34. Needs the `cellxgene-census` package, which is not part of the project
environment; install it in a separate virtual environment (`python -m venv census_env &&
census_env/bin/pip install cellxgene-census`) and run this script with that interpreter:

    census_env/bin/python analysis/scripts/fetch_census_pt_segments.py --data-root DATA --results-root RESULTS

It reads the notebook 13 ortholog map and the notebook 12 gene universe from the results root and
writes per-donor x sex x segment raw-count sums to DATA/external/census_pt_segments/.
"""
import argparse
import os
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
import cellxgene_census

parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
parser.add_argument('--data-root', type=Path, default=os.environ.get('PSEUDOSPACE_DATA_ROOT'))
parser.add_argument('--results-root', type=Path, default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT'))
args = parser.parse_args()
if args.data_root is None or args.results_root is None:
    parser.error('Pass --data-root and --results-root (or set PSEUDOSPACE_DATA_ROOT / PSEUDOSPACE_RESULTS_ROOT).')
DATA, RESULTS = args.data_root.expanduser(), args.results_root.expanduser()
OUT = DATA / 'external' / 'census_pt_segments'
OUT.mkdir(parents=True, exist_ok=True)
VERSION = '2025-11-08'
orth = pd.read_csv(RESULTS / 'minimal_pt_scfates' / 'ortholog_map_used.csv')
universe = pd.read_csv(RESULTS / 'pt_pathway_remodeling_scfates' / 'gene_statistics.csv', index_col=0).index
orth = orth[orth.mouse_symbol.isin(universe)]
SETS = {
    'mouse': dict(org='mus_musculus', datasets=['25818bf7-e2a7-41ec-8ff2-bc369c0ff4f5'],
                  genes=sorted(orth.mouse_symbol), to_mouse=None),
    'human': dict(org='homo_sapiens', datasets=['09b518f9-da64-44cc-aec8-70a89d55611f'],
                  genes=sorted(orth.human_symbol), to_mouse=dict(zip(orth.human_symbol, orth.mouse_symbol))),
}
with cellxgene_census.open_soma(census_version=VERSION) as census:
    for label, cfg in SETS.items():
        ds = ' or '.join(f"dataset_id == '{d}'" for d in cfg['datasets'])
        obs_filter = (f"({ds}) and is_primary_data == True and cell_type in ['epithelial cell of proximal tubule segment 1', "
                      "'epithelial cell of proximal tubule segment 2', 'epithelial cell of proximal tubule segment 3', "
                      "'kidney proximal convoluted tubule epithelial cell']")
        adata = cellxgene_census.get_anndata(census, organism=cfg['org'], obs_value_filter=obs_filter,
                                             var_value_filter=f"feature_name in {cfg['genes']}", X_name='raw',
                                             obs_column_names=['dataset_id', 'donor_id', 'sex', 'cell_type',
                                                               'development_stage', 'raw_sum'])
        print(label, adata.shape)
        obs = adata.obs.copy()
        obs['group'] = obs.donor_id.astype(str) + '|' + obs.sex.astype(str) + '|' + obs.cell_type.astype(str) + '|' + obs.development_stage.astype(str)
        codes, groups = pd.factorize(obs.group)
        design = sp.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))), shape=(len(groups), len(codes)))
        counts = design @ sp.csr_matrix(adata.X)
        totals = np.asarray(design @ obs.raw_sum.to_numpy(float)).ravel()
        cells = np.asarray(design.sum(axis=1)).ravel()
        frame = pd.DataFrame(counts.toarray(), index=groups, columns=adata.var.feature_name.astype(str))
        if cfg['to_mouse']:
            frame = frame.rename(columns=cfg['to_mouse'])
        frame = frame.T.groupby(level=0).sum().T  # duplicated symbols collapse
        meta = pd.DataFrame([g.split('|') for g in groups], index=groups,
                            columns=['donor_id', 'sex', 'cell_type', 'development_stage'])
        meta['n_cells'] = cells
        meta['total_counts_all_genes'] = totals
        frame.to_csv(OUT / f'{label}_pt_segment_pseudobulk_counts.csv.gz')
        meta.to_csv(OUT / f'{label}_pt_segment_pseudobulk_meta.csv')
        print(meta.groupby(['sex', 'cell_type']).n_cells.agg(['size', 'sum']).to_string())
(OUT / 'SOURCE.txt').write_text(
    f'CELLxGENE Census {VERSION}; public data.\n'
    'mouse: dataset 25818bf7-e2a7-41ec-8ff2-bc369c0ff4f5 (Multi-omic and spatial analysis of mouse kidneys highlights sex-specific ...)\n'
    'human: dataset 09b518f9-da64-44cc-aec8-70a89d55611f (Multimodal benchmarking dataset for renal cortex characterization)\n'
    'Counts are raw sums per donor x sex x cell type x stage over the PT ortholog universe; total_counts_all_genes is the sum of raw_sum over cells.\n')
