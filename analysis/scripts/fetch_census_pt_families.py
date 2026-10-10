"""Per-donor PT segment pseudobulk from CELLxGENE Census (public data) for mouse paralog-family members.

Input to notebook 68 (family sums of the species-level positive control). Needs the `cellxgene-census`
package, which is not part of the project environment; run it with a separate interpreter, as for
`fetch_census_pt_segments.py`:

    census_env/bin/python analysis/scripts/fetch_census_pt_families.py --data-root DATA --results-root RESULTS

Genes: every mouse-panel gene that notebook 62 places in a non-one-to-one HCOP orthogroup (support >= 3; class
1:many, many:1 or many:many), read from RESULTS/pt_non_one_to_one_orthologs/tables/gene_inventory_mouse.csv and
matched to Census by Ensembl id (feature_id), not by symbol. Cells, dataset, version and grouping are those of
`fetch_census_pt_segments.py` (mouse dataset 25818bf7, primary data, PT segment 1/2/3 and convoluted PT, one group
per donor x sex x cell type x stage), so the rows align with DATA/external/census_pt_segments/mouse_*.
Writes next to those files: raw-count sums per group (columns = our mouse panel symbols), the group metadata, the
gene match table and SOURCE_families.txt. Genes already in mouse_pt_segment_pseudobulk_counts.csv.gz are fetched
again as controls; notebook 68 checks that they are identical.
"""
import argparse
import os
from pathlib import Path

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
OUT = DATA / 'external' / 'census_pt_segments'
OUT.mkdir(parents=True, exist_ok=True)
VERSION = '2025-11-08'
DATASET = '25818bf7-e2a7-41ec-8ff2-bc369c0ff4f5'
CELL_TYPES = ['epithelial cell of proximal tubule segment 1', 'epithelial cell of proximal tubule segment 2',
              'epithelial cell of proximal tubule segment 3', 'kidney proximal convoluted tubule epithelial cell']
FAMILY_CLASSES = ('1:many', 'many:1', 'many:many')

inventory = pd.read_csv(RESULTS / 'pt_non_one_to_one_orthologs' / 'tables' / 'gene_inventory_mouse.csv',
                        usecols=['gene', 'gene_id', 'component', 'ortholog_class'])
wanted = inventory[inventory.ortholog_class.isin(FAMILY_CLASSES)].reset_index(drop=True)
if wanted.gene_id.isna().any() or wanted.gene_id.duplicated().any() or wanted.gene.duplicated().any():
    raise ValueError('Expected one Ensembl id per family member.')
print(f'{len(wanted):,} mouse family members in {wanted.component.nunique():,} orthogroups')

with cellxgene_census.open_soma(census_version=VERSION) as census:
    var = census['census_data']['mus_musculus'].ms['RNA'].var.read(
        column_names=['feature_id', 'feature_name']).concat().to_pandas()
    present = sorted(set(wanted.gene_id) & set(var.feature_id))
    adata = cellxgene_census.get_anndata(
        census, organism='mus_musculus',
        obs_value_filter=f"dataset_id == '{DATASET}' and is_primary_data == True and cell_type in {CELL_TYPES}",
        var_value_filter=f'feature_id in {present}', X_name='raw',
        obs_column_names=['dataset_id', 'donor_id', 'sex', 'cell_type', 'development_stage', 'raw_sum'])
print('Census cells x genes:', adata.shape)

obs = adata.obs.copy()
obs['group'] = (obs.donor_id.astype(str) + '|' + obs.sex.astype(str) + '|' + obs.cell_type.astype(str) + '|'
                + obs.development_stage.astype(str))
codes, groups = pd.factorize(obs.group)
design = sp.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))), shape=(len(groups), len(codes)))
sums = design @ sp.csr_matrix(adata.X)
symbol_of = dict(zip(wanted.gene_id, wanted.gene))
columns = [symbol_of[i] for i in adata.var.feature_id.astype(str)]
frame = pd.DataFrame(sums.toarray(), index=groups, columns=columns)
frame = frame.reindex(columns=[g for g in wanted.gene if g in set(columns)])
meta = pd.DataFrame([g.split('|') for g in groups], index=groups,
                    columns=['donor_id', 'sex', 'cell_type', 'development_stage'])
meta['n_cells'] = np.asarray(design.sum(axis=1)).ravel()
meta['total_counts_all_genes'] = np.asarray(design @ obs.raw_sum.to_numpy(float)).ravel()

census_name = dict(zip(var.feature_id, var.feature_name))
genes = wanted.assign(found=wanted.gene_id.isin(set(present)),
                      census_feature_name=wanted.gene_id.map(census_name))


def symbol_audit(columns, our_ids, census_var, species):
    """Check a symbol-matched pseudobulk: does each column's symbol name our panel's Ensembl id in Census?

    status: match (one Census feature with that symbol, our id), several_ids (our id among several Census features
    with that symbol, which the symbol fetch summed), other_id (Census gives the symbol to another feature: the column
    measured another gene), absent (no Census feature has the symbol).
    """
    ids = census_var.groupby('feature_name').feature_id.apply(list)
    rows = []
    for column, ours in zip(columns, our_ids):
        found = ids.get(ours[0], []) if ours[0] is not None else []
        status = ('absent' if not found else 'other_id' if ours[1] not in found
                  else 'match' if len(found) == 1 else 'several_ids')
        rows.append({'species': species, 'column': column, 'symbol': ours[0], 'our_gene_id': ours[1],
                     'census_ids_for_symbol': ';'.join(found), 'status': status})
    return pd.DataFrame(rows)


# Audit of the symbol-matched pseudobulks of fetch_census_pt_segments.py (metadata only): our Ensembl ids come from
# the tubule_by_gene matrices (the panels), symbols in the human file are mapped back through the ortholog map.
import anndata as ad   # noqa: E402

with cellxgene_census.open_soma(census_version=VERSION) as census:
    human_var = census['census_data']['homo_sapiens'].ms['RNA'].var.read(
        column_names=['feature_id', 'feature_name']).concat().to_pandas()
panel_ids = {}
for species, sample in (('mouse', 'Ctrl1A2'), ('human', 'HUK1_COR1')):
    handle = ad.read_h5ad(DATA / 'tubule_by_gene' / f'{sample}_tubule_by_gene_caleb.h5ad', backed='r')
    panel_ids[species] = handle.var['gene_ids'].astype(str).to_dict()
    handle.file.close()
orth = pd.read_csv(RESULTS / 'minimal_pt_scfates' / 'ortholog_map_used.csv')
to_human = dict(zip(orth.mouse_symbol, orth.human_symbol))
mouse_columns = pd.read_csv(OUT / 'mouse_pt_segment_pseudobulk_counts.csv.gz', index_col=0, nrows=0).columns
human_columns = pd.read_csv(OUT / 'human_pt_segment_pseudobulk_counts.csv.gz', index_col=0, nrows=0).columns
audit = pd.concat([
    symbol_audit(mouse_columns, [(c, panel_ids['mouse'].get(c)) for c in mouse_columns], var, 'mouse'),
    symbol_audit(human_columns, [(to_human.get(c), panel_ids['human'].get(to_human.get(c))) for c in human_columns],
                 human_var, 'human')], ignore_index=True)
audit.to_csv(OUT / 'census_symbol_audit.csv', index=False)
print(audit.groupby(['species', 'status']).size().to_string())
print(audit[audit.status.ne('match')].to_string())
frame.to_csv(OUT / 'mouse_pt_segment_pseudobulk_counts_families.csv.gz')
meta.to_csv(OUT / 'mouse_pt_segment_pseudobulk_meta_families.csv')
genes.to_csv(OUT / 'mouse_family_genes_fetched.csv', index=False)
(OUT / 'SOURCE_families.txt').write_text(
    f'CELLxGENE Census {VERSION}; public data. Mouse dataset {DATASET} only.\n'
    'Written by analysis/scripts/fetch_census_pt_families.py (notebook 68 input).\n'
    'mouse_pt_segment_pseudobulk_counts_families.csv.gz: raw-count sums per donor x sex x cell type x stage for every\n'
    '  mouse-panel gene in a non-one-to-one HCOP orthogroup of notebook 62 (support >= 3; 1:many, many:1, many:many),\n'
    '  matched by Ensembl id; columns are our mouse panel symbols. Same cells and groups as\n'
    '  mouse_pt_segment_pseudobulk_counts.csv.gz (primary data; PT segment 1/2/3 and convoluted PT).\n'
    'mouse_pt_segment_pseudobulk_meta_families.csv: group metadata of this fetch (n_cells, total_counts_all_genes =\n'
    '  sum of raw_sum over cells); notebook 68 checks it equals mouse_pt_segment_pseudobulk_meta.csv.\n'
    'mouse_family_genes_fetched.csv: requested genes, whether Census has the Ensembl id, and the Census symbol.\n'
    f'{int(genes.found.sum()):,} of {len(genes):,} requested genes found.\n'
    'census_symbol_audit.csv: for every column of the symbol-matched mouse_/human_pt_segment_pseudobulk_counts.csv.gz,\n'
    '  whether Census gives that symbol to our panel Ensembl id (match), to several features including ours (several_ids,\n'
    '  summed by the symbol fetch), to another feature only (other_id: the column measured another gene) or to none.\n'
    f'  Columns by status: {audit.groupby(["species", "status"]).size().to_dict()}.\n')
print(meta.groupby(['sex', 'cell_type']).n_cells.agg(['size', 'sum']).to_string())
print(f'found {int(genes.found.sum()):,} of {len(genes):,}; missing: {sorted(genes.loc[~genes.found, "gene"])}')
