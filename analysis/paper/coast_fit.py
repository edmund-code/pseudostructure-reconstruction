"""The paper's two fits built with COAST (github.com/edmundtsou/COAST), the finalized method, in the run layout the
figure scripts read, plus `refit` for the reruns of Fig. 3 (gene halves) and Fig. 6 (injury-blind axis).

  python analysis/paper/coast_fit.py --run aki      # control vs AKI mouse
  python analysis/paper/coast_fit.py --run hvm      # control mouse vs human, on shared 1:1 orthologs

Steps (the settings of the frozen package run, analysis/package_configs/<config>):
  units      stage-01 units of the package (segmentation polygon x gene counts; results/package_pseudospace/runs/
             <run>/01_units); hvm genes restricted to the shared 1:1 orthologs of 00_ortholog_space, in mouse symbols
  QC         coast.pp.filter_units (>= 500 counts, bin coverage >= 0.9, MAD 3 per sample x segmentation class on log
             counts, log counts per bin and the complexity residual) and coast.pp.doublets (tiers 2, 3, 4 discard;
             glomeruli exempt; regions and adjacency of marker_panels.yaml qc_regions)
  pass 1     coast.pp.filter_genes (5 % of units, 20 counts, panel genes kept), normalize, integrate (Harmony on
             sample, variable genes in every condition), Leiden (resolution 2), coast.tl.label_clusters
  pass 2     the nephron tubule labels only (drop_classes, rejected and Unassigned removed), the same steps again
  curve      coast.tl.pseudostructure on the 50 Harmony dimensions, oriented by the axis_anchors (early -> late)

Writes <results-root>/coast_runs/<run>/: 01_units (link to the package's stage 01), 02_qc/{units_qc.h5ad, qc_flags.parquet},
03_integration/{pass1_labelled,pass2}.h5ad and cluster scores, 04_pseudospace/{pseudospace.parquet,
curve_nodes.npz, order_metrics.csv, manifest.json}.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

REPO = Path(__file__).resolve().parents[2]
CONFIGS = {'aki': 'analysis/package_configs/aki', 'hvm': 'analysis/package_configs/human_vs_mouse'}
PACKAGE_RUNS = 'package_pseudospace/runs'
OUT = 'coast_runs'
SEED = 0


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--run', choices=sorted(CONFIGS), required=True)
    p.add_argument('--data-root', default=os.environ.get('PSEUDOSPACE_DATA_ROOT', str(REPO / 'data')))
    p.add_argument('--results-root', default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT', str(REPO / 'results')))
    p.add_argument('--visium-root', default=os.environ.get('PSEUDOSPACE_VISIUM_ROOT'))
    return p.parse_args(argv)


def panels(run: str) -> dict:
    return yaml.safe_load((REPO / CONFIGS[run] / 'marker_panels.yaml').read_text())


def samples(run: str, data_root: Path, visium_root: str | None) -> pd.DataFrame:
    s = pd.read_csv(REPO / CONFIGS[run] / 'samples.csv', comment='#')
    s['bins_dir'] = [Path(b.replace('$PSEUDOSPACE_VISIUM_ROOT', visium_root or '')) if b.startswith('$')
                     else data_root / b for b in s.bins_dir]
    s['segmentation_path'] = [data_root / p for p in s.segmentation_path]
    return s.set_index('sample')


def protected_genes(P: dict) -> list:
    return sorted({g for gs in P['label_panels'].values() for g in gs} | set(P['axis_anchors']['early'])
                  | set(P['axis_anchors']['late']) | set(P.get('extra_marker_genes', []))
                  | set(P.get('extra_protected_genes', [])) | {g for gs in P['qc_regions']['markers'].values() for g in gs})


def exclude_regex(genes=()) -> str:
    """COAST's default mitochondrial / ribosomal exclusion plus a list of gene names withdrawn from the features."""
    import coast

    if not len(genes):
        return coast.pp.DEFAULT_EXCLUDE
    return coast.pp.DEFAULT_EXCLUDE + '|^(?:' + '|'.join(map(re.escape, sorted(genes))) + ')$'


# ------------------------------------------------------------------------------------------------- units and QC
def load_units(run: str, results_root: Path, S: pd.DataFrame):
    import anndata as ad

    stage01 = results_root / PACKAGE_RUNS / run / '01_units'
    binning = pd.read_csv(stage01 / 'binning_summary.csv').set_index('sample')
    shared = (pd.read_csv(results_root / PACKAGE_RUNS / run / '00_ortholog_space' / 'shared_orthologs.csv')
              if run == 'hvm' else None)
    parts = {}
    for sample, row in S.iterrows():
        a = ad.read_h5ad(stage01 / f'{sample}_units.h5ad')
        a.obs_names = sample + '_' + a.obs_names.astype(str)
        a.obs['unit_id'] = a.obs_names.str.removeprefix(sample + '_')
        a.obs['condition'] = row['group']
        a.obs['bin_coverage'] = a.obs['n_spots'] * binning.loc[sample, 'bin_area_px'] / a.obs['polygon_area_px']
        if shared is not None:   # one gene space: the shared 1:1 orthologs, named by their mouse symbol
            human = row['group'] == 'Human'
            a = a[:, (shared.human_source if human else shared.mouse_symbol).to_numpy()].copy()
            a.var_names = shared.mouse_symbol.to_numpy()
        a.var = a.var[['gene_ids']]   # Ensembl ids (fig56_annotate.py); differ by species, so dropped for hvm
        parts[sample] = a
    adata = ad.concat(parts, merge='same')
    adata.obs['sample'] = pd.Categorical(adata.obs['sample'], categories=list(S.index))
    adata.obs['condition'] = adata.obs['condition'].astype('category')
    return adata, stage01


def bins_for_doublets(adata, stage01: Path, S: pd.DataFrame, markers: list, run: str, results_root: Path):
    """Bins x marker genes (raw counts) of the units in `adata`, with obs unit / x / y / total_counts."""
    import anndata as ad
    import scanpy as sc

    parts = []
    for sample, row in S.iterrows():
        h5 = (results_root / PACKAGE_RUNS / run / '00_ortholog_space' / 'bins' / sample if run == 'hvm'
              else row['bins_dir']) / 'filtered_feature_bc_matrix.h5'
        m = sc.read_10x_h5(h5)
        m.var_names_make_unique()
        total = np.asarray(m.X.sum(1)).ravel()
        assign = pd.read_parquet(stage01 / f'{sample}_bin_to_unit.parquet')
        assign['unit'] = sample + '_' + assign['polygon_id'].astype(str)
        assign = assign[assign['unit'].isin(adata.obs_names)]
        idx = m.obs_names.get_indexer(assign['bin_id'])
        assert (idx >= 0).all(), f'{sample}: bins of the assignment missing from the matrix'
        genes = [g for g in markers if g in m.var_names]
        b = ad.AnnData(m.X[idx][:, m.var_names.get_indexer(genes)].tocsr(),
                       obs=pd.DataFrame({'unit': assign['unit'].to_numpy(), 'x': assign['x'].to_numpy(float),
                                         'y': assign['y'].to_numpy(float), 'total_counts': total[idx]},
                                        index=sample + ':' + assign['bin_id'].to_numpy()),
                       var=pd.DataFrame(index=genes))
        parts.append(b)
        del m
    return ad.concat(parts, join='outer', fill_value=0)


def polygons_for(adata, S: pd.DataFrame) -> list:
    from shapely.geometry import shape

    out = pd.Series(index=adata.obs_names, dtype=object)
    for sample, row in S.iterrows():
        features = json.loads(Path(row['segmentation_path']).read_text())['features']
        obs = adata.obs[adata.obs['sample'] == sample]
        out[obs.index] = [shape(features[int(i)]['geometry']) for i in obs['feature_index']]
    return out.tolist()


def quality_control(adata, stage01: Path, S: pd.DataFrame, P: dict, run: str, results_root: Path):
    import coast

    floors = coast.pp.filter_units(adata, min_counts=500, min_coverage=0.9, n_mads=None)
    regions, adjacency = P['qc_regions']['markers'], P['qc_regions']['adjacency']
    markers = sorted({g for gs in regions.values() for g in gs})
    bins = bins_for_doublets(floors, stage01, S, markers, run, results_root)
    coast.pp.doublets(floors, bins, polygons_for(floors, S), regions, adjacency, class_key='seg_class',
                      exempt=('glomerulus',), discard=('2', '3', '4'), seed=SEED)
    mad = coast.pp.filter_units(floors, min_counts=500, min_coverage=0.9, n_mads=3.0,
                                groupby=['sample', 'seg_class'], bins_key='n_spots', min_group_size=50)
    keep = floors.obs_names.isin(mad.obs_names) & ~floors.obs['doublet'].to_numpy(bool)
    flags = floors.obs[['sample', 'seg_class', 'doublet', *[c for c in floors.obs if c.startswith('doublet_tier')]]].copy()
    flags['mad_outlier'] = ~floors.obs_names.isin(mad.obs_names)
    flags['qc_retained'] = keep
    qc = floors[keep].copy()
    counts = {'polygons': adata.obs['sample'].value_counts(sort=False),
              'with bins': adata.obs.loc[adata.obs['n_spots'] > 0, 'sample'].value_counts(sort=False),
              'hard floors': floors.obs['sample'].value_counts(sort=False),
              'MAD retained': mad.obs['sample'].value_counts(sort=False),
              'doublets': floors.obs.loc[floors.obs['doublet'], 'sample'].value_counts(sort=False),
              'after QC': qc.obs['sample'].value_counts(sort=False)}
    coast.pp.normalize(qc)
    qc.layers['lognorm'] = qc.X.copy()
    qc.var['mt'] = qc.var_names.str.match(r'^(?:mt-|MT-)')
    qc.var['ribo'] = qc.var_names.str.match(r'^(?:Rp[ls]\d|RP[LS]\d|Mrp[ls]\d|MRP[LS]\d)')
    qc.var['exclude_from_harmony_hvg'] = qc.var_names.str.contains(coast.pp.DEFAULT_EXCLUDE)
    qc.obs['qc_retained'] = True
    return qc, pd.DataFrame(counts).fillna(0).astype(int), flags


# ------------------------------------------------------------------------------------------- two passes and curve
def _pass(units, P: dict, group_key: str, exclude: str, key: str):
    """filter_genes -> normalize -> integrate -> Leiden on QC units (raw counts in layers['counts'])."""
    import coast
    import scanpy as sc

    a = coast.pp.filter_genes(units, keep=protected_genes(P))
    coast.pp.normalize(a)
    coast.pp.integrate(a, batch_key='sample', group_key=group_key, exclude=exclude, random_state=SEED)
    sc.tl.leiden(a, resolution=2, flavor='igraph', n_iterations=2, random_state=SEED, key_added=key)
    return a


def two_pass(qc, P: dict, exclude_genes=(), pre_umap: bool = False):
    """Both passes on QC units. Returns (pass 1, pass 2 with obs['pseudostructure'], cluster score table)."""
    import coast
    import scanpy as sc

    exclude = exclude_regex(exclude_genes)
    units = qc.copy()
    units.X = units.layers['counts'].copy()
    p1 = _pass(units, P, 'condition', exclude, 'leiden_coarse')
    if pre_umap:   # Fig. 1: the same units before Harmony
        pre = sc.AnnData(obs=p1.obs[[]], obsm={'X_pca': p1.obsm['X_pca']})
        sc.pp.neighbors(pre, use_rep='X_pca', n_neighbors=25, transformer='sklearn', random_state=SEED)
        sc.tl.umap(pre, min_dist=0.3, random_state=SEED)
        p1.obsm['X_umap_pre_harmony'] = pre.obsm['X_umap']
    scores = coast.tl.label_clusters(p1, P['label_panels'], cluster_key='leiden_coarse', key_added='segment_class')
    rejected = P.get('rejected_labels', {})
    family = P['label_to_family']
    label = p1.obs['segment_class'].astype(str).replace(rejected)
    p1.obs['segment_class'] = label.astype('category')
    p1.obs['coarse_class'] = label.map(lambda l: family.get(l, l)).astype('category')
    p1.obs['broad_tubule_marker_call'] = p1.obs['coarse_class']
    drop = set(P['drop_classes']) | set(rejected.values()) | {'Unassigned'}
    tubule = ~p1.obs['coarse_class'].astype(str).isin(drop) & ~label.isin(drop)
    p2 = _pass(units[tubule.to_numpy()].copy(), P, 'condition', exclude, 'leiden_pass2')
    for col in ('segment_class', 'coarse_class', 'broad_tubule_marker_call', 'leiden_coarse'):
        p2.obs[col] = p1.obs.loc[p2.obs_names, col].astype(str).astype('category')
    coast.tl.pseudostructure(p2, early_genes=P['axis_anchors']['early'], late_genes=P['axis_anchors']['late'],
                             seed=SEED)
    return p1, p2, scores


def order_metrics(table: pd.DataFrame, order: list) -> pd.DataFrame:
    """Spearman rho of the coordinate against the label's place in `order`, and inverted adjacent label medians."""
    from scipy.stats import spearmanr

    rows = []
    for name, part in [('all', table)] + [(c, g) for c, g in table.groupby('condition', observed=True)]:
        rank = part['segment_class'].map({l: i for i, l in enumerate(order)})
        ok = rank.notna()
        med = part[ok].groupby('segment_class', observed=True)['pseudospace'].median().reindex(order).dropna()
        rows.append({'coordinate': 'pseudospace', 'units': name, 'n_units': int(len(part)),
                     'segment_order_rho': float(spearmanr(med.index.map({l: i for i, l in enumerate(order)}),
                                                          med.to_numpy()).statistic),
                     'n_inversions': int((np.diff(med.to_numpy()) < 0).sum()),
                     'unit_rank_rho': float(spearmanr(rank[ok], part.loc[ok, 'pseudospace']).statistic)})
    return pd.DataFrame(rows)


def stage04_table(p2) -> pd.DataFrame:
    t = pd.DataFrame({'sample': p2.obs['sample'].astype(str), 'condition': p2.obs['condition'].astype(str),
                      'segment_class': p2.obs['segment_class'].astype(str),
                      'coarse_class': p2.obs['coarse_class'].astype(str),
                      'pseudospace': p2.obs['pseudostructure'].astype(float)}, index=p2.obs_names.copy())
    t.index.name = 'obs_name'
    return t


def write_stage04(p2, P: dict, folder: Path, extra: dict | None = None) -> pd.DataFrame:
    folder.mkdir(parents=True, exist_ok=True)
    table = stage04_table(p2)
    table.to_parquet(folder / 'pseudospace.parquet')
    np.savez(folder / 'curve_nodes.npz', nodes=p2.uns['coast']['nodes'],
             pseudospace=p2.uns['coast']['node_pseudostructure'])
    order_metrics(table, P['segment_order']).to_csv(folder / 'order_metrics.csv', index=False)
    root = table['pseudospace'].idxmin()
    (folder / 'manifest.json').write_text(json.dumps({'stage': '04', 'method': 'COAST', 'summary': {
        'n_units': int(len(table)), 'root_unit': str(root), 'root_rule': 'unit at pseudostructure 0',
        'n_nodes': int(len(p2.uns['coast']['nodes'])), **(extra or {})}}, indent=1))
    return table


def h5_safe(a):
    """Drop what h5ad cannot store (graph params with None values are fine; uns objects from COAST are arrays)."""
    a.uns.pop('coast', None)
    return a


def build(args) -> None:
    import coast

    t0 = time.time()
    run, results_root, data_root = args.run, Path(args.results_root), Path(args.data_root)
    out = results_root / OUT / run
    for stage in ('02_qc', '03_integration', '04_pseudospace'):
        (out / stage).mkdir(parents=True, exist_ok=True)
    link = out / '01_units'
    if not link.exists():
        link.symlink_to(results_root / PACKAGE_RUNS / run / '01_units', target_is_directory=True)
    P, S = panels(run), samples(run, data_root, args.visium_root)
    adata, stage01 = load_units(run, results_root, S)
    qc, counts, flags = quality_control(adata, stage01, S, P, run, results_root)
    flags.to_parquet(out / '02_qc' / 'qc_flags.parquet')   # the hard-floor units: doublet tiers, MAD, retained
    print(counts.to_string())
    qc.write_h5ad(out / '02_qc' / 'units_qc.h5ad')
    p1, p2, scores = two_pass(qc, P, pre_umap=True)
    scores.to_csv(out / '03_integration' / 'pass1_cluster_scores.csv')
    nodes = dict(p2.uns['coast'])
    h5_safe(p1).write_h5ad(out / '03_integration' / 'pass1_labelled.h5ad')
    table = write_stage04(p2, P, out / '04_pseudospace')
    p2.uns['coast'] = nodes
    h5_safe(p2).write_h5ad(out / '03_integration' / 'pass2.h5ad')
    counts['pass 2 (curve)'] = table['sample'].value_counts().reindex(counts.index).fillna(0).astype(int)
    pt = table[table['coarse_class'] == 'PT']
    for seg in ('PT-S1', 'PT-S2', 'PT-S3'):
        counts[seg] = pt[pt['segment_class'] == seg]['sample'].value_counts().reindex(counts.index).fillna(0).astype(int)
    counts.loc['total'] = counts.sum()
    counts.to_csv(out / 'stage_counts.csv')
    (out / 'record.json').write_text(json.dumps({
        'coast_version': coast.__version__, 'coast_path': str(Path(coast.__file__).parent),
        'run': run, 'seconds': round(time.time() - t0),
        'pt_medians': pt.groupby('segment_class')['pseudospace'].median().round(4).to_dict()}, indent=1))
    print(counts.to_string())
    print(pt.groupby('segment_class')['pseudospace'].median().round(4).to_dict(), f'{time.time() - t0:.0f} s')


if __name__ == '__main__':
    build(parse_args())
