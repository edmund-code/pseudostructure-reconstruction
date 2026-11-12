"""Is the human PT less zonated than the mouse PT in external data? (post hoc, 2026-11-11)

S3-vs-S1 log2 fold changes per gene, human against mouse, in our sections and in external sources: Cheval 2012 SAGE
on microdissected S1/S3 (both species, one lab), Lake 2023 human snRNA (normal donors; cortex and medulla S3),
CELLxGENE census mouse snRNA, and Chen 2021 mouse microdissection. Same-species pairs calibrate how much two datasets
disagree on steepness. Usage: python attenuation_external_check.py <repo> <lake_region_pseudobulk.csv>
Note: docs/results/human-gradient-attenuation-external.md."""
import sys
from pathlib import Path
import numpy as np
import pandas as pd


def lake_region_pseudobulk(h5ad: Path, out: Path):
    """Normal-donor Lake 2023 PT nuclei summed by segment and region (raw counts), written to out."""
    import anndata as ad
    a = ad.read_h5ad(h5ad, backed='r')
    o = a.obs
    pt = o['subclass.l2'].astype(str).isin(['PT-S1', 'PT-S2', 'PT-S3']) & o.disease.astype(str).eq('normal')
    grp = (o.loc[pt, 'subclass.l2'].astype(str) + '|' + o.loc[pt, 'region'].astype(str)).to_numpy()
    idx = np.flatnonzero(pt.to_numpy())
    X = a.raw.X if a.raw is not None else a.X
    sums = {}
    for name in ['PT-S1|C', 'PT-S3|C', 'PT-S3|M', 'PT-S3|C-M']:
        k = np.sort(idx[grp == name])
        sums[name] = sum(np.asarray(X[k[i:i + 4000]].sum(0)).ravel() for i in range(0, len(k), 4000))
    var = a.raw.var if a.raw is not None else a.var
    pd.DataFrame(sums, index=var['feature_name'].astype(str).values).groupby(level=0).sum().to_csv(out)


repo = Path(sys.argv[1])
# argv[2]: Lake region pseudobulk CSV; built with lake_region_pseudobulk() when missing
ext = repo / 'data/external'
g = pd.read_csv(repo / 'results/paper_figures/fig5/screen/genes.csv', index_col=0)
fc = {}
fc['ours mouse'] = np.mean([g[f'cluster_S3_{m}'] - g[f'cluster_S1_{m}'] for m in ('Ctrl1A2', 'Ctrl1A4')], 0)
fc['ours human'] = np.mean([g[f'cluster_S3_{m}'] - g[f'cluster_S1_{m}'] for m in ('HUK1_COR1', 'HUK1_MED1')], 0)
fc = pd.DataFrame(fc, index=g.index)
assess = g.index[g.assessable.astype(bool)]

# Cheval 2012: same lab and method (SAGE on microdissected S1 and S3) in both species, tags per 10,000
c = pd.read_csv(ext / 'cheval2012/cheval2012_table_s1.csv')
ms = pd.read_csv(ext / 'cheval2012/cheval2012_mouse_symbols.csv')
c = c.merge(ms, on='human_symbol', how='left') if 'mouse_symbol' in ms.columns else c
key = 'mouse_symbol' if 'mouse_symbol' in c.columns else None
c['gene'] = c[key] if key else c.human_symbol.str.capitalize()
c = c.drop_duplicates('gene').set_index('gene')
for sp in ('human', 'mouse'):
    c[f'chev_{sp}'] = np.log2(c[f'{sp}_S3'] + 1) - np.log2(c[f'{sp}_S1'] + 1)
c['chev_level'] = c[['human_S1', 'human_S3', 'mouse_S1', 'mouse_S3']].min(1).clip(lower=0) + 0
c['chev_ok'] = (c[['human_S1', 'human_S3']].max(1) >= 5) & (c[['mouse_S1', 'mouse_S3']].max(1) >= 5)

# Lake 2023 human snRNA, cortex, per-segment means (log2 CPM); mouse microdissection (Chen 2021, TPM)
lake = pd.read_csv(repo / 'results/pt_literature_deep/lake_cortex_segment_means.csv', index_col=0)
lake.index = lake.index.str.capitalize()
lake_fc = lake.lake_log2cpm_S3 - lake.lake_log2cpm_S1
lake_lvl = lake[['lake_log2cpm_S1', 'lake_log2cpm_S3']].max(1)
mb = pd.read_csv(ext / 'segment_microdissection/GSE150338_Mouse_TPM_Replicates.txt.gz', sep='\t').drop_duplicates('Gene_symbol').set_index('Gene_symbol')
seg = {s: np.log2(mb[[k for k in mb.columns if f'_PT{s}_' in k]].apply(pd.to_numeric, errors='coerce').clip(lower=0) + 1).mean(1) for s in ('S1', 'S3')}
chen_fc, chen_lvl = seg['S3'] - seg['S1'], pd.concat(seg, axis=1).max(1)

# census mouse snRNA (adult), pooled pseudobulk S1 and S3
cm = pd.read_csv(ext / 'census_pt_segments/mouse_pt_segment_pseudobulk_counts.csv.gz', index_col=0)
meta = pd.read_csv(ext / 'census_pt_segments/mouse_pt_segment_pseudobulk_meta.csv', index_col=0)
meta = meta[~meta.development_stage.str.contains('embryo|E16|Theiler', case=False)]
def pb(label):
    cols = meta.index[meta.cell_type.str.contains(f'segment {label}')]
    cols = [k for k in cols if k in cm.columns] or [k for k in cols if k in cm.index]
    sub = cm[cols] if cols and cols[0] in cm.columns else cm.loc[cols].T
    return np.log2(1e6 * sub.sum(1) / sub.sum(1).sum() + 1)
try:
    cs1, cs3 = pb('1'), pb('3')
    sn_fc, sn_lvl = cs3 - cs1, pd.concat([cs1, cs3], axis=1).max(1)
except Exception as e:
    print('census mouse failed:', e); sn_fc = sn_lvl = None


ROWS = []


def amp(h, m, name, n_min=30, group=None, key=None):
    ok = h.notna() & m.notna()
    h, m = h[ok], m[ok]
    sd = h.std() / m.std()
    # total-least-squares slope of human on mouse (symmetric in noise)
    cov = np.cov(m, h)
    lam = np.linalg.eigh(cov)
    v = lam[1][:, -1]
    tls = v[1] / v[0]
    from scipy.stats import spearmanr
    print(f'{name:58s} n={len(h):5d}  SD ratio h/m={sd:5.2f}  TLS slope={tls:5.2f}  rho={spearmanr(h, m).statistic:5.2f}')
    if key:
        ROWS.append({'key': key, 'group': group, 'label': name.strip(), 'n': len(h), 'sd_ratio': sd})


print('S3 - S1 log2 fold change, human against mouse, on genes in our assessable universe\n')
base = pd.Index(assess)
ourok = fc.loc[base]
amp(ourok['ours human'], ourok['ours mouse'], 'OURS: human sections vs mouse sections', group='human vs mouse', key='ours')

common = base.intersection(c.index[c.chev_ok])
amp(c.loc[common, 'chev_human'], c.loc[common, 'chev_mouse'], 'CHEVAL 2012 (SAGE, microdissection): human vs mouse')
amp(fc.loc[common, 'ours human'], fc.loc[common, 'ours mouse'], '  ours on the same genes')

if sn_fc is not None:
    common = base.intersection(lake.index[lake_lvl >= 4]).intersection(sn_lvl.index[sn_lvl >= 4])
    amp(lake_fc[common], sn_fc[common], 'SNRNA: Lake human cortex vs census mouse')
    amp(fc.loc[common, 'ours human'], fc.loc[common, 'ours mouse'], '  ours on the same genes')
common = base.intersection(lake.index[lake_lvl >= 4]).intersection(chen_lvl.index[chen_lvl >= 3])
amp(lake_fc[common], chen_fc[common], 'MIXED: Lake human snRNA vs Chen mouse microdissection')
amp(fc.loc[common, 'ours human'], fc.loc[common, 'ours mouse'], '  ours on the same genes')

# genes zonated in mouse by an independent source (Chen microdissection |FC| >= 2), ratio of folds per gene
print('\nGenes zonated in Chen mouse microdissection (|S3-S1| >= 2 log2), median |human| / |mouse| fold per source:')
z = base.intersection(chen_fc.index[(chen_fc.abs() >= 2) & (chen_lvl >= 3)])
for name, h, m in (('ours', fc['ours human'], fc['ours mouse']),
                   ('Cheval', c['chev_human'].where(c.chev_ok), c['chev_mouse'].where(c.chev_ok)),
                   ('Lake vs census mouse', lake_fc.where(lake_lvl >= 4), sn_fc.where(sn_lvl >= 4) if sn_fc is not None else None)):
    if m is None:
        continue
    k = z.intersection(h.dropna().index).intersection(m.dropna().index)
    same = np.sign(h[k]) == np.sign(m[k])
    r = (h[k].abs() / m[k].abs().clip(lower=0.25))[same]
    print(f'  {name:22s} genes={len(k):4d} same direction={same.mean():.2f}  median |h|/|m|={r.median():.2f}')

# Lake by region: is the S3 outside the cortex steeper? (tests whether cortex-only tissue misses the PT's end)
lake_pb = Path(sys.argv[2])
if not lake_pb.exists():
    lake_region_pseudobulk(ext / 'kpmp_lake2023/lake2023_snRNA_v1.5.h5ad', lake_pb)
pbr = pd.read_csv(lake_pb, index_col=0)
pbr.index = pbr.index.str.capitalize()
cpm = np.log2(1e6 * pbr / pbr.sum() + 1)
lvl = cpm[['PT-S1|C', 'PT-S3|C', 'PT-S3|M']].max(1)
print('\nLake 2023 by region (normal donors), against census mouse snRNA S3 - S1:')
common = base.intersection(cpm.index[lvl >= 4]).intersection(sn_lvl.index[sn_lvl >= 4])
for late in ('PT-S3|C', 'PT-S3|C-M', 'PT-S3|M'):
    amp(cpm.loc[common, late] - cpm.loc[common, 'PT-S1|C'], sn_fc[common], f'  human {late:10s} - S1 cortex vs mouse',
        group='human vs mouse', key={'PT-S3|C': 'atlas_cortex', 'PT-S3|M': 'atlas_medulla'}.get(late))
amp(fc.loc[common, 'ours human'], fc.loc[common, 'ours mouse'], '  ours on the same genes')
k = z.intersection(common)
for late in ('PT-S3|C', 'PT-S3|M'):
    h = cpm.loc[k, late] - cpm.loc[k, 'PT-S1|C']; m = sn_fc[k]
    same = np.sign(h) == np.sign(m)
    print(f'  Chen-zonated genes, {late}: n={len(k)} same dir={same.mean():.2f} median |h|/|m|={(h.abs()/m.abs().clip(lower=.25))[same].median():.2f}')

print('\nCalibration, same species in different datasets (should be ~1 if datasets agree on steepness):')
common = base.intersection(sn_lvl.index[sn_lvl >= 4]).intersection(chen_lvl.index[chen_lvl >= 3])
amp(sn_fc[common], chen_fc[common], '  census mouse snRNA vs Chen mouse microdissection', group='same species',
    key='mouse_atlas_vs_micro')
amp(fc.loc[common, 'ours mouse'], sn_fc[common], '  ours mouse vs census mouse snRNA', group='same species',
    key='ours_mouse_vs_atlas')
amp(fc.loc[common, 'ours mouse'], chen_fc[common], '  ours mouse vs Chen mouse microdissection')
common = base.intersection(lake.index[lake_lvl >= 4])
amp(fc.loc[common, 'ours human'], lake_fc[common], '  ours human vs Lake human cortex', group='same species',
    key='ours_human_vs_atlas')
k = z.intersection(common).intersection(sn_lvl.index[sn_lvl >= 4])
for name, h, m in (('census mouse / Chen', sn_fc, chen_fc), ('ours mouse / census mouse', fc['ours mouse'], sn_fc),
                   ('ours human / Lake cortex', fc['ours human'], lake_fc)):
    hh, mm = h[k], m[k]; same = np.sign(hh) == np.sign(mm)
    print(f'  Chen-zonated genes, {name:26s}: n={len(k)} same dir={same.mean():.2f} median ratio={(hh.abs()/mm.abs().clip(lower=.25))[same].median():.2f}')

out = repo / 'results/human_gradient_attenuation'
out.mkdir(parents=True, exist_ok=True)
pd.DataFrame(ROWS).to_csv(out / 'amplitude_ratios.csv', index=False)
