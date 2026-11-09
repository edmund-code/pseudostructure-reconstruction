# Usage: python analysis/paper/fig56_rank.py <out_dir>  (run from the repo root; inputs under results/ and data/external/)
import numpy as np, pandas as pd, sys
from pathlib import Path
pd.set_option('display.width', 250, 'display.max_colwidth', 60, 'display.max_rows', 200)
S = Path(sys.argv[1]); R = Path('results/paper_figures')
ann = pd.read_csv(S / 'annotation.csv').drop_duplicates('mouse').set_index('mouse')
ann['clin'] = (3 * (ann.max_phase >= 4) + 1.5 * ((ann.max_phase > 0) & (ann.max_phase < 4)) + 2 * ann.tubulopathy_panel
               + 1 * (ann.gwas_eGFR | ann.gwas_urate | ann.gwas_urine_metab) + np.log10(ann.lit_count.fillna(0) + 1) / 2)
A = ['human', 'clin', 'max_phase', 'top_drugs', 'tubulopathy_panel', 'gwas_eGFR', 'gwas_urate', 'gwas_urine_metab', 'lit_count']
def a(df, col='gene'): return df.join(ann[A], on=col)

# ---------------- Fig 5
g5 = pd.read_csv(R / 'fig5/screen/genes.csv').set_index('gene')
M, H = ['Ctrl1A2', 'Ctrl1A4'], ['HUK1_COR1', 'HUK1_MED1']
g5['E_mouse'] = g5[[f'E_{m}' for m in M]].mean(1); g5['E_human'] = g5[[f'E_{h}' for h in H]].mean(1)
calls = g5[g5.catR | g5.catH].copy(); calls['cat'] = np.where(calls.catR, 'reversed', 'steeper in human')
c5 = a(calls.reset_index())[['gene', 'cat', 'E_mouse', 'E_human'] + A].sort_values('clin', ascending=False)
c5.to_csv(S / 'fig5_calls_annotated.csv', index=False)
print('=== Fig 5 calibrated gene calls, by clinical score'); print(c5.round(2).head(35).to_string(index=False))
# descriptive class for every clinically relevant screened gene
ok = g5.screened & g5.assessable
e = g5.loc[ok, [f'E_{s}' for s in M + H]]
cls = pd.Series('mixed', index=e.index)
sm, sh = np.sign(e[[f'E_{m}' for m in M]]), np.sign(e[[f'E_{h}' for h in H]])
big = (e.abs() >= .5)
cons = (sm.nunique(1) == 1) & (sh.nunique(1) == 1) & (sm.iloc[:, 0] == sh.iloc[:, 0]) & big.all(1)
cls[cons] = 'conserved'
cls[g5.loc[ok, 'catR']] = 'reversed (call)'
mouse_only = big[[f'E_{m}' for m in M]].all(1) & (sm.nunique(1) == 1) & (e[[f'E_{h}' for h in H]].abs() < .25).all(1)
cls[mouse_only & (cls == 'mixed')] = 'mouse gradient, human flat'
human_only = big[[f'E_{h}' for h in H]].all(1) & (sh.nunique(1) == 1) & (e[[f'E_{m}' for m in M]].abs() < .25).all(1)
cls[human_only & (cls == 'mixed')] = 'human gradient, mouse flat'
d5 = a(pd.DataFrame({'gene': e.index, 'class': cls.values, 'E_mouse': g5.loc[ok, 'E_mouse'].values, 'E_human': g5.loc[ok, 'E_human'].values}))
d5c = d5[(d5.max_phase > 0) | d5.tubulopathy_panel | d5.gwas_urate | d5.gwas_eGFR | d5.gwas_urine_metab]
d5c.sort_values(['class', 'clin'], ascending=[True, False]).to_csv(S / 'fig5_clinical_genes_by_class.csv', index=False)
print('\n=== Fig 5: clinically annotated screened genes by class'); print(d5c['class'].value_counts().to_dict())
for k in ['human gradient, mouse flat', 'mouse gradient, human flat', 'conserved']:
    print('--', k); print(d5c[d5c['class'] == k].sort_values('clin', ascending=False).head(15)[['gene', 'E_mouse', 'E_human', 'human', 'max_phase', 'top_drugs', 'tubulopathy_panel', 'gwas_urate', 'lit_count']].round(2).to_string(index=False))
h5 = pd.read_csv(R / 'fig5/screen/interpret_handoffs.csv')
h5 = h5.join(ann[['human', 'max_phase', 'top_drugs']].add_suffix('_early'), on='early').join(ann[['human', 'max_phase', 'top_drugs']].add_suffix('_late'), on='late')
print('\n=== Fig 5 shifted handoffs'); print(h5[['early', 'late', 'human_early', 'human_late', 'shift', 'beyond_warp', 'max_phase_early', 'top_drugs_early', 'max_phase_late', 'top_drugs_late']].round(2).to_string(index=False))
h5.to_csv(S / 'fig5_handoffs_annotated.csv', index=False)

# ---------------- Fig 6
g6 = pd.read_csv(R / 'fig6/screen/genes.csv').set_index('gene')
C, I = ['Ctrl1A2', 'Ctrl1A4'], ['IR2A2', 'IR2A4']
ok6 = g6.universe & g6.neighbour_ok & g6.converged if 'universe' in g6 else g6.neighbour_ok
rows = []
for gname, r in g6[ok6].iterrows():
    cl = {m: np.array([r[f'cluster_S{s}_{m}'] for s in (1, 2, 3)]) for m in C + I}
    if not all(np.isfinite(v).all() for v in cl.values()):
        continue
    pairs = [cl[i] - cl[c] for i in I for c in C]            # IR minus control, 4 pairs, per segment
    cond = np.mean(pairs, 0)
    rel1 = (cl['IR2A2'] + cl['Ctrl1A2']) / 2 - (cl['IR2A4'] + cl['Ctrl1A4']) / 2
    rel2 = (cl['IR2A2'] + cl['Ctrl1A4']) / 2 - (cl['IR2A4'] + cl['Ctrl1A2']) / 2
    dep = cond[2] - cond[0]                                  # S3 response minus S1 response
    pdep = [p[2] - p[0] for p in pairs]
    rdep = [rel1[2] - rel1[0], rel2[2] - rel2[0]]
    rows.append({'gene': gname, 'resp_S1': cond[0], 'resp_S2': cond[1], 'resp_S3': cond[2], 'S3_minus_S1': dep,
                 'pairs_agree': all(np.sign(x) == np.sign(dep) for x in pdep) and min(abs(x) for x in pdep) >= .5,
                 'specific': all(abs(x) <= .25 * abs(dep) for x in rdep), 'catA': bool(r.catA)})
d6 = a(pd.DataFrame(rows))
d6['segment_dependent'] = d6.pairs_agree & d6.specific & (d6.S3_minus_S1.abs() >= 1)
d6.to_csv(S / 'fig6_segment_response.csv', index=False)
print('\n=== Fig 6: injury response differs S3 vs S1 (|Δ|>=1, 4 pairs agree >=0.5, relabelings <=0.25x):', d6.segment_dependent.sum(), 'of', len(d6))
clin6 = d6[d6.segment_dependent & ((d6.max_phase > 0) | d6.tubulopathy_panel | d6.gwas_eGFR | d6.gwas_urate | d6.gwas_urine_metab)]
print(clin6.sort_values('clin', ascending=False)[['gene', 'human', 'resp_S1', 'resp_S2', 'resp_S3', 'S3_minus_S1', 'max_phase', 'top_drugs', 'tubulopathy_panel', 'lit_count']].round(2).head(40).to_string(index=False))
h6 = pd.read_csv(R / 'fig6/screen/handoffs.csv'); h6 = h6[h6.shifted]
h6 = h6.join(ann[['human', 'max_phase', 'top_drugs']].add_suffix('_early'), on='early').join(ann[['human', 'max_phase', 'top_drugs']].add_suffix('_late'), on='late')
print('\n=== Fig 6 shifted handoffs'); print(h6[['early', 'late', 'human_early', 'human_late', 'shift', 'max_phase_early', 'top_drugs_early', 'max_phase_late', 'top_drugs_late']].round(2).to_string(index=False))
h6.to_csv(S / 'fig6_handoffs_annotated.csv', index=False)
ig = a(pd.read_csv(R / 'fig6/screen/interpret_genes.csv'))
print('\n=== Fig 6 named within-segment genes'); print(ig.sort_values('clin', ascending=False)[['gene', 'segment', 'W', 'W_control', 'W_IR', 'human', 'max_phase', 'top_drugs', 'tubulopathy_panel', 'lit_count']].round(2).to_string(index=False))
