# Usage: python analysis/paper/fig56_annotate.py <out_dir>  (run from the repo root; inputs under results/ and data/external/)
"""Annotate every gene of the Fig. 5 and Fig. 6 screens: human ortholog, Open Targets drugs and literature,
PanelApp tubulopathy genes (green), kidney GWAS loci (eGFR, urate, urine metabolites)."""
import json, re, time, gzip
from pathlib import Path
import numpy as np, pandas as pd, requests, anndata as ad
R = Path('results'); OUT = Path(__import__('sys').argv[1])
EXT = Path('data/external')
f5 = pd.read_csv(R / 'paper_figures/fig5/screen/genes.csv'); f6 = pd.read_csv(R / 'paper_figures/fig6/screen/genes.csv')
h5 = pd.read_csv(R / 'paper_figures/fig5/screen/handoffs.csv'); h6 = pd.read_csv(R / 'paper_figures/fig6/screen/handoffs.csv')
mouse = sorted(set(f5.gene) | set(f6.gene) | set(h5.early) | set(h5.late) | set(h6.early) | set(h6.late))
var = ad.read_h5ad(R / 'package_pseudospace/runs/aki/02_qc/units_qc.h5ad', backed='r').var
mid = var['gene_ids'].reindex(mouse)
comp = pd.read_csv(EXT / 'ensembl_compara_human_mouse/Compara.98.protein_default.homologies.human_mouse.tsv.gz', sep='\t',
                   usecols=['gene_stable_id', 'homology_type', 'homology_gene_stable_id', 'identity'])
comp['rank'] = comp.homology_type.map({'ortholog_one2one': 0, 'ortholog_one2many': 1, 'ortholog_many2many': 2}).fillna(3)
comp = comp.sort_values(['homology_gene_stable_id', 'rank', 'identity'], ascending=[True, True, False]).drop_duplicates('homology_gene_stable_id')
comp = comp.set_index('homology_gene_stable_id')
ann = pd.DataFrame({'mouse': mouse, 'mouse_id': mid.to_numpy()})
ann['human_id'] = comp.gene_stable_id.reindex(ann.mouse_id).to_numpy()
ann['homology'] = comp.homology_type.reindex(ann.mouse_id).to_numpy()
ids = ann.human_id.dropna().unique().tolist()
Q = '''{ targets(ensemblIds: %s) { id approvedSymbol approvedName
  drugAndClinicalCandidates { count rows { maxClinicalStage drug { name } } }
  literatureOcurrences { count } } }'''
STAGE = {'APPROVAL': 4, 'PHASE_4': 4, 'PHASE_3': 3, 'PHASE_2_3': 2.5, 'PHASE_2': 2, 'PHASE_1_2': 1.5, 'PHASE_1': 1,
         'EARLY_PHASE_1': 0.5, 'PRECLINICAL': 0.1}
rows = {}
for i in range(0, len(ids), 40):
    for attempt in range(4):
        try:
            r = requests.post('https://api.platform.opentargets.org/api/v4/graphql', json={'query': Q % json.dumps(ids[i:i + 40])}, timeout=60)
            r.raise_for_status(); data = r.json()['data']['targets']; break
        except Exception as e:
            time.sleep(3 * (attempt + 1)); data = None
    for t in data or []:
        dc = t['drugAndClinicalCandidates']
        stages = [(STAGE.get(x['maxClinicalStage'], 0), x['drug']['name']) for x in dc['rows']]
        stages.sort(reverse=True)
        rows[t['id']] = {'human': t['approvedSymbol'], 'human_name': t['approvedName'], 'n_drugs': dc['count'],
                         'max_phase': stages[0][0] if stages else 0,
                         'top_drugs': '; '.join(dict.fromkeys(n for _, n in stages[:6])),
                         'lit_count': t['literatureOcurrences']['count']}
    print(i, flush=True)
ot = pd.DataFrame.from_dict(rows, orient='index')
ann = ann.join(ot, on='human_id')
# PanelApp green tubulopathy / stone genes
panel = set()
for pid in (292, 149, 3993):
    for g in json.loads((EXT / f'panelapp_2026-10-18/panel_{pid}.json').read_text())['genes']:
        if str(g.get('confidence_level')) == '3':
            panel.add(g['entity_name'])
def gwas(pmids):
    out = set()
    for pmid in pmids:
        for f in sorted((EXT / 'gwas_catalog_v2_2026-10-18').glob(f'assoc_{pmid}_p*.json')):
            for a in json.loads(f.read_text()).get('_embedded', {}).get('associations', []):
                for mg in a.get('mapped_genes') or []:
                    out |= {x.strip() for x in re.split(r'[,;]| - ', str(mg)) if x.strip()}
    return out
G = {'gwas_eGFR': gwas([34272381]), 'gwas_urate': gwas([31578528]), 'gwas_urine_metab': gwas([31959995, 37277652])}
sym = ann.human.fillna(ann.mouse.str.upper())
ann['tubulopathy_panel'] = sym.isin(panel)
for k, v in G.items():
    ann[k] = sym.isin(v)
ann.to_csv(OUT / 'annotation.csv', index=False)
print(ann.notna().sum().to_dict()); print('drug targets (any phase):', (ann.n_drugs > 0).sum(), 'approved:', (ann.max_phase >= 4).sum())
