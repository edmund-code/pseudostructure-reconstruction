"""Fig. 6 placement check: are IR tubules placed on the axis by genes that do not respond to injury?

Injury-responsive genes (whole-PT pseudobulk: all four IR x control differences share a sign and are >= 0.5 log2) are
withdrawn from the integration HVG pool, and stages 03-04 are rerun once with the Fig. 3d machinery. Each section's
PT positions on that injury-blind axis are compared with the reference, against the gene-half A vs B agreement as
the noise of a rerun. Protocol (frozen before this ran): docs/paper/fig6-aki-along-pt.md, section 2.

  OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 python analysis/paper/fig6_placement.py
Writes <results-root>/paper_figures/fig6/placement/: responsive_genes.csv, injury_blind.parquet (+ .json, the
rerun), placement.csv (rho per section and label).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = next(p for p in Path(__file__).resolve().parents if (p / 'pseudospace').is_dir())
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from fig3_common import RUNS  # noqa: E402

CONTROLS, INJURED = ['Ctrl1A2', 'Ctrl1A4'], ['IR2A2', 'IR2A4']
PT = ['PT-S1', 'PT-S2', 'PT-S3']
MIN_DETECTED, RESPONSE_LOG2, TOLERANCE = 0.05, 0.5, 0.05
NAME = 'injury_blind'


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--data-root', default=os.environ.get('PSEUDOSPACE_DATA_ROOT', str(REPO / 'data')))
    p.add_argument('--results-root', default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT', str(REPO / 'results')))
    p.add_argument('--threads', type=int, default=4)
    p.add_argument('--tmp', default='/dev/shm' if Path('/dev/shm').is_dir() else '')
    return p.parse_args(argv)


def responsive_genes(root: Path) -> pd.DataFrame:
    """Whole-PT pseudobulk log2 counts per 10k per section, and the four IR - control differences."""
    import anndata as ad
    from scipy import sparse

    qc = ad.read_h5ad(root / RUNS['aki'] / '02_qc/units_qc.h5ad')
    ref = pd.read_parquet(root / RUNS['aki'] / '04_pseudospace/pseudospace.parquet')
    pt = ref[ref.segment_class.astype(str).isin(PT)]
    counts = sparse.csr_matrix(qc[pt.index.to_numpy()].layers['counts'])
    mt = qc.var['mt'].to_numpy(bool)
    sample = pt['sample'].astype(str).to_numpy()
    table = pd.DataFrame(index=qc.var_names)
    detected = np.zeros(qc.n_vars, bool)
    for s in CONTROLS + INJURED:
        c = counts[sample == s]
        lib = float(c[:, ~mt].sum())
        table[s] = np.log2(1e4 * (np.asarray(c.sum(0)).ravel() + 0.5) / lib)
        detected |= np.asarray((c > 0).mean(0)).ravel() >= MIN_DETECTED
    pairs = np.stack([table[i] - table[c] for i in INJURED for c in CONTROLS], 1)
    table['detected'] = detected & ~mt
    table['min_abs_pair'] = np.abs(pairs).min(1)
    table['condition'] = pairs.mean(1)
    same = (np.sign(pairs) == np.sign(pairs[:, :1])).all(1)
    table['responsive'] = table.detected & same & (table.min_abs_pair >= RESPONSE_LOG2)
    return table


def scores(root: Path, blind: pd.Series) -> pd.DataFrame:
    from scipy.stats import spearmanr

    ref = pd.read_parquet(root / RUNS['aki'] / '04_pseudospace/pseudospace.parquet')
    pt = ref[ref.segment_class.astype(str).isin(PT)]
    halves = root / 'paper_figures/fig3/robustness/aki/reruns'
    rows = []
    for s in CONTROLS + INJURED:
        for lab in ['all PT', *PT]:
            k = (pt['sample'].astype(str) == s) & (pt.segment_class.astype(str).isin(PT if lab == 'all PT' else [lab]))
            units = pt.index[k]
            both = units.intersection(blind.dropna().index)
            rho = spearmanr(pt.loc[both, 'pseudospace'], blind[both]).statistic
            hr = []
            for k5 in range(5):
                a = pd.read_parquet(halves / f'gene_half_s{k5}_A.parquet')['t'].reindex(units)
                b = pd.read_parquet(halves / f'gene_half_s{k5}_B.parquet')['t'].reindex(units)
                ok = a.notna() & b.notna()
                hr.append(spearmanr(a[ok], b[ok]).statistic)
            rows.append({'sample': s, 'units': lab, 'n': int(len(both)), 'rho_reference_vs_blind': float(rho),
                         'rho_half_A_vs_B_mean': float(np.mean(hr))})
    return pd.DataFrame(rows)


def main(argv=None):
    import fig3_robustness as fr

    args = parse_args(argv)
    root = Path(args.results_root)
    out = root / 'paper_figures/fig6/placement'
    out.mkdir(parents=True, exist_ok=True)
    table = responsive_genes(root)
    table.to_csv(out / 'responsive_genes.csv')
    genes = table.index[table.responsive].tolist()
    print(f'{table.detected.sum()} genes assessed; {len(genes)} injury-responsive '
          f'({int((table.condition[table.responsive] > 0).sum())} up, {int((table.condition[table.responsive] < 0).sum())} down)',
          flush=True)

    if not (out / f'{NAME}.parquet').exists():
        rerun = argparse.Namespace(data_root=args.data_root, results_root=args.results_root, run='aki',
                                   threads=args.threads, tmp=args.tmp)
        fr.run_job(rerun, NAME, out, spec={'condition': 'exclude', 'fraction': 1.0, 'seed': 0, 'half': '',
                                           'genes': genes})
    blind = pd.read_parquet(out / f'{NAME}.parquet')['t']
    sc = scores(root, blind)
    sc.to_csv(out / 'placement.csv', index=False)
    print(sc.round(3).to_string(index=False))
    verdict = {}
    for lab in ('all PT', 'PT-S3'):
        r = sc[sc.units == lab].set_index('sample')['rho_reference_vs_blind']
        floor = min(r[c] for c in CONTROLS) - TOLERANCE
        verdict[lab] = {'floor': floor, **{s: bool(r[s] >= floor) for s in INJURED}}
    passed = all(v[s] for v in verdict.values() for s in INJURED)
    (out / 'verdict.json').write_text(json.dumps({'n_responsive': len(genes), 'checks': verdict, 'pass': passed}, indent=1))
    print('placement check', 'PASSES' if passed else 'FAILS', json.dumps(verdict))


if __name__ == '__main__':
    main()
