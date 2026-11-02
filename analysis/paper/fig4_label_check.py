"""Post-hoc check of the Fig. 4 calls (not in the frozen protocol): do they hold within a single label?

Near a boundary, S1- and S2-labelled units overlap on x. A gene that only steps between labels can then look graded
on x. This re-measures each call on units of one label only, from raw counts (log2 of summed counts per 10k):
 category 2  change between the first and last third of the segment core, units with that segment's label
 category 1  (peak at b23) S2 units: last third of the S2 core minus its first third (should rise);
             S3 units: first third of the S3 core minus its last third (should fall)

  python analysis/paper/fig4_label_check.py
Writes <results-root>/paper_figures/fig4/screen/label_check.csv. Protocol: docs/paper/fig4-continuous-vs-discrete.md.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fig4_screen as fs  # noqa: E402


def third_rates(y, lib, x, lo, hi):
    third = (hi - lo) / 3
    a, b = (x >= lo) & (x <= lo + third), (x >= hi - third) & (x <= hi)
    rate = lambda k: np.log2(1e4 * (y[k].sum() + 0.5) / lib[k].sum())  # noqa: E731
    return rate(a), rate(b)


def main(argv=None):
    import anndata as ad
    from scipy import sparse

    args = fs.parse_args(argv)
    root = Path(args.results_root)
    out = root / args.out
    genes = pd.read_csv(out / 'genes.csv')
    calls = genes[genes.cat1 | genes.cat2].copy()

    qc = ad.read_h5ad(root / fs.RUNS['aki'] / '02_qc/units_qc.h5ad')
    ref = pd.read_parquet(root / fs.RUNS['aki'] / '04_pseudospace/pseudospace.parquet')
    ref = ref.set_index('obs_name') if 'obs_name' in ref.columns else ref
    pt = ref[ref['sample'].astype(str).isin(fs.MICE) & ref.segment_class.astype(str).isin(fs.PT)]
    units = pt.index.to_numpy()
    counts = sparse.csc_matrix(qc[units].layers['counts'])
    lib = np.asarray(counts[:, ~qc.var['mt'].to_numpy(bool)].sum(1)).ravel()
    sample = pt['sample'].astype(str).to_numpy()
    label = pt.segment_class.astype(str).map(fs.PT).to_numpy()
    axes = {(k, h): pd.read_parquet(root / 'paper_figures/fig3/robustness/aki/reruns' /
                                    f'gene_half_s{k}_{h}.parquet')['t'].reindex(units).to_numpy(float)
            for k in range(fs.N_SPLITS) for h in 'AB'}

    rows = []
    for pat, grp in calls.groupby('pattern'):
        pos = np.stack([axes[(k, pat[k])] for k in range(fs.N_SPLITS)], 1)
        ok = np.isfinite(pos).sum(1) >= fs.MIN_AXES
        x = np.full(len(units), np.nan)
        r = pd.Series(np.nanmean(pos[ok], 1)).rank(method='average').to_numpy()
        x[ok] = (r - 0.5) / ok.sum()
        for _, g in grp.iterrows():
            b12, b23 = g.b12, g.b23
            cores = [(0.0, b12 - fs.MARGIN), (b12 + fs.MARGIN, b23 - fs.MARGIN), (b23 + fs.MARGIN, 1.0)]
            y = counts[:, qc.var_names.get_loc(g.gene)].toarray().ravel()
            row = {'gene': g.gene, 'category': 1 if g.cat1 else 2}
            for m in fs.MICE:
                k = ok & (sample == m)
                if g.cat1:
                    s2 = k & (label == 1)
                    first, last = third_rates(y[s2], lib[s2], x[s2], *cores[1])
                    s3 = k & (label == 2)
                    f3, l3 = third_rates(y[s3], lib[s3], x[s3], *cores[2])
                    row.update({f'rise_in_S2_{m}': last - first, f'fall_in_S3_{m}': f3 - l3})
                else:
                    seg = int(g[f'change_segment_{m}'])
                    kk = k & (label == seg)
                    first, last = third_rates(y[kk], lib[kk], x[kk], *cores[seg])
                    row.update({f'segment_{m}': seg + 1, f'within_label_change_{m}': last - first,
                                f'curve_change_{m}': g[f'change_{m}']})
            rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(out / 'label_check.csv', index=False)
    a, b = fs.MICE
    c2 = df[df.category == 2]
    same = np.sign(c2[f'within_label_change_{a}']) == np.sign(c2[f'curve_change_{a}'])
    same &= np.sign(c2[f'within_label_change_{b}']) == np.sign(c2[f'curve_change_{b}'])
    size = c2[[f'within_label_change_{a}', f'within_label_change_{b}']].abs().min(1)
    print('category 2:', len(c2), 'calls; same sign within label in both mice:', int(same.sum()),
          '| of those >= 1 log2:', int((same & (size >= 1)).sum()), '| >= 0.5 log2:', int((same & (size >= 0.5)).sum()))
    print('median within-label / curve change:',
          float(np.median(c2[f'within_label_change_{a}'] / c2[f'curve_change_{a}'])),
          float(np.median(c2[f'within_label_change_{b}'] / c2[f'curve_change_{b}'])))
    c1 = df[df.category == 1]
    print(c1.round(2).to_string(index=False))


if __name__ == '__main__':
    main()
