"""Fig. 4b ranking: how much more of each gene's deviance the curve explains than straight lines through the clusters.

Same units, genes, axes and NB model as fig4_screen.py. Per gene and control mouse:
- straight lines: NB on [1, x, (x - c2)+], c2 = S2 centre on the gene's axis ((b12 + b23) / 2), i.e. a piecewise-linear
  curve through the three segment centres, extended linearly at both ends;
- curve: the screen's cubic B-spline plus (x - c2)+, so the straight lines are nested in it;
- dispersion per gene under the curve design; gain = (D_lines - D_curve) / D_intercept.
The rank uses the smaller gain of the two mice. Protocol and result: docs/paper/fig4-continuous-vs-discrete.md.

  python analysis/paper/fig4_linear_gain.py [--results-root <dir>] [--jobs 8]

Writes <results-root>/paper_figures/fig4/screen/linear_gain.csv (one row per screened gene, highest gain first).
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = next(p for p in Path(__file__).resolve().parents if (p / 'pseudospace').is_dir())
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from fig3_common import RUNS  # noqa: E402
from fig4_screen import BASIS_DF, MICE, MIN_DETECTED, N_SPLITS, PT, boundary, gene_halves, pattern_axis  # noqa: E402


def designs(x: np.ndarray, c2: float, knots) -> dict:
    from pseudospace.stats_gam import make_gam_design
    hinge = np.maximum(x - c2, 0)[:, None]
    return {'intercept': np.ones((len(x), 1)), 'lines': np.column_stack([np.ones_like(x), x, hinge]),
            'curve': np.column_stack([make_gam_design(x, knots), hinge])}


def main(argv=None):
    import anndata as ad
    from scipy import sparse
    from pseudospace.nb_gam import nb_dispersion, nb_fit
    from pseudospace.stats_gam import gam_internal_knots

    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--results-root', default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT', str(REPO / 'results')))
    p.add_argument('--jobs', type=int, default=8)
    args = p.parse_args(argv)
    root = Path(args.results_root)

    qc = ad.read_h5ad(root / RUNS['aki'] / '02_qc/units_qc.h5ad')
    ref = pd.read_parquet(root / RUNS['aki'] / '04_pseudospace/pseudospace.parquet')
    ref = ref.set_index('obs_name') if 'obs_name' in ref.columns else ref
    pt = ref[ref['sample'].astype(str).isin(MICE) & ref.segment_class.astype(str).isin(PT)]
    units = pt.index.to_numpy()
    counts = sparse.csr_matrix(qc[units].layers['counts'])
    mt = qc.var['mt'].to_numpy(bool)
    lib = np.asarray(counts[:, ~mt].sum(1)).ravel()
    sample = pt['sample'].astype(str).to_numpy()
    label = pt.segment_class.astype(str).map(PT).to_numpy()
    detected = np.ones(qc.n_vars, bool)
    for m in MICE:
        detected &= np.asarray((counts[sample == m] > 0).mean(0)).ravel() >= MIN_DETECTED
    keep = detected & ~mt
    genes = qc.var_names[keep]
    halves = gene_halves(qc.var)
    axes = {(k, h): pd.read_parquet(root / 'paper_figures/fig3/robustness/aki/reruns' / f'gene_half_s{k}_{h}.parquet')['t']
            .reindex(units).to_numpy(float) for k in range(N_SPLITS) for h in 'AB'}
    pattern = np.array([''.join('B' if g in halves[k] else 'A' for k in range(N_SPLITS)) for g in genes])

    rows = []
    for pat in sorted(set(pattern)):
        x, ok = pattern_axis(axes, pat)
        c2 = (boundary(x[ok], label[ok], 0) + boundary(x[ok], label[ok], 1)) / 2
        knots = gam_internal_knots(x[ok], BASIS_DF)
        gi = np.flatnonzero(pattern == pat)
        sub = counts[:, np.flatnonzero(keep)[gi]]
        fits = {}
        for m in MICE:
            k = ok & (sample == m)
            y, d = sub[k].toarray(), designs(x[k], c2, knots)
            w, off = np.ones(k.sum()), np.log(lib[k])
            disp = nb_dispersion(y, d['curve'], w, off, n_jobs=args.jobs)
            fits[m] = {n: nb_fit(y, d[n], w, off, disp, n_jobs=args.jobs) for n in d}
        for j, g in enumerate(genes[gi]):
            row = {'gene': g, 'pattern': pat, 'c2': c2}
            for m in MICE:
                f = fits[m]
                d0, dl, dc = (f[n]['deviance'][j] for n in ('intercept', 'lines', 'curve'))
                conv = bool(f['lines']['converged'][j] and f['curve']['converged'][j])
                row.update({f'gain_{m}': (dl - dc) / d0 if conv else np.nan, f'lines_explained_{m}': 1 - dl / d0,
                            **{f'lines_b{i}_{m}': f['lines']['beta'][j, i] for i in range(3)}})
            rows.append(row)
        print(pat, len(gi), 'genes', flush=True)

    df = pd.DataFrame(rows)
    df['gain'] = df[[f'gain_{m}' for m in MICE]].min(1)
    # nested fits: the curve can only lower the deviance; a clearly negative gain means a fit failed
    assert (df['gain'].dropna() > -1e-6).all(), 'curve deviance above the straight lines: nested fit failed'
    df.sort_values('gain', ascending=False).to_csv(root / 'paper_figures/fig4/screen/linear_gain.csv', index=False)
    print(df.sort_values('gain', ascending=False).head(15)[['gene', 'gain']].to_string(index=False))


if __name__ == '__main__':
    main()
