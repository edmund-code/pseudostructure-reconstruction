"""Fig. 5 follow-up (post hoc, descriptive): the named checks of the protocol, run on every call of fig5_screen.py.

Genes (category R and H calls):
 admixture  per section, Spearman rho over the 20 equal-count bins between the gene's rate and a distal-identity
            rate (AL, DCT and CD markers) on the gene's own axis (descriptive: the distal signal itself rises along
            the PT in both species, so rho mostly restates the gene's direction). The flag (post hoc) is a bound:
            the distal admixture fraction f of the S1 and S3 cores (marker rate in PT units over the rate in that
            neighbour's own units, written to interpret_admixture_fraction.csv) caps the S3-minus-S1 change that
            admixture can add to a gene whose neighbour rate is r times its PT rate at log2(1 + f3 r) - log2(1 + f1 r).
            Admixture can only raise late-PT expression, so it can never produce a fall.
 low mito   raw R of each human section recomputed on its lowest-mito tertile (stored by the screen): same sign?
 cortical   (post hoc) both human sections are cortex, while much mouse S3 lies in the outer stripe: each mouse's raw R
            recomputed with its S3-core units restricted to cortical ones (within the glomerulus distance that holds
            90% of that section's S1-core units, as the 2026-10-21 check); same sign?
Pathways (calls):
 raw        the set's mean beyond-flattening residual computed on raw R instead of curves (same OLS construction)
 low mito   the same with the human sections' low-mito raw R
 drop3      the set's d without its three largest same-sign members, and its p against random sets of that size
 leading    the five members furthest in the call's direction, and the top-3 share
Handoffs (shifted, mouse-only and human-only pairs): what moves, read from each partner's E per section.

  python analysis/paper/fig5_interpret.py
Reads <results-root>/paper_figures/fig5/screen/ and writes interpret_{genes,pathways,handoffs}.csv there.
Protocol: docs/paper/fig5-human-vs-mouse-along-pt.md.
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
import fig4_screen as fs  # noqa: E402
import fig5_screen as f5  # noqa: E402
from fig3_common import RUNS  # noqa: E402

N_NULL = 2000
DISTAL = ['Slc12a1', 'Umod', 'Slc12a3', 'Calb1', 'Pvalb', 'Aqp2', 'Kcnj1', 'Cldn16']
NEIGHBOUR_MARKERS = {'AL': ['Umod', 'Slc12a1'], 'DCT': ['Slc12a3'], 'CD': ['Aqp2']}


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--results-root', default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT', str(REPO / 'results')))
    return p.parse_args(argv)


def residual(E: dict, valid: np.ndarray) -> np.ndarray:
    """Mean over the four human x mouse pairs of the OLS residual of E(human) on E(mouse) (screen's construction)."""
    res = []
    for a in f5.HUMAN:
        for b in f5.MOUSE:
            k = valid & np.isfinite(E[a]) & np.isfinite(E[b])
            slope, icpt = np.polyfit(E[b][k], E[a][k], 1)
            r = np.full(len(valid), np.nan)
            r[k] = E[a][k] - (icpt + slope * E[b][k])
            res.append(r)
    return np.mean(res, 0)


def main(argv=None):
    import anndata as ad
    from scipy import sparse
    from scipy.spatial import cKDTree
    from scipy.stats import spearmanr

    args = parse_args(argv)
    root = Path(args.results_root)
    screen = root / 'paper_figures/fig5/screen'
    genes = pd.read_csv(screen / 'genes.csv', index_col=0)
    ok = genes.assessable.to_numpy(bool)
    rng = np.random.default_rng(5)

    # ---------------------------------------------------------------------------------------------- genes
    qc = ad.read_h5ad(root / RUNS['hvm'] / '02_qc/units_qc.h5ad')
    ref = pd.read_parquet(root / RUNS['hvm'] / '04_pseudospace/pseudospace.parquet')
    ref = ref.set_index('obs_name') if 'obs_name' in ref.columns else ref
    pt = ref[ref['sample'].astype(str).isin(f5.SECTIONS) & ref.segment_class.astype(str).isin(fs.PT)]
    units = pt.index.to_numpy()
    counts = sparse.csc_matrix(qc[units].layers['counts'])
    mt = qc.var['mt'].to_numpy(bool)
    lib = np.asarray(counts[:, ~mt].sum(1)).ravel()
    sample = pt['sample'].astype(str).to_numpy()
    xy = qc.obs.loc[units, ['centroid_x', 'centroid_y']].to_numpy(float)
    p1 = ad.read_h5ad(root / RUNS['hvm'] / '03_integration/pass1_labelled.h5ad', backed='r').obs
    glom = np.full(len(units), np.nan)
    for m in f5.SECTIONS:
        gl = p1[(p1['sample'].astype(str) == m) & (p1.segment_class.astype(str) == 'Podocyte')]
        tree = cKDTree(gl[['centroid_x', 'centroid_y']].to_numpy(float))
        glom[sample == m] = tree.query(xy[sample == m])[0]
    label = pt.segment_class.astype(str).map(fs.PT).to_numpy()
    distal = np.asarray(counts[:, [qc.var_names.get_loc(g) for g in DISTAL if g in qc.var_names]].sum(1)).ravel()
    axes = {(k, h): pd.read_parquet(root / 'paper_figures/fig3/robustness/hvm/reruns' /
                                    f'gene_half_s{k}_{h}.parquet')['t'].reindex(units).to_numpy(float)
            for k in range(fs.N_SPLITS) for h in 'AB'}
    called = genes.index[genes.catR | genes.catH]
    grows = []
    for g in called:
        r = genes.loc[g]
        pat = r.pattern
        pos = np.stack([axes[(k, pat[k])] for k in range(fs.N_SPLITS)], 1)
        okp = np.isfinite(pos).sum(1) >= fs.MIN_AXES
        mean = np.where(okp, np.nanmean(np.where(np.isfinite(pos), pos, np.nan), 1), np.nan)
        x = f5.mouse_quantile(mean, sample, okp)
        y = np.asarray(counts[:, qc.var_names.get_loc(g)].todense()).ravel()
        row = {'gene': g, 'category': 'R' if r.catR else 'H', **{f'E_{m}': r[f'E_{m}'] for m in f5.SECTIONS},
               **{f'R_{m}': r[f'R_{m}'] for m in f5.SECTIONS},
               'neighbour_ratio_mouse': r.neighbour_ratio_mouse, 'neighbour_ratio_human': r.neighbour_ratio_human}
        for m in f5.SECTIONS:
            k = okp & (sample == m)
            bins = np.array_split(np.argsort(x[k]), fs.N_BINS)
            gr = [y[k][b].sum() / lib[k][b].sum() for b in bins]
            dr = [distal[k][b].sum() / lib[k][b].sum() for b in bins]
            row[f'admixture_rho_{m}'] = float(spearmanr(gr, dr).statistic)
        for m in f5.MOUSE:
            k = okp & (sample == m)
            b12, b23 = r.b12, r.b23
            s1 = k & (label == 0) & (x <= b12 - fs.MARGIN)
            s3 = k & (label == 2) & (x >= b23 + fs.MARGIN)
            cort = s3 & (glom <= np.nanpercentile(glom[s1], 90))
            rate = lambda u: np.log2(1e4 * (y[u].sum() + 0.5) / max(lib[u].sum(), 1.0))  # noqa: E731
            row[f'R_cortical_{m}'] = rate(cort) - rate(s1)
            row[f'cortical_share_{m}'] = cort.sum() / max(s3.sum(), 1)
            row[f'cortical_same_sign_{m}'] = bool(np.sign(row[f'R_cortical_{m}']) == np.sign(r[f'E_{m}']))
        for m in f5.HUMAN:
            row[f'lowmito_same_sign_{m}'] = bool(np.sign(r[f'R_lowmito_{m}']) == np.sign(r[f'E_{m}']))
            row[f'R_lowmito_{m}'] = r[f'R_lowmito_{m}']
        grows.append(row)
    gdf = pd.DataFrame(grows)
    # distal admixture fraction of the S1 and S3 cores, per section and neighbour (on the reference x)
    nb_all = ref[ref['sample'].astype(str).isin(f5.SECTIONS) & ref.segment_class.astype(str).isin(fs.NEIGHBOURS)]
    nbc = sparse.csr_matrix(qc[nb_all.index.to_numpy()].layers['counts'])
    nblib = np.asarray(nbc[:, ~mt].sum(1)).ravel()
    nbs, nbg = nb_all['sample'].astype(str).to_numpy(), nb_all.segment_class.astype(str).map(fs.NEIGHBOURS).to_numpy()
    ms = np.isin(sample, f5.MOUSE)
    ref_t = pt['pseudospace'].to_numpy(float)
    c = np.sort(ref_t[ms])
    xref = np.interp(ref_t, c, (np.arange(len(c)) + 0.5) / len(c))
    b12m, b23m = genes.b12.median(), genes.b23.median()
    frows = []
    for m in f5.SECTIONS:
        for core, (lo, hi) in (('S1', (0, b12m - fs.MARGIN)), ('S3', (b23m + fs.MARGIN, 1))):
            k = (sample == m) & (xref >= lo) & (xref <= hi)
            row = {'section': m, 'core': core, 'n': int(k.sum())}
            for grp, mk in NEIGHBOUR_MARKERS.items():
                cols = [qc.var_names.get_loc(z) for z in mk]
                kn = (nbs == m) & (nbg == grp)
                row[f'f_{grp}'] = (counts[k][:, cols].sum() / lib[k].sum()) / (nbc[kn][:, cols].sum() / nblib[kn].sum())
            frows.append(row)
    fr = pd.DataFrame(frows)
    fr.to_csv(screen / 'interpret_admixture_fraction.csv', index=False)
    fmax = fr.set_index(['section', 'core'])[[f'f_{g}' for g in NEIGHBOUR_MARKERS]].max(axis=1)
    if len(gdf):
        def bound(r, m):
            ratio = r.neighbour_ratio_mouse if m in f5.MOUSE else r.neighbour_ratio_human
            return np.log2(1 + fmax[(m, 'S3')] * ratio) - np.log2(1 + fmax[(m, 'S1')] * ratio)
        for m in f5.SECTIONS:
            gdf[f'admixture_bound_{m}'] = [bound(r, m) for _, r in gdf.iterrows()]
        # flagged when admixture could produce the gene's rise in some section (falls cannot come from admixture)
        gdf['admixture_flag'] = np.any([(gdf[f'E_{m}'] > 0) & (gdf[f'admixture_bound_{m}'] >= gdf[f'E_{m}'].abs())
                                        for m in f5.SECTIONS], 0)
        gdf['admixture_bound_max'] = gdf[[f'admixture_bound_{m}' for m in f5.SECTIONS]].max(axis=1)
        gdf['lowmito_ok'] = np.all([gdf[f'lowmito_same_sign_{m}'] for m in f5.HUMAN], 0)
        gdf['cortical_ok'] = np.all([gdf[f'cortical_same_sign_{m}'] for m in f5.MOUSE], 0)
    gdf.to_csv(screen / 'interpret_genes.csv', index=False)
    print(fr.round(4).to_string(index=False))

    # ---------------------------------------------------------------------------------------------- pathways
    pw = pd.read_csv(screen / 'pathways.csv')
    calls = pw[(pw.split == 'species') & pw.call].copy()
    if not len(calls):     # nothing called: read the ten best-ranked sets instead (descriptive, marked uncalled)
        calls = pw[pw.split == 'species'].nsmallest(10, 'p').copy()
    sets = json.loads((root / 'pt_package_aki/tables/gene_sets.json').read_text())
    E = {m: genes[f'E_{m}'].to_numpy() for m in f5.SECTIONS}
    R = {m: genes[f'R_{m}'].to_numpy() for m in f5.SECTIONS}
    Rlow = {**{m: R[m] for m in f5.MOUSE}, **{m: genes[f'R_lowmito_{m}'].to_numpy() for m in f5.HUMAN}}
    # cortical mouse raw R for every assessable gene (post hoc), one pass per axis
    R_cort = {m: np.full(len(genes), np.nan) for m in f5.MOUSE}
    gpos = {g: i for i, g in enumerate(genes.index)}
    for pat, grp in genes[ok].groupby('pattern'):
        pos = np.stack([axes[(k, pat[k])] for k in range(fs.N_SPLITS)], 1)
        okp = np.isfinite(pos).sum(1) >= fs.MIN_AXES
        x = f5.mouse_quantile(np.where(okp, np.nanmean(np.where(np.isfinite(pos), pos, np.nan), 1), np.nan),
                              sample, okp)
        b12, b23 = grp.b12.iloc[0], grp.b23.iloc[0]
        sub = counts[:, [qc.var_names.get_loc(g) for g in grp.index]]
        for m in f5.MOUSE:
            k = okp & (sample == m)
            s1 = k & (label == 0) & (x <= b12 - fs.MARGIN)
            cort = k & (label == 2) & (x >= b23 + fs.MARGIN) & (glom <= np.nanpercentile(glom[s1], 90))
            r1, r3 = (np.log2(1e4 * (np.asarray(sub[u].sum(0)).ravel() + 0.5) / lib[u].sum()) for u in (s1, cort))
            R_cort[m][[gpos[g] for g in grp.index]] = r3 - r1
    Rcort = {**R_cort, **{m: R[m] for m in f5.HUMAN}}
    d = pd.Series(genes['d_species'].to_numpy(), index=genes.index)
    d_raw = pd.Series(residual(R, ok), index=genes.index)
    d_low = pd.Series(residual(Rlow, ok), index=genes.index)
    d_cort = pd.Series(residual(Rcort, ok), index=genes.index)
    pool = d[ok].to_numpy()
    rows = []
    for _, c in calls.iterrows():
        mem = [g for g in sets[c.set] if g in genes.index and ok[genes.index.get_loc(g)]]
        w = d.loc[mem]
        sign = np.sign(c.d)
        order = (w * sign).sort_values(ascending=False)
        rest = order.index[3:]
        null = np.array([pool[rng.choice(len(pool), len(rest), replace=False)].mean() for _ in range(N_NULL)])
        d3 = w.loc[rest].mean()
        p3 = (1 + (np.abs(null - null.mean()) >= abs(d3 - null.mean())).sum()) / (1 + N_NULL)
        rows.append({'set': c.set, 'called': bool(c.call), 'n': len(mem), 'd': c.d, 'p': c.p, 'q': c.q, 'raw_d': d_raw.loc[mem].mean(),
                     'raw_same_sign': bool(np.sign(d_raw.loc[mem].mean()) == sign),
                     'lowmito_d': d_low.loc[mem].mean(), 'lowmito_same_sign': bool(np.sign(d_low.loc[mem].mean()) == sign),
                     'cortical_d': d_cort.loc[mem].mean(),
                     'cortical_same_sign': bool(np.sign(d_cort.loc[mem].mean()) == sign),
                     'drop3_d': d3, 'drop3_p': p3,
                     'top3_share': float(order.iloc[:3].sum() / max(order[order > 0].sum(), 1e-9)),
                     'mouse_E': np.mean([genes.loc[mem, f'E_{m}'].mean() for m in f5.MOUSE]),
                     'human_E': np.mean([genes.loc[mem, f'E_{m}'].mean() for m in f5.HUMAN]),
                     'leading_genes': ', '.join(f'{g} {w[g]:+.2f}' for g in order.index[:6])})
    out = pd.DataFrame(rows)
    if len(out):
        out['robust'] = out.raw_same_sign & out.lowmito_same_sign & (out.drop3_p <= 0.05)
    out.to_csv(screen / 'interpret_pathways.csv', index=False)

    # ---------------------------------------------------------------------------------------------- handoffs
    hp = pd.read_csv(screen / 'handoffs.csv')
    sel = hp[hp.shifted | hp.mouse_only | hp.human_only]
    cols = ['early', 'late', 'defined_in', 'shift', 'shift_lo', 'shift_hi', 'beyond_warp', 'mouse_only', 'human_only',
            'quantile_shift', 'neighbour_ok_early', 'neighbour_ok_late'] + \
        [f'x_{m}' for m in f5.SECTIONS] + [f'E_{r}_{m}' for r in ('early', 'late') for m in f5.SECTIONS]
    sel[cols].to_csv(screen / 'interpret_handoffs.csv', index=False)

    pd.set_option('display.width', 250)
    pd.set_option('display.max_colwidth', 90)
    print(gdf.round(2).to_string(index=False))
    print(out.round(3).to_string(index=False))
    print(sel[cols].round(2).to_string(index=False))


if __name__ == '__main__':
    main()
