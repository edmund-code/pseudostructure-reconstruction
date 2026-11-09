"""Fig. 4 screen: genes whose pattern along PT the S1/S2/S3 means and the bulk mean cannot express.

Control PT units of the control/AKI fit. Each gene is placed on the mean of the five Fig. 3d gene-half axes built
without it, and fitted in each control mouse by NB GLMs: bulk, clusters (S1/S2/S3), clusters + curve, curve.
Category 1 is a peak between segments, category 2 a change within a segment; both need the curve to beat the
clusters in both mice at a decoy-calibrated 5% threshold. Protocol (frozen before this ran):
docs/paper/fig4-continuous-vs-discrete.md.

  python analysis/paper/fig4_screen.py [--results-root <dir>] [--jobs 8]

Writes <results-root>/paper_figures/fig4/screen/: genes.csv (one row per screened gene), curves.npz (curves,
cluster means, bulk and bin means of every screened gene in each mouse), axes.csv (boundaries of each axis).
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

MICE = ['Ctrl1A2', 'Ctrl1A4']
PT = {'PT-S1': 0, 'PT-S2': 1, 'PT-S3': 2}
NEIGHBOURS = {'mTAL': 'AL', 'cTAL': 'AL', 'Macula-densa': 'AL', 'DCT1': 'DCT', 'DCT2': 'DCT',
              'CCD': 'CD', 'OMCD': 'CD', 'IMCD': 'CD'}
N_SPLITS, MIN_AXES = 5, 3
MIN_DETECTED = 0.10        # share of PT units with a count, in each mouse
MARGIN = 0.06              # cores lie further than this from a boundary; the transition window lies within it
PEAK_LOG2, CHANGE_LOG2 = 0.5, 1.0
FDR = 0.05
NEIGHBOUR_RATIO = 2.0
BASIS_DF, GRID, N_BINS = 6, np.linspace(0, 1, 201), 20


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--results-root', default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT', str(REPO / 'results')))
    p.add_argument('--out', default='paper_figures/fig4/screen', help='relative to the results root')
    p.add_argument('--jobs', type=int, default=8)
    return p.parse_args(argv)


def gene_halves(var: pd.DataFrame) -> list[set]:
    """Half A of each Fig. 3d split, drawn by the robustness script's own code."""
    import anndata as ad
    from fig3_robustness import EXCLUDE_COL, perturb

    out = []
    for k in range(N_SPLITS):
        tmp, info = perturb(ad.AnnData(var=var.copy()), {'condition': 'gene_half', 'fraction': 0.5, 'seed': k,
                                                         'half': 'A'}, 1e4)
        out.append(set(tmp.var_names[~tmp.var[EXCLUDE_COL].to_numpy(bool)]))
    return out


def pattern_axis(axes: dict, pat: str):
    """Rank-scaled mean of the half-axes that exclude a gene with this pattern; ok = units on >= MIN_AXES axes."""
    pos = np.stack([axes[(k, pat[k])] for k in range(N_SPLITS)], 1)
    ok = np.isfinite(pos).sum(1) >= MIN_AXES
    mean = np.where(ok, np.nanmean(np.where(np.isfinite(pos), pos, np.nan), 1), np.nan)
    x = np.full(len(pos), np.nan)
    r = pd.Series(mean[ok]).rank(method='average').to_numpy()
    x[ok] = (r - 0.5) / ok.sum()
    return x, ok


def boundary(x: np.ndarray, label: np.ndarray, lo: int) -> float:
    """Threshold on x misclassifying the fewest units between labels lo and lo + 1."""
    k = (label == lo) | (label == lo + 1)
    xs, hi = x[k], label[k] == lo + 1
    order = np.argsort(xs)
    xs, hi = xs[order], hi[order]
    # errors at a cut after position i: lo-units above it + hi-units at or below it
    err = (~hi[::-1]).cumsum()[::-1] - ~hi + hi.cumsum()
    i = int(np.argmin(err))
    return float((xs[i] + xs[min(i + 1, len(xs) - 1)]) / 2)


def shuffle_within(x: np.ndarray, label: np.ndarray, rng) -> np.ndarray:
    out = x.copy()
    for g in np.unique(label):
        idx = np.flatnonzero(label == g)
        out[idx] = x[rng.permutation(idx)]
    return out


def log2_cp10k(eta: np.ndarray) -> np.ndarray:
    return (eta + np.log(1e4)) / np.log(2)


def fit_mouse(y, lib, x, label, knots, rng, jobs):
    """NB fits of one mouse on one axis. y: units x genes (dense)."""
    from pseudospace.nb_gam import nb_dispersion, nb_fit, nested_nb_lr
    from pseudospace.stats_gam import make_gam_design

    w, off = np.ones(len(x)), np.log(lib)
    onehot = np.eye(3)[label]
    spline = make_gam_design(x, knots)
    designs = {'clusters': onehot, 'full': np.column_stack([onehot, spline[:, 1:]])}
    disp = nb_dispersion(y, designs['full'], w, off, n_jobs=jobs)
    lr = nested_nb_lr(y, designs, w, off, disp, {'beyond': ('clusters', 'full')}, n_jobs=jobs, fill=None)

    xs = shuffle_within(x, label, rng)
    decoy = {'clusters': onehot, 'full': np.column_stack([onehot, make_gam_design(xs, knots)[:, 1:]])}
    ddisp = nb_dispersion(y, decoy['full'], w, off, n_jobs=jobs)
    dlr = nested_nb_lr(y, decoy, w, off, ddisp, {'beyond': ('clusters', 'full')}, n_jobs=jobs, fill=None)

    curve = nb_fit(y, spline, w, off, disp, n_jobs=jobs)
    clusters = nb_fit(y, onehot, w, off, disp, n_jobs=jobs)
    bulk = nb_fit(y, np.ones((len(x), 1)), w, off, disp, n_jobs=jobs)
    order = np.argsort(x)
    bins = np.array_split(order, N_BINS)
    with np.errstate(divide='ignore'):
        bin_means = np.stack([np.log2(1e4 * y[b].sum(0) / lib[b].sum()) for b in bins], 1)
    return {'lr': lr['beyond'], 'decoy_lr': dlr['beyond'],
            'converged': lr['converged'] & curve['converged'] & clusters['converged'],
            'curve': log2_cp10k(curve['beta'] @ make_gam_design(GRID, knots).T),
            'clusters': log2_cp10k(clusters['beta']), 'bulk': log2_cp10k(bulk['beta'][:, 0]),
            'bins': bin_means, 'bin_x': np.array([x[b].mean() for b in bins])}


def decoy_threshold(real: np.ndarray, decoy: np.ndarray, fdr: float) -> float:
    """Smallest c with #(decoy >= c) / #(real >= c) <= fdr."""
    real, decoy = real[np.isfinite(real)], np.sort(decoy[np.isfinite(decoy)])
    for c in np.sort(real):
        n_real = (real >= c).sum()
        n_decoy = len(decoy) - np.searchsorted(decoy, c, 'left')
        if n_decoy / n_real <= fdr:
            return float(c)
    return np.inf


def shape_calls(curve: np.ndarray, b12: float, b23: float) -> dict:
    """Category 1 and 2 shape criteria of one curve (log2 on GRID)."""
    cores = [(0.0, b12 - MARGIN), (b12 + MARGIN, b23 - MARGIN), (b23 + MARGIN, 1.0)]
    core_mean = [curve[(GRID >= lo) & (GRID <= hi)].mean() for lo, hi in cores]
    i = int(np.argmax(curve))
    xp = GRID[i]
    peak_at = 12 if abs(xp - b12) <= MARGIN else 23 if abs(xp - b23) <= MARGIN else 0
    flank = {12: (0, 1), 23: (1, 2)}.get(peak_at)
    height = curve[i] - max(core_mean[flank[0]], core_mean[flank[1]]) if flank else np.nan
    changes = []
    for lo, hi in cores:
        third = (hi - lo) / 3
        first = curve[(GRID >= lo) & (GRID <= lo + third)].mean()
        last = curve[(GRID >= hi - third) & (GRID <= hi)].mean()
        changes.append(last - first)
    seg = int(np.argmax(np.abs(changes)))
    return {'peak_x': xp, 'peak_at': peak_at, 'peak_height': height, 'change_segment': seg,
            'change': changes[seg], 'changes': changes, 'core_mean': core_mean}


def main(argv=None):
    import anndata as ad
    from scipy import sparse
    from pseudospace.stats_gam import gam_internal_knots

    args = parse_args(argv)
    root = Path(args.results_root)
    out = root / args.out
    out.mkdir(parents=True, exist_ok=True)

    qc = ad.read_h5ad(root / RUNS['aki'] / '02_qc/units_qc.h5ad')
    ref = pd.read_parquet(root / RUNS['aki'] / '04_pseudospace/pseudospace.parquet')
    ref = ref.set_index('obs_name') if 'obs_name' in ref.columns else ref
    ctrl = ref[ref['sample'].astype(str).isin(MICE)]
    pt = ctrl[ctrl.segment_class.astype(str).isin(PT)]
    units = pt.index.to_numpy()

    counts = sparse.csr_matrix(qc[units].layers['counts'])
    mt = qc.var['mt'].to_numpy(bool)
    lib = np.asarray(counts[:, ~mt].sum(1)).ravel()
    sample = pt['sample'].astype(str).to_numpy()
    label = pt.segment_class.astype(str).map(PT).to_numpy()

    detected = np.ones(qc.n_vars, bool)
    for m in MICE:
        k = sample == m
        detected &= np.asarray((counts[k] > 0).mean(0)).ravel() >= MIN_DETECTED
    keep = detected & ~mt
    genes = qc.var_names[keep]

    # neighbour filter: counts per 10k in AL, DCT and CD units of the control sections against PT
    nb = ctrl[ctrl.segment_class.astype(str).isin(NEIGHBOURS)]
    nbc = sparse.csr_matrix(qc[nb.index.to_numpy()].layers['counts'])[:, keep]
    nblib = np.asarray(sparse.csr_matrix(qc[nb.index.to_numpy()].layers['counts'])[:, ~mt].sum(1)).ravel()
    pt_rate = np.asarray(counts[:, keep].sum(0)).ravel() / lib.sum()
    group = nb.segment_class.astype(str).map(NEIGHBOURS).to_numpy()
    nb_rate = {g: np.asarray(nbc[group == g].sum(0)).ravel() / nblib[group == g].sum() for g in ('AL', 'DCT', 'CD')}
    neighbour_ratio = np.max(np.stack(list(nb_rate.values())), 0) / np.maximum(pt_rate, 1e-12)

    # axes: each gene gets the mean of the five half-axes built without it
    halves = gene_halves(qc.var)
    axes = {}
    for k in range(N_SPLITS):
        for h in 'AB':
            t = pd.read_parquet(root / 'paper_figures/fig3/robustness/aki/reruns' / f'gene_half_s{k}_{h}.parquet')['t']
            axes[(k, h)] = t.reindex(units).to_numpy(float)
    pattern = np.array([''.join('B' if g in halves[k] else 'A' for k in range(N_SPLITS)) for g in genes])
    # pattern letter = the half whose axis excludes the gene

    rows, store = [], {m: {} for m in MICE}
    axis_rows = []
    rng = np.random.default_rng(0)
    for pi, pat in enumerate(sorted(set(pattern))):
        x, ok = pattern_axis(axes, pat)
        b12, b23 = boundary(x[ok], label[ok], 0), boundary(x[ok], label[ok], 1)
        knots = gam_internal_knots(x[ok], BASIS_DF)
        axis_rows.append({'pattern': pat, 'n_units': int(ok.sum()), 'b12': b12, 'b23': b23})
        gi = np.flatnonzero(pattern == pat)
        sub = counts[:, np.flatnonzero(keep)[gi]]
        fits = {}
        for m in MICE:
            k = ok & (sample == m)
            fits[m] = fit_mouse(sub[k].toarray(), lib[k], x[k], label[k], knots, rng, args.jobs)
        for j, g in enumerate(genes[gi]):
            row = {'gene': g, 'pattern': pat, 'b12': b12, 'b23': b23,
                   'neighbour_ratio': float(neighbour_ratio[gi[j]])}
            for m in MICE:
                f = fits[m]
                s = shape_calls(f['curve'][j], b12, b23)
                row.update({f'lr_{m}': f['lr'][j], f'decoy_lr_{m}': f['decoy_lr'][j],
                            f'converged_{m}': bool(f['converged'][j]), f'peak_x_{m}': s['peak_x'],
                            f'peak_at_{m}': s['peak_at'], f'peak_height_{m}': s['peak_height'],
                            f'change_segment_{m}': s['change_segment'], f'change_{m}': s['change'],
                            f'bulk_{m}': f['bulk'][j]})
                row.update({f'cluster_S{i + 1}_{m}': f['clusters'][j, i] for i in range(3)})
                store[m][g] = (f['curve'][j], f['clusters'][j], f['bulk'][j], f['bins'][j], f['bin_x'])
            rows.append(row)
        print(f'axis {pi + 1}: {pat}  {len(gi)} genes  b12 {b12:.3f} b23 {b23:.3f}', flush=True)

    df = pd.DataFrame(rows)
    thresholds = {}
    for m in MICE:
        conv = df[f'converged_{m}']
        thresholds[m] = decoy_threshold(df.loc[conv, f'lr_{m}'].to_numpy(), df.loc[conv, f'decoy_lr_{m}'].to_numpy(), FDR)
        df[f'beyond_{m}'] = conv & (df[f'lr_{m}'] >= thresholds[m])
    a, b = MICE
    df['beyond_both'] = df[f'beyond_{a}'] & df[f'beyond_{b}']
    df['neighbour_ok'] = df['neighbour_ratio'] <= NEIGHBOUR_RATIO
    df['cat1'] = (df.beyond_both & df.neighbour_ok & (df[f'peak_at_{a}'] > 0) & (df[f'peak_at_{a}'] == df[f'peak_at_{b}'])
                  & (df[f'peak_height_{a}'] >= PEAK_LOG2) & (df[f'peak_height_{b}'] >= PEAK_LOG2))
    df['cat2'] = (df.beyond_both & df.neighbour_ok & (df[f'change_segment_{a}'] == df[f'change_segment_{b}'])
                  & (np.sign(df[f'change_{a}']) == np.sign(df[f'change_{b}']))
                  & (df[f'change_{a}'].abs() >= CHANGE_LOG2) & (df[f'change_{b}'].abs() >= CHANGE_LOG2))
    df['cat1_score'] = df[[f'peak_height_{a}', f'peak_height_{b}']].min(1)
    df['cat2_score'] = df[[f'change_{a}', f'change_{b}']].abs().min(1)
    df.to_csv(out / 'genes.csv', index=False)
    pd.DataFrame(axis_rows).to_csv(out / 'axes.csv', index=False)

    keep_genes = df.loc[df.cat1 | df.cat2 | df.beyond_both, 'gene'].tolist()
    np.savez_compressed(out / 'curves.npz', grid=GRID, genes=np.array(keep_genes), **{
        f'{part}|{m}': np.stack([store[m][g][i] for g in keep_genes]) if keep_genes else np.empty(0)
        for m in MICE for i, part in enumerate(('curve', 'clusters', 'bulk', 'bins', 'bin_x'))})
    summary = {'n_units': {m: int((sample == m).sum()) for m in MICE}, 'n_genes_screened': int(len(df)),
               'n_axes': len(axis_rows), 'thresholds': thresholds,
               'n_beyond': {m: int(df[f'beyond_{m}'].sum()) for m in MICE}, 'n_beyond_both': int(df.beyond_both.sum()),
               'n_neighbour_dropped': int((df.beyond_both & ~df.neighbour_ok).sum()),
               'n_cat1': int(df.cat1.sum()), 'n_cat2': int((df.cat2 & ~df.cat1).sum()),
               'n_cat1_by_boundary': df.loc[df.cat1, f'peak_at_{a}'].value_counts().to_dict(),
               'n_cat2_by_segment': (df.loc[df.cat2 & ~df.cat1, f'change_segment_{a}'] + 1).value_counts().to_dict(),
               'top_cat1': df[df.cat1].nlargest(10, 'cat1_score').gene.tolist(),
               'top_cat2': df[df.cat2 & ~df.cat1].nlargest(10, 'cat2_score').gene.tolist()}
    (out / 'summary.json').write_text(json.dumps(summary, indent=1, default=str))
    print(json.dumps(summary, indent=1, default=str))


if __name__ == '__main__':
    main()
