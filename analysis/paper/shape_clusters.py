"""Shape-difference screen for Fig. 5 (mouse vs human) and Fig. 6 (control vs AKI): genes, then pathways.

Each section's curve (log2(CP10k + 0.5) over pseudostructure 0.05-0.95) is centred and scaled to unit variance, so
only its shape is compared, never its level or amplitude (the two species are measured on different probe panels).
A difference curve is the mean shape of one condition's two sections minus the other's.
 1  eligible: the gene changes at least 1.4-fold along the PT in both sections of at least one condition, so that it
    has a shape somewhere (a z-scored flat curve is noise)
 2  genes: a gene changes shape when its difference curve has RMS >= TAU (z units; 2 is a full reversal, sqrt(2)
    unrelated shapes) and both balanced section swaps give at most RHO times that. The rule applied with each swap as
    the contrast gives the null count. Changed genes' difference curves, rescaled to unit variance, are clustered by
    Ward linkage into K shape groups.
 3  pathways: a Reactome / KEGG / Hallmark set's difference curve is the mean of its eligible members' curves, so it
    is large only when members change the same way. A set changes when that RMS >= SET_TAU and both swaps give at
    most SET_RHO times it (self-contained; the swaps again give the null). Each set also gets a competitive p against
    N_DRAW size-matched random eligible genes (BH). Changed sets are clustered by their difference curves (Ward,
    SET_K groups).
Writes genes.csv, pathways.csv, shapes.npz and summary.json to <results-root>/paper_figures/<fig>/shape_clusters/.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from statsmodels.stats.multitest import multipletests

REPO = next(p for p in Path(__file__).resolve().parents if (p / 'pseudospace').is_dir())
RUNS = {'hvm': ('fig5', ['Ctrl1A2', 'Ctrl1A4'], ['HUK1_COR1', 'HUK1_MED1']),
        'aki': ('fig6', ['Ctrl1A2', 'Ctrl1A4'], ['IR2A2', 'IR2A4'])}
TRIM, AMP_LOG2 = 0.05, 0.5
TAU, RHO, K = 0.5, 0.5, 4
SET_MIN, SET_MAX, SET_TAU, SET_RHO, SET_K, N_DRAW = 10, 300, 0.5, 0.5, 4, 10000


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--run', choices=RUNS, required=True)
    p.add_argument('--results-root', default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT', str(REPO / 'results')))
    return p.parse_args(argv)


def zshape(c):
    c = c - c.mean(-1, keepdims=True)
    return c / c.std(-1, keepdims=True)


def partitions(a, b):
    """The condition split and the two balanced section swaps, each as (group 1, group 2)."""
    return {'condition': (a, b), 'swap 1': ([a[0], b[1]], [a[1], b[0]]), 'swap 2': ([a[0], b[0]], [a[1], b[1]])}


def rms(v):
    return np.sqrt((v ** 2).mean(-1))


def specific(r, tau, rho):
    """r: {partition: RMS}. A partition passes when its RMS >= tau and every other partition's is <= rho times it."""
    return {k: (r[k] >= tau) & np.all([r[o] <= rho * r[k] for o in r if o != k], 0) for k in r}


def ward(curves, k):
    """Ward clusters of unit-variance curves, numbered by where the cluster's mean difference peaks."""
    u = zshape(curves)
    lab = fcluster(linkage(u, 'ward'), k, 'maxclust')
    peak = {c: np.argmax(np.abs(u[lab == c].mean(0))) for c in np.unique(lab)}
    order = {c: i + 1 for i, c in enumerate(sorted(peak, key=peak.get))}
    return np.array([order[c] for c in lab])


def main(argv=None):
    args = parse_args(argv)
    rroot = Path(args.results_root)
    fig, a, b = RUNS[args.run]
    out = rroot / 'paper_figures' / fig / 'shape_clusters'
    out.mkdir(parents=True, exist_ok=True)
    screen = rroot / 'paper_figures' / fig / 'screen'
    zf = np.load(screen / 'all_curves.npz', allow_pickle=True)
    grid = zf['grid']
    w = (grid >= TRIM) & (grid <= 1 - TRIM)
    info = pd.read_csv(screen / 'genes.csv').set_index('gene')
    keep = info.screened if 'screened' in info else info.universe & info.neighbour_ok
    names = np.asarray(zf['genes']).astype(str)
    sel = keep.reindex(names).fillna(False).to_numpy(bool)
    logc = {m: np.log2(2 ** zf[f'curve|{m}'][sel][:, w].astype(float) + 0.5) for m in a + b}
    graded = {k: np.all([np.ptp(logc[m], 1) >= AMP_LOG2 for m in g], 0) for k, g in (('A', a), ('B', b))}
    elig = graded['A'] | graded['B']
    genes = names[sel][elig]
    z = {m: zshape(v[elig]) for m, v in logc.items()}
    parts = partitions(a, b)
    d = {k: np.mean([z[m] for m in g2], 0) - np.mean([z[m] for m in g1], 0) for k, (g1, g2) in parts.items()}

    # ------------------------------------------------------------------------------------------------------ genes
    r = {k: rms(v) for k, v in d.items()}
    call = specific(r, TAU, RHO)
    hit = call['condition']
    lab = np.zeros(len(genes), int)
    lab[hit] = ward(d['condition'][hit], K)
    mean_a, mean_b = (np.mean([z[m] for m in g], 0) for g in (a, b))
    tab = pd.DataFrame({'gene': genes, 'graded_A': graded['A'][elig], 'graded_B': graded['B'][elig],
                        **{f'rms_{k}': v for k, v in r.items()}, 'changed': hit, 'cluster': lab,
                        'r_conditions': (zshape(mean_a) * zshape(mean_b)).mean(1)})
    tab.to_csv(out / 'genes.csv', index=False)

    # --------------------------------------------------------------------------------------------------- pathways
    sets = json.loads((rroot / 'pt_package_aki/tables/gene_sets.json').read_text())
    gi = {g: i for i, g in enumerate(genes)}
    mem = {s: [gi[g] for g in v if g in gi] for s, v in sets.items()}
    mem = {s: i for s, i in mem.items() if SET_MIN <= len(i) <= SET_MAX}
    rng = np.random.default_rng(0)
    null = {}
    rows = []
    for s, idx in mem.items():
        n = len(idx)
        if n not in null:
            draws = np.array([rng.choice(len(genes), n, replace=False) for _ in range(N_DRAW)])
            null[n] = rms(d['condition'][draws].mean(1))
        row = {'set': s, 'n': n, **{f'rms_{k}': float(rms(v[idx].mean(0))) for k, v in d.items()}}
        row['p'] = (1 + (null[n] >= row['rms_condition']).sum()) / (N_DRAW + 1)
        row['genes'] = ','.join(genes[idx])
        rows.append(row)
    pw = pd.DataFrame(rows)
    pw['q'] = multipletests(pw.p, method='fdr_bh')[1]
    pcall = specific({k: pw[f'rms_{k}'].to_numpy() for k in d}, SET_TAU, SET_RHO)
    pw['changed'] = pcall['condition']
    sets_hit = pw.set[pw.changed].tolist()
    curve = {s: (mean_a[mem[s]].mean(0), mean_b[mem[s]].mean(0)) for s in sets_hit}
    pw['group'] = 0
    if len(sets_hit) >= SET_K:
        pw.loc[pw.changed, 'group'] = ward(np.array([curve[s][1] - curve[s][0] for s in sets_hit]), SET_K)
    pw = pw.sort_values(['changed', 'rms_condition'], ascending=False)
    pw.to_csv(out / 'pathways.csv', index=False)
    np.savez(out / 'shapes.npz', x=grid[w], genes=genes, mean_a=mean_a, mean_b=mean_b,
             **{f'z|{m}': z[m] for m in a + b})
    summary = {'run': args.run, 'conditions': [a, b], 'amp_log2': AMP_LOG2, 'tau': TAU, 'rho': RHO, 'k': K,
               'set_tau': SET_TAU, 'set_rho': SET_RHO, 'set_k': SET_K, 'n_draw': N_DRAW,
               'n_screened': int(sel.sum()), 'n_eligible': int(elig.sum()),
               'genes_changed': {k: int(v.sum()) for k, v in call.items()},
               'gene_clusters': {int(c): int((lab == c).sum()) for c in range(1, K + 1)},
               'n_sets_tested': len(pw), 'sets_changed': {k: int(v.sum()) for k, v in pcall.items()},
               'sets_changed_and_q05': int((pw.changed & (pw.q <= 0.05)).sum())}
    (out / 'summary.json').write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))
    show = pw[pw.changed].assign(set=lambda x: x.set.str.slice(0, 72))
    print(show.groupby('group').head(8)[['group', 'set', 'n', 'rms_condition', 'rms_swap 1', 'rms_swap 2', 'q']]
          .sort_values(['group', 'rms_condition'], ascending=[True, False]).round(3).to_string(index=False))


if __name__ == '__main__':
    main()
