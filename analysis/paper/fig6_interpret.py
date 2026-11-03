"""Fig. 6 follow-up (post hoc, descriptive): how robust is each pathway call and each handoff shift, and what drives it?

For every pathway x segment call of fig6_screen.py:
 raw      the set's within-segment change measured on raw counts of units carrying that segment's label (mean over
          members of the screen's rule-4 change), its condition value and whether all four pairs share its sign
 trimmed  the set's within-segment change with the PT ends trimmed (S1 core from x = 0.05, S3 core to x = 0.95), so a
          spline's end behaviour cannot carry it
 drop3    the set's change with its three largest same-sign members removed, and its p against random sets of that
          size (does it rest on a few genes?)
 leading  the five members with the largest change in the call's direction, and their share of the set's total
Both: the within-segment change in control and in IR separately (does injury create a gradient, or erase one?).
For every single-gene within-segment call (category A; its list missed calibration, so these are checked examples, not a
list): the same raw, trimmed and per-condition readouts, and the neighbour ratios.
Gradient retention: genes with a within-segment change of at least 0.75 log2 in one control mouse, and how much of it
the other control and each IR mouse keep (median signed change / median selected change). Selecting in one control
and reading the other avoids regression to the mean. Split by injury-responsive genes (fig6_placement.py) against genes
whose whole-PT level barely changes (|condition| < 0.25 log2), which bounds what noisier IR positions could explain.
For every handoff pair that shifts: each partner's within-S3 and within-S2 raw-count change in control and in IR
(does the late partner lose its S3 rise, does the early partner stay on?).

  python analysis/paper/fig6_interpret.py
Reads <results-root>/paper_figures/fig6/screen/ and writes interpret_pathways.csv, interpret_handoffs.csv there.
Protocol and the reason for this follow-up: docs/paper/fig6-aki-along-pt.md.
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
sys.path.insert(0, str(Path(__file__).resolve().parent))
import fig6_screen as f6  # noqa: E402

END_TRIM = 0.05
N_NULL = 2000


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--results-root', default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT', str(REPO / 'results')))
    return p.parse_args(argv)


def trimmed_w(curves: dict, b12: np.ndarray, b23: np.ndarray, seg: int) -> np.ndarray:
    """Per-gene within-core change (condition contrast) with the PT ends trimmed; curves {mouse: (G, grid)}."""
    grid = f6.GRID
    out = {}
    for m, c in curves.items():
        w = np.empty(len(c))
        for i in range(len(c)):
            lo, hi = f6.cores(b12[i], b23[i])[seg]
            lo, hi = (max(lo, END_TRIM), hi) if seg == 0 else (lo, min(hi, 1 - END_TRIM)) if seg == 2 else (lo, hi)
            (a0, a1), (z0, z1) = f6.thirds(lo, hi)
            w[i] = c[i, (grid >= z0) & (grid <= z1)].mean() - c[i, (grid >= a0) & (grid <= a1)].mean()
        out[m] = w
    return f6.split_value(out, 'condition')


def main(argv=None):
    args = parse_args(argv)
    root = Path(args.results_root)
    screen = root / 'paper_figures/fig6/screen'
    genes = pd.read_csv(screen / 'genes.csv', index_col=0)
    pw = pd.read_csv(screen / 'pathways.csv')
    calls = pw[(pw.split == 'condition') & pw.call].copy()
    sets = json.loads((root / 'pt_package_aki/tables/gene_sets.json').read_text())
    z = np.load(screen / 'all_curves.npz', allow_pickle=True)
    names = list(z['genes'])
    pos = {g: i for i, g in enumerate(names)}
    screened = genes.universe & genes.neighbour_ok & genes.converged
    rng = np.random.default_rng(6)

    trimmed = {}
    for seg in (0, 2):
        idx = np.array([pos[g] for g in genes.index])
        trimmed[seg] = pd.Series(trimmed_w({m: z[f'curve|{m}'][idx] for m in f6.MICE}, genes.b12.to_numpy(),
                                           genes.b23.to_numpy(), seg), index=genes.index)

    rows = []
    for _, c in calls.iterrows():
        s = int(c.segment)
        ok = screened & ~genes[f'unsupported_S{s}']
        mem = [g for g in sets[c.set] if ok.get(g, False)]
        w = genes.loc[mem, f'W_condition_S{s}']
        sign = np.sign(c.W)
        lab = {m: genes.loc[mem, f'label_change_S{s}_{m}'].mean() for m in f6.MICE}
        raw = f6.split_value(lab, 'condition')
        raw_pairs = f6.split_pairs(lab, 'condition')
        trim = trimmed[s - 1].loc[mem].mean() if s in (1, 3) else w.mean()
        order = (w * sign).sort_values(ascending=False)
        rest = order.index[3:]
        pool = genes.loc[ok, f'W_condition_S{s}'].to_numpy()
        null = np.array([pool[rng.choice(len(pool), len(rest), replace=False)].mean() for _ in range(N_NULL)])
        d3 = w.loc[rest].mean()
        p3 = (1 + (np.abs(null - null.mean()) >= abs(d3 - null.mean())).sum()) / (1 + N_NULL)
        top = order.index[:5]
        rows.append({'set': c.set, 'segment': s, 'n': len(mem), 'W': c.W, 'q': c.q,
                     'min_abs_pair': np.abs([c[f'pair_{k}'] for k in range(4)]).min(),
                     'relabel_ratio': max(abs(c.other_0), abs(c.other_1)) / abs(c.W),
                     'raw_within_label': raw, 'raw_same_sign': bool(np.sign(raw) == sign),
                     'raw_pairs_agree': bool((np.sign(raw_pairs) == sign).all()),
                     'trimmed_W': trim, 'trimmed_same_sign': bool(np.sign(trim) == sign),
                     'drop3_W': d3, 'drop3_p': p3, 'top3_share': float(order.iloc[:3].sum() / max(order.sum(), 1e-9)),
                     'leading_genes': ', '.join(f'{g} {w[g]:+.2f}' for g in top),
                     'W_control': genes.loc[mem, [f'W_S{s}_{m}' for m in f6.CONTROLS]].mean(axis=1).mean(),
                     'W_IR': genes.loc[mem, [f'W_S{s}_{m}' for m in f6.INJURED]].mean(axis=1).mean(),
                     **{f'segment_view_S{k}': genes.loc[mem, f'segment_view_S{k}'].mean() for k in (1, 2, 3)}})
    out = pd.DataFrame(rows)
    out['robust'] = out.raw_same_sign & out.raw_pairs_agree & out.trimmed_same_sign & (out.drop3_p <= 0.05)
    out.to_csv(screen / 'interpret_pathways.csv', index=False)

    grows = []
    for s in (1, 2, 3):
        for g, r in genes[genes[f'catA_S{s}']].iterrows():
            lab = {m: r[f'label_change_S{s}_{m}'] for m in f6.MICE}
            grows.append({'gene': g, 'segment': s, 'W': r[f'W_condition_S{s}'], 'raw_within_label': r[f'W_label_condition_S{s}'],
                          'trimmed_W': trimmed[s - 1][g] if s in (1, 3) else r[f'W_condition_S{s}'],
                          'W_control': np.mean([r[f'W_S{s}_{m}'] for m in f6.CONTROLS]),
                          'W_IR': np.mean([r[f'W_S{s}_{m}'] for m in f6.INJURED]),
                          'raw_control': np.mean([lab[m] for m in f6.CONTROLS]), 'raw_IR': np.mean([lab[m] for m in f6.INJURED]),
                          'min_abs_pair': r[f'W_min_abs_pair_S{s}'],
                          'neighbour_ratio_max': max(r.neighbour_ratio_control, r.neighbour_ratio_IR),
                          **{f'segment_view_S{k}': r[f'segment_view_S{k}'] for k in (1, 2, 3)},
                          'bulk_control': np.mean([r[f'bulk_{m}'] for m in f6.CONTROLS]),
                          'bulk_IR': np.mean([r[f'bulk_{m}'] for m in f6.INJURED])})
    gdf = pd.DataFrame(grows)
    gdf['trimmed_same_sign'] = np.sign(gdf.trimmed_W) == np.sign(gdf.W)
    gdf.to_csv(screen / 'interpret_genes.csv', index=False)

    resp = pd.read_csv(root / 'paper_figures/fig6/placement/responsive_genes.csv', index_col=0).reindex(genes.index)
    groups = {'all': pd.Series(True, index=genes.index), 'injury-responsive': resp.responsive.fillna(False).astype(bool),
              'level unchanged': resp.condition.abs() < 0.25}
    rrows = []
    for sel, held in (('Ctrl1A2', 'Ctrl1A4'), ('Ctrl1A4', 'Ctrl1A2')):
        for s in (1, 2, 3):
            ok = screened & ~genes[f'unsupported_S{s}']
            w = genes.loc[ok, [f'W_S{s}_{m}' for m in f6.MICE]]
            w.columns = f6.MICE
            pick = w[w[sel].abs() >= 0.75]
            for name, mask in groups.items():
                p_ = pick[mask.reindex(pick.index).fillna(False).to_numpy(bool)]
                sign, den = np.sign(p_[sel]), p_[sel].abs().median()
                rrows.append({'selected_in': sel, 'segment': s, 'genes': name, 'n': len(p_),
                              **{f'kept_{m}': (p_[m] * sign).median() / den for m in f6.MICE if m != sel},
                              **{f'same_sign_{m}': (np.sign(p_[m]) == sign).mean() for m in f6.MICE if m != sel}})
    pd.DataFrame(rrows).to_csv(screen / 'interpret_retention.csv', index=False)

    hp = pd.read_csv(screen / 'handoffs.csv')
    hrows = []
    for _, r in hp[hp.shifted | hp.early.isin(['Slc7a9', 'Slc5a2', 'Slc5a12'])].iterrows():
        row = {'early': r.early, 'late': r.late, 'shift': r['shift'], 'shifted': r.shifted,
               'level_early': r.level_early, 'level_late': r.level_late}
        for role, g in (('early', r.early), ('late', r.late)):
            if g not in genes.index:
                continue
            for s in (2, 3):
                for cond, mice in (('control', f6.CONTROLS), ('IR', f6.INJURED)):
                    row[f'{role}_raw_S{s}_{cond}'] = np.mean([genes.loc[g, f'label_change_S{s}_{m}'] for m in mice])
        hrows.append(row)
    pd.DataFrame(hrows).to_csv(screen / 'interpret_handoffs.csv', index=False)

    pd.set_option('display.width', 250)
    pd.set_option('display.max_colwidth', 70)
    print(f'{len(out)} calls; robust on all four checks: {int(out.robust.sum())}')
    print(out.groupby('segment').robust.agg(['size', 'sum']).to_string())
    cols = ['set', 'segment', 'n', 'W', 'raw_within_label', 'raw_pairs_agree', 'trimmed_W', 'drop3_p', 'top3_share']
    print(out.sort_values(['robust', 'segment', 'W'], ascending=[False, True, True])[cols].round(3).to_string(index=False))
    print(pd.DataFrame(hrows).round(2).to_string(index=False))
    print(gdf.sort_values(['segment', 'W']).round(2).to_string(index=False))


if __name__ == '__main__':
    main()
