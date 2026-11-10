"""Question bank (post hoc, descriptive): per-section gradient readouts for the external question list.

PT questions read the Fig. 5 / Fig. 6 screen curves (gene-held-out NB curves, control-quantile x). Distal and
whole-nephron questions read binned profiles on each run's package pseudostructure t (not gene-held-out; see the
protocol for that ceiling). Each readout is contrasted between groups and against the two balanced relabelings.
Protocol (frozen before any readout): docs/paper/question-bank.md.

  python analysis/paper/question_bank.py [--results-root <dir>]
Post hoc additions (not in the frozen protocol; see the doc): injury_origin.csv classifies every gene induced in a PT
label (all four IR - control pairs >= 1 log2, label pseudobulk) as within the healthy PT range, at a level found in a
healthy distal segment, or above every healthy segment, with its healthy and IR S3 - S1 slopes; the same rule is run
on both relabelings.

Writes <results-root>/paper_figures/question_bank/: readouts.csv (one row per readout), profiles.npz,
injury_origin.csv, injury_origin_relabelings.json.
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
from pseudospace.handoffs import crossing, scaled_profile  # noqa: E402

GROUPS = {'hvm': (['HUK1_COR1', 'HUK1_MED1'], ['Ctrl1A2', 'Ctrl1A4']),
          'aki': (['IR2A2', 'IR2A4'], ['Ctrl1A2', 'Ctrl1A4'])}
N_BINS, MIN_UNITS = 60, 15
END_TRIM = 0.05
NEPHRON_GRID = (np.arange(N_BINS) + 0.5) / N_BINS
# label windows for distal questions: from the 5th percentile of t in the first label to the 95th in the last,
# pooled over the reference (mouse control) sections
WINDOWS = {'TAL': ['mTAL', 'cTAL'], 'DCT-CD': ['DCT1', 'DCT2', 'CNT', 'CCD'], 'nephron': None}

# (question id, kind, window, genes). kinds: peak / domain per gene; cross = first->second switch point;
# overlap of two genes; module_cross = mean scaled profile of the first list crossing that of the second
HVM = [
    ('H1 phosphate', 'pt', [('peak', 'Slc34a1'), ('peak', 'Slc34a3'), ('domain', 'Slc34a1'), ('domain', 'Slc34a3'),
                            ('cross', 'Slc34a1', 'Slc34a3')]),
    ('H2 drug secretion', 'pt', [(k, g) for g in ['Slc22a6', 'Slc22a8', 'Slc22a2', 'Slc47a1', 'Abcc2', 'Abcc4']
                                 for k in ('peak', 'domain')] + [('cross', 'Slc22a6', 'Slc22a2')]),
    ('H3 SGLT2 zone', 'pt', [('domain', 'Slc5a2'), ('domain', 'Slc5a1'), ('domain', 'Slc2a2'),
                             ('cross', 'Slc5a2', 'Slc5a1')]),
    ('H4 urate', 'pt', [(k, g) for g in ['Slc22a12', 'Pdzk1'] for k in ('peak', 'domain')]),
    ('H5 nutrient to secretion', 'pt', [('module_cross', ['Slc5a2', 'Slc34a1'],
                                         ['Slc22a6', 'Slc22a8', 'Slc22a2', 'Slc47a1'])]),
    ('H6 TAL transition', 'TAL', [(k, g) for g in ['Slc12a1', 'Umod', 'Kcnj1', 'Bsnd', 'Casr'] for k in ('peak', 'domain')]),
    ('H7 TAL claudins', 'TAL', [(k, g) for g in ['Cldn10', 'Cldn14', 'Cldn16', 'Cldn19'] for k in ('peak', 'domain')]
     + [('cross', 'Cldn10', 'Cldn16')]),
    ('H8 DCT1 to DCT2', 'DCT-CD', [(k, g) for g in ['Slc12a3', 'Pvalb', 'Trpm6', 'Calb1', 'Trpv5']
                                   for k in ('peak', 'domain')] + [('cross', 'Pvalb', 'Calb1')]),
    ('H9 NCC-ENaC overlap', 'DCT-CD', [('overlap', 'Slc12a3', 'Scnn1g'), ('overlap', 'Slc12a3', 'Scnn1b'),
                                       ('cross', 'Slc12a3', 'Scnn1g')]),
    ('H10 AQP2 onset', 'DCT-CD', [('domain', 'Aqp2'), ('domain', 'Avpr2'), ('domain', 'Hsd11b2'),
                                  ('cross', 'Slc12a3', 'Aqp2')]),
    ('H11 ENaC and K secretion', 'DCT-CD', [(k, g) for g in ['Scnn1a', 'Scnn1b', 'Scnn1g', 'Nr3c2', 'Kcnj1', 'Kcnj10']
                                            for k in ('peak', 'domain')]),
    ('H12 intercalated cells', 'DCT-CD', [(k, g) for g in ['Slc4a1', 'Slc26a4', 'Atp6v1b1', 'Foxi1']
                                          for k in ('peak', 'domain')] + [('cross', 'Slc4a1', 'Slc26a4')]),
    ('H14 boundary genes', 'nephron', [('domain', g) for g in ['Calb1', 'Pvalb', 'Aqp2', 'Slc8a1', 'Tmem213']]),
]
AKI_PT = ['Havcr1', 'Krt8', 'Krt18', 'Krt20', 'Sox9', 'Vim', 'Lcn2',                      # injury
          'Slc5a2', 'Slc34a1', 'Lrp2', 'Cubn', 'Slc22a6', 'Slc22a8', 'Slc22a2', 'Slc7a13',  # mature transport
          'Mki67', 'Pcna', 'Top2a', 'Acsm2', 'Vcam1', 'Il34', 'Icam1', 'Spp1', 'Clu',      # repair, failed repair
          'Ppara', 'Cpt1a', 'Acox1', 'Acadm']                                                # fatty-acid oxidation
AKI_DISTAL = ['Slc12a1', 'Umod', 'Kcnj1', 'Lcn2', 'Spp1', 'Clu', 'Krt8', 'Krt18', 'Sox9', 'Havcr1',
              'Slc12a3', 'Pvalb', 'Calb1', 'Trpv5', 'Scnn1a', 'Scnn1b', 'Scnn1g', 'Aqp2', 'Hsd11b2', 'Avpr2']
AKI_PAIRS = [('Pvalb', 'Calb1'), ('Slc12a3', 'Aqp2'), ('Slc12a3', 'Scnn1g')]


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--results-root', default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT', str(REPO / 'results')))
    return p.parse_args(argv)


# ------------------------------------------------------------------------------------------------- readouts
def peak(profile, grid):
    return float(grid[np.nanargmax(profile)]) if np.isfinite(profile).any() else np.nan


def domain(profile, grid):
    s = scaled_profile(profile)
    on = np.flatnonzero(s >= 0.5)
    return (float(grid[on[0]]), float(grid[on[-1]])) if len(on) else (np.nan, np.nan)


def overlap(a, b):
    sa, sb = scaled_profile(a), scaled_profile(b)
    ok = np.isfinite(sa) & np.isfinite(sb)
    return float(((sa >= 0.5) & (sb >= 0.5))[ok].mean()) if ok.any() else np.nan


def module(profiles):
    return np.nanmean([scaled_profile(p) for p in profiles], axis=0)


def contrast(q: dict, test, ref) -> dict:
    """Group contrast, the two balanced relabelings, the four cross pairs and the call rule of the protocol."""
    v = lambda a, b: np.mean([q[s] for s in a]) - np.mean([q[s] for s in b])  # noqa: E731
    c = v(test, ref)
    r1 = v([test[0], ref[0]], [test[1], ref[1]])
    r2 = v([test[0], ref[1]], [test[1], ref[0]])
    pairs = [q[t] - q[r] for t in test for r in ref]
    ok = np.isfinite(c) and all(np.isfinite(pairs))
    same = ok and (all(p > 0 for p in pairs) or all(p < 0 for p in pairs))
    return {'contrast': c, 'relabel_1': r1, 'relabel_2': r2, 'pairs_same_sign': bool(same),
            'call': bool(same and abs(c) > max(abs(r1), abs(r2)))}


# ------------------------------------------------------------------------------------------------- inputs
def pt_curves(root: Path, fig: str):
    z = np.load(root / f'paper_figures/{fig}/screen/all_curves.npz', allow_pickle=True)
    ax = pd.read_csv(root / f'paper_figures/{fig}/screen/axes.csv')
    genes = list(z['genes'])
    secs = [k.split('|')[1] for k in z.files if k.startswith('curve|')]
    # the screen's curves are log2 counts per 10k with no floor, so a gene absent in a section runs to -inf;
    # floor them as the nephron profiles are floored (+0.5 per 10k)
    get = lambda g, s: (np.log2(2 ** z[f'curve|{s}'][genes.index(g)].astype(float) + 0.5)  # noqa: E731
                        if g in genes else None)
    # the spline ends overshoot (Fig. 5 clips them); read only x in [END_TRIM, 1 - END_TRIM], as fig6_interpret.py
    k = (z['grid'] >= END_TRIM) & (z['grid'] <= 1 - END_TRIM)
    trimmed = lambda g, s: None if get(g, s) is None else get(g, s)[k]  # noqa: E731
    return z['grid'][k], trimmed, float(ax.b12.median()), float(ax.b23.median()), secs


def nephron_profiles(root: Path, run: str, genes):
    """Per section: 60 equal-width bins of t, log2 counts per 10k of pooled counts, 3-bin running mean."""
    import anndata as ad
    from scipy import sparse
    qc = ad.read_h5ad(root / RUNS[run] / '02_qc/units_qc.h5ad')
    ref = pd.read_parquet(root / RUNS[run] / '04_pseudospace/pseudospace.parquet')
    ref = ref.set_index('obs_name') if 'obs_name' in ref.columns else ref
    ref = ref[np.isfinite(ref['pseudospace'].astype(float))]
    units = ref.index.intersection(qc.obs_names)
    counts = sparse.csr_matrix(qc[units].layers['counts'])
    lib = np.asarray(counts[:, ~qc.var['mt'].to_numpy(bool)].sum(1)).ravel()
    genes = [g for g in genes if g in qc.var_names]
    sub = counts[:, [qc.var_names.get_loc(g) for g in genes]].toarray()
    t = ref.loc[units, 'pseudospace'].astype(float).to_numpy()
    sample = ref.loc[units, 'sample'].astype(str).to_numpy()
    label = ref.loc[units, 'segment_class'].astype(str).to_numpy()
    t = register(t, sample, label, GROUPS[run][1])
    b = np.clip((t * N_BINS).astype(int), 0, N_BINS - 1)
    prof = {}
    for s in np.unique(sample):
        k = sample == s
        n = np.bincount(b[k], minlength=N_BINS)
        L = np.bincount(b[k], weights=lib[k], minlength=N_BINS)
        for j, g in enumerate(genes):
            y = np.bincount(b[k], weights=sub[k, j], minlength=N_BINS)
            r = np.where(n >= MIN_UNITS, np.log2(1e4 * y / np.maximum(L, 1) + 0.5), np.nan)
            prof[(g, s)] = pd.Series(r).rolling(3, center=True, min_periods=2).mean().to_numpy()
    return prof, t, sample, label


def register(t, sample, label, ref_secs):
    """Map each section's t onto the reference sections' t by its own segment-label medians (piecewise linear,
    ends fixed at 0 and 1), so a gene's position is read relative to the section's morphological segments."""
    k = np.isin(sample, ref_secs)
    labels = sorted(np.unique(label[k]), key=lambda l: np.median(t[k & (label == l)]))
    target = [0.0] + [float(np.median(t[k & (label == l)])) for l in labels] + [1.0]
    out = t.copy()
    for s in np.unique(sample):
        m = sample == s
        own = [0.0] + [float(np.median(t[m & (label == l)])) if (m & (label == l)).sum() >= 20 else np.nan
                       for l in labels] + [1.0]
        ok = np.isfinite(own)
        out[m] = np.interp(t[m], np.array(own)[ok], np.array(target)[ok])
    return out


def window(t, sample, label, ref_secs, labels):
    if labels is None:
        return 0.0, 1.0
    k = np.isin(sample, ref_secs)
    lo = np.quantile(t[k & (label == next(l for l in labels if (label[k] == l).any()))], 0.05)
    hi = np.quantile(t[k & (label == [l for l in labels if (label[k] == l).any()][-1])], 0.95)
    return float(lo), float(hi)


# ------------------------------------------------------------------------------------------------- main
def measure(spec, get, grid, secs, test, ref, question, scope):
    """One row per readout of one question: per-section values plus contrast and call."""
    rows = []
    for item in spec:
        kind, args = item[0], item[1:]
        vals = {}
        for s in secs:
            if kind in ('peak', 'domain'):
                p = get(args[0], s)
                if p is None:
                    break
                vals[s] = peak(p, grid) if kind == 'peak' else domain(p, grid)
            elif kind == 'cross':
                a, b = get(args[0], s), get(args[1], s)
                if a is None or b is None:
                    break
                vals[s] = crossing(scaled_profile(a), scaled_profile(b), grid)
            elif kind == 'overlap':
                a, b = get(args[0], s), get(args[1], s)
                if a is None or b is None:
                    break
                vals[s] = overlap(a, b)
            elif kind == 'module_cross':
                ps = [[get(g, s) for g in lst] for lst in args]
                ps = [[p for p in lst if p is not None] for lst in ps]
                vals[s] = crossing(module(ps[0]), module(ps[1]), grid)
        if len(vals) < len(secs):
            rows.append({'question': question, 'scope': scope, 'readout': kind, 'genes': ' > '.join(map(str, args)),
                         'note': 'gene not measured'})
            continue
        parts = ({'start': {s: v[0] for s, v in vals.items()}, 'end': {s: v[1] for s, v in vals.items()},
                  'width': {s: v[1] - v[0] for s, v in vals.items()}} if kind == 'domain' else {'': vals})
        for part, q in parts.items():
            rows.append({'question': question, 'scope': scope, 'readout': f'{kind} {part}'.strip(),
                         'genes': ' > '.join(a if isinstance(a, str) else '+'.join(a) for a in args),
                         **{f'v_{s}': q[s] for s in secs}, **contrast(q, test, ref)})
    return rows


def aki_response_rows(get, grid, b12, b23, genes, test, ref, question, scope):
    """Per IR section: response = IR curve - mean control curve; its peak position, peak value and proximal front.
    Contrasted against the relabeling of the same construction (each 'IR' section against the mean of the other two)."""
    rows = []
    allsecs = test + ref
    for g in genes:
        if get(g, test[0]) is None:
            rows.append({'question': question, 'scope': scope, 'readout': 'response', 'genes': g, 'note': 'gene not measured'})
            continue
        cur = {s: get(g, s) for s in allsecs}
        out = {}
        for split, (a, b) in {'condition': (test, ref), 'relabel_1': ([test[0], ref[0]], [test[1], ref[1]]),
                              'relabel_2': ([test[0], ref[1]], [test[1], ref[0]])}.items():
            resp = np.mean([cur[s] for s in a], 0) - np.mean([cur[s] for s in b], 0)
            i = int(np.nanargmax(np.abs(resp)))
            sign = np.sign(resp[i])
            hit = np.flatnonzero(sign * resp >= 0.5 * abs(resp[i]))
            seg = [np.nanmean(resp[(grid >= lo) & (grid < hi)]) for lo, hi in ((0, b12), (b12, b23), (b23, 1.01))]
            out[split] = {'peak_x': grid[i], 'peak': resp[i], 'front_x': grid[hit[0]], 'tail_x': grid[hit[-1]],
                          'S1': seg[0], 'S2': seg[1], 'S3': seg[2]}
        pair_S = {s: {k: np.nanmean((cur[s] - np.mean([cur[r] for r in ref], 0))[(grid >= lo) & (grid < hi)])
                      for k, (lo, hi) in zip(('S1', 'S2', 'S3'), ((0, b12), (b12, b23), (b23, 1.01)))} for s in test}
        i = int(np.argmin(np.abs(grid - out['condition']['peak_x'])))
        pairs = [cur[t][i] - cur[r][i] for t in test for r in ref]
        same = all(np.sign(p) == np.sign(out['condition']['peak']) for p in pairs)
        # where the gene sits in IR itself (scaled IR level per IR section): the injury domain
        lvl = {s: domain(cur[s], grid) + (peak(cur[s], grid),) for s in test}
        ctl_max = max(float(np.nanmax(cur[r])) for r in ref)
        rows.append({'question': question, 'scope': scope, 'readout': 'response', 'genes': g,
                     **{f'ir_{k}_{s}': v[j] for s, v in lvl.items() for j, k in enumerate(('start', 'end', 'peak'))},
                     'control_max_log2': ctl_max, 'ir_max_log2': max(float(np.nanmax(cur[s])) for s in test),
                     **{f'{k}': v for k, v in out['condition'].items()}, 'pairs_same_sign': same,
                     **{f'{k}_IR2A2': v for k, v in pair_S['IR2A2'].items()},
                     **{f'{k}_IR2A4': v for k, v in pair_S['IR2A4'].items()},
                     'relabel_1_peak': out['relabel_1']['peak'], 'relabel_2_peak': out['relabel_2']['peak'],
                     'call': bool(same and abs(out['condition']['peak']) > max(abs(out['relabel_1']['peak']),
                                                                       abs(out['relabel_2']['peak'])))})
    return rows


def injury_origin(root: Path):
    """Per PT label, genes induced or repressed in all four IR - control pairs (>= 1 log2, label pseudobulk), and
    where their IR level sits against the healthy nephron (margin 0.5 log2). Run for the condition and both relabelings."""
    import anndata as ad
    from scipy import sparse
    qc = ad.read_h5ad(root / RUNS['aki'] / '02_qc/units_qc.h5ad')
    ref = pd.read_parquet(root / RUNS['aki'] / '04_pseudospace/pseudospace.parquet')
    ref = ref.set_index('obs_name') if 'obs_name' in ref.columns else ref
    units = ref.index.intersection(qc.obs_names)
    counts = sparse.csr_matrix(qc[units].layers['counts'])
    mt = qc.var['mt'].to_numpy(bool)
    lib = np.asarray(counts[:, ~mt].sum(1)).ravel()
    sample = ref.loc[units, 'sample'].astype(str).to_numpy()
    label = ref.loc[units, 'segment_class'].astype(str).to_numpy()
    pt, distal = ['PT-S1', 'PT-S2', 'PT-S3'], ['mTAL', 'cTAL', 'DCT2', 'CCD', 'IMCD']
    lev = {}
    for s in np.unique(sample):
        for l in pt + distal:
            k = (sample == s) & (label == l)
            lev[(s, l)] = np.log2(1e4 * np.asarray(counts[k].sum(0)).ravel() / lib[k].sum() + 0.5)
    genes = qc.var_names.to_numpy()

    def run(test, refs):
        mean = lambda secs, l: np.mean([lev[(s, l)] for s in secs], 0)  # noqa: E731
        pt_max = np.max([mean(refs, l) for l in pt], 0)
        pt_min = np.min([mean(refs, l) for l in pt], 0)
        d_lev = np.array([mean(refs, l) for l in distal])
        rows = []
        for l in pt:
            pairs = np.array([lev[(t, l)] - lev[(r, l)] for t in test for r in refs])
            ir = mean(test, l)
            for direction, hit in (('up', (pairs >= 1).all(0) & ~mt), ('down', (pairs <= -1).all(0) & ~mt)):
                for i in np.flatnonzero(hit):
                    if direction == 'up':
                        cls = ('within healthy PT range' if ir[i] <= pt_max[i] + 0.5 else
                               'level found in healthy distal nephron' if ir[i] <= d_lev[:, i].max() + 0.5 else
                               'above every healthy segment')
                    else:
                        cls = 'within healthy PT range' if ir[i] >= pt_min[i] - 0.5 else 'below every healthy PT label'
                    rows.append({'segment': l, 'direction': direction, 'gene': genes[i], 'class': cls,
                                 'ir_level': ir[i], 'healthy_pt_max': pt_max[i], 'healthy_distal_max': d_lev[:, i].max(),
                                 'healthy_distal_where': distal[int(np.argmax(d_lev[:, i]))],
                                 'healthy_slope': (mean(refs, 'PT-S3') - mean(refs, 'PT-S1'))[i],
                                 **{f'slope_{t}': (lev[(t, 'PT-S3')] - lev[(t, 'PT-S1')])[i] for t in test}})
        return pd.DataFrame(rows)

    test, refs = GROUPS['aki']
    relabel = {'relabel_1': run([test[0], refs[0]], [test[1], refs[1]]),
               'relabel_2': run([test[0], refs[1]], [test[1], refs[0]])}
    return run(test, refs), {k: len(v) for k, v in relabel.items()}


def main(argv=None):
    args = parse_args(argv)
    root = Path(args.results_root)
    out = root / 'paper_figures/question_bank'
    out.mkdir(parents=True, exist_ok=True)
    rows, saved = [], {}

    # human vs mouse
    test, ref = GROUPS['hvm']
    grid, get_pt, _, _, secs = pt_curves(root, 'fig5')
    distal_genes = sorted({g for _, w, spec in HVM if w != 'pt' for it in spec for g in it[1:]})
    prof, t, sample, label = nephron_profiles(root, 'hvm', distal_genes)
    saved.update({f'hvm|{g}|{s}': v for (g, s), v in prof.items()})
    for q, w, spec in HVM:
        if w == 'pt':
            rows += measure(spec, get_pt, grid, secs, test, ref, q, 'PT')
            continue
        lo, hi = window(t, sample, label, ref, WINDOWS[w])
        k = (NEPHRON_GRID >= lo) & (NEPHRON_GRID <= hi)
        get_w = lambda g, s: prof[(g, s)][k] if (g, s) in prof else None  # noqa: E731
        rows += measure(spec, get_w, NEPHRON_GRID[k], test + ref, test, ref, q, f'{w} t {lo:.2f}-{hi:.2f}')

    # AKI
    test, ref = GROUPS['aki']
    grid, get_pt, b12, b23, secs = pt_curves(root, 'fig6')
    rows += aki_response_rows(get_pt, grid, b12, b23, AKI_PT, test, ref, 'A PT responses', 'PT')
    prof, t, sample, label = nephron_profiles(root, 'aki', AKI_DISTAL)
    saved.update({f'aki|{g}|{s}': v for (g, s), v in prof.items()})
    get_n = lambda g, s: prof.get((g, s))  # noqa: E731
    lab = {l: np.median(t[np.isin(sample, ref) & (label == l)]) for l in np.unique(label)}
    # whole-nephron responses: per label window, IR minus control, per IR section
    for g in AKI_DISTAL:
        if (g, ref[0]) not in prof:
            continue
        for l in sorted(lab, key=lab.get):
            kk = np.isin(sample, ref) & (label == l)
            lo, hi = np.quantile(t[kk], [0.1, 0.9])
            k = (NEPHRON_GRID >= lo) & (NEPHRON_GRID <= hi)
            q = {s: np.nanmean(prof[(g, s)][k]) for s in test + ref}
            rows.append({'question': 'A nephron responses', 'scope': l, 'readout': 'level', 'genes': g,
                         **{f'v_{s}': q[s] for s in test + ref}, **contrast(q, test, ref)})
    lo, hi = window(t, sample, label, ref, ['DCT2', 'CCD'])
    k = (NEPHRON_GRID >= lo) & (NEPHRON_GRID <= hi)
    get_w = lambda g, s: prof[(g, s)][k] if (g, s) in prof else None  # noqa: E731
    spec = [('cross', a, b) for a, b in AKI_PAIRS] + [('domain', g) for g in ['Aqp2', 'Scnn1g', 'Slc12a3', 'Calb1']]
    rows += measure(spec, get_w, NEPHRON_GRID[k], test + ref, test, ref, 'A distal boundaries', f'DCT-CD t {lo:.2f}-{hi:.2f}')
    lo, hi = window(t, sample, label, ref, ['mTAL', 'cTAL'])
    k = (NEPHRON_GRID >= lo) & (NEPHRON_GRID <= hi)
    rows += aki_response_rows(lambda g, s: prof[(g, s)][k] if (g, s) in prof else None, NEPHRON_GRID[k],
                              lo + (hi - lo) / 3, lo + 2 * (hi - lo) / 3, ['Slc12a1', 'Umod', 'Kcnj1', 'Lcn2', 'Spp1'],
                              test, ref, 'A TAL responses', f'TAL t {lo:.2f}-{hi:.2f} (thirds)')

    origin, null = injury_origin(root)
    origin.to_csv(out / 'injury_origin.csv', index=False)
    (out / 'injury_origin_relabelings.json').write_text(json.dumps(null))
    print('injury origin: condition', len(origin), 'gene x label calls; relabelings', null)

    df = pd.DataFrame(rows)
    df.to_csv(out / 'readouts.csv', index=False)
    np.savez_compressed(out / 'profiles.npz', grid=NEPHRON_GRID, **saved)
    print(df.groupby('question').call.agg(['size', 'sum']).to_string())


if __name__ == '__main__':
    main()
