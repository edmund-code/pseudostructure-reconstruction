"""Fig. 5 screen: mouse against human along the PT, on the Fig. 4 axes of the human vs mouse fit.

Each gene is placed on the mean of the five Fig. 3d gene-half axes of the human vs mouse fit built without it; x is
the mouse quantile (human units are mapped through the mouse distribution). Per gene and section, NB GLMs give the
curve, the label means and the bulk mean. From those:
 gradients    E = S3-core minus S1-core mean of each section's curve, and the raw label-in-core contrast R
 category R   gradients that reverse between species; category H: steeper in human (relabeling-calibrated)
 pathways     sets whose gradient departs from the genome-wide flattening line, against random sets
 handoffs     notebook 101's pairs (mouse and human): crossings per section, shifts with bootstrap intervals and a
              warp check
 conservation how much of a gradient chosen in one section the other sections keep
Protocol (frozen before this ran): docs/paper/fig5-human-vs-mouse-along-pt.md.

  OMP_NUM_THREADS=1 python analysis/paper/fig5_screen.py --jobs 24
Writes <results-root>/paper_figures/fig5/screen/.
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
import fig6_screen as f6  # noqa: E402
from fig3_common import RUNS  # noqa: E402

MOUSE, HUMAN = ['Ctrl1A2', 'Ctrl1A4'], ['HUK1_COR1', 'HUK1_MED1']
SECTIONS = MOUSE + HUMAN
# the three balanced 2-vs-2 splits: (plus, minus); each one's "pairs" are its four cross differences
SPLITS = {'species': (HUMAN, MOUSE), 'relabeling 1': (['HUK1_COR1', 'Ctrl1A2'], ['HUK1_MED1', 'Ctrl1A4']),
          'relabeling 2': (['HUK1_COR1', 'Ctrl1A4'], ['HUK1_MED1', 'Ctrl1A2'])}
GRAD_LOG2, HUMAN_STEEP_LOG2, STEEPER_BY = 0.5, 1.0, 0.5
SPECIFIC = 0.25
MIN_SUPPORT = 25
MAX_FALSE_SHARE = 0.10
SHIFT_MIN, N_BOOT = 0.05, 200
MIN_SET, N_RANDOM_SETS, SET_Q = 10, 1000, 0.05
SELECT_LOG2 = 0.75
GRID = fs.GRID


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--results-root', default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT', str(REPO / 'results')))
    p.add_argument('--out', default='paper_figures/fig5/screen', help='relative to the results root')
    p.add_argument('--jobs', type=int, default=8)
    p.add_argument('--smoke', action='store_true', help='two axes, few bootstraps and random sets (code check only)')
    return p.parse_args(argv)


def split_value(q: dict, split: str):
    plus, minus = SPLITS[split]
    return np.mean([q[m] for m in plus], 0) - np.mean([q[m] for m in minus], 0)


def split_pairs(q: dict, split: str):
    plus, minus = SPLITS[split]
    return np.stack([q[a] - q[b] for a in plus for b in minus])


def mouse_quantile(mean, sample, ok):
    """Mouse units: rank among the mouse PT units (Fig. 4's x); human units: mapped through that distribution."""
    x = np.full(len(mean), np.nan)
    ms = ok & np.isin(sample, MOUSE)
    r = pd.Series(mean[ms]).rank(method='average').to_numpy()
    x[ms] = (r - 0.5) / ms.sum()
    order = np.argsort(mean[ms])
    hu = ok & np.isin(sample, HUMAN)
    x[hu] = np.interp(mean[hu], mean[ms][order], x[ms][order])
    return x


def end_contrast(curves: np.ndarray, b12, b23, xr) -> np.ndarray:
    """S3-core mean minus S1-core mean of curves on GRID (..., grid), over grid points inside the section's x range."""
    inside = (GRID >= xr[0]) & (GRID <= xr[1])
    s1 = inside & (GRID <= b12 - fs.MARGIN)
    s3 = inside & (GRID >= b23 + fs.MARGIN)
    if not s1.any() or not s3.any():
        return np.full(curves.shape[:-1], np.nan)
    return curves[..., s3].mean(-1) - curves[..., s1].mean(-1)


def raw_contrast(counts, lib, x, label, b12, b23):
    """log2 rate (summed counts + 0.5, per 10k): S3-labelled units in the S3 core minus S1-labelled units in the S1
    core; and the two unit counts."""
    rates, n = [], []
    for k in ((label == 0) & (x <= b12 - fs.MARGIN), (label == 2) & (x >= b23 + fs.MARGIN)):
        rates.append(np.log2(1e4 * (np.asarray(counts[k].sum(0)).ravel() + 0.5) / max(lib[k].sum(), 1.0)))
        n.append(int(k.sum()))
    return rates[1] - rates[0], n


def fit_section(y, lib, x, label, knots, jobs):
    from pseudospace.nb_gam import nb_fit

    w, off = np.ones(len(x)), np.log(lib)
    disp, curve, grid_curve = f6.fit_curves(y, lib, x, label, knots, jobs)
    clusters = nb_fit(y, np.eye(3)[label], w, off, disp, n_jobs=jobs)
    bulk = nb_fit(y, np.ones((len(x), 1)), w, off, disp, n_jobs=jobs)
    order = np.argsort(x)
    bins = np.array_split(order, fs.N_BINS)
    with np.errstate(divide='ignore'):
        bin_means = np.stack([np.log2(1e4 * y[b].sum(0) / lib[b].sum()) for b in bins], 1)
    return {'alpha': disp['alpha'], 'curve': grid_curve, 'converged': curve['converged'] & clusters['converged'],
            'clusters': fs.log2_cp10k(clusters['beta']), 'bulk': fs.log2_cp10k(bulk['beta'][:, 0]),
            'bins': bin_means, 'bin_x': np.array([x[b].mean() for b in bins])}


# ------------------------------------------------------------------------------------------------- rules
def category_r(E: dict, R: dict, ok: np.ndarray, plus, minus) -> np.ndarray:
    """The gradient reverses between the two arms: one sign in both minus sections, the other in both plus
    sections, |E| >= GRAD_LOG2 in all four, raw R agreeing with each section's E."""
    s = np.sign(E[minus[0]])
    secs = list(plus) + list(minus)
    r1 = (s != 0) & (np.sign(E[minus[1]]) == s) & np.all([np.sign(E[p]) == -s for p in plus], 0)
    r1 &= np.all([np.abs(E[m]) >= GRAD_LOG2 for m in secs], 0)
    r2 = np.all([np.sign(R[m]) == np.sign(E[m]) for m in secs], 0)
    return r1 & r2 & ok


def category_h(E: dict, R: dict, ok: np.ndarray, plus, minus) -> np.ndarray:
    """Steeper in the plus arm: one sign everywhere, |E| >= HUMAN_STEEP_LOG2 in both plus sections, steeper by
    STEEPER_BY in all four cross pairs, raw R with that sign everywhere and steeper in all four cross pairs."""
    s = np.sign(E[plus[0]])
    secs = list(plus) + list(minus)
    same = (s != 0) & np.all([np.sign(E[m]) == s for m in secs], 0)
    big = np.all([np.abs(E[p]) >= HUMAN_STEEP_LOG2 for p in plus], 0)
    steeper = np.all([np.abs(E[p]) - np.abs(E[q]) >= STEEPER_BY for p in plus for q in minus], 0)
    raw = np.all([np.sign(R[m]) == s for m in secs], 0) & \
        np.all([np.abs(R[p]) > np.abs(R[q]) for p in plus for q in minus], 0)
    return same & big & steeper & raw & ok


def flattening_residuals(E: dict, valid: np.ndarray, split: str):
    """For each cross pair (a in plus, b in minus): OLS of E[a] on E[b] over the valid genes; residuals (NaN where
    not valid) and (intercept, slope). Returns the split's mean residual, the four pair residuals and the lines."""
    plus, minus = SPLITS[split]
    res, lines = [], []
    for a in plus:
        for b in minus:
            k = valid & np.isfinite(E[a]) & np.isfinite(E[b])
            slope, icpt = np.polyfit(E[b][k], E[a][k], 1)
            r = np.full(len(valid), np.nan)
            r[k] = E[a][k] - (icpt + slope * E[b][k])
            res.append(r)
            lines.append((float(icpt), float(slope)))
    res = np.stack(res)
    return res.mean(0), res, lines


def set_screen(gene_d: dict, members: dict, rng) -> pd.DataFrame:
    """gene_d: {split or 'split|k': (G,) residuals over the screened genes, NaN where unsupported}; members: {set:
    index array}. A set uses its valid members (>= MIN_SET) against random sets of that size from the valid genes;
    two-sided empirical p, BH within each split."""
    rows = []
    for split in SPLITS:
        d = gene_d[split]
        valid = np.isfinite(d)
        pool = np.flatnonzero(valid)
        null = {}
        for name, idx in members.items():
            ix = idx[valid[idx]]
            if len(ix) < MIN_SET:
                continue
            if len(ix) not in null:
                null[len(ix)] = np.array([d[rng.choice(pool, len(ix), replace=False)].mean()
                                          for _ in range(N_RANDOM_SETS)])
            obs, nl = d[ix].mean(), null[len(ix)]
            centre = nl.mean()
            p = (1 + (np.abs(nl - centre) >= abs(obs - centre)).sum()) / (1 + N_RANDOM_SETS)
            pairs = np.array([gene_d[f'{split}|{k}'][ix].mean() for k in range(4)])
            others = np.array([gene_d[o][ix].mean() for o in SPLITS if o != split])
            rows.append({'split': split, 'set': name, 'n': len(ix), 'd': obs, 'p': p,
                         'pairs_agree': bool((np.sign(pairs) == np.sign(obs)).all()),
                         'specific': bool((np.abs(others) <= SPECIFIC * abs(obs)).all()),
                         **{f'pair_{k}': pairs[k] for k in range(4)}, **{f'other_{j}': others[j] for j in range(2)}})
    df = pd.DataFrame(rows)
    df['q'] = np.nan
    for split in SPLITS:
        k = df.split == split
        df.loc[k, 'q'] = f6.bh(df.loc[k, 'p'].to_numpy())
    df['call'] = (df.q <= SET_Q) & df.pairs_agree & df.specific
    return df


def switch_point(profile: np.ndarray, rising: bool) -> float:
    """Where a scaled profile first crosses 0.5 in its gradient's direction."""
    from pseudospace.handoffs import crossing

    half = np.full_like(profile, 0.5)
    return crossing(half, profile, GRID) if rising else crossing(profile, half, GRID)


# ------------------------------------------------------------------------------------------------- main
def main(argv=None):
    import anndata as ad
    from scipy import sparse
    from pseudospace.handoffs import crossing, scaled_profile
    from pseudospace.stats_gam import gam_internal_knots

    global N_BOOT, N_RANDOM_SETS
    args = parse_args(argv)
    root = Path(args.results_root)
    out = root / args.out
    if args.smoke:
        out, N_BOOT, N_RANDOM_SETS = out.with_name(out.name + '_smoke'), 3, 50
    out.mkdir(parents=True, exist_ok=True)

    qc = ad.read_h5ad(root / RUNS['hvm'] / '02_qc/units_qc.h5ad')
    ref = pd.read_parquet(root / RUNS['hvm'] / '04_pseudospace/pseudospace.parquet')
    ref = ref.set_index('obs_name') if 'obs_name' in ref.columns else ref
    pt = ref[ref['sample'].astype(str).isin(SECTIONS) & ref.segment_class.astype(str).isin(fs.PT)]
    units = pt.index.to_numpy()
    counts = sparse.csc_matrix(qc[units].layers['counts'])
    mt = qc.var['mt'].to_numpy(bool)
    lib = np.asarray(counts[:, ~mt].sum(1)).ravel()
    sample = pt['sample'].astype(str).to_numpy()
    label = pt.segment_class.astype(str).map(fs.PT).to_numpy()
    mito = np.asarray(counts[:, mt].sum(1)).ravel() / np.asarray(counts.sum(1)).ravel()

    det = {m: np.asarray((counts[sample == m] > 0).mean(0)).ravel() for m in SECTIONS}
    universe = np.all([det[m] >= fs.MIN_DETECTED for m in SECTIONS], 0) & ~mt

    ratio = {}
    for sp, secs in (('mouse', MOUSE), ('human', HUMAN)):
        k = np.isin(sample, secs)
        pt_rate = np.asarray(counts[k].sum(0)).ravel() / lib[k].sum()
        nb = ref[ref['sample'].astype(str).isin(secs) & ref.segment_class.astype(str).isin(fs.NEIGHBOURS)]
        nbc = sparse.csr_matrix(qc[nb.index.to_numpy()].layers['counts'])
        nblib = np.asarray(nbc[:, ~mt].sum(1)).ravel()
        group = nb.segment_class.astype(str).map(fs.NEIGHBOURS).to_numpy()
        rates = [np.asarray(nbc[group == g].sum(0)).ravel() / nblib[group == g].sum() for g in ('AL', 'DCT', 'CD')]
        ratio[sp] = np.max(np.stack(rates), 0) / np.maximum(pt_rate, 1e-12)
    neighbour_ok = (ratio['mouse'] <= fs.NEIGHBOUR_RATIO) & (ratio['human'] <= fs.NEIGHBOUR_RATIO)

    p101 = pd.read_csv(root / 'pt_transporter_handoffs/tables/pairs.csv')
    p101 = p101[p101.confirmed.astype(bool)]
    src = p101.groupby(['early', 'late']).species.agg(lambda v: ' and '.join(sorted(set(v)))).reset_index()
    fam = p101.drop_duplicates(['early', 'late'])[['early', 'late', 'family']]
    pairs_101 = src.merge(fam, on=['early', 'late']).rename(columns={'species': 'defined_in'})
    poster = pd.DataFrame({'early': ['Slc5a12', 'Slc5a2'], 'late': ['Slc5a8', 'Slc5a1'], 'defined_in': 'poster',
                           'family': 'poster'})
    handoff_pairs = pd.concat([pairs_101, poster], ignore_index=True)
    handoff_genes = sorted(set(handoff_pairs.early) | set(handoff_pairs.late))
    handoff_genes = [g for g in handoff_genes if g in qc.var_names]
    fit_mask = universe.copy()
    fit_mask[[qc.var_names.get_loc(g) for g in handoff_genes]] = True

    halves = fs.gene_halves(qc.var)
    axes = {(k, h): pd.read_parquet(root / 'paper_figures/fig3/robustness/hvm/reruns' /
                                    f'gene_half_s{k}_{h}.parquet')['t'].reindex(units).to_numpy(float)
            for k in range(fs.N_SPLITS) for h in 'AB'}
    names = qc.var_names.to_numpy()
    pattern = np.array([''.join('B' if g in halves[k] else 'A' for k in range(fs.N_SPLITS)) for g in names])

    rows, store, axis_rows, boot_inputs, xranges = [], {}, [], {}, {}
    patterns = sorted(set(pattern[fit_mask]))
    if args.smoke:
        patterns = patterns[:2]
        fit_mask &= np.isin(pattern, patterns)
        handoff_genes = [g for g in handoff_genes if fit_mask[qc.var_names.get_loc(g)]]
    for pi, pat in enumerate(patterns):
        pos = np.stack([axes[(k, pat[k])] for k in range(fs.N_SPLITS)], 1)
        ok = np.isfinite(pos).sum(1) >= fs.MIN_AXES
        mean = np.where(ok, np.nanmean(np.where(np.isfinite(pos), pos, np.nan), 1), np.nan)
        x = mouse_quantile(mean, sample, ok)
        ms = ok & np.isin(sample, MOUSE)
        b12, b23 = fs.boundary(x[ms], label[ms], 0), fs.boundary(x[ms], label[ms], 1)
        knots = gam_internal_knots(x[ms], fs.BASIS_DF)
        lo_m, hi_m = mean[ms].min(), mean[ms].max()
        xr = {m: tuple(np.quantile(x[ok & (sample == m)], [0.025, 0.975])) for m in SECTIONS}
        xranges[pat] = xr
        axis_rows.append({'pattern': pat, 'b12': b12, 'b23': b23,
                          **{f'n_{m}': int((ok & (sample == m)).sum()) for m in SECTIONS},
                          **{f'x_lo_{m}': xr[m][0] for m in SECTIONS}, **{f'x_hi_{m}': xr[m][1] for m in SECTIONS},
                          **{f'clipped_low_{m}': float((mean[ok & (sample == m)] < lo_m).mean()) for m in HUMAN},
                          **{f'clipped_high_{m}': float((mean[ok & (sample == m)] > hi_m).mean()) for m in HUMAN}})
        gi = np.flatnonzero(fit_mask & (pattern == pat))
        sub = counts[:, gi]
        fits, raw, nraw, zero = {}, {}, {}, {}
        for m in SECTIONS:
            k = ok & (sample == m)
            y = sub[k].toarray()
            fits[m] = fit_section(y, lib[k], x[k], label[k], knots, args.jobs)
            raw[m], nraw[m] = raw_contrast(sub[k], lib[k], x[k], label[k], b12, b23)
            zero[m] = f6.core_zero(sub[k], x[k], b12, b23)[:, [0, 2]].any(1)
            # the lowest-mito tertile of each human section (named cell-state check)
            if m in HUMAN:
                lowm = mito[k] <= np.quantile(mito[k], 1 / 3)
                raw[f'{m}|lowmito'], _ = raw_contrast(sub[k][lowm], lib[k][lowm], x[k][lowm], label[k][lowm], b12, b23)
            hk = np.isin(names[gi], handoff_genes)
            if hk.any():
                boot_inputs.setdefault(pat, {})[m] = {'x': x[k], 'lib': lib[k], 'y': y[:, hk], 'genes': names[gi][hk],
                                                      'alpha': fits[m]['alpha'][hk], 'knots': knots}
        support = all(min(nraw[m]) >= MIN_SUPPORT for m in SECTIONS)
        for j, g in enumerate(names[gi]):
            store[g] = {m: (fits[m]['curve'][j], fits[m]['clusters'][j], fits[m]['bulk'][j], fits[m]['bins'][j],
                            fits[m]['bin_x']) for m in SECTIONS}
            row = {'gene': g, 'pattern': pat, 'b12': b12, 'b23': b23, 'universe': bool(universe[gi[j]]),
                   'neighbour_ok': bool(neighbour_ok[gi[j]]), 'neighbour_ratio_mouse': ratio['mouse'][gi[j]],
                   'neighbour_ratio_human': ratio['human'][gi[j]],
                   'converged': bool(np.all([fits[m]['converged'][j] for m in SECTIONS])), 'support': support,
                   'unsupported': bool(np.any([zero[m][j] for m in SECTIONS])),
                   **{f'detected_{m}': det[m][gi[j]] for m in SECTIONS}}
            for m in SECTIONS:
                row[f'E_{m}'] = float(end_contrast(fits[m]['curve'][j], b12, b23, xr[m]))
                row[f'R_{m}'] = raw[m][j]
                row[f'bulk_{m}'] = fits[m]['bulk'][j]
                row.update({f'cluster_S{s + 1}_{m}': fits[m]['clusters'][j, s] for s in range(3)})
                if m in HUMAN:
                    row[f'R_lowmito_{m}'] = raw[f'{m}|lowmito'][j]
            rows.append(row)
        print(f'axis {pi + 1}/{len(patterns)}: {pat}  {len(gi)} genes  b12 {b12:.3f} b23 {b23:.3f}', flush=True)

    df = pd.DataFrame(rows).set_index('gene')
    genes = df.index.to_numpy()
    E = {m: df[f'E_{m}'].to_numpy() for m in SECTIONS}
    R = {m: df[f'R_{m}'].to_numpy() for m in SECTIONS}
    screened = (df.universe & df.neighbour_ok & df.converged).to_numpy()
    assessable = screened & df.support.to_numpy() & ~df.unsupported.to_numpy() & \
        np.all([np.isfinite(E[m]) for m in SECTIONS], 0)
    df['screened'], df['assessable'] = screened, assessable

    # ---------------------------------------------------------------------------------------------- categories
    cat = {}
    for split, (plus, minus) in SPLITS.items():
        cat[('R', split)] = category_r(E, R, assessable, plus, minus)
        cat[('H', split)] = category_h(E, R, assessable, plus, minus)
        cat[('H reverse', split)] = category_h(E, R, assessable, minus, plus)
    df['catR'], df['catH'] = cat[('R', 'species')], cat[('H', 'species')]
    df['steeper_in_mouse'] = cat[('H reverse', 'species')]
    for split in SPLITS:
        df[f'E_{split}'] = split_value(E, split)
    df['E_min_abs'] = np.min([np.abs(E[m]) for m in SECTIONS], 0)
    df['H_min_pair'] = np.min([np.abs(E[p]) - np.abs(E[q]) for p in HUMAN for q in MOUSE], 0)

    def calibration(c):
        real = int(cat[(c, 'species')].sum())
        rel = [int(cat[(c, s)].sum()) for s in ('relabeling 1', 'relabeling 2')]
        share = np.mean(rel) / max(real, 1)
        return {'calls': real, 'relabeling_calls': rel, 'false_share': share, 'calibrated': share <= MAX_FALSE_SHARE}

    # ---------------------------------------------------------------------------------------------- flattening, sets
    gene_d, lines = {}, {}
    for split in SPLITS:
        d, pr, ln = flattening_residuals(E, assessable, split)
        gene_d[split] = d
        lines[split] = ln
        for k in range(4):
            gene_d[f'{split}|{k}'] = pr[k]
        df[f'd_{split}'] = d
    for k, (a, b) in enumerate((a, b) for a in HUMAN for b in MOUSE):
        df[f'd_pair_{a}_{b}'] = gene_d[f'species|{k}']

    sets = json.loads((root / 'pt_package_aki/tables/gene_sets.json').read_text())
    pos_of = {g: i for i, g in enumerate(genes)}
    members = {}
    for name, mem in sets.items():
        idx = np.array([pos_of[g] for g in mem if g in pos_of and assessable[pos_of[g]]], int)
        if len(idx) >= MIN_SET:
            members[name] = idx
    sets_df = set_screen(gene_d, members, np.random.default_rng(55))
    sets_df.to_csv(out / 'pathways.csv', index=False)
    n_set_calls = {s: int(sets_df[sets_df.split == s].call.sum()) for s in SPLITS}
    set_share = np.mean([n_set_calls['relabeling 1'], n_set_calls['relabeling 2']]) / max(n_set_calls['species'], 1)

    # ---------------------------------------------------------------------------------------------- handoffs
    pat0 = pattern[qc.var_names.get_loc(handoff_genes[0])]
    xrange_ = xranges[pat0]

    def prof(c, m):
        return scaled_profile(np.where((GRID >= xrange_[m][0]) & (GRID <= xrange_[m][1]), c, np.nan))

    def pair_positions(get_curve):
        o = np.full((len(handoff_pairs), len(SECTIONS)), np.nan)
        for i, (e, l) in enumerate(zip(handoff_pairs.early, handoff_pairs.late)):
            if e not in store or l not in store:
                continue
            for j, m in enumerate(SECTIONS):
                o[i, j] = crossing(prof(get_curve(e, m), m), prof(get_curve(l, m), m), GRID)
        return o

    xs_obs = pair_positions(lambda g, m: store[g][m][0])
    from pseudospace.nb_gam import nb_fit
    from pseudospace.stats_gam import make_gam_design
    brng = np.random.default_rng(5)
    boot = np.full((N_BOOT, len(handoff_pairs), len(SECTIONS)), np.nan)
    for b in range(N_BOOT):
        bcurve = {}
        draws = {m: brng.integers(0, 1 << 31) for m in SECTIONS}
        for per in boot_inputs.values():
            for m, d in per.items():
                idx = np.random.default_rng(draws[m]).integers(0, len(d['x']), len(d['x']))
                fit = nb_fit(d['y'][idx], make_gam_design(d['x'][idx], d['knots']), np.ones(len(idx)),
                             np.log(d['lib'][idx]), d['alpha'])
                cv = fs.log2_cp10k(fit['beta'] @ make_gam_design(GRID, d['knots']).T)
                for j, g in enumerate(d['genes']):
                    bcurve[(g, m)] = cv[j]
        boot[b] = pair_positions(lambda g, m: bcurve[(g, m)])
        if (b + 1) % 50 == 0:
            print(f'bootstrap {b + 1}/{N_BOOT}', flush=True)

    # warp check: each section's switch points of the assessable genes graded there
    switch = {}
    for m in SECTIONS:
        sp = []
        for i in np.flatnonzero(assessable & (np.abs(E[m]) >= GRAD_LOG2)):
            sp.append(switch_point(prof(store[genes[i]][m][0], m), E[m][i] > 0))
        sp = np.array(sp)
        switch[m] = np.sort(sp[np.isfinite(sp)])

    hp = handoff_pairs.copy()
    q = {m: xs_obs[:, j] for j, m in enumerate(SECTIONS)}
    for j, m in enumerate(SECTIONS):
        hp[f'x_{m}'] = xs_obs[:, j]
        hp[f'x_lo_{m}'] = np.nanpercentile(boot[:, :, j], 2.5, 0)
        hp[f'x_hi_{m}'] = np.nanpercentile(boot[:, :, j], 97.5, 0)
        hp[f'boot_present_{m}'] = np.isfinite(boot[:, :, j]).mean(0)
        hp[f'quantile_{m}'] = [np.searchsorted(switch[m], v, 'right') / len(switch[m]) if np.isfinite(v) else np.nan
                               for v in xs_obs[:, j]]
    in_mouse = np.isfinite(q['Ctrl1A2']) & np.isfinite(q['Ctrl1A4'])
    in_human = np.isfinite(q['HUK1_COR1']) & np.isfinite(q['HUK1_MED1'])
    hp['in_mouse'], hp['in_human'] = in_mouse, in_human
    hp['eligible'] = in_mouse & in_human
    hp['shift'] = split_value(q, 'species')
    pr = split_pairs(q, 'species')
    hp['pairs_agree'] = (np.sign(pr) == np.sign(hp['shift'].to_numpy())[None]).all(0)
    for split in ('relabeling 1', 'relabeling 2'):
        hp[f'shift_{split}'] = split_value(q, split)
    hp['specific'] = (hp['shift_relabeling 1'].abs() <= SPECIFIC * hp['shift'].abs()) & \
                     (hp['shift_relabeling 2'].abs() <= SPECIFIC * hp['shift'].abs())
    bshift = split_value({m: boot[:, :, j] for j, m in enumerate(SECTIONS)}, 'species')
    hp['shift_lo'] = np.nanpercentile(bshift, 2.5, 0)
    hp['shift_hi'] = np.nanpercentile(bshift, 97.5, 0)
    hp['shifted'] = hp.eligible & hp.pairs_agree & (hp['shift'].abs() >= SHIFT_MIN) & hp.specific & \
        ((hp.shift_lo > 0) | (hp.shift_hi < 0))
    qq = {m: hp[f'quantile_{m}'].to_numpy() for m in SECTIONS}
    hp['quantile_shift'] = split_value(qq, 'species')
    qpr = split_pairs(qq, 'species')
    hp['beyond_warp'] = hp.shifted & (np.sign(qpr) == np.sign(hp['shift'].to_numpy())[None]).all(0)
    hp['mouse_only'] = in_mouse & ~np.isfinite(q['HUK1_COR1']) & ~np.isfinite(q['HUK1_MED1'])
    hp['human_only'] = in_human & ~np.isfinite(q['Ctrl1A2']) & ~np.isfinite(q['Ctrl1A4'])
    for role in ('early', 'late'):
        for m in SECTIONS:
            hp[f'E_{role}_{m}'] = [df.loc[g, f'E_{m}'] if g in df.index else np.nan for g in hp[role]]
        hp[f'neighbour_ok_{role}'] = [bool(df.loc[g, 'neighbour_ok']) if g in df.index else False for g in hp[role]]
    hp.to_csv(out / 'handoffs.csv', index=False)
    np.save(out / 'handoff_bootstrap.npy', boot)

    # ---------------------------------------------------------------------------------------------- conservation
    crow = []
    for sel in SECTIONS:
        pick = assessable & (np.abs(E[sel]) >= SELECT_LOG2)
        sign, den = np.sign(E[sel][pick]), np.median(np.abs(E[sel][pick]))
        crow.append({'selected_in': sel, 'n': int(pick.sum()),
                     **{f'kept_{m}': float(np.median(E[m][pick] * sign) / den) for m in SECTIONS if m != sel},
                     **{f'same_sign_{m}': float((np.sign(E[m][pick]) == sign).mean()) for m in SECTIONS if m != sel}})
    pd.DataFrame(crow).to_csv(out / 'conservation.csv', index=False)

    # ---------------------------------------------------------------------------------------------- save
    df.to_csv(out / 'genes.csv')
    pd.DataFrame(axis_rows).to_csv(out / 'axes.csv', index=False)
    np.savez_compressed(out / 'all_curves.npz', grid=GRID, genes=genes,
                        **{f'{part}|{m}': np.stack([store[g][m][i] for g in genes]).astype(np.float32)
                           for m in SECTIONS for i, part in enumerate(('curve', 'clusters', 'bulk', 'bins', 'bin_x'))})
    summary = {'n_units': {m: int((sample == m).sum()) for m in SECTIONS}, 'n_screened': int(screened.sum()),
               'n_assessable': int(assessable.sum()), 'n_neighbour_dropped': int((df.universe & ~df.neighbour_ok).sum()),
               'n_axes': len(axis_rows), 'mito_fraction_median': {m: float(np.median(mito[sample == m])) for m in SECTIONS},
               'category_R': calibration('R'), 'category_H': calibration('H'),
               'steeper_in_mouse_count': int(cat[('H reverse', 'species')].sum()),
               'steeper_reverse_relabeling_calls': [int(cat[('H reverse', s)].sum()) for s in ('relabeling 1', 'relabeling 2')],
               'flattening_lines': lines,
               'pathways': {'sets': len(members), 'calls': n_set_calls, 'calibration_share': set_share,
                            'calibrated': set_share <= MAX_FALSE_SHARE},
               'handoffs': {'pairs': int(len(hp)), 'in_mouse': int(in_mouse.sum()), 'in_human': int(in_human.sum()),
                            'eligible': int(hp.eligible.sum()), 'shifted': int(hp.shifted.sum()),
                            'beyond_warp': int(hp.beyond_warp.sum()), 'mouse_only': int(hp.mouse_only.sum()),
                            'human_only': int(hp.human_only.sum()),
                            'shifted_pairs': hp.loc[hp.shifted, ['early', 'late', 'shift']].round(3).values.tolist(),
                            'mouse_only_pairs': hp.loc[hp.mouse_only, ['early', 'late']].values.tolist(),
                            'human_only_pairs': hp.loc[hp.human_only, ['early', 'late']].values.tolist()},
               'switch_points': {m: int(len(v)) for m, v in switch.items()}}
    (out / 'summary.json').write_text(json.dumps(summary, indent=1, default=str))
    print(json.dumps(summary, indent=1, default=str))


if __name__ == '__main__':
    main()
