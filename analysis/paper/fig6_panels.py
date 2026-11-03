"""Fig. 6 panels: control against AKI along the PT, drawn from fig6_placement.py and fig6_screen.py.

 a  one control and one IR section, PT tubules coloured by position on the shared axis
 b  where each section's PT tubules sit on the axis (density per section)
 c  placement check: positions on the injury-blind axis against the reference, per section and label
 d  injury markers along the PT, control and IR, each mouse
 e  programs along the PT (cell cycle, transport, injury markers): the AKI response against the segment view
 f  every gene with a large AKI response: response along the PT beside its segment view
 g  genes whose AKI response changes inside a segment (category A)
 h  transporter handoffs: crossing positions in control and IR
 i  pathways whose AKI response changes inside a segment

Reads <results-root>/paper_figures/fig6/{placement,screen}/ and writes PNGs to <results-root>/paper_figures/fig6/panels/.
Protocol: docs/paper/fig6-aki-along-pt.md.
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
from fig3_panels import COLOURS, INK, MUTED, SECTION_COLOURS, style  # noqa: E402

CONTROLS, INJURED = ['Ctrl1A2', 'Ctrl1A4'], ['IR2A2', 'IR2A4']
MICE = CONTROLS + INJURED
SEGS = ['PTS1', 'PTS2', 'PTS3']
PT = {'PT-S1': 0, 'PT-S2': 1, 'PT-S3': 2}
CTRL, IR = '#2F6DB5', '#D9692B'
STEP = '#555555'
UP, DOWN = '#B2182B', '#2166AC'     # direction of an AKI response, as in panel f's colour map
# panel boxes of make_fig6.py in SVG units (1 unit = 0.01 inch here)
SIZE = {'a': (360, 242), 'b': (380, 242), 'c': (380, 242), 'd': (580, 232), 'e': (560, 232),
        'f': (480, 432), 'g': (660, 192), 'h': (660, 222),
        'i': (1160, 212)}


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--results-root', default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT', str(REPO / 'results')))
    p.add_argument('--out', default='paper_figures/fig6', help='relative to the results root')
    p.add_argument('--panels', default='abcdefghi')
    p.add_argument('--screen', default='screen', help='screen folder under --out')
    return p.parse_args(argv)


def figure(plt, key):
    w, h = SIZE[key]
    return plt.figure(figsize=(w / 100, h / 100))


def segment_band(ax, b12, b23, y, h, alpha=1.0):
    from matplotlib.patches import Rectangle

    for (lo, hi), seg in zip(((0, b12), (b12, b23), (b23, 1)), SEGS):
        ax.add_patch(Rectangle((lo, y), hi - lo, h, color=COLOURS[seg], lw=0, clip_on=False, alpha=alpha))


def reference_x(root: Path) -> pd.DataFrame:
    """PT units of the reference fit with x = control quantile of the package position (display only)."""
    ref = pd.read_parquet(root / RUNS['aki'] / '04_pseudospace/pseudospace.parquet')
    pt = ref[ref.segment_class.astype(str).isin(PT)].copy()
    pt['sample'] = pt['sample'].astype(str)
    ctrl = pt['sample'].isin(CONTROLS)
    c = np.sort(pt.loc[ctrl, 'pseudospace'].to_numpy())
    q = (np.arange(len(c)) + 0.5) / len(c)
    pt['x'] = np.interp(pt['pseudospace'], c, q)
    pt['label'] = pt.segment_class.astype(str).map(PT)
    return pt


# ---------------------------------------------------------------------------------------------------- a
def panel_a(plt, root, out):
    import anndata as ad

    pt = reference_x(root)
    obs = ad.read_h5ad(root / RUNS['aki'] / '02_qc/units_qc.h5ad', backed='r').obs
    ref = pd.read_parquet(root / RUNS['aki'] / '04_pseudospace/pseudospace.parquet')
    fig = figure(plt, 'a')
    cmap = plt.get_cmap('viridis')
    for i, s in enumerate(('Ctrl1A2', 'IR2A2')):
        ax = fig.add_axes([0.01 + i * 0.47, 0.13, 0.46, 0.80])
        units = ref.index[ref['sample'].astype(str) == s]
        xy = obs.loc[units, ['x_centroid', 'y_centroid']].to_numpy(float)
        other = ~units.isin(pt.index)
        ax.scatter(xy[other, 0], xy[other, 1], s=0.15, color='#DDDDDD', lw=0, rasterized=True)
        p = pt[pt['sample'] == s]
        pxy = obs.loc[p.index, ['x_centroid', 'y_centroid']].to_numpy(float)
        ax.scatter(pxy[:, 0], pxy[:, 1], s=0.35, c=p['x'], cmap=cmap, vmin=0, vmax=1, lw=0, rasterized=True)
        ax.set_aspect('equal')
        ax.invert_yaxis()
        ax.axis('off')
        fig.text(0.01 + i * 0.47 + 0.23, 0.95, 'control' if i == 0 else 'AKI', fontsize=7.5, ha='center', va='top',
                 color=CTRL if i == 0 else IR)
    cax = fig.add_axes([0.30, 0.07, 0.40, 0.025])
    grad = np.linspace(0, 1, 256)[None]
    cax.imshow(grad, aspect='auto', cmap=cmap, extent=(0, 1, 0, 1))
    cax.set_yticks([])
    cax.set_xticks([0, 1])
    cax.set_xticklabels(['start of PT', 'end of PT'], fontsize=6)
    cax.tick_params(length=0, pad=1.5)
    for sp in cax.spines.values():
        sp.set_visible(False)
    fig.savefig(out / 'f6_a.png', facecolor='white', dpi=300)
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------- b
def panel_b(plt, root, out, b12, b23):
    from scipy.stats import gaussian_kde

    pt = reference_x(root)
    fig = figure(plt, 'b')
    ax = fig.add_axes([0.16, 0.16, 0.80, 0.80])
    grid = np.linspace(0, 1, 301)
    step = 1.0
    for i, s in enumerate(MICE[::-1]):
        x = pt.loc[pt['sample'] == s, 'x'].to_numpy()
        d = gaussian_kde(x, bw_method=0.06)(grid)
        d = d / d.max() * 0.9 * step
        y0 = i * step
        c = SECTION_COLOURS[s]
        ax.fill_between(grid, y0, y0 + d, color=c, alpha=0.35, lw=0)
        ax.plot(grid, y0 + d, color=c, lw=1.0)
        ax.text(-0.02, y0 + 0.25 * step, s, ha='right', va='center', fontsize=6.5, color=c)
    for b in (b12, b23):
        ax.axvline(b, color='#BBBBBB', lw=0.5, zorder=0)
    ax.set_ylim(-0.32 * step, len(MICE) * step)
    segment_band(ax, b12, b23, -0.28 * step, 0.12 * step)
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['start of PT', 'end'])
    ax.set_xlabel('pseudostructure (control quantile)', labelpad=1)
    fig.savefig(out / 'f6_b.png', facecolor='white')
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------- c
def panel_c(plt, root, out):
    pl = pd.read_csv(root / 'placement/placement.csv')
    fig = figure(plt, 'c')
    ax = fig.add_axes([0.15, 0.20, 0.82, 0.74])
    groups = ['all PT', 'PT-S1', 'PT-S2', 'PT-S3']
    names = {'all PT': 'all PT', 'PT-S1': 'S1', 'PT-S2': 'S2', 'PT-S3': 'S3'}
    for gi, g in enumerate(groups):
        for si, s in enumerate(MICE):
            r = pl[(pl['sample'] == s) & (pl.units == g)].iloc[0]
            xp = gi + (si - 1.5) * 0.17
            c = SECTION_COLOURS[s]
            ax.plot([xp, xp], [r.rho_half_A_vs_B_mean, r.rho_reference_vs_blind], color=c, lw=0.6, alpha=0.6)
            ax.scatter([xp], [r.rho_reference_vs_blind], s=16, color=c, zorder=3, lw=0)
            ax.scatter([xp], [r.rho_half_A_vs_B_mean], s=14, facecolor='white', edgecolor=c, lw=0.8, zorder=3)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([names[g] for g in groups])
    ax.set_xlim(-0.5, len(groups) - 0.5)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel('Spearman ρ')
    ax.scatter([0.02], [0.13], s=16, color=MUTED, lw=0, transform=ax.transAxes, clip_on=False)
    ax.text(0.05, 0.13, 'injury-blind axis vs reference', transform=ax.transAxes, va='center', fontsize=6, color=MUTED)
    ax.scatter([0.02], [0.05], s=14, facecolor='white', edgecolor=MUTED, lw=0.8, transform=ax.transAxes, clip_on=False)
    ax.text(0.05, 0.05, 'gene half A vs half B', transform=ax.transAxes, va='center', fontsize=6, color=MUTED)
    fig.savefig(out / 'f6_c.png', facecolor='white')
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------- data
def load(root: Path, screen='screen'):
    genes = pd.read_csv(root / screen / 'genes.csv').set_index('gene')
    z = np.load(root / screen / 'curves.npz')
    names = list(z['genes'])
    curves = {m: {part: dict(zip(names, z[f'{part}|{m}'])) for part in ('curve', 'clusters', 'bulk', 'bins', 'bin_x')}
              for m in MICE}
    programs = {k.split('|', 1)[1]: z[k] for k in z.files if k.startswith('program|')}
    return genes, z['grid'], curves, programs


def covered(grid, bin_x):
    """Grid points inside a mouse's sampled range (first to last equal-count bin centre): no curve is drawn
    where that mouse has no tubules."""
    return (grid >= bin_x[0]) & (grid <= bin_x[-1])


MIN_FINITE_BINS = 16     # of 20 equal-count bins with a count, to draw a mouse's curve


def clipped(curves, m, gene):
    """A mouse's curve clipped to its bin-mean range (+-0.25), as in Fig. 4: a spline that dives where counts run
    out, or wanders where a gene is barely detected, is not drawn beyond the data."""
    b = curves[m]['bins'][gene]
    b = b[np.isfinite(b)]
    c = curves[m]['curve'][gene]
    return np.clip(c, b.min() - 0.25, b.max() + 0.25) if len(b) else np.full_like(c, np.nan)


def gene_levels(ax, gene, row, grid, curves, ylim=None):
    """One gene in all four mice: curves (control blue, IR orange) where each mouse has tubules, bin means."""
    for m in MICE:
        c = SECTION_COLOURS[m]
        bx, bm = curves[m]['bin_x'][gene], curves[m]['bins'][gene]
        ok = np.isfinite(bm)
        ax.scatter(bx[ok], bm[ok], s=3, color=c, lw=0, alpha=0.7, zorder=3)
        if ok.sum() < MIN_FINITE_BINS:     # barely detected in this mouse: its bins only, no curve
            continue
        k = covered(grid, bx)
        ax.plot(grid[k], clipped(curves, m, gene)[k], color=c, lw=1.2, zorder=4)
    lo, hi = ylim or ax.get_ylim()
    ax.set_ylim(lo, hi)
    pad = (hi - lo) * 0.06
    segment_band(ax, row.b12, row.b23, lo - 2.2 * pad, pad)
    for b in (row.b12, row.b23):
        ax.axvline(b, color='#CCCCCC', lw=0.5, zorder=0)
    ax.set_xlim(0, 1)
    ax.set_xticks([])


def level_ylim(gene, curves, pad=0.6):
    v = np.concatenate([curves[m]['bins'][gene] for m in MICE])
    v = v[np.isfinite(v)]
    return v.min() - pad, v.max() + pad


# ---------------------------------------------------------------------------------------------------- d
D_GENES = ['Havcr1', 'Krt20', 'Mki67', 'Slc34a1']


def panel_d(plt, genes, grid, curves, out):
    fig = figure(plt, 'd')
    for i, g in enumerate(D_GENES):
        ax = fig.add_axes([0.07 + i * 0.235, 0.12, 0.19, 0.76])
        gene_levels(ax, g, genes.loc[g], grid, curves, level_ylim(g, curves))
        ax.set_title(g, fontsize=7.5, style='italic', color=INK, pad=2)
        if i == 0:
            ax.set_ylabel('log2 counts per 10k')
    fig.savefig(out / 'f6_d.png', facecolor='white')
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------- e
E_PROGRAMS = [('C', 'cell-cycle program'), ('T_PT', 'PT transporters')]


def staircase(b12, b23, values):
    return np.array([0, b12, b12, b23, b23, 1]), np.repeat(values, 2)


def panel_e(plt, grid, programs, out, b12, b23):
    fig = figure(plt, 'e')
    for i, (key, title) in enumerate(E_PROGRAMS):
        ax = fig.add_axes([0.09 + i * 0.48, 0.12, 0.38, 0.76])
        for r in ('relabeling 1', 'relabeling 2'):
            ax.plot(grid, programs[f'{key}|{r}'], color='#9A9A9A', lw=0.8, ls=(0, (2, 1.5)))
        for a in INJURED:
            for c in CONTROLS:
                ax.plot(grid, programs[f'{key}|{a}−{c}'], color=IR, lw=0.5, alpha=0.5)
        sx, sy = staircase(b12, b23, programs[f'{key}|segment_view'])
        ax.plot(sx, sy, color=STEP, lw=1.4, solid_capstyle='butt')
        ax.plot(grid, programs[f'{key}|condition'], color=IR, lw=1.8)
        ax.axhline(0, color='#BBBBBB', lw=0.5, ls=(0, (2, 2)), zorder=0)
        lo, hi = ax.get_ylim()
        pad = (hi - lo) * 0.06
        segment_band(ax, b12, b23, lo - 2.2 * pad, pad)
        ax.set_ylim(lo, hi)
        ax.set_xlim(0, 1)
        ax.set_xticks([])
        ax.set_title(title, fontsize=7.5, color=INK, pad=2)
        if i == 0:
            ax.set_ylabel('log2 AKI / control')
    fig.savefig(out / 'f6_e.png', facecolor='white')
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------- f
def response(gene, curves):
    """Condition contrast of a gene's clipped curves: a spline dive is not drawn as a response."""
    out = {m: clipped(curves, m, gene) for m in MICE}
    return np.mean([out[m] for m in INJURED], 0) - np.mean([out[m] for m in CONTROLS], 0)


def panel_f(plt, genes, grid, curves, out, b12, b23):
    from matplotlib.colors import TwoSlopeNorm
    from matplotlib.patches import Rectangle
    from scipy.cluster.hierarchy import leaves_list, linkage

    rows = genes[genes.in_map].index.tolist()
    resp = np.stack([response(g, curves) for g in rows])
    seg = genes.loc[rows, [f'segment_view_S{s + 1}' for s in range(3)]].to_numpy()
    lim = 3.0
    order = leaves_list(linkage(np.clip(resp, -lim, lim)[:, ::4], 'average'))
    resp, seg = resp[order], seg[order]
    fig = figure(plt, 'f')
    norm = TwoSlopeNorm(0, -lim, lim)
    n = len(rows)
    ax = fig.add_axes([0.07, 0.07, 0.66, 0.89])
    ax.imshow(np.clip(resp, -lim, lim), aspect='auto', cmap='RdBu_r', norm=norm, extent=(0, 1, n, 0),
              interpolation='nearest', rasterized=True)
    for b in (b12, b23):
        ax.axvline(b, color='black', lw=0.4, alpha=0.4)
    ax.set_yticks([])
    ax.set_xticks([])
    segment_band(ax, b12, b23, n * 1.01, n * 0.018)
    ax.set_ylabel(f'{n} genes', labelpad=2)
    ax.set_title('along the PT', fontsize=7, color=MUTED, pad=2)
    ax2 = fig.add_axes([0.76, 0.07, 0.12, 0.89])
    ax2.imshow(np.clip(seg, -lim, lim), aspect='auto', cmap='RdBu_r', norm=norm, extent=(0, 3, n, 0),
               interpolation='nearest', rasterized=True)
    ax2.set_yticks([])
    ax2.set_xticks([])
    for k, name in enumerate(SEGS):
        ax2.add_patch(Rectangle((k, n * 1.01), 1, n * 0.018, color=COLOURS[name], lw=0, clip_on=False))
    ax2.set_title('S1 S2 S3', fontsize=7, color=MUTED, pad=2)
    cax = fig.add_axes([0.92, 0.35, 0.025, 0.30])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap='RdBu_r'), cax=cax)
    cb.set_ticks([-lim, 0, lim])
    cb.set_ticklabels([f'-{lim:g}', '0', f'+{lim:g}'])
    cb.outline.set_visible(False)
    cb.ax.tick_params(length=0, labelsize=6)
    cb.set_label('log2 AKI / control', fontsize=6, color=MUTED, labelpad=2)
    fig.savefig(out / 'f6_f.png', facecolor='white')
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------- g
G_PAIRS = [('Slc6a19', 'Slc6a13'), ('Slc39a5', 'Slc39a14'), ('Slc7a9', 'Slc7a13')]


def scaled(gene, curves, m, grid):
    """A mouse's curve min-max scaled over its sampled range, as the handoff position is defined (notebook 101)."""
    from pseudospace.handoffs import scaled_profile

    bx = curves[m]['bin_x'][gene]
    return scaled_profile(np.where(covered(grid, bx), curves[m]['curve'][gene], np.nan))


def panel_g(plt, grid, curves, hp, out, b12, b23):
    fig = figure(plt, 'g')
    for i, (e, l) in enumerate(G_PAIRS):
        ax = fig.add_axes([0.07 + i * 0.32, 0.13, 0.26, 0.70])
        r = hp[(hp.early == e) & (hp.late == l)].iloc[0]
        for mice, col in ((CONTROLS, CTRL), (INJURED, IR)):
            for g, ls in ((e, '-'), (l, (0, (3, 1.5)))):
                ax.plot(grid, np.nanmean([scaled(g, curves, m, grid) for m in mice], 0), color=col, lw=1.3, ls=ls)
            xs = np.nanmean([r[f'x_{m}'] for m in mice])
            ax.plot([xs], [-0.04], marker='^', ms=4.5, color=col, clip_on=False, zorder=5)
        ax.set_ylim(-0.06, 1.06)
        segment_band(ax, b12, b23, -0.19, 0.06)
        ax.set_xlim(0, 1)
        ax.set_xticks([])
        ax.set_yticks([0, 1])
        ax.set_title(f'{e} → {l}', fontsize=7, style='italic', color=INK, pad=2)
        if i == 0:
            ax.set_ylabel('scaled expression', labelpad=1)
    for j, (lab, ls) in enumerate((('early partner', '-'), ('late partner', (0, (3, 1.5))))):
        x0 = 0.66 + j * 0.17
        fig.lines.append(plt.Line2D([x0, x0 + 0.035], [0.965, 0.965], transform=fig.transFigure, color=MUTED, lw=1.1,
                                    ls=ls))
        fig.text(x0 + 0.042, 0.965, lab, fontsize=6, color=MUTED, va='center')
    fig.savefig(out / 'f6_g.png', facecolor='white')
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------- h
def panel_h(plt, hp, out, b12, b23):
    """Left: every handoff, crossing in control against AKI. Right: the shifted pairs, control to AKI."""
    fig = figure(plt, 'h')
    e = hp[hp.eligible].copy()
    e['cx'] = e[[f'x_{m}' for m in CONTROLS]].mean(1)
    e['ix'] = e[[f'x_{m}' for m in INJURED]].mean(1)
    sh = e.shifted.to_numpy(bool)
    ax = fig.add_axes([0.07, 0.17, 0.27, 0.78])
    for b in (b12, b23):
        ax.axvline(b, color='#E6E6E6', lw=0.5, zorder=0)
        ax.axhline(b, color='#E6E6E6', lw=0.5, zorder=0)
    ax.plot([0, 1], [0, 1], color='#BBBBBB', lw=0.6, ls=(0, (2, 2)), zorder=1)
    ax.scatter(e.cx[~sh], e.ix[~sh], s=8, color='#B8B8B8', lw=0, zorder=3)
    ax.scatter(e.cx[sh], e.ix[sh], s=13, color=IR, lw=0, zorder=4)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([0, 0.5, 1])
    ax.set_yticks([0, 0.5, 1])
    ax.set_aspect('equal')
    ax.set_xlabel('control', labelpad=1)
    ax.set_ylabel('AKI', labelpad=1)

    d = e[sh].sort_values('shift')
    ax2 = fig.add_axes([0.58, 0.17, 0.40, 0.78])
    for y, (_, r) in enumerate(d.iterrows()):
        ax2.plot([r.cx, r.ix], [y, y], color='#BBBBBB', lw=0.9, zorder=1)
        for m in CONTROLS + INJURED:
            ax2.scatter([r[f'x_{m}']], [y], s=6, facecolor='white', edgecolor=SECTION_COLOURS[m], lw=0.6, zorder=2)
        ax2.scatter([r.cx], [y], s=16, color=CTRL, lw=0, zorder=3)
        ax2.scatter([r.ix], [y], s=16, color=IR, lw=0, zorder=3)
    ax2.set_yticks(range(len(d)))
    ax2.set_yticklabels([f'{a} → {b}' for a, b in zip(d.early, d.late)], fontsize=6.3, style='italic')
    ax2.tick_params(axis='y', length=0, pad=2)
    ax2.set_ylim(-0.8, len(d) - 0.4)
    ax2.invert_yaxis()
    ax2.spines['left'].set_visible(False)
    for b in (b12, b23):
        ax2.axvline(b, color='#E6E6E6', lw=0.5, zorder=0)
    ax2.set_xlim(0, 1)
    ax2.set_xticks([0, 1])
    ax2.set_xticklabels(['start of PT', 'end'])
    segment_band(ax2, b12, b23, len(d) - 0.25, 0.32)
    fig.savefig(out / 'f6_h.png', facecolor='white')
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------- i
def set_profiles(root: Path, screen: str, calls: pd.DataFrame, genes: pd.DataFrame):
    """Condition response along the PT of each called set: the mean of its members' Δ(x), over the screened genes
    whose curves are supported in that segment (as in the screen)."""
    sets = json.loads((root.parents[1] / 'pt_package_aki/tables/gene_sets.json').read_text())
    z = np.load(root / screen / 'all_curves.npz', allow_pickle=True)
    names = list(z['genes'])
    delta = np.mean([z[f'curve|{m}'] for m in INJURED], 0) - np.mean([z[f'curve|{m}'] for m in CONTROLS], 0)
    pos = {g: i for i, g in enumerate(names)}
    screened = genes.universe & genes.neighbour_ok & genes.converged
    out = []
    for _, r in calls.iterrows():
        ok = screened & ~genes[f'unsupported_S{r.segment}']
        idx = [pos[g] for g in sets[r.set] if g in pos and ok.get(g, False)]
        out.append(delta[idx].mean(0))
    return z['grid'], np.array(out)


def panel_i(plt, root, screen, genes, out, b12, b23):
    pw = pd.read_csv(root / screen / 'pathways.csv')
    calls = pw[(pw.split == 'condition') & pw.call].reset_index(drop=True)
    grid, prof = set_profiles(root, screen, calls, genes)
    fig = figure(plt, 'i')
    core = [(0, b12 - 0.06), (b12 + 0.06, b23 - 0.06), (b23 + 0.06, 1)]
    for s in range(3):
        ax = fig.add_axes([0.06 + s * 0.325, 0.16, 0.27, 0.72])
        lo, hi = core[s]
        k = (grid >= lo) & (grid <= hi)
        rows = calls.index[calls.segment == s + 1]
        for sign, col in ((1, UP), (-1, DOWN)):
            sel = [i for i in rows if np.sign(calls.W[i]) == sign]
            if not sel:
                continue
            cen = prof[sel][:, k] - prof[sel][:, k].mean(1, keepdims=True)    # each set centred on its core mean
            for c in cen:
                ax.plot(grid[k], c, color=col, lw=0.4, alpha=0.35)
            ax.plot(grid[k], np.median(cen, 0), color=col, lw=1.8)
            ax.text(0.03 if sign > 0 else 0.97, 0.95 if sign > 0 else 0.05, str(len(sel)), transform=ax.transAxes,
                    fontsize=6.5, color=col, ha='left' if sign > 0 else 'right', va='top' if sign > 0 else 'bottom')
        ax.axhline(0, color=STEP, lw=1.4, zorder=0)
        ax.set_xlim(lo, hi)
        ax.set_ylim(-0.55, 0.55)
        ax.set_xticks([])
        ax.set_title(f'S{s + 1}', fontsize=7.5, color=COLOURS[SEGS[s]], pad=2)
        from matplotlib.patches import Rectangle
        ax.add_patch(Rectangle((lo, -0.66), hi - lo, 0.05, color=COLOURS[SEGS[s]], lw=0, clip_on=False))
        if s == 0:
            ax.set_ylabel('log2 AKI / control,\ncentred in the segment', labelpad=1)
        else:
            ax.set_yticklabels([])
    fig.savefig(out / 'f6_i.png', facecolor='white')
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------- main
def main(argv=None):
    args = parse_args(argv)
    root = Path(args.results_root) / args.out
    rroot = Path(args.results_root)
    out = root / 'panels'
    out.mkdir(parents=True, exist_ok=True)
    plt = style()
    have_screen = (root / args.screen / 'genes.csv').exists()
    if have_screen:
        genes, grid, curves, programs = load(root, args.screen)
        hp = pd.read_csv(root / args.screen / 'handoffs.csv')
        b12, b23 = float(genes.b12.median()), float(genes.b23.median())
    else:
        b12, b23 = 0.29, 0.675
    if 'a' in args.panels:
        panel_a(plt, rroot, out)
    if 'b' in args.panels:
        panel_b(plt, rroot, out, b12, b23)
    if 'c' in args.panels:
        panel_c(plt, root, out)
    if not have_screen:
        return
    if 'd' in args.panels:
        panel_d(plt, genes, grid, curves, out)
    if 'e' in args.panels:
        panel_e(plt, grid, programs, out, b12, b23)
    if 'f' in args.panels:
        panel_f(plt, genes, grid, curves, out, b12, b23)
    if 'g' in args.panels:
        panel_g(plt, grid, curves, hp, out, b12, b23)
    if 'h' in args.panels:
        panel_h(plt, hp, out, b12, b23)
    if 'i' in args.panels:
        panel_i(plt, root, args.screen, genes, out, b12, b23)


if __name__ == '__main__':
    main()
