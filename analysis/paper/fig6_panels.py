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
        'f': (560, 312), 'g': (580, 312), 'h': (1160, 412), 'i': (1160, 322)}


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


# ---------------------------------------------------------------------------------------------------- h
# h, poster layout: one row per transporter pair, control and AKI side by side, the switch marked in both
G_PAIRS = [('Slc16a1', 'Slc16a4', 'Monocarboxylates'), ('Slc6a19', 'Slc6a13', 'Amino acids'),
           ('Slc3a2', 'Slc3a1', 'Amino-acid carrier'), ('Slc39a5', 'Slc39a14', 'Zinc'), ('Slc7a9', 'Slc7a13', 'Cystine')]
PROTEIN = {'Slc16a1': 'MCT1', 'Slc16a4': 'MCT4', 'Slc6a19': 'B0AT1', 'Slc6a13': 'GAT2', 'Slc3a2': '4F2hc',
           'Slc3a1': 'rBAT', 'Slc39a5': 'ZIP5', 'Slc39a14': 'ZIP14', 'Slc7a9': 'b0,+AT', 'Slc7a13': 'AGT1',
           'Slc39a8': 'ZIP8', 'Slc16a9': 'MCT9', 'Slc16a12': 'MCT12', 'Slc22a6': 'OAT1', 'Slc22a2': 'OCT2',
           'Slc46a3': 'SLC46A3', 'Slc46a1': 'PCFT', 'Slc6a8': 'CRT1', 'Slc5a9': 'SGLT4', 'Slc5a3': 'SMIT1'}
EARLY = {'control': '#3A8DD0', 'AKI': '#F08A2C'}
LATE = {'control': '#1B3A5C', 'AKI': '#8A3A12'}
BAND = {'PTS1': '#DCEAF6', 'PTS2': '#E4EFDF', 'PTS3': '#F2EADF'}


def scaled(gene, curves, m, grid):
    """A mouse's curve min-max scaled over its sampled range, as the handoff position is defined (notebook 101)."""
    from pseudospace.handoffs import scaled_profile

    bx = curves[m]['bin_x'][gene]
    return scaled_profile(np.where(covered(grid, bx), curves[m]['curve'][gene], np.nan))


def where_in(x, b12, b23):
    seg, lo, hi = (1, 0, b12) if x < b12 else (2, b12, b23) if x < b23 else (3, b23, 1)
    part = ('early', 'middle', 'late')[min(int(3 * (x - lo) / (hi - lo)), 2)]
    return f'{part} S{seg}'


def panel_h(plt, grid, curves, hp, out, b12, b23):
    from matplotlib.patches import FancyBboxPatch

    fig = figure(plt, 'h')
    n = len(G_PAIRS)
    top, bottom, gap = 0.86, 0.09, 0.025
    h = (top - bottom - gap * (n - 1)) / n
    cols = {'control': (0.20, 0.33), 'AKI': (0.57, 0.33)}
    edges = ((0, b12), (b12, b23), (b23, 1))
    for name, (x0, w) in cols.items():
        fig.text(x0 + w / 2, 0.975, name, ha='center', va='top', fontsize=9, weight='bold',
                 color=CTRL if name == 'control' else IR)
        for (lo, hi), seg in zip(edges, SEGS):
            fig.text(x0 + w * (lo + hi) / 2, 0.895, f' {seg} ', ha='center', va='bottom', fontsize=6.5, weight='bold',
                     color='white', bbox=dict(boxstyle='round,pad=0.25', fc=COLOURS[seg], ec='none'))
    fig.lines.append(plt.Line2D([0.005, 0.03], [0.965, 0.965], transform=fig.transFigure, color=CTRL, lw=1.2,
                                ls=(0, (3, 1.5))))
    fig.text(0.035, 0.965, 'switch in control', fontsize=6.5, color=MUTED, va='center')
    fig.lines.append(plt.Line2D([0.005, 0.03], [0.925, 0.925], transform=fig.transFigure, color=IR, lw=1.2,
                                ls=(0, (1, 1.2))))
    fig.text(0.035, 0.925, 'switch in AKI', fontsize=6.5, color=MUTED, va='center')
    fig.text(0.955, 0.895, 'level in AKI\n(log2 / control)', ha='center', va='bottom', fontsize=6, color=MUTED)
    for r, (e, l, substrate) in enumerate(G_PAIRS):
        y0 = top - (r + 1) * h - r * gap
        row = hp[(hp.early == e) & (hp.late == l)].iloc[0]
        xc = np.nanmean([row[f'x_{m}'] for m in CONTROLS])
        xi = np.nanmean([row[f'x_{m}'] for m in INJURED])
        fig.text(0.005, y0 + h * 0.80, substrate, fontsize=7.5, weight='bold', color=INK, va='center')
        t = fig.text(0.005, y0 + h * 0.50, e, fontsize=7, style='italic', weight='bold', color=EARLY['control'],
                     va='center')
        right = t.get_window_extent(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted()).x1
        fig.text(right + 0.003, y0 + h * 0.50, f'→ {l}', fontsize=7, style='italic', weight='bold',
                 color=LATE['control'], va='center')
        fig.text(0.005, y0 + h * 0.20, f'{PROTEIN.get(e, e)} → {PROTEIN.get(l, l)}', fontsize=6.3, color=MUTED,
                 va='center')
        for name, (x0, w) in cols.items():
            ax = fig.add_axes([x0, y0, w, h])
            for (lo, hi), seg in zip(edges, SEGS):
                ax.axvspan(lo, hi, color=BAND[seg], lw=0, zorder=0)
            mice = CONTROLS if name == 'control' else INJURED
            for gene, col in ((e, EARLY[name]), (l, LATE[name])):
                per = [scaled(gene, curves, m, grid) for m in mice]
                for c in per:
                    ax.plot(grid, c, color=col, lw=0.5, alpha=0.6, zorder=2)
                mean = np.nanmean(per, 0)
                ax.fill_between(grid, 0, mean, color=col, alpha=0.12, lw=0, zorder=1)
                ax.plot(grid, mean, color=col, lw=1.6, zorder=3)
            ax.axvline(xc, color=CTRL, lw=1.1, ls=(0, (3, 1.5)), zorder=4)
            if name == 'AKI':
                ax.axvline(xi, color=IR, lw=1.3, ls=(0, (1, 1.2)), zorder=4)
            else:
                ax.text(xc + 0.012, 1.0, where_in(xc, b12, b23), fontsize=6, weight='bold', color=INK, va='top',
                        bbox=dict(boxstyle='square,pad=0.1', fc='white', ec='none', alpha=0.8), zorder=5)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1.08)
            ax.set_xticks([])
            ax.set_yticks([0, 1] if name == 'control' else [])
            for sp in ('top', 'right'):
                ax.spines[sp].set_visible(False)
            if r == n - 1:
                ax.set_xlabel('start → end of PT', labelpad=2)
            if name == 'control' and r == n // 2:
                ax.set_ylabel('scaled expression', labelpad=2)
        bx = fig.add_axes([0.925, y0 + h * 0.15, 0.06, h * 0.7])
        lv = [row.level_early, row.level_late]
        bx.barh([1, 0], lv, color=[EARLY['AKI'], LATE['AKI']], height=0.7)
        bx.axvline(0, color='#999999', lw=0.6)
        bx.set_xlim(-5, 3)
        bx.set_yticks([])
        bx.set_xticks([-4, 0, 2] if r == n - 1 else [])
        bx.tick_params(labelsize=6)
        for sp in ('top', 'right', 'left'):
            bx.spines[sp].set_visible(False)
    fig.savefig(out / 'f6_h.png', facecolor='white')
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------- g
def panel_g(plt, hp, out, b12, b23):
    """The pairs that move: switch position in control (blue) and AKI (orange), each mouse as an open dot."""
    fig = figure(plt, 'g')
    e = hp[hp.eligible & hp.shifted].copy()
    e['cx'] = e[[f'x_{m}' for m in CONTROLS]].mean(1)
    e['ix'] = e[[f'x_{m}' for m in INJURED]].mean(1)
    d = e.sort_values('shift')
    ax = fig.add_axes([0.40, 0.12, 0.57, 0.84])
    for (lo, hi), seg in zip(((0, b12), (b12, b23), (b23, 1)), SEGS):
        ax.axvspan(lo, hi, color=BAND[seg], lw=0, zorder=0)
    for y, (_, r) in enumerate(d.iterrows()):
        ax.annotate('', xy=(r.ix, y), xytext=(r.cx, y), zorder=2,
                    arrowprops=dict(arrowstyle='-|>', color='#888888', lw=1.0, shrinkA=4, shrinkB=4, mutation_scale=7))
        for m in CONTROLS + INJURED:
            ax.scatter([r[f'x_{m}']], [y], s=9, facecolor='white', edgecolor=SECTION_COLOURS[m], lw=0.7, zorder=3)
        ax.scatter([r.cx], [y], s=26, color=CTRL, lw=0, zorder=4)
        ax.scatter([r.ix], [y], s=26, color=IR, lw=0, zorder=4)
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels([f'{a} → {b}  ({PROTEIN.get(a, a)} → {PROTEIN.get(b, b)})' for a, b in zip(d.early, d.late)],
                       fontsize=6.3)
    ax.tick_params(axis='y', length=0, pad=3)
    ax.set_ylim(len(d) - 0.4, -0.6)
    ax.spines['left'].set_visible(False)
    ax.set_xlim(0, 1)
    ax.set_xticks([0, b12, b23, 1])
    ax.set_xticklabels(['start', 'S1 | S2', 'S2 | S3', 'end of PT'], fontsize=6)
    fig.savefig(out / 'f6_g.png', facecolor='white')
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------- i
# poster layout: per pathway, its S1/S2/S3 means beside the PT drawn as a tube coloured by the response along it
I_SETS = [('Reactome_2022::S Phase R-HSA-69242', 'S phase'),
          ('Reactome_2022::Integrin Cell Surface Interactions R-HSA-216083', 'Integrin cell-surface interactions'),
          ('Reactome_2022::Cholesterol Biosynthesis R-HSA-191273', 'Cholesterol biosynthesis'),
          ('KEGG_2019_Mouse::Oxidative phosphorylation', 'Oxidative phosphorylation')]
# one per theme of the calls, chosen after the screen (presentation choice, recorded in the protocol doc)


def nephron_path(b12, b23, n=400):
    """Points along a drawn PT (glomerulus top right; S1 convoluted, a U-turn, S2 convoluted, S3 straight down into
    the outer medulla) and the PT position x of each point."""
    t = np.linspace(0, 1, n)
    s1 = np.c_[0.80 - 0.56 * t, 0.76 + 0.035 * np.sin(2 * np.pi * 3 * t)]
    a = np.linspace(np.pi / 2, 3 * np.pi / 2, n // 3)
    turn = np.c_[0.24 + 0.12 * np.cos(a), 0.64 + 0.12 * np.sin(a)]
    s2 = np.c_[0.24 + 0.46 * t, 0.52 + 0.03 * np.sin(2 * np.pi * 3 * t)]
    diag = np.c_[np.linspace(0.70, 0.76, n // 4), np.linspace(0.52, 0.42, n // 4)]
    s3 = np.c_[np.full(n // 2, 0.76), np.linspace(0.42, 0.05, n // 2)]
    out, xs = [], []
    for pieces, (lo, hi) in (((s1, turn), (0, b12)), ((s2, diag), (b12, b23)), ((s3,), (b23, 1))):
        pts = np.vstack(pieces)
        d = np.r_[0, np.cumsum(np.hypot(*np.diff(pts, axis=0).T))]
        out.append(pts)
        xs.append(lo + (hi - lo) * d / d[-1])
    return np.vstack(out), np.concatenate(xs)


def draw_nephron(ax, prof, grid, b12, b23, cmap):
    from matplotlib.collections import LineCollection
    from matplotlib.patches import Circle, Rectangle

    ax.add_patch(Rectangle((0, 0.30), 1, 0.72, color='#FBF6EC', lw=0, zorder=0))
    ax.add_patch(Rectangle((0, 0), 1, 0.30, color='#F1E7DA', lw=0, zorder=0))
    ax.add_patch(Circle((0.88, 0.84), 0.075, fc='none', ec='#999999', lw=1.0, zorder=1))
    ax.add_patch(Circle((0.88, 0.84), 0.05, fc='#D5D8DD', ec='#999999', lw=0.6, zorder=1))
    pts, xs = nephron_path(b12, b23)
    seg = np.stack([pts[:-1], pts[1:]], 1)
    mid = (xs[:-1] + xs[1:]) / 2
    outline = [COLOURS['PTS1'] if v < b12 else COLOURS['PTS2'] if v < b23 else COLOURS['PTS3'] for v in mid]
    ax.add_collection(LineCollection(seg, colors=outline, linewidths=7.5, capstyle='round', zorder=2))
    v = np.interp(mid, grid, prof)
    a = np.abs(v)
    lo, hi = a.min(), a.max()
    ax.add_collection(LineCollection(seg, colors=cmap(0.12 + 0.88 * (a - lo) / max(hi - lo, 1e-9)), linewidths=5,
                                     capstyle='round', zorder=3))
    for b in (b12, b23):          # segment boundaries across the tube
        k = np.argmin(np.abs(xs - b))
        d = pts[min(k + 1, len(pts) - 1)] - pts[max(k - 1, 0)]
        nrm = np.array([-d[1], d[0]]) / np.hypot(*d) * 0.045
        ax.plot(*np.c_[pts[k] - nrm, pts[k] + nrm], color=INK, lw=1.0, zorder=4)
    k = np.argmax(np.abs(np.interp(xs, grid, prof)))
    ax.plot(*pts[k], marker='o', ms=5.5, mfc='white', mec=INK, mew=1.2, zorder=5)
    for name, xy in (('PTS1', (0.47, 0.88)), ('PTS2', (0.47, 0.40)), ('PTS3', (0.82, 0.20))):
        ax.text(*xy, name, fontsize=6.5, weight='bold', color=INK, ha='left' if name == 'PTS3' else 'center',
                va='center', zorder=6)
    ax.text(0.02, 0.95, 'cortex', fontsize=6, style='italic', color=MUTED, va='top')
    ax.text(0.02, 0.12, 'outer\nmedulla', fontsize=6, style='italic', color=MUTED)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.0)
    ax.axis('off')
    return v.min(), v.max()


def panel_i(plt, root, screen, genes, out, b12, b23):
    from matplotlib.colors import TwoSlopeNorm

    pw = pd.read_csv(root / screen / 'pathways.csv')
    calls = pw[(pw.split == 'condition') & pw.call]
    sets = json.loads((root.parents[1] / 'pt_package_aki/tables/gene_sets.json').read_text())
    z = np.load(root / screen / 'all_curves.npz', allow_pickle=True)
    grid = z['grid']
    pos = {g: i for i, g in enumerate(z['genes'])}
    delta = np.mean([z[f'curve|{m}'] for m in INJURED], 0) - np.mean([z[f'curve|{m}'] for m in CONTROLS], 0)
    screened = genes.universe & genes.neighbour_ok & genes.converged
    supported = screened & ~genes[[f'unsupported_S{s + 1}' for s in range(3)]].any(axis=1)
    norm = TwoSlopeNorm(0, -1.5, 1.5)
    fig = figure(plt, 'i')
    n = len(I_SETS)
    for k, (key, title) in enumerate(I_SETS):
        x0, w = 0.01 + k / n, 1 / n - 0.02
        mem = [g for g in sets[key] if supported.get(g, False)]
        prof = delta[[pos[g] for g in mem]].mean(0)
        segv = genes.loc[mem, [f'segment_view_S{s + 1}' for s in range(3)]].mean().to_numpy()
        up = segv.mean() > 0
        mine = calls[calls.set == key].sort_values('p').iloc[0]
        called = set(calls.loc[calls.set == key, 'segment'])
        wcol = genes.loc[mem, f'W_condition_S{mine.segment}']
        top3 = (wcol * np.sign(mine.W)).sort_values(ascending=False).index[:3]
        fig.text(x0, 0.965, f'{title} {"↑" if up else "↓"}', fontsize=8, weight='bold', color=INK, va='top')
        fig.text(x0, 0.875, ', '.join(top3), fontsize=6.5, style='italic', color=MUTED, va='top')
        fig.text(x0, 0.79, 'S1/S2/S3 means', fontsize=6.3, weight='bold', color=INK, va='center')
        for s in range(3):
            c = plt.get_cmap('RdBu_r')(norm(segv[s]))
            dark = abs(segv[s]) > 0.9
            fig.text(x0 + w * (0.42 + 0.20 * s), 0.79, f'S{s + 1} {segv[s]:+.2f}', fontsize=6.3, va='center',
                     color='white' if dark else INK, weight='bold' if s + 1 in called else 'normal',
                     bbox=dict(boxstyle='round,pad=0.25', fc=c, ec=INK if s + 1 in called else 'none', lw=0.8))
        ax = fig.add_axes([x0, 0.02, w, 0.66])
        lo, hi = draw_nephron(ax, prof, grid, b12, b23, plt.get_cmap('Reds' if up else 'Blues'))
        fig.text(x0, 0.715, f'pseudostructure  {lo:+.2f} to {hi:+.2f}', fontsize=6.3, color=INK, va='center')
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
        panel_g(plt, hp, out, b12, b23)
    if 'h' in args.panels:
        panel_h(plt, grid, curves, hp, out, b12, b23)
    if 'i' in args.panels:
        panel_i(plt, root, args.screen, genes, out, b12, b23)


if __name__ == '__main__':
    main()
