"""Fig. 6 panels (main) and their Extended Data companions, drawn from fig6_placement.py and fig6_screen.py.

Main figure (make_fig6.py):
 programs       cell-cycle program and PT transporters: AKI response along the PT against the segment means
 pathways       distinct pathways whose AKI response changes inside a segment: response along the PT, S1/S2/S3
                means, and the difference
 handoff_rows   five handoffs in the poster layout: early and late partner, control and IR, switches, AKI levels
Extended Data (make_ed_fig_aki.py):
 maps           one control and one IR section, PT tubules coloured by position on the shared axis
 coverage       where each section's PT tubules sit on the axis (density per section)
 placement      placement check: positions on the injury-blind axis against the reference, per section and label
 genes          injury and transport genes along the PT, control and IR, each mouse
 response_map   every gene with a large AKI response: response along the PT beside its segment view
 handoff_moves  the transporter handoffs that move: switch position in control and IR, per mouse

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
SIZE = {'programs': (420, 452), 'pathways': (720, 452), 'handoff_rows': (1160, 412),
        'maps': (360, 242), 'coverage': (380, 242), 'placement': (380, 242), 'genes': (580, 312),
        'response_map': (560, 312), 'handoff_moves': (580, 312)}
MAIN = ['programs', 'pathways', 'handoff_rows']
EXTENDED = ['maps', 'coverage', 'placement', 'genes', 'response_map', 'handoff_moves']


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--results-root', default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT', str(REPO / 'results')))
    p.add_argument('--out', default='paper_figures/fig6', help='relative to the results root')
    p.add_argument('--panels', nargs='+', default=None, help='panel names (default: all)')
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
    fig = figure(plt, 'maps')
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
    fig.savefig(out / 'f6_maps.png', facecolor='white', dpi=300)
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------- b
def panel_b(plt, root, out, b12, b23):
    from scipy.stats import gaussian_kde

    pt = reference_x(root)
    fig = figure(plt, 'coverage')
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
    fig.savefig(out / 'f6_coverage.png', facecolor='white')
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------- c
def panel_c(plt, root, out):
    pl = pd.read_csv(root / 'placement/placement.csv')
    fig = figure(plt, 'placement')
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
    fig.savefig(out / 'f6_placement.png', facecolor='white')
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
    fig = figure(plt, 'genes')
    for i, g in enumerate(D_GENES):
        ax = fig.add_axes([0.07 + i * 0.235, 0.12, 0.19, 0.76])
        gene_levels(ax, g, genes.loc[g], grid, curves, level_ylim(g, curves))
        ax.set_title(g, fontsize=7.5, style='italic', color=INK, pad=2)
        if i == 0:
            ax.set_ylabel('log2 counts per 10k')
    fig.savefig(out / 'f6_genes.png', facecolor='white')
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------- e
E_PROGRAMS = [('C', 'cell-cycle program'), ('T_PT', 'PT transporters')]


def staircase(b12, b23, values):
    return np.array([0, b12, b12, b23, b23, 1]), np.repeat(values, 2)


def panel_e(plt, grid, programs, out, b12, b23):
    fig = figure(plt, 'programs')
    for i, (key, title) in enumerate(E_PROGRAMS):
        ax = fig.add_axes([0.16, 0.56 - i * 0.48, 0.80, 0.36])
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
        ax.set_ylabel('log2 AKI / control')
    fig.savefig(out / 'f6_programs.png', facecolor='white')
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
    fig = figure(plt, 'response_map')
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
    fig.savefig(out / 'f6_response_map.png', facecolor='white')
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

    fig = figure(plt, 'handoff_rows')
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
    fig.savefig(out / 'f6_handoff_rows.png', facecolor='white')
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------- g
def panel_g(plt, hp, out, b12, b23):
    """The pairs that move: switch position in control (blue) and AKI (orange), each mouse as an open dot."""
    fig = figure(plt, 'handoff_moves')
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
    fig.savefig(out / 'f6_handoff_moves.png', facecolor='white')
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------- pathways
LIBRARY_ORDER = ('MSigDB_Hallmark_2020', 'Reactome_2022', 'KEGG_2019_Mouse')
JACCARD = 0.3      # calls in one segment sharing more of their screened genes than this are drawn once


DISEASE_WORDS = ('disease', 'cancer', 'carcinoma', 'leukemia', 'infection', 'tuberculosis', 'diabetic', 'hepatitis')


def pathway_name(key, width=48):
    import re

    lib, name = key.split('::', 1)
    name = ('Hallmark ' if lib.startswith('MSigDB') else '') + re.sub(r'\s+R-HSA-\d+$', '', name)
    return name if len(name) <= width else name[:name.rfind(' ', 0, width)].rstrip(',') + ' …'


def distinct_calls(calls, sets, screened):
    """One row per family of overlapping calls within a segment: candidates in library order (Hallmark, Reactome,
    KEGG), then by p; a call is kept unless it shares >= JACCARD of its genes with one already kept."""
    # KEGG "Human Diseases" sets (Alzheimer disease, small cell lung cancer, ...) are not drawn: in kidney PT they
    # stand in for their respiratory-chain or cell-cycle members (presentation choice; they stay in pathways.csv)
    calls = calls[~calls.set.str.lower().str.contains('|'.join(DISEASE_WORDS))]
    calls = calls.assign(lib=calls.set.str.split('::').str[0].map({l: i for i, l in enumerate(LIBRARY_ORDER)}))
    kept = []
    for seg, grp in calls.groupby('segment'):
        mine = []
        for _, r in grp.sort_values(['lib', 'p', 'n'], ascending=[True, True, False]).iterrows():
            m = {g for g in sets[r.set] if screened.get(g, False)}
            if all(len(m & o) / len(m | o) < JACCARD for o in mine):
                mine.append(m)
                kept.append(r)
    kept = pd.DataFrame(kept)
    # a pathway called in two segments is one row, ordered by its first called segment
    out = kept.groupby('set').agg(segments=('segment', lambda v: sorted(set(v))), first=('segment', 'min'),
                                  W=('W', 'first')).reset_index()
    return out.sort_values(['first', 'W']).reset_index(drop=True)


def panel_pathways(plt, root, screen, genes, out, b12, b23):
    """Distinct pathway calls x position: the AKI response along the PT, the same pathways' S1/S2/S3 means, and the
    response centred within each segment (the change the segment means cannot show). Called segments are outlined."""
    from matplotlib.colors import TwoSlopeNorm
    from matplotlib.patches import Rectangle

    pw = pd.read_csv(root / screen / 'pathways.csv')
    calls = pw[(pw.split == 'condition') & pw.call]
    sets = json.loads((root.parents[1] / 'pt_package_aki/tables/gene_sets.json').read_text())
    screened = genes.universe & genes.neighbour_ok & genes.converged
    rows = distinct_calls(calls, sets, screened)
    z = np.load(root / screen / 'all_curves.npz', allow_pickle=True)
    grid = z['grid']
    pos = {g: i for i, g in enumerate(z['genes'])}
    delta = np.mean([z[f'curve|{m}'] for m in INJURED], 0) - np.mean([z[f'curve|{m}'] for m in CONTROLS], 0)
    resp, segv = [], []
    for _, r in rows.iterrows():
        ok = screened & ~genes[[f'unsupported_S{s}' for s in r.segments]].any(axis=1)
        mem = [g for g in sets[r.set] if ok.get(g, False)]
        resp.append(delta[[pos[g] for g in mem]].mean(0))
        segv.append(genes.loc[mem, [f'segment_view_S{s + 1}' for s in range(3)]].mean().to_numpy())
    resp, segv = np.array(resp), np.array(segv)
    # within each segment core, the response minus its own mean over that core: only the change inside segments
    # (the transition windows, where S1/S2 or S2/S3 tubules mix, are left blank)
    within = np.full_like(resp, np.nan)
    for lo, hi in ((0, b12 - 0.06), (b12 + 0.06, b23 - 0.06), (b23 + 0.06, 1.0)):
        k = (grid >= lo) & (grid <= hi)
        within[:, k] = resp[:, k] - resp[:, k].mean(1, keepdims=True)
    n = len(rows)
    fig = figure(plt, 'pathways')
    top, bottom = 0.90, 0.07
    lim, mlim = 1.5, 0.5
    blocks = [('along the PT', resp, TwoSlopeNorm(0, -lim, lim), 'RdBu_r', 0.36, 0.30),
              ('S1/S2/S3 means', segv, TwoSlopeNorm(0, -lim, lim), 'RdBu_r', 0.675, 0.075),
              ('within each segment', within, TwoSlopeNorm(0, -mlim, mlim), 'PuOr_r', 0.77, 0.20)]
    for title, data, norm, cmap, x0, w in blocks:
        ax = fig.add_axes([x0, bottom, w, top - bottom])
        three = data.shape[1] == 3
        cm = plt.get_cmap(cmap).copy()
        cm.set_bad('#F3F3F3')
        ax.imshow(np.ma.masked_invalid(data), aspect='auto', cmap=cm, norm=norm, interpolation='nearest',
                  extent=(0, 3, n, 0) if three else (0, 1, n, 0))
        ax.set_xticks([])
        ax.set_yticks([])
        if three:
            for s in range(3):
                ax.add_patch(Rectangle((s, n + 0.25), 1, 0.45, color=COLOURS[SEGS[s]], lw=0, clip_on=False))
        else:
            for b in (b12, b23):
                ax.axvline(b, color='#333333', lw=0.5)
            segment_band(ax, b12, b23, n + 0.25, 0.45)
            cores = [(0, b12 - 0.06), (b12 + 0.06, b23 - 0.06), (b23 + 0.06, 1)]
            for y, r in rows.iterrows():
                for seg in r.segments:
                    lo, hi = cores[int(seg) - 1]
                    ax.add_patch(Rectangle((lo, y + 0.06), hi - lo, 0.88, fc='none', ec=INK, lw=0.7))
        ax.set_ylim(n + 0.75, 0)
        ax.set_title(title, fontsize=7.5, color=MUTED, pad=2)
        for sp in ax.spines.values():
            sp.set_visible(False)
        cax = fig.add_axes([x0 + w * 0.15, 0.955, w * 0.7, 0.012]) if not three else None
        if cax is not None:
            cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation='horizontal')
            cb.set_ticks([norm.vmin, 0, norm.vmax])
            cb.set_ticklabels([f'{norm.vmin:+g}', '0', f'{norm.vmax:+g}'])
            cb.ax.tick_params(length=0, labelsize=5.5, pad=1)
            cb.ax.xaxis.set_ticks_position('top')
            cb.outline.set_visible(False)
    for y, r in rows.iterrows():
        yy = top - (y + 0.5) / (n + 0.75) * (top - bottom)
        fig.text(0.345, yy, pathway_name(r.set), fontsize=6.2, ha='right', va='center', color=INK)
    for seg in (1, 2, 3):
        idx = np.flatnonzero(rows['first'].to_numpy() == seg)
        if len(idx):
            y0 = top - idx[0] / (n + 0.75) * (top - bottom)
            y1 = top - (idx[-1] + 1) / (n + 0.75) * (top - bottom)
            fig.lines.append(plt.Line2D([0.352, 0.352], [y1 + 0.003, y0 - 0.003], transform=fig.transFigure,
                                        color=COLOURS[SEGS[seg - 1]], lw=3))
    rows.assign(name=rows.set.map(pathway_name)).to_csv(out / 'f6_pathways_rows.csv', index=False)
    fig.savefig(out / 'f6_pathways.png', facecolor='white')
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------- main
def main(argv=None):
    args = parse_args(argv)
    root = Path(args.results_root) / args.out
    rroot = Path(args.results_root)
    out = root / 'panels'
    out.mkdir(parents=True, exist_ok=True)
    plt = style()
    args.panels = args.panels or MAIN + EXTENDED
    have_screen = (root / args.screen / 'genes.csv').exists()
    if have_screen:
        genes, grid, curves, programs = load(root, args.screen)
        hp = pd.read_csv(root / args.screen / 'handoffs.csv')
        b12, b23 = float(genes.b12.median()), float(genes.b23.median())
    else:
        b12, b23 = 0.29, 0.675
    if 'maps' in args.panels:
        panel_a(plt, rroot, out)
    if 'coverage' in args.panels:
        panel_b(plt, rroot, out, b12, b23)
    if 'placement' in args.panels:
        panel_c(plt, root, out)
    if not have_screen:
        return
    if 'genes' in args.panels:
        panel_d(plt, genes, grid, curves, out)
    if 'programs' in args.panels:
        panel_e(plt, grid, programs, out, b12, b23)
    if 'response_map' in args.panels:
        panel_f(plt, genes, grid, curves, out, b12, b23)
    if 'handoff_moves' in args.panels:
        panel_g(plt, hp, out, b12, b23)
    if 'handoff_rows' in args.panels:
        panel_h(plt, grid, curves, hp, out, b12, b23)
    if 'pathways' in args.panels:
        panel_pathways(plt, root, args.screen, genes, out, b12, b23)


if __name__ == '__main__':
    main()
