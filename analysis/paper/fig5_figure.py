"""Fig. 5 drawn as one journal-style canvas (180 mm wide), restyled 2026-11-11 after reading published continuous-axis
methods figures (GASTON, Genes2Genes, STcompare, STalign, Halpern 2017 liver zonation):
 a  one mouse and one human section, PT tubules coloured by pseudostructure
 b  Slc5a2 (SGLT2) expression in the same tubules
 c  Slc5a2 against pseudostructure: binned means (points) and fitted curves (lines) in log2(CP10k + 0.5)
 d  the same for five more transporters (Slc2a2, GLUT2, is the unmoved comparator); the sixth slot is the section key
 e  transporter pairs in the poster style: mean log2(1 + CP10k) curves per species, axis from 0, no switch lines
 f  each shifted switch as an arrow from its mouse to its human position, from the curves and from straight lines
    through the S1/S2/S3 cluster means (the orthology panel is Extended Data, fig5_panels.panel_orthology)
Panel letters only; claims live in the legend and text. Gene symbols are italic (mouse Slc5a2, human SLC5A2).
Writes fig5.png, fig5.pdf and fig5_wireframe.png to <results-root>/paper_figures/fig5/.
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
from fig5_panels import (EARLY, H_PAIRS, HUMAN, HUMAN_C, LATE, MOUSE, MOUSE_C, PROTEIN, SECTIONS,  # noqa: E402
                         XRANGE, handoff_range, load, reference_x)
from fig6_panels import BAND, SEGS  # noqa: E402

MM = 1 / 25.4
INK, GREY = '#222222', '#8A8A8A'
TRIM = 0.05                      # spline ends overshoot (question-bank.md); curves and bins shown on [0.05, 0.95]
# Abcc2 (MRP2) is a strong call but its two human curves disagree in shape; it is shown in Extended Data
GENES_D = ['Slc2a2', 'Slc22a6', 'Slc22a8', 'Slc22a2', 'Slc22a12']
EXPR_CMAP = None
# main e: the shifted pairs whose handover is visible in counts; urate (SLC22A23 -> URAT1) is in the Extended Data
# scaled rows only, because URAT1 is the more abundant transcript along the whole PT in both species
E_PAIRS = H_PAIRS[:3]
SHADE = {'Ctrl1A2': MOUSE_C, 'Ctrl1A4': '#7FA8D8', 'HUK1_COR1': HUMAN_C, 'HUK1_MED1': '#E58BBE'}


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--results-root', default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT', str(REPO / 'results')))
    return p.parse_args(argv)


def style():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family': ['Nimbus Sans', 'Liberation Sans', 'DejaVu Sans'], 'font.size': 6,
                         'axes.linewidth': 0.5, 'axes.edgecolor': INK, 'axes.labelcolor': INK, 'axes.labelsize': 6,
                         'xtick.color': INK, 'ytick.color': INK, 'xtick.labelsize': 5.5, 'ytick.labelsize': 5.5,
                         'xtick.major.width': 0.5, 'ytick.major.width': 0.5, 'xtick.major.size': 2,
                         'ytick.major.size': 2, 'axes.spines.top': False, 'axes.spines.right': False,
                         'pdf.fonttype': 42, 'savefig.dpi': 400})
    return plt


def sym(g, sp):
    return g if sp == 'mouse' else g.upper()


def segments(ax, b12, b23, top=True):
    for b in (b12, b23):
        ax.axvline(b, color=GREY, lw=0.5, ls=(0, (2, 2)), zorder=0)
    if top:
        for x, s in zip((b12 / 2, (b12 + b23) / 2, (b23 + 1) / 2), ('S1', 'S2', 'S3')):
            ax.text(x, 1.0, s, transform=ax.get_xaxis_transform(), ha='center', va='bottom', fontsize=5.5, color=GREY)


def log_with_bins(gene, m, curves, pos, grid, lo_x=TRIM, hi_x=1 - TRIM):
    """Curve and bin means of one section in log2(counts per 10k + 0.5), the floor the question bank uses (the
    screen's curves are unfloored log2 per 10k), on [lo_x, hi_x] and the section's bin range."""
    i = pos[gene]
    c, bx, bm = (curves[m][k][i].astype(float) for k in ('curve', 'bin_x', 'bins'))
    fin = np.isfinite(bm)
    k = (grid >= max(lo_x, bx[fin].min())) & (grid <= min(hi_x, bx[fin].max()))
    keep = fin & (bx >= lo_x) & (bx <= hi_x)
    floor = lambda v: np.log2(2 ** v + 0.5)  # noqa: E731
    return np.where(k, floor(c), np.nan), bx[keep], floor(bm[keep])


def curve_axis(ax, gene, curves, pos, grid, b12, b23, title=True, ylabel=False, xlabel=False, seglabels=True, inside=False):
    """Expression in log2(CP10k + 0.5) on one axis for all four sections. Levels are not comparable across species
    (different probe panels); positions are."""
    vals = []
    for m in SECTIONS:
        c, bx, bm = log_with_bins(gene, m, curves, pos, grid)
        ax.scatter(bx, bm, s=4, color=SHADE[m], lw=0, zorder=3)
        ax.plot(grid, c, color=SHADE[m], lw=0.9, zorder=4)
        vals += [bm, c[np.isfinite(c)]]
    v = np.concatenate(vals)
    lo, hi = np.floor(v.min()), np.ceil(v.max())
    segments(ax, b12, b23, top=seglabels)
    ax.set_xlim(0, 1)
    ax.set_ylim(lo - 0.05 * (hi - lo), hi + 0.22 * (hi - lo))
    from matplotlib.ticker import MaxNLocator
    ax.yaxis.set_major_locator(MaxNLocator(4, integer=True))
    ax.set_xticks([0, 0.5, 1])
    ax.set_xticklabels(['0', '0.5', '1'])
    if title:
        ax.text(0.03, 0.97 if inside else 1.04, gene, transform=ax.transAxes, fontsize=6.5, style='italic',
                va='top' if inside else 'bottom', color=INK)
    if ylabel:
        ax.set_ylabel('log2(CP10k + 0.5)', labelpad=1.5)
    if xlabel:
        ax.set_xlabel('Pseudostructure', labelpad=1.5)
    else:
        ax.set_xticklabels([])


def panel_maps(fig, root, rect_a, rect_b):
    import anndata as ad

    pt = reference_x(root)
    a = ad.read_h5ad(root / RUNS['hvm'] / '02_qc/units_qc.h5ad', backed='r')
    obs = a.obs
    ref = pd.read_parquet(root / RUNS['hvm'] / '04_pseudospace/pseudospace.parquet')
    gi = list(a.var_names).index('Slc5a2')
    for rect, what in ((rect_a, 'x'), (rect_b, 'expr')):
        x0, y0, w, h = rect
        for i, (s, sp) in enumerate((('Ctrl1A2', 'mouse'), ('HUK1_COR1', 'human'))):
            ax = fig.add_axes([x0 + i * w / 2, y0, w / 2 * 0.96, h])
            units = ref.index[ref['sample'].astype(str) == s]
            xy = obs.loc[units, ['x_centroid', 'y_centroid']].to_numpy(float)
            ax.scatter(xy[:, 0], xy[:, 1], s=0.12, color='#E3E3E3', lw=0, rasterized=True)
            p = pt[pt['sample'] == s]
            pxy = obs.loc[p.index, ['x_centroid', 'y_centroid']].to_numpy(float)
            if what == 'x':
                v, cmap, vmax = p['x'].to_numpy(), 'viridis', 1
            else:
                idx = np.flatnonzero(obs.index.isin(p.index))
                v = np.asarray(a.layers['lognorm'][idx, gi].toarray()).ravel() if hasattr(a.layers['lognorm'][idx, gi], 'toarray') \
                    else np.asarray(a.layers['lognorm'][idx, gi]).ravel()
                v = pd.Series(v, index=obs.index[idx]).loc[p.index].to_numpy()
                cmap, vmax = EXPR_CMAP, np.quantile(v, 0.98)
            order = np.argsort(v) if what == 'expr' else np.arange(len(v))
            sc = ax.scatter(pxy[order, 0], pxy[order, 1], s=0.35, c=v[order], cmap=cmap, vmin=0, vmax=vmax, lw=0,
                            rasterized=True)
            ax.set_aspect('equal')
            ax.invert_yaxis()
            ax.axis('off')
            fig.text(x0 + i * w / 2 + w / 4, y0 + h + 0.004, sp if what == 'x' else sym('Slc5a2', sp), fontsize=6.5,
                     ha='center', va='bottom', color=INK, style='normal' if what == 'x' else 'italic')
        cax = fig.add_axes([x0 + w * 0.3, y0 - 0.012, w * 0.4, 0.006])
        cb = fig.colorbar(sc, cax=cax, orientation='horizontal')
        cb.outline.set_visible(False)
        cb.ax.tick_params(length=1.5, pad=1, labelsize=5)
        if what == 'x':
            cb.set_ticks([0, 1])
            cb.set_ticklabels(['0', '1'])
            cb.set_label('Pseudostructure', fontsize=5.5, labelpad=1)
        else:
            cb.set_ticks([0, vmax])
            cb.set_ticklabels(['0', 'max'])
            cb.set_label('Expression (log)', fontsize=5.5, labelpad=1)


def draw_handoffs(fig, R, top, curves, pos, grid, hp, b12, b23):
    """Transporter pairs in the poster layout: per pair, mouse and human side by side, each species' mean curve (thick)
    over its two sections (thin) in log2(1 + CP10k), one y axis from 0 per species and row, shared by the two partners.
    No switch lines: the switch is defined on min-max scaled curves (Extended Data), which this scale does not show.
    Heights compare the two partners within a species, never mouse with human (different probe panels)."""
    edges = ((0, b12), (b12, b23), (b23, 1))
    cols = {'mouse': (40, 66), 'human': (112, 66)}         # x (mm) and width of the two curve columns
    rh, gap = 9.5, 2.2
    for name, (x0, w) in cols.items():
        fig.text(*R(x0 + w / 2, top), name, ha='center', va='top', fontsize=6.5, weight='bold',
                 color=MOUSE_C if name == 'mouse' else HUMAN_C)
        for (lo, hi), seg in zip(edges, ('S1', 'S2', 'S3')):
            fig.text(*R(x0 + w * (lo + hi) / 2, top + 3.6), seg, ha='center', va='top', fontsize=5.5, color=GREY)
    y0 = top + 6.5
    for r, (e, l, substrate) in enumerate(E_PAIRS):
        yt = y0 + r * (rh + gap)
        fig.text(*R(1.5, yt + 1.0), substrate, fontsize=6, weight='bold', color=INK, va='center')
        t = fig.text(*R(1.5, yt + 4.2), e, fontsize=5.8, style='italic', color=EARLY['mouse'], va='center')
        fig.canvas.draw_idle()
        right = t.get_window_extent(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted()).x1
        fig.text(right + 0.003, R(0, yt + 4.2)[1], f' → {l}', fontsize=5.8, style='italic', color=LATE['mouse'],
                 va='center')
        fig.text(*R(1.5, yt + 7.3), f'{PROTEIN.get(e, e)} → {PROTEIN.get(l, l)}', fontsize=5.2, color=GREY,
                 va='center')
        for name, (x0, w) in cols.items():
            ax = fig.add_axes(R(x0, yt, w, rh, rect=True))
            for (lo, hi), seg in zip(edges, SEGS):
                ax.axvspan(lo, hi, color=BAND[seg], lw=0, zorder=0)
            secs = MOUSE if name == 'mouse' else HUMAN
            hi_y = 0.0
            for gene, col in ((e, EARLY[name]), (l, LATE[name])):
                per = []
                for m in secs:
                    lo_x, hi_x = XRANGE[m]
                    c = np.where((grid >= lo_x) & (grid <= hi_x), curves[m]['curve'][pos[gene]], np.nan)
                    per.append(np.log2(1 + 2.0 ** c))         # stored curves are log2 CP10k
                mean = np.mean(per, 0)          # drawn only where both sections of the species have tubules
                for c in per:
                    ax.plot(grid, np.where(np.isfinite(mean), c, np.nan), color=col, lw=0.4, alpha=0.6, zorder=2)
                ax.plot(grid, mean, color=col, lw=1.1, zorder=3)
                hi_y = max(hi_y, *(np.nanmax(c) for c in per))
            ax.set_ylim(0, hi_y * 1.08)
            ax.set_xlim(0, 1)
            from matplotlib.ticker import MaxNLocator
            ax.yaxis.set_major_locator(MaxNLocator(3, integer=True))
            ax.set_xticks([0, 0.5, 1])
            ax.set_xticklabels(['0', '0.5', '1'] if r == len(E_PAIRS) - 1 else [])
            if r == len(E_PAIRS) - 1:
                ax.set_xlabel('Pseudostructure', labelpad=1)
            if name == 'mouse' and r == len(E_PAIRS) // 2:
                ax.set_ylabel('log$_2$(1 + CP10k)', labelpad=1)
    return y0 + len(E_PAIRS) * (rh + gap)


def panel_shift(fig, rect, genes, hp, b12, b23):
    """Each shifted pair (both partners in the screened set) on the PT axis: an arrow from the mouse to the human switch
    (mean of the two sections) from the COAST curves (ink), and the same from straight lines through each section's
    S1/S2/S3 cluster means (grey, dashed; post hoc), as fig6_panels.cluster_switch_shift locates them."""
    from matplotlib.lines import Line2D
    from pseudospace.handoffs import crossing, scaled_profile

    inside = set(genes.index[genes.universe.astype(bool)])
    d = hp[hp.shifted & hp.early.isin(inside) & hp.late.isin(inside)].sort_values('shift').reset_index(drop=True)
    centres = [b12 / 2, (b12 + b23) / 2, (b23 + 1) / 2]
    g = np.linspace(0, 1, 201)

    def lines(gene, m):
        return scaled_profile(np.interp(g, centres, [genes.loc[gene, f'cluster_S{k}_{m}'] for k in (1, 2, 3)]))

    ax = fig.add_axes(rect)
    for (lo, hi), seg in zip(((0, b12), (b12, b23), (b23, 1)), SEGS):
        ax.axvspan(lo, hi, color=BAND[seg], lw=0, zorder=0)
    for x, s in zip(centres, ('S1', 'S2', 'S3')):
        ax.text(x, 1.01, s, transform=ax.get_xaxis_transform(), ha='center', va='bottom', fontsize=5.5, color=GREY)
    for y, r in d.iterrows():
        xm, xh = (np.mean([r[f'x_{m}'] for m in secs]) for secs in (MOUSE, HUMAN))
        cm, ch = (np.mean([crossing(lines(r.early, m), lines(r.late, m), g) for m in secs]) for secs in (MOUSE, HUMAN))
        for (x0, x1, dy, col, lw, ls) in ((xm, xh, -0.13, INK, 1.0, '-'), (cm, ch, 0.2, '#9A9A9A', 0.8, (0, (2, 1)))):
            if not (np.isfinite(x0) and np.isfinite(x1)):
                ax.text(0.99, y + dy, 'cluster means: no switch', fontsize=4.8, color='#7A7A7A', ha='right',
                        va='center')
                continue
            ax.annotate('', xy=(x1, y + dy), xytext=(x0, y + dy), arrowprops=dict(
                arrowstyle='-|>', color=col, lw=lw, ls=ls, mutation_scale=6, shrinkA=0, shrinkB=0))
            ax.plot([x0, x0], [y + dy - 0.12, y + dy + 0.12], color=col, lw=lw)
        ax.text(-0.015, y, f'{PROTEIN.get(r.early, r.early)} → {PROTEIN.get(r.late, r.late)}', fontsize=5.5,
                ha='right', va='center', transform=ax.get_yaxis_transform())
    ax.set_ylim(len(d) - 0.45, -0.6)
    ax.set_xlim(0, 1)
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.set_xticks([0, 0.5, 1])
    ax.set_xticklabels(['0', '0.5', '1'])
    ax.set_xlabel('Pseudostructure', labelpad=1)
    ax.legend(handles=[Line2D([], [], color=INK, lw=1.0, marker='>', ms=3, label='COAST curves'),
                       Line2D([], [], color='#9A9A9A', lw=0.8, ls=(0, (2, 1)), marker='>', ms=3,
                              label='Cluster means')],
              title='Switch point, mouse → human', title_fontsize=5.5, frameon=False, fontsize=5.5,
              loc='center left', bbox_to_anchor=(1.02, 0.5), handlelength=2.0, alignment='left')
    return d


def main(argv=None):
    args = parse_args(argv)
    rroot = Path(args.results_root)
    root = rroot / 'paper_figures/fig5'
    screen = root / 'screen'
    plt = style()
    global EXPR_CMAP
    from matplotlib.colors import LinearSegmentedColormap
    EXPR_CMAP = LinearSegmentedColormap.from_list('expr', ['#F3E9F0', '#C51B7D', '#3D0425'])
    genes, grid, curves, pos = load(screen)
    hp = pd.read_csv(screen / 'handoffs.csv')
    b12, b23 = float(genes.b12.median()), float(genes.b23.median())
    XRANGE.update(handoff_range(screen, hp, genes))

    W, H = 180, 204                                       # mm
    fig = plt.figure(figsize=(W * MM, H * MM))

    def R(x, y, w=0, h=0, rect=False):
        """mm from the top-left corner to figure fractions; rect=True gives [x0, y0, w, h] for add_axes."""
        if rect:
            return [x / W, 1 - (y + h) / H, w / W, h / H]
        return x / W, 1 - y / H

    boxes = {}

    def panel(letter, x0, y0, x1, y1):
        fig.text(*R(x0, y0), letter, fontsize=8, weight='bold', ha='left', va='top', color=INK)
        boxes[letter] = (x0, y0, x1, y1)

    # row 1: tissue
    panel('a', 1, 1, 89, 50)
    panel('b', 91, 1, 179, 50)
    panel_maps(fig, rroot, R(8, 7, 78, 34, rect=True), R(98, 7, 78, 34, rect=True))

    # row 2: curves
    panel('c', 1, 52, 62, 104)
    ax = fig.add_axes(R(12, 60, 46, 34, rect=True))
    curve_axis(ax, 'Slc5a2', curves, pos, grid, b12, b23, ylabel=True, xlabel=True, inside=True)
    panel('d', 64, 52, 179, 104)
    gw, gh = 30, 15
    for k, g in enumerate(GENES_D):
        r, c = divmod(k, 3)
        ax = fig.add_axes(R(75 + c * (gw + 6.5), 59 + r * (gh + 6), gw, gh, rect=True))
        curve_axis(ax, g, curves, pos, grid, b12, b23, ylabel=c == 0, xlabel=r == 1 or c == 2, seglabels=False)
    key = fig.add_axes(R(75 + 2 * (gw + 6.5), 59 + gh + 9, gw, gh, rect=True))
    key.axis('off')
    for m, nm in zip(SECTIONS, ('Mouse section 1', 'Mouse section 2', 'Human section 1', 'Human section 2')):
        key.plot([], [], color=SHADE[m], lw=0.9, marker='o', ms=2, label=nm)
    key.plot([], [], color=GREY, lw=0.5, ls=(0, (2, 2)), label='S1/S2/S3 boundaries')
    key.legend(frameon=False, fontsize=5.3, loc='center left', handlelength=1.6, handletextpad=0.4, labelspacing=0.4)

    # row 3: transporter switches (earlier handoff-rows style)
    panel('e', 1, 106, 179, 154)
    draw_handoffs(fig, R, 107, curves, pos, grid, hp, b12, b23)

    # row 4: where each shifted switch sits in mouse and human, and what clustering recovers
    panel('f', 1, 157, 179, 203)
    panel_shift(fig, R(40, 164, 95, 32, rect=True), genes, hp, b12, b23)

    fig.savefig(root / 'fig5.png', facecolor='white', dpi=400)
    fig.savefig(root / 'fig5.pdf', facecolor='white')
    plt.close(fig)

    from matplotlib.patches import Rectangle
    wf = plt.figure(figsize=(W * MM, H * MM))
    for k, (x0, y0, x1, y1) in boxes.items():
        wf.patches.append(Rectangle((x0 / W, 1 - y1 / H), (x1 - x0) / W, (y1 - y0) / H, transform=wf.transFigure,
                                    fc='#E6E6E6', ec='#9A9A9A', lw=0.6))
        wf.text(x0 / W + 0.01, 1 - y0 / H - 0.01, k, fontsize=9, weight='bold', va='top')
    wf.savefig(root / 'fig5_wireframe.png', facecolor='white', dpi=200)
    plt.close(wf)


if __name__ == '__main__':
    main()
