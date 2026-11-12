"""Fig. 5 drawn as one journal-style canvas (180 mm wide), restyled 2026-11-11 after reading published continuous-axis
methods figures (GASTON, Genes2Genes, STcompare, STalign, Halpern 2017 liver zonation):
 a  one mouse and one human section, PT tubules coloured by pseudostructure
 b  Slc5a2 (SGLT2) expression in the same tubules
 c  Slc5a2 against pseudostructure: binned means (points) and fitted curves (lines) in log2(CP10k + 0.5)
 d  the same for five more transporters (Slc2a2, GLUT2, is the unmoved comparator); the sixth slot is the section key
 e  transporter switches in the earlier handoff-rows style: mean log2(CP10k + 0.5) curves per species and switch lines
 f  each shifted switch's species shift, from the curves and from straight lines through the S1/S2/S3 cluster means
 g  share of PT genes strongly zonated in each orthology class relative to the same species' one-to-one genes
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
from fig5_panels import (EARLY, H_PAIRS, HUMAN, HUMAN_C, LATE, MOUSE, MOUSE_C, ORTHO_CLASSES, PROTEIN, SECTIONS,  # noqa: E402
                         XRANGE, handoff_range, load, reference_x)
from fig6_panels import BAND, SEGS  # noqa: E402
from fig6_panels import cluster_switch_shift  # noqa: E402

MM = 1 / 25.4
INK, GREY = '#222222', '#8A8A8A'
TRIM = 0.05                      # spline ends overshoot (question-bank.md); curves and bins shown on [0.05, 0.95]
# Abcc2 (MRP2) is a strong call but its two human curves disagree in shape; it is shown in Extended Data
GENES_D = ['Slc2a2', 'Slc22a6', 'Slc22a8', 'Slc22a2', 'Slc22a12']
EXPR_CMAP = None
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
    """The earlier handoff-rows panel (fig5_panels.panel_handoff_rows) on the journal canvas: per pair, mouse and
    human side by side, each species' mean scaled curve (thick, shaded) over its two sections (thin), the mouse switch
    (dashed blue) and human switch (dotted magenta). Expression is log2(CP10k + 0.5), as in c and d, with a y axis per species and
    row; the switch lines come from the normalized crossing, so they need not fall where the drawn curves cross."""
    edges = ((0, b12), (b12, b23), (b23, 1))
    cols = {'mouse': (40, 66), 'human': (112, 66)}         # x (mm) and width of the two curve columns
    rh, gap = 9.5, 2.2
    for name, (x0, w) in cols.items():
        fig.text(*R(x0 + w / 2, top), name, ha='center', va='top', fontsize=6.5, weight='bold',
                 color=MOUSE_C if name == 'mouse' else HUMAN_C)
        for (lo, hi), seg in zip(edges, ('S1', 'S2', 'S3')):
            fig.text(*R(x0 + w * (lo + hi) / 2, top + 3.6), seg, ha='center', va='top', fontsize=5.5, color=GREY)
    for k, (txt, ls, col) in enumerate((('Mouse switch', (0, (3, 1.5)), MOUSE_C), ('Human switch', (0, (1, 1.2)), HUMAN_C))):
        y = top + 0.6 + k * 2.6
        fig.lines.append(plt_line(fig, R(6, y), R(11, y), col, ls))
        fig.text(*R(12, y), txt, fontsize=5.5, va='center', color=INK)
    y0 = top + 6.5
    for r, (e, l, substrate) in enumerate(H_PAIRS):
        yt = y0 + r * (rh + gap)
        row = hp[(hp.early == e) & (hp.late == l)].iloc[0]
        xs = {sp: np.nanmean([row[f'x_{m}'] for m in secs]) for sp, secs in (('mouse', MOUSE), ('human', HUMAN))}
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
            drawn = []
            for gene, col in ((e, EARLY[name]), (l, LATE[name])):
                per = []
                for m in secs:
                    lo_x, hi_x = XRANGE[m]
                    c = np.where((grid >= lo_x) & (grid <= hi_x), curves[m]['curve'][pos[gene]], np.nan)
                    per.append(np.log2(2.0 ** c + 0.5))       # log2(CP10k + 0.5), as in c and d
                mean = np.mean(per, 0)          # drawn only where both sections of the species have tubules
                for c in per:
                    ax.plot(grid, np.where(np.isfinite(mean), c, np.nan), color=col, lw=0.4, alpha=0.6, zorder=2)
                ax.plot(grid, mean, color=col, lw=1.1, zorder=3)
                drawn += [mean, *per]
            ax.axvline(xs['mouse'], color=MOUSE_C, lw=0.8, ls=(0, (3, 1.5)), zorder=4)
            if name == 'human':
                ax.axvline(xs['human'], color=HUMAN_C, lw=0.9, ls=(0, (1, 1.2)), zorder=4)
            v = np.concatenate([d[np.isfinite(d)] for d in drawn])
            lo_y, hi_y = v.min(), v.max()
            pad = 0.08 * max(hi_y - lo_y, 0.5)
            ax.set_ylim(lo_y - pad, hi_y + pad)
            ax.set_xlim(0, 1)
            from matplotlib.ticker import MaxNLocator
            ax.yaxis.set_major_locator(MaxNLocator(3, integer=True))
            ax.set_xticks([0, 0.5, 1])
            ax.set_xticklabels(['0', '0.5', '1'] if r == len(H_PAIRS) - 1 else [])
            if r == len(H_PAIRS) - 1:
                ax.set_xlabel('Pseudostructure', labelpad=1)
            if name == 'mouse' and r == len(H_PAIRS) // 2:
                ax.set_ylabel('log$_2$(CP10k + 0.5)', labelpad=1)
    return y0 + len(H_PAIRS) * (rh + gap)


def plt_line(fig, p0, p1, col, ls):
    from matplotlib.lines import Line2D
    return Line2D([p0[0], p1[0]], [p0[1], p1[1]], transform=fig.transFigure, color=col, lw=0.9, ls=ls)


def panel_shift(fig, rect, genes, hp, b12, b23):
    d = cluster_switch_shift(genes, hp, b12, b23, MOUSE, HUMAN)
    inside = set(genes.index[genes.universe.astype(bool)])
    d = d[d.early.isin(inside) & d.late.isin(inside)].reset_index(drop=True)
    ax = fig.add_axes(rect)
    bh = 0.34
    for y, r in d.iterrows():
        ax.barh(y - bh / 2 - 0.02, r.curve, height=bh, color=INK, lw=0, zorder=2)
        if np.isfinite(r.clusters):
            ax.barh(y + bh / 2 + 0.02, r.clusters, height=bh, facecolor='white', edgecolor=INK, lw=0.5, zorder=2)
        else:
            ax.text(0.008, y + bh / 2 + 0.02, 'no switch', fontsize=4.6, ha='left', va='center', color=GREY, zorder=3)
        ax.text(-0.03, y, f'{r.early} → {r.late}', fontsize=5, style='italic', ha='right', va='center',
                transform=ax.get_yaxis_transform())
    ax.axvline(0, color=INK, lw=0.5, zorder=3)
    ax.set_ylim(len(d) - 0.4, -0.6)
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.set_xlabel('Switch shift, human − mouse', labelpad=1.5)
    ax.bar(0, 0, color=INK, lw=0, label='COAST curve')
    ax.bar(0, 0, facecolor='white', edgecolor=INK, lw=0.5, label='Cluster means')
    ax.legend(frameon=False, fontsize=5, loc='lower left', bbox_to_anchor=(-0.02, 1.0), ncol=2, handletextpad=0.3,
              borderaxespad=0.1, columnspacing=0.8, handlelength=1.1)
    ax.set_xticks([-0.5, -0.25, 0, 0.25])
    ax.set_xticklabels(['−0.5', '', '0', '0.25'])
    return d


def panel_orthology(fig, rect, rroot):
    """Each orthology class's share of strongly zonated PT genes relative to the one-to-one share of the same species
    (so the flatter human gradients do not enter the comparison), with Katz log 95% intervals (0.5 added to empty
    cells)."""
    t = pd.read_csv(rroot / 'pt_non_one_to_one_package/tables/h1_by_class.csv').set_index(['species', 'class_support3'])
    ax = fig.add_axes(rect)
    classes = [c for c in ORTHO_CLASSES if c != '1:1']
    for i, sp in enumerate(('mouse', 'human')):
        col = MOUSE_C if sp == 'mouse' else HUMAN_C
        n0, s0 = t.loc[(sp, '1:1'), 'expressed'], t.loc[(sp, '1:1'), 'strong']
        for j, c in enumerate(classes):
            n, s = t.loc[(sp, c), 'expressed'], t.loc[(sp, c), 'strong']
            ss, ss0 = (s + 0.5, s0 + 0.5) if s == 0 else (s, s0)
            rr = (s / n) / (s0 / n0)
            se = np.sqrt(1 / ss - 1 / n + 1 / ss0 - 1 / n0)
            lo, hi = (ss / n) / (ss0 / n0) * np.exp(-1.96 * se), (ss / n) / (ss0 / n0) * np.exp(1.96 * se)
            x = j + (i - 0.5) * 0.36
            if s == 0:      # no strongly zonated gene: open bar at the 0.5-corrected estimate
                rr = (ss / n) / (ss0 / n0)
                ax.bar(x, rr - 1, width=0.3, bottom=1, facecolor='white', edgecolor=col, lw=0.6, zorder=2)
                hi = min(hi, 12)
            else:           # bars grow from the 1:1 reference, so the dashed line is the baseline
                ax.bar(x, rr - 1, width=0.3, bottom=1, color=col, lw=0, zorder=2,
                       label=sp if j == 0 else None)
            ax.plot([x, x], [lo, hi], color=INK, lw=0.6, zorder=4)
            ax.plot([x - 0.07, x + 0.07], [lo, lo], color=INK, lw=0.6, zorder=4)
            ax.plot([x - 0.07, x + 0.07], [hi, hi], color=INK, lw=0.6, zorder=4)
    ax.axhline(1, color=GREY, lw=0.5, ls=(0, (2, 2)), zorder=0)
    ax.set_yscale('log')
    ax.set_ylim(0.1, 12)
    ax.set_yticks([0.1, 0.3, 1, 3, 10])
    ax.set_yticklabels(['0.1', '0.3', '1', '3', '10'])
    ax.set_xticks(range(len(classes)))
    ax.set_xticklabels(classes)
    ax.set_xlim(-0.6, len(classes) - 0.4)
    ax.set_xlabel('Orthology (human : mouse copies)', labelpad=1.5)
    ax.set_ylabel('Strongly zonated, relative to 1:1', labelpad=1.5)
    ax.legend(frameon=False, fontsize=5.5, loc='upper left', handletextpad=0.2)


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

    W, H = 180, 228                                       # mm
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
    panel('e', 1, 106, 179, 172)
    draw_handoffs(fig, R, 107, curves, pos, grid, hp, b12, b23)

    # row 4: what clustering recovers, orthology
    panel('f', 1, 175, 99, 227)
    panel_shift(fig, R(33, 183, 50, 36, rect=True), genes, hp, b12, b23)
    panel('g', 101, 175, 179, 227)
    panel_orthology(fig, R(115, 181, 62, 36, rect=True), rroot)

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
