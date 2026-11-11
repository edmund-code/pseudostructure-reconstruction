"""Fig. 5 drawn as one journal-style canvas (180 mm wide), restyled 2026-11-11 after reading published continuous-axis
methods figures (GASTON, Genes2Genes, STcompare, STalign, Halpern 2017 liver zonation):
 a  one mouse and one human section, PT tubules coloured by pseudostructure
 b  Slc5a2 (SGLT2) expression in the same tubules
 c  Slc5a2 against pseudostructure: binned means (points) and fitted curves (lines) in log2(CP10k + 0.5)
 d  the same for five more transporters (Slc2a2, GLUT2, is the unmoved comparator); the sixth slot is the section key
 e  two transporter switches, mouse and human, with each species' switch point
 f  each shifted switch's species shift from the curves against straight lines through the S1/S2/S3 cluster means
 g  share of PT genes strongly zonated, by orthology class
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
from fig5_panels import (HUMAN, HUMAN_C, MOUSE, MOUSE_C, ORTHO_CLASSES, PT, SECTIONS, XRANGE, handoff_range,  # noqa: E402
                         load, reference_x)
from fig6_panels import cluster_switch_shift  # noqa: E402

MM = 1 / 25.4
INK, GREY = '#222222', '#8A8A8A'
TRIM = 0.05                      # spline ends overshoot (question-bank.md); curves and bins shown on [0.05, 0.95]
# Abcc2 (MRP2) is a strong call but its two human curves disagree in shape; it is shown in Extended Data
GENES_D = ['Slc2a2', 'Slc22a6', 'Slc22a8', 'Slc22a2', 'Slc22a12']
SWITCHES = [('Slc5a2', 'Slc5a10'), ('Slc22a23', 'Slc22a12')]
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


def label(fig, x, y, s):
    fig.text(x, y, s, fontsize=8, weight='bold', ha='left', va='top', color=INK)


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


def panel_switch(fig, rect, curves, pos, grid, hp, b12, b23):
    from pseudospace.handoffs import scaled_profile

    x0, y0, w, h = rect
    cols = {'mouse': (MOUSE, ('#7FB3E6', '#173A63')), 'human': (HUMAN, ('#F0A3CB', '#7A0F4D'))}
    gw, gh = w / 2 - 0.012, h / 2 - 0.03
    for r, (e, l) in enumerate(SWITCHES):
        row = hp[(hp.early == e) & (hp.late == l)].iloc[0]
        for c, (sp, (secs, (ce, cl))) in enumerate(cols.items()):
            ax = fig.add_axes([x0 + c * (gw + 0.024), y0 + (1 - r) * (gh + 0.06), gw, gh])
            for gene, col in ((e, ce), (l, cl)):
                for m in secs:
                    lo_x, hi_x = XRANGE[m]
                    raw = np.where((grid >= lo_x) & (grid <= hi_x), curves[m]['curve'][pos[gene]], np.nan)
                    prof = scaled_profile(raw)
                    lo, hi = np.nanmin(raw), np.nanmax(raw)
                    bx, bm = curves[m]['bin_x'][pos[gene]], curves[m]['bins'][pos[gene]].astype(float)
                    k = np.isfinite(bm) & (bx >= lo_x) & (bx <= hi_x)
                    ax.scatter(bx[k], (bm[k] - lo) / (hi - lo), s=3, color=col, lw=0, zorder=3)
                    ax.plot(grid, prof, color=col, lw=0.8, zorder=4)
            xs = np.nanmean([row[f'x_{m}'] for m in secs])
            ax.axvline(xs, color=MOUSE_C if sp == 'mouse' else HUMAN_C, lw=0.8, zorder=5)
            if sp == 'human':
                ax.axvline(np.nanmean([row[f'x_{m}'] for m in MOUSE]), color=MOUSE_C, lw=0.8, ls=(0, (2, 1.5)), zorder=5)
            segments(ax, b12, b23, top=r == 0)
            ax.set_xlim(0, 1)
            ax.set_ylim(-0.15, 1.2)
            ax.set_xticks([0, 0.5, 1])
            ax.set_xticklabels(['0', '0.5', '1'])
            ax.set_yticks([0, 1])
            if c:
                ax.set_yticklabels([])
            else:
                ax.set_ylabel('Scaled expr.', labelpad=1)
            if r == 0:
                ax.set_xticklabels([])
            else:
                ax.set_xlabel('Pseudostructure', labelpad=1.5)
            t = ax.text(0.02, 1.13, sym(e, sp), transform=ax.transAxes, fontsize=6, style='italic', color=ce,
                        va='bottom')
            ax.annotate(f' → {sym(l, sp)}', xycoords=t, xy=(1, 0), va='bottom', fontsize=6, style='italic', color=cl)


def panel_shift(fig, rect, genes, hp, b12, b23):
    d = cluster_switch_shift(genes, hp, b12, b23, MOUSE, HUMAN)
    inside = set(genes.index[genes.universe.astype(bool)])
    d = d[d.early.isin(inside) & d.late.isin(inside)].reset_index(drop=True)
    ax = fig.add_axes(rect)
    for y, r in d.iterrows():
        if np.isfinite(r.clusters):
            ax.plot([r.clusters, r.curve], [y, y], color='#C8C8C8', lw=0.8, zorder=1)
            ax.scatter([r.clusters], [y], s=9, facecolor='white', edgecolor=INK, lw=0.6, zorder=3)
        else:
            ax.text(0, y, '×', fontsize=6, ha='center', va='center', color=INK, zorder=3)
        ax.scatter([r.curve], [y], s=9, color=INK, lw=0, zorder=4)
        ax.text(-0.03, y, f'{r.early} → {r.late}', fontsize=5, style='italic', ha='right', va='center',
                transform=ax.get_yaxis_transform())
    ax.axvline(0, color=GREY, lw=0.5, zorder=0)
    ax.set_ylim(-0.8, len(d) - 0.2)
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.set_xlabel('Switch shift, human − mouse', labelpad=1.5)
    ax.scatter([], [], s=9, color=INK, label='COAST curve')
    ax.scatter([], [], s=9, facecolor='white', edgecolor=INK, lw=0.6, label='Cluster means')
    ax.legend(frameon=False, fontsize=5, loc='lower left', bbox_to_anchor=(-0.02, 1.0), ncol=2, handletextpad=0.2,
              borderaxespad=0.1, columnspacing=0.8)
    ax.set_xticks([-0.5, -0.25, 0, 0.25])
    ax.set_xticklabels(['−0.5', '', '0', '0.25'])
    return d


def panel_orthology(fig, rect, rroot):
    from statsmodels.stats.proportion import proportion_confint

    t = pd.read_csv(rroot / 'pt_non_one_to_one_package/tables/h1_by_class.csv').set_index(['species', 'class_support3'])
    ax = fig.add_axes(rect)
    for i, sp in enumerate(('mouse', 'human')):
        col = MOUSE_C if sp == 'mouse' else HUMAN_C
        for j, c in enumerate(ORTHO_CLASSES):
            n, s = t.loc[(sp, c), 'expressed'], t.loc[(sp, c), 'strong']
            lo, hi = proportion_confint(s, n, method='wilson')
            x = j + (i - 0.5) * 0.38
            ax.bar(x, s / n * 100, width=0.36, color=col, lw=0, label=sp if j == 0 else None)
            ax.plot([x, x], [lo * 100, hi * 100], color=INK, lw=0.5)
    ax.set_xticks(range(len(ORTHO_CLASSES)))
    ax.set_xticklabels(ORTHO_CLASSES)
    ax.set_xlabel('Orthology (human : mouse copies)', labelpad=1.5)
    ax.set_ylabel('Strongly zonated (%)', labelpad=1.5)
    ax.set_yticks([0, 5, 10, 15, 20])
    ax.legend(frameon=False, fontsize=5.5, loc='upper left', handlelength=0.8, handletextpad=0.3)


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

    W, H = 180 * MM, 196 * MM
    fig = plt.figure(figsize=(W, H))
    boxes = {}

    # row 1: tissue
    ra, rb = [0.045, 0.775, 0.43, 0.19], [0.545, 0.775, 0.43, 0.19]
    panel_maps(fig, rroot, ra, rb)
    label(fig, 0.01, 0.99, 'a')
    label(fig, 0.51, 0.99, 'b')
    boxes['a'], boxes['b'] = (0.01, 0.745, 0.49, 0.99), (0.51, 0.745, 0.99, 0.99)

    # row 2: curves
    ax = fig.add_axes([0.065, 0.455, 0.25, 0.215])
    curve_axis(ax, 'Slc5a2', curves, pos, grid, b12, b23, ylabel=True, xlabel=True, inside=True)
    label(fig, 0.01, 0.715, 'c')
    boxes['c'] = (0.01, 0.42, 0.34, 0.715)
    gx0, gw, gh = 0.405, 0.17, 0.095
    for k, g in enumerate(GENES_D):
        r, c = divmod(k, 3)
        ax = fig.add_axes([gx0 + c * (gw + 0.025), 0.565 - r * (gh + 0.025), gw, gh])
        curve_axis(ax, g, curves, pos, grid, b12, b23, ylabel=c == 0, xlabel=r == 1 or c == 2, seglabels=False)
    key = fig.add_axes([gx0 + 2 * (gw + 0.025), 0.565 - (gh + 0.025), gw, gh])
    key.axis('off')
    for m, nm in zip(SECTIONS, ('Mouse section 1', 'Mouse section 2', 'Human section 1', 'Human section 2')):
        key.plot([], [], color=SHADE[m], lw=0.9, marker='o', ms=2, label=nm)
    key.plot([], [], color=GREY, lw=0.5, ls=(0, (2, 2)), label='S1/S2/S3 boundaries')
    key.legend(frameon=False, fontsize=5.5, loc='center left', handlelength=1.6, handletextpad=0.4, labelspacing=0.5)
    label(fig, 0.36, 0.715, 'd')
    boxes['d'] = (0.36, 0.42, 0.99, 0.715)

    # row 3: switches, shift, orthology
    panel_switch(fig, [0.065, 0.075, 0.30, 0.29], curves, pos, grid, hp, b12, b23)
    label(fig, 0.01, 0.395, 'e')
    boxes['e'] = (0.01, 0.01, 0.39, 0.395)
    panel_shift(fig, [0.53, 0.075, 0.16, 0.29], genes, hp, b12, b23)
    label(fig, 0.40, 0.395, 'f')
    boxes['f'] = (0.40, 0.01, 0.70, 0.395)
    panel_orthology(fig, [0.765, 0.075, 0.225, 0.27], rroot)
    label(fig, 0.71, 0.395, 'g')
    boxes['g'] = (0.71, 0.01, 0.99, 0.395)

    fig.savefig(root / 'fig5.png', facecolor='white', dpi=400)
    fig.savefig(root / 'fig5.pdf', facecolor='white')
    plt.close(fig)

    from matplotlib.patches import Rectangle
    wf = plt.figure(figsize=(W, H))
    for k, (x0, y0, x1, y1) in boxes.items():
        wf.patches.append(Rectangle((x0, y0), x1 - x0, y1 - y0, transform=wf.transFigure, fc='#E6E6E6', ec='#9A9A9A',
                                    lw=0.6))
        wf.text(x0 + 0.01, y1 - 0.01, k, fontsize=9, weight='bold', va='top')
    wf.savefig(root / 'fig5_wireframe.png', facecolor='white', dpi=200)
    plt.close(wf)


if __name__ == '__main__':
    main()
