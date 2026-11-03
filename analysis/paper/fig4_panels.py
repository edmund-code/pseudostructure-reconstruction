"""Fig. 4 panels drawn from the screen (fig4_screen.py), sized for the figure layout (make_fig4.py).

 a  one gene three ways: bulk, S1/S2/S3 clusters, pseudostructure (the top category 1 gene)
 b  every gene that passes either category: curves next to the same genes painted with their cluster means, and bulk
 c  the top three category 1 genes (peak between segments)
 d  the top three category 2 genes (change within a segment)

Reads <results-root>/paper_figures/fig4/screen/ and writes PNGs to <results-root>/paper_figures/fig4/panels/.
Protocol: docs/paper/fig4-continuous-vs-discrete.md.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = next(p for p in Path(__file__).resolve().parents if (p / 'pseudospace').is_dir())
sys.path.insert(0, str(Path(__file__).resolve().parent))
from fig3_panels import COLOURS, INK, MUTED, SECTION_COLOURS, style  # noqa: E402

MICE = ['Ctrl1A2', 'Ctrl1A4']
SEGS = ['PTS1', 'PTS2', 'PTS3']
CLUSTER, BULK = '#555555', '#AAAAAA'
LABEL_LOG2 = 0.5
# panel boxes of make_fig4.py in SVG units (1 unit = 0.01 inch here)
SIZE = {'a': (540, 292), 'b': (580, 967), 'c': (540, 282), 'd': (540, 282)}


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--results-root', default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT', str(REPO / 'results')))
    p.add_argument('--out', default='paper_figures/fig4', help='relative to the results root')
    return p.parse_args(argv)


def load(root: Path):
    genes = pd.read_csv(root / 'screen/genes.csv').set_index('gene')
    z = np.load(root / 'screen/curves.npz')
    names = list(z['genes'])
    curves = {m: {part: dict(zip(names, z[f'{part}|{m}'])) for part in ('curve', 'clusters', 'bulk', 'bins', 'bin_x')}
              for m in MICE}
    # Post-hoc rule (after the screen, before the figure was final; recorded in the protocol doc): a call is shown
    # only if it holds on raw counts within single labels, by >= LABEL_LOG2 in the same direction in both mice.
    # Category 2 examples are ranked by that within-label change: the curve's change can be inflated where a
    # spline dives into near-zero counts (Dpys in late S3).
    chk = pd.read_csv(root / 'screen/label_check.csv').set_index('gene')
    a, b = MICE
    within = chk[[f'within_label_change_{a}', f'within_label_change_{b}']]
    same = (np.sign(within[f'within_label_change_{a}']) == np.sign(chk[f'curve_change_{a}'])) & \
           (np.sign(within[f'within_label_change_{b}']) == np.sign(chk[f'curve_change_{b}']))
    ok2 = same & (within.abs().min(1) >= LABEL_LOG2)
    peak = chk[[f'rise_in_S2_{a}', f'rise_in_S2_{b}', f'fall_in_S3_{a}', f'fall_in_S3_{b}']].min(1)
    genes['label_ok'] = False
    genes.loc[chk.index, 'label_ok'] = np.where(chk.category == 1, peak >= LABEL_LOG2, ok2)
    genes['label_score'] = np.nan
    genes.loc[chk.index, 'label_score'] = np.where(chk.category == 1, peak, within.abs().min(1))
    genes['cat1'] &= genes.label_ok
    genes['cat2'] &= genes.label_ok
    return genes, z['grid'], curves


def segment_band(ax, b12, b23, y, h):
    from matplotlib.patches import Rectangle

    for (lo, hi), seg in zip(((0, b12), (b12, b23), (b23, 1)), SEGS):
        ax.add_patch(Rectangle((lo, y), hi - lo, h, color=COLOURS[seg], lw=0, clip_on=False))


def gene_axes(ax, gene, row, grid, curves, show=('bins', 'bulk', 'clusters', 'curve'), ylim=None):
    """One gene in both mice: bin means, bulk, cluster steps and the curve."""
    b12, b23 = row.b12, row.b23
    edges = ((0, b12), (b12, b23), (b23, 1))
    for m in MICE:
        c = SECTION_COLOURS[m]
        if 'bins' in show:
            ax.scatter(curves[m]['bin_x'][gene], curves[m]['bins'][gene], s=5, color=c, lw=0, alpha=0.8, zorder=3)
        if 'bulk' in show:
            ax.plot([0, 1], [curves[m]['bulk'][gene]] * 2, color=BULK, lw=1.0, ls=(0, (2, 1.5)), zorder=1)
        if 'clusters' in show:
            for (lo, hi), v in zip(edges, curves[m]['clusters'][gene]):
                ax.plot([lo, hi], [v, v], color=CLUSTER, lw=1.3, zorder=2, solid_capstyle='butt')
        if 'curve' in show:
            ax.plot(grid, curves[m]['curve'][gene], color=c, lw=1.4, zorder=4)
    lo, hi = ylim or ax.get_ylim()
    ax.set_ylim(lo, hi)
    pad = (hi - lo) * 0.06
    segment_band(ax, b12, b23, lo - 2.2 * pad, pad)
    for b in (b12, b23):
        ax.axvline(b, color='#CCCCCC', lw=0.5, zorder=0)
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['', ''])
    ax.tick_params(axis='x', length=0)


def ylim_of(gene, curves, pad=0.5):
    """The data's range (bin means and cluster means); a spline that dives where counts run out is clipped."""
    v = np.concatenate([np.r_[curves[m]['bins'][gene], curves[m]['clusters'][gene]] for m in MICE])
    v = v[np.isfinite(v)]
    return v.min() - pad, v.max() + pad


def examples(genes):
    c1 = genes[genes.cat1].sort_values('cat1_score', ascending=False).index.tolist()
    c2 = genes[genes.cat2 & ~genes.cat1].sort_values('label_score', ascending=False).index.tolist()
    return c1, c2


# ---------------------------------------------------------------------------------------------------- a
def panel_a(plt, genes, grid, curves, out):
    c1, c2 = examples(genes)
    gene = (c1 or c2)[0]
    fig = plt.figure(figsize=(SIZE['a'][0] / 100, SIZE['a'][1] / 100))
    lim = ylim_of(gene, curves)
    titles = ('bulk', 'S1 / S2 / S3 clusters', 'pseudostructure')
    for i, show in enumerate((('bins', 'bulk'), ('bins', 'clusters'), ('bins', 'curve'))):
        ax = fig.add_axes([0.08 + i * 0.31, 0.16, 0.27, 0.70])
        gene_axes(ax, gene, genes.loc[gene], grid, curves, show, lim)
        ax.set_title(titles[i], fontsize=7.5, color=INK, pad=3)
        if i:
            ax.set_yticklabels([])
        else:
            ax.set_ylabel('log2 counts per 10k')
            ax.text(0.02, 0.97, gene, transform=ax.transAxes, fontsize=8, style='italic', color=INK, va='top')
        if i == 2:
            for j, m in enumerate(MICE):
                ax.text(0.03, 0.96 - j * 0.08, m, transform=ax.transAxes, fontsize=6.5, color=SECTION_COLOURS[m],
                        va='top')
    fig.savefig(out / 'f4_a.png', facecolor='white')
    plt.close(fig)
    return gene


# ---------------------------------------------------------------------------------------------------- b
def panel_b(plt, genes, grid, curves, out, mark):
    from matplotlib.colors import LinearSegmentedColormap

    c1 = genes[genes.cat1].copy()
    c2 = genes[genes.cat2 & ~genes.cat1].copy()
    c1['key'] = c1['peak_at_Ctrl1A2'] * 10 + c1[['peak_x_Ctrl1A2', 'peak_x_Ctrl1A4']].mean(1)
    sign = np.sign(c2['change_Ctrl1A2'])
    c2['key'] = c2['change_segment_Ctrl1A2'] * 10 + (sign < 0) * 5 + c2['cat2_score'].rank(ascending=False) / 1e4
    rows = c1.sort_values('key').index.tolist() + c2.sort_values('key').index.tolist()
    n1 = len(c1)
    curve = np.stack([np.mean([curves[m]['curve'][g] for m in MICE], 0) for g in rows])
    steps = np.zeros_like(curve)
    bulk = np.zeros(len(rows))
    for i, g in enumerate(rows):
        b12, b23 = genes.loc[g, ['b12', 'b23']]
        seg = np.digitize(grid, [b12, b23])
        steps[i] = np.mean([curves[m]['clusters'][g] for m in MICE], 0)[seg]
        bulk[i] = np.mean([curves[m]['bulk'][g] for m in MICE])
    lo, hi = curve.min(1, keepdims=True), curve.max(1, keepdims=True)
    scale = lambda v: (v - lo) / (hi - lo)  # noqa: E731
    cmap = LinearSegmentedColormap.from_list('seq', ['#FFFFFF', '#C9D3EE', '#5B7DC8', '#1D2F6F'])
    fig = plt.figure(figsize=(SIZE['b'][0] / 100, SIZE['b'][1] / 100))
    top, bottom = 0.955, 0.04
    widths = [(0.12, 0.38), (0.54, 0.38), (0.95, 0.035)]
    titles = ('pseudostructure', 'S1 / S2 / S3 clusters', 'bulk')
    for (x0, w), mat, title in zip(widths, (scale(curve), scale(steps), scale(bulk[:, None])), titles):
        ax = fig.add_axes([x0, bottom, w, top - bottom])
        ax.set_title(title, fontsize=7.5, color=INK, pad=3)
        ax.imshow(np.clip(mat, 0, 1), aspect='auto', interpolation='nearest', cmap=cmap, vmin=0, vmax=1,
                  extent=(0, 1, len(rows), 0))
        ax.axhline(n1, color=INK, lw=0.8)
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        if w > 0.1:
            band = fig.add_axes([x0, bottom - 0.014, w, 0.009])
            b12, b23 = genes.loc[rows, 'b12'].mean(), genes.loc[rows, 'b23'].mean()
            segment_band(band, b12, b23, 0, 1)
            band.set_xlim(0, 1)
            band.set_ylim(0, 1)
            band.axis('off')
    lab = fig.add_axes([0.0, bottom, 0.115, top - bottom])
    lab.set_xlim(0, 1)
    lab.set_ylim(len(rows), 0)
    lab.axis('off')
    for g in mark:
        if g in rows:
            lab.text(0.95, rows.index(g) + 0.5, g, ha='right', va='center', fontsize=6, style='italic', color=INK)
    fig.savefig(out / 'f4_b.png', facecolor='white')
    plt.close(fig)
    return {'n_rows': len(rows), 'n_cat1': n1}


# ---------------------------------------------------------------------------------------------------- b (lines)
GROUPS = [('peak at S2/S3', lambda g: g.cat1),
          ('rises within S2', lambda g: ~g.cat1 & g.cat2 & (g.change_segment_Ctrl1A2 == 1) & (g.change_Ctrl1A2 > 0)),
          ('falls within S2', lambda g: ~g.cat1 & g.cat2 & (g.change_segment_Ctrl1A2 == 1) & (g.change_Ctrl1A2 < 0)),
          ('falls within S3', lambda g: ~g.cat1 & g.cat2 & (g.change_segment_Ctrl1A2 == 2) & (g.change_Ctrl1A2 < 0))]
MISSED = '#F4B6B0'


def centred(gene, genes, grid, curves):
    """Curve and cluster staircase of one gene, both mice averaged, as log2 fold from its bulk PT mean."""
    b12, b23 = genes.loc[gene, ['b12', 'b23']]
    seg = np.digitize(grid, [b12, b23])
    bulk = np.mean([curves[m]['bulk'][gene] for m in MICE])
    curve = np.mean([curves[m]['curve'][gene] for m in MICE], 0) - bulk
    steps = np.mean([curves[m]['clusters'][gene] for m in MICE], 0)[seg] - bulk
    # the spline can dive where counts run out; keep it inside the range of this gene's own bin means
    bins = np.concatenate([curves[m]['bins'][gene] for m in MICE]) - bulk
    bins = bins[np.isfinite(bins)]
    return np.clip(curve, bins.min() - 0.25, bins.max() + 0.25), steps


def panel_b_lines(plt, genes, grid, curves, out):
    fig = plt.figure(figsize=(SIZE['b'][0] / 100, SIZE['b'][1] / 100))
    rows = [(name, genes.index[rule(genes)].tolist()) for name, rule in GROUPS]
    rows = [(n, g) for n, g in rows if g]
    h, gap, top = 0.78 / len(rows), 0.20 / len(rows), 0.96
    counts = {}
    for i, (name, members) in enumerate(rows):
        ax = fig.add_axes([0.14, top - (i + 1) * h - i * gap + gap * 0.35, 0.82, h])
        pairs = [centred(g, genes, grid, curves) for g in members]
        cur, stp = np.stack([c for c, _ in pairs]), np.stack([s for _, s in pairs])
        for c in cur:
            ax.plot(grid, c, color=SECTION_COLOURS['Ctrl1A4'], lw=0.5, alpha=0.35, zorder=1)
        mc, ms = np.median(cur, 0), np.median(stp, 0)
        ax.fill_between(grid, mc, ms, color=MISSED, lw=0, zorder=2)
        ax.plot(grid, ms, color=CLUSTER, lw=2.0, zorder=3, solid_capstyle='butt')
        ax.plot(grid, mc, color=SECTION_COLOURS['Ctrl1A2'], lw=2.2, zorder=4)
        ax.axhline(0, color=BULK, lw=0.8, ls=(0, (2, 1.5)), zorder=0)
        both = np.r_[mc, ms]
        b12, b23 = genes.loc[members, 'b12'].mean(), genes.loc[members, 'b23'].mean()
        gene_band(ax, b12, b23, min(-2, np.floor(both.min() - 1)), max(2, np.ceil(both.max() + 1)))
        ax.set_title(f'{name}  ·  {len(members)} genes', fontsize=7.5, color=INK, pad=3, loc='left')
        ax.set_ylabel('log2 fold from PT mean')
        counts[name] = len(members)
    fig.savefig(out / 'f4_b.png', facecolor='white')
    plt.close(fig)
    return counts


def gene_band(ax, b12, b23, lo, hi):
    ax.set_ylim(lo, hi)
    pad = (hi - lo) * 0.05
    segment_band(ax, b12, b23, lo - 2.4 * pad, pad)
    for b in (b12, b23):
        ax.axvline(b, color='#CCCCCC', lw=0.5, zorder=0)
    ax.set_xlim(0, 1)
    ax.set_xticks([])


# ---------------------------------------------------------------------------------------------------- c, d
def panel_examples(plt, key, gene_list, genes, grid, curves, out):
    fig = plt.figure(figsize=(SIZE[key][0] / 100, SIZE[key][1] / 100))
    for i, g in enumerate(gene_list[:3]):
        ax = fig.add_axes([0.08 + i * 0.315, 0.13, 0.25, 0.74])
        gene_axes(ax, g, genes.loc[g], grid, curves, ('bins', 'clusters', 'curve'), ylim_of(g, curves))
        ax.set_title(g, fontsize=8, style='italic', color=INK, pad=3)
        if i == 0:
            ax.set_ylabel('log2 counts per 10k')
    fig.savefig(out / f'f4_{key}.png', facecolor='white')
    plt.close(fig)


def main(argv=None):
    args = parse_args(argv)
    root = Path(args.results_root) / args.out
    out = root / 'panels'
    out.mkdir(parents=True, exist_ok=True)
    plt = style()
    genes, grid, curves = load(root)
    c1, c2 = examples(genes)
    a = panel_a(plt, genes, grid, curves, out)
    # presentation choice made after the screen: panel a shows the top gene, so c and d skip it (no duplicate)
    c1, c2 = [g for g in c1 if g != a][:3], [g for g in c2 if g != a][:3]
    b = panel_b_lines(plt, genes, grid, curves, out)
    if c1:
        panel_examples(plt, 'c', c1, genes, grid, curves, out)
    if c2:
        panel_examples(plt, 'd', c2, genes, grid, curves, out)
    pd.Series({'panel_a': a, 'cat1_examples': ','.join(c1), 'cat2_examples': ','.join(c2), **b}).to_csv(
        out / 'panels.csv', header=False)
    print('a:', a, '| c:', c1, '| d:', c2, '| b rows:', b)


if __name__ == '__main__':
    main()
