"""Fig. 6 on the same canvas and in the same order as Fig. 5 (analysis/paper/fig5_figure.py), with control and AKI
in place of mouse and human:
 a  one control and one IR section, PT tubules coloured by pseudostructure
 b  Havcr1 (KIM-1) expression in the same tubules
 c  Havcr1 against pseudostructure: binned means (points) and fitted curves (lines) in log2(CP10k + 0.5)
 d  four transporters and enzymes whose loss is positional (Slc34a1, Slc7a13, Slc22a6, Acadm)
 e  three genes whose healthy within-segment gradient injury erases (Cyp4a14, Agt, Cyp24a1)
 f  transporter pairs in the poster style: mean log2(1 + CP10k) curves per condition, axis from 0, no switch lines
 g  each shifted switch as an arrow from its control to its AKI position, from the curves and from straight lines
    through the S1/S2/S3 cluster means
Both conditions are measured on the one mouse probe panel, so unlike Fig. 5 the two columns of f share a y axis.
Writes fig6.png, fig6.pdf and fig6_wireframe.png to <results-root>/paper_figures/fig6/.
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
import fig5_figure as F  # noqa: E402
from fig5_panels import XRANGE  # noqa: E402
from fig6_panels import CONTROLS, CTRL, INJURED, IR, MICE, PROTEIN, reference_x  # noqa: E402

LEAD = 'Havcr1'
# d: each loses its own specialised transport or metabolism, at different positions (fig6-aki-along-pt.md)
GENES_D = ['Slc34a1', 'Slc7a13', 'Slc22a6', 'Acadm']
# e: healthy within-segment gradients injury erases (Category A examples; the controls carry the gradient)
GENES_E = ['Cyp4a14', 'Agt', 'Cyp24a1']
# f: shifted pairs whose partners both pass the neighbour rule, one per substrate class
G_PAIRS = [('Slc6a19', 'Slc6a13', 'Amino acids'), ('Slc22a6', 'Slc22a2', 'Drug secretion'),
           ('Slc39a5', 'Slc39a8', 'Zinc')]
SHADE = {'Ctrl1A2': CTRL, 'Ctrl1A4': '#7FA8D8', 'IR2A2': IR, 'IR2A4': '#F0A978'}
KEY = {'Ctrl1A2': 'Control 1', 'Ctrl1A4': 'Control 2', 'IR2A2': 'AKI 1', 'IR2A4': 'AKI 2'}
EARLY = {'control': '#3A8DD0', 'AKI': '#F0A978'}
LATE = {'control': '#1B3A5C', 'AKI': '#8A3B08'}
GROUPS = {'control': (CONTROLS, CTRL, EARLY['control'], LATE['control']),
          'AKI': (INJURED, IR, EARLY['AKI'], LATE['AKI'])}


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--results-root', default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT', str(REPO / 'results')))
    p.add_argument('--screen', default='screen', help='screen folder under paper_figures/fig6')
    return p.parse_args(argv)


def load(screen: Path):
    """curves.npz holds the genes the screen kept, with the bins Fig. 5's curve panels draw."""
    z = np.load(screen / 'curves.npz', allow_pickle=True)
    names = list(z['genes'])
    pos = {g: i for i, g in enumerate(names)}
    curves = {m: {part: z[f'{part}|{m}'] for part in ('curve', 'bins', 'bin_x')} for m in MICE}
    return z['grid'], curves, pos


def main(argv=None):
    args = parse_args(argv)
    rroot = Path(args.results_root)
    root = rroot / 'paper_figures/fig6'
    screen = root / args.screen
    plt = F.style()
    from matplotlib.colors import LinearSegmentedColormap
    F.EXPR_CMAP = LinearSegmentedColormap.from_list('expr', ['#FBEDE2', '#D9692B', '#4A1D02'])
    genes = pd.read_csv(screen / 'genes.csv').set_index('gene')
    grid, curves, pos = load(screen)
    hp = pd.read_csv(screen / 'handoffs.csv')
    axes = pd.read_csv(screen / 'axes.csv')
    b12, b23 = float(axes.b12.median()), float(axes.b23.median())
    # each mouse's drawn range: the sampled range of the handoff genes' axis (axes.csv carries no x range here)
    XRANGE.update({m: (float(curves[m]['bin_x'][pos[LEAD]].min()), float(curves[m]['bin_x'][pos[LEAD]].max()))
                   for m in MICE})

    W, H = 180, 238
    fig, R, panel, boxes = F.canvas(plt, W, H)
    panel('a', 1, 1, 89, 50)
    panel('b', 91, 1, 179, 50)
    F.panel_maps(fig, rroot, R(8, 7, 78, 34, rect=True), R(98, 7, 78, 34, rect=True), run='aki',
                 shown=(('Ctrl1A2', 'control'), ('IR2A2', 'AKI')), gene=LEAD, pt=reference_x(rroot))
    F.curve_rows(fig, R, panel, curves, pos, grid, b12, b23, LEAD, GENES_D, GENES_E, MICE, SHADE, KEY)
    panel('f', 1, 140, 179, 188)
    F.draw_handoffs(fig, R, 141, curves, pos, grid, b12, b23, pairs=G_PAIRS, groups=GROUPS, protein=PROTEIN,
                    shared_y=True)
    panel('g', 1, 191, 179, 237)
    # the handoff screen applies no neighbour filter: Slc5a3 (SMIT1) is 24x higher in thick limb than in PT, so
    # that pair is thick-limb signal and is dropped here (post hoc, as fig6_panels.neighbour_free does)
    keep = genes.index[genes.universe.astype(bool) & genes.neighbour_ok.astype(bool)]
    d = F.panel_shift(fig, R(40, 198, 95, 32, rect=True), genes, hp, b12, b23, groups=GROUPS, protein=PROTEIN,
                      title='Switch point, control → AKI', keep=keep)
    d.to_csv(root / 'fig6_switches.csv', index=False)
    fig.savefig(root / 'fig6.png', facecolor='white', dpi=400)
    fig.savefig(root / 'fig6.pdf', facecolor='white')
    plt.close(fig)
    F.wireframe(plt, W, H, boxes, root / 'fig6_wireframe.png')
    print(d.round(3).to_string())


if __name__ == '__main__':
    main()
