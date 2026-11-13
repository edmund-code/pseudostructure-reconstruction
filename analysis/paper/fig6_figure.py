"""Fig. 6 on the same canvas and in the same order as Fig. 5 (analysis/paper/fig5_figure.py), with control and AKI
in place of mouse and human:
 a  one control and one IR section, PT tubules coloured by pseudostructure
 b  Havcr1 (KIM-1) expression in the same tubules
 c  Havcr1 against pseudostructure: binned means (points) and fitted curves (lines) in log2(CP10k + 0.5)
 d  one gene from each gene shape cluster (Hmox1, Slc6a13, Slc22a2, Slc6a8)
 e  three genes whose healthy within-segment gradient injury erases (Cyp4a14, Agt, Cyp24a1)
 f  the four gene shape clusters of shape_clusters.py: mean z-scored shape per section and per condition
 g  one called pathway per pathway group, drawn as in f
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
from fig6_panels import CONTROLS, CTRL, INJURED, IR, MICE, reference_x  # noqa: E402

LEAD = 'Havcr1'
# d: one gene from each gene shape cluster (shape_clusters.py), in cluster order
# (COAST fit: the earlier gene of a cluster where it is still there, otherwise its transporter with the largest change)
GENES_D = ['Slc26a1', 'Hmox1', 'Slc6a13', 'Slc6a8']
# e: healthy within-segment gradients injury erases (Category A examples; the controls carry the gradient)
GENES_E = ['Cyp4a14', 'Agt', 'Cyp24a1']
# f: genes named on each shape cluster; g: per pathway group, the called set with the smallest competitive p
# (COAST fit: f keeps the earlier named genes that are still in a cluster, otherwise its two largest changes)
EXAMPLES = {1: ['Khk', 'Aamdc'], 2: ['Hmox1', 'Tk1'], 3: ['Slc6a13', 'Slc22a12'], 4: ['Acadm', 'Havcr1']}
HIGHLIGHT = ['MSigDB_Hallmark_2020::E2F Targets', 'MSigDB_Hallmark_2020::Myc Targets V1',
             'Reactome_2022::Citric Acid (TCA) Cycle And Respiratory Electron Transport R-HSA-1428517',
             'KEGG_2019_Mouse::Lysosome']
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
    grid, curves, pos = load(screen)
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
    F.shape_rows(fig, R, panel, root / 'shape_clusters', GROUPS, SHADE, b12, b23, HIGHLIGHT, EXAMPLES)
    fig.savefig(root / 'fig6.png', facecolor='white', dpi=400)
    fig.savefig(root / 'fig6.pdf', facecolor='white')
    plt.close(fig)
    F.wireframe(plt, W, H, boxes, root / 'fig6_wireframe.png')


if __name__ == '__main__':
    main()
