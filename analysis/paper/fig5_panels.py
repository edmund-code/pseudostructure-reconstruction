"""Fig. 5 panels (main) and their Extended Data companions, drawn from fig5_screen.py and fig5_interpret.py.

Main figure (make_fig5.py):
 reversals      reversed gradients, both limbs externally confirmed, ranked by clinical relevance (fig56_rank.py)
 contrast       every gene's S1->S3 gradient, mouse against human: drug targets that keep their direction beside the
                reversals
 pathway_curves members' mean curves of clinically relevant pathways (display only: the set screen called none)
 handoff_rows   transporter handoffs in the poster layout, mouse and human side by side, with switch positions
Extended Data (make_ed_fig_hvm.py):
 maps           one mouse and one human section, PT tubules coloured by position on the shared axis
 coverage       where each section's PT tubules sit on the axis (density per section)
 conservation   how much of a gradient chosen in one section each other section keeps
 flattening     every assessable gene's gradient, human against mouse, with the reversed and steeper calls
 all_calls      every reversed and steeper call along the PT, mouse beside human
 pathway_null   the set screen's sorted p against both relabelings (it called no set)
 handoff_moves  every handoff found in both species: switch position in mouse against human
 checks         the admixture bound and cortical check of every called gene
(relocated and pathways draw set calls; the frozen screen made none, so they are skipped.)

Curves are drawn centred on each section's own mean (levels are not comparable across the two probe panels).
Reads <results-root>/paper_figures/fig5/screen/ and writes PNGs to <results-root>/paper_figures/fig5/panels/.
Protocol: docs/paper/fig5-human-vs-mouse-along-pt.md.
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
from fig6_panels import BAND, LIBRARY_ORDER, PROTEIN, SEGS, covered, pathway_name, segment_band  # noqa: E402

MOUSE, HUMAN = ['Ctrl1A2', 'Ctrl1A4'], ['HUK1_COR1', 'HUK1_MED1']
SECTIONS = MOUSE + HUMAN
PT = {'PT-S1': 0, 'PT-S2': 1, 'PT-S3': 2}
MOUSE_C, HUMAN_C = '#2F6DB5', '#C51B7D'
EARLY = {'mouse': '#3A8DD0', 'human': '#E58BBE'}
LATE = {'mouse': '#1B3A5C', 'human': '#7A0F4D'}
REVERSED_C = '#D95F02'
SIZE = {'reversals': (1160, 232), 'reversals_grid': (760, 352), 'contrast': (390, 352), 'pathway_curves': (1160, 232), 'rises': (1160, 232), 'falls': (1160, 232), 'relocated': (870, 232), 'handoff_rows': (1160, 412), 'steeper': (1160, 232),
        'maps': (360, 242), 'coverage': (380, 242), 'conservation': (380, 242), 'flattening': (380, 312),
        'pathways': (780, 312), 'handoff_moves': (580, 252), 'checks': (560, 252), 'all_calls': (760, 612),
        'pathway_null': (380, 242)}
MAIN = ['reversals_grid', 'contrast', 'pathway_curves', 'handoff_rows']
EXTENDED = ['maps', 'coverage', 'conservation', 'flattening', 'all_calls', 'pathway_null', 'handoff_moves', 'checks']
MIN_FINITE_BINS = 16

# examples are chosen after the screen and the follow-up (presentation choice, recorded in the protocol doc)
# reversals whose two limbs both reproduce externally (mouse microdissection in both sexes, human atlas donors;
# docs/paper/fig5-human-vs-mouse-along-pt.md): rising in mouse and falling in human, and the reverse
RISES = [('Mme', 'Mme'), ('Slc17a3', 'Slc17a3'), ('Hpgd', 'Hpgd'), ('Lpar6', 'Lpar6')]
FALLS = [('Dcxr', 'Dcxr'), ('Slc5a11', 'Slc5a11'), ('Pitpnc1', 'Pitpnc1'), ('Sorbs1', 'Sorbs1')]
# main panel a (2026-11-09): the reversals with a drug, GWAS or disease hook, all both-limbs confirmed
# (docs/paper/fig5-6-highlights.md); Lpar6 and Pitpnc1 stay in Extended Data (no kidney literature)
REVERSED = [('Mme', 'Mme'), ('Slc17a3', 'Slc17a3'), ('Hpgd', 'Hpgd'), ('Slc5a11', 'Slc5a11'), ('Dcxr', 'Dcxr')]
# contrast panel: approved-drug or disease-panel genes whose gradient keeps its sign, |E| >= 0.5, in all four sections
CONSERVED_TARGETS = ['Slc5a2', 'Slc22a12', 'Slc22a8', 'Ace', 'Folh1', 'Dpep1', 'Pah', 'Slc2a2', 'Slc4a4', 'Gatm']
# pathway display (no call: the frozen set screen called none; human gradients are attenuated, so a flat human curve
# is not evidence of a missing gradient)
PATHWAY_CURVES = [('Fatty-acid oxidation', 'Reactome_2022::Mitochondrial Fatty Acid Beta-Oxidation R-HSA-77289'),
                  ('Gluconeogenesis', 'Reactome_2022::Gluconeogenesis R-HSA-70263'),
                  ('Lysosome', 'KEGG_2019_Mouse::Lysosome'),
                  ('MHC class II', 'Reactome_2022::MHC Class II Antigen Presentation R-HSA-2132295')]
RELOCATED: list = []          # (label, gene-set key); the frozen set screen called none
H_PAIRS = [('Slc5a2', 'Slc5a1', 'Glucose'), ('Slc6a19', 'Slc6a18', 'Neutral amino acids'),
           ('Slc7a7', 'Slc7a13', 'Cationic amino acids'), ('Slc22a23', 'Slc22a12', 'Urate'),
           ('Slc22a23', 'Slc22a2', 'Organic cations')]
PROTEIN = {**PROTEIN, 'Slc5a12': 'SMCT2', 'Slc5a8': 'SMCT1', 'Slc5a2': 'SGLT2', 'Slc5a1': 'SGLT1',
           'Slc5a10': 'SGLT5', 'Slc5a11': 'SMIT2', 'Slc6a18': 'B0AT3', 'Slc7a8': 'LAT2', 'Slc7a7': 'y+LAT1',
           'Slc22a8': 'OAT3', 'Slc22a7': 'OAT2', 'Slc22a12': 'URAT1', 'Slc22a13': 'OAT10', 'Slc22a5': 'OCTN2',
           'Slc22a23': 'SLC22A23', 'Slc16a13': 'MCT13', 'Slc16a10': 'MCT10', 'Slc6a6': 'TauT', 'Slc17a3': 'NPT4',
           'Slc10a2': 'ASBT', 'Mme': 'MME', 'Hpgd': '15-PGDH', 'Dcxr': 'DCXR', 'Ace': 'ACE', 'Folh1': 'PSMA',
           'Dpep1': 'DPEP1', 'Pah': 'PAH', 'Slc2a2': 'GLUT2', 'Slc4a4': 'NBCe1', 'Gatm': 'GATM', 'Slc22a2': 'OCT2', 'Slc17a1': 'NPT1', 'Abcc4': 'MRP4', 'Abcc2': 'MRP2', 'Slc12a1': 'NKCC2'}


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--results-root', default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT', str(REPO / 'results')))
    p.add_argument('--out', default='paper_figures/fig5', help='relative to the results root')
    p.add_argument('--panels', nargs='+', default=None, help='panel names (default: all)')
    return p.parse_args(argv)


def figure(plt, key):
    w, h = SIZE[key]
    return plt.figure(figsize=(w / 100, h / 100))


def species_of(m):
    return 'mouse' if m in MOUSE else 'human'


# ---------------------------------------------------------------------------------------------------- data
def load(screen: Path):
    genes = pd.read_csv(screen / 'genes.csv', index_col=0)
    z = np.load(screen / 'all_curves.npz', allow_pickle=True)
    names = list(z['genes'])
    pos = {g: i for i, g in enumerate(names)}
    curves = {m: {part: z[f'{part}|{m}'] for part in ('curve', 'bins', 'bin_x')} for m in SECTIONS}
    return genes, z['grid'], curves, pos


def centred(gene, m, curves, pos, grid):
    """A section's curve and bin means centred on the mean of its curve over the section's sampled range, curve
    clipped to the bin range (+-0.25); NaN outside the range."""
    i = pos[gene]
    bx, bm, c = curves[m]['bin_x'][i], curves[m]['bins'][i].astype(float), curves[m]['curve'][i].astype(float)
    k = covered(grid, bx)
    fin = np.isfinite(bm)
    centre = c[k].mean()
    if fin.any():
        c = np.clip(c, bm[fin].min() - 0.25, bm[fin].max() + 0.25)
    c = np.where(k, c - centre, np.nan)
    if fin.sum() < MIN_FINITE_BINS:
        c = np.full_like(c, np.nan)
    return c, bx, bm - centre


def set_curves(members, genes, curves, pos, grid):
    """Per section, the mean over the set's assessable members of their centred curves."""
    ok = [g for g in members if g in genes.index and genes.loc[g, 'assessable'] and g in pos]
    out = {}
    for m in SECTIONS:
        c = np.array([centred(g, m, curves, pos, grid)[0] for g in ok])
        # drawn only where >= 80% of members have a curve: past that, members' coverage ends one by one and the
        # mean jumps
        full = np.isfinite(c).mean(0) >= 0.8
        out[m] = np.where(full, np.nanmean(np.where(np.isfinite(c), c, np.nan), 0) if full.any() else np.nan, np.nan)
    return out, len(ok)


def decorate(ax, b12, b23, lo, hi):
    ax.set_ylim(lo, hi)
    pad = (hi - lo) * 0.06
    segment_band(ax, b12, b23, lo - 2.2 * pad, pad)
    for b in (b12, b23):
        ax.axvline(b, color='#CCCCCC', lw=0.5, zorder=0)
    ax.set_xlim(0, 1)
    ax.set_xticks([])


# ---------------------------------------------------------------------------------------------------- curves
def panel_examples(plt, genes, grid, curves, pos, out, b12, b23, key, examples, sets=None, ncols=None):
    """Small multiples: each example's centred curves (mouse blues, human magentas), bin means for single genes;
    one row, or a grid of ncols columns."""
    fig = figure(plt, key)
    n = max(len(examples), 1)
    ncols = ncols or n
    nrows = -(-n // ncols)
    rh = 0.87 / nrows
    for i, (label, what) in enumerate(examples):
        r, ci = divmod(i, ncols)
        ax = fig.add_axes([0.07 + ci * (0.93 / ncols), 0.13 / nrows + (nrows - 1 - r) * rh + (0.02 if nrows > 1 else 0),
                           0.93 / ncols - 0.06, rh - 0.13])
        if isinstance(what, str) and what in pos:
            vals = []
            for m in SECTIONS:
                c, bx, bm = centred(what, m, curves, pos, grid)
                ok = np.isfinite(bm)
                ax.scatter(bx[ok], bm[ok], s=3, color=SECTION_COLOURS[m], lw=0, alpha=0.7, zorder=3)
                ax.plot(grid, c, color=SECTION_COLOURS[m], lw=1.2, zorder=4)
                vals.append(bm[ok])
            v = np.concatenate(vals)
            lo, hi = v.min() - 0.4, v.max() + 0.4
            ax.set_title(label, fontsize=7.5, style='italic', color=INK, pad=2)
        else:
            mc, k = set_curves(sets[what], genes, curves, pos, grid)
            inner = (grid >= 0.03) & (grid <= 0.97)     # the free spline ends hook (fig4 doc, PT-end uptick)
            mc = {m: np.where(inner, mc[m], np.nan) for m in SECTIONS}
            for m in SECTIONS:
                ax.plot(grid, mc[m], color=SECTION_COLOURS[m], lw=1.2)
            v = np.concatenate([mc[m][np.isfinite(mc[m])] for m in SECTIONS])
            lo, hi = v.min() - 0.1, v.max() + 0.1
            ax.set_title(f'{label} ({k} genes)', fontsize=7.5, color=INK, pad=2)
        ax.axhline(0, color='#E2E2E2', lw=0.5, zorder=0)
        decorate(ax, b12, b23, lo, hi)
        if ci == 0:
            ax.set_ylabel('log2, centred')
    fig.savefig(out / f'f5_{key}.png', facecolor='white')
    plt.close(fig)


def panel_contrast(plt, genes, out):
    """Every screened, assessable gene's S1->S3 gradient E, mouse mean against human mean; named drug targets that keep
    their direction (ink) beside the reversed calls (orange). Axes differ: human gradients are about 0.6x mouse."""
    fig = figure(plt, 'contrast')
    ax = fig.add_axes([0.15, 0.14, 0.80, 0.80])
    ok = genes.screened & genes.assessable
    em = genes.loc[ok, [f'E_{m}' for m in MOUSE]].mean(1)
    eh = genes.loc[ok, [f'E_{m}' for m in HUMAN]].mean(1)
    ax.axhline(0, color='#CCCCCC', lw=0.6, zorder=0)
    ax.axvline(0, color='#CCCCCC', lw=0.6, zorder=0)
    ax.scatter(em, eh, s=3, color='#CFCFCF', lw=0, zorder=1, rasterized=True)
    named = [(g, INK) for g in CONSERVED_TARGETS] + [(g, REVERSED_C) for g, _ in REVERSED]
    for g, col in named:
        assert g in em.index, f'{g} not screened and assessable'
        x, y = em[g], eh[g]
        ax.scatter([x], [y], s=16, color=col, lw=0, zorder=3)
        dx, dy, ha = LABEL_SIDE.get(g, (3, 2, 'left'))
        ax.annotate(PROTEIN.get(g, g.upper()), (x, y), xytext=(dx, dy), textcoords='offset points', fontsize=5.8,
                    color=col, ha=ha, zorder=4)
    ax.set_xlim(-7, 7)
    ax.set_ylim(-2, 2)
    ax.set_xlabel('mouse gradient (log2, S1 → S3)', labelpad=1)
    ax.set_ylabel('human gradient', labelpad=1)
    fig.savefig(out / 'f5_contrast.png', facecolor='white')
    plt.close(fig)


LABEL_SIDE = {'Pah': (-3, 2, 'right'), 'Slc22a8': (-3, -2, 'right'), 'Folh1': (3, -8, 'left'),
              'Slc2a2': (3, -8, 'left'), 'Slc5a2': (3, -8, 'left'), 'Slc5a11': (3, -6, 'left')}


# ---------------------------------------------------------------------------------------------------- handoffs
def handoff_range(screen: Path, hp: pd.DataFrame, genes: pd.DataFrame) -> dict:
    """Each section's 2.5-97.5% x range on the axis the screen used for every pair (the first handoff gene's)."""
    first = sorted(g for g in set(hp.early) | set(hp.late) if g in genes.index)[0]
    ax = pd.read_csv(screen / 'axes.csv', dtype={'pattern': str}).set_index('pattern').loc[genes.loc[first, 'pattern']]
    return {m: (ax[f'x_lo_{m}'], ax[f'x_hi_{m}']) for m in SECTIONS}


XRANGE: dict = {}


def scaled(gene, m, curves, pos, grid):
    """A section's curve min-max scaled over the x range the switch was computed on (fig5_screen.py)."""
    from pseudospace.handoffs import scaled_profile

    lo, hi = XRANGE[m]
    return scaled_profile(np.where((grid >= lo) & (grid <= hi), curves[m]['curve'][pos[gene]], np.nan))


def panel_handoff_rows(plt, grid, curves, pos, hp, out, b12, b23):
    fig = figure(plt, 'handoff_rows')
    n = len(H_PAIRS)
    top, bottom, gap = 0.86, 0.09, 0.025
    h = (top - bottom - gap * (n - 1)) / n
    cols = {'mouse': (0.20, 0.31), 'human': (0.545, 0.31)}
    edges = ((0, b12), (b12, b23), (b23, 1))
    for name, (x0, w) in cols.items():
        fig.text(x0 + w / 2, 0.975, name, ha='center', va='top', fontsize=9, weight='bold',
                 color=MOUSE_C if name == 'mouse' else HUMAN_C)
        for (lo, hi), seg in zip(edges, SEGS):
            fig.text(x0 + w * (lo + hi) / 2, 0.895, f' {seg} ', ha='center', va='bottom', fontsize=6.5, weight='bold',
                     color='white', bbox=dict(boxstyle='round,pad=0.25', fc=COLOURS[seg], ec='none'))
    fig.lines.append(plt.Line2D([0.005, 0.03], [0.965, 0.965], transform=fig.transFigure, color=MOUSE_C, lw=1.2,
                                ls=(0, (3, 1.5))))
    fig.text(0.035, 0.965, 'switch in mouse', fontsize=6.5, color=MUTED, va='center')
    fig.lines.append(plt.Line2D([0.005, 0.03], [0.925, 0.925], transform=fig.transFigure, color=HUMAN_C, lw=1.2,
                                ls=(0, (1, 1.2))))
    fig.text(0.035, 0.925, 'switch in human', fontsize=6.5, color=MUTED, va='center')
    sx0, sw = 0.885, 0.105
    fig.text(sx0 + sw / 2, 0.895, 'switch\n(95% interval)', ha='center', va='bottom', fontsize=6, color=MUTED)
    for r, (e, l, substrate) in enumerate(H_PAIRS):
        y0 = top - (r + 1) * h - r * gap
        row = hp[(hp.early == e) & (hp.late == l)].iloc[0]
        xs = {sp: np.nanmean([row[f'x_{m}'] for m in secs]) for sp, secs in (('mouse', MOUSE), ('human', HUMAN))}
        fig.text(0.005, y0 + h * 0.80, substrate, fontsize=7.5, weight='bold', color=INK, va='center')
        t = fig.text(0.005, y0 + h * 0.50, e, fontsize=7, style='italic', weight='bold', color=EARLY['mouse'],
                     va='center')
        right = t.get_window_extent(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted()).x1
        fig.text(right + 0.003, y0 + h * 0.50, f'→ {l}', fontsize=7, style='italic', weight='bold',
                 color=LATE['mouse'], va='center')
        fig.text(0.005, y0 + h * 0.20, f'{PROTEIN.get(e, e)} → {PROTEIN.get(l, l)}', fontsize=6.3, color=MUTED,
                 va='center')
        for name, (x0, w) in cols.items():
            ax = fig.add_axes([x0, y0, w, h])
            for (lo, hi), seg in zip(edges, SEGS):
                ax.axvspan(lo, hi, color=BAND[seg], lw=0, zorder=0)
            secs = MOUSE if name == 'mouse' else HUMAN
            for gene, col in ((e, EARLY[name]), (l, LATE[name])):
                per = [scaled(gene, m, curves, pos, grid) for m in secs]
                mean = np.mean(per, 0)          # drawn only where both sections of the species have tubules
                for c in per:
                    ax.plot(grid, np.where(np.isfinite(mean), c, np.nan), color=col, lw=0.5, alpha=0.6, zorder=2)
                ax.fill_between(grid, 0, mean, color=col, alpha=0.12, lw=0, zorder=1)
                ax.plot(grid, mean, color=col, lw=1.6, zorder=3)
            if np.isfinite(xs['mouse']):
                ax.axvline(xs['mouse'], color=MOUSE_C, lw=1.1, ls=(0, (3, 1.5)), zorder=4)
            if name == 'human' and np.isfinite(xs['human']):
                ax.axvline(xs['human'], color=HUMAN_C, lw=1.3, ls=(0, (1, 1.2)), zorder=4)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1.08)
            ax.set_xticks([])
            ax.set_yticks([0, 1] if name == 'mouse' else [])
            for sp in ('top', 'right'):
                ax.spines[sp].set_visible(False)
            if r == n - 1:
                ax.set_xlabel('start → end of PT', labelpad=2)
            if name == 'mouse' and r == n // 2:
                ax.set_ylabel('scaled expression', labelpad=2)
        # switch positions: each section's crossing with its bootstrap interval
        sx = fig.add_axes([sx0, y0 + h * 0.1, sw, h * 0.8])
        for (lo, hi), seg in zip(edges, SEGS):
            sx.axvspan(lo, hi, color=BAND[seg], lw=0, zorder=0)
        for j, m in enumerate(SECTIONS):
            yy = 3 - j
            v, a, b = row[f'x_{m}'], row[f'x_lo_{m}'], row[f'x_hi_{m}']
            col = MOUSE_C if m in MOUSE else HUMAN_C
            if np.isfinite(v):
                sx.plot([a, b], [yy, yy], color=col, lw=0.9, zorder=2)
                sx.scatter([v], [yy], s=10, color=col, lw=0, zorder=3)
            else:
                sx.text(0.5, yy, 'no switch', fontsize=5.2, color=col, ha='center', va='center')
        sx.set_xlim(0, 1)
        sx.set_ylim(-0.7, 3.7)
        sx.set_yticks([])
        sx.set_xticks([0, 1] if r == n - 1 else [])
        if r == n - 1:
            sx.set_xticklabels(['start', 'end'], fontsize=5.5)
        for sp in ('top', 'right', 'left'):
            sx.spines[sp].set_visible(False)
    fig.savefig(out / 'f5_handoff_rows.png', facecolor='white')
    plt.close(fig)


def panel_handoff_moves(plt, hp, out, b12, b23):
    """Every pair with a switch in all four sections: mouse switch against human switch (lines: each section's
    spread); shifted pairs coloured by direction and listed beside, with the pairs that switch in mouse only."""
    later, earlier = '#6A4BB0', '#2F6DB5'
    fig = figure(plt, 'handoff_moves')
    ax = fig.add_axes([0.10, 0.15, 0.40, 0.80])
    e = hp[hp.eligible].copy()
    e['mx'] = e[[f'x_{m}' for m in MOUSE]].mean(1)
    e['hx'] = e[[f'x_{m}' for m in HUMAN]].mean(1)
    for bnd in (b12, b23):
        ax.axvline(bnd, color='#E2E2E2', lw=0.5, zorder=0)
        ax.axhline(bnd, color='#E2E2E2', lw=0.5, zorder=0)
    ax.plot([0, 1], [0, 1], color='#BBBBBB', lw=0.6, zorder=0)
    for _, r in e.sort_values('shifted').iterrows():
        col = (later if r['shift'] > 0 else earlier) if r.shifted else '#BDBDBD'
        ax.plot([r.mx, r.mx], [r[f'x_{HUMAN[0]}'], r[f'x_{HUMAN[1]}']], color=col, lw=0.6, alpha=0.7)
        ax.plot([r[f'x_{MOUSE[0]}'], r[f'x_{MOUSE[1]}']], [r.hx, r.hx], color=col, lw=0.6, alpha=0.7)
        ax.scatter([r.mx], [r.hx], s=18 if r.shifted else 8, color=col, lw=0, zorder=3)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel('switch in mouse', labelpad=1)
    ax.set_ylabel('switch in human', labelpad=1)
    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    name = lambda r: f'{PROTEIN.get(r.early, r.early)} → {PROTEIN.get(r.late, r.late)}'  # noqa: E731
    columns = ((0.55, (('later in human', later, e[e.shifted & (e['shift'] > 0)].sort_values('shift', ascending=False)),)),
               (0.78, (('earlier in human', earlier, e[e.shifted & (e['shift'] < 0)].sort_values('shift')),
                       ('mouse only', MUTED, hp[hp.mouse_only]))))
    for x0, groups in columns:
        y = 0.95
        for title, col, d in groups:
            fig.text(x0, y, title, fontsize=6.8, weight='bold', color=col, va='top')
            y -= 0.07
            for _, r in d.iterrows():
                shift = f'  {r["shift"]:+.2f}' if np.isfinite(r['shift']) else ''
                fig.text(x0 + 0.005, y, f'{name(r)}{shift}', fontsize=5.8, color=INK, va='top')
                y -= 0.058
            y -= 0.04
    fig.savefig(out / 'f5_handoff_moves.png', facecolor='white')
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------- context
def reference_x(root: Path) -> pd.DataFrame:
    ref = pd.read_parquet(root / RUNS['hvm'] / '04_pseudospace/pseudospace.parquet')
    pt = ref[ref.segment_class.astype(str).isin(PT) & ref['sample'].astype(str).isin(SECTIONS)].copy()
    pt['sample'] = pt['sample'].astype(str)
    ms = pt['sample'].isin(MOUSE)
    c = np.sort(pt.loc[ms, 'pseudospace'].to_numpy())
    q = (np.arange(len(c)) + 0.5) / len(c)
    pt['x'] = np.interp(pt['pseudospace'], c, q)
    return pt


def panel_maps(plt, root, out):
    import anndata as ad

    pt = reference_x(root)
    obs = ad.read_h5ad(root / RUNS['hvm'] / '02_qc/units_qc.h5ad', backed='r').obs
    ref = pd.read_parquet(root / RUNS['hvm'] / '04_pseudospace/pseudospace.parquet')
    fig = figure(plt, 'maps')
    cmap = plt.get_cmap('viridis')
    for i, s in enumerate(('Ctrl1A2', 'HUK1_COR1')):
        ax = fig.add_axes([0.01 + i * 0.49, 0.13, 0.47, 0.80])
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
        fig.text(0.01 + i * 0.49 + 0.235, 0.95, 'mouse' if i == 0 else 'human', fontsize=7.5, ha='center', va='top',
                 color=MOUSE_C if i == 0 else HUMAN_C)
    cax = fig.add_axes([0.30, 0.07, 0.40, 0.025])
    cax.imshow(np.linspace(0, 1, 256)[None], aspect='auto', cmap=cmap, extent=(0, 1, 0, 1))
    cax.set_yticks([])
    cax.set_xticks([0, 1])
    cax.set_xticklabels(['start of PT', 'end of PT'], fontsize=6)
    cax.tick_params(length=0, pad=1.5)
    for sp in cax.spines.values():
        sp.set_visible(False)
    fig.savefig(out / 'f5_maps.png', facecolor='white', dpi=300)
    plt.close(fig)


def panel_coverage(plt, root, out, b12, b23):
    from scipy.stats import gaussian_kde

    pt = reference_x(root)
    fig = figure(plt, 'coverage')
    ax = fig.add_axes([0.20, 0.16, 0.76, 0.80])
    grid = np.linspace(0, 1, 301)
    for i, s in enumerate(SECTIONS[::-1]):
        x = pt.loc[pt['sample'] == s, 'x'].to_numpy()
        d = gaussian_kde(x, bw_method=0.06)(grid)
        d = d / d.max() * 0.9
        c = SECTION_COLOURS[s]
        ax.fill_between(grid, i, i + d, color=c, alpha=0.35, lw=0)
        ax.plot(grid, i + d, color=c, lw=1.0)
        ax.text(-0.02, i + 0.25, s, ha='right', va='center', fontsize=6.5, color=c)
    for b in (b12, b23):
        ax.axvline(b, color='#BBBBBB', lw=0.5, zorder=0)
    ax.set_ylim(-0.32, len(SECTIONS))
    segment_band(ax, b12, b23, -0.28, 0.12)
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['start of PT', 'end'])
    ax.set_xlabel('pseudostructure (mouse quantile)', labelpad=1)
    fig.savefig(out / 'f5_coverage.png', facecolor='white')
    plt.close(fig)


def panel_conservation(plt, screen, out):
    """Median share of a gradient chosen in one section that each other section keeps (signed)."""
    d = pd.read_csv(screen / 'conservation.csv')
    fig = figure(plt, 'conservation')
    ax = fig.add_axes([0.17, 0.20, 0.80, 0.72])
    for i, sel in enumerate(SECTIONS):
        r = d[d.selected_in == sel].iloc[0]
        for j, m in enumerate(s for s in SECTIONS if s != sel):
            ax.scatter([i + (j - 1) * 0.2], [r[f'kept_{m}']], s=22, color=SECTION_COLOURS[m], lw=0, zorder=3,
                       marker='o' if species_of(m) == species_of(sel) else 'D')
    ax.axhline(0, color='#BBBBBB', lw=0.5)
    ax.axhline(1, color='#E2E2E2', lw=0.5)
    ax.set_xticks(range(len(SECTIONS)))
    ax.set_xticklabels(['mouse 1', 'mouse 2', 'human 1', 'human 2'])
    for t, s in zip(ax.get_xticklabels(), SECTIONS):
        t.set_color(SECTION_COLOURS[s])
    ax.set_xlim(-0.6, len(SECTIONS) - 0.4)
    ax.set_ylabel('gradient kept', labelpad=1)
    ax.set_xlabel('gradient chosen in', labelpad=1)
    fig.savefig(out / 'f5_conservation.png', facecolor='white')
    plt.close(fig)


def panel_flattening(plt, genes, summary, out):
    """Every assessable gene: mean human E against mean mouse E, the genome-wide line, the calls."""
    fig = figure(plt, 'flattening')
    ax = fig.add_axes([0.17, 0.14, 0.78, 0.80])
    g = genes[genes.assessable]
    mx = g[[f'E_{m}' for m in MOUSE]].mean(1)
    hx = g[[f'E_{m}' for m in HUMAN]].mean(1)
    ax.scatter(mx, hx, s=2, color='#C8C8C8', lw=0, rasterized=True)
    for flag, col in (('catR', '#B2182B'), ('catH', HUMAN_C)):
        k = g[flag].to_numpy(bool)
        ax.scatter(mx[k], hx[k], s=9, color=col, lw=0, zorder=3)
    slope = np.mean([s for _, s in summary['flattening_lines']['species']])
    icpt = np.mean([a for a, _ in summary['flattening_lines']['species']])
    lim = np.nanmax(np.abs(np.r_[mx, hx])) * 1.05
    xx = np.array([-lim, lim])
    ax.plot(xx, xx, color='#BBBBBB', lw=0.6, ls=(0, (3, 2)))
    ax.plot(xx, icpt + slope * xx, color=INK, lw=0.9)
    ax.axhline(0, color='#E5E5E5', lw=0.5, zorder=0)
    ax.axvline(0, color='#E5E5E5', lw=0.5, zorder=0)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel('mouse gradient, S3 − S1 (log2)', labelpad=1)
    ax.set_ylabel('human gradient, S3 − S1 (log2)', labelpad=1)
    ax.set_aspect('equal')
    fig.savefig(out / 'f5_flattening.png', facecolor='white')
    plt.close(fig)


def distinct(calls, sets, ok, jaccard=0.3):
    calls = calls.assign(lib=calls.set.str.split('::').str[0].map({l: i for i, l in enumerate(LIBRARY_ORDER)}))
    kept, mem = [], []
    for _, r in calls.sort_values(['lib', 'p', 'n'], ascending=[True, True, False]).iterrows():
        m = {g for g in sets[r.set] if ok.get(g, False)}
        if all(len(m & o) / len(m | o) < jaccard for o in mem):
            mem.append(m)
            kept.append(r)
    return pd.DataFrame(kept).sort_values('d')


def panel_pathways(plt, screen, sets, genes, out):
    """Distinct pathway calls: the species offset d and its four human x mouse pairs; relabelings as open marks."""
    pw = pd.read_csv(screen / 'pathways.csv')
    calls = pw[(pw.split == 'species') & pw.call]
    ok = genes.assessable.to_dict()
    d = distinct(calls, sets, ok) if len(calls) else calls
    fig = figure(plt, 'pathways')
    ax = fig.add_axes([0.46, 0.08, 0.52, 0.90])
    for y, (_, r) in enumerate(d.iterrows()):
        ax.scatter([r[f'pair_{k}'] for k in range(4)], [y] * 4, s=6, color='#999999', lw=0, zorder=2)
        ax.scatter([r.other_0, r.other_1], [y, y], s=10, facecolor='white', edgecolor='#999999', lw=0.6, zorder=2)
        ax.scatter([r.d], [y], s=22, color='#B2182B' if r.d > 0 else '#2166AC', lw=0, zorder=3)
    ax.axvline(0, color='#BBBBBB', lw=0.6)
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels([f'{pathway_name(s)} ({n})' for s, n in zip(d.set, d.n)], fontsize=5.6)
    ax.tick_params(axis='y', length=0)
    ax.set_ylim(len(d) - 0.4, -0.6)
    ax.spines['left'].set_visible(False)
    ax.set_xlabel('offset from the genome-wide trend (log2): earlier in human ← → later in human', labelpad=1)
    fig.savefig(out / 'f5_pathways.png', facecolor='white')
    plt.close(fig)


def panel_checks(plt, screen, out):
    """Per call: the largest gradient admixture could add (x) against the call's smallest |E| (y); open marks fail
    the cortical-S3 check (post hoc). Every call sits far above the diagonal."""
    g = pd.read_csv(screen / 'interpret_genes.csv')
    g['emin'] = g[[f'E_{m}' for m in SECTIONS]].abs().min(1)
    fig = figure(plt, 'checks')
    ax = fig.add_axes([0.17, 0.15, 0.78, 0.80])
    for _, r in g.iterrows():
        col = '#B2182B' if r.category == 'R' else HUMAN_C
        ax.scatter([r.admixture_bound_max], [r.emin], s=16, lw=0.8, zorder=3,
                   facecolor=col if r.cortical_ok else 'white', edgecolor=col)
    lim = max(g.emin.max(), 1) * 1.05
    ax.plot([0, lim], [0, lim], color='#BBBBBB', lw=0.6, ls=(0, (3, 2)))
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel('largest gradient admixture could add (log2)', labelpad=1)
    ax.set_ylabel('smallest gradient of the call (log2)', labelpad=1)
    fig.savefig(out / 'f5_checks.png', facecolor='white')
    plt.close(fig)


def panel_pathway_null(plt, screen, out):
    """Every set's p for the species contrast against the same set's p under each relabeling (sorted): the species
    curve does not rise above the relabelings, so no set is called."""
    pw = pd.read_csv(screen / 'pathways.csv')
    fig = figure(plt, 'pathway_null')
    ax = fig.add_axes([0.17, 0.15, 0.78, 0.80])
    styles = {'species': (INK, 1.4, '-'), 'relabeling 1': ('#999999', 0.9, (0, (3, 1.5))),
              'relabeling 2': ('#999999', 0.9, (0, (1, 1.2)))}
    for split, (col, lw, ls) in styles.items():
        p = np.sort(pw.loc[pw.split == split, 'p'].to_numpy())
        exp = (np.arange(1, len(p) + 1)) / (len(p) + 1)
        ax.plot(-np.log10(exp), -np.log10(p), color=col, lw=lw, ls=ls, label=split)
    lim = ax.get_xlim()[1]
    ax.plot([0, lim], [0, lim], color='#DDDDDD', lw=0.6, zorder=0)
    ax.set_xlabel('expected −log10 p', labelpad=1)
    ax.set_ylabel('observed −log10 p', labelpad=1)
    ax.legend(frameon=False, fontsize=6, loc='upper left')
    fig.savefig(out / 'f5_pathway_null.png', facecolor='white')
    plt.close(fig)


def panel_all_calls(plt, genes, grid, curves, pos, out, b12, b23):
    """Every reversed and steeper-in-human call: each species' mean centred curve, scaled per row to the larger
    species range, mouse block beside human block; rows grouped by category and ordered by the mouse curve's peak."""
    rows = [(g, 'R') for g in genes.index[genes.catR]] + [(g, 'H') for g in genes.index[genes.catH]]
    data = []
    for g, c in rows:
        mc = np.nanmean([centred(g, m, curves, pos, grid)[0] for m in MOUSE], 0)
        hc = np.nanmean([centred(g, m, curves, pos, grid)[0] for m in HUMAN], 0)
        scale = max(np.nanmax(np.abs(mc)), np.nanmax(np.abs(hc)), 1e-9)
        data.append((g, c, mc / scale, hc / scale, float(grid[np.nanargmax(mc)])))
    data.sort(key=lambda d: (d[1] != 'R', d[4]))
    M = np.array([d[2] for d in data])
    H = np.array([d[3] for d in data])
    fig = figure(plt, 'all_calls')
    for j, (mat, name, col) in enumerate(((M, 'mouse', MOUSE_C), (H, 'human', HUMAN_C))):
        ax = fig.add_axes([0.14 + j * 0.42, 0.06, 0.40, 0.88])
        ax.imshow(mat, aspect='auto', cmap='RdBu_r', vmin=-1, vmax=1, extent=(0, 1, len(mat), 0),
                  interpolation='nearest')
        segment_band(ax, b12, b23, len(mat) + 0.2, 0.6)
        ax.set_ylim(len(mat) + 0.9, 0)
        ax.set_xticks([])
        ax.set_yticks(np.arange(len(mat)) + 0.5 if j == 0 else [])
        if j == 0:
            ax.set_yticklabels([f'{d[0]}' for d in data], fontsize=5.6, style='italic')
        ax.tick_params(axis='y', length=0)
        ax.set_title(name, fontsize=8, color=col, pad=3)
        n_r = sum(d[1] == 'R' for d in data)
        ax.axhline(n_r, color='white', lw=2)
        for sp in ax.spines.values():
            sp.set_visible(False)
    fig.savefig(out / 'f5_all_calls.png', facecolor='white')
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------- main
def main(argv=None):
    args = parse_args(argv)
    rroot = Path(args.results_root)
    root = rroot / args.out
    screen = root / 'screen'
    out = root / 'panels'
    out.mkdir(parents=True, exist_ok=True)
    plt = style()
    panels = args.panels or MAIN + EXTENDED
    genes, grid, curves, pos = load(screen)
    summary = json.loads((screen / 'summary.json').read_text())
    sets = json.loads((rroot / 'pt_package_aki/tables/gene_sets.json').read_text())
    hp = pd.read_csv(screen / 'handoffs.csv')
    b12, b23 = float(genes.b12.median()), float(genes.b23.median())
    if 'maps' in panels:
        panel_maps(plt, rroot, out)
    if 'coverage' in panels:
        panel_coverage(plt, rroot, out, b12, b23)
    if 'conservation' in panels:
        panel_conservation(plt, screen, out)
    if 'flattening' in panels:
        panel_flattening(plt, genes, summary, out)
    if 'reversals' in panels:
        panel_examples(plt, genes, grid, curves, pos, out, b12, b23, 'reversals', REVERSED)
    if 'reversals_grid' in panels:   # rising in mouse on top, falling in mouse below
        panel_examples(plt, genes, grid, curves, pos, out, b12, b23, 'reversals_grid', REVERSED, ncols=3)
    if 'contrast' in panels:
        panel_contrast(plt, genes, out)
    if 'pathway_curves' in panels:
        panel_examples(plt, genes, grid, curves, pos, out, b12, b23, 'pathway_curves', PATHWAY_CURVES, sets)
    if 'rises' in panels:
        panel_examples(plt, genes, grid, curves, pos, out, b12, b23, 'rises', RISES)
    if 'falls' in panels:
        panel_examples(plt, genes, grid, curves, pos, out, b12, b23, 'falls', FALLS)
    if 'relocated' in panels and RELOCATED:
        panel_examples(plt, genes, grid, curves, pos, out, b12, b23, 'relocated', RELOCATED, sets)
    if 'pathways' in panels and summary['pathways']['calls']['species']:
        panel_pathways(plt, screen, sets, genes, out)
    if 'pathway_null' in panels:
        panel_pathway_null(plt, screen, out)
    if 'all_calls' in panels:
        panel_all_calls(plt, genes, grid, curves, pos, out, b12, b23)
    XRANGE.update(handoff_range(screen, hp, genes))
    if 'handoff_rows' in panels and H_PAIRS:
        panel_handoff_rows(plt, grid, curves, pos, hp, out, b12, b23)
    if 'handoff_moves' in panels:
        panel_handoff_moves(plt, hp, out, b12, b23)
    if 'checks' in panels and (screen / 'interpret_genes.csv').exists():
        panel_checks(plt, screen, out)


if __name__ == '__main__':
    main()
