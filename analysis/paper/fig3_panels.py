"""Fig. 3 panels drawn from data, sized for the figure layout (make_fig3.py), plus the numbers behind them.

 a  UMAP of the control/AKI fit with the elastic curve, coloured by pseudostructure and by segment label
 b  segment labels along pseudostructure (mouse filled, human outline) and rho per section with a permutation null
 c  one marker per segment, mean expression along pseudostructure: control, AKI, human
 d  robustness reruns (needs fig3_robustness.py output)
 e  trajectory baselines (needs fig3_baselines.py output)
 f  PT pseudostructure against distance from the capsule (needs fig3_depth.py output)

Reads private data under --results-root and writes PNGs and CSVs to <results-root>/paper_figures/fig3/panels/.
Protocol: docs/paper/fig3-validation-experiments.md.
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
from fig3_common import GROUP_ORDER, RUNS, SECTIONS, load_positions, section_rho  # noqa: E402

COLOURS = {'PTS1': '#4C8DC9', 'PTS2': '#5BA55B', 'PTS3': '#C9A46B', 'TL': '#999999', 'AL': '#EE8A2E',
           'DCT': '#E06AA6', 'CNT': '#6FC2B4', 'CD': '#3F5FA8'}
# one marker per segment, the same genes in both species (open decision 5 of the handoff, settled on the
# human-poster set because every gene is measured in both fits)
MARKERS = [('Slc5a2', 'PTS1'), ('Slc13a3', 'PTS2'), ('Slc7a13', 'PTS3'), ('Slc12a1', 'AL'), ('Slc12a3', 'DCT'),
           ('Slc8a1', 'CNT'), ('Aqp2', 'CD')]
SECTION_NAMES = {'Ctrl1A2': 'Ctrl1A2', 'Ctrl1A4': 'Ctrl1A4', 'IR2A2': 'IR2A2', 'IR2A4': 'IR2A4',
                 'HUK1_COR1': 'human 1', 'HUK1_MED1': 'human 2'}
SECTION_COLOURS = {'Ctrl1A2': '#2F6DB5', 'Ctrl1A4': '#7FA7D9', 'IR2A2': '#D9692B', 'IR2A4': '#EFA36E',
                   'HUK1_COR1': '#6A4BB0', 'HUK1_MED1': '#A893D6'}
INK, MUTED = '#222222', '#777777'
# panel boxes of make_fig3.py in SVG units (1 unit = 0.01 inch here)
SIZE = {'a': (760, 345), 'b': (360, 345), 'c': (640, 610), 'd': (480, 170), 'e': (480, 170), 'f': (480, 170)}


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--results-root', default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT', str(REPO / 'results')))
    p.add_argument('--out', default='paper_figures/fig3', help='relative to the results root')
    p.add_argument('--panels', default='abcdef')
    p.add_argument('--permutations', type=int, default=1000)
    return p.parse_args(argv)


def style():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family': ['Nimbus Sans', 'Liberation Sans', 'DejaVu Sans'], 'font.size': 7, 'axes.linewidth': 0.6,
                         'axes.edgecolor': '#999999', 'axes.labelcolor': MUTED, 'xtick.color': MUTED,
                         'ytick.color': MUTED, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
                         'axes.spines.top': False, 'axes.spines.right': False, 'savefig.dpi': 300})
    return plt


def figure(plt, key):
    w, h = SIZE[key]
    return plt.figure(figsize=(w / 100, h / 100))


def pass2(results_root, run):
    import anndata as ad
    return ad.read_h5ad(Path(results_root) / RUNS[run] / '03_integration/pass2.h5ad', backed='r')


# ---------------------------------------------------------------------------------------------------- a
def panel_a(args, plt, out):
    pos = load_positions(args.results_root, 'aki')
    a = pass2(args.results_root, 'aki')
    umap = pd.DataFrame(a.obsm['X_umap'], index=a.obs_names).reindex(pos.index).to_numpy()
    emb = pd.DataFrame(a.obsm['X_harmony'][:, :50], index=a.obs_names).reindex(pos.index).to_numpy()
    nodes = np.load(Path(args.results_root) / RUNS['aki'] / '04_pseudospace/curve_nodes.npz')
    order = np.argsort(nodes['pseudospace'])
    centres = nodes['nodes'][order]
    nearest = np.argmin(((emb[:, None, :] - centres[None]) ** 2).sum(-1), axis=1)
    node_xy = np.array([np.median(umap[nearest == k], 0) if (nearest == k).any() else [np.nan, np.nan]
                        for k in range(len(centres))])
    ok = np.isfinite(node_xy[:, 0])
    node_xy = node_xy[ok]
    rng = np.random.default_rng(0)
    draw = rng.permutation(len(pos))
    fig = figure(plt, 'a')
    for i, mode in enumerate(('t', 'group')):
        ax = fig.add_axes([0.01 + 0.5 * i, 0.04, 0.42, 0.92])
        if mode == 't':
            sc = ax.scatter(*umap[draw].T, c=pos.t.to_numpy()[draw], cmap='plasma', vmin=0, vmax=1, s=0.6, lw=0,
                            rasterized=True)
            cax = fig.add_axes([0.445, 0.25, 0.009, 0.5])
            cb = fig.colorbar(sc, cax=cax, ticks=[0, 0.5, 1])
            cb.outline.set_linewidth(0)
            cb.set_label('pseudostructure', color=MUTED)
        else:
            ax.scatter(*umap[draw].T, c=pos.group.map(COLOURS).to_numpy()[draw], s=0.6, lw=0, rasterized=True)
            for g in GROUP_ORDER:
                k = (pos.group == g).to_numpy()
                if k.any():
                    x, y = np.median(umap[k], 0)
                    ax.text(x, y, g, ha='center', va='center', fontsize=7, weight='bold', color=COLOURS[g],
                            bbox=dict(boxstyle='round,pad=0.25', fc='white', ec=COLOURS[g], lw=0.8), zorder=6)
        ax.plot(*node_xy.T, color='k', lw=1.1, zorder=3)
        ax.scatter(*node_xy.T, s=7, facecolor='white', edgecolor='k', lw=0.6, zorder=4)
        for j, lab in ((0, '0'), (-1, '1')):
            ax.scatter(*node_xy[j], marker='*', s=90, color='#D7263D', edgecolor='k', lw=0.4, zorder=5)
            ax.annotate(lab, node_xy[j], xytext=(6, -2), textcoords='offset points', fontsize=8, weight='bold')
        ax.set_axis_off()
    fig.savefig(out / 'f3_a.png', transparent=False, facecolor='white')
    plt.close(fig)
    return {'n_units': int(len(pos))}


# ---------------------------------------------------------------------------------------------------- b
def label_statistics(args):
    """Unit-level rho of t with the label rank per section, and a within-section permutation null."""
    rng = np.random.default_rng(0)
    rows = []
    for run in ('aki', 'hvm'):
        pos = load_positions(args.results_root, run)
        sections = [s for s, r in SECTIONS.items() if r == run]
        obs = section_rho(pos.t, pos, sections).set_index('sample')
        for s in sections:
            d = pos[(pos['sample'] == s) & pos['rank'].notna()]
            ranks = d['rank'].to_numpy(float)
            t = d.t.to_numpy()
            null = np.array([np.corrcoef(pd.Series(rng.permutation(t)).rank(), pd.Series(ranks).rank())[0, 1]
                             for _ in range(args.permutations)])
            rho = obs.loc[s, 'rho']
            rows.append({'sample': s, 'run': run, 'n': int(obs.loc[s, 'n']), 'rho': rho,
                         'null_max_abs': float(np.abs(null).max()),
                         'p_perm': (1 + int((null >= rho).sum())) / (1 + args.permutations)})
    return pd.DataFrame(rows)


def panel_b(args, plt, out):
    from scipy.stats import gaussian_kde

    stats = label_statistics(args)
    stats.to_csv(out / 'f3_b_rho.csv', index=False)
    mouse = load_positions(args.results_root, 'aki')
    human = load_positions(args.results_root, 'hvm')
    human = human[human['sample'].isin(['HUK1_COR1', 'HUK1_MED1'])]
    groups = [g for g in GROUP_ORDER if (mouse.group == g).any() or (human.group == g).any()]
    fig = figure(plt, 'b')
    ax = fig.add_axes([0.16, 0.30, 0.80, 0.68])
    grid = np.linspace(0, 1, 400)
    for i, g in enumerate(groups):
        base = len(groups) - 1 - i
        for d, filled in ((mouse, True), (human, False)):
            t = d.t[d.group == g].to_numpy()
            if len(t) < 20:
                continue
            y = gaussian_kde(t)(grid)
            y = 0.85 * y / y.max()
            if filled:
                ax.fill_between(grid, base, base + y, color=COLOURS[g], alpha=0.45, lw=0)
                ax.plot(grid, base + y, color=COLOURS[g], lw=0.9)
            else:
                ax.plot(grid, base + y, color=COLOURS[g], lw=0.9, ls=(0, (2, 1.2)))
        ax.text(-0.02, base + 0.15, g, ha='right', va='bottom', fontsize=7, color=INK)
        ax.axhline(base, color='#DDDDDD', lw=0.5, zorder=0)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, len(groups))
    ax.fill_between([0.70, 0.76], len(groups) - 0.55, len(groups) - 0.35, color='#999999', alpha=0.45, lw=0)
    ax.text(0.78, len(groups) - 0.45, 'mouse', va='center', fontsize=6, color=MUTED)
    ax.plot([0.70, 0.76], [len(groups) - 0.85] * 2, color='#999999', lw=0.9, ls=(0, (2, 1.2)))
    ax.text(0.78, len(groups) - 0.85, 'human', va='center', fontsize=6, color=MUTED)
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.set_xticks([0, 0.5, 1])
    ax.set_xticklabels(['0', '0.5', '1'])
    ax.set_xlabel('pseudostructure', labelpad=1)
    # rho per section
    rax = fig.add_axes([0.16, 0.02, 0.80, 0.14])
    rax.set_axis_off()
    rax.text(-0.02, 0.62, 'ρ', ha='right', va='center', fontsize=7, weight='bold', color=INK,
             transform=rax.transAxes)
    for i, r in enumerate(stats.itertuples()):
        x = (i + 0.5) / len(stats)
        rax.text(x, 0.62, f'{r.rho:.2f}', ha='center', va='center', fontsize=8, weight='bold',
                 color=SECTION_COLOURS[r.sample], transform=rax.transAxes)
        rax.text(x, 0.08, SECTION_NAMES[r.sample], ha='center', va='center', fontsize=5.5, color=MUTED,
                 transform=rax.transAxes)
    fig.savefig(out / 'f3_b.png', facecolor='white')
    plt.close(fig)
    return stats


# ---------------------------------------------------------------------------------------------------- c
def binned_markers(results_root, run, samples, bins=100, min_units=5):
    pos = load_positions(results_root, run)
    pos = pos[pos['sample'].isin(samples)]
    a = pass2(results_root, run)
    genes = [g for g, _ in MARKERS]
    lognorm = a.layers['lognorm']        # backed mode reads layers into memory; only X stays on disk
    cols = [a.var_names.get_loc(g) for g in genes]
    rows = a.obs_names.get_indexer(pos.index)
    sub = lognorm[rows][:, cols]
    expr = pd.DataFrame(sub.toarray() if hasattr(sub, 'toarray') else np.asarray(sub), index=pos.index, columns=genes)
    b = np.minimum((pos.t.to_numpy() * bins).astype(int), bins - 1)
    count = np.bincount(b, minlength=bins)
    mean = np.vstack([np.bincount(b, weights=expr[g].to_numpy(), minlength=bins) / np.maximum(count, 1) for g in genes])
    mean[:, count < min_units] = np.nan
    z = (mean - np.nanmean(mean, 1, keepdims=True)) / np.nanstd(mean, 1, keepdims=True)
    label = pd.Series(pos.group.to_numpy()).groupby(b).agg(lambda s: s.value_counts().index[0]).reindex(range(bins))
    label[count < min_units] = None
    return z, label


def panel_c(args, plt, out):
    from matplotlib.colors import ListedColormap, LinearSegmentedColormap

    blocks = [('Control', 'aki', ['Ctrl1A2', 'Ctrl1A4'], '#2F6DB5', False),
              ('AKI', 'aki', ['IR2A2', 'IR2A4'], '#D9692B', False),
              ('Human', 'hvm', ['HUK1_COR1', 'HUK1_MED1'], '#6A4BB0', True)]
    cmap = LinearSegmentedColormap.from_list('bwr', ['#3B4CC0', '#C9D3EE', '#FFFFFF', '#F4B6B0', '#D7263D'])
    cmap.set_bad('#F2F2F2')
    fig = figure(plt, 'c')
    top, gap, n = 0.985, 0.045, len(MARKERS)
    unit = (top - 0.06 - 2 * gap) / (3 * (n + 1))
    for k, (name, run, samples, colour, upper) in enumerate(blocks):
        z, label = binned_markers(args.results_root, run, samples)
        y1 = top - k * ((n + 1) * unit + gap)
        lab_ax = fig.add_axes([0.20, y1 - unit, 0.72, unit * 0.9])
        codes = label.map({g: i for i, g in enumerate(GROUP_ORDER)}).to_numpy(float)
        lab_ax.imshow(codes[None], aspect='auto', interpolation='nearest',
                      cmap=ListedColormap([COLOURS[g] for g in GROUP_ORDER]), vmin=-0.5, vmax=len(GROUP_ORDER) - 0.5)
        lab_ax.set_axis_off()
        for g in GROUP_ORDER:
            idx = np.flatnonzero(label.to_numpy() == g)
            if len(idx) > 3:
                lab_ax.text(np.median(idx), 0, g, ha='center', va='center', fontsize=5.5, color='white', weight='bold')
        hax = fig.add_axes([0.20, y1 - (n + 1) * unit, 0.72, n * unit])
        im = hax.imshow(np.clip(z, -3, 3), aspect='auto', interpolation='nearest', cmap=cmap, vmin=-3, vmax=3,
                        extent=(0, 1, n, 0))
        hax.set_yticks(np.arange(n) + 0.5)
        hax.set_yticklabels([(g.upper() if upper else g) + f' ({s})' for g, s in MARKERS], fontsize=6, style='italic',
                            color=INK)
        hax.tick_params(axis='y', length=0)
        for sp in hax.spines.values():
            sp.set_visible(False)
        if k in (1, 2):
            hax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
        if k == 2:
            hax.set_xlabel('pseudostructure', labelpad=1)
        else:
            hax.set_xticks([])
        fig.text(0.015, y1 - (n + 1) * unit / 2, name, rotation=90, ha='left', va='center', fontsize=9, weight='bold',
                 color=colour)
    cax = fig.add_axes([0.945, 0.45, 0.012, 0.2])
    cb = fig.colorbar(im, cax=cax, ticks=[-3, 0, 3])
    cb.outline.set_linewidth(0)
    fig.savefig(out / 'f3_c.png', facecolor='white')
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------- d
def panel_d(args, plt, out):
    root = Path(args.results_root) / 'paper_figures/fig3/robustness'
    fits = {run: pd.read_csv(root / run / 'robustness.csv') for run in ('aki', 'hvm') if (root / run / 'robustness.csv').exists()}
    if not fits:
        print('d: no robustness results yet')
        return None
    from matplotlib.lines import Line2D
    fig = figure(plt, 'd')
    ax = fig.add_axes([0.13, 0.26, 0.60, 0.68])
    arms = {'tubules': '#2F6DB5', 'reads': '#D9692B', 'genes': '#5BA55B'}
    dashes = {'aki': '-', 'hvm': (0, (2, 1.2))}           # median over seeds; every rerun is a dot (mouse filled)
    for run, r in fits.items():
        for k, (cond, c) in enumerate(arms.items()):
            d = r[r.condition == cond]
            g = d.groupby('fraction').rho_vs_reference.median()
            ax.plot(100 * g.index, g, color=c, ls=dashes[run], lw=1.0)
            jitter = (k - 1) * 0.9 + (0.45 if run == 'hvm' else -0.45)
            ax.scatter(100 * d.fraction + jitter, d.rho_vs_reference, s=3.5, lw=0.5, zorder=3,
                       color=c if run == 'aki' else 'white', edgecolor=c)
    lo = min(float(r.rho_vs_reference.min()) for r in fits.values())
    ax.set_xlim(0, 103)
    ax.set_ylim(np.floor(lo * 10) / 10 - 0.02, 1.02)
    ax.set_xticks([10, 30, 50, 70, 100])
    ax.set_xlabel('fraction kept (%)', labelpad=1)
    ax.set_ylabel('ρ with full-data axis')
    keys = [Line2D([], [], color=c, lw=1.1) for c in arms.values()]
    keys += [Line2D([], [], color='#999999', lw=1.1), Line2D([], [], color='#999999', lw=1.1, ls=dashes['hvm'])]
    ax.legend(keys, [*arms, 'mouse', 'human + mouse'], loc='lower right', fontsize=6, frameon=False,
              handlelength=1.6, labelspacing=0.25, borderaxespad=0.2)
    hx = fig.add_axes([0.80, 0.26, 0.17, 0.68])           # five complementary gene splits: rho between the halves
    for i, (run, filled) in enumerate((('aki', True), ('hvm', False))):
        p = root / run / 'gene_halves.csv'
        if not p.exists():
            continue
        h = pd.read_csv(p)['rho_A_vs_B']
        hx.scatter(i + np.linspace(-0.12, 0.12, len(h)), h, s=9, color=arms['genes'] if filled else 'white',
                   edgecolor=arms['genes'], lw=0.8, zorder=3)
    hx.set_xlim(-0.6, 1.6)
    hx.set_ylim(ax.get_ylim())
    hx.set_xticks([0, 1])
    hx.set_xticklabels(['mouse', 'human\n+ mouse'], fontsize=6)
    hx.set_yticklabels([])
    hx.set_title('gene half A vs B', fontsize=6, color=MUTED, pad=2)
    fig.savefig(out / 'f3_d.png', facecolor='white')
    plt.close(fig)
    return fits


# ---------------------------------------------------------------------------------------------------- e
def panel_e(args, plt, out):
    path = Path(args.results_root) / 'paper_figures/fig3/baselines/baseline_rho.csv'
    if not path.exists():
        print('e: no baseline results yet')
        return None
    r = pd.read_csv(path)
    r = r[r.apply(lambda x: SECTIONS.get(x['sample']) == x['run'], axis=1)]
    names = {'t': 'pseudostructure', 'pseudostructure': 'pseudostructure', 'dpt': 'DPT', 'slingshot': 'Slingshot',
             'scfates': 'scFates', 'gaston': 'GASTON'}
    methods = [m for m in ('t', 'pseudostructure', 'dpt', 'slingshot', 'scfates', 'gaston') if m in set(r.method)]
    fig = figure(plt, 'e')
    ax = fig.add_axes([0.25, 0.24, 0.70, 0.74])
    for i, m in enumerate(methods):
        d = r[r.method == m]
        y = len(methods) - 1 - i
        ax.barh(y, d.rho.mean(), height=0.6, color='#DCE7F3' if i == 0 else '#EEEEEE', edgecolor='#999999', lw=0.5)
        jitter = np.linspace(-0.18, 0.18, len(d))
        ax.scatter(d.rho, y + jitter, s=9, c=d['sample'].map(SECTION_COLOURS), lw=0, zorder=3)
        ax.text(-0.02, y, names[m], ha='right', va='center', fontsize=7, weight='bold' if i == 0 else 'normal',
                color=INK, transform=ax.get_yaxis_transform())
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    lo = min(0.0, float(r.rho.min()) - 0.05)
    ax.set_xlim(lo, 1)
    ax.axvline(0, color='#BBBBBB', lw=0.5)
    ax.set_xlabel('Spearman ρ with segment order', labelpad=1)
    fig.savefig(out / 'f3_e.png', facecolor='white')
    plt.close(fig)
    return r


# ---------------------------------------------------------------------------------------------------- f
def panel_f(args, plt, out):
    from scipy.stats import spearmanr

    d = pd.read_parquet(Path(args.results_root) / 'paper_figures/fig3/depth/depth.parquet')
    pt = d[d.group.isin(['PTS1', 'PTS2', 'PTS3'])]
    rows = []
    for s, g in pt.groupby('sample'):
        row = {'sample': s, 'n_pt': len(g), 'rho_t_depth': spearmanr(g.t, g.depth_um).statistic,
               'rho_label_rank_depth': spearmanr(g['rank'], g.depth_um).statistic}
        for lab, x in g.groupby('group'):
            row[f'rho_within_{lab}'] = spearmanr(x.t, x.depth_um).statistic
            row[f'median_depth_{lab}'] = x.depth_um.median()
        rows.append(row)
    stats = pd.DataFrame(rows)
    stats.to_csv(out / 'f3_f_depth_rho.csv', index=False)
    fig = figure(plt, 'f')
    ax = fig.add_axes([0.14, 0.26, 0.52, 0.70])
    ctrl = pt[pt['sample'].isin(['Ctrl1A2', 'Ctrl1A4'])]
    ax.scatter(ctrl.depth_um, ctrl.t, s=0.8, c=ctrl.group.map(COLOURS), alpha=0.35, lw=0, rasterized=True)
    edges = np.arange(0, 2001, 100)
    for s, g in pt.groupby('sample'):
        b = np.digitize(g.depth_um, edges)
        med = g.groupby(b).t.median()
        n = g.groupby(b).size()
        med = med[n >= 20]
        ax.plot((edges[med.index - 1] + 50), med, color=SECTION_COLOURS[s], lw=1.2)
    ax.set_xlim(0, 2000)
    ax.set_ylim(0, pt.t.quantile(0.995) + 0.03)
    ax.set_xlabel('distance from capsule (µm)', labelpad=1)
    ax.set_ylabel('pseudostructure')
    for i, r in enumerate(stats.sort_values('sample').itertuples()):
        fig.text(0.72, 0.78 - i * 0.13, f'{SECTION_NAMES[r.sample]}', fontsize=6.5, color=SECTION_COLOURS[r.sample])
        fig.text(0.90, 0.78 - i * 0.13, f'{r.rho_t_depth:.2f}', fontsize=7, weight='bold', ha='right',
                 color=SECTION_COLOURS[r.sample])
    fig.text(0.90, 0.90, 'ρ (PT)', fontsize=6.5, ha='right', color=MUTED)
    fig.savefig(out / 'f3_f.png', facecolor='white')
    plt.close(fig)
    return stats


def main(argv=None):
    args = parse_args(argv)
    plt = style()
    out = Path(args.results_root) / args.out / 'panels'
    out.mkdir(parents=True, exist_ok=True)
    draw = {'a': panel_a, 'b': panel_b, 'c': panel_c, 'd': panel_d, 'e': panel_e, 'f': panel_f}
    for k in args.panels:
        res = draw[k](args, plt, out)
        if isinstance(res, pd.DataFrame):
            print(f'--- {k}\n{res.to_string()}')
        else:
            print(f'--- {k} {res or ""}')


if __name__ == '__main__':
    main()
