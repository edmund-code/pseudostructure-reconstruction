"""Fig. 1 (overview) panels drawn from data, sized for the layout of the figure composer (make_fig1).

 slides   all six sections: H&E thumbnails and Visium HD total counts per bin, in the same pixel frame
 seg      one control-mouse cortex window: H&E, then the segmentation's instances over it
 agg      the same window: 2 um bins coloured by total counts inside each mask, then a tubule x gene matrix
 harm     Harmony, before and after: control vs AKI mouse (AKI fit), mouse vs human (human-mouse fit)
 graph    the elastic curve on each fit's integrated UMAP, units coloured by pseudostructure
 plane    a larger window of the same section with every segmented cross-section outlined, the numbered
          cross-sections the composer carries onto the nephron, and the segment spans of the bar (JSON)
 ds       downstream examples: control PT gradients (Fig. 4 screen) and the SGLT2 -> SGLT1 handoff in mouse
          and human (Fig. 5 screen)

Reads private data (the package runs under --results-root, Visium HD bins, segmentations under --data-root, H&E
under --he-root) and writes PNGs and JSON to <results-root>/paper_figures/fig1/panels/.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = next(p for p in Path(__file__).resolve().parents if (p / 'pseudospace').is_dir())
sys.path[:0] = [str(REPO), str(Path(__file__).resolve().parent)]
from fig3_common import RUNS, load_positions  # noqa: E402
from fig3_panels import COLOURS, MUTED, style  # noqa: E402

# sample -> (package config that lists its bins, H&E file stem, card label)
SAMPLES = {'Ctrl1A2': ('aki', 'Ctrl_1A2', 'Control'), 'Ctrl1A4': ('aki', 'Ctrl_1A4', 'Control'),
           'IR2A2': ('aki', 'IR_2A2', 'AKI'), 'IR2A4': ('aki', 'IR_2A4', 'AKI'),
           'HUK1_COR1': ('human_vs_mouse', 'HUK1_COR_1', 'Human'), 'HUK1_MED1': ('human_vs_mouse', 'HUK1_MED_1', 'Human')}
SEGMENTATION = {'Ctrl1A2': 'segmentations/v4fix/Ctrl1A2_v4fix.geojson'}
CTRL, AKI, MOUSE, HUMAN = '#2F6DB5', '#D9692B', '#2F6DB5', '#C51B7D'   # human magenta: purple blurred into blue in the UMAPs
GROUPS = ['PTS1', 'PTS2', 'PTS3', 'AL', 'DCT', 'CNT', 'CD']
MARKERS = ['Slc5a2', 'Slc13a3', 'Slc7a13', 'Slc12a1', 'Slc12a3', 'Slc8a1', 'Aqp2']
UM_PER_PX = 0.4408725121969895          # Ctrl1A2 full-resolution frame (scalefactors_json.json)
SECTION = 'Ctrl1A2'                     # the section shown in seg, agg and plane
# windows in full-resolution pixels (x0, y0, side): plane = 450 um, six or more balanced groups, then most units;
# seg/agg = a 200 um window inside it; None = choose with window_search() (deterministic)
PLANE_WINDOW = None
SEG_WINDOW = None
COUNT_CMAP = 'inferno'


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--results-root', default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT', str(REPO / 'results')))
    p.add_argument('--data-root', default=os.environ.get('PSEUDOSPACE_DATA_ROOT', str(REPO / 'data')))
    p.add_argument('--visium-root', default=os.environ.get('PSEUDOSPACE_VISIUM_ROOT'),
                   help='folder holding the mouse Visium HD bins (data/README.md)')
    p.add_argument('--he-root', default=os.environ.get('PSEUDOSPACE_HE_ROOT'),
                   help='folder holding the six full-resolution H&E TIFFs')
    p.add_argument('--out', default='paper_figures/fig1', help='relative to the results root')
    p.add_argument('--panels', default='slides,seg,agg,harm,graph,plane,ds')
    return p.parse_args(argv)


# ------------------------------------------------------------------------------------------------- inputs
def bins_dir(args, sample) -> Path:
    config = SAMPLES[sample][0]
    sheet = pd.read_csv(REPO / 'analysis/package_configs' / config / 'samples.csv', comment='#').set_index('sample')
    if args.visium_root:
        os.environ['PSEUDOSPACE_VISIUM_ROOT'] = args.visium_root
    path = Path(os.path.expandvars(sheet.loc[sample, 'bins_dir']))
    return path if path.is_absolute() else Path(args.data_root) / path


def bin_totals(args, sample, cache: Path) -> pd.DataFrame:
    """Total counts of every in-tissue 2 um bin with its full-resolution pixel position (cached)."""
    import h5py

    f = cache / f'{sample}_bin_totals.parquet'
    if f.exists():
        return pd.read_parquet(f)
    d = bins_dir(args, sample)
    with h5py.File(d / 'filtered_feature_bc_matrix.h5', 'r') as h:
        m = h['matrix']
        indptr = m['indptr'][:]
        cs = np.concatenate([[0], np.cumsum(m['data'][:], dtype=np.int64)])
        barcodes = m['barcodes'][:].astype(str)
    pos = pd.read_parquet(d / 'spatial/tissue_positions.parquet').set_index('barcode')
    capture = pos[['pxl_col_in_fullres', 'pxl_row_in_fullres']].agg(['min', 'max'])
    pos = pos.loc[barcodes]
    out = pd.DataFrame({'barcode': barcodes, 'x': pos.pxl_col_in_fullres.to_numpy(np.float32),
                        'y': pos.pxl_row_in_fullres.to_numpy(np.float32),
                        'total': (cs[indptr[1:]] - cs[indptr[:-1]]).astype(np.int32)})
    out.attrs['capture'] = capture.to_dict()
    cache.mkdir(parents=True, exist_ok=True)
    out.to_parquet(f)
    (cache / f'{sample}_capture.json').write_text(json.dumps(
        {'x0': float(capture.loc['min', 'pxl_col_in_fullres']), 'x1': float(capture.loc['max', 'pxl_col_in_fullres']),
         'y0': float(capture.loc['min', 'pxl_row_in_fullres']), 'y1': float(capture.loc['max', 'pxl_row_in_fullres'])}))
    return out


def he_path(args, sample) -> Path:
    return Path(args.he_root) / f'{SAMPLES[sample][1]}.tif'


def he_size(path) -> tuple[int, int]:
    w = int(subprocess.check_output(['vipsheader', '-f', 'width', str(path)]))
    h = int(subprocess.check_output(['vipsheader', '-f', 'height', str(path)]))
    return w, h


def he_crop(args, sample, x0, y0, side, tmp: Path) -> np.ndarray:
    from PIL import Image

    tmp.mkdir(parents=True, exist_ok=True)
    f = tmp / f'{sample}_{int(x0)}_{int(y0)}_{int(side)}.png'
    if not f.exists():
        subprocess.run(['vips', 'crop', str(he_path(args, sample)), str(f), str(int(x0)), str(int(y0)),
                        str(int(side)), str(int(side))], check=True)
    return np.asarray(Image.open(f).convert('RGB'))


def section_obs(args, run, sample) -> pd.DataFrame:
    """Units of one section on a fit's curve: centroid, feature index, group and pseudostructure t."""
    import h5py

    with h5py.File(Path(args.results_root) / RUNS[run] / '03_integration/pass2.h5ad', 'r') as h:
        o = h['obs']
        idx = o['_index'][:].astype(str)
        cols = {k: o[k][:] for k in ('x_centroid', 'y_centroid', 'feature_index')}
    obs = pd.DataFrame(cols, index=idx)
    pos = load_positions(args.results_root, run)
    d = pos.join(obs, how='inner')
    return d[d['sample'] == sample]


def polygons(args, sample, x0, y0, side):
    """Every segmented polygon of the section that intersects a window: [(feature_index, Nx2 array)]."""
    import shapely
    from pseudospace_reconstruction.io.segmentation import load_geometries

    _, geoms = load_geometries(Path(args.data_root) / SEGMENTATION[sample])
    box = shapely.box(x0, y0, x0 + side, y0 + side)
    out = []
    for i, g in enumerate(geoms):
        if g is None or not g.intersects(box):
            continue
        g = g if g.geom_type == 'Polygon' else max(g.geoms, key=lambda p: p.area)
        out.append((i, np.asarray(g.exterior.coords)))
    return out


def window_search(d: pd.DataFrame, side_um: float, inside=None, step=60):
    """Window with at least six groups, then the best balance (up to 5 units per group), then the most units.

    Mouse labels cluster in space (CD and CNT are rare in the labyrinth), so requiring all seven groups lands on the
    hilum; six balanced groups give a typical cortical view."""
    w = side_um / UM_PER_PX
    x, y = d.x_centroid.to_numpy(), d.y_centroid.to_numpy()
    code = pd.Categorical(d.group, categories=GROUPS).codes
    lo_x, lo_y = (inside[0], inside[1]) if inside else (x.min(), y.min())
    hi_x, hi_y = (inside[0] + inside[2] - w, inside[1] + inside[2] - w) if inside else (x.max() - w, y.max() - w)
    best = None
    for wx in np.arange(lo_x, hi_x + 1, step):
        kx = (x >= wx) & (x <= wx + w)
        for wy in np.arange(lo_y, hi_y + 1, step):
            k = kx & (y >= wy) & (y <= wy + w)
            n = np.bincount(code[k & (code >= 0)], minlength=len(GROUPS))
            key = (min(int((n > 0).sum()), 6), int(np.minimum(n, 5).sum()), int(k.sum()))
            if best is None or key > best[0]:
                best = (key, (float(wx), float(wy), float(w)))
    return best[1]


def save(fig, path, dark=False):
    fig.savefig(path, facecolor='black' if dark else 'white', dpi=300)


# ------------------------------------------------------------------------------------------------- slides
def panel_slides(args, plt, out, cache):
    """Six H&E thumbnails and six total-count maps in the same frame (one PNG each), plus the capture boxes."""
    from matplotlib import colormaps
    from PIL import Image
    from scipy.ndimage import gaussian_filter

    thumbs = cache / 'he'
    thumbs.mkdir(parents=True, exist_ok=True)
    meta = {}
    cmap = colormaps[COUNT_CMAP]
    for s in SAMPLES:
        src = he_path(args, s)
        full_w, full_h = he_size(src)
        f = thumbs / f'{s}_thumb.png'
        if not f.exists():
            subprocess.run(['vips', 'thumbnail', str(src), str(f), '900'], check=True)
        he = np.asarray(Image.open(f).convert('RGB')).copy()
        h, w = he.shape[:2]
        sc = w / full_w
        trim = max(2, int(0.006 * h))                      # scanner edge rows/columns
        he[:trim], he[-trim:], he[:, :trim], he[:, -trim:] = 247, 247, 247, 247
        Image.fromarray(he).save(out / f'f1_he_{s}.png')
        b = bin_totals(args, s, cache / 'bins')
        px = np.clip((b.x.to_numpy() * sc).astype(int), 0, w - 1)
        py = np.clip((b.y.to_numpy() * sc).astype(int), 0, h - 1)
        flat = py * w + px
        tot = np.bincount(flat, weights=b.total.to_numpy(float), minlength=h * w).reshape(h, w)
        n = np.bincount(flat, minlength=h * w).reshape(h, w).astype(float)
        # light smoothing (about 25 um) hides the row/column banding of the bin grid at thumbnail scale
        tot_s, n_s = gaussian_filter(tot, 0.8), gaussian_filter(n, 0.8)
        v = np.log1p(np.divide(tot_s, n_s, out=np.zeros_like(tot_s), where=n_s > 0))
        lo, hi = np.quantile(v[n > 0], [0.02, 0.995])
        rgb = (cmap(np.clip((v - lo) / (hi - lo), 0, 1))[..., :3] * 255).astype(np.uint8)
        rgb[n == 0] = (8, 8, 12)
        Image.fromarray(rgb).save(out / f'f1_counts_{s}.png')
        cap = json.loads((cache / 'bins' / f'{s}_capture.json').read_text())
        meta[s] = {'label': SAMPLES[s][2], 'median_counts_per_bin': float(np.median(b.total)), 'n_bins': int(len(b)),
                   # capture area as fractions of the thumbnail
                   'capture': {'x0': cap['x0'] * sc / w, 'x1': cap['x1'] * sc / w,
                               'y0': cap['y0'] * sc / h, 'y1': cap['y1'] * sc / h}}
    (out / 'f1_slides.json').write_text(json.dumps(meta, indent=1))
    # colour bar for the count maps
    fig = plt.figure(figsize=(0.9, 0.12))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(np.linspace(0, 1, 256)[None], aspect='auto', cmap=cmap)
    ax.set_axis_off()
    save(fig, out / 'f1_counts_bar.png')
    plt.close(fig)
    return meta


# ------------------------------------------------------------------------------------------------- windows
def windows(args):
    global PLANE_WINDOW, SEG_WINDOW
    d = section_obs(args, 'hvm', SECTION)
    if PLANE_WINDOW is None:
        PLANE_WINDOW = window_search(d, 450)
    if SEG_WINDOW is None:
        SEG_WINDOW = window_search(d, 200, inside=PLANE_WINDOW, step=40)
    return d


def instance_colours(n, seed=3):
    rng = np.random.default_rng(seed)
    from matplotlib.colors import hsv_to_rgb
    return hsv_to_rgb(np.c_[rng.uniform(0, 1, n), rng.uniform(0.45, 0.8, n), rng.uniform(0.85, 1, n)])


def square_axes(plt, px=420):
    fig = plt.figure(figsize=(px / 300, px / 300))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    return fig, ax


def panel_seg(args, plt, out, cache):
    from matplotlib.patches import Polygon

    x0, y0, side = SEG_WINDOW
    img = he_crop(args, SECTION, x0, y0, side, cache / 'crops')
    polys = polygons(args, SECTION, x0, y0, side)
    for k, with_masks in enumerate((False, True)):
        fig, ax = square_axes(plt)
        ax.imshow(img, extent=(x0, x0 + side, y0 + side, y0), interpolation='lanczos')
        if with_masks:
            cols = instance_colours(len(polys))
            for (i, xy), c in zip(polys, cols):
                ax.add_patch(Polygon(xy, closed=True, fc=(*c, 0.5), ec=c * 0.6, lw=0.7))
        ax.set_xlim(x0, x0 + side)
        ax.set_ylim(y0 + side, y0)
        save(fig, out / f"f1_seg_{'masks' if with_masks else 'he'}.png")
        plt.close(fig)
    return {'window_px': SEG_WINDOW, 'window_um': side * UM_PER_PX, 'n_polygons': len(polys)}


def panel_agg(args, plt, out, cache, d):
    import h5py
    from matplotlib.patches import Polygon

    x0, y0, side = SEG_WINDOW
    b = bin_totals(args, SECTION, cache / 'bins')
    k = b.x.between(x0, x0 + side) & b.y.between(y0, y0 + side)
    b = b[k]
    assign = pd.read_parquet(Path(args.results_root) / RUNS['hvm'] / f'01_units/{SECTION}_bin_to_unit.parquet',
                             columns=['bin_id', 'polygon_id'])
    in_mask = b.barcode.isin(assign.bin_id[assign.polygon_id.notna()]).to_numpy()
    v = np.log1p(b.total.to_numpy(float))
    lo, hi = np.quantile(v, [0.02, 0.99])
    fig, ax = square_axes(plt)
    fig.set_facecolor('black')
    side_pts = 2 / UM_PER_PX / side * (420 / 300 * 72)          # one 2 um bin in points
    for mask, alpha in ((~in_mask, 0.28), (in_mask, 1.0)):
        ax.scatter(b.x[mask], b.y[mask], c=v[mask], cmap=COUNT_CMAP, vmin=lo, vmax=hi, marker='s',
                   s=side_pts ** 2, lw=0, alpha=alpha, rasterized=True)
    for i, xy in polygons(args, SECTION, x0, y0, side):
        ax.add_patch(Polygon(xy, closed=True, fc='none', ec='white', lw=0.6, alpha=0.9))
    ax.set_xlim(x0, x0 + side)
    ax.set_ylim(y0 + side, y0)
    ax.set_facecolor('black')
    save(fig, out / 'f1_agg_bins.png', dark=True)
    plt.close(fig)

    # tubule x gene matrix: the window's units on the curve, marker genes, log-normalized, scaled per gene
    w = d[d.x_centroid.between(x0, x0 + side) & d.y_centroid.between(y0, y0 + side)]
    w = w.sort_values(['y_centroid', 'x_centroid']).head(16)
    with h5py.File(Path(args.results_root) / RUNS['hvm'] / '03_integration/pass2.h5ad', 'r') as h:
        names = h['obs']['_index'][:].astype(str)
        genes = h['var']['_index'][:].astype(str)
        rows = np.flatnonzero(np.isin(names, w.index))
        order = pd.Series(np.arange(len(rows)), index=names[rows]).reindex(w.index).to_numpy()
        lay = h['layers']['lognorm']
        cols = [int(np.flatnonzero(genes == g)[0]) for g in MARKERS]
        mat = np.zeros((len(rows), len(cols)))
        for r, i in enumerate(rows):
            a, z = lay['indptr'][i], lay['indptr'][i + 1]
            ind, val = lay['indices'][a:z], lay['data'][a:z]
            for c, j in enumerate(cols):
                hit = np.flatnonzero(ind == j)
                mat[r, c] = val[hit[0]] if hit.size else 0
    mat = mat[order]
    mat = mat / np.maximum(mat.max(0, keepdims=True), 1e-9)
    fig = plt.figure(figsize=(600 / 300, 420 / 300))
    ax = fig.add_axes([0.0, 0.0, 1, 1])
    ax.imshow(mat, cmap='Blues', vmin=0, vmax=1, aspect='auto', interpolation='nearest')
    ax.set_xticks(np.arange(len(cols)) - 0.5, minor=True)
    ax.set_yticks(np.arange(len(rows)) - 0.5, minor=True)
    ax.grid(which='minor', color='white', lw=0.6)
    ax.tick_params(which='both', length=0)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    save(fig, out / 'f1_agg_matrix.png')
    plt.close(fig)
    return {'n_bins': int(len(b)), 'bins_in_masks': int(in_mask.sum()), 'matrix_units': int(len(rows)),
            'matrix_genes': MARKERS}


# ------------------------------------------------------------------------------------------------- harmony
def panel_harm(args, plt, out):
    import h5py

    stats = {}
    for run, groups in (('aki', {'Control': CTRL, 'AKI': AKI}), ('hvm', {'mouse': MOUSE, 'human': HUMAN})):
        with h5py.File(Path(args.results_root) / RUNS[run] / '03_integration/pass1_labelled.h5ad', 'r') as h:
            sample = h['obs']['sample']
            if isinstance(sample, h5py.Group):
                sample = sample['categories'][:].astype(str)[sample['codes'][:]]
            else:
                sample = sample[:].astype(str)
            emb = {'pre': h['obsm']['X_umap_pre_harmony'][:], 'post': h['obsm']['X_umap'][:]}
        if run == 'aki':
            key = np.where(pd.Series(sample).str.startswith('Ctrl'), 'Control', 'AKI')
        else:
            key = np.where(pd.Series(sample).str.startswith('HUK'), 'human', 'mouse')
        rng = np.random.default_rng(0)
        draw = rng.permutation(len(key))
        for when, xy in emb.items():
            fig, ax = square_axes(plt, 330)
            ax.scatter(*xy[draw].T, c=[groups[k] for k in key[draw]], s=0.25, lw=0, alpha=0.6, rasterized=True)
            ax.set_aspect('equal', 'datalim')
            save(fig, out / f'f1_harm_{run}_{when}.png')
            plt.close(fig)
        stats[run] = {g: int((key == g).sum()) for g in groups}
    return stats


# ------------------------------------------------------------------------------------------------- curve
def panel_graph(args, plt, out):
    import h5py

    stats = {}
    for run in ('aki', 'hvm'):
        pos = load_positions(args.results_root, run)
        with h5py.File(Path(args.results_root) / RUNS[run] / '03_integration/pass2.h5ad', 'r') as h:
            names = h['obs']['_index'][:].astype(str)
            umap = pd.DataFrame(h['obsm']['X_umap'][:], index=names).reindex(pos.index).to_numpy()
            emb = pd.DataFrame(h['obsm']['X_harmony'][:, :50], index=names).reindex(pos.index).to_numpy()
        nodes = np.load(Path(args.results_root) / RUNS[run] / '04_pseudospace/curve_nodes.npz')
        centres = nodes['nodes'][np.argsort(nodes['pseudospace'])]
        nearest = np.argmin(((emb[:, None, :] - centres[None]) ** 2).sum(-1), axis=1)
        node_xy = np.array([np.median(umap[nearest == k], 0) if (nearest == k).any() else [np.nan, np.nan]
                            for k in range(len(centres))])
        node_xy = node_xy[np.isfinite(node_xy[:, 0])]
        draw = np.random.default_rng(0).permutation(len(pos))
        fig, ax = square_axes(plt)
        ax.scatter(*umap[draw].T, c=pos.t.to_numpy()[draw], cmap='plasma', vmin=0, vmax=1, s=0.35, lw=0,
                   alpha=0.7, rasterized=True)
        ax.plot(*node_xy.T, color='k', lw=1.0, zorder=3)
        ax.scatter(*node_xy.T, s=5, facecolor='white', edgecolor='k', lw=0.5, zorder=4)
        for j in (0, -1):
            ax.scatter(*node_xy[j], marker='o', s=26, color='#D7263D', edgecolor='k', lw=0.4, zorder=5)
        ax.set_aspect('equal', 'datalim')
        save(fig, out / f'f1_graph_{run}.png')
        plt.close(fig)
        stats[run] = {'n_units': int(len(pos)), 'n_nodes': int(len(centres))}
    fig = plt.figure(figsize=(0.9, 0.12))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(np.linspace(0, 1, 256)[None], aspect='auto', cmap='plasma')
    ax.set_axis_off()
    save(fig, out / 'f1_t_bar.png')
    plt.close(fig)
    return stats


# ------------------------------------------------------------------------------------------------- output plane
def group_spans(pos: pd.DataFrame) -> dict:
    """Bar spans: between consecutive groups, the threshold on t that misplaces the fewest units of the two."""
    present = [g for g in GROUPS if (pos.group == g).any()]
    cuts = []
    for a, b in zip(present, present[1:]):
        ta, tb = np.sort(pos.t[pos.group == a].to_numpy()), np.sort(pos.t[pos.group == b].to_numpy())
        grid = np.linspace(np.median(ta), np.median(tb), 400)
        wrong = [(ta > c).sum() + (tb <= c).sum() for c in grid]
        cuts.append(float(grid[int(np.argmin(wrong))]))
    edges = [0.0] + cuts + [1.0]
    return {g: [edges[i], edges[i + 1]] for i, g in enumerate(present)}


def panel_plane(args, plt, out, cache, d):
    from matplotlib.patches import Polygon

    x0, y0, side = PLANE_WINDOW
    img = he_crop(args, SECTION, x0, y0, side, cache / 'crops')
    polys = polygons(args, SECTION, x0, y0, side)
    w = d[d.x_centroid.between(x0, x0 + side) & d.y_centroid.between(y0, y0 + side)].copy()
    pos = load_positions(args.results_root, 'hvm')
    spans = group_spans(pos[pos.group.isin(GROUPS)])
    med = pos.groupby('group').t.median()
    # one numbered cross-section per group: near the group's median t, away from the edges and each other
    chosen = []
    margin = 0.08 * side
    for g in GROUPS:
        c = w[(w.group == g) & w.x_centroid.between(x0 + margin, x0 + side - margin)
              & w.y_centroid.between(y0 + margin, y0 + side - margin)].copy()
        c['score'] = (c.t - med[g]).abs()
        for name, r in c.sort_values('score').iterrows():
            if all(np.hypot(r.x_centroid - q.x_centroid, r.y_centroid - q.y_centroid) > 0.13 * side for q in chosen):
                chosen.append(r.rename(name))
                break
    sel = pd.DataFrame(chosen)
    # number them by fitted position, so the bar reads 1..n while the section shows them scattered
    sel = sel.sort_values('t')
    sel['number'] = np.arange(1, len(sel) + 1)
    fi = set(sel.feature_index.astype(int))
    fig, ax = square_axes(plt, 900)
    ax.imshow(img, extent=(x0, x0 + side, y0 + side, y0), interpolation='lanczos')
    for i, xy in polys:
        hit = i in fi
        ax.add_patch(Polygon(xy, closed=True, fc='none', ec='#111111' if hit else (0.40, 0.12, 0.28, 0.85),
                             lw=2.6 if hit else 1.1, zorder=3 if hit else 2))
    ax.set_xlim(x0, x0 + side)
    ax.set_ylim(y0 + side, y0)
    save(fig, out / 'f1_plane.png')
    plt.close(fig)
    numbered = [{'number': int(r.number), 'group': r.group, 't': float(r.t), 'unit': str(n),
                 'fx': float((r.x_centroid - x0) / side), 'fy': float((r.y_centroid - y0) / side)}
                for n, r in sel.iterrows()]
    meta = {'window_px': PLANE_WINDOW, 'window_um': side * UM_PER_PX, 'n_polygons': len(polys),
            'n_units_on_curve': int(len(w)), 'spans': spans, 'numbered': numbered,
            'group_medians': {g: float(med[g]) for g in GROUPS if g in med}}
    (out / 'f1_plane.json').write_text(json.dumps(meta, indent=1))
    return meta


# ------------------------------------------------------------------------------------------------- downstream
def bands(ax, b12, b23, top):
    for (lo, hi), g in zip(((0, b12), (b12, b23), (b23, 1)), ('PTS1', 'PTS2', 'PTS3')):
        ax.axvspan(lo, hi, color=COLOURS[g], alpha=0.10, lw=0, zorder=0)
        ax.text((lo + hi) / 2, top, g.replace('PT', ''), ha='center', va='bottom', fontsize=7, weight='bold',
                color=COLOURS[g], transform=ax.get_xaxis_transform())


def ds_axes(plt):
    fig = plt.figure(figsize=(780 / 300, 700 / 300))
    ax = fig.add_axes([0.09, 0.21, 0.88, 0.70])
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.03, 1.08)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['PT start', 'PT end'])
    for t, ha in zip(ax.get_xticklabels(), ('left', 'right')):
        t.set_ha(ha)
    ax.set_yticks([])
    ax.set_ylabel('expression (scaled)', labelpad=2)
    return fig, ax


def legend_row(fig, items):
    """One row of legend entries under the axes: (label, colour, linestyle or None, italic)."""
    from matplotlib.lines import Line2D

    renderer = fig.canvas.get_renderer()
    x = 0.09
    for label, col, ls, italic in items:
        if ls:
            fig.lines.append(Line2D([x, x + 0.05], [0.06, 0.06], transform=fig.transFigure, color=col, lw=1.4, ls=ls))
            x += 0.065
        t = fig.text(x, 0.035, label, fontsize=7.5, color=col, style='italic' if italic else 'normal', va='bottom')
        x += t.get_window_extent(renderer).width / fig.bbox.width + 0.04


def masked(curve, lo, hi, grid):
    return np.where((grid >= lo) & (grid <= hi), curve, np.nan)


def panel_ds(args, plt, out):
    from pseudospace.handoffs import crossing, scaled_profile

    root = Path(args.results_root) / 'paper_figures'
    # c: a sharp and a gradual control-PT gene (Fig. 4 screen); each mouse inside its own covered range,
    # scaled to its own maximum
    z = np.load(root / 'fig4/screen/curves.npz')
    grid, genes = z['grid'], list(z['genes'])
    b = pd.read_csv(root / 'fig4/screen/axes.csv').set_index('pattern').loc['AAAAA']
    fig, ax = ds_axes(plt)
    bands(ax, b.b12, b.b23, 1.01)
    items = []
    for gene, col in (('Slc5a2', '#3F5FA8'), ('Kap', '#A9844B')):
        i = genes.index(gene)
        for m in ('Ctrl1A2', 'Ctrl1A4'):
            bx = z[f'bin_x|{m}'][i]
            c = masked(np.exp(z[f'curve|{m}'][i]), np.nanmin(bx), np.nanmax(bx), grid)
            ax.plot(grid, c / np.nanmax(c), color=col, lw=1.4, alpha=0.9)
        items.append((gene, col, None, True))
    legend_row(fig, items)
    save(fig, out / 'f1_ds_gradients.png')
    plt.close(fig)

    # d: SGLT2 (Slc5a2) -> SGLT1 (Slc5a1) in mouse and human (Fig. 5 screen), scaled as the screen scales them;
    # species curves are drawn only where both sections are covered
    z = np.load(root / 'fig5/screen/all_curves.npz', allow_pickle=True)
    grid, genes = z['grid'], [str(g) for g in z['genes']]
    g5 = pd.read_csv(root / 'fig5/screen/genes.csv').set_index('gene')
    a = pd.read_csv(root / 'fig5/screen/axes.csv').set_index('pattern').loc[g5.loc['Slc5a2', 'pattern']]
    ho = pd.read_csv(root / 'fig5/screen/handoffs.csv')
    pair = ho[(ho.early == 'Slc5a2') & (ho.late == 'Slc5a1')].iloc[0]
    fig, ax = ds_axes(plt)
    bands(ax, a.b12, a.b23, 1.01)
    switch = {}
    for species, members, col in (('mouse', ('Ctrl1A2', 'Ctrl1A4'), MOUSE), ('human', ('HUK1_COR1', 'HUK1_MED1'), HUMAN)):
        prof = {}
        for gene in ('Slc5a2', 'Slc5a1'):
            i = genes.index(gene)
            prof[gene] = np.mean([scaled_profile(masked(z[f'curve|{m}'][i], a[f'x_lo_{m}'], a[f'x_hi_{m}'], grid))
                                  for m in members], 0)
        ax.plot(grid, prof['Slc5a2'], color=col, lw=1.5)
        ax.plot(grid, prof['Slc5a1'], color=col, lw=1.5, ls=(0, (3, 1.6)))
        x = float(np.mean([pair[f'x_{m}'] for m in members]))
        switch[species] = {'mean_of_sections': x, 'crossing_of_means': crossing(prof['Slc5a2'], prof['Slc5a1'], grid)}
        ax.axvline(x, color=col, lw=0.9, ls=(0, (2, 2)), ymax=0.93)
    legend_row(fig, [('Slc5a2', '#444444', '-', True), ('Slc5a1', '#444444', (0, (3, 1.6)), True),
                     ('mouse', MOUSE, None, False), ('human', HUMAN, None, False)])
    save(fig, out / 'f1_ds_species.png')
    plt.close(fig)
    return {'gradients': ['Slc5a2', 'Kap'], 'handoff': 'Slc5a2 -> Slc5a1', 'switch': switch}


def main(argv=None):
    args = parse_args(argv)
    plt = style()
    out = Path(args.results_root) / args.out / 'panels'
    cache = Path(args.results_root) / args.out / 'cache'
    out.mkdir(parents=True, exist_ok=True)
    want = set(args.panels.split(','))
    record = {}
    d = windows(args) if want & {'seg', 'agg', 'plane'} else None
    if 'slides' in want:
        record['slides'] = panel_slides(args, plt, out, cache)
    if 'seg' in want:
        record['seg'] = panel_seg(args, plt, out, cache)
    if 'agg' in want:
        record['agg'] = panel_agg(args, plt, out, cache, d)
    if 'harm' in want:
        record['harm'] = panel_harm(args, plt, out)
    if 'graph' in want:
        record['graph'] = panel_graph(args, plt, out)
    if 'plane' in want:
        record['plane'] = panel_plane(args, plt, out, cache, d)
    if 'ds' in want:
        record['ds'] = panel_ds(args, plt, out)
    old = json.loads((out / 'record.json').read_text()) if (out / 'record.json').exists() else {}
    (out / 'record.json').write_text(json.dumps({**old, **record}, indent=1, default=str))
    print(json.dumps(record, indent=1, default=str)[:3000])


if __name__ == '__main__':
    main()
