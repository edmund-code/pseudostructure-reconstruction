"""Fig. 2 panels coloured by segmentation class, so every outline and fill uses one palette.

 heads    panel a thumbnails from one held-out Ctrl1A2 patch: H&E, class map, boundary map, glomerulus centres
 tiles    panel c gallery: outlines coloured by the predicted class
            held-out    Ctrl1A2 validation block, manual outlines in white under the model's, tile PQ
            unseen      Ctrl1A4 (no section of it was used in training)
            injured     IR2A4 window of fig2_gallery.py (most PT-S3 units)
            human       HUK1_MED1 window of fig2_gallery.py (the human model separates tubule from glomerulus only)
 maps     panel e: the four mouse sections, QC-retained cross-sections filled by class, the area outside the
          Visium HD capture area in grey; class counts for the bars

Classes: mouse model background/proximal/distal/collecting/glomerulus; human model background/tubule/glomerulus.
Manual labels Tubules/Tubules2/Tubules3/Glomeruli are proximal/distal/collecting/glomerulus (kidney_panoptic README).

Reads the segmentation GeoJSONs under --data-root, the manual labels (--manual), H&E under --he-root and the package
runs under --results-root; writes PNGs and JSON to <results-root>/paper_figures/fig2/classes/.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

REPO = next(p for p in Path(__file__).resolve().parents if (p / 'pseudospace').is_dir())
CLS = {'tubule_proximal': '#2F6FAE', 'tubule_distal': '#E3A21A', 'tubule_collecting': '#1E9E7A',
       'glomerulus': '#D98BB6', 'tubule': '#3A3A3A'}
MANUAL = {'Tubules': 'tubule_proximal', 'Tubules2': 'tubule_distal', 'Tubules3': 'tubule_collecting',
          'Glomeruli': 'glomerulus'}
SEGMENTATION = {'Ctrl1A2': 'Ctrl_1A2_v4.geojson', 'Ctrl1A4': 'Ctrl_1A4_v4.geojson', 'IR2A2': 'IR_2A2_v4.geojson',
                'IR2A4': 'IR_2A4_v4.geojson', 'HUK1_MED1': 'HUK1_MED1_v2.geojson'}
HE = {'Ctrl1A2': 'Ctrl_1A2', 'Ctrl1A4': 'Ctrl_1A4', 'IR2A2': 'IR_2A2', 'IR2A4': 'IR_2A4', 'HUK1_MED1': 'HUK1_MED_1'}
# validation blocks of Ctrl1A2 (kidney_panoptic manifest_v4.csv, roi of Ctrl_1A2_val), full-resolution pixels
VAL_BLOCKS = [(5613, 0, 7661, 2048), (9709, 2048, 11757, 4096), (5613, 4096, 7661, 6144),
              (9709, 6144, 11757, 8192), (5613, 8192, 7661, 10240), (9709, 10240, 11757, 11574)]
UM_PER_PX = 0.4409
TILE_UM, PATCH_PX, OUT_PX = 250.0, 512, 900


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--results-root', default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT', str(REPO / 'results')))
    p.add_argument('--data-root', default=os.environ.get('PSEUDOSPACE_DATA_ROOT', str(REPO / 'data')))
    p.add_argument('--he-root', default=os.environ.get('PSEUDOSPACE_HE_ROOT'), help='folder of the full-resolution H&E TIFFs')
    p.add_argument('--manual', default=os.environ.get('PSEUDOSPACE_CTRL1A2_MANUAL'),
                   help='manual outlines of the Ctrl1A2 validation blocks (Ctrl_1A2_final.geojson)')
    p.add_argument('--out', default='paper_figures/fig2/classes', help='relative to the results root')
    p.add_argument('--panels', default='heads,tiles,maps')
    return p.parse_args(argv)


# ------------------------------------------------------------------------------------------------- inputs
def features(path):
    """[(feature_index, class name, shapely polygon)] of a GeoJSON; the largest part of a multipolygon."""
    import shapely.geometry as sg

    g = json.loads(Path(path).read_text())
    out = []
    for i, f in enumerate(g['features'] if isinstance(g, dict) else g):
        geom = sg.shape(f['geometry'])
        if geom.is_empty:
            continue
        geom = geom if geom.geom_type == 'Polygon' else max(geom.geoms, key=lambda q: q.area)
        c = (f.get('properties') or {}).get('classification')
        name = c.get('name') if isinstance(c, dict) else c
        out.append((i, MANUAL.get(name, name), geom))
    return out


def in_window(feats, x0, y0, side):
    import shapely

    box = shapely.box(x0, y0, x0 + side, y0 + side)
    return [f for f in feats if f[2].intersects(box)]


def he_crop(args, sample, x0, y0, side, cache: Path):
    from PIL import Image

    cache.mkdir(parents=True, exist_ok=True)
    f = cache / f'{sample}_{int(x0)}_{int(y0)}_{int(side)}.png'
    if not f.exists():
        subprocess.run(['vips', 'crop', str(Path(args.he_root) / f'{HE[sample]}.tif'), str(f), str(int(x0)),
                        str(int(y0)), str(int(side)), str(int(side))], check=True)
    return Image.open(f).convert('RGB')


def rgb(hexc):
    return tuple(int(hexc[i:i + 2], 16) for i in (1, 3, 5))


def outline_tile(img, feats, x0, y0, side, px, width=6, manual=None):
    """H&E window resized to px, manual outlines (white) under the model's class-coloured outlines."""
    from PIL import Image, ImageDraw

    s = px / side
    img = img.resize((px, px), Image.LANCZOS)
    draw = ImageDraw.Draw(img)

    def ring(geom):
        return [((x - x0) * s, (y - y0) * s) for x, y in geom.exterior.coords]

    for _, _, geom in manual or []:
        draw.line(ring(geom), fill=(255, 255, 255), width=width + 4, joint='curve')
    for _, cls, geom in feats:
        draw.line(ring(geom), fill=rgb(CLS[cls]), width=width, joint='curve')
    return img


def pq(pred, gt):
    """Class-agnostic panoptic quality of two polygon lists (match at IoU > 0.5) and class agreement of matches."""
    tp, ious, same = 0, [], 0
    used = set()
    for _, gc, g in gt:
        best, bj = 0.0, None
        for j, (_, pc, p) in enumerate(pred):
            if j in used or not g.intersects(p):
                continue
            iou = g.intersection(p).area / g.union(p).area
            if iou > best:
                best, bj = iou, j
        if best > 0.5:
            tp += 1
            ious.append(best)
            used.add(bj)
            same += pred[bj][1] == gc
    fp, fn = len(pred) - tp, len(gt) - tp
    return {'pq': float(np.sum(ious) / (tp + 0.5 * fp + 0.5 * fn)) if tp else 0.0, 'tp': tp, 'fp': fp, 'fn': fn,
            'class_agreement': same / tp if tp else np.nan}


def best_window(feats, side, bounds, step, classes):
    """Window inside bounds with every class present, then the best balance (up to 8 per class), then the most."""
    cx = np.array([f[2].centroid.x for f in feats])
    cy = np.array([f[2].centroid.y for f in feats])
    code = np.array([classes.index(f[1]) if f[1] in classes else -1 for f in feats])
    best = None
    for (bx0, by0, bx1, by1) in bounds:
        for x in np.arange(bx0, bx1 - side + 1, step):
            kx = (cx >= x) & (cx <= x + side)
            for y in np.arange(by0, by1 - side + 1, step):
                k = kx & (cy >= y) & (cy <= y + side) & (code >= 0)
                n = np.bincount(code[k], minlength=len(classes))
                key = (int((n > 0).sum()), int(np.minimum(n, 8).sum()), int(k.sum()))
                if best is None or key > best[0]:
                    best = (key, (float(x), float(y)))
    return best[1]


# ------------------------------------------------------------------------------------------------- panels
def held_out_window(args, feats_manual):
    side = TILE_UM / UM_PER_PX
    return best_window(feats_manual, side, VAL_BLOCKS, 60,
                       ['tubule_proximal', 'tubule_distal', 'tubule_collecting', 'glomerulus']), side


def panel_heads(args, out, cache):
    """Thumbnails for the U-Net diagram: one held-out patch (512 px) and the maps the three heads predict."""
    from PIL import Image, ImageDraw, ImageFilter

    manual = features(args.manual)
    (wx, wy), side = held_out_window(args, manual)
    cx, cy = wx + side / 2, wy + side / 2
    x0, y0 = cx - PATCH_PX / 2, cy - PATCH_PX / 2
    pred = in_window(features(Path(args.data_root) / SEGMENTATION['Ctrl1A2']), x0, y0, PATCH_PX)
    he = he_crop(args, 'Ctrl1A2', x0, y0, PATCH_PX, cache)
    px = 360
    s = px / PATCH_PX
    he.resize((px, px), Image.LANCZOS).save(out / 'head_he.png')

    def poly(geom):
        return [((x - x0) * s, (y - y0) * s) for x, y in geom.exterior.coords]

    cls_img = Image.new('RGB', (px, px), 'white')
    d = ImageDraw.Draw(cls_img)
    for _, c, g in pred:
        d.polygon(poly(g), fill=rgb(CLS[c]))
    cls_img.save(out / 'head_classes.png')
    bnd = Image.new('RGB', (px, px), 'white')
    d = ImageDraw.Draw(bnd)
    for _, c, g in pred:
        d.polygon(poly(g), fill=(203, 211, 225))
        d.line(poly(g), fill=(74, 90, 122), width=3, joint='curve')
    bnd.save(out / 'head_boundary.png')
    grey = he.resize((px, px), Image.LANCZOS).convert('L').convert('RGB')
    heat = np.zeros((px, px))
    yy, xx = np.mgrid[0:px, 0:px]
    for _, c, g in pred:
        if c == 'glomerulus':
            gx, gy = (g.centroid.x - x0) * s, (g.centroid.y - y0) * s
            heat = np.maximum(heat, np.exp(-((xx - gx) ** 2 + (yy - gy) ** 2) / (2 * 14 ** 2)))
    base = np.asarray(grey, float) * 0.75
    mag = np.array(rgb(CLS['glomerulus']), float)
    blend = base * (1 - heat[..., None]) + mag * heat[..., None]
    Image.fromarray(blend.astype(np.uint8)).save(out / 'head_centre.png')
    return {'patch_px': [x0, y0, PATCH_PX], 'n_polygons': len(pred),
            'n_glomeruli': int(sum(c == 'glomerulus' for _, c, _ in pred))}


def panel_tiles(args, out, cache):
    manual_all = features(args.manual)
    (wx, wy), side = held_out_window(args, manual_all)
    meta = {}
    gallery = pd.read_csv(Path(args.results_root) / 'paper_figures/fig2/gallery/tiles.csv').set_index('sample')
    windows = {'held_out': ('Ctrl1A2', wx, wy)}
    ctrl4 = features(Path(args.data_root) / SEGMENTATION['Ctrl1A4'])
    xs = [f[2].centroid.x for f in ctrl4]
    ys = [f[2].centroid.y for f in ctrl4]
    windows['unseen'] = ('Ctrl1A4', *best_window(ctrl4, side, [(min(xs), min(ys), max(xs), max(ys))], 400,
                                                ['tubule_proximal', 'tubule_distal', 'tubule_collecting', 'glomerulus']))
    windows['injured'] = ('IR2A4', float(gallery.loc['IR2A4', 'x0_px']), float(gallery.loc['IR2A4', 'y0_px']))
    windows['human'] = ('HUK1_MED1', float(gallery.loc['HUK1_MED1', 'x0_px']), float(gallery.loc['HUK1_MED1', 'y0_px']))
    for key, (sample, x0, y0) in windows.items():
        feats = ctrl4 if sample == 'Ctrl1A4' else features(Path(args.data_root) / SEGMENTATION[sample])
        pred = in_window(feats, x0, y0, side)
        manual = in_window(manual_all, x0, y0, side) if key == 'held_out' else None
        img = outline_tile(he_crop(args, sample, x0, y0, side, cache), pred, x0, y0, side, OUT_PX, manual=manual)
        img.save(out / f'tile_{key}.png')
        m = {'sample': sample, 'x0_px': x0, 'y0_px': y0, 'tile_um': TILE_UM, 'n_predicted': len(pred),
             'classes': pd.Series([c for _, c, _ in pred]).value_counts().to_dict()}
        if manual:
            # score only objects whose centroid lies in the window, so cut-off edge objects do not count
            import shapely
            box = shapely.box(x0, y0, x0 + side, y0 + side)
            m['pq'] = pq([f for f in pred if box.contains(f[2].centroid)], [f for f in manual if box.contains(f[2].centroid)])
        meta[key] = m
    return meta


def panel_maps(args, out, cache):
    """Whole-section maps: QC-retained cross-sections filled by class; outside the capture area in grey."""
    import h5py
    from PIL import Image, ImageDraw

    with h5py.File(Path(args.results_root) / 'package_pseudospace/runs/aki/02_qc/units_qc.h5ad', 'r') as h:
        o = h['obs']

        def col(k):
            v = o[k]
            return v['categories'][:].astype(str)[v['codes'][:]] if isinstance(v, h5py.Group) else v[:]
        qc = pd.DataFrame({k: col(k) for k in ('sample', 'seg_class', 'qc_retained', 'feature_index')})
    qc = qc[qc.qc_retained.astype(bool)]
    meta = {}
    fig1 = Path(args.results_root) / 'paper_figures/fig1/cache'
    for sample in ('Ctrl1A2', 'Ctrl1A4', 'IR2A2', 'IR2A4'):
        tif = Path(args.he_root) / f'{HE[sample]}.tif'
        full_w = int(subprocess.check_output(['vipsheader', '-f', 'width', str(tif)]))
        th = cache / f'{sample}_map_he.png'
        if not th.exists():
            subprocess.run(['vips', 'thumbnail', str(tif), str(th), '1600'], check=True)
        img = Image.open(th).convert('RGB')
        w, h = img.size
        s = w / full_w
        a = np.asarray(img, float)
        e = max(2, int(0.006 * h))                       # scanner edge rows/columns
        a[:e], a[-e:], a[:, :e], a[:, -e:] = 247, 247, 247, 247
        cap = json.loads((fig1 / 'bins' / f'{sample}_capture.json').read_text())
        inside = np.zeros((h, w), bool)
        inside[max(int(cap['y0'] * s), 0):int(cap['y1'] * s), max(int(cap['x0'] * s), 0):int(cap['x1'] * s)] = True
        tissue = a.mean(2) < 225
        out_cap = tissue & ~inside                        # grey, desaturated: pink would read as glomerulus
        a[out_cap] = (0.35 * a[out_cap].mean(1) + 0.65 * 214)[:, None]
        a[tissue & inside] = 0.35 * a[tissue & inside] + 0.65 * 235      # faint H&E under the fills
        img = Image.fromarray(a.astype(np.uint8))
        draw = ImageDraw.Draw(img)
        keep = qc[qc['sample'] == sample]
        feats = {i: (c, g) for i, c, g in features(Path(args.data_root) / SEGMENTATION[sample])}
        for fi, cls in zip(keep.feature_index.astype(int), keep.seg_class):
            if fi in feats:
                _, g = feats[fi]
                draw.polygon([(x * s, y * s) for x, y in g.exterior.coords], fill=rgb(CLS[cls]))
        ys, xs = np.where(tissue)
        img = img.crop((max(xs.min() - 20, 0), max(ys.min() - 20, 0), min(xs.max() + 20, w), min(ys.max() + 20, h)))
        img.save(out / f'map_{sample}.png')
        meta[sample] = {'retained': int(len(keep)), 'classes': keep.seg_class.value_counts().to_dict()}
    return meta


def main(argv=None):
    args = parse_args(argv)
    out = Path(args.results_root) / args.out
    cache = out / 'cache'
    out.mkdir(parents=True, exist_ok=True)
    want = set(args.panels.split(','))
    rec = json.loads((out / 'record.json').read_text()) if (out / 'record.json').exists() else {}
    if 'heads' in want:
        rec['heads'] = panel_heads(args, out, cache)
    if 'tiles' in want:
        rec['tiles'] = panel_tiles(args, out, cache)
    if 'maps' in want:
        rec['maps'] = panel_maps(args, out, cache)
    rec['palette'] = CLS
    (out / 'record.json').write_text(json.dumps(rec, indent=1, default=str))
    print(json.dumps(rec, indent=1, default=str))


if __name__ == '__main__':
    main()
