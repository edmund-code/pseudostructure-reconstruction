"""Fig. 2c gallery tiles: H&E with the segmentation model's predicted outlines, on sections it never trained on.

 IR2A4      the --tile-um window holding the most PT-S3 units of the control/AKI fit (injured outer stripe)
 HUK1_MED1  the --tile-um window holding the most PT units of the human-vs-mouse fit (the human model trained on
            HUK1_COR1 blocks only)

The window is chosen from the unit centroids alone, before any image is looked at. PQ needs manual outlines,
which do not exist for these tiles yet. Protocol: docs/paper/fig3-validation-experiments.md.

  python analysis/paper/fig2_gallery.py --he-dir <dir with IR_2A4.tif ...>
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
sys.path.insert(0, str(Path(__file__).resolve().parent))
from fig3_common import RUNS  # noqa: E402

# sample -> (fit, H&E file, predicted-outline GeoJSON under the data root, labels whose units pick the window)
TILES = {'IR2A4': ('aki', 'IR_2A4.tif', 'segmentations/v4fix/IR2A4_v4fix.geojson', ['PT-S3']),
         'HUK1_MED1': ('hvm', 'HUK1_MED_1.tif', 'HUK1_MED1_v2.geojson', ['PT-S1', 'PT-S2', 'PT-S3'])}
UM_PER_PX = 0.4409   # full-resolution H&E of every section (Space Ranger scalefactors: 0.44074-0.44089)
OUTLINE = (255, 212, 0)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--data-root', default=os.environ.get('PSEUDOSPACE_DATA_ROOT', str(REPO / 'data')))
    p.add_argument('--results-root', default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT', str(REPO / 'results')))
    p.add_argument('--he-dir', default=os.environ.get('PSEUDOSPACE_HE_DIR'), required='PSEUDOSPACE_HE_DIR' not in os.environ)
    p.add_argument('--out', default='paper_figures/fig2/gallery', help='relative to the results root')
    p.add_argument('--tile-um', type=float, default=250.0)
    p.add_argument('--pixels', type=int, default=900, help='output tile width and height')
    return p.parse_args(argv)


def densest_window(xy: np.ndarray, size: float, step: float) -> tuple[float, float]:
    """Top-left corner of the size x size window (on a step grid) holding the most points."""
    x0, y0 = xy.min(0)
    gx = np.floor((xy[:, 0] - x0) / step).astype(int)
    gy = np.floor((xy[:, 1] - y0) / step).astype(int)
    grid = np.zeros((gy.max() + 1, gx.max() + 1))
    np.add.at(grid, (gy, gx), 1)
    k = int(round(size / step))
    c = np.pad(grid, ((1, 0), (1, 0))).cumsum(0).cumsum(1)
    win = c[k:, k:] - c[:-k, k:] - c[k:, :-k] + c[:-k, :-k]
    iy, ix = np.unravel_index(np.argmax(win), win.shape)
    return x0 + ix * step, y0 + iy * step


def rings(geometry) -> list[np.ndarray]:
    polys = geometry['coordinates'] if geometry['type'] == 'MultiPolygon' else [geometry['coordinates']]
    return [np.asarray(ring, float) for poly in polys for ring in poly]


def main(argv=None):
    import anndata as ad
    from PIL import Image, ImageDraw

    args = parse_args(argv)
    Image.MAX_IMAGE_PIXELS = None
    out = Path(args.results_root) / args.out
    out.mkdir(parents=True, exist_ok=True)
    size = args.tile_um / UM_PER_PX
    meta = []
    for sample, (run, he, geojson, labels) in TILES.items():
        obs = ad.read_h5ad(Path(args.results_root) / RUNS[run] / '03_integration/pass2.h5ad', backed='r').obs
        pick = obs[(obs['sample'] == sample) & obs['segment_class'].astype(str).isin(labels)]
        x, y = densest_window(pick[['centroid_x', 'centroid_y']].to_numpy(float), size, size / 10)
        box = (int(x), int(y), int(x + size), int(y + size))
        tile = Image.open(Path(args.he_dir) / he).crop(box).convert('RGB')
        scale = args.pixels / tile.width
        tile = tile.resize((args.pixels, args.pixels), Image.LANCZOS)
        draw = ImageDraw.Draw(tile)
        features = json.loads((Path(args.data_root) / geojson).read_text())
        features = features['features'] if isinstance(features, dict) else features
        n = 0
        for f in features:
            for ring in rings(f['geometry']):
                if ring[:, 0].max() < box[0] or ring[:, 0].min() > box[2] or ring[:, 1].max() < box[1] or ring[:, 1].min() > box[3]:
                    continue
                pts = [((px - box[0]) * scale, (py - box[1]) * scale) for px, py in ring]
                draw.line(pts + pts[:1], fill=OUTLINE, width=3, joint='curve')
                n += 1
        tile.save(out / f'tile_{sample}.png')
        meta.append({'sample': sample, 'fit': run, 'labels': '+'.join(labels), 'x0_px': box[0], 'y0_px': box[1],
                     'tile_um': args.tile_um, 'units_in_window': int(((pick.centroid_x.between(box[0], box[2])) &
                                                                      (pick.centroid_y.between(box[1], box[3]))).sum()),
                     'outlines_drawn': n})
    pd.DataFrame(meta).to_csv(out / 'tiles.csv', index=False)
    print(pd.DataFrame(meta).to_string())


if __name__ == '__main__':
    main()
