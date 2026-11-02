"""Fig. 3f: distance of every tubule cross-section from the renal capsule, measured on H&E, against pseudostructure.

Tissue mask: Otsu threshold on the saturation of the H&E at 1/16 scale, closed, holes filled, largest component.
Capsule: the part of the mask boundary within --capsule-um of the mask's convex hull; boundary deeper inside the
hull (hilum, cracks, cut faces) is not capsule. Depth: Euclidean distance from the unit centroid to the nearest
capsule pixel. The method never sees depth; nothing is selected on it. Protocol: docs/paper/fig3-validation-experiments.md.

Reads private data (H&E slides, package runs) and writes under --results-root (outside Git).

  python analysis/paper/fig3_depth.py --he-dir <dir with Ctrl_1A2.tif ...> --visium-root <mouse Visium HD root>
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import ndimage
from scipy.spatial import ConvexHull

REPO = next(p for p in Path(__file__).resolve().parents if (p / 'pseudospace').is_dir())
sys.path.insert(0, str(Path(__file__).resolve().parent))
from fig3_common import SECTIONS, load_positions  # noqa: E402

HE_FILES = {'Ctrl1A2': 'Ctrl_1A2.tif', 'Ctrl1A4': 'Ctrl_1A4.tif', 'IR2A2': 'IR_2A2.tif', 'IR2A4': 'IR_2A4.tif',
            'HUK1_COR1': 'HUK1_COR_1.tif', 'HUK1_MED1': 'HUK1_MED_1.tif'}
CONFIGS = {'aki': 'analysis/package_configs/aki', 'hvm': 'analysis/package_configs/human_vs_mouse'}
# Human pieces are cut blocks of cortex: on the H&E no outer edge is capsule (checked by eye, 2026-11-02), so
# depth is not defined for them and they are skipped rather than given a substitute edge.
NO_CAPSULE = {'HUK1_COR1', 'HUK1_MED1'}


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--data-root', default=os.environ.get('PSEUDOSPACE_DATA_ROOT', str(REPO / 'data')))
    p.add_argument('--results-root', default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT', str(REPO / 'results')))
    p.add_argument('--visium-root', default=os.environ.get('PSEUDOSPACE_VISIUM_ROOT'))
    p.add_argument('--he-dir', default=os.environ.get('PSEUDOSPACE_HE_DIR'), required='PSEUDOSPACE_HE_DIR' not in os.environ)
    p.add_argument('--out', default='paper_figures/fig3/depth', help='relative to the results root')
    p.add_argument('--factor', type=int, default=16, help='H&E downsampling for the mask (16 -> about 7 um/px)')
    p.add_argument('--capsule-um', type=float, default=100.0)
    return p.parse_args(argv)


def microns_per_pixel(args, run: str, sample: str) -> float:
    sheet = pd.read_csv(REPO / CONFIGS[run] / 'samples.csv', comment='#').set_index('sample')
    bins = sheet.loc[sample, 'bins_dir'].replace('$PSEUDOSPACE_VISIUM_ROOT', str(args.visium_root))
    bins = Path(bins) if Path(bins).is_absolute() else Path(args.data_root) / bins
    return float(json.loads((bins / 'spatial/scalefactors_json.json').read_text())['microns_per_pixel'])


def otsu(values: np.ndarray) -> float:
    hist, edges = np.histogram(values, 256)
    centres = (edges[:-1] + edges[1:]) / 2
    w0 = np.cumsum(hist)
    w1 = w0[-1] - w0
    m0 = np.cumsum(hist * centres) / np.maximum(w0, 1)
    m1 = (np.sum(hist * centres) - np.cumsum(hist * centres)) / np.maximum(w1, 1)
    return float(centres[np.argmax(w0 * w1 * (m0 - m1) ** 2)])


def tissue_mask(rgb: np.ndarray) -> np.ndarray:
    rgb = rgb.astype(float)
    sat = (rgb.max(2) - rgb.min(2)) / np.maximum(rgb.max(2), 1)
    mask = ndimage.binary_closing(sat > otsu(sat), iterations=3)
    mask = ndimage.binary_fill_holes(mask)
    lab, n = ndimage.label(mask)
    return lab == (np.argmax(np.bincount(lab.ravel())[1:]) + 1)


def capsule_and_depth(mask: np.ndarray, px_um: float, capsule_um: float):
    """Capsule pixels and the distance (um) of every pixel from the nearest one."""
    from matplotlib.path import Path as Polygon

    boundary = mask & ~ndimage.binary_erosion(mask)
    yx = np.argwhere(boundary)
    hull = Polygon(yx[ConvexHull(yx).vertices][:, ::-1])
    gy, gx = np.mgrid[:mask.shape[0], :mask.shape[1]]
    inside = hull.contains_points(np.c_[gx.ravel(), gy.ravel()], radius=1).reshape(mask.shape) | boundary
    capsule = boundary & (ndimage.distance_transform_edt(inside) * px_um <= capsule_um)
    return capsule, ndimage.distance_transform_edt(~capsule) * px_um


def main(argv=None):
    from PIL import Image
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    args = parse_args(argv)
    Image.MAX_IMAGE_PIXELS = None
    out = Path(args.results_root) / args.out
    out.mkdir(parents=True, exist_ok=True)
    positions = {run: load_positions(args.results_root, run) for run in set(SECTIONS.values())}
    tables, summary = [], []
    for sample, run in SECTIONS.items():
        if sample in NO_CAPSULE:
            summary.append({'sample': sample, 'status': 'no capsule on the section'})
            continue
        mpp = microns_per_pixel(args, run, sample)
        rgb = np.asarray(Image.open(Path(args.he_dir) / HE_FILES[sample]).reduce(args.factor))
        px_um = mpp * args.factor
        mask = tissue_mask(rgb)
        capsule, depth_map = capsule_and_depth(mask, px_um, args.capsule_um)
        import anndata as ad
        obs = ad.read_h5ad(Path(args.results_root) / f'package_pseudospace/runs/{run}/03_integration/pass2.h5ad',
                           backed='r').obs
        units = positions[run][positions[run]['sample'] == sample].join(obs[['centroid_x', 'centroid_y']])
        iy = np.clip((units.centroid_y / args.factor).astype(int), 0, mask.shape[0] - 1)
        ix = np.clip((units.centroid_x / args.factor).astype(int), 0, mask.shape[1] - 1)
        units['depth_um'] = depth_map[iy, ix]
        units['in_mask'] = mask[iy, ix]
        tables.append(units.assign(run=run))
        summary.append({'sample': sample, 'status': 'measured', 'um_per_mask_px': px_um,
                        'capsule_fraction_of_boundary': float(capsule.sum() / (mask & ~ndimage.binary_erosion(mask)).sum()),
                        'units_outside_mask': int((~units.in_mask).sum()), 'max_depth_um': float(depth_map[mask].max())})
        fig, ax = plt.subplots(figsize=(10, 10 * mask.shape[0] / mask.shape[1]))
        ax.imshow(rgb)
        ax.contour(mask, [0.5], colors='k', linewidths=0.6)
        cy, cx = np.nonzero(capsule)
        ax.scatter(cx, cy, s=0.4, c='red')
        sc = ax.scatter(ix, iy, s=0.3, c=units.depth_um, cmap='viridis')
        fig.colorbar(sc, ax=ax, shrink=0.6, label='depth (um)')
        ax.set_axis_off()
        fig.savefig(out / f'qc_{sample}.png', dpi=110, bbox_inches='tight')
        plt.close(fig)
    pd.concat(tables).to_parquet(out / 'depth.parquet')
    pd.DataFrame(summary).to_csv(out / 'sections.csv', index=False)
    print(pd.DataFrame(summary).to_string())


if __name__ == '__main__':
    main()
