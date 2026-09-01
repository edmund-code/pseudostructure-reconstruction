#!/usr/bin/env python
"""Render full-resolution detail crops from a finished WSI run, for visual QC.

    python scripts/render_crops.py --run-dir results/v2_ctrl_1a2 \
      --slide /path/to/slide.tif --n 4 --size 1024

The slide-wide overlay is downsampled ~16x and is only good for spotting gross failures
(instances on glass, whole regions missed). Judging whether two touching tubules were
actually separated needs native resolution, which is what this produces: side-by-side raw
image and instance outlines at 1:1.

Crops are chosen at the densest tissue locations unless ``--at x,y`` is given.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import _bootstrap  # noqa: F401
import cv2
import numpy as np

from kidney_panoptic.data.wsi import WSIReader

PALETTE = np.array([[0, 0, 0], [220, 50, 50], [50, 200, 80], [60, 90, 230],
                    [235, 180, 40], [190, 70, 210], [40, 200, 210]], np.uint8)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run-dir", required=True, help="dir containing instance_map.npz + summary.json")
    p.add_argument("--slide", required=True)
    p.add_argument("--n", type=int, default=4)
    p.add_argument("--size", type=int, default=1024)
    p.add_argument("--at", action="append", default=None, help="x,y (work px); repeatable")
    p.add_argument("--out", default=None)
    args = p.parse_args()

    run = Path(args.run_dir)
    with open(run / "summary.json") as f:
        summary = json.load(f)
    labels = np.load(run / "instance_map.npz")["instance_map"]
    cls_pairs = np.load(run / "instance_map.npz")["class_of_id"]
    class_of_id = {int(a): int(b) for a, b in cls_pairs} if len(cls_pairs) else {}
    names = summary["config"]["classes"]

    x0, y0 = summary["region"][0], summary["region"][1]
    reader = WSIReader(args.slide, level=summary["level"],
                       manual_mpp=summary["working_mpp"])

    spots = ([tuple(int(v) for v in s.split(",")) for s in args.at] if args.at
             else densest_spots(labels, args.size, args.n))
    out_dir = Path(args.out or run / "crops")
    out_dir.mkdir(parents=True, exist_ok=True)

    cls_lut = np.zeros(int(labels.max()) + 1, np.int32)
    for i, c in class_of_id.items():
        if i < len(cls_lut):
            cls_lut[i] = c

    for (cx, cy) in spots:
        s = args.size
        sub = labels[cy:cy + s, cx:cx + s]
        rgb = reader.read_region_work(x0 + cx, y0 + cy, sub.shape[1], sub.shape[0])

        edges = np.zeros(sub.shape, bool)
        edges[:-1] |= sub[:-1] != sub[1:]
        edges[:, :-1] |= sub[:, :-1] != sub[:, 1:]
        edges[1:] |= sub[1:] != sub[:-1]
        edges[:, 1:] |= sub[:, 1:] != sub[:, :-1]
        edges &= sub > 0

        colors = PALETTE[np.clip(cls_lut[np.clip(sub, 0, len(cls_lut) - 1)],
                                 0, len(PALETTE) - 1)]
        over = rgb.copy()
        fill = sub > 0
        over[fill] = (0.80 * over[fill] + 0.20 * colors[fill]).astype(np.uint8)
        over[edges] = colors[edges]

        panel = np.concatenate([rgb, over], axis=1)
        name = f"crop_x{cx}_y{cy}.png"
        cv2.imwrite(str(out_dir / name), cv2.cvtColor(panel, cv2.COLOR_RGB2BGR))
        n_here = len(set(np.unique(sub)) - {0})
        print(f"  {name}: {n_here} instances")

    reader.close()
    print(f"-> {out_dir}")


def densest_spots(labels: np.ndarray, size: int, n: int) -> list[tuple[int, int]]:
    """Pick n well-separated windows with the most distinct instances."""
    H, W = labels.shape
    step = size
    scored = []
    for y in range(0, max(1, H - size), step):
        for x in range(0, max(1, W - size), step):
            sub = labels[y:y + size:4, x:x + size:4]
            scored.append((len(set(np.unique(sub)) - {0}), x, y))
    scored.sort(reverse=True)
    picks: list[tuple[int, int]] = []
    for _, x, y in scored:
        if all(abs(x - px) >= size or abs(y - py) >= size for px, py in picks):
            picks.append((x, y))
        if len(picks) >= n:
            break
    return picks


if __name__ == "__main__":
    sys.exit(main())
