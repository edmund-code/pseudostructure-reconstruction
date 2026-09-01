#!/usr/bin/env python
"""Build and inspect the tissue mask for a slide.

    python scripts/make_tissue_mask.py --slide /path/to/slide.tif \
      --config configs/default.yaml --out results/ctrl_1a2_tissue

Writes ``tissue_mask.png`` (downsampled), ``thumbnail.png`` and ``tissue.json``. Use it to
check the saturation/gray thresholds before training or inference on a new stain or
scanner: calling tissue "glass" trains the model to suppress real structures, so the
defaults are deliberately conservative in the other direction.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import _bootstrap  # noqa: F401
import cv2
import numpy as np

from kidney_panoptic.data.tissue import build_tissue_mask
from kidney_panoptic.data.wsi import WSIReader
from kidney_panoptic.utils.config import load_config


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--slide", required=True)
    p.add_argument("--config", default=None)
    p.add_argument("--out", required=True)
    p.add_argument("--downsample", type=int, default=32)
    p.add_argument("--level", type=int, default=None)
    p.add_argument("--set", nargs="*", dest="overrides", default=[])
    args = p.parse_args()

    cfg = load_config(args.config, args.overrides) if args.config else {"wsi": {}, "tissue": None}
    wcfg = dict(cfg.get("wsi", {}) or {})
    level = args.level if args.level is not None else wcfg.get("level", 1)

    try:
        reader = WSIReader(args.slide, level=level, manual_mpp=wcfg.get("manual_mpp"),
                           target_mpp=float(wcfg.get("target_mpp", 0.44068)))
    except ValueError:
        print(f"level {level} unavailable — falling back to level 0")
        reader = WSIReader(args.slide, level=0, manual_mpp=wcfg.get("manual_mpp"),
                           target_mpp=float(wcfg.get("target_mpp", 0.44068)))

    W, H = reader.level_dimensions_work
    ds = args.downsample
    sw, sh = max(1, W // ds), max(1, H // ds)
    thumb = np.zeros((sh, sw, 3), np.uint8)
    band = max(ds, 4096)
    for y in range(0, H, band):
        h = min(band, H - y)
        strip = reader.read_region_work(0, y, W, h)
        th = max(1, h // ds)
        thumb[y // ds:y // ds + th] = cv2.resize(strip, (sw, th),
                                                 interpolation=cv2.INTER_AREA)

    tcfg = dict(cfg.get("tissue") or {})
    for k in ("open_px", "close_px", "dilate_px"):
        if k in tcfg:
            tcfg[k] = max(1, int(round(tcfg[k] / ds)))
    for k in ("min_object_px", "fill_holes_px"):
        if k in tcfg:
            tcfg[k] = max(1, int(tcfg[k] / (ds * ds)))

    mask = build_tissue_mask(thumb, tcfg or None)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    cv2.imwrite(str(out / "thumbnail.png"), cv2.cvtColor(thumb, cv2.COLOR_RGB2BGR))
    cv2.imwrite(str(out / "tissue_mask.png"), (mask * 255).astype(np.uint8))

    overlay = thumb.copy()
    overlay[~mask] = (0.5 * overlay[~mask] + 0.5 * np.array([255, 0, 0])).astype(np.uint8)
    cv2.imwrite(str(out / "overlay.png"), cv2.cvtColor(overlay, cv2.COLOR_RGB2BGR))

    mpp2 = reader.working_mpp ** 2
    info = {
        "slide": str(args.slide),
        "level": reader.level,
        "working_mpp": reader.working_mpp,
        "size_work": [W, H],
        "downsample": ds,
        "tissue_fraction": float(mask.mean()),
        "tissue_area_mm2": float(mask.sum() * (ds ** 2) * mpp2 / 1e6),
        "params": tcfg,
    }
    with open(out / "tissue.json", "w") as f:
        json.dump(info, f, indent=2)
    print(json.dumps(info, indent=2))
    print(f"-> {out}")
    reader.close()


if __name__ == "__main__":
    sys.exit(main())
