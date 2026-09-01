#!/usr/bin/env python
"""Whole-slide inference: tiled prediction -> windowed decode -> stitched GeoJSON.

    python scripts/wsi_infer.py \
      --checkpoint runs/v2_panoptic/best_pq.pt \
      --slide /path/to/slide.tif \
      --out-dir results/v2_ctrl_1a2 --save-viz --gpu 1

Outputs under ``--out-dir``:
    instances.geojson    QuPath-ready polygons, full-resolution coords, class + score
    instance_map.npz     int32 slide-wide instance ids (``--no-save-map`` to skip)
    tissue_mask.png      the tissue gate that was applied
    summary.json         counts per class, areas, glass rate, timings, config used
    overlay.png          downsampled visualisation (``--save-viz``)

``--dry-run`` prints the plan (slide geometry, tissue coverage, window and tile counts)
without loading the model, which is the cheap way to sanity-check a new slide's mpp.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from collections import Counter
from pathlib import Path

import _bootstrap  # noqa: F401
import cv2
import numpy as np
import torch

from kidney_panoptic.data.tissue import build_tissue_mask
from kidney_panoptic.data.wsi import WSIReader
from kidney_panoptic.eval.diagnostics import glass_rate
from kidney_panoptic.infer.tiled import TiledPredictor, decode_windows, reader_fetcher
from kidney_panoptic.postprocess.geojson_out import instance_polygons, write_geojson
from kidney_panoptic.postprocess.stitch import InstanceStitcher
from kidney_panoptic.utils.checkpoint import build_from_checkpoint
from kidney_panoptic.utils.config import apply_overrides


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--checkpoint", required=True)
    p.add_argument("--slide", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--gpu", type=int, default=None)
    p.add_argument("--save-viz", action="store_true")
    p.add_argument("--no-save-map", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--region", default=None,
                   help="x,y,w,h in working pixels — infer only this region")
    p.add_argument("--tta", choices=["none", "rot4", "d4"], default=None,
                   help="test-time augmentation; d4 is 8x the forward passes "
                        "(~65 min instead of ~8 on a full slide) for sharper walls")
    p.add_argument("--set", nargs="*", dest="overrides", default=[])
    return p.parse_args()


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    device = "cpu"
    if args.gpu is not None:
        device = f"cuda:{args.gpu}"
    elif torch.cuda.is_available():
        device = "cuda"

    # Config comes from the checkpoint so the class vocabulary matches the weights.
    if args.dry_run:
        from kidney_panoptic.utils.checkpoint import load_checkpoint
        cfg = dict(load_checkpoint(args.checkpoint)["config"])
        model = spec = None
    else:
        model, cfg, spec = build_from_checkpoint(args.checkpoint, device=device)
    if args.overrides:
        apply_overrides(cfg, args.overrides)
    if args.tta is not None:
        cfg.setdefault("infer", {})["tta"] = args.tta

    wcfg, icfg = cfg["wsi"], cfg["infer"]
    reader = open_slide(args.slide, wcfg)
    W, H = reader.level_dimensions_work
    print(f"slide {Path(args.slide).name}: {W}x{H} @ level {reader.level} "
          f"({reader.working_mpp:.5f} µm/px, scale_to_full {reader.scale_to_full})")

    x0, y0, rw, rh = (0, 0, W, H)
    if args.region:
        x0, y0, rw, rh = (int(v) for v in args.region.split(","))
        print(f"restricted to region ({x0},{y0}) {rw}x{rh}")

    t_start = time.time()
    tissue = build_slide_tissue(reader, cfg, W, H)
    tissue_frac = float(tissue.mean())
    print(f"tissue mask: {tissue_frac * 100:.1f}% of the slide "
          f"({time.time() - t_start:.0f}s)")

    window = int(icfg.get("decode_window", 2048))
    overlap = int(icfg.get("decode_overlap", 256))
    windows = [w for w in decode_windows(rw, rh, window, overlap)]
    live = [w for w in windows
            if tissue[y0 + w[1]:y0 + w[1] + w[3], x0 + w[0]:x0 + w[0] + w[2]].any()]
    print(f"decode windows: {len(live)} with tissue of {len(windows)} total")

    if args.dry_run:
        tile, ov = int(icfg.get("tile", 512)), int(icfg.get("overlap", 128))
        per = ((window + 2 * tile) // max(tile - ov, 1) + 1) ** 2
        print(f"dry run: ~{len(live) * per} tile forwards at tile={tile} overlap={ov}")
        reader.close()
        return

    predictor = TiledPredictor(model, cfg, spec, device=device,
                               amp=bool(cfg["train"].get("amp", True)))
    fetch = reader_fetcher(reader)

    def tissue_fetch(x, y, w, h):
        return tissue[max(0, y):y + h, max(0, x):x + w]

    stitcher = InstanceStitcher(rh, rw, float(icfg.get("stitch_iou", 0.5)))
    t0 = time.time()
    for i, (wx, wy, ww, wh) in enumerate(live):
        res, _ = predictor.predict_and_decode(
            fetch, x0 + wx, y0 + wy, ww, wh, W, H,
            tissue=tissue[y0 + wy:y0 + wy + wh, x0 + wx:x0 + wx + ww],
            tissue_fetch=tissue_fetch if icfg.get("skip_empty_tiles", True) else None,
        )
        stitcher.add_window(wx, wy, res.id_map, res.class_of_id, res.score_of_id)
        done = i + 1
        if done % 5 == 0 or done == len(live):
            rate = (time.time() - t0) / done
            print(f"  window {done}/{len(live)}  {res.n_instances} inst  "
                  f"{rate:.1f}s/window  eta {(len(live) - done) * rate / 60:.1f} min",
                  flush=True)

    labels, class_of_id, score_of_id = stitcher.finalize()
    infer_secs = time.time() - t0
    print(f"decoded {len(class_of_id)} instances in {infer_secs / 60:.1f} min")

    # ------------------------------------------------------------------ export
    counts = Counter(spec.names[c] for c in class_of_id.values())
    region_tissue = tissue[y0:y0 + rh, x0:x0 + rw]
    gr = glass_rate(labels, region_tissue,
                    min_tissue_frac=float(icfg.get("tissue_min_instance_frac", 0.5)))

    out_scale = output_scale(cfg, reader)
    feats = instance_polygons(
        labels, class_of_id, list(spec.names), score_of_id=score_of_id,
        offset_xy=(x0, y0), scale=out_scale,
        simplify_px=float(icfg.get("polygon_simplify_px", 1.0)),
        max_vertices=int(icfg.get("geojson_max_vertices", 400)),
        mpp=float(reader.working_mpp),
    )
    n_feats = write_geojson(feats, out_dir / "instances.geojson")
    print(f"wrote {n_feats} polygons -> {out_dir / 'instances.geojson'}")

    if not args.no_save_map and icfg.get("save_instance_map", True):
        np.savez_compressed(
            out_dir / "instance_map.npz", instance_map=labels,
            class_of_id=np.array(sorted(class_of_id.items()), dtype=np.int64),
            # Scores travel with the map so the GeoJSON can be re-exported (different
            # simplification, different coordinate space) without re-running inference.
            score_of_id=np.array([[i, score_of_id.get(i, 0.0)]
                                  for i in sorted(class_of_id)], dtype=np.float64),
        )
        print(f"wrote instance map -> {out_dir / 'instance_map.npz'}")

    ds = int(icfg.get("viz_downsample", 16))
    cv2.imwrite(str(out_dir / "tissue_mask.png"),
                (tissue[::ds, ::ds] * 255).astype(np.uint8))
    if args.save_viz:
        save_overlay(reader, labels, class_of_id, spec, x0, y0, rw, rh, ds,
                     out_dir / "overlay.png")
        print(f"wrote overlay -> {out_dir / 'overlay.png'}")

    mpp2 = reader.working_mpp ** 2
    areas = Counter({n: 0 for n in spec.names[1:]})
    area_px = np.bincount(labels.ravel())
    for i, c in class_of_id.items():
        if i < len(area_px):
            areas[spec.names[c]] += int(area_px[i])

    summary = {
        "slide": str(args.slide),
        "checkpoint": str(args.checkpoint),
        "level": reader.level,
        "working_mpp": reader.working_mpp,
        "scale_to_full": reader.scale_to_full,
        "slide_size_work": [W, H],
        "region": [x0, y0, rw, rh],
        "tissue_fraction": tissue_frac,
        "tissue_area_mm2": float(region_tissue.sum() * mpp2 / 1e6),
        "n_instances": len(class_of_id),
        "counts_by_class": dict(counts),
        "area_px_by_class": dict(areas),
        "area_mm2_by_class": {k: float(v * mpp2 / 1e6) for k, v in areas.items()},
        "density_per_mm2_tissue": {
            k: float(v / max(region_tissue.sum() * mpp2 / 1e6, 1e-9))
            for k, v in counts.items()
        },
        "glass": gr,
        "n_decode_windows": len(live),
        "infer_minutes": infer_secs / 60.0,
        "config": {k: v for k, v in cfg.items() if k != "_config_dir"},
    }
    with open(out_dir / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps({k: summary[k] for k in
                      ("n_instances", "counts_by_class", "density_per_mm2_tissue")}, indent=2))
    print(f"summary -> {out_dir / 'summary.json'}")
    reader.close()


def output_scale(cfg: dict, reader: WSIReader) -> float:
    """Work-pixel -> exported-coordinate scale factor.

    Getting this wrong is silent and confusing: predictions land at exactly the level
    downsample factor away from the annotations they should sit on top of, which looks like
    a catastrophic model failure in QuPath rather than a units bug.
    """
    space = str(cfg["infer"].get("output_space", "auto")).lower()
    if space == "auto":
        space = "level" if cfg["patch"].get("annotations_space") == "level" else "full"
    if space == "level":
        return 1.0
    if space == "full":
        return float(reader.scale_to_full)
    raise ValueError(f"infer.output_space must be auto|level|full, got '{space}'")


def open_slide(path: str, wcfg: dict) -> WSIReader:
    """Open the slide, with a clear message when a flat TIFF carries no real mpp."""
    try:
        return WSIReader(path, level=wcfg.get("level", 1),
                         target_mpp=float(wcfg.get("target_mpp", 0.44068)),
                         manual_mpp=wcfg.get("manual_mpp"),
                         anisotropy_tol=float(wcfg.get("anisotropy_tol", 0.02)),
                         prefer=wcfg.get("prefer", "finer"))
    except ValueError as e:
        if "level" in str(e):
            # A flat (single-level) export has no level 1; read level 0 instead.
            return WSIReader(path, level=0,
                             target_mpp=float(wcfg.get("target_mpp", 0.44068)),
                             manual_mpp=wcfg.get("manual_mpp"),
                             prefer=wcfg.get("prefer", "finer"))
        raise


def build_slide_tissue(reader: WSIReader, cfg: dict, W: int, H: int) -> np.ndarray:
    """Tissue mask at 1/N resolution, upsampled — full-res would be needlessly slow."""
    ds = int(cfg["infer"].get("tissue_downsample", 32))
    small_w, small_h = max(1, W // ds), max(1, H // ds)
    thumb = np.zeros((small_h, small_w, 3), np.uint8)
    band = max(ds, 4096 // ds * ds)
    for y in range(0, H, band):
        h = min(band, H - y)
        strip = reader.read_region_work(0, y, W, h)
        sh = max(1, h // ds)
        thumb[y // ds:y // ds + sh] = cv2.resize(strip, (small_w, sh),
                                                 interpolation=cv2.INTER_AREA)
    small = build_tissue_mask(thumb, _scaled_tissue_cfg(cfg.get("tissue"), ds))
    return cv2.resize(small.astype(np.uint8), (W, H),
                      interpolation=cv2.INTER_NEAREST).astype(bool)


def _scaled_tissue_cfg(tcfg: dict | None, ds: int) -> dict:
    """Morphology radii are in pixels — rescale them for the downsampled thumbnail."""
    out = dict(tcfg or {})
    for k in ("open_px", "close_px", "dilate_px"):
        if k in out:
            out[k] = max(1, int(round(out[k] / ds)))
    for k in ("min_object_px", "fill_holes_px"):
        if k in out:
            out[k] = max(1, int(out[k] / (ds * ds)))
    return out


def save_overlay(reader, labels, class_of_id, spec, x0, y0, rw, rh, ds, path):
    """Downsampled RGB with per-class instance boundaries drawn over it."""
    thumb = cv2.resize(reader.read_region_work(x0, y0, rw, rh),
                       (max(1, rw // ds), max(1, rh // ds)), interpolation=cv2.INTER_AREA)
    small = labels[::ds, ::ds]
    small = small[:thumb.shape[0], :thumb.shape[1]]
    palette = np.array([[0, 0, 0], [220, 50, 50], [50, 200, 80], [60, 90, 230],
                        [235, 180, 40], [190, 70, 210], [40, 200, 210]], np.uint8)

    cls_map = np.zeros(int(max(class_of_id, default=0)) + 1, np.int32)
    for i, c in class_of_id.items():
        cls_map[i] = c
    cls_img = cls_map[np.clip(small, 0, len(cls_map) - 1)]

    edges = np.zeros(small.shape, bool)
    edges[:-1] |= small[:-1] != small[1:]
    edges[:, :-1] |= small[:, :-1] != small[:, 1:]
    edges &= small > 0

    out = thumb.copy()
    colors = palette[np.clip(cls_img, 0, len(palette) - 1)]
    fill = small > 0
    out[fill] = (0.72 * out[fill] + 0.28 * colors[fill]).astype(np.uint8)
    out[edges] = colors[edges]
    cv2.imwrite(str(path), cv2.cvtColor(out, cv2.COLOR_RGB2BGR))


if __name__ == "__main__":
    sys.exit(main())
