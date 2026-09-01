#!/usr/bin/env python
"""Region evaluation: PQ/SQ/DQ, AJI, semantic Dice, split/merge, glass, shift, canary.

    python scripts/eval_regions.py --checkpoint runs/v4_panoptic/best_pq.pt \
      --config configs/data.yaml --manifest data/manifest_eval_nx050.csv \
      --split val --size 4096 --gpu 1 --out results/v2_region_eval.json

Regions are taken from the annotated area of each evaluation slide. Ground truth is
rasterized from the same GeoJSON the training set came from, and predictions falling in the
IGNORED part of the mask are dropped before scoring — on a partially annotated slide a
prediction on an undrawn but real tubule is not a false positive. Predictions on TRUSTED
background (glass, the ring around each annotation) are kept and DO count against the score.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import _bootstrap  # noqa: F401
import numpy as np
import torch
from shapely.affinity import scale as shapely_scale
from shapely.geometry import box

from kidney_panoptic.data.classes import ClassSpec
from kidney_panoptic.data.geojson_io import load_geojson_instances
from kidney_panoptic.data.ignore import build_valid_mask
from kidney_panoptic.data.manifest import load_manifest
from kidney_panoptic.data.rasterize import rasterize_instances
from kidney_panoptic.data.targets import semantic_target
from kidney_panoptic.data.tissue import build_tissue_mask
from kidney_panoptic.data.wsi import WSIReader
from kidney_panoptic.eval.diagnostics import constant_input_canary, glass_rate, shift_stability
from kidney_panoptic.eval.instance_metrics import (
    evaluate_instances, evaluate_per_class, filter_by_valid, subtype_merge_breakdown,
)
from kidney_panoptic.eval.semantic_metrics import semantic_scores
from kidney_panoptic.infer.tiled import TiledPredictor, array_fetcher
from kidney_panoptic.utils.checkpoint import build_from_checkpoint
from kidney_panoptic.utils.config import apply_overrides, load_config, resolve_path

# Fraction of ADJACENT tubule pairs that differ in subtype, measured on Nx_4wk_050,
# Nx_4wk_001 and Nx_24h_1 (same/cross = 99/38, 188/63, 190/79). This is the null for the
# merge breakdown below: a cross-subtype merge fraction near this value means merging is
# subtype-independent, and a class-aware split cannot be the fix.
CROSS_SUBTYPE_NULL = 0.274


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--checkpoint", required=True)
    p.add_argument("--config", default=None, help="only for the manifest; model cfg comes "
                                                  "from the checkpoint")
    p.add_argument("--manifest", default=None)
    p.add_argument("--split", default="val")
    p.add_argument("--region", action="append", default=None,
                   help="x,y in working px; repeatable. Default: auto-pick annotated regions")
    p.add_argument("--size", type=int, default=4096)
    p.add_argument("--n-regions", type=int, default=3)
    p.add_argument("--gpu", type=int, default=None)
    p.add_argument("--out", default="results/v2_region_eval.json")
    p.add_argument("--skip-shift", action="store_true")
    p.add_argument("--set", nargs="*", dest="overrides", default=[])
    return p.parse_args()


def main():
    args = parse_args()
    device = f"cuda:{args.gpu}" if args.gpu is not None else (
        "cuda" if torch.cuda.is_available() else "cpu")

    model, cfg, spec = build_from_checkpoint(args.checkpoint, device=device)
    if args.overrides:
        apply_overrides(cfg, args.overrides)

    manifest_path = args.manifest
    if manifest_path is None:
        if not args.config:
            raise SystemExit("Pass --manifest, or --config so it can be read from there.")
        fcfg = load_config(args.config)
        manifest_path = resolve_path(fcfg, fcfg["manifest"])
    rows = load_manifest(manifest_path, splits=(args.split,))
    print(f"evaluating {len(rows)} {args.split} slide(s): "
          f"{[r['slide_id'] for r in rows]}")

    predictor = TiledPredictor(model, cfg, spec, device=device)

    results = {"checkpoint": str(args.checkpoint), "size": args.size, "slides": {}}
    pooled = {"TP": 0, "FP": 0, "FN": 0, "sq_sum": 0.0, "n_gt": 0, "n_pred": 0,
              "aji_sum": 0.0, "n_split": 0, "n_merged_gt": 0}

    for row in rows:
        res = eval_slide(row, cfg, spec, predictor, args)
        results["slides"][row["slide_id"]] = res
        for r in res["regions"]:
            m = r["instance"]["overall"]
            pooled["TP"] += m["TP"]; pooled["FP"] += m["FP"]; pooled["FN"] += m["FN"]
            pooled["sq_sum"] += m["SQ"] * m["TP"]
            pooled["aji_sum"] += m["AJI"] * m["n_gt"]
            pooled["n_gt"] += m["n_gt"]; pooled["n_pred"] += m["n_pred"]
            pooled["n_split"] += m["n_split"]; pooled["n_merged_gt"] += m["n_merged_gt"]

    tp, fp, fn = pooled["TP"], pooled["FP"], pooled["FN"]
    sq = pooled["sq_sum"] / tp if tp else 0.0
    dq = tp / (tp + 0.5 * fp + 0.5 * fn) if (tp + fp + fn) else 0.0
    results["pooled"] = {
        "PQ": sq * dq, "SQ": sq, "DQ": dq,
        "AJI": pooled["aji_sum"] / max(pooled["n_gt"], 1),
        "split_rate": pooled["n_split"] / max(pooled["n_gt"], 1),
        "merge_rate": pooled["n_merged_gt"] / max(pooled["n_gt"], 1),
        "TP": tp, "FP": fp, "FN": fn,
        "n_gt": pooled["n_gt"], "n_pred": pooled["n_pred"],
    }

    # Pool the merge breakdown across every region and compare with the adjacency null.
    same = cross = 0
    all_groups = []
    for s in results["slides"].values():
        for r in s["regions"]:
            sm = r.get("subtype_merge") or {}
            same += sm.get("n_same_subtype", 0)
            cross += sm.get("n_cross_subtype", 0)
            all_groups += sm.get("groups", [])
    tot = same + cross
    results["subtype_merge_pooled"] = {
        "n_merging_preds": tot,
        "n_same_subtype": same,
        "n_cross_subtype": cross,
        "cross_subtype_fraction": (cross / tot) if tot else None,
        "null_cross_subtype_adjacency": CROSS_SUBTYPE_NULL,
        "groups": all_groups,
    }
    if tot:
        print(f"\nmerge breakdown by raw label: {tot} merging predictions — "
              f"{cross} cross-subtype, {same} same-subtype "
              f"({100 * cross / tot:.0f}% cross vs {100 * CROSS_SUBTYPE_NULL:.0f}% null)")
        for g in all_groups[:10]:
            print(f"    pred {g['pred_id']:>4} absorbed {g['labels']}")
    else:
        print("\nmerge breakdown: no merging predictions found")

    print("\nconstant-input canary (tissue gate disabled)...")
    results["canary"] = constant_input_canary(predictor, spec)
    print(f"  instances on constant input: "
          f"{results['canary']['per_level']}  -> "
          f"{'PASS' if results['canary']['passed'] else 'FAIL'}")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        json.dump(results, f, indent=2)

    p = results["pooled"]
    print(f"\n{'=' * 72}")
    print(f"POOLED over {p['n_gt']} GT objects: PQ {p['PQ']:.3f}  SQ {p['SQ']:.3f}  "
          f"DQ {p['DQ']:.3f}  AJI {p['AJI']:.3f}")
    print(f"  split {p['split_rate'] * 100:.1f}%  merge {p['merge_rate'] * 100:.1f}%  "
          f"TP {tp} FP {fp} FN {fn}")
    print(f"written -> {out}")


def eval_slide(row, cfg, spec: ClassSpec, predictor, args) -> dict:
    wcfg = cfg["wsi"]
    reader = WSIReader(row["image_path"], level=wcfg.get("level", 1),
                       target_mpp=float(wcfg.get("target_mpp", 0.44068)),
                       manual_mpp=wcfg.get("manual_mpp"), prefer=wcfg.get("prefer", "finer"))
    W, H = reader.level_dimensions_work

    raw = load_geojson_instances(row["geojson_path"])
    to_work = 1.0 if cfg["patch"].get("annotations_space", "full") == "level" \
        else reader.scale_full_to_work
    polys, names, bg_polys, roi_polys = [], [], [], []
    for inst in raw:
        role = spec.role(inst.class_name)
        if role == "ignore":
            continue
        g = inst.polygon if to_work == 1.0 else shapely_scale(
            inst.polygon, xfact=to_work, yfact=to_work, origin=(0, 0))
        if role == "instance":
            polys.append(g); names.append(inst.class_name)
        elif role == "background":
            bg_polys.append(g)
        else:
            roi_polys.append(g)

    regions = pick_regions(polys, args, W, H)
    print(f"\n[{row['slide_id']}] {len(regions)} region(s) of {args.size}px: {regions}")

    out_regions = []
    for (rx, ry) in regions:
        out_regions.append(eval_region(reader, rx, ry, args.size, polys, names,
                                       bg_polys, roi_polys, cfg, spec, predictor,
                                       skip_shift=args.skip_shift))
    reader.close()
    return {"slide_id": row["slide_id"], "regions": out_regions}


def pick_regions(polys, args, W, H) -> list[tuple[int, int]]:
    """Explicit --region wins; otherwise take the densest annotated windows."""
    if args.region:
        return [tuple(int(v) for v in r.split(",")[:2]) for r in args.region]
    size = args.size
    cents = np.array([[p.centroid.x, p.centroid.y] for p in polys])
    if len(cents) == 0:
        return [(0, 0)]
    picks, used = [], np.zeros(len(cents), bool)
    for _ in range(args.n_regions):
        best, best_n = None, -1
        for i in range(len(cents)):
            if used[i]:
                continue
            x = int(np.clip(cents[i, 0] - size / 2, 0, max(0, W - size)))
            y = int(np.clip(cents[i, 1] - size / 2, 0, max(0, H - size)))
            inside = ((cents[:, 0] >= x) & (cents[:, 0] < x + size) &
                      (cents[:, 1] >= y) & (cents[:, 1] < y + size) & ~used)
            if inside.sum() > best_n:
                best, best_n, best_inside = (x, y), int(inside.sum()), inside
        if best is None or best_n <= 0:
            break
        picks.append(best)
        used |= best_inside
    return picks or [(0, 0)]


def eval_region(reader, x, y, size, polys, names, bg_polys, roi_polys, cfg, spec,
                predictor, skip_shift: bool) -> dict:
    rgb = reader.read_region_work(x, y, size, size)
    tissue = build_tissue_mask(rgb, cfg.get("tissue"))
    window = box(x, y, x + size, y + size)

    hits = [i for i, p in enumerate(polys) if p.intersects(window)]
    gt, name_of_id = rasterize_instances([polys[i] for i in hits],
                                         [names[i] for i in hits], size, size,
                                         offset_xy=(x, y))
    gt_class = {i: spec.canonical_index(n) for i, n in name_of_id.items()}

    valid = build_valid_mask(
        gt, strategy=cfg.get("ignore_strategy", "dilated_gt_two_radius_plus_tissue"),
        tissue_mask=tissue,
        background_mask=_poly_mask(bg_polys, x, y, size),
        dense_roi_mask=_poly_mask(roi_polys, x, y, size),
        dilate_px_inner=int(cfg.get("annotation", {}).get("dilate_px_inner", 5)),
        dilate_px=int(cfg.get("annotation", {}).get("dilate_px", 30)),
    )

    res, maps = predictor.predict_and_decode(array_fetcher(rgb), 0, 0, size, size,
                                             size, size, tissue=tissue)

    inst = evaluate_per_class(gt, res.id_map, gt_class, res.class_of_id,
                              list(spec.names), valid=valid)
    sem_pred = np.argmax(maps["semantic_prob"], axis=0)
    sem = semantic_scores(semantic_target(gt, gt_class), sem_pred, spec.n_classes,
                          list(spec.names), valid=valid)
    gl = glass_rate(res.id_map, tissue)

    # Merging is the dominant failure mode; record WHICH objects get fused, keyed by their
    # RAW annotation label, so same-subtype and cross-subtype merges can be told apart.
    # Predictions are filtered by the valid mask first so this matches the scored numbers.
    scored_pred = filter_by_valid(res.id_map, valid, 0.5)
    subtype_merge = subtype_merge_breakdown(gt, scored_pred, name_of_id)

    out = {"origin": [x, y], "n_gt": int(len(name_of_id)),
           "instance": inst, "semantic": sem, "glass": gl,
           "subtype_merge": subtype_merge,
           "valid_fraction": float(valid.mean())}
    if not skip_shift:
        out["shift_stability"] = shift_stability(predictor, rgb, shift=128,
                                                 tissue_cfg=cfg.get("tissue"))
        print(f"  region ({x},{y}): PQ {inst['overall']['PQ']:.3f} "
              f"AJI {inst['overall']['AJI']:.3f} "
              f"split {inst['overall']['split_rate'] * 100:.0f}% "
              f"glass {gl['glass_instance_rate'] * 100:.1f}% "
              f"shift {out['shift_stability']['matched_fraction'] * 100:.0f}%")
    else:
        print(f"  region ({x},{y}): PQ {inst['overall']['PQ']:.3f} "
              f"AJI {inst['overall']['AJI']:.3f}")
    return out


def _poly_mask(polys, x, y, size):
    if not polys:
        return np.zeros((size, size), bool)
    hits = [p for p in polys if p.intersects(box(x, y, x + size, y + size))]
    if not hits:
        return np.zeros((size, size), bool)
    m, _ = rasterize_instances(hits, ["_"] * len(hits), size, size, offset_xy=(x, y))
    return m > 0


if __name__ == "__main__":
    sys.exit(main())
