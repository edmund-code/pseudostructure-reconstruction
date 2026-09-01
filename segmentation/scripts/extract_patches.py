#!/usr/bin/env python
"""Extract training patches from the slides listed in a manifest.

    python scripts/extract_patches.py --config configs/data.yaml

For each slide: read the chosen pyramid level natively, load its GeoJSON, rasterize
instances into an id map, build a tissue mask, and write patches with their auxiliary
masks. Semantic / boundary / center targets are NOT stored — they are regenerated from the
instance map after augmentation (see data/targets.py).
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import _bootstrap  # noqa: F401
import numpy as np
from shapely.affinity import scale as shapely_scale
from shapely.geometry import box

from kidney_panoptic.data.classes import ClassSpec
from kidney_panoptic.data.geojson_io import load_geojson_instances
from kidney_panoptic.data.manifest import load_manifest, manifest_summary
from kidney_panoptic.data.patching import iter_tiles, parse_rois, patch_in_rois, save_patch
from kidney_panoptic.data.rasterize import rasterize_instances
from kidney_panoptic.data.tissue import build_tissue_mask
from kidney_panoptic.data.wsi import WSIReader
from kidney_panoptic.utils.config import load_config, resolve_path
from kidney_panoptic.utils.seed import set_seed


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", required=True)
    p.add_argument("--manifest", default=None, help="overrides `manifest` in the config")
    p.add_argument("--out", default=None, help="overrides `paths.patches_dir`")
    p.add_argument("--slides", nargs="*", default=None, help="restrict to these slide_ids")
    p.add_argument("--set", nargs="*", dest="overrides", default=[])
    p.add_argument("--dry-run", action="store_true", help="report the plan, write nothing")
    return p.parse_args()


def main():
    args = parse_args()
    cfg = load_config(args.config, args.overrides)
    set_seed(int(cfg.get("seed", 1337)))
    spec = ClassSpec.from_config(cfg)

    manifest_path = args.manifest or cfg.get("manifest")
    if not manifest_path:
        raise SystemExit("No manifest: pass --manifest or set `manifest:` in the config.")
    rows = load_manifest(resolve_path(cfg, manifest_path))
    if args.slides:
        rows = [r for r in rows if r["slide_id"] in set(args.slides)]
        if not rows:
            raise SystemExit(f"No manifest rows match --slides {args.slides}")
    print(f"Manifest: {manifest_summary(rows)}")

    out_dir = Path(args.out or resolve_path(cfg, cfg["paths"]["patches_dir"]))
    pcfg = cfg["patch"]
    size, stride = int(pcfg["size"]), int(pcfg["stride"])
    rng = np.random.default_rng(int(cfg.get("seed", 1337)))

    summary = {}
    for row in rows:
        summary[row["slide_id"]] = extract_slide(row, cfg, spec, out_dir, size, stride,
                                                 rng, args.dry_run)

    total = sum(s["n_patches"] for s in summary.values())
    print(f"\n{'=' * 72}\nTotal patches: {total}  ->  {out_dir}")
    counts = Counter()
    for s in summary.values():
        counts.update(s["class_counts"])
    print("Instances by class: " + json.dumps(dict(counts.most_common())))
    if not args.dry_run:
        out_dir.mkdir(parents=True, exist_ok=True)
        with open(out_dir / "extraction_summary.json", "w") as f:
            json.dump({"summary": summary, "class_counts": dict(counts),
                       "classes": list(spec.names), "patch_size": size,
                       "stride": stride}, f, indent=2)


def extract_slide(row, cfg, spec: ClassSpec, out_dir, size, stride, rng, dry_run) -> dict:
    slide_id = row["slide_id"]
    wcfg = cfg["wsi"]
    # Per-slide overrides. The cohort is not homogeneous: the Nx .ndpi slides carry the
    # 0.44 um/px plane at level 1, while Ctrl_1A2.tif is a single-level export that is
    # already 0.44 um/px at level 0 and whose only resolution tag is a 300-DPI print
    # artefact. A manifest `level` / `manual_mpp` column pins each slide to the plane its
    # annotations were drawn on; blank falls back to the config.
    level = int(row["level"]) if str(row.get("level", "")).strip() else wcfg.get("level", 1)
    mmpp = (float(row["manual_mpp"]) if str(row.get("manual_mpp", "")).strip()
            else wcfg.get("manual_mpp"))
    reader = WSIReader(
        row["image_path"],
        level=level,
        target_mpp=float(wcfg.get("target_mpp", 0.44068)),
        manual_mpp=mmpp,
        anisotropy_tol=float(wcfg.get("anisotropy_tol", 0.02)),
        prefer=wcfg.get("prefer", "finer"),
    )
    W, H = reader.level_dimensions_work
    pcfg = cfg["patch"]

    raw = load_geojson_instances(
        row["geojson_path"],
        each_part_separate=(pcfg.get("multipolygon") == "each_part_separate"),
        multipolygon_min_area_frac=float(pcfg.get("multipolygon_min_area_frac", 0.05)),
    )

    # Coordinate space: 'level' annotations are already at the working level; 'full' ones
    # are level-0 QuPath coordinates and must be scaled down.
    to_work = 1.0 if pcfg.get("annotations_space", "full") == "level" \
        else reader.scale_full_to_work
    polys, names, roles = [], [], []
    for inst in raw:
        role = spec.role(inst.class_name)
        if role == "ignore":
            continue
        g = inst.polygon if to_work == 1.0 else shapely_scale(
            inst.polygon, xfact=to_work, yfact=to_work, origin=(0, 0)
        )
        polys.append(g)
        names.append(inst.class_name)
        roles.append(role)

    inst_polys = [p for p, r in zip(polys, roles) if r == "instance"]
    inst_names = [n for n, r in zip(names, roles) if r == "instance"]
    bg_polys = [p for p, r in zip(polys, roles) if r == "background"]
    roi_polys = [p for p, r in zip(polys, roles) if r == "dense_roi"]

    print(f"\n[{slide_id}] {W}x{H} @ {reader.working_mpp:.5f} µm/px  "
          f"instances={len(inst_polys)} bg_polys={len(bg_polys)} rois={len(roi_polys)}")
    if not inst_polys:
        print(f"[{slide_id}] no instance annotations after class mapping — skipped")
        return {"n_patches": 0, "class_counts": {}}

    # Restrict tiling to the annotated area: outside it nothing is labelled, so mining
    # background there would train the model on structures nobody drew.
    xs = [p.bounds for p in inst_polys]
    if pcfg.get("restrict_to_annotations", True):
        m = int(pcfg.get("annotation_margin", 512))
        x_min = max(0, int(min(b[0] for b in xs)) - m)
        y_min = max(0, int(min(b[1] for b in xs)) - m)
        x_max = min(W, int(max(b[2] for b in xs)) + m)
        y_max = min(H, int(max(b[3] for b in xs)) + m)
    else:
        x_min, y_min, x_max, y_max = 0, 0, W, H

    # Optional manifest `roi` column: one or more "x0,y0,x1,y1" rectangles in WORKING
    # coordinates, separated by ";". This is how ONE slide contributes to two splits --
    # list it twice with complementary rectangle sets.
    #
    # A patch is kept only if it lies ENTIRELY inside one rectangle. That containment rule
    # is what creates the buffer: a patch straddling a rectangle border fails the test for
    # both splits and is simply dropped, so no validation patch can share a pixel with a
    # training patch even at stride 256. Nothing else enforces that.
    #
    # Use INTERLEAVED BLOCKS, not one contiguous band. Kidney anatomy is spatially
    # organised: a band split of Ctrl_1A2 put 85% of the validation objects in
    # tubule_collecting because the medulla sits at one end of the tissue, which measures
    # the medulla and nothing else. Interleaving keeps both splits anatomically
    # representative. The cost is honest and worth stating: adjacent blocks share staining
    # and section artefacts, so this band is NOT an independent test of cross-slide
    # generalisation -- Nx_4wk_050 remains the fully held-out slide for that.
    rects = parse_rois(str(row.get("roi", "") or ""))
    if rects:
        rx0 = min(r[0] for r in rects); ry0 = min(r[1] for r in rects)
        rx1 = max(r[2] for r in rects); ry1 = max(r[3] for r in rects)
        x_min, y_min = max(x_min, rx0), max(y_min, ry0)
        x_max, y_max = min(x_max, rx1), min(y_max, ry1)
        if x_max <= x_min or y_max <= y_min:
            print(f"[{slide_id}] roi does not intersect the annotated area — skipped")
            return {"n_patches": 0, "class_counts": {}}
        print(f"[{slide_id}] {len(rects)} roi rect(s) -> ({x_min},{y_min})-({x_max},{y_max})")
    region_w, region_h = x_max - x_min, y_max - y_min
    print(f"[{slide_id}] tiling region ({x_min},{y_min}) {region_w}x{region_h}")

    min_fg = float(pcfg.get("min_fg_fraction", 0.02))
    bg_keep = float(pcfg.get("bg_keep_ratio", 0.0))
    bg_glass_keep = float(pcfg.get("bg_glass_keep_ratio", bg_keep))
    bg_min_glass = float(pcfg.get("bg_min_glass_frac", 0.20))
    class_counts = Counter()
    n_written = n_bg = n_bg_glass = 0

    for tx, ty in iter_tiles(region_w, region_h, size, stride):
        x, y = x_min + tx, y_min + ty
        if x + size > W or y + size > H:
            continue
        if not patch_in_rois(x, y, size, rects):
            continue
        window = box(x, y, x + size, y + size)
        hits = [i for i, p in enumerate(inst_polys) if p.intersects(window)]

        id_map = np.zeros((size, size), np.int32)
        class_idx_of_id: dict[int, int] = {}
        subtype_of_id: dict[int, str] = {}
        if hits:
            sub_polys = [inst_polys[i] for i in hits]
            sub_names = [inst_names[i] for i in hits]
            id_map, name_of_id = rasterize_instances(sub_polys, sub_names, size, size,
                                                     offset_xy=(x, y))
            for iid, raw_name in name_of_id.items():
                class_idx_of_id[iid] = spec.canonical_index(raw_name)
                subtype_of_id[iid] = raw_name

        frac = float((id_map > 0).mean())
        rgb = reader.read_region_work(x, y, size, size)
        tissue = build_tissue_mask(rgb, cfg.get("tissue"))

        if frac < min_fg:
            # Keep glassy background patches preferentially: off-tissue pixels are trusted
            # background, whereas an unannotated all-tissue tile is almost entirely ignored.
            glass = 1.0 - float(tissue.mean())
            keep_p = bg_glass_keep if glass >= bg_min_glass else bg_keep
            if keep_p <= 0 or rng.random() > keep_p:
                continue
            n_bg += 1
            n_bg_glass += glass >= bg_min_glass
        bg_mask = _mask_from_polys(bg_polys, x, y, size)
        roi_mask = _mask_from_polys(roi_polys, x, y, size)

        for ci in class_idx_of_id.values():
            class_counts[spec.names[ci]] += 1

        if not dry_run:
            save_patch(out_dir, slide_id, x, y, rgb, id_map, class_idx_of_id,
                       subtype_of_id=subtype_of_id, tissue=tissue, bg_poly=bg_mask,
                       dense_roi=roi_mask, class_names=list(spec.names),
                       extra_meta={"working_mpp": reader.working_mpp,
                                   "scale_to_full": reader.scale_to_full,
                                   "split": row.get("split", "train")})
        n_written += 1
        if n_written % 50 == 0:
            print(f"[{slide_id}]   {n_written} patches...", flush=True)

    reader.close()
    print(f"[{slide_id}] wrote {n_written} patches "
          f"({n_bg} background-only, {n_bg_glass} of them glassy) | {dict(class_counts)}")
    return {"n_patches": n_written, "n_background": n_bg, "n_background_glass": int(n_bg_glass),
            "class_counts": dict(class_counts)}


def _mask_from_polys(polys, x, y, size) -> np.ndarray:
    if not polys:
        return np.zeros((size, size), np.uint8)
    from shapely.geometry import box as _box

    window = _box(x, y, x + size, y + size)
    hits = [p for p in polys if p.intersects(window)]
    if not hits:
        return np.zeros((size, size), np.uint8)
    m, _ = rasterize_instances(hits, ["_"] * len(hits), size, size, offset_xy=(x, y))
    return (m > 0).astype(np.uint8)


if __name__ == "__main__":
    sys.exit(main())
