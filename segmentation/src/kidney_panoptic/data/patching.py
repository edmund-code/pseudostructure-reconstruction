"""Patch tiling and the on-disk patch format.

On-disk layout (per patch, under ``patches_dir/<slide_id>/``)::

    <pid>.png    RGB at working resolution
    <pid>.npz    id_map (int32), tissue (uint8), bg_poly (uint8), dense_roi (uint8)
    <pid>.json   metadata: slide_id, origin_xy, size, class_idx_of_id, subtype_of_id,
                 per-class instance counts, fg_fraction, tissue_fraction

Only the instance map and the three auxiliary masks are stored. Semantic, boundary and
center targets are REGENERATED from ``id_map`` after augmentation (see data/targets.py),
which is what keeps geometric augmentation correct by construction.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Iterator, Optional

import cv2
import numpy as np


Rect = tuple[int, int, int, int]


def parse_rois(spec: str) -> list[Rect]:
    """Parse a manifest ``roi`` cell: ``"x0,y0,x1,y1;x0,y0,x1,y1;..."`` -> list of rects.

    Empty / blank returns ``[]``, meaning "no restriction".
    """
    rects: list[Rect] = []
    for part in (spec or "").split(";"):
        part = part.strip()
        if not part:
            continue
        vals = [int(float(v)) for v in part.split(",")]
        if len(vals) != 4:
            raise ValueError(f"roi rect needs 4 numbers 'x0,y0,x1,y1', got {part!r}")
        x0, y0, x1, y1 = vals
        if x1 <= x0 or y1 <= y0:
            raise ValueError(f"roi rect {part!r} is empty or inverted")
        rects.append((x0, y0, x1, y1))
    return rects


def patch_in_rois(x: int, y: int, size: int, rects: list[Rect]) -> bool:
    """True if the patch [x,x+size) x [y,y+size) lies ENTIRELY inside some rect.

    Full containment, not intersection, and that is the whole point. When one slide feeds
    two splits via complementary rect sets, a patch straddling the interface fails this
    test for BOTH sets and is dropped — which is what keeps a validation patch from sharing
    pixels with a training patch. Relax this to `intersects` and the two splits silently
    overlap by up to one patch width, inflating validation scores with memorised pixels.
    """
    if not rects:
        return True
    return any(x >= a and y >= b and x + size <= c and y + size <= d
               for a, b, c, d in rects)


def iter_tiles(width: int, height: int, size: int, stride: int) -> Iterator[tuple[int, int]]:
    """Yield (x, y) top-left origins covering [0,width) x [0,height); last row/col clamp."""
    def axis(extent: int) -> list[int]:
        last = max(0, extent - size)
        vals = list(range(0, last + 1, stride))
        if not vals or vals[-1] != last:
            vals.append(last)
        return sorted(set(vals))

    for y in axis(height):
        for x in axis(width):
            yield x, y


def fg_fraction(id_map: np.ndarray) -> float:
    return float((id_map > 0).mean())


def patch_id(slide_id: str, x: int, y: int) -> str:
    return f"{slide_id}_x{x}_y{y}"


def save_patch(
    out_dir: str | Path,
    slide_id: str,
    x: int,
    y: int,
    image: np.ndarray,
    id_map: np.ndarray,
    class_idx_of_id: dict[int, int],
    *,
    subtype_of_id: Optional[dict[int, str]] = None,
    tissue: Optional[np.ndarray] = None,
    bg_poly: Optional[np.ndarray] = None,
    dense_roi: Optional[np.ndarray] = None,
    class_names: Optional[list[str]] = None,
    extra_meta: Optional[dict] = None,
) -> str:
    out_dir = Path(out_dir) / slide_id
    out_dir.mkdir(parents=True, exist_ok=True)
    pid = patch_id(slide_id, x, y)
    h, w = id_map.shape

    cv2.imwrite(str(out_dir / f"{pid}.png"), cv2.cvtColor(image, cv2.COLOR_RGB2BGR))

    def _u8(m):
        return (np.zeros((h, w), np.uint8) if m is None else m.astype(np.uint8))

    np.savez_compressed(
        out_dir / f"{pid}.npz",
        id_map=id_map.astype(np.int32),
        tissue=_u8(tissue),
        bg_poly=_u8(bg_poly),
        dense_roi=_u8(dense_roi),
    )

    counts: dict[str, int] = {}
    if class_names:
        for ci in class_idx_of_id.values():
            counts[class_names[ci]] = counts.get(class_names[ci], 0) + 1

    meta = {
        "patch_id": pid,
        "slide_id": slide_id,
        "origin_xy": [int(x), int(y)],
        "size": [int(w), int(h)],
        "class_idx_of_id": {str(k): int(v) for k, v in class_idx_of_id.items()},
        "subtype_of_id": {str(k): v for k, v in (subtype_of_id or {}).items()},
        "class_counts": counts,
        "n_instances": len(class_idx_of_id),
        "fg_fraction": fg_fraction(id_map),
        "tissue_fraction": float(_u8(tissue).mean()) if tissue is not None else None,
        **(extra_meta or {}),
    }
    with open(out_dir / f"{pid}.json", "w") as f:
        json.dump(meta, f, indent=2)
    return pid


def load_patch(meta_path: str | Path) -> dict:
    """Load a saved patch into image / id_map / masks / meta."""
    meta_path = Path(meta_path)
    with open(meta_path) as f:
        meta = json.load(f)
    base = str(meta_path.with_suffix(""))
    bgr = cv2.imread(base + ".png", cv2.IMREAD_COLOR)
    if bgr is None:
        raise FileNotFoundError(f"Missing patch image for {meta_path}")
    arrs = np.load(base + ".npz")
    return {
        "image": cv2.cvtColor(bgr, cv2.COLOR_BGR2RGB),
        "id_map": arrs["id_map"],
        "tissue": arrs["tissue"].astype(bool),
        "bg_poly": arrs["bg_poly"].astype(bool),
        "dense_roi": arrs["dense_roi"].astype(bool),
        "class_idx_of_id": {int(k): int(v) for k, v in meta["class_idx_of_id"].items()},
        "meta": meta,
    }


def find_patch_metas(patches_dir: str | Path,
                     slide_ids: Optional[list[str]] = None) -> list[Path]:
    patches_dir = Path(patches_dir)
    if not patches_dir.exists():
        raise FileNotFoundError(
            f"Patch directory not found: {patches_dir}. Run scripts/extract_patches.py first."
        )
    subdirs = ([patches_dir / s for s in slide_ids] if slide_ids
               else sorted(p for p in patches_dir.iterdir() if p.is_dir()))
    metas: list[Path] = []
    for d in subdirs:
        if d.is_dir():
            metas.extend(sorted(d.glob("*.json")))
    return metas
