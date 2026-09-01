"""Rasterize instance polygons into an instance-ID label map + per-instance class lookup.

Given polygons already expressed in a patch's local pixel frame, produce:
* ``id_map``     : int32 array, 0 = background, 1..N = instance ids.
* ``class_of_id``: dict instance_id -> class_name.

Touching instances keep distinct ids (no merging), which is what the HoVer targets
and watershed evaluation rely on. Painting is done largest-area-first so that when a
smaller instance is drawn it overwrites overlap from a larger one (the smaller, more
specific object wins its pixels), mirroring typical annotation intent.

(Vendored unchanged from the v1 `kidneyseg` package, which validated it across the
NDPI cohort; kept self-contained here so kidney_panoptic has no cross-tree import.)
"""
from __future__ import annotations

from typing import Iterable

import numpy as np
from shapely.affinity import translate
from shapely.geometry import Polygon

try:
    from rasterio.features import rasterize as _rio_rasterize
    _HAVE_RIO = True
except Exception:  # pragma: no cover
    _HAVE_RIO = False

import cv2


def rasterize_instances(
    polygons: list[Polygon],
    class_names: list[str],
    height: int,
    width: int,
    offset_xy: tuple[float, float] = (0.0, 0.0),
) -> tuple[np.ndarray, dict[int, str]]:
    """Rasterize polygons (in working-res coords) into an instance-id map.

    Parameters
    ----------
    polygons, class_names : aligned lists; one class per polygon.
    height, width         : output patch size (working-res pixels).
    offset_xy             : subtract this (patch origin) from polygon coords first.
    """
    id_map = np.zeros((height, width), dtype=np.int32)
    class_of_id: dict[int, str] = {}
    ox, oy = offset_xy

    # Largest first so smaller instances overwrite overlap.
    order = sorted(range(len(polygons)), key=lambda i: polygons[i].area, reverse=True)
    next_id = 1
    for idx in order:
        poly = translate(polygons[idx], xoff=-ox, yoff=-oy)
        if poly.is_empty:
            continue
        mask = _rasterize_one(poly, height, width)
        if mask.sum() == 0:
            continue
        id_map[mask] = next_id
        class_of_id[next_id] = class_names[idx]
        next_id += 1

    return id_map, class_of_id


def _rasterize_one(poly: Polygon, height: int, width: int) -> np.ndarray:
    if _HAVE_RIO:
        arr = _rio_rasterize(
            [(poly, 1)], out_shape=(height, width), fill=0,
            all_touched=False, dtype="uint8",
        )
        return arr.astype(bool)
    # OpenCV fallback.
    mask = np.zeros((height, width), np.uint8)
    exterior = np.array(poly.exterior.coords, dtype=np.int32)
    cv2.fillPoly(mask, [exterior], 1)
    for interior in poly.interiors:
        cv2.fillPoly(mask, [np.array(interior.coords, dtype=np.int32)], 0)
    return mask.astype(bool)


def relabel_sequential(id_map: np.ndarray) -> np.ndarray:
    """Relabel an arbitrary-id map to 0..N contiguous ids (background stays 0)."""
    out = np.zeros_like(id_map)
    for new_id, old_id in enumerate(sorted(i for i in np.unique(id_map) if i != 0), start=1):
        out[id_map == old_id] = new_id
    return out
