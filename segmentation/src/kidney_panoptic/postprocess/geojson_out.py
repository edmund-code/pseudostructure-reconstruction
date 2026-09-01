"""Instance label map -> QuPath-compatible GeoJSON polygons.

Contours are traced per instance inside its bounding box (never over the whole slide), so
memory stays bounded regardless of slide size.

Two properties are non-negotiable if QuPath is to open the file at all, and both were
learned the hard way:

**Exterior rings only.** QuPath's GeoJSON import does not round-trip polygons with interior
rings; a file with holes fails to load rather than loading without them. Holes are therefore
stripped by default (``keep_holes=False``). The cost is real: a tubule lumen is a genuine
hole, so the drawn POLYGON encloses more than the object does. The ``Area px`` measurement
is NOT affected — it is computed from the instance mask, which still excludes the lumen — so
the number attached to each object stays correct while its outline is the filled exterior.
Anyone measuring area from polygon geometry rather than from the ``Area px`` measurement (or
from ``instance_map.npz``) will overcount by the lumen.

**Valid geometry.** Contour simplification with ``cv2.approxPolyDP`` readily produces
self-intersecting rings — measured at 5.4 % of instances on a real slide — and QuPath's JTS
geometry engine throws on those. Every ring is therefore repaired with ``shapely.make_valid``
(``buffer(0)`` fallback) and simplified with shapely's topology-preserving simplifier, which
cannot introduce a self-intersection, rather than with ``approxPolyDP``.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable, Iterator, Optional

import cv2
import numpy as np
from scipy import ndimage as ndi
from shapely.geometry import Polygon
from shapely.validation import make_valid

# QuPath signed-int ARGB colours.
_PALETTE = [-65536, -16711936, -16776961, -256, -65281, -16711681,
            -23296, -8388608, -16744448, -10496]


def class_color(class_name: str, class_list: Optional[list[str]] = None) -> int:
    if class_list and class_name in class_list:
        return _PALETTE[class_list.index(class_name) % len(_PALETTE)]
    return _PALETTE[abs(hash(class_name)) % len(_PALETTE)]


def instance_polygons(
    labels: np.ndarray,
    class_of_id: dict[int, int],
    class_names: list[str],
    *,
    score_of_id: Optional[dict[int, float]] = None,
    offset_xy: tuple[float, float] = (0.0, 0.0),
    scale: float = 1.0,
    keep_holes: bool = False,
    simplify_px: float = 1.0,
    max_vertices: int = 400,
    mpp: Optional[float] = None,
) -> Iterator[dict]:
    """Yield one GeoJSON Feature per instance, each a valid single-ring Polygon."""
    ox, oy = offset_xy
    objects = ndi.find_objects(labels)
    emitted = 0
    for idx, sl in enumerate(objects, start=1):
        if sl is None or idx not in class_of_id:
            continue
        mask = (labels[sl] == idx).astype(np.uint8)
        if not mask.any():
            continue
        area_px = int(mask.sum())
        y0, x0 = sl[0].start, sl[1].start

        # Pad by 1 so contours of instances touching the bbox edge close properly.
        padded = np.pad(mask, 1)
        contours, _ = cv2.findContours(padded, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        if not contours:
            continue
        ring = max(contours, key=cv2.contourArea).reshape(-1, 2)
        if len(ring) < 4:
            continue

        # Move to output coordinates BEFORE repair, so simplification tolerances are in
        # the same units as the emitted geometry.
        xy = ring.astype(np.float64)
        xy[:, 0] = (xy[:, 0] - 1 + x0 + ox) * scale
        xy[:, 1] = (xy[:, 1] - 1 + y0 + oy) * scale

        poly = _repair(xy)
        if poly is None:
            continue
        poly = _simplify_to_cap(poly, simplify_px, max_vertices)
        if poly is None:
            continue

        coords = [[round(float(a), 2), round(float(b), 2)] for a, b in poly.exterior.coords]
        if coords[0] != coords[-1]:
            coords.append(coords[0])
        if len(coords) < 4:
            continue

        name = class_names[class_of_id[idx]]
        measurements = [{"name": "Area px", "value": float(area_px)}]
        if mpp:
            measurements.append({"name": "Area um^2",
                                 "value": round(float(area_px * mpp * mpp), 2)})
        if score_of_id and idx in score_of_id:
            measurements.append({"name": "Score", "value": round(float(score_of_id[idx]), 4)})

        yield {
            "type": "Feature",
            "id": str(emitted),
            "geometry": {"type": "Polygon", "coordinates": [coords]},
            "properties": {
                "objectType": "annotation",
                "classification": {"name": name, "colorRGB": class_color(name, class_names)},
                "isLocked": False,
                "measurements": measurements,
                "instance_id": int(idx),
            },
        }
        emitted += 1


def _repair(xy: np.ndarray) -> Optional[Polygon]:
    """Exterior-only, valid, non-empty Polygon from a raw contour ring.

    ``make_valid`` on a self-intersecting ring returns a MultiPolygon (or a
    GeometryCollection); the largest polygonal piece is the instance and the rest are
    slivers thrown off by the self-intersection.
    """
    try:
        poly = Polygon(xy)
    except Exception:
        return None
    if not poly.is_valid:
        repaired = make_valid(poly)
        poly = _largest_polygon(repaired)
        if poly is None:
            poly = Polygon(xy).buffer(0)
            poly = _largest_polygon(poly)
    if poly is None or poly.is_empty:
        return None
    poly = Polygon(poly.exterior)          # strip holes — QuPath cannot import them
    if not poly.is_valid:
        poly = _largest_polygon(make_valid(poly))
    if poly is None or poly.is_empty or poly.area <= 0:
        return None
    return poly


def _largest_polygon(geom) -> Optional[Polygon]:
    if geom is None or geom.is_empty:
        return None
    if geom.geom_type == "Polygon":
        return geom
    parts = [g for g in getattr(geom, "geoms", []) if g.geom_type == "Polygon" and not g.is_empty]
    return max(parts, key=lambda g: g.area) if parts else None


def _simplify_to_cap(poly: Polygon, tol: float, max_vertices: int) -> Optional[Polygon]:
    """Topology-preserving simplification, escalating until the vertex cap is met.

    shapely's simplifier preserves validity, unlike ``cv2.approxPolyDP``, which is what
    was producing self-intersecting rings that QuPath refused to load.
    """
    tol = max(float(tol), 0.0)
    if tol > 0:
        poly = _largest_polygon(poly.simplify(tol, preserve_topology=True)) or poly
    guard = max(tol, 0.5)
    while len(poly.exterior.coords) > max_vertices and guard < 512:
        guard *= 2
        simpler = _largest_polygon(poly.simplify(guard, preserve_topology=True))
        if simpler is None or simpler.is_empty:
            break
        poly = simpler
    if not poly.is_valid:
        poly = _largest_polygon(make_valid(poly))
    if poly is None or poly.is_empty or len(poly.exterior.coords) < 4:
        return None
    return Polygon(poly.exterior)


def write_geojson(features: Iterable[dict], path: str | Path) -> int:
    """Stream features into a FeatureCollection. Returns the count written."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(path, "w", encoding="utf-8") as f:
        f.write('{"type":"FeatureCollection","features":[\n')
        for feat in features:
            if feat is None:
                continue
            f.write(("," if n else "") + json.dumps(feat))
            n += 1
        f.write("\n]}")
    return n
