"""QuPath GeoJSON loading (fail-loud class resolution) and writing.

Loading
-------
Parses a FeatureCollection of Polygon/MultiPolygon features. Each feature is ONE
instance. The class name is resolved by trying, in order:
    properties.classification.name -> properties.classification (str)
    -> properties.name -> properties.class
If none are present, we raise listing the actual property keys found (fail loud).

MultiPolygon default (plan R5): the LARGEST part is the instance; smaller parts are
dropped below ``multipolygon_min_area_frac`` of the largest. Set
``each_part_separate=True`` to instead emit one instance per part.

Writing
-------
QuPath-compatible: exterior-only polygons (holes stripped), geometry repaired with
shapely.make_valid (fallback buffer(0)), adaptive simplification to a vertex cap,
coords rounded, unique string ids, classification.name + colorRGB, streamed to disk.

(Vendored unchanged from the v1 `kidneyseg` package, which validated it across the
NDPI cohort; kept self-contained here so kidney_panoptic has no cross-tree import.)
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Optional

import numpy as np
from shapely.geometry import MultiPolygon, Polygon, shape as shapely_shape
from shapely.validation import make_valid


@dataclass
class Instance:
    """One loaded annotation instance in full-resolution coordinates."""
    polygon: Polygon
    class_name: str
    raw_properties: dict = field(default_factory=dict)


def _resolve_class_name(props: dict) -> str:
    classification = props.get("classification")
    if isinstance(classification, dict) and classification.get("name"):
        return str(classification["name"])
    if isinstance(classification, str) and classification:
        return classification
    if props.get("name"):
        return str(props["name"])
    if props.get("class"):
        return str(props["class"])
    raise KeyError(
        "Could not resolve class name from feature properties. Tried "
        "classification.name, classification, name, class. "
        f"Actual property keys present: {sorted(props.keys())}"
    )


def _geom_to_polygons(geom) -> list[Polygon]:
    if geom.geom_type == "Polygon":
        return [geom]
    if geom.geom_type == "MultiPolygon":
        return list(geom.geoms)
    return []


def load_geojson_instances(
    path: str | Path,
    *,
    each_part_separate: bool = False,
    multipolygon_min_area_frac: float = 0.05,
    ignore_names: Optional[Iterable[str]] = None,
) -> list[Instance]:
    """Load instances from a QuPath GeoJSON FeatureCollection (full-res coords)."""
    path = Path(path)
    with open(path) as f:
        data = json.load(f)

    feats = data.get("features", data if isinstance(data, list) else [])
    ignore_set = set(ignore_names or [])
    instances: list[Instance] = []

    for feat in feats:
        geom_dict = feat.get("geometry")
        if not geom_dict:
            continue
        props = feat.get("properties", {}) or {}
        class_name = _resolve_class_name(props)
        if class_name in ignore_set:
            continue

        geom = shapely_shape(geom_dict)
        if not geom.is_valid:
            geom = make_valid(geom)
        parts = [p for p in _geom_to_polygons(geom) if p.area > 0]
        if not parts:
            continue

        if each_part_separate:
            for p in parts:
                instances.append(Instance(_clean_polygon(p), class_name, props))
        else:
            largest = max(parts, key=lambda p: p.area)
            # Keep only parts above the area-fraction threshold; merge them.
            thresh = largest.area * multipolygon_min_area_frac
            kept = [p for p in parts if p.area >= thresh]
            merged = kept[0]
            for p in kept[1:]:
                merged = merged.union(p)
            if merged.geom_type == "MultiPolygon":
                merged = max(merged.geoms, key=lambda p: p.area)
            instances.append(Instance(_clean_polygon(merged), class_name, props))

    return instances


def _clean_polygon(poly: Polygon) -> Polygon:
    """Strip holes (QuPath limitation) and repair geometry."""
    if poly.geom_type != "Polygon":
        poly = max(_geom_to_polygons(poly), key=lambda p: p.area)
    exterior = Polygon(poly.exterior)
    if not exterior.is_valid:
        repaired = make_valid(exterior)
        if repaired.geom_type == "MultiPolygon":
            repaired = max(repaired.geoms, key=lambda p: p.area)
        if repaired.geom_type == "Polygon":
            exterior = repaired
        else:
            exterior = exterior.buffer(0)
    return exterior


# --------------------------------------------------------------------- writing

# QuPath signed-int ARGB colors (a few defaults; classes beyond this cycle a palette).
_PALETTE = [-65536, -16711936, -16776961, -256, -65281, -16711681, -23296, -8388608, -16744448, -10496]


def class_color(class_name: str, class_list: Optional[list[str]] = None) -> int:
    if class_list and class_name in class_list:
        return _PALETTE[class_list.index(class_name) % len(_PALETTE)]
    return _PALETTE[abs(hash(class_name)) % len(_PALETTE)]


def repair_geometry(poly: Polygon) -> Optional[Polygon]:
    """make_valid first, buffer(0) fallback; return a single exterior-only Polygon."""
    if poly.is_valid and poly.geom_type == "Polygon":
        g = poly
    else:
        g = make_valid(poly)
        if g.geom_type == "MultiPolygon":
            g = max(g.geoms, key=lambda p: p.area)
        if g.geom_type != "Polygon":
            g = poly.buffer(0)
            if g.geom_type == "MultiPolygon":
                g = max(g.geoms, key=lambda p: p.area)
    if g.is_empty or g.geom_type != "Polygon":
        return None
    return Polygon(g.exterior)


def polygon_to_feature(
    poly: Polygon,
    class_name: str,
    feature_id: str,
    *,
    class_list: Optional[list[str]] = None,
    max_vertices: int = 500,
    round_decimals: int = 2,
    measurements: Optional[list[dict]] = None,
) -> Optional[dict]:
    """Convert a polygon (full-res coords) to a QuPath Feature dict, or None if invalid."""
    g = repair_geometry(poly)
    if g is None:
        return None
    # Adaptive simplification to cap vertex count.
    tol = 0.5
    while len(g.exterior.coords) > max_vertices and tol < 100:
        g = g.simplify(tol, preserve_topology=True)
        if g.geom_type == "MultiPolygon":
            g = max(g.geoms, key=lambda p: p.area)
        if g.is_empty:
            return None
        tol *= 2
    coords = [[round(float(x), round_decimals), round(float(y), round_decimals)]
              for x, y in g.exterior.coords]
    if len(coords) < 4:
        return None
    return {
        "type": "Feature",
        "id": feature_id,
        "geometry": {"type": "Polygon", "coordinates": [coords]},
        "properties": {
            "objectType": "annotation",
            "classification": {"name": class_name, "colorRGB": class_color(class_name, class_list)},
            "isLocked": False,
            "measurements": measurements or [],
        },
    }


def write_geojson_features(features: Iterable[dict], path: str | Path) -> int:
    """Stream features to a FeatureCollection file. Returns count written."""
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
