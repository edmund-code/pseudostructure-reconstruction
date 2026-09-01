"""GeoJSON export must be loadable by QuPath.

Two properties are load-or-fail, both learned from a file QuPath refused:
  * every Polygon has exactly ONE ring (QuPath does not import interior rings), and
  * every geometry is VALID (JTS throws on self-intersections, which approxPolyDP makes).
"""
import json

import numpy as np
import pytest
from shapely.geometry import shape

from kidney_panoptic.postprocess.geojson_out import instance_polygons, write_geojson

NAMES = ["background", "tubule", "glomerulus", "vessel", "blood_cell", "nerve"]


def ring_with_hole():
    """An annulus — a tubule cross-section with a lumen."""
    lbl = np.zeros((80, 80), np.int32)
    yy, xx = np.ogrid[:80, :80]
    r2 = (xx - 40) ** 2 + (yy - 40) ** 2
    lbl[(r2 < 30 ** 2) & (r2 > 12 ** 2)] = 1
    return lbl


def spiky_blob(seed=0):
    """A jagged shape whose contour simplification is prone to self-intersection."""
    rng = np.random.default_rng(seed)
    lbl = np.zeros((120, 120), np.int32)
    yy, xx = np.ogrid[:120, :120]
    ang = np.arctan2(yy - 60, xx - 60)
    rad = 35 + 14 * np.sin(7 * ang) + rng.normal(0, 2.5, (120, 120))
    lbl[((xx - 60) ** 2 + (yy - 60) ** 2) < rad ** 2] = 1
    return lbl


def features_for(lbl, **kw):
    return list(instance_polygons(lbl, {1: 1}, NAMES, score_of_id={1: 0.9}, **kw))


def test_polygons_have_exactly_one_ring_even_for_an_annulus():
    feats = features_for(ring_with_hole())
    assert feats
    for f in feats:
        assert len(f["geometry"]["coordinates"]) == 1, "QuPath cannot import interior rings"


def test_every_emitted_geometry_is_valid():
    for seed in range(12):
        for f in features_for(spiky_blob(seed)):
            g = shape(f["geometry"])
            assert g.is_valid, f"invalid geometry from seed {seed}"
            assert not g.is_empty


def test_geometry_stays_valid_under_aggressive_simplification():
    """The vertex cap is where approxPolyDP used to introduce self-intersections."""
    for seed in range(8):
        for f in features_for(spiky_blob(seed), simplify_px=2.0, max_vertices=12):
            g = shape(f["geometry"])
            assert g.is_valid
            assert len(f["geometry"]["coordinates"][0]) <= 14   # cap + closing vertex


def test_rings_are_closed_and_have_at_least_four_points():
    for f in features_for(spiky_blob(3)):
        ring = f["geometry"]["coordinates"][0]
        assert ring[0] == ring[-1], "GeoJSON rings must be explicitly closed"
        assert len(ring) >= 4


def test_feature_schema_matches_what_qupath_expects():
    f = features_for(spiky_blob(1))[0]
    assert f["type"] == "Feature"
    assert isinstance(f["id"], str)
    p = f["properties"]
    assert p["objectType"] == "annotation"
    assert p["classification"]["name"] == "tubule"
    assert isinstance(p["classification"]["colorRGB"], int)
    assert p["isLocked"] is False
    assert all(isinstance(m["value"], float) for m in p["measurements"])
    assert {m["name"] for m in p["measurements"]} >= {"Area px", "Score"}


def test_ids_are_sequential_strings_from_zero():
    lbl = np.zeros((60, 160), np.int32)
    lbl[10:50, 10:50] = 1
    lbl[10:50, 60:100] = 2
    lbl[10:50, 110:150] = 3
    feats = list(instance_polygons(lbl, {1: 1, 2: 2, 3: 3}, NAMES))
    assert [f["id"] for f in feats] == ["0", "1", "2"]
    assert [f["properties"]["instance_id"] for f in feats] == [1, 2, 3]


def test_coordinates_are_offset_and_scaled():
    lbl = np.zeros((40, 40), np.int32)
    lbl[10:30, 10:30] = 1
    f = list(instance_polygons(lbl, {1: 1}, NAMES, offset_xy=(100, 200), scale=2.0))[0]
    xs = [c[0] for c in f["geometry"]["coordinates"][0]]
    ys = [c[1] for c in f["geometry"]["coordinates"][0]]
    assert min(xs) == pytest.approx((100 + 10) * 2, abs=4)
    assert min(ys) == pytest.approx((200 + 10) * 2, abs=4)


def test_written_file_parses_as_a_feature_collection(tmp_path):
    path = tmp_path / "out.geojson"
    n = write_geojson(features_for(spiky_blob(5)), path)
    d = json.loads(path.read_text())
    assert d["type"] == "FeatureCollection"
    assert len(d["features"]) == n >= 1
    for f in d["features"]:
        assert shape(f["geometry"]).is_valid


def test_empty_label_map_writes_a_valid_empty_collection(tmp_path):
    path = tmp_path / "empty.geojson"
    n = write_geojson(instance_polygons(np.zeros((20, 20), np.int32), {}, NAMES), path)
    assert n == 0
    assert json.loads(path.read_text())["features"] == []
