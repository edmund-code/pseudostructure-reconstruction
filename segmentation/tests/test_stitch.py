import numpy as np
import pytest

from kidney_panoptic.infer.tiled import decode_windows, gaussian_window, tile_origins
from kidney_panoptic.postprocess.stitch import InstanceStitcher, UnionFind


def test_union_find_collapses_chains_to_the_smallest_id():
    uf = UnionFind()
    uf.union(3, 7); uf.union(7, 11); uf.union(20, 21)
    assert uf.find(11) == uf.find(3) == 3
    assert uf.find(21) == 20
    assert uf.find(3) != uf.find(20)


def _bar(h, w, y0, y1, x0, x1, label=1):
    a = np.zeros((h, w), np.int32)
    a[y0:y1, x0:x1] = label
    return a


def test_object_split_across_two_windows_becomes_one_instance():
    """A bar spanning the seam: each window sees part of it, output must be one object."""
    s = InstanceStitcher(64, 200, stitch_iou=0.5)
    # window A covers x 0..128, window B covers x 100..200 (28 px overlap)
    s.add_window(0, 0, _bar(64, 128, 20, 40, 60, 128), {1: 1}, {1: 0.9})
    s.add_window(100, 0, _bar(64, 100, 20, 40, 0, 50), {1: 1}, {1: 0.9})
    labels, cls, score = s.finalize()
    ids = [i for i in np.unique(labels) if i]
    assert len(ids) == 1, f"the bar must stitch into one instance, got {len(ids)}"
    assert labels[30, 70] == labels[30, 140]
    assert cls[ids[0]] == 1


def test_distinct_objects_in_the_overlap_are_not_fused():
    s = InstanceStitcher(64, 200, stitch_iou=0.5)
    a = np.zeros((64, 128), np.int32); a[5:15, 100:120] = 1
    b = np.zeros((64, 100), np.int32); b[40:55, 10:30] = 1     # far from a
    s.add_window(0, 0, a, {1: 1})
    s.add_window(100, 0, b, {1: 2})
    labels, cls, _ = s.finalize()
    assert len([i for i in np.unique(labels) if i]) == 2
    assert sorted(cls.values()) == [1, 2]


def test_three_windows_crossing_one_object_collapse_to_one_id():
    s = InstanceStitcher(64, 300, stitch_iou=0.5)
    s.add_window(0, 0, _bar(64, 128, 20, 40, 40, 128), {1: 3})
    s.add_window(100, 0, _bar(64, 128, 20, 40, 0, 128), {1: 3})
    s.add_window(200, 0, _bar(64, 100, 20, 40, 0, 60), {1: 3})
    labels, cls, _ = s.finalize()
    assert len([i for i in np.unique(labels) if i]) == 1
    assert labels[30, 50] == labels[30, 150] == labels[30, 250]


def test_class_of_a_merged_instance_comes_from_the_largest_fragment():
    s = InstanceStitcher(64, 200, stitch_iou=0.5)
    s.add_window(0, 0, _bar(64, 128, 20, 40, 20, 128), {1: 1})     # big fragment, tubule
    s.add_window(100, 0, _bar(64, 100, 20, 40, 0, 40), {1: 3})     # small, vessel
    labels, cls, _ = s.finalize()
    ids = [i for i in np.unique(labels) if i]
    assert len(ids) == 1 and cls[ids[0]] == 1


def test_finalize_produces_contiguous_ids_matching_the_class_map():
    s = InstanceStitcher(64, 200)
    a = np.zeros((64, 100), np.int32)
    a[5:15, 5:15] = 1; a[30:40, 30:40] = 2
    s.add_window(0, 0, a, {1: 1, 2: 2})
    labels, cls, score = s.finalize()
    ids = sorted(i for i in np.unique(labels) if i)
    assert ids == list(range(1, len(ids) + 1))
    assert set(ids) == set(cls) == set(score)


def test_empty_stitcher_finalizes_cleanly():
    labels, cls, score = InstanceStitcher(32, 32).finalize()
    assert labels.max() == 0 and cls == {} and score == {}


def test_gaussian_window_is_strictly_positive_and_peaks_at_the_centre():
    w = gaussian_window(64)
    assert (w > 0).all(), "a zero weight would leave a pixel unnormalised"
    assert np.unravel_index(np.argmax(w), w.shape) in [(31, 31), (32, 32), (31, 32), (32, 31)]


def test_tile_origins_cover_the_requested_span_on_the_global_grid():
    origins = tile_origins(1000, 1500, size=512, step=384, limit=4000)
    assert all(o % 384 == 0 or o == 4000 - 512 for o in origins)
    # every pixel of [1000, 1500) must be covered by some tile
    covered = np.zeros(1500 - 1000, bool)
    for o in origins:
        lo, hi = max(o, 1000), min(o + 512, 1500)
        if hi > lo:
            covered[lo - 1000:hi - 1000] = True
    assert covered.all()


def test_decode_windows_tile_the_slide_completely():
    W, H = 900, 700
    covered = np.zeros((H, W), bool)
    for x, y, w, h in decode_windows(W, H, 512, 128):
        covered[y:y + h, x:x + w] = True
        assert x + w <= W and y + h <= H
    assert covered.all()


def test_output_scale_follows_the_annotation_space():
    """Predictions must land on the same coordinates as the annotations they trained on."""
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    from wsi_infer import output_scale

    class FakeReader:
        scale_to_full = 2

    level_cfg = {"infer": {"output_space": "auto"}, "patch": {"annotations_space": "level"}}
    full_cfg = {"infer": {"output_space": "auto"}, "patch": {"annotations_space": "full"}}
    assert output_scale(level_cfg, FakeReader()) == 1.0
    assert output_scale(full_cfg, FakeReader()) == 2.0
    assert output_scale({"infer": {"output_space": "full"},
                         "patch": {"annotations_space": "level"}}, FakeReader()) == 2.0
    with pytest.raises(ValueError, match="output_space"):
        output_scale({"infer": {"output_space": "nonsense"}, "patch": {}}, FakeReader())
