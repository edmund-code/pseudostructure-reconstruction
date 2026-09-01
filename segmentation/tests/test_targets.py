import numpy as np
import pytest

from kidney_panoptic.data.targets import (
    BG, BOUNDARY, INTERIOR, boundary_target, build_targets, center_heatmaps,
    semantic_target,
)


def test_semantic_maps_ids_to_classes():
    lbl = np.zeros((10, 10), np.int32)
    lbl[2:5, 2:5] = 1
    lbl[6:9, 6:9] = 2
    sem = semantic_target(lbl, {1: 1, 2: 3})
    assert sem[3, 3] == 1 and sem[7, 7] == 3 and sem[0, 0] == 0
    assert sem.dtype == np.int64


def test_semantic_unlisted_id_falls_back_to_background():
    lbl = np.zeros((6, 6), np.int32)
    lbl[1:3, 1:3] = 7
    assert semantic_target(lbl, {}).max() == 0


def test_boundary_three_classes_and_interior_survives(touching_pair):
    bnd, touch = boundary_target(touching_pair, boundary_width_px=3)
    assert set(np.unique(bnd)) <= {BG, INTERIOR, BOUNDARY}
    assert (bnd == INTERIOR).any(), "objects must keep an interior at width 3"
    # background stays background — the band is carved out of the objects, not grown
    assert bnd[0, 0] == BG
    assert (bnd[touching_pair == 0] == BG).all()


def test_touching_pixels_are_only_at_the_instance_interface(touching_pair):
    _, touch = boundary_target(touching_pair, boundary_width_px=3)
    assert touch.any(), "a zero-gap pair must produce touching pixels"
    xs = np.nonzero(touch)[1]
    assert xs.min() >= 44 and xs.max() <= 51, "touching band must hug the x=48 seam"
    # An isolated object has no touching pixels at all.
    solo = np.zeros((40, 40), np.int32)
    solo[10:30, 10:30] = 1
    assert not boundary_target(solo, 3)[1].any()


def test_boundary_band_separates_touching_interiors(touching_pair):
    """The point of the representation: interiors of two touching objects disconnect."""
    from scipy import ndimage as ndi

    bnd, _ = boundary_target(touching_pair, boundary_width_px=3)
    n = ndi.label(bnd == INTERIOR)[1]
    assert n == 2, f"expected 2 disconnected interiors, got {n}"


def test_wider_band_eats_more_interior(touching_pair):
    a = (boundary_target(touching_pair, 2)[0] == INTERIOR).sum()
    b = (boundary_target(touching_pair, 5)[0] == INTERIOR).sum()
    assert b < a


def test_center_heatmap_peaks_at_centroid_of_compact_classes():
    lbl = np.zeros((60, 60), np.int32)
    lbl[10:30, 10:30] = 1     # glomerulus (compact, class 2 -> channel 0)
    lbl[40:50, 40:50] = 2     # tubule (elongated, no channel)
    heat = center_heatmaps(lbl, {1: 2, 2: 1}, compact_indices=[2, 4])
    assert heat.shape == (2, 60, 60)
    peak = np.unravel_index(np.argmax(heat[0]), heat[0].shape)
    assert abs(peak[0] - 19.5) <= 1 and abs(peak[1] - 19.5) <= 1
    # The centroid is at (19.5, 19.5) — between pixels — so the discrete maximum sits
    # just under 1.0. That is the sub-pixel centroid being encoded honestly, not a bug.
    assert 0.95 <= heat[0].max() <= 1.0
    assert heat[1].max() == 0.0, "no blood cells present -> empty channel"


def test_center_sigma_scales_with_object_size():
    small = np.zeros((80, 80), np.int32); small[38:42, 38:42] = 1
    big = np.zeros((80, 80), np.int32); big[20:60, 20:60] = 1
    hs = center_heatmaps(small, {1: 2}, [2])[0]
    hb = center_heatmaps(big, {1: 2}, [2])[0]
    assert (hb > 0.5).sum() > (hs > 0.5).sum()


def test_targets_regenerate_identically_under_flip(touching_pair):
    """Flip-invariance by construction — nothing in the targets carries orientation."""
    cmap = {1: 1, 2: 1}
    a = build_targets(touching_pair, cmap, [2, 4])
    b = build_targets(np.ascontiguousarray(touching_pair[:, ::-1]), cmap, [2, 4])
    assert np.array_equal(a.boundary[:, ::-1], b.boundary)
    assert np.array_equal(a.semantic[:, ::-1], b.semantic)
    assert np.allclose(a.centers[:, :, ::-1], b.centers)


def test_empty_input_is_all_background():
    lbl = np.zeros((16, 16), np.int32)
    t = build_targets(lbl, {}, [2, 4])
    assert (t.semantic == 0).all() and (t.boundary == BG).all()
    assert t.centers.shape == (2, 16, 16) and t.centers.max() == 0


# ------------------------------------------------- 4-class boundary (ADR 0006)
def test_four_class_boundary_promotes_the_instance_instance_wall(touching_pair):
    from kidney_panoptic.data.targets import TOUCHING

    b3, touch = boundary_target(touching_pair, 3, n_classes=3)
    b4, touch4 = boundary_target(touching_pair, 3, n_classes=4)

    assert set(np.unique(b3)) <= {BG, INTERIOR, BOUNDARY}
    assert TOUCHING in np.unique(b4), "the wall between two objects must get its own class"
    assert np.array_equal(touch, touch4)
    # The promoted pixels are exactly the touching mask, and nothing else moved.
    assert np.array_equal(b4 == TOUCHING, touch)
    assert np.array_equal(b4 == INTERIOR, b3 == INTERIOR)
    assert np.array_equal(b4 == BG, b3 == BG)
    assert np.array_equal((b4 == BOUNDARY) | (b4 == TOUCHING), b3 == BOUNDARY)


def test_four_class_boundary_on_an_isolated_object_has_no_touching_pixels():
    from kidney_panoptic.data.targets import TOUCHING

    solo = np.zeros((40, 40), np.int32)
    solo[10:30, 10:30] = 1
    b4, _ = boundary_target(solo, 3, n_classes=4)
    assert TOUCHING not in np.unique(b4)
    assert BOUNDARY in np.unique(b4), "object-to-background contour is still class 2"


def test_boundary_target_rejects_an_unsupported_class_count(touching_pair):
    with pytest.raises(ValueError, match="3 or 4"):
        boundary_target(touching_pair, 3, n_classes=5)


def test_build_targets_threads_boundary_classes(touching_pair):
    from kidney_panoptic.data.targets import TOUCHING, targets_from_config

    t3 = targets_from_config(touching_pair, {1: 1, 2: 1}, [4], {"model": {"boundary_classes": 3}})
    t4 = targets_from_config(touching_pair, {1: 1, 2: 1}, [4], {"model": {"boundary_classes": 4}})
    assert TOUCHING not in np.unique(t3.boundary)
    assert TOUCHING in np.unique(t4.boundary)
