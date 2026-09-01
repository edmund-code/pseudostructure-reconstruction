"""ROI rectangles: the mechanism that lets ONE slide feed two splits.

The failure this guards against is silent and flattering: if a validation patch shares
pixels with a training patch, validation PQ goes up and nothing looks wrong. At stride 256
with 512 px patches, an off-by-one in the containment rule leaks half a patch.
"""
import pytest

from kidney_panoptic.data.patching import parse_rois, patch_in_rois


def test_blank_roi_means_no_restriction():
    assert parse_rois("") == []
    assert parse_rois("   ") == []
    assert patch_in_rois(0, 0, 512, []) is True


def test_parses_multiple_rects():
    assert parse_rois("0,0,100,200;300,400,500,600") == [(0, 0, 100, 200),
                                                         (300, 400, 500, 600)]
    assert parse_rois("0,0,10,10;") == [(0, 0, 10, 10)]


@pytest.mark.parametrize("bad", ["1,2,3", "1,2,3,4,5", "10,0,5,20", "0,0,0,10"])
def test_malformed_rects_raise(bad):
    """Silently skipping a bad rect would silently change the split."""
    with pytest.raises(ValueError):
        parse_rois(bad)


def test_patch_must_be_fully_inside_not_merely_overlapping():
    rects = [(0, 0, 1024, 1024)]
    assert patch_in_rois(0, 0, 512, rects)
    assert patch_in_rois(512, 512, 512, rects)      # exactly flush with the far corner
    assert not patch_in_rois(513, 0, 512, rects)    # one pixel over the right edge
    assert not patch_in_rois(0, 513, 512, rects)    # one pixel over the bottom edge
    assert not patch_in_rois(-1, 0, 512, rects)


def test_a_patch_may_satisfy_any_one_rect():
    rects = [(0, 0, 600, 600), (2000, 2000, 2600, 2600)]
    assert patch_in_rois(2000, 2000, 512, rects)
    assert not patch_in_rois(1000, 1000, 512, rects)


def test_complementary_rect_sets_cannot_both_accept_a_patch():
    """The real invariant. Two adjacent blocks assigned to different splits must never
    both admit the same patch — otherwise train and val overlap."""
    train = [(0, 0, 2048, 2048)]
    val = [(2048, 0, 4096, 2048)]
    size, stride = 512, 256
    shared = [(x, y) for x in range(0, 4096 - size + 1, stride)
              for y in range(0, 2048 - size + 1, stride)
              if patch_in_rois(x, y, size, train) and patch_in_rois(x, y, size, val)]
    assert shared == [], f"{len(shared)} patches accepted by both splits"


def test_patches_straddling_the_interface_are_dropped_by_both_splits():
    """Containment implies a buffer: an origin whose patch crosses the block border is
    rejected by both rect sets. The invariant worth asserting is not the buffer's width but
    its consequence — the last training patch must END at or before the first val patch
    BEGINS, so the two splits share no pixel."""
    train = [(0, 0, 2048, 2048)]
    val = [(2048, 0, 4096, 2048)]
    size, stride = 512, 256
    origins = range(0, 4096 - size + 1, stride)
    dropped = [x for x in origins if not patch_in_rois(x, 0, size, train)
               and not patch_in_rois(x, 0, size, val)]
    assert dropped, "no buffer at all — every origin landed in one split or the other"

    last_train = max(x for x in origins if patch_in_rois(x, 0, size, train))
    first_val = min(x for x in origins if patch_in_rois(x, 0, size, val))
    assert last_train + size <= first_val, "train and val patches share pixels"
