import numpy as np
import pytest

from kidney_panoptic.data.ignore import (
    STRATEGIES, apply_ignore_to_semantic, build_valid_mask, trusted_bg_fraction,
)


@pytest.fixture
def scene():
    """40x80: one annotated instance on tissue (left), glass on the right."""
    id_map = np.zeros((40, 80), np.int32)
    id_map[10:30, 5:25] = 1
    tissue = np.zeros((40, 80), bool)
    tissue[:, :40] = True
    return id_map, tissue


def test_plus_tissue_makes_glass_trusted_background(scene):
    id_map, tissue = scene
    v = build_valid_mask(id_map, strategy="dilated_gt_two_radius_plus_tissue",
                         tissue_mask=tissue, dilate_px_inner=5)
    assert v[0, 70], "off-tissue glass must be trusted background"
    assert v[20, 15], "the annotated instance must be valid"
    assert not v[35, 35], "unannotated TISSUE must stay ignored"


def test_two_radius_without_tissue_leaves_glass_ignored(scene):
    id_map, tissue = scene
    v = build_valid_mask(id_map, strategy="dilated_gt_two_radius",
                         tissue_mask=tissue, dilate_px_inner=5)
    assert not v[0, 70], "without the tissue term, glass is ignored"
    assert v[20, 15]


def test_plus_tissue_requires_a_tissue_mask(scene):
    id_map, _ = scene
    with pytest.raises(ValueError, match="requires tissue_mask"):
        build_valid_mask(id_map, strategy="dilated_gt_two_radius_plus_tissue")


def test_ring_width_controls_trusted_background_around_instances(scene):
    id_map, tissue = scene
    thin = build_valid_mask(id_map, strategy="dilated_gt_two_radius",
                            tissue_mask=tissue, dilate_px_inner=3)
    thick = build_valid_mask(id_map, strategy="dilated_gt_two_radius",
                             tissue_mask=tissue, dilate_px_inner=12)
    assert trusted_bg_fraction(thick, id_map) > trusted_bg_fraction(thin, id_map)


def test_instances_always_win_over_every_other_term(scene):
    """An annotation on a tissue-mask miss must not be relabelled background."""
    id_map, tissue = scene
    tissue[:] = False                      # pathological: tissue detector saw nothing
    v = build_valid_mask(id_map, strategy="dilated_gt_two_radius_plus_tissue",
                         tissue_mask=tissue)
    assert v[(id_map > 0)].all()


def test_annotated_background_polygons_become_trusted(scene):
    id_map, tissue = scene
    bg = np.zeros_like(tissue)
    bg[32:38, 30:38] = True                # a drawn "Background" region on tissue
    v = build_valid_mask(id_map, strategy="dilated_gt_two_radius_plus_tissue",
                         tissue_mask=tissue, background_mask=bg)
    assert v[35, 34]


def test_dense_roi_makes_everything_inside_valid(scene):
    id_map, tissue = scene
    roi = np.zeros_like(tissue)
    roi[:, 30:40] = True
    v = build_valid_mask(id_map, strategy="dense_roi", dense_roi_mask=roi)
    assert v[20, 35] and not v[20, 60]
    assert v[(id_map > 0)].all()


def test_none_strategy_is_all_valid(scene):
    id_map, _ = scene
    assert build_valid_mask(id_map, strategy="none").all()


def test_unknown_strategy_fails_loudly(scene):
    with pytest.raises(ValueError, match="Unknown ignore strategy"):
        build_valid_mask(scene[0], strategy="whatever")


def test_all_documented_strategies_are_reachable(scene):
    id_map, tissue = scene
    roi = np.ones_like(tissue)
    for s in STRATEGIES:
        v = build_valid_mask(id_map, strategy=s, tissue_mask=tissue, dense_roi_mask=roi)
        assert v.shape == id_map.shape and v.dtype == bool


def test_apply_ignore_stamps_255(scene):
    id_map, tissue = scene
    v = build_valid_mask(id_map, strategy="dilated_gt_two_radius", tissue_mask=tissue)
    sem = np.ones_like(id_map, dtype=np.int64)
    out = apply_ignore_to_semantic(sem, v)
    assert (out[~v] == 255).all() and (out[v] == 1).all()
