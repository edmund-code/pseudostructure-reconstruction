"""Decode tests driven by SYNTHETIC PERFECT PREDICTIONS.

Each test builds the probability maps a perfectly trained model would emit for a known
scene, then checks the decode recovers the scene. That isolates the decode from model
quality: a failure here is a decode bug, never a training problem.
"""
import numpy as np
import pytest
from scipy import ndimage as ndi

from kidney_panoptic.data.targets import BOUNDARY, INTERIOR, boundary_target, center_heatmaps
from kidney_panoptic.postprocess.decode import decode_panoptic
from kidney_panoptic.eval.instance_metrics import panoptic_quality


def perfect_maps(id_map, class_idx_of_id, spec, boundary_width=3, sharpness=1.0):
    """The probability maps an ideal model would produce for this scene."""
    bnd, _ = boundary_target(id_map, boundary_width)
    C = spec.n_classes
    H, W = id_map.shape

    boundary_prob = np.full((3, H, W), (1 - sharpness) / 3, np.float32)
    for c in (0, 1, 2):
        boundary_prob[c][bnd == c] += sharpness
    boundary_prob /= boundary_prob.sum(0, keepdims=True)

    sem = np.zeros((H, W), np.int64)
    for i, c in class_idx_of_id.items():
        sem[id_map == i] = c
    semantic_prob = np.full((C, H, W), (1 - sharpness) / C, np.float32)
    for c in range(C):
        semantic_prob[c][sem == c] += sharpness
    semantic_prob /= semantic_prob.sum(0, keepdims=True)

    centers = center_heatmaps(id_map, class_idx_of_id, spec.compact_indices)
    return semantic_prob, boundary_prob, centers.astype(np.float32)


DECODE_CFG = {
    "decode": {"fg_threshold": 0.5, "interior_threshold": 0.5, "min_seed_area": 10,
               "min_area": {"default": 20}, "center_threshold": 0.35,
               "center_min_distance": 8},
    "infer": {"tissue_gate": False},
}


def test_touching_objects_are_recovered_as_two(touching_pair, spec, cfg):
    sem, bnd, cen = perfect_maps(touching_pair, {1: 1, 2: 1}, spec)
    res = decode_panoptic(sem, bnd, cen, DECODE_CFG, spec)
    assert res.n_instances == 2, "a zero-gap pair must decode to two instances"
    assert panoptic_quality(touching_pair, res.id_map)["PQ"] > 0.75


def test_class_vote_assigns_the_right_class(spec):
    lbl = np.zeros((64, 96), np.int32)
    lbl[12:52, 8:48] = 1
    lbl[12:52, 56:88] = 2
    cmap = {1: 1, 2: 4}          # tubule_proximal, glomerulus
    sem, bnd, cen = perfect_maps(lbl, cmap, spec)
    res = decode_panoptic(sem, bnd, cen, DECODE_CFG, spec)
    got = sorted(res.class_of_id.values())
    assert got == [1, 4], f"expected tubule_proximal+glomerulus, got {got}"


def test_empty_input_yields_no_instances(spec):
    H = W = 48
    bnd = np.zeros((3, H, W), np.float32); bnd[0] = 1.0
    sem = np.zeros((spec.n_classes, H, W), np.float32); sem[0] = 1.0
    cen = np.zeros((len(spec.compact_indices), H, W), np.float32)
    res = decode_panoptic(sem, bnd, cen, DECODE_CFG, spec)
    assert res.n_instances == 0 and res.id_map.max() == 0


def test_tissue_gate_removes_off_tissue_instances(spec):
    lbl = np.zeros((40, 80), np.int32)
    lbl[10:30, 5:35] = 1        # on tissue
    lbl[10:30, 45:75] = 2       # on glass
    sem, bnd, cen = perfect_maps(lbl, {1: 1, 2: 1}, spec)
    tissue = np.zeros((40, 80), bool); tissue[:, :40] = True
    cfg = {**DECODE_CFG, "infer": {"tissue_gate": True, "tissue_min_instance_frac": 0.5}}
    res = decode_panoptic(sem, bnd, cen, cfg, spec, tissue_mask=tissue)
    assert res.n_instances == 1


def test_min_area_is_applied_per_class(spec):
    lbl = np.zeros((60, 60), np.int32)
    lbl[5:11, 5:11] = 1         # 36 px "blood cell"
    lbl[20:50, 20:50] = 2       # 900 px "tubule"
    cmap = {1: 4, 2: 1}   # glomerulus (tiny), tubule_proximal (large)
    sem, bnd, cen = perfect_maps(lbl, cmap, spec)
    cfg = {"decode": {**DECODE_CFG["decode"],
                      "min_area": {"default": 10, "glomerulus": 1000, "tubule_proximal": 10}},
           "infer": {"tissue_gate": False}}
    res = decode_panoptic(sem, bnd, cen, cfg, spec)
    assert list(res.class_of_id.values()) == [1], "the small glomerulus must be dropped, not the tubule"


def test_center_peaks_split_two_glomeruli_with_a_weak_boundary(spec):
    """The compact path exists for exactly this: no boundary evidence, two centers."""
    lbl = np.zeros((80, 140), np.int32)
    yy, xx = np.ogrid[:80, :140]
    lbl[((xx - 40) ** 2 + (yy - 40) ** 2) < 30 ** 2] = 1
    lbl[((xx - 98) ** 2 + (yy - 40) ** 2) < 30 ** 2] = 2

    sem, bnd, cen = perfect_maps(lbl, {1: 4, 2: 4}, spec)   # 4 = glomerulus
    # Erase the boundary evidence between them: the only remaining cue is the centers.
    seam = slice(60, 80)
    bnd[2][:, seam] = 0.02
    bnd[1][:, seam] = np.where(lbl[:, seam] > 0, 0.95, 0.02)
    bnd[0][:, seam] = np.where(lbl[:, seam] > 0, 0.03, 0.96)
    bnd /= bnd.sum(0, keepdims=True)

    cfg = {"decode": {**DECODE_CFG["decode"], "center_min_distance": 20,
                      "min_area": {"default": 50}},
           "infer": {"tissue_gate": False}}
    res = decode_panoptic(sem, bnd, cen, cfg, spec)
    assert res.n_instances == 2, f"centers should split the fused blob, got {res.n_instances}"


def test_decode_is_translation_equivariant(touching_pair, spec):
    """The decode reads only local fields, so shifting the scene shifts the output."""
    sem, bnd, cen = perfect_maps(touching_pair, {1: 1, 2: 1}, spec)
    pad = ((0, 0), (7, 0), (11, 0))
    shifted = decode_panoptic(np.pad(sem, pad, constant_values=0),
                              np.pad(bnd, pad, constant_values=0),
                              np.pad(cen, pad, constant_values=0), DECODE_CFG, spec)
    base = decode_panoptic(sem, bnd, cen, DECODE_CFG, spec)
    assert shifted.n_instances == base.n_instances


def test_decode_on_a_window_matches_decode_on_the_whole(spec):
    """Chunking cannot manufacture or destroy an instance away from the cut."""
    lbl = np.zeros((64, 256), np.int32)
    for k in range(4):
        lbl[12:52, 8 + k * 60: 56 + k * 60] = k + 1
    sem, bnd, cen = perfect_maps(lbl, {k + 1: 1 for k in range(4)}, spec)

    whole = decode_panoptic(sem, bnd, cen, DECODE_CFG, spec)
    half = decode_panoptic(sem[:, :, :128], bnd[:, :, :128], cen[:, :, :128],
                           DECODE_CFG, spec)
    # The first two objects lie entirely within the first half.
    assert whole.n_instances == 4
    assert half.n_instances >= 2


def test_merge_fragments_rejects_parallel_glue(spec):
    """Two side-by-side tubules must NOT be fused even with a weak interface."""
    lbl = np.zeros((80, 60), np.int32)
    lbl[10:70, 8:26] = 1
    lbl[10:70, 30:48] = 2
    sem, bnd, cen = perfect_maps(lbl, {1: 1, 2: 1}, spec)

    # Make them adjacent with almost no boundary signal between them.
    lbl2 = lbl.copy(); lbl2[10:70, 26:30] = 1
    sem, bnd, cen = perfect_maps(lbl2, {1: 1, 2: 1}, spec)
    bnd[2] *= 0.05
    bnd /= bnd.sum(0, keepdims=True)

    cfg = {"decode": {**DECODE_CFG["decode"], "merge_fragments": True,
                      "merge_max_boundary_prob": 0.9,
                      "merge_max_minor_axis_growth": 1.15,
                      "min_area": {"default": 50}},
           "infer": {"tissue_gate": False}}
    res = decode_panoptic(sem, bnd, cen, cfg, spec)
    assert res.n_instances >= 2, "minor-axis growth must veto parallel glue"


def test_scores_are_probabilities(touching_pair, spec):
    sem, bnd, cen = perfect_maps(touching_pair, {1: 1, 2: 1}, spec)
    res = decode_panoptic(sem, bnd, cen, DECODE_CFG, spec)
    assert all(0.0 <= s <= 1.0 for s in res.score_of_id.values())
    assert set(res.class_of_id) == set(res.score_of_id)


def test_ids_are_contiguous_from_one(touching_pair, spec):
    sem, bnd, cen = perfect_maps(touching_pair, {1: 1, 2: 1}, spec)
    res = decode_panoptic(sem, bnd, cen, DECODE_CFG, spec)
    ids = sorted(i for i in np.unique(res.id_map) if i)
    assert ids == list(range(1, len(ids) + 1))
    assert set(ids) == set(res.class_of_id)


def test_relabel_drops_instances_with_no_class_rather_than_calling_them_background(spec):
    """A classless instance must not be exported as a background polygon."""
    from kidney_panoptic.postprocess.decode import relabel

    labels = np.zeros((10, 10), np.int32)
    labels[1:4, 1:4] = 1
    labels[6:9, 6:9] = 2          # deliberately absent from the class map
    out, cls, score = relabel(labels, {1: 3}, {1: 0.5})
    assert set(cls.values()) == {3}
    assert sorted(np.unique(out)) == [0, 1]
    assert 0 not in cls.values(), "background is never a valid instance class"


# --------------------------------------------------------------- semantic split
SPLIT_CFG = {
    "decode": {**DECODE_CFG["decode"],
               "min_area": {"default": 200},
               "semantic_split": {"enabled": True, "min_class_prob": 0.50, "open_px": 2,
                                  "min_component_area": 300, "min_minor_frac": 0.15,
                                  "max_parts": 4, "elevation": "boundary_plus_ambiguity",
                                  "ambiguity_weight": 1.0, "pairs": None}},
    "infer": {"tissue_gate": False},
}


def _two_segment_blob(spec, erase_wall=True):
    """One connected blob whose left half is tubule_proximal and right half is tubule_distal."""
    lbl = np.zeros((80, 160), np.int32)
    lbl[15:65, 10:150] = 1                       # a SINGLE instance in the GT sense
    sem, bnd, cen = perfect_maps(lbl, {1: 1}, spec)
    # Repaint semantics: left half tubule_proximal (1), right half tubule_distal (2).
    sem[:] = 0.02
    sem[0][lbl == 0] = 0.96
    sem[1][:, :80][lbl[:, :80] > 0] = 0.95
    sem[2][:, 80:][lbl[:, 80:] > 0] = 0.95
    sem /= sem.sum(0, keepdims=True)
    if erase_wall:
        bnd[2][:, 70:90] = 0.02                  # boundary head "missed" the wall
        bnd /= bnd.sum(0, keepdims=True)
    return sem, bnd, cen


def test_semantic_split_separates_two_segments_with_no_boundary_evidence(spec):
    """The mechanism's entire purpose: split where the boundary head failed."""
    sem, bnd, cen = _two_segment_blob(spec)
    assert decode_panoptic(sem, bnd, cen, DECODE_CFG, spec).n_instances == 1
    res = decode_panoptic(sem, bnd, cen, SPLIT_CFG, spec)
    assert res.n_instances == 2, f"expected a split into 2, got {res.n_instances}"
    assert sorted(res.class_of_id.values()) == [1, 2]


def test_semantic_split_does_not_punch_holes(spec):
    """Every pixel of the parent must survive in exactly one child."""
    sem, bnd, cen = _two_segment_blob(spec)
    before = decode_panoptic(sem, bnd, cen, DECODE_CFG, spec).id_map > 0
    after = decode_panoptic(sem, bnd, cen, SPLIT_CFG, spec).id_map > 0
    assert np.array_equal(before, after), "splitting must repartition, never delete pixels"


def test_semantic_split_ignores_speckle(spec):
    """A few stray pixels of another class must not manufacture an instance."""
    lbl = np.zeros((80, 160), np.int32)
    lbl[15:65, 10:150] = 1
    sem, bnd, cen = perfect_maps(lbl, {1: 1}, spec)
    rng = np.random.default_rng(0)
    ys, xs = np.nonzero(lbl > 0)
    pick = rng.permutation(len(ys))[: int(0.03 * len(ys))]
    sem[2][ys[pick], xs[pick]] = 0.95
    sem[1][ys[pick], xs[pick]] = 0.03
    sem /= sem.sum(0, keepdims=True)
    assert decode_panoptic(sem, bnd, cen, SPLIT_CFG, spec).n_instances == 1


def test_semantic_split_respects_min_minor_frac(spec):
    """A tiny minority region is not a second object."""
    lbl = np.zeros((80, 160), np.int32)
    lbl[15:65, 10:150] = 1
    sem, bnd, cen = perfect_maps(lbl, {1: 1}, spec)
    sem[2][:, 140:][lbl[:, 140:] > 0] = 0.95     # ~7% of the instance
    sem[1][:, 140:][lbl[:, 140:] > 0] = 0.03
    sem /= sem.sum(0, keepdims=True)
    cfg = {**SPLIT_CFG, "decode": {**SPLIT_CFG["decode"],
                                   "semantic_split": {**SPLIT_CFG["decode"]["semantic_split"],
                                                      "min_minor_frac": 0.15}}}
    assert decode_panoptic(sem, bnd, cen, cfg, spec).n_instances == 1


def test_semantic_split_respects_the_pair_allow_list(spec):
    sem, bnd, cen = _two_segment_blob(spec)
    cfg = {**SPLIT_CFG, "decode": {**SPLIT_CFG["decode"],
                                   "semantic_split": {**SPLIT_CFG["decode"]["semantic_split"],
                                                      "pairs": [["tubule_proximal", "tubule_collecting"]]}}}
    assert decode_panoptic(sem, bnd, cen, cfg, spec).n_instances == 1


def test_semantic_split_does_not_fire_on_a_single_class_blob(spec, touching_pair):
    """No regression on the ordinary case: same-class neighbours are unaffected."""
    sem, bnd, cen = perfect_maps(touching_pair, {1: 1, 2: 1}, spec)
    assert decode_panoptic(sem, bnd, cen, SPLIT_CFG, spec).n_instances == \
           decode_panoptic(sem, bnd, cen, DECODE_CFG, spec).n_instances


def test_semantic_split_is_translation_equivariant(spec):
    sem, bnd, cen = _two_segment_blob(spec)
    pad = ((0, 0), (9, 0), (13, 0))
    base = decode_panoptic(sem, bnd, cen, SPLIT_CFG, spec)
    moved = decode_panoptic(np.pad(sem, pad), np.pad(bnd, pad), np.pad(cen, pad),
                            SPLIT_CFG, spec)
    assert moved.n_instances == base.n_instances


def test_semantic_split_rejects_an_unknown_elevation_mode(spec):
    sem, bnd, cen = _two_segment_blob(spec)
    cfg = {**SPLIT_CFG, "decode": {**SPLIT_CFG["decode"],
                                   "semantic_split": {**SPLIT_CFG["decode"]["semantic_split"],
                                                      "elevation": "nonsense"}}}
    with pytest.raises(ValueError, match="elevation"):
        decode_panoptic(sem, bnd, cen, cfg, spec)


def test_decode_accepts_a_four_channel_boundary_head(spec, touching_pair):
    """bg / interior / boundary-to-bg / boundary-to-instance: both walls are ridges."""
    sem, bnd3, cen = perfect_maps(touching_pair, {1: 1, 2: 1}, spec)
    from kidney_panoptic.data.targets import boundary_target
    b4, _ = boundary_target(touching_pair, 3, n_classes=4)
    bnd4 = np.zeros((4,) + touching_pair.shape, np.float32)
    for c in range(4):
        bnd4[c][b4 == c] = 1.0
    res = decode_panoptic(sem, bnd4, cen, DECODE_CFG, spec)
    assert res.n_instances == 2
