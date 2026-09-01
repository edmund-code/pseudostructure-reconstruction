import numpy as np
import pytest

from kidney_panoptic.eval.instance_metrics import (
    aggregated_jaccard_index, evaluate_instances, f1_at_iou, filter_by_valid,
    panoptic_quality, split_merge_rates,
)
from kidney_panoptic.eval.semantic_metrics import semantic_scores


def _two_boxes():
    lbl = np.zeros((40, 80), np.int32)
    lbl[10:30, 5:35] = 1
    lbl[10:30, 45:75] = 2
    return lbl


def test_pq_perfect_match_is_one():
    a = _two_boxes()
    m = panoptic_quality(a, a.copy())
    assert m["PQ"] == pytest.approx(1.0)
    assert m["SQ"] == pytest.approx(1.0) and m["DQ"] == pytest.approx(1.0)
    assert (m["TP"], m["FP"], m["FN"]) == (2, 0, 0)


def test_pq_dq_alias_matches():
    a = _two_boxes()
    m = panoptic_quality(a, a.copy())
    assert m["DQ"] == m["RQ"]


def test_pq_empty_vs_empty_is_one_and_empty_vs_full_is_zero():
    z = np.zeros((20, 20), np.int32)
    assert panoptic_quality(z, z)["PQ"] == 1.0
    a = _two_boxes()
    assert panoptic_quality(a, np.zeros_like(a))["PQ"] == 0.0
    assert panoptic_quality(np.zeros_like(a), a)["PQ"] == 0.0


def test_pq_missed_object_halves_dq():
    gt = _two_boxes()
    pred = gt.copy()
    pred[pred == 2] = 0
    m = panoptic_quality(gt, pred)
    assert (m["TP"], m["FP"], m["FN"]) == (1, 0, 1)
    assert m["DQ"] == pytest.approx(1 / 1.5)


def test_aji_perfect_and_disjoint():
    a = _two_boxes()
    assert aggregated_jaccard_index(a, a.copy()) == pytest.approx(1.0)
    b = np.zeros_like(a)
    b[35:39, 0:4] = 1
    assert aggregated_jaccard_index(a, b) < 0.05


def test_aji_uses_its_own_matching_not_pqs():
    """A prediction below the PQ threshold still contributes to AJI."""
    gt = np.zeros((20, 20), np.int32); gt[2:18, 2:18] = 1
    pred = np.zeros((20, 20), np.int32); pred[2:8, 2:18] = 1     # ~37% IoU
    assert panoptic_quality(gt, pred)["TP"] == 0
    assert aggregated_jaccard_index(gt, pred) > 0.3


def test_split_detected_at_low_overlap_where_pq_reports_only_a_miss():
    gt = np.zeros((40, 40), np.int32)
    gt[5:35, 5:35] = 1
    pred = np.zeros((40, 40), np.int32)
    pred[5:19, 5:35] = 1        # two halves, neither above IoU 0.5
    pred[21:35, 5:35] = 2
    assert panoptic_quality(gt, pred)["TP"] == 0
    sm = split_merge_rates(gt, pred, overlap_frac=0.20)
    assert sm["n_split"] == 1 and sm["split_rate"] == pytest.approx(1.0)


def test_merge_rate_counts_fused_gt_objects():
    gt = _two_boxes()
    pred = np.zeros_like(gt)
    pred[10:30, 5:75] = 1       # one prediction swallowing both
    sm = split_merge_rates(gt, pred, overlap_frac=0.20)
    assert sm["n_merged_gt"] == 2 and sm["merge_rate"] == pytest.approx(1.0)
    assert sm["n_split"] == 0


def test_f1_thresholds_are_monotonic():
    gt = np.zeros((40, 40), np.int32); gt[5:35, 5:35] = 1
    pred = np.zeros((40, 40), np.int32); pred[7:35, 7:35] = 1
    assert f1_at_iou(gt, pred, 0.5)["F1"] >= f1_at_iou(gt, pred, 0.95)["F1"]


def test_filter_by_valid_drops_ignored_predictions_but_keeps_trusted_ones():
    pred = np.zeros((20, 40), np.int32)
    pred[5:15, 2:12] = 1        # in the valid half
    pred[5:15, 22:32] = 2       # in the ignored half
    valid = np.zeros((20, 40), bool)
    valid[:, :20] = True
    out = filter_by_valid(pred, valid, 0.5)
    assert set(np.unique(out)) == {0, 1}


def test_evaluate_instances_bundle_has_every_reported_key():
    gt = _two_boxes()
    m = evaluate_instances(gt, gt.copy())
    for k in ("PQ", "SQ", "DQ", "AJI", "split_rate", "merge_rate", "mean_dice",
              "n_gt", "n_pred", "count_error", "F1@0.5"):
        assert k in m, k


def test_semantic_scores_ignore_and_absent_classes():
    gt = np.zeros((10, 10), np.int64)
    gt[2:6, 2:6] = 1
    gt[8:, :] = 255                       # ignored band
    pred = np.zeros((10, 10), np.int64)
    pred[2:6, 2:6] = 1
    pred[8:, :] = 3                       # wrong, but ignored
    s = semantic_scores(gt, pred, 4, ["background", "tubule", "glomerulus", "vessel"])
    assert s["per_class"]["tubule"]["dice"] == pytest.approx(1.0)
    assert s["per_class"]["glomerulus"]["dice"] is None, "absent class must not score 1.0"
    assert s["mean_dice"] == pytest.approx(1.0)
    assert s["n_valid_px"] == 80


def test_glass_rate_reports_count_and_area_separately():
    """They diverge: one huge blob on glass is a small count and a large area."""
    from kidney_panoptic.eval.diagnostics import glass_rate

    id_map = np.zeros((40, 80), np.int32)
    id_map[10:30, 5:25] = 1        # on tissue, 400 px
    id_map[2:38, 42:78] = 2        # on glass, 1296 px — one instance, most of the area
    tissue = np.zeros((40, 80), bool); tissue[:, :40] = True

    g = glass_rate(id_map, tissue)
    assert g["n_instances"] == 2 and g["n_glass_instances"] == 1
    assert g["glass_instance_rate"] == pytest.approx(0.5)
    assert g["glass_area_rate"] > 0.7, "by area the failure is much larger than by count"
    assert g["predicted_area_px"] == 400 + 1296


def test_glass_rate_is_zero_when_everything_is_on_tissue():
    from kidney_panoptic.eval.diagnostics import glass_rate

    id_map = np.zeros((20, 20), np.int32); id_map[5:15, 5:15] = 1
    g = glass_rate(id_map, np.ones((20, 20), bool))
    assert g["n_glass_instances"] == 0 and g["glass_area_rate"] == 0.0


def test_merged_groups_reports_which_gt_objects_were_absorbed():
    from kidney_panoptic.eval.instance_metrics import merged_groups

    gt = _two_boxes()
    pred = np.zeros_like(gt)
    pred[10:30, 5:75] = 1              # one prediction swallowing both GT objects
    groups = merged_groups(gt, pred)
    assert len(groups) == 1
    assert groups[0]["pred_id"] == 1
    assert sorted(groups[0]["gt_ids"]) == [1, 2]


def test_merged_groups_empty_when_nothing_is_fused():
    from kidney_panoptic.eval.instance_metrics import merged_groups

    gt = _two_boxes()
    assert merged_groups(gt, gt.copy()) == []


def test_subtype_breakdown_separates_same_from_cross_label_merges():
    from kidney_panoptic.eval.instance_metrics import subtype_merge_breakdown

    gt = np.zeros((40, 160), np.int32)
    gt[10:30, 5:35] = 1
    gt[10:30, 40:70] = 2               # merged with 1, DIFFERENT label
    gt[10:30, 85:110] = 3
    gt[10:30, 115:140] = 4             # merged with 3, SAME label
    pred = np.zeros_like(gt)
    pred[10:30, 5:70] = 1
    pred[10:30, 85:140] = 2

    labels = {1: "Tubules", 2: "Tubules2", 3: "Tubules", 4: "Tubules"}
    b = subtype_merge_breakdown(gt, pred, labels)
    assert b["n_merging_preds"] == 2
    assert b["n_cross_subtype"] == 1 and b["n_same_subtype"] == 1
    assert b["cross_subtype_fraction"] == pytest.approx(0.5)
    assert {tuple(sorted(g["labels"])) for g in b["groups"]} == {
        ("Tubules", "Tubules2"), ("Tubules", "Tubules")}


def test_subtype_breakdown_is_none_when_there_are_no_merges():
    from kidney_panoptic.eval.instance_metrics import subtype_merge_breakdown

    gt = _two_boxes()
    b = subtype_merge_breakdown(gt, gt.copy(), {1: "Tubules", 2: "Tubules2"})
    assert b["n_merging_preds"] == 0 and b["cross_subtype_fraction"] is None
