"""Instance metrics: PQ / SQ / DQ, AJI, F1@IoU, mean Dice, split and merge rates.

Inputs are integer instance-id maps (0 = background). Two separate matchings are used and
they must not be shared:

* **PQ** uses unique IoU > 0.5 matching. Above 0.5 the match is provably unique, so the
  Hungarian assignment only resolves ties in the aggregate-IoU sense.
* **AJI** uses its own greedy per-GT argmax matching, replicating Kumar et al. 2017.
  Reusing PQ's matching here silently changes the number.

``DQ`` is the detection/recognition quality ``TP / (TP + 0.5 FP + 0.5 FN)``; it appears in
the literature under both names and is reported as ``DQ`` with ``RQ`` as an alias.

Split / merge rates are computed at a LOWER overlap threshold than PQ matching on purpose:
a GT object broken into two halves may have no prediction above IoU 0.5 at all, so counting
fragments at IoU 0.5 would report the failure as a plain miss and hide the mechanism.
"""
from __future__ import annotations

from typing import Optional

import numpy as np
from scipy.optimize import linear_sum_assignment


def instance_ids(lbl: np.ndarray) -> list[int]:
    return [int(i) for i in np.unique(lbl) if i != 0]


def _contingency(gt: np.ndarray, pred: np.ndarray):
    """Sparse-ish intersection matrix via a single bincount over paired labels."""
    gt_ids, pred_ids = instance_ids(gt), instance_ids(pred)
    G, P = len(gt_ids), len(pred_ids)
    inter = np.zeros((G, P), dtype=np.int64)
    if G == 0 or P == 0:
        return inter, gt_ids, pred_ids, np.zeros(G, np.int64), np.zeros(P, np.int64)

    gt_pos = {v: i for i, v in enumerate(gt_ids)}
    pred_pos = {v: i for i, v in enumerate(pred_ids)}
    gmax, pmax = max(gt_ids), max(pred_ids)
    g_lut = np.zeros(gmax + 1, np.int64)
    p_lut = np.zeros(pmax + 1, np.int64)
    for v, i in gt_pos.items():
        g_lut[v] = i + 1
    for v, i in pred_pos.items():
        p_lut[v] = i + 1

    gi = g_lut[np.clip(gt, 0, gmax)]
    pi = p_lut[np.clip(pred, 0, pmax)]
    both = (gi > 0) & (pi > 0)
    if both.any():
        flat = (gi[both] - 1) * P + (pi[both] - 1)
        counts = np.bincount(flat, minlength=G * P)
        inter = counts.reshape(G, P)

    gt_areas = np.array([(gt == v).sum() for v in gt_ids], dtype=np.int64)
    pred_areas = np.array([(pred == v).sum() for v in pred_ids], dtype=np.int64)
    return inter, gt_ids, pred_ids, gt_areas, pred_areas


def pairwise_iou(gt: np.ndarray, pred: np.ndarray):
    inter, gt_ids, pred_ids, ga, pa = _contingency(gt, pred)
    if inter.size == 0:
        return np.zeros(inter.shape), gt_ids, pred_ids, inter, ga, pa
    union = ga[:, None] + pa[None, :] - inter
    iou = np.where(union > 0, inter / np.maximum(union, 1), 0.0)
    return iou, gt_ids, pred_ids, inter, ga, pa


# ------------------------------------------------------------------------- PQ
def panoptic_quality(gt: np.ndarray, pred: np.ndarray, iou_threshold: float = 0.5) -> dict:
    iou, gt_ids, pred_ids, *_ = pairwise_iou(gt, pred)
    n_gt, n_pred = len(gt_ids), len(pred_ids)
    if n_gt == 0 and n_pred == 0:
        return {"PQ": 1.0, "SQ": 1.0, "DQ": 1.0, "RQ": 1.0, "TP": 0, "FP": 0, "FN": 0,
                "n_gt": 0, "n_pred": 0}
    pairs = []
    if n_gt and n_pred:
        rows, cols = linear_sum_assignment(-(iou * (iou > iou_threshold)))
        pairs = [(r, c) for r, c in zip(rows, cols) if iou[r, c] > iou_threshold]
    tp = len(pairs)
    fp, fn = n_pred - tp, n_gt - tp
    sq = float(sum(iou[r, c] for r, c in pairs) / tp) if tp else 0.0
    dq = tp / (tp + 0.5 * fp + 0.5 * fn) if (tp + fp + fn) else 0.0
    return {"PQ": sq * dq, "SQ": sq, "DQ": dq, "RQ": dq, "TP": tp, "FP": fp, "FN": fn,
            "n_gt": n_gt, "n_pred": n_pred}


def aggregated_jaccard_index(gt: np.ndarray, pred: np.ndarray) -> float:
    """AJI with its own greedy matching (Kumar et al. 2017)."""
    iou, gt_ids, pred_ids, inter, ga, pa = pairwise_iou(gt, pred)
    if not gt_ids:
        return 1.0 if not pred_ids else 0.0
    if not pred_ids:
        return 0.0
    used = np.zeros(len(pred_ids), bool)
    num = den = 0
    for gi in range(len(gt_ids)):
        if iou[gi].max(initial=0) == 0:
            den += int(ga[gi])
            continue
        pj = int(np.argmax(iou[gi]))
        num += int(inter[gi, pj])
        den += int(ga[gi] + pa[pj] - inter[gi, pj])
        used[pj] = True
    den += int(pa[~used].sum())
    return float(num / den) if den else 0.0


def f1_at_iou(gt: np.ndarray, pred: np.ndarray, threshold: float = 0.5) -> dict:
    iou, gt_ids, pred_ids, *_ = pairwise_iou(gt, pred)
    n_gt, n_pred = len(gt_ids), len(pred_ids)
    tp = 0
    if n_gt and n_pred:
        rows, cols = linear_sum_assignment(-(iou * (iou >= threshold)))
        tp = sum(1 for r, c in zip(rows, cols) if iou[r, c] >= threshold)
    fp, fn = n_pred - tp, n_gt - tp
    prec = tp / (tp + fp) if (tp + fp) else (1.0 if n_gt == 0 else 0.0)
    rec = tp / (tp + fn) if (tp + fn) else (1.0 if n_pred == 0 else 0.0)
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
    return {"F1": f1, "precision": prec, "recall": rec, "TP": tp, "FP": fp, "FN": fn}


def mean_matched_dice(gt: np.ndarray, pred: np.ndarray, iou_threshold: float = 0.5) -> dict:
    iou, gt_ids, pred_ids, *_ = pairwise_iou(gt, pred)
    if not gt_ids or not pred_ids:
        return {"mean_iou": 0.0, "mean_dice": 0.0}
    rows, cols = linear_sum_assignment(-(iou * (iou > iou_threshold)))
    matched = [iou[r, c] for r, c in zip(rows, cols) if iou[r, c] > iou_threshold]
    if not matched:
        return {"mean_iou": 0.0, "mean_dice": 0.0}
    return {"mean_iou": float(np.mean(matched)),
            "mean_dice": float(np.mean([2 * v / (1 + v) for v in matched]))}


# --------------------------------------------------------------- split / merge
def split_merge_rates(gt: np.ndarray, pred: np.ndarray,
                      overlap_frac: float = 0.20) -> dict:
    """Fragmentation and fusion rates.

    A GT object is SPLIT when >= 2 predictions each cover at least ``overlap_frac`` of it.
    A prediction MERGES when it covers at least ``overlap_frac`` of >= 2 GT objects; the
    merge rate is reported over GT objects (how many GT objects were fused into something).

    ``overlap_frac`` is intentionally well below the PQ threshold: two halves of a split
    tubule may each sit at IoU ~0.45 and would otherwise be invisible to the metric.
    """
    inter, gt_ids, pred_ids, ga, pa = _contingency(gt, pred)
    n_gt, n_pred = len(gt_ids), len(pred_ids)
    if n_gt == 0:
        return {"split_rate": 0.0, "merge_rate": 0.0, "n_split": 0, "n_merged_gt": 0,
                "n_gt": 0, "n_merging_preds": 0}

    gt_cov = inter / np.maximum(ga[:, None], 1)               # frac of GT covered by pred
    n_split = int((gt_cov >= overlap_frac).sum(axis=1).__ge__(2).sum()) if n_pred else 0

    n_merged_gt = n_merging = 0
    if n_pred:
        hits = gt_cov >= overlap_frac                          # (G, P)
        per_pred = hits.sum(axis=0)
        merging = per_pred >= 2
        n_merging = int(merging.sum())
        n_merged_gt = int(hits[:, merging].any(axis=1).sum()) if n_merging else 0

    return {
        "split_rate": n_split / n_gt,
        "merge_rate": n_merged_gt / n_gt,
        "n_split": n_split,
        "n_merged_gt": n_merged_gt,
        "n_merging_preds": n_merging,
        "n_gt": n_gt,
    }


def merged_groups(gt: np.ndarray, pred: np.ndarray,
                  overlap_frac: float = 0.20) -> list[dict]:
    """For each merging prediction, which GT ids it absorbed.

    ``split_merge_rates`` reduces this to a rate; this returns the groups themselves so a
    caller can ask *what kind* of objects get fused — e.g. whether the two halves of a merge
    are the same tubule subtype or different ones. That question decides whether a
    class-aware split can help at all, so it needs the identities, not just the count.
    """
    inter, gt_ids, pred_ids, ga, _ = _contingency(gt, pred)
    if len(gt_ids) == 0 or len(pred_ids) == 0:
        return []
    hits = (inter / np.maximum(ga[:, None], 1)) >= overlap_frac
    groups = []
    for pj, pid in enumerate(pred_ids):
        absorbed = [int(gt_ids[gi]) for gi in np.nonzero(hits[:, pj])[0]]
        if len(absorbed) >= 2:
            groups.append({"pred_id": int(pid), "gt_ids": absorbed})
    return groups


def subtype_merge_breakdown(
    gt: np.ndarray,
    pred: np.ndarray,
    label_of_gt_id: dict[int, str],
    *,
    overlap_frac: float = 0.20,
) -> dict:
    """Split merged groups into same-label vs cross-label.

    Interpretation needs the NULL — the rate at which adjacent objects differ in label
    anyway. Measured on this cohort, only 27.4 % of adjacent tubule pairs are cross-subtype,
    so a cross-subtype merge fraction near 27 % means merging is subtype-INDEPENDENT and a
    class-aware split will not help much; well above it means the subtype border really is
    where the model fails.
    """
    groups = merged_groups(gt, pred, overlap_frac)
    same = cross = 0
    detail = []
    for g in groups:
        labels = [label_of_gt_id.get(i, "?") for i in g["gt_ids"]]
        is_cross = len(set(labels)) > 1
        cross += is_cross
        same += not is_cross
        detail.append({**g, "labels": labels, "cross_subtype": bool(is_cross)})
    total = same + cross
    return {
        "n_merging_preds": total,
        "n_same_subtype": same,
        "n_cross_subtype": cross,
        "cross_subtype_fraction": (cross / total) if total else None,
        "groups": detail,
    }


def count_error(gt: np.ndarray, pred: np.ndarray) -> dict:
    n_gt, n_pred = len(instance_ids(gt)), len(instance_ids(pred))
    return {"n_gt": n_gt, "n_pred": n_pred, "count_error": n_pred - n_gt,
            "count_error_pct": (100.0 * (n_pred - n_gt) / n_gt) if n_gt else float("nan")}


# ------------------------------------------------------------------- ignore
def filter_by_valid(pred: np.ndarray, valid: np.ndarray,
                    min_valid_frac: float = 0.5) -> np.ndarray:
    """Drop predicted instances that lie mostly in the ignored region.

    On partially annotated slides a prediction sitting on an unannotated (but real)
    structure is not a false positive — there is simply no label there. Predictions on
    TRUSTED background (glass, the ring around annotations) stay, because those genuinely
    are false positives.
    """
    out = pred.copy()
    valid = valid.astype(bool)
    for pid in instance_ids(pred):
        m = pred == pid
        if (valid & m).sum() / max(m.sum(), 1) < min_valid_frac:
            out[m] = 0
    return out


def evaluate_instances(
    gt: np.ndarray,
    pred: np.ndarray,
    *,
    valid: Optional[np.ndarray] = None,
    iou_threshold: float = 0.5,
    f1_thresholds: tuple[float, ...] = (0.5, 0.75),
    split_overlap_frac: float = 0.20,
    min_valid_frac: float = 0.5,
) -> dict:
    """Full instance-metric bundle for one region."""
    if valid is not None:
        pred = filter_by_valid(pred, valid, min_valid_frac)
    out: dict = {}
    out.update(panoptic_quality(gt, pred, iou_threshold))
    out["AJI"] = aggregated_jaccard_index(gt, pred)
    out.update(mean_matched_dice(gt, pred, iou_threshold))
    out.update(split_merge_rates(gt, pred, split_overlap_frac))
    out.update(count_error(gt, pred))
    for t in f1_thresholds:
        r = f1_at_iou(gt, pred, t)
        out[f"F1@{t}"] = r["F1"]
    return out


def evaluate_per_class(
    gt: np.ndarray,
    pred: np.ndarray,
    gt_class: dict[int, int],
    pred_class: dict[int, int],
    class_names: list[str],
    **kwargs,
) -> dict:
    """Per-class instance metrics plus the class-agnostic overall pass.

    Each class is evaluated on id maps masked to that class, so a tubule predicted where a
    vessel is counts as both a vessel FN and a tubule FP — which is the honest accounting.
    """
    res: dict = {"overall": evaluate_instances(gt, pred, **kwargs)}
    for ci, name in enumerate(class_names):
        if ci == 0:
            continue
        g = np.where(np.isin(gt, [i for i, c in gt_class.items() if c == ci]), gt, 0)
        p = np.where(np.isin(pred, [i for i, c in pred_class.items() if c == ci]), pred, 0)
        if len(instance_ids(g)) == 0 and len(instance_ids(p)) == 0:
            continue
        res[name] = evaluate_instances(g, p, **kwargs)
    thing = [v for k, v in res.items() if k != "overall"]
    if thing:
        res["mean_class"] = {
            "PQ": float(np.mean([v["PQ"] for v in thing])),
            "SQ": float(np.mean([v["SQ"] for v in thing])),
            "DQ": float(np.mean([v["DQ"] for v in thing])),
            "AJI": float(np.mean([v["AJI"] for v in thing])),
        }
    return res
