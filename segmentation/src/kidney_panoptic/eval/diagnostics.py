"""Failure-mode probes that ordinary segmentation metrics do not catch.

PQ on a hand-picked tissue-dense region can look healthy while the model invents hundreds
of objects on empty glass, fragments differently every time the tile grid moves, or
responds to a featureless input with a fixed positional prior. Each probe below targets one
of those, and each is cheap enough to run on every checkpoint.
"""
from __future__ import annotations

from typing import Optional, Sequence

import numpy as np

from ..data.classes import ClassSpec
from ..data.tissue import build_tissue_mask
from ..eval.instance_metrics import pairwise_iou
from ..infer.tiled import TiledPredictor, array_fetcher


def constant_input_canary(
    predictor: TiledPredictor,
    spec: ClassSpec,
    *,
    size: int = 1024,
    levels: Sequence[int] = (200, 220, 235, 255),
) -> dict:
    """Feed featureless gray images and count what comes out. Should be nothing.

    The tissue gate is DISABLED here on purpose. Gating would zero this metric by
    construction and hide the thing it is meant to measure: whether the network itself has
    learned a position-locked foreground prior that fires without any image evidence.
    """
    cfg = dict(predictor.cfg)
    cfg["infer"] = {**dict(cfg.get("infer", {}) or {}), "tissue_gate": False}
    saved, predictor.cfg = predictor.cfg, cfg
    try:
        per_level = {}
        for lv in levels:
            img = np.full((size, size, 3), int(lv), np.uint8)
            res, maps = predictor.predict_and_decode(
                array_fetcher(img), 0, 0, size, size, size, size, tissue=None
            )
            fg_prob = 1.0 - maps["boundary_prob"][0]
            per_level[int(lv)] = {
                "n_instances": int(res.n_instances),
                "instance_area_frac": float((res.id_map > 0).mean()),
                "mean_fg_prob": float(fg_prob.mean()),
                "max_fg_prob": float(fg_prob.max()),
            }
        return {
            "per_level": per_level,
            "mean_instances": float(np.mean([v["n_instances"] for v in per_level.values()])),
            "max_instances": int(max(v["n_instances"] for v in per_level.values())),
            "max_area_frac": float(max(v["instance_area_frac"] for v in per_level.values())),
            "passed": all(v["n_instances"] == 0 for v in per_level.values()),
        }
    finally:
        predictor.cfg = saved


def glass_rate(
    id_map: np.ndarray,
    tissue: np.ndarray,
    *,
    min_tissue_frac: float = 0.5,
) -> dict:
    """How much of the output sits on slide background, by COUNT and by AREA.

    Both are reported because they diverge: a handful of enormous hallucinated blobs on
    glass is a small count and a huge area, and quoting only the count makes that look fine.
    """
    tissue = tissue.astype(bool)
    # Per-instance work is done inside each bounding box. Scanning the whole array once
    # per instance is fine on a 4096 px region and completely infeasible slide-wide, where
    # there are ~10^5 instances over ~10^9 pixels.
    from scipy import ndimage as ndi

    n_glass = 0
    n_instances = 0
    for idx, sl in enumerate(ndi.find_objects(id_map), start=1):
        if sl is None:
            continue
        m = id_map[sl] == idx
        area = int(m.sum())
        if area == 0:
            continue
        n_instances += 1
        if int((tissue[sl] & m).sum()) / area < min_tissue_frac:
            n_glass += 1
    fg = id_map > 0
    glass_area = int((fg & ~tissue).sum())
    return {
        "n_instances": n_instances,
        "n_glass_instances": n_glass,
        "glass_instance_rate": n_glass / n_instances if n_instances else 0.0,
        "predicted_area_px": int(fg.sum()),
        "glass_area_px": glass_area,
        "glass_area_rate": glass_area / max(int(fg.sum()), 1),
        "glass_area_frac_of_glass": glass_area / max(int((~tissue).sum()), 1),
        "tissue_frac": float(tissue.mean()),
    }


def shift_stability(
    predictor: TiledPredictor,
    rgb: np.ndarray,
    *,
    shift: int = 128,
    size: Optional[int] = None,
    iou_threshold: float = 0.5,
    tissue_cfg: Optional[dict] = None,
) -> dict:
    """Re-infer the same content with the tile grid moved, and see if it agrees.

    The two runs read the SAME pixels; only the phase of the tiling grid relative to the
    content differs. Any instance that appears in one and not the other is an artefact of
    where the tile edges happened to fall. Reported as the fraction of instances matched at
    IoU >= threshold in the common area, from both directions.
    """
    H, W = rgb.shape[:2]
    size = int(size or min(H, W) - shift)
    if size <= shift:
        raise ValueError(f"Region ({H}x{W}) too small for shift {shift}.")

    a_img = rgb[0:size, 0:size]
    b_img = rgb[shift:shift + size, shift:shift + size]

    def run(img):
        tissue = build_tissue_mask(img, tissue_cfg) if tissue_cfg is not None else None
        res, _ = predictor.predict_and_decode(
            array_fetcher(img), 0, 0, img.shape[1], img.shape[0],
            img.shape[1], img.shape[0], tissue=tissue,
        )
        return res.id_map

    a = run(a_img)
    b = run(b_img)

    # Common content: a[shift:size, shift:size] == b[0:size-shift, 0:size-shift]
    ca = a[shift:size, shift:size]
    cb = b[0:size - shift, 0:size - shift]
    ca = _drop_border_touching(ca)
    cb = _drop_border_touching(cb)

    iou, a_ids, b_ids, *_ = pairwise_iou(ca, cb)
    n_a, n_b = len(a_ids), len(b_ids)
    if n_a == 0 and n_b == 0:
        return {"matched_fraction": 1.0, "n_a": 0, "n_b": 0, "mean_matched_iou": 1.0,
                "count_delta": 0}
    if n_a == 0 or n_b == 0:
        return {"matched_fraction": 0.0, "n_a": n_a, "n_b": n_b, "mean_matched_iou": 0.0,
                "count_delta": n_b - n_a}

    best_a = iou.max(axis=1)
    best_b = iou.max(axis=0)
    matched_a = int((best_a >= iou_threshold).sum())
    matched_b = int((best_b >= iou_threshold).sum())
    return {
        "matched_fraction": (matched_a + matched_b) / (n_a + n_b),
        "matched_fraction_a": matched_a / n_a,
        "matched_fraction_b": matched_b / n_b,
        "mean_matched_iou": float(best_a[best_a >= iou_threshold].mean())
        if matched_a else 0.0,
        "n_a": n_a,
        "n_b": n_b,
        "count_delta": n_b - n_a,
    }


def _drop_border_touching(labels: np.ndarray) -> np.ndarray:
    """Remove instances clipped by the crop — they are not comparable between runs."""
    border = np.unique(np.concatenate([
        labels[0], labels[-1], labels[:, 0], labels[:, -1]
    ]))
    out = labels.copy()
    for i in border:
        if i != 0:
            out[out == i] = 0
    return out
