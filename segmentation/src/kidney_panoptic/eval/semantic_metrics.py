"""Per-class semantic Dice / IoU with ignore support."""
from __future__ import annotations

from typing import Optional

import numpy as np


def semantic_scores(
    gt: np.ndarray,
    pred: np.ndarray,
    n_classes: int,
    class_names: Optional[list[str]] = None,
    valid: Optional[np.ndarray] = None,
    ignore_index: int = 255,
) -> dict:
    """Dice and IoU per class, plus macro means over classes PRESENT in the ground truth.

    Absent classes are reported as ``null`` rather than 1.0 — a class with no GT and no
    prediction is undefined, and scoring it as perfect inflates the macro mean.
    """
    gt = np.asarray(gt)
    pred = np.asarray(pred)
    m = gt != ignore_index
    if valid is not None:
        m &= valid.astype(bool)
    g, p = gt[m], pred[m]

    names = class_names or [str(i) for i in range(n_classes)]
    per: dict[str, dict] = {}
    dices, ious = [], []
    for c in range(n_classes):
        gc, pc = (g == c), (p == c)
        inter = int((gc & pc).sum())
        gs, ps = int(gc.sum()), int(pc.sum())
        if gs == 0 and ps == 0:
            per[names[c]] = {"dice": None, "iou": None, "gt_px": 0, "pred_px": 0}
            continue
        dice = 2 * inter / (gs + ps) if (gs + ps) else 0.0
        iou = inter / (gs + ps - inter) if (gs + ps - inter) else 0.0
        per[names[c]] = {"dice": float(dice), "iou": float(iou), "gt_px": gs, "pred_px": ps}
        if gs > 0:
            dices.append(dice)
            ious.append(iou)

    fg = (g > 0)
    fg_p = (p > 0)
    fg_inter = int((fg & fg_p).sum())
    fg_dice = 2 * fg_inter / max(int(fg.sum()) + int(fg_p.sum()), 1)

    return {
        "per_class": per,
        "mean_dice": float(np.mean(dices)) if dices else 0.0,
        "mean_iou": float(np.mean(ious)) if ious else 0.0,
        "foreground_dice": float(fg_dice),
        "pixel_accuracy": float((g == p).mean()) if g.size else 0.0,
        "n_valid_px": int(m.sum()),
    }
