"""Composite panoptic loss: semantic + boundary/interior + centers.

Per-head terms are returned separately so a training log can show which head is actually
moving. The boundary head carries an extra class weight on channel 2 (the separator) and
an extra per-pixel weight on instance-instance contact pixels — those two together are the
only supervision that decides whether two touching tubules come out as one object or two,
and they are a tiny fraction of the pixels.
"""
from __future__ import annotations

from typing import Optional

import torch
import torch.nn as nn

from .seg_losses import (
    IGNORE_INDEX,
    masked_cross_entropy,
    masked_dice,
    masked_mse,
    masked_tversky,
)


class PanopticLoss(nn.Module):
    def __init__(self, cfg: dict, n_classes: int, semantic_class_weights: Optional[torch.Tensor] = None):
        super().__init__()
        lc = dict(cfg.get("loss", {}) or {})
        self.sem = dict(lc.get("semantic", {}) or {})
        self.bnd = dict(lc.get("boundary", {}) or {})
        self.cen = dict(lc.get("center", {}) or {})
        self.tversky = dict(lc.get("tversky", {}) or {})
        self.n_classes = n_classes

        if semantic_class_weights is not None:
            self.register_buffer("sem_weights", semantic_class_weights.float())
        else:
            self.sem_weights = None

        bw = float(self.bnd.get("boundary_class_weight", 3.0))
        tw = float(self.bnd.get("touching_class_weight", 6.0))
        self.boundary_classes = int(cfg.get("model", {}).get("boundary_classes", 3))
        if self.boundary_classes == 4:
            # (background, interior, boundary-to-background, boundary-to-INSTANCE)
            w = [1.0, 1.0, bw, tw]
        else:
            # (background, interior, boundary)
            w = [1.0, 1.0, bw]
        self.register_buffer("bnd_weights", torch.tensor(w, dtype=torch.float32))

    def forward(self, out: dict, batch: dict) -> tuple[torch.Tensor, dict]:
        valid = batch["valid"]
        pixel_w = batch.get("weight")

        # ---- semantic ----
        sem_ce = masked_cross_entropy(
            out["semantic_logits"], batch["semantic"],
            weight=self.sem_weights, ignore_index=IGNORE_INDEX,
        )
        sem_dice = masked_dice(out["semantic_logits"],
                               batch["semantic"].clamp(max=self.n_classes - 1), valid)
        l_sem = float(self.sem.get("ce", 1.0)) * sem_ce + float(self.sem.get("dice", 1.0)) * sem_dice

        # ---- boundary / interior ----
        bnd_ce = masked_cross_entropy(
            out["boundary_logits"], batch["boundary"],
            weight=self.bnd_weights, pixel_weight=pixel_w, ignore_index=IGNORE_INDEX,
        )
        bnd_tgt = batch["boundary"].clamp(max=self.boundary_classes - 1)
        if bool(self.tversky.get("enabled", False)):
            bnd_region = masked_tversky(
                out["boundary_logits"], bnd_tgt, valid,
                alpha=float(self.tversky.get("alpha", 0.45)),
                beta=float(self.tversky.get("beta", 0.55)),
            )
        else:
            bnd_region = masked_dice(out["boundary_logits"], bnd_tgt, valid)
        l_bnd = float(self.bnd.get("ce", 1.0)) * bnd_ce + float(self.bnd.get("dice", 1.0)) * bnd_region

        # ---- centers ----
        l_cen = masked_mse(out["center_heatmaps"], batch["centers"], valid)
        l_cen = float(self.cen.get("mse", 25.0)) * l_cen

        total = l_sem + l_bnd + l_cen
        parts = {
            "loss": float(total.detach()),
            "sem_ce": float(sem_ce.detach()),
            "sem_dice": float(sem_dice.detach()),
            "bnd_ce": float(bnd_ce.detach()),
            "bnd_region": float(bnd_region.detach()),
            "center_mse": float(l_cen.detach()),
            "loss_semantic": float(l_sem.detach()),
            "loss_boundary": float(l_bnd.detach()),
        }
        return total, parts
