"""Masked segmentation losses.

Every loss here takes an explicit ``valid`` mask and normalises over it. On partially
annotated slides most of the tile carries no usable label, so a loss that averages over
all pixels silently rescales itself with the annotation density and stops being comparable
between patches.

One consequence worth stating: because each configuration normalises over its OWN valid
region, total loss values are **not comparable across different ignore strategies**. That
is exactly why checkpoint selection uses PQ (docs/adr/0004).
"""
from __future__ import annotations

from typing import Optional

import torch
import torch.nn.functional as F

IGNORE_INDEX = 255


def masked_cross_entropy(
    logits: torch.Tensor,          # (B, C, H, W)
    target: torch.Tensor,          # (B, H, W) int64, IGNORE_INDEX where invalid
    weight: Optional[torch.Tensor] = None,       # (C,) class weights
    pixel_weight: Optional[torch.Tensor] = None,  # (B, H, W) per-pixel multiplier
    ignore_index: int = IGNORE_INDEX,
) -> torch.Tensor:
    ce = F.cross_entropy(logits, target.clamp(max=logits.shape[1] - 1).long(),
                         weight=weight, reduction="none")
    valid = (target != ignore_index).float()
    if pixel_weight is not None:
        ce = ce * pixel_weight
        denom = (valid * pixel_weight).sum()
    else:
        denom = valid.sum()
    return (ce * valid).sum() / denom.clamp(min=1.0)


def masked_dice(
    logits: torch.Tensor,          # (B, C, H, W)
    target: torch.Tensor,          # (B, H, W) int64
    valid: torch.Tensor,           # (B, H, W) float {0,1}
    *,
    exclude_background: bool = True,
    eps: float = 1.0,
) -> torch.Tensor:
    """Soft multi-class Dice over valid pixels; returns 1 - mean Dice.

    Only classes PRESENT in the target (or predicted with meaningful mass) contribute, so
    a patch containing no vessels is not penalised for correctly predicting no vessels.
    """
    B, C, H, W = logits.shape
    probs = torch.softmax(logits, dim=1)
    v = valid.unsqueeze(1).float()
    tgt = target.clamp(max=C - 1).long()
    onehot = F.one_hot(tgt, C).permute(0, 3, 1, 2).float()

    probs, onehot = probs * v, onehot * v
    start = 1 if exclude_background else 0
    dims = (0, 2, 3)
    inter = (probs[:, start:] * onehot[:, start:]).sum(dims)
    denom = probs[:, start:].sum(dims) + onehot[:, start:].sum(dims)
    dice = (2 * inter + eps) / (denom + eps)

    present = onehot[:, start:].sum(dims) > 0
    if present.any():
        return 1.0 - dice[present].mean()
    return 1.0 - dice.mean()


def masked_tversky(
    logits: torch.Tensor,
    target: torch.Tensor,
    valid: torch.Tensor,
    alpha: float = 0.45,
    beta: float = 0.55,
    *,
    exclude_background: bool = True,
    eps: float = 1.0,
) -> torch.Tensor:
    """Tversky loss (alpha weights FP, beta weights FN).

    Kept near-balanced by default. A strong recall bias only makes sense when background
    is trusted everywhere; with an ignore mask it mostly inflates predictions into the
    unlabelled region where nothing pushes back.
    """
    B, C, H, W = logits.shape
    probs = torch.softmax(logits, dim=1)
    v = valid.unsqueeze(1).float()
    onehot = F.one_hot(target.clamp(max=C - 1).long(), C).permute(0, 3, 1, 2).float()
    probs, onehot = probs * v, onehot * v

    start = 1 if exclude_background else 0
    dims = (0, 2, 3)
    tp = (probs[:, start:] * onehot[:, start:]).sum(dims)
    fp = (probs[:, start:] * (1 - onehot[:, start:])).sum(dims)
    fn = ((1 - probs[:, start:]) * onehot[:, start:]).sum(dims)
    tv = (tp + eps) / (tp + alpha * fp + beta * fn + eps)

    present = onehot[:, start:].sum(dims) > 0
    return 1.0 - (tv[present].mean() if present.any() else tv.mean())


def masked_mse(
    pred: torch.Tensor,            # (B, K, H, W)
    target: torch.Tensor,          # (B, K, H, W)
    valid: torch.Tensor,           # (B, H, W)
) -> torch.Tensor:
    if pred.shape[1] == 0:
        return pred.sum() * 0.0
    v = valid.unsqueeze(1).float()
    se = (pred - target) ** 2 * v
    return se.sum() / (v.sum() * pred.shape[1]).clamp(min=1.0)


def class_weights_from_counts(
    counts: torch.Tensor,
    mode: str = "inverse_sqrt_freq",
    clip: float = 12.0,
) -> torch.Tensor:
    """Class weights from pixel or instance counts, normalised to mean 1.

    ``inverse_sqrt_freq`` is the default: full inverse frequency on a 5-class kidney
    problem puts a weight of several hundred on nerve, which destabilises training long
    before it improves nerve recall.
    """
    counts = counts.float().clamp(min=1.0)
    if mode == "none":
        return torch.ones_like(counts)
    freq = counts / counts.sum()
    w = 1.0 / freq if mode == "inverse_freq" else 1.0 / freq.sqrt()
    w = w.clamp(max=float(clip) * w.min())
    return w / w.mean()
