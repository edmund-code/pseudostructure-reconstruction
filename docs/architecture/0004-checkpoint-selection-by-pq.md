# ADR 0004 — Select checkpoints by validation PQ, not validation loss

**Status:** accepted
**Date:** 2026-07-29

## Decision

`scripts/train.py` decodes a sample of validation patches every `train.eval_every` epochs,
computes pooled instance metrics, and saves `best_pq.pt` on improvement in **PQ**.
Validation loss is logged and never used for selection.

## Why

**1. Loss is not comparable across configurations.** Each ignore strategy normalises over
its own valid region. Two arms differing only in `ignore_strategy` produce loss values on
different denominators — one can be numerically lower while being the worse model. This was
observed directly on the predecessor, where the two Gate-1 arms' val losses (0.6397 vs
0.6249) were explicitly documented as incommensurable.

**2. Loss and instance quality diverge in a specific, predictable way.** The pixel losses
are dominated by interior and background — the overwhelming majority of pixels. A model can
reduce total loss by smoothing its boundary predictions, and smoothing is exactly what
fuses touching tubules. Per-pixel accuracy improves while the object count collapses. PQ is
computed *after* the watershed, so it sees the fusion; the loss cannot.

**3. What the system is for is counts and morphometrics.** Those are object-level
quantities. The selection metric should be the one whose failure modes are the ones that
matter.

## How it is computed

- Pooled, not averaged per patch: TP/FP/FN are summed across patches and SQ is
  TP-weighted. A 512 px patch may contain three objects; a per-patch mean lets a
  three-object patch swing the number as hard as a thirty-object one.
- Predictions in the ignored region are filtered out first (ADR 0003).
- `--max-val-patches` caps the count because the decode, not the forward pass, is the slow
  part of validation.

## Consequences

- Validation is ~10× more expensive than a loss-only pass. Hence `eval_every: 2`.
- PQ on a small validation set is noisy. With ~100 GT objects, PQ differences below roughly
  0.02 are not resolvable and must not be read as improvements — see the predecessor's
  Region 3 case, where three objects moving from matched to unmatched moved PQ by 0.06 on a
  25-object region.
- `last.pt` is saved every epoch alongside `best_pq.pt`, so a selection judged wrong in
  hindsight can be revisited without retraining.
