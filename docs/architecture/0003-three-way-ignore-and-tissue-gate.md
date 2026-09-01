# ADR 0003 — Three-way ignore mask and the inference tissue gate

**Status:** accepted
**Date:** 2026-07-29

## Context

Annotation on these slides is **partial**: a few hundred objects are drawn per slide out of
tens of thousands present. Two naive readings both fail:

- *Everything unannotated is background.* Trains the model to suppress real tubules; the
  majority of its "negative" pixels are positives.
- *Everything unannotated is unknown.* Then nothing anywhere pushes the foreground head
  toward zero. On the predecessor this was measured directly: fed a **constant gray image**
  the model still emitted 42 saturated instances covering 26 % of the tile, with the pattern
  essentially identical (r = 0.993–1.000) from gray 200 to 255 — a pure decoder prior with
  no input dependence. On production slides 17–26 % of predicted instances sat on blank
  glass, and slide-wide sampling put it at **54.7 %**.

The way out is that one kind of background needs no annotation at all. **Glass is
identifiable without a label.** H&E stain is chromatic; slide background is not.

## Decision

A three-way mask (`data/ignore.py`, default
`ignore_strategy: dilated_gt_two_radius_plus_tissue`):

| region | source | treated as |
|---|---|---|
| annotated instance pixels | GeoJSON | trusted foreground |
| a 5 px ring around each instance | dilation | trusted background |
| everything off tissue | saturation-based tissue mask | trusted background |
| annotated `Background` polygons | GeoJSON | trusted background |
| inside a dense ROI | GeoJSON | trusted (fully labelled) |
| **tissue without annotation** | — | **ignored** |

Every loss normalises over the valid region only. At inference the same tissue mask is
applied as a **gate**: instances lying more than half off tissue are dropped
(`infer.tissue_gate`, `tissue_min_instance_frac`).

## Why the gate as well as the supervision

They fail differently and the cost of each error is asymmetric. Supervision fixes the cause
(the model learns glass is background); the gate is a cheap, deterministic backstop for the
residue and for stain/scanner shifts the training set never saw. Calling tissue "glass"
would suppress real structures, so the tissue parameters are deliberately conservative —
a low saturation threshold plus a 16 px outward dilation, so the tissue *edge* never lands
in trusted background.

## Evaluation policy

Evaluation uses **the same mask**. Predictions lying mostly in the *ignored* region are
dropped before scoring — on a partially annotated slide a prediction on an undrawn but real
tubule is not a false positive, and counting it as one measures annotation density rather
than model quality. Predictions on *trusted background* (glass, the ring) are **kept and do
count against the score**.

This is a deliberate departure from the predecessor, which stripped the `_plus_tissue`
term at scoring time for back-compatibility with an older results table. Greenfield has no
such obligation, and the honest accounting is worth more: with it, PQ and the glass rate
cannot tell contradictory stories about the same output.

## Consequences

- `patch.bg_glass_keep_ratio` must stay non-zero. If no glass-containing patch enters the
  training set, the trusted-background term has nothing to act on. Extraction therefore
  keeps glassy background patches preferentially (60 %) over tissue-only ones (3 %) — the
  latter are nearly all ignored pixels and buy almost no gradient for a full forward pass.
- The tissue thresholds are calibrated for H&E at this scanner. A new stain needs
  `scripts/make_tissue_mask.py` run and the numbers re-checked before training.
- Reported metrics are conditional on annotation coverage and are **not** comparable to
  numbers computed on densely labelled data. Where dense ROIs exist, `ignore_strategy:
  dense_roi` gives an unbiased estimate; without them the bias direction should be stated
  whenever a number is quoted.
- `constant_input_canary` deliberately runs with the gate **disabled**. Gating would zero
  the metric by construction and hide the very prior it exists to detect.
