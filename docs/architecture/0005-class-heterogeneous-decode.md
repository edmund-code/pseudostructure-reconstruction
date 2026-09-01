# ADR 0005 — Class-heterogeneous instance decode

**Status:** accepted
**Date:** 2026-07-29

## Context

The six classes are not variations on one object. Measured on this cohort's ground truth,
they differ by more than an order of magnitude in size and by a large factor in shape:

| class | typical form | separation from neighbours |
|---|---|---|
| tubule | elongated cross-section, irregular | touching, gap often 0 px |
| vessel | large, elongated, can span > 1000 px | usually isolated |
| glomerulus | compact, round, several hundred px across | isolated, occasionally adjacent |
| blood cell | tiny compact blob, tens of px | clustered |
| nerve | small, elongated | rare, isolated |

A single decode tuned for tubules under-segments blood cells and over-splits vessels; one
tuned for blood cells shatters everything else. Uniform parameters are a compromise that is
wrong for every class.

## Decision

One shared field prediction, then **per-class cleanup** in `postprocess/decode.py`:

- **All classes** — watershed of boundary probability seeded by interior components; class
  by majority semantic vote inside the instance with boundary pixels excluded.
- **Compact classes** (`compact_classes: [glomerulus, blood_cell]`) — center heatmap peaks.
  An instance containing ≥ 2 peaks is re-split by flooding its **negative distance
  transform** from those peaks. Geometric, not boundary-driven, because for round objects
  the distance transform is the reliable cue and the capsule between two adjacent glomeruli
  is often invisible to the boundary head. An instance with **no** peak is kept, not
  dropped: the center head is auxiliary, and treating its silence as a veto costs recall
  for nothing.
- **Elongated classes** — optional conservative re-merge, `merge_fragments`, **off by
  default**. Both a weak shared interface *and* small minor-axis growth are required.
- **Every class** — its own minimum area (`decode.min_area`), from 40 px² for blood cells to
  600 px² for glomeruli.

Vessels are explicitly excluded from the tubule compactness priors: a legitimate vessel can
exceed 1300 px in extent, and a size or minor-axis cap tuned for tubules would truncate it.

## The minor-axis veto

Two tubules lying **parallel and touching** also have a weak interface, so an
interface-only merge rule fuses them. The distinguishing signal is geometric: fusing two
parallel objects roughly doubles the minor axis, whereas rejoining two fragments of one
object barely changes it. `merge_max_minor_axis_growth: 1.15` is a hard reject.

Merging is nevertheless **off by default**. A wrong merge destroys two objects and is
invisible in per-pixel metrics; a wrong split is visible and recoverable. The predecessor's
equivalent interface rule measured AUC 0.823 at ~16–21 % precision — not good enough to run
unattended. `tests/test_decode.py::test_merge_fragments_rejects_parallel_glue` pins the
veto so the rule cannot silently regress if someone does enable it.

## Consequences

- More knobs. They are all in one config block with per-class keys, and the defaults are
  derived from GT geometry rather than tuned on the metric.
- The class vote happens *after* the split, so a mis-split object can end up with two
  different class labels. Rare in practice; would need an instance-level classifier head to
  fix properly.
- Center supervision only exists for compact classes, so `split_compact_by_centers` cannot
  help a fused pair of tubules. That is by design — for tubules the boundary head is the
  right instrument, and adding centers for elongated objects would reintroduce the
  centroid-is-non-local problem of ADR 0001.
