# ADR 0001 — Instance representation: interior/boundary, not centroid-HV

**Status:** accepted
**Date:** 2026-07-29

## Context

The instances to separate are mostly **tubule cross-sections**: elongated, irregular, and
touching with an inter-instance gap that is frequently zero pixels. A minority — glomeruli,
blood cells — are compact and well separated. Both must come out of the same model.

The predecessor system used the HoVer-Net representation: per-pixel horizontal and vertical
offsets to the owning object's centroid, decoded by taking a signed Sobel of the offset
fields to build an energy landscape, then watershedding it.

## Decision

The default instance representation is a **3-class per-pixel field — background / instance
interior / instance boundary** — decoded by connected components on the interior and a
marker-controlled watershed of the boundary probability. Compact classes additionally get a
**center heatmap**. Centroid-HV is not used; it survives only behind `model.legacy_hv`
(default `false`) so a comparison can be run without forking the code.

## Why

**1. Centroid offset is a non-local quantity, and everything downstream inherits that.**
The value at a pixel depends on where the object's far end is, which may be outside the
tile. Two adjacent tiles that disagree about a pixel's owner produce a step discontinuity
in the field, and the signed-Sobel decode reads any step edge as a boundary. This is not
hypothetical: on the v1 system it produced a measured 3–5× excess of straight cuts exactly
on tile seams, and it is a property of the representation, not of any particular
implementation of it.

Interior/boundary has no such coupling. A pixel's label depends only on what is within
`boundary_width_px` of it. The decode is translation-equivariant and window-invariant by
construction (`tests/test_decode.py::test_decode_is_translation_equivariant`,
`::test_decode_on_a_window_matches_decode_on_the_whole`), so no amount of blending, halo
sizing or seam repair is needed to make tiling safe — the problem does not arise.

**2. Centroid offset is ill-conditioned for exactly the objects we care about.** For a long
curved tubule the centroid may not lie inside the object at all, so the offset field points
at a location that is not part of the instance; and the gradient magnitude that separates
two objects scales with object size, so one threshold cannot serve tubules and blood cells
at once. A boundary band is the same 3 px whatever the object's size or aspect ratio.

**3. The evidence for the split is local anyway.** What actually distinguishes two touching
tubules in H&E is a thin wall — a local intensity/texture feature a few pixels wide. A
representation whose supervision sits exactly on that wall uses the available evidence
directly. Predicting a centroid vector asks the network to solve a harder, global problem
first and then recovers the boundary by differentiating the answer.

**4. Regional-vs-local seeding stops being a dilemma.** The v1 decode had to choose between
`h_maxima` (regional, correct for elongated objects, far too slow) and `peak_local_max`
(local, fast, double-seeds elongated objects along their distance-transform ridge), and
needed a `dist_min_marker_distance` fudge to paper over the difference. Interior components
are seeded by connected components, which is regional *and* linear-time: one connected
interior is one seed no matter how elongated it is.

## Consequences

- Two touching objects are separated iff the boundary head fires between them. When it does
  not, they merge — a merge failure rather than the split failure HV tended toward. The
  boundary class carries an extra loss weight (`boundary_class_weight: 3.0`) and
  instance-instance contact pixels carry another (`touching_weight: 2.0`) precisely because
  those pixels are a small fraction of the image and decide the outcome.
- The boundary band is carved out of the objects, so predicted instances are eroded by
  roughly `boundary_width_px`. This is a systematic area bias affecting morphometrics.
  Watershed flooding restores most of it (basins expand to fill the foreground mask), but
  the residual bias should be measured before areas are reported as absolute.
- Compact classes needed a second mechanism: two adjacent glomeruli are separated by a thin
  capsule the boundary head often misses. That is what the center head is for
  (ADR 0005), and it is why the decode is class-heterogeneous rather than uniform.
- A boundary-band representation cannot express nested or overlapping instances. Nothing in
  kidney morphometry needs that here.

## Alternatives considered

- **StarDist.** Strong for convex/star-shaped objects; tubule cross-sections are frequently
  not star-convex, and the radial representation is also object-global.
- **Per-pixel embeddings + clustering.** Handles arbitrary shape, but clustering is not
  window-invariant — the same failure class as HV, plus an unstable inference-time
  hyperparameter. Left as an optional extra head, not the default.
- **Mask2Former-style query-based panoptic.** Queries are global objects; tiling a gigapixel
  slide reintroduces cross-tile identity, and it needs far more labelled data than the few
  hundred annotated instances available here.
