"""Supervision targets derived from an instance-id map.

Everything the network is trained on is a pure function of ``id_map`` (int32, 0 =
background, 1..N = instances) plus ``class_of_id``. Nothing is stored pre-computed, so
any geometric augmentation is applied to ``id_map`` and the targets are REGENERATED from
the warped result. Flip/rotation sign handling is then correct by construction — there is
no orientation-carrying field (that is one of the practical wins of dropping centroid-HV;
see docs/adr/0001).

Three target families
---------------------
``semantic``  (H, W) int64 in [0, n_classes)   — 0 = background.
``boundary``  (H, W) int64 in {0, 1, 2}        — 0 = background, 1 = instance interior,
                                                 2 = instance boundary band.
``centers``   (K, H, W) float32 in [0, 1]      — one Gaussian-peak channel per COMPACT
                                                 class (glomerulus, blood cell).

Why interior/boundary
---------------------
Adjacent tubules share a contour with zero gap between them. A boundary band carved out
of BOTH neighbours turns "two touching objects" into "two interior components separated by
a ridge", which connected components can then split — a purely LOCAL operation, so it
survives tiling. The band is taken *inside* the instances (never grown into background),
which keeps background pixels honestly background.

The band is computed from the LABEL map, not per-instance contours, so an instance-instance
contact and an instance-background contact both produce a band, and the pixels where two
DIFFERENT instances meet are additionally returned as ``touching`` for loss upweighting —
those are the pixels that actually decide whether two tubules become one object or two.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Optional, Sequence

import numpy as np
from scipy.ndimage import maximum_filter, minimum_filter

BG, INTERIOR, BOUNDARY, TOUCHING = 0, 1, 2, 3

_BIG = np.int64(np.iinfo(np.int32).max)


def _disk(radius: int) -> np.ndarray:
    """Boolean disk footprint of the given radius (radius 0 -> single pixel)."""
    r = int(radius)
    if r <= 0:
        return np.ones((1, 1), dtype=bool)
    y, x = np.ogrid[-r:r + 1, -r:r + 1]
    return (x * x + y * y) <= r * r


# --------------------------------------------------------------------- semantic
def semantic_target(
    id_map: np.ndarray,
    class_idx_of_id: dict[int, int],
) -> np.ndarray:
    """(H, W) int64 semantic map; unlisted ids fall back to background."""
    out = np.zeros(id_map.shape, dtype=np.int64)
    if not class_idx_of_id:
        return out
    ids = np.asarray(sorted(class_idx_of_id), dtype=np.int64)
    lut = np.zeros(int(max(ids.max(), id_map.max())) + 1, dtype=np.int64)
    for i in ids:
        lut[i] = class_idx_of_id[int(i)]
    np.clip(id_map, 0, len(lut) - 1, out=out)
    return lut[out]


# -------------------------------------------------------------- boundary/interior
def boundary_target(
    id_map: np.ndarray,
    boundary_width_px: int = 3,
    interior_erode_px: int = 0,
    n_classes: int = 3,
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(boundary_map, touching_mask)``.

    ``boundary_map`` is (H, W) int64. A pixel of an instance is BOUNDARY when a *different*
    label (another instance OR background) lies within ``boundary_width_px``.
    ``touching_mask`` is the subset of BOUNDARY pixels where the competing label is another
    INSTANCE rather than background.

    ``n_classes``:
      3 -> BG / INTERIOR / BOUNDARY. Touching pixels are BOUNDARY; the mask is returned
           separately for per-pixel loss upweighting.
      4 -> BG / INTERIOR / BOUNDARY / TOUCHING, promoting the instance-instance wall to its
           own class. Those pixels are the only thing standing between two touching tubules
           and a merge, and they are a small minority of an already-small boundary class; at
           3 classes they are diluted by the far more numerous object-to-background
           contours, which are easy and which the head can satisfy while ignoring the walls.
    """
    lab = np.ascontiguousarray(id_map.astype(np.int64))
    fg = lab > 0
    out = np.where(fg, INTERIOR, BG).astype(np.int64)
    touching = np.zeros(lab.shape, dtype=bool)
    if not fg.any():
        return out, touching

    fp = _disk(boundary_width_px)
    lab_max = maximum_filter(lab, footprint=fp, mode="nearest")
    lab_min = minimum_filter(lab, footprint=fp, mode="nearest")
    # Any label transition within the footprint, restricted to instance pixels: this is
    # the band lying INSIDE each instance along its contour.
    band = fg & (lab_max != lab_min)

    # Instance-instance contact: ignore background when taking the min, so a pixel is
    # "touching" only when two distinct NONZERO labels are within the footprint.
    lab_nz = np.where(fg, lab, _BIG)
    min_nz = minimum_filter(lab_nz, footprint=fp, mode="nearest")
    touching = band & (min_nz != _BIG) & (lab_max != min_nz)

    out[band] = BOUNDARY

    if interior_erode_px > 0:
        interior = out == INTERIOR
        eroded = minimum_filter(interior.astype(np.uint8), footprint=_disk(interior_erode_px),
                                mode="constant", cval=0).astype(bool)
        out[interior & ~eroded] = BOUNDARY

    if n_classes == 4:
        # Promote the instance-instance wall AFTER the erosion step, so a pixel that became
        # BOUNDARY only through erosion is not mislabelled as a wall between two objects.
        out[touching] = TOUCHING
    elif n_classes != 3:
        raise ValueError(f"boundary target n_classes must be 3 or 4, got {n_classes}")

    return out, touching


def interior_mask(id_map: np.ndarray, boundary_width_px: int = 3) -> np.ndarray:
    """Convenience: boolean interior (instance minus its boundary band)."""
    bmap, _ = boundary_target(id_map, boundary_width_px)
    return bmap == INTERIOR


# ---------------------------------------------------------------------- centers
def center_heatmaps(
    id_map: np.ndarray,
    class_idx_of_id: dict[int, int],
    compact_indices: Sequence[int],
    sigma_scale: float = 0.125,
    sigma_min: float = 2.0,
    sigma_max: float = 24.0,
) -> np.ndarray:
    """(K, H, W) float32 Gaussian center heatmaps, one channel per compact class.

    Sigma scales with object size (``sigma_scale * sqrt(area)``) and is clamped, so a
    glomerulus and a blood cell each get a peak proportioned to themselves rather than one
    global sigma that is too fat for one and too thin for the other. Channels are combined
    with ``max`` so overlapping peaks do not sum past 1.
    """
    H, W = id_map.shape
    K = len(compact_indices)
    heat = np.zeros((K, H, W), dtype=np.float32)
    if K == 0 or not class_idx_of_id:
        return heat
    chan_of_class = {c: k for k, c in enumerate(compact_indices)}

    for inst_id, cls_idx in class_idx_of_id.items():
        k = chan_of_class.get(int(cls_idx))
        if k is None:
            continue
        mask = id_map == inst_id
        area = int(mask.sum())
        if area == 0:
            continue
        ys, xs = np.nonzero(mask)
        cy, cx = float(ys.mean()), float(xs.mean())
        sigma = float(np.clip(sigma_scale * np.sqrt(area), sigma_min, sigma_max))
        _draw_gaussian(heat[k], cy, cx, sigma)
    return heat


def _draw_gaussian(canvas: np.ndarray, cy: float, cx: float, sigma: float) -> None:
    """Max-composite an unnormalized Gaussian (peak 1.0) into ``canvas`` in place."""
    H, W = canvas.shape
    r = int(np.ceil(3.0 * sigma))
    y0, y1 = max(0, int(cy) - r), min(H, int(cy) + r + 1)
    x0, x1 = max(0, int(cx) - r), min(W, int(cx) + r + 1)
    if y0 >= y1 or x0 >= x1:
        return
    yy = np.arange(y0, y1, dtype=np.float32)[:, None]
    xx = np.arange(x0, x1, dtype=np.float32)[None, :]
    g = np.exp(-(((yy - cy) ** 2 + (xx - cx) ** 2) / (2.0 * sigma * sigma)))
    np.maximum(canvas[y0:y1, x0:x1], g, out=canvas[y0:y1, x0:x1])


# ------------------------------------------------------------------- bundle
@dataclass
class PanopticTargets:
    semantic: np.ndarray      # (H, W) int64
    boundary: np.ndarray      # (H, W) int64
    touching: np.ndarray      # (H, W) bool
    centers: np.ndarray       # (K, H, W) float32
    id_map: np.ndarray        # (H, W) int32


def build_targets(
    id_map: np.ndarray,
    class_idx_of_id: dict[int, int],
    compact_indices: Sequence[int],
    boundary_width_px: int = 3,
    interior_erode_px: int = 0,
    sigma_scale: float = 0.125,
    sigma_min: float = 2.0,
    sigma_max: float = 24.0,
    boundary_classes: int = 3,
) -> PanopticTargets:
    """Generate every training target from an instance map in one pass."""
    sem = semantic_target(id_map, class_idx_of_id)
    bnd, touch = boundary_target(id_map, boundary_width_px, interior_erode_px,
                                 n_classes=boundary_classes)
    cen = center_heatmaps(id_map, class_idx_of_id, compact_indices,
                          sigma_scale, sigma_min, sigma_max)
    return PanopticTargets(semantic=sem, boundary=bnd, touching=touch, centers=cen,
                           id_map=id_map.astype(np.int32))


def targets_from_config(id_map: np.ndarray, class_idx_of_id: dict[int, int],
                        compact_indices: Sequence[int], cfg: dict) -> PanopticTargets:
    t = (cfg or {}).get("targets", {}) or {}
    return build_targets(
        id_map, class_idx_of_id, compact_indices,
        boundary_width_px=int(t.get("boundary_width_px", 3)),
        interior_erode_px=int(t.get("interior_erode_px", 0)),
        sigma_scale=float(t.get("center_sigma_scale", 0.125)),
        sigma_min=float(t.get("center_sigma_min", 2.0)),
        sigma_max=float(t.get("center_sigma_max", 24.0)),
        boundary_classes=int((cfg or {}).get("model", {}).get("boundary_classes", 3)),
    )
