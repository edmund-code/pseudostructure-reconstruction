"""Tissue vs. slide-background (glass) segmentation for a working-resolution RGB patch.

Why this exists
---------------
Under ``dilated_gt_two_radius`` the trusted region is ~2% of the tissue (measured on three
4096^2 regions of Nx_4wk_050: 1.2-2.3%). Empty glass falls in the ignored 98%, so nothing
ever pushes the FG head toward zero there. The measured consequence: fed a CONSTANT image
the head still emits 42-44 saturated instances covering ~26% of the tile, and the pattern is
identical (Pearson 0.993-1.000) across input levels from gray 200 to 255 — a pure decoder
prior with no input dependence (``scripts/diag_decoder_prior.py``). On production slides
17-26% of all predicted instances sit on blank glass.

Glass is trivially and reliably identifiable without any annotation, so it can be supervised
as background for free. That is what this module enables.

Method
------
H&E stain is chromatic; glass is not. Saturation separates them almost perfectly:

    tissue = (S > sat_threshold) OR (gray < gray_threshold)

The gray term catches dark, low-saturation material (debris, folds, pen) that should not be
called background. Morphological closing then opening removes speckle, and small holes are
filled so lumina inside tissue are not mislabelled as glass.

Everything is configurable under ``annotation.tissue`` and the defaults are deliberately
CONSERVATIVE: it is much worse to call tissue "glass" (which would train the FG head to
suppress real structures) than to call glass "tissue" (which merely reverts to the current
ignore behaviour). Hence a low saturation threshold and a dilation of the tissue region.

(Vendored unchanged from the v1 `kidneyseg` package, which validated it across the
NDPI cohort; kept self-contained here so kidney_panoptic has no cross-tree import.)
"""
from __future__ import annotations

import cv2
import numpy as np

# Calibrated on Nx_4wk_050 level 1 (0.44068 µm/px), 1024² crops of pure glass and of
# cortex (scripts/diag_ignore_coverage.py records the same statistics per run):
#
#            saturation                     gray
#   glass    p50 0.035  p95 0.059  p99 0.071    p05 231  p50 234  p95 237
#   tissue   p05 0.176  p50 0.447  p95 0.537    p05  79  p50 136  p95 205
#
# The two populations are cleanly separated on both channels, so the thresholds sit at the
# midpoint of each gap: saturation 0.12 (between 0.071 and 0.176), gray 215 (between 205
# and 231). Re-derive these if the stain or scanner changes.
DEFAULTS = {
    "method": "saturation",   # saturation | otsu
    "sat_threshold": 0.12,    # HSV S in [0, 1]; above this is tissue
    "gray_threshold": 215,    # anything darker than this is tissue regardless of hue
    # OPEN BEFORE CLOSE, and the order is load-bearing. Glass carries ~5% of pixels above
    # the saturation threshold as JPEG/scanner colour noise. Closing first bridges that
    # speckle into a solid blob and swallows the whole glass region as "tissue" — measured
    # on ablation region 0 as 0.7% glass detected against ~28% actual. Opening first
    # deletes the speckle, and only then does closing fill genuine gaps inside tissue.
    "open_px": 5,
    "close_px": 9,
    "min_object_px": 4096,    # drop tissue islands smaller than this (px)
    "fill_holes_px": 65536,   # fill ENCLOSED glass-labelled holes inside tissue up to this
    "dilate_px": 16,          # grow tissue outward; keeps the tissue EDGE out of "trusted bg"
}


def tissue_params(cfg: dict | None) -> dict:
    p = dict(DEFAULTS)
    p.update({k: v for k, v in (cfg or {}).items() if k in DEFAULTS})
    return p


def _disk(r: int) -> np.ndarray:
    return cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * int(r) + 1, 2 * int(r) + 1))


def build_tissue_mask(rgb: np.ndarray, cfg: dict | None = None) -> np.ndarray:
    """Return (H, W) bool — True = tissue, False = glass / slide background."""
    p = tissue_params(cfg)
    if rgb.ndim != 3 or rgb.shape[2] < 3:
        raise ValueError(f"build_tissue_mask expects an (H, W, 3) RGB array, got {rgb.shape}")
    img = np.ascontiguousarray(rgb[:, :, :3].astype(np.uint8))

    hsv = cv2.cvtColor(img, cv2.COLOR_RGB2HSV)
    sat = hsv[:, :, 1].astype(np.float32) / 255.0
    gray = cv2.cvtColor(img, cv2.COLOR_RGB2GRAY)

    if p["method"] == "otsu":
        thr, _ = cv2.threshold((sat * 255).astype(np.uint8), 0, 255,
                               cv2.THRESH_BINARY + cv2.THRESH_OTSU)
        sat_hit = sat > (thr / 255.0)
    else:
        sat_hit = sat > float(p["sat_threshold"])

    mask = (sat_hit | (gray < int(p["gray_threshold"]))).astype(np.uint8)

    # Open before close — see the note on DEFAULTS; reversing these silently loses glass.
    if p["open_px"] > 0:
        mask = cv2.morphologyEx(mask, cv2.MORPH_OPEN, _disk(p["open_px"]))
    if p["close_px"] > 0:
        mask = cv2.morphologyEx(mask, cv2.MORPH_CLOSE, _disk(p["close_px"]))

    m = mask.astype(bool)
    if p["min_object_px"] > 0:
        m = _drop_small(m, int(p["min_object_px"]))
    if p["fill_holes_px"] > 0:
        m = _fill_enclosed_holes(m, int(p["fill_holes_px"]))
    if p["dilate_px"] > 0:
        m = cv2.dilate(m.astype(np.uint8), _disk(p["dilate_px"])).astype(bool)
    return m


def _drop_small(mask: np.ndarray, min_px: int) -> np.ndarray:
    """Remove connected components smaller than ``min_px``."""
    if not mask.any():
        return mask
    n, lab, stats, _ = cv2.connectedComponentsWithStats(mask.astype(np.uint8), 8)
    if n <= 1:
        return mask
    keep = np.zeros(n, dtype=bool)
    keep[1:] = stats[1:, cv2.CC_STAT_AREA] >= min_px
    return keep[lab]


def _fill_enclosed_holes(mask: np.ndarray, max_px: int) -> np.ndarray:
    """Fill background components fully ENCLOSED by tissue, up to ``max_px``.

    Components touching the patch border are not holes — they are the slide background
    continuing off-frame. Filling those would swallow the entire glass region of any
    patch straddling the tissue edge, which is precisely the patch this whole mechanism
    exists to supervise.
    """
    if not mask.any():
        return mask
    bg = ~mask
    n, lab, stats, _ = cv2.connectedComponentsWithStats(bg.astype(np.uint8), 8)
    if n <= 1:
        return mask
    h, w = mask.shape
    border = np.unique(np.concatenate([lab[0], lab[-1], lab[:, 0], lab[:, -1]]))
    fill = np.zeros(n, dtype=bool)
    for i in range(1, n):
        if i in border:
            continue
        if stats[i, cv2.CC_STAT_AREA] <= max_px:
            fill[i] = True
    return mask | fill[lab]


def glass_fraction(rgb: np.ndarray, cfg: dict | None = None) -> float:
    """Fraction of the patch that is slide background."""
    return float(1.0 - build_tissue_mask(rgb, cfg).mean())
