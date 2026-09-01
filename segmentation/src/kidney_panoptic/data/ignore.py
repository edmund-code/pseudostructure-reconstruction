"""Three-way supervision mask for partially annotated slides.

On a WSI only a fraction of the real objects are drawn. Treating every unannotated pixel
as background teaches the model to suppress real structures; treating every unannotated
pixel as "unknown" leaves the foreground head with nothing pushing it toward zero on empty
glass, and it then hallucinates instances on blank background.

So the mask is three-way:

    TRUSTED FOREGROUND  annotated instance pixels
    TRUSTED BACKGROUND  a thin ring around each annotated instance
                        + everything off tissue (glass needs no annotation to be known
                          background)
                        + explicitly annotated `background_names` polygons
                        + everything inside a dense ROI that is not an instance
    IGNORE              tissue without annotation — real structures may live there

``build_valid_mask`` returns the boolean VALID region (trusted FG ∪ trusted BG); the
complement is ignored in every loss.

Strategies
----------
``dilated_gt_two_radius_plus_tissue``  (default) all three trusted sources above.
``dilated_gt_two_radius``              same, minus the off-tissue term.
``dilated_gt``                         legacy single wide collar around instances.
``dense_roi``                          only dense ROIs are supervised (fully-labelled data).
``none``                               every pixel valid — correct ONLY for dense labels.

Evaluation uses the SAME mask. A prediction on glass is a genuine false positive and is
scored as one; that is the honest reading and it is what makes the glass-rate metric and
PQ tell a consistent story.
"""
from __future__ import annotations

from typing import Optional

import cv2
import numpy as np

STRATEGIES = (
    "dilated_gt_two_radius_plus_tissue",
    "dilated_gt_two_radius",
    "dilated_gt",
    "dense_roi",
    "none",
)


def _disk(r: int) -> np.ndarray:
    return cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * int(r) + 1, 2 * int(r) + 1))


def _dilate(mask: np.ndarray, r: int) -> np.ndarray:
    if r <= 0:
        return mask.copy()
    return cv2.dilate(mask.astype(np.uint8), _disk(r)).astype(bool)


def build_valid_mask(
    id_map: np.ndarray,
    *,
    strategy: str = "dilated_gt_two_radius_plus_tissue",
    tissue_mask: Optional[np.ndarray] = None,
    background_mask: Optional[np.ndarray] = None,
    dense_roi_mask: Optional[np.ndarray] = None,
    dilate_px_inner: int = 5,
    dilate_px: int = 30,
) -> np.ndarray:
    """Boolean (H, W) VALID mask; False = ignored in the loss.

    Parameters
    ----------
    id_map          instance ids (0 = background).
    tissue_mask     True = tissue. Required by ``*_plus_tissue``.
    background_mask True inside annotated background polygons (trusted background).
    dense_roi_mask  True inside fully-annotated ROIs (everything there is trusted).
    """
    if strategy not in STRATEGIES:
        raise ValueError(f"Unknown ignore strategy '{strategy}'. Known: {list(STRATEGIES)}")

    shape = id_map.shape
    fg = id_map > 0

    if strategy == "none":
        return np.ones(shape, dtype=bool)

    if strategy == "dense_roi":
        if dense_roi_mask is None:
            raise ValueError("strategy 'dense_roi' requires dense_roi_mask")
        return dense_roi_mask.astype(bool) | fg

    if strategy == "dilated_gt":
        valid = _dilate(fg, dilate_px)
    else:  # both two_radius variants
        valid = _dilate(fg, dilate_px_inner)

    if strategy == "dilated_gt_two_radius_plus_tissue":
        if tissue_mask is None:
            raise ValueError(
                "strategy 'dilated_gt_two_radius_plus_tissue' requires tissue_mask. "
                "Build one with data.tissue.build_tissue_mask(rgb, cfg['tissue'])."
            )
        # Off-tissue is trusted BACKGROUND. Annotated instances always win, so an
        # instance sitting on a tissue-mask miss is never relabelled to background.
        valid = valid | (~tissue_mask.astype(bool))

    if background_mask is not None:
        valid = valid | background_mask.astype(bool)
    if dense_roi_mask is not None:
        valid = valid | dense_roi_mask.astype(bool)

    return valid | fg


def valid_fraction(valid: np.ndarray) -> float:
    return float(valid.mean())


def trusted_bg_fraction(valid: np.ndarray, id_map: np.ndarray) -> float:
    """Fraction of the patch that is trusted BACKGROUND (valid and not an instance)."""
    return float((valid & (id_map == 0)).mean())


def apply_ignore_to_semantic(semantic: np.ndarray, valid: np.ndarray,
                             ignore_index: int = 255) -> np.ndarray:
    """Stamp ``ignore_index`` into a semantic target wherever the mask is invalid."""
    out = semantic.astype(np.int64).copy()
    out[~valid.astype(bool)] = ignore_index
    return out
