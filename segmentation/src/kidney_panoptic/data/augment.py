"""Online augmentation: identical geometry for image and every mask.

Design rules
------------
1. **One geometric transform, applied to everything.** Image uses bilinear, every mask
   uses nearest. The instance map is warped, and the boundary / center / semantic targets
   are regenerated from it afterwards (data/targets.py) — never warped themselves.

2. **Out-of-frame pixels become IGNORED, not background.** A rotation pulls content in
   from outside the extracted patch. The image border is reflected so the network never
   sees an implausible black wedge, but an ``inbounds`` mask records which pixels are real
   and those that are not are stripped from the valid mask. Reflecting the image while
   calling the reflected corner "background" would be straightforward label noise.

3. **Photometric jitter touches the image only.** HED stain jitter is the H&E-appropriate
   one (perturb haematoxylin/eosin concentrations independently) and is applied before the
   coarser HSV/brightness jitter.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import cv2
import numpy as np

try:
    from skimage.color import hed2rgb, rgb2hed
    _HAVE_HED = True
except Exception:  # pragma: no cover
    _HAVE_HED = False


@dataclass
class Sample:
    """A patch and every mask that must follow its geometry."""
    image: np.ndarray        # (H, W, 3) uint8
    id_map: np.ndarray       # (H, W) int32
    tissue: np.ndarray       # (H, W) bool
    bg_poly: np.ndarray      # (H, W) bool
    dense_roi: np.ndarray    # (H, W) bool
    inbounds: np.ndarray     # (H, W) bool — False where content came from outside the patch
    class_idx_of_id: dict


def _warp_masks(sample: Sample, warp) -> Sample:
    """Apply a callable ``warp(arr, interp, border_value)`` to every field."""
    return Sample(
        image=warp(sample.image, cv2.INTER_LINEAR, None),
        id_map=warp(sample.id_map.astype(np.int32), cv2.INTER_NEAREST, 0),
        tissue=warp(sample.tissue.astype(np.uint8), cv2.INTER_NEAREST, 0).astype(bool),
        bg_poly=warp(sample.bg_poly.astype(np.uint8), cv2.INTER_NEAREST, 0).astype(bool),
        dense_roi=warp(sample.dense_roi.astype(np.uint8), cv2.INTER_NEAREST, 0).astype(bool),
        inbounds=warp(sample.inbounds.astype(np.uint8), cv2.INTER_NEAREST, 0).astype(bool),
        class_idx_of_id=sample.class_idx_of_id,
    )


# ------------------------------------------------------------------- geometry
def random_flip_rot90(s: Sample, rng: np.random.Generator, cfg: dict) -> Sample:
    def apply(fn):
        return Sample(fn(s.image), fn(s.id_map), fn(s.tissue), fn(s.bg_poly),
                      fn(s.dense_roi), fn(s.inbounds), s.class_idx_of_id)

    if rng.random() < cfg.get("hflip", 0.5):
        s = apply(lambda a: np.ascontiguousarray(a[:, ::-1]))
    if rng.random() < cfg.get("vflip", 0.5):
        s = apply(lambda a: np.ascontiguousarray(a[::-1]))
    if rng.random() < cfg.get("rot90", 0.5):
        k = int(rng.integers(1, 4))
        s = apply(lambda a, k=k: np.ascontiguousarray(np.rot90(a, k)))
    return s


def random_affine(s: Sample, rng: np.random.Generator, cfg: dict) -> Sample:
    """Rotation + anisotropic scale + shear + translation about the patch centre.

    Anisotropy and shear are anatomically defensible here: a tubule is a tube sectioned at
    an arbitrary obliquity, so its cross-section is an ellipse of varying eccentricity.
    Jittering the aspect ratio therefore generates plausible sections rather than distorted
    ones — which would not be true for, say, nuclei.
    """
    h, w = s.id_map.shape
    max_deg = float(cfg.get("rotate_deg", 180))
    lo, hi = cfg.get("scale", [1.0, 1.0])
    a_lo, a_hi = cfg.get("aspect", [1.0, 1.0])
    shear_deg = float(cfg.get("shear_deg", 0.0))
    tmax = float(cfg.get("translate_px", 0))

    angle = float(rng.uniform(-max_deg, max_deg)) if max_deg else 0.0
    scale = float(rng.uniform(lo, hi))
    aspect = float(np.sqrt(rng.uniform(a_lo, a_hi)))
    shear = np.deg2rad(float(rng.uniform(-shear_deg, shear_deg))) if shear_deg else 0.0
    if angle == 0.0 and scale == 1.0 and aspect == 1.0 and shear == 0.0 and tmax == 0:
        return s

    th = np.deg2rad(angle)
    R = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    S = np.diag([scale * aspect, scale / aspect])
    Sh = np.array([[1.0, np.tan(shear)], [0.0, 1.0]])
    A = R @ S @ Sh
    cx, cy = w / 2.0, h / 2.0
    M = np.zeros((2, 3), np.float64)
    M[:, :2] = A
    M[:, 2] = [cx - (A[0, 0] * cx + A[0, 1] * cy), cy - (A[1, 0] * cx + A[1, 1] * cy)]
    if tmax:
        M[0, 2] += float(rng.uniform(-tmax, tmax))
        M[1, 2] += float(rng.uniform(-tmax, tmax))

    def warp(arr, interp, border_value):
        if border_value is None:  # the image: reflect so no implausible flat wedge appears
            return cv2.warpAffine(arr, M, (w, h), flags=interp,
                                  borderMode=cv2.BORDER_REFLECT_101)
        return cv2.warpAffine(arr, M, (w, h), flags=interp,
                              borderMode=cv2.BORDER_CONSTANT, borderValue=border_value)

    return _warp_masks(s, warp)


def random_elastic(s: Sample, rng: np.random.Generator, cfg: dict) -> Sample:
    """Low-magnitude elastic warp — enough to vary contours, not enough to bend anatomy."""
    h, w = s.id_map.shape
    alpha = float(cfg.get("elastic_alpha", 24.0))
    sigma = float(cfg.get("elastic_sigma", 6.0))
    if alpha <= 0:
        return s
    dx = cv2.GaussianBlur(rng.random((h, w), dtype=np.float32) * 2 - 1, (0, 0), sigma) * alpha
    dy = cv2.GaussianBlur(rng.random((h, w), dtype=np.float32) * 2 - 1, (0, 0), sigma) * alpha
    gx, gy = np.meshgrid(np.arange(w, dtype=np.float32), np.arange(h, dtype=np.float32))
    mx, my = (gx + dx).astype(np.float32), (gy + dy).astype(np.float32)

    def warp(arr, interp, border_value):
        if border_value is None:
            return cv2.remap(arr, mx, my, interp, borderMode=cv2.BORDER_REFLECT_101)
        return cv2.remap(arr, mx, my, interp, borderMode=cv2.BORDER_CONSTANT,
                         borderValue=border_value)

    return _warp_masks(s, warp)


# ---------------------------------------------------------------- photometric
def hed_jitter(image: np.ndarray, rng: np.random.Generator,
               alpha_sigma: float, beta_od: float) -> np.ndarray:
    """Perturb haematoxylin / eosin / DAB concentrations independently (Tellez et al.).

    This is the stain-aware augmentation for H&E: it varies how strongly each dye reads
    without changing tissue geometry, which is what actually differs between scanners and
    staining runs.

    TWO terms on DIFFERENT scales, hence two parameters:

        alpha   multiplicative, dimensionless   ~ U(1 - alpha_sigma, 1 + alpha_sigma)
        beta    ADDITIVE, in skimage OPTICAL-DENSITY units ~ U(-beta_od, +beta_od)

    They were previously a single ``sigma`` with an undocumented ``* 0.05`` on beta, which
    made the additive term look negligible next to the OD channel ranges. It is not:
    ``hed2rgb`` is exponential, so a uniform OD offset produces a systematic RGB shift far
    larger than its size in OD space suggests. Measured mean |dRGB| over six real patches:

        alpha only (0.05)          3.00 / 255
        beta as shipped (0.0025)   3.25 / 255      <- already the larger of the two
        beta without the 0.05     37.07 / 255      <- would destroy the image

    So keep beta in OD units and tune it there. Values around 0.002-0.006 are sane;
    0.05 is not.
    """
    if not _HAVE_HED or (alpha_sigma <= 0 and beta_od <= 0):
        return image
    hed = rgb2hed(image.astype(np.float32) / 255.0)
    alpha = rng.uniform(1 - alpha_sigma, 1 + alpha_sigma, size=3).astype(np.float32)
    beta = rng.uniform(-beta_od, beta_od, size=3).astype(np.float32)
    rgb = np.clip(hed2rgb(hed * alpha[None, None, :] + beta[None, None, :]), 0, 1)
    return (rgb * 255.0).astype(np.uint8)


def _hed_params(cfg: dict) -> tuple[float, float]:
    """Read the split alpha/beta keys, falling back to the legacy single ``hed_sigma``.

    The legacy fallback reproduces the old behaviour exactly (beta was internally scaled by
    0.05), so an old config keeps its meaning instead of silently becoming 20x stronger.
    """
    if "hed_alpha_sigma" in cfg or "hed_beta_od" in cfg:
        return float(cfg.get("hed_alpha_sigma", 0.0)), float(cfg.get("hed_beta_od", 0.0))
    legacy = float(cfg.get("hed_sigma", 0.0))
    return legacy, 0.05 * legacy


def photometric(image: np.ndarray, rng: np.random.Generator, cfg: dict) -> np.ndarray:
    a_sigma, b_od = _hed_params(cfg)
    img = image
    if rng.random() < float(cfg.get("hed_p", 1.0)):
        img = hed_jitter(img, rng, a_sigma, b_od)

    if rng.random() < cfg.get("hsv_p", 0.0):
        hsv = cv2.cvtColor(img, cv2.COLOR_RGB2HSV).astype(np.int16)
        hsv[..., 0] = (hsv[..., 0] + int(rng.integers(-8, 9))) % 180
        hsv[..., 1] = np.clip(hsv[..., 1] + int(rng.integers(-20, 21)), 0, 255)
        hsv[..., 2] = np.clip(hsv[..., 2] + int(rng.integers(-20, 21)), 0, 255)
        img = cv2.cvtColor(hsv.astype(np.uint8), cv2.COLOR_HSV2RGB)

    if rng.random() < cfg.get("brightness_contrast_p", 0.0):
        a = float(rng.uniform(0.85, 1.15))     # contrast
        b = float(rng.uniform(-20, 20))        # brightness
        img = np.clip(img.astype(np.float32) * a + b, 0, 255).astype(np.uint8)

    if rng.random() < cfg.get("blur_p", 0.0):
        lo, hi = cfg.get("blur_sigma", [0.4, 1.2])
        img = cv2.GaussianBlur(img, (0, 0), float(rng.uniform(lo, hi)))

    if rng.random() < cfg.get("noise_p", 0.0):
        lo, hi = cfg.get("noise_sigma", [2.0, 8.0])
        sigma = float(rng.uniform(lo, hi))
        img = np.clip(img.astype(np.float32) + rng.normal(0, sigma, img.shape),
                      0, 255).astype(np.uint8)

    # JPEG is the scanner-realistic corruption: NDPI/SVS tiles are JPEG-compressed, and
    # Ctrl_1A2 arrived through a different export path entirely. Cheap insurance for the
    # out-of-distribution gap rather than a generic "more augmentation" knob.
    if rng.random() < cfg.get("jpeg_p", 0.0):
        lo, hi = cfg.get("jpeg_quality", [55, 90])
        ok, buf = cv2.imencode(".jpg", cv2.cvtColor(img, cv2.COLOR_RGB2BGR),
                               [int(cv2.IMWRITE_JPEG_QUALITY), int(rng.integers(lo, hi + 1))])
        if ok:
            img = cv2.cvtColor(cv2.imdecode(buf, cv2.IMREAD_COLOR), cv2.COLOR_BGR2RGB)

    return img


# ---------------------------------------------------------------- copy-paste
def copy_paste(
    s: Sample,
    donor_image: np.ndarray,
    donor_id_map: np.ndarray,
    donor_inst_id: int,
    donor_class_idx: int,
    rng: np.random.Generator,
    blend_px: int = 3,
    max_overlap_frac: float = 0.15,
) -> Sample:
    """Paste one donor instance onto ``s`` at a random tissue location.

    The paste is rejected (returned unchanged) if it would bury an existing annotated
    instance — destroying a labelled object to add a synthetic one trades a real training
    signal for a fake one. Pasted pixels become a NEW instance id, so the valid mask picks
    them up as trusted foreground automatically.
    """
    src = donor_id_map == donor_inst_id
    if not src.any():
        return s
    ys, xs = np.nonzero(src)
    y0, y1, x0, x1 = ys.min(), ys.max() + 1, xs.min(), xs.max() + 1
    obj_mask = src[y0:y1, x0:x1]
    obj_rgb = donor_image[y0:y1, x0:x1]
    oh, ow = obj_mask.shape
    H, W = s.id_map.shape
    if oh >= H or ow >= W:
        return s

    # Prefer a landing site on tissue that is not already occupied.
    for _ in range(8):
        ty = int(rng.integers(0, H - oh + 1))
        tx = int(rng.integers(0, W - ow + 1))
        dst_ids = s.id_map[ty:ty + oh, tx:tx + ow]
        dst_tissue = s.tissue[ty:ty + oh, tx:tx + ow]
        occupied = ((dst_ids > 0) & obj_mask).sum()
        if occupied <= max_overlap_frac * obj_mask.sum() and dst_tissue[obj_mask].mean() > 0.5:
            break
    else:
        return s

    # Feathered alpha so the seam is not a learnable straight edge.
    alpha = obj_mask.astype(np.float32)
    if blend_px > 0:
        alpha = cv2.GaussianBlur(alpha, (0, 0), float(blend_px)) * obj_mask
    alpha = alpha[..., None]

    image = s.image.copy()
    region = image[ty:ty + oh, tx:tx + ow].astype(np.float32)
    image[ty:ty + oh, tx:tx + ow] = (
        region * (1 - alpha) + obj_rgb.astype(np.float32) * alpha
    ).astype(np.uint8)

    id_map = s.id_map.copy()
    new_id = int(id_map.max()) + 1
    sub = id_map[ty:ty + oh, tx:tx + ow]
    sub[obj_mask] = new_id
    id_map[ty:ty + oh, tx:tx + ow] = sub

    classes = dict(s.class_idx_of_id)
    classes[new_id] = int(donor_class_idx)
    return Sample(image, id_map, s.tissue, s.bg_poly, s.dense_roi, s.inbounds, classes)


def random_cutout(s: Sample, rng: np.random.Generator, cfg: dict) -> Sample:
    """Erase rectangles from the IMAGE and mark them IGNORED via ``inbounds``.

    The masks are deliberately left untouched. Painting a flat rectangle and leaving the
    instance map alone would teach the network that a featureless patch is *background* —
    exactly the label noise design rule 2 exists to prevent, and precisely the wrong lesson
    for a model whose failure mode is merging. Marking the hole ignored means cutout
    regularises the surrounding context without ever supervising the hole itself.
    """
    h, w = s.id_map.shape
    n_lo, n_hi = cfg.get("cutout_n", [1, 3])
    f_lo, f_hi = cfg.get("cutout_frac", [0.02, 0.08])
    image = s.image.copy()
    inbounds = s.inbounds.copy()
    fill = image.reshape(-1, image.shape[-1]).mean(axis=0).astype(np.uint8)

    for _ in range(int(rng.integers(n_lo, n_hi + 1))):
        frac = float(rng.uniform(f_lo, f_hi))
        aspect = float(rng.uniform(0.5, 2.0))
        rh = int(np.clip(np.sqrt(frac * h * w * aspect), 4, h - 1))
        rw = int(np.clip(np.sqrt(frac * h * w / aspect), 4, w - 1))
        y0 = int(rng.integers(0, h - rh + 1))
        x0 = int(rng.integers(0, w - rw + 1))
        image[y0:y0 + rh, x0:x0 + rw] = fill
        inbounds[y0:y0 + rh, x0:x0 + rw] = False

    return Sample(image, s.id_map, s.tissue, s.bg_poly, s.dense_roi, inbounds,
                  s.class_idx_of_id)


# ------------------------------------------------------------------ pipeline
def augment_sample(s: Sample, rng: np.random.Generator, cfg: dict) -> Sample:
    """Full geometric + photometric pipeline (copy-paste is applied by the dataset)."""
    s = random_flip_rot90(s, rng, cfg)
    s = random_affine(s, rng, cfg)
    if rng.random() < cfg.get("elastic_p", 0.0):
        s = random_elastic(s, rng, cfg)
    s = Sample(photometric(s.image, rng, cfg), s.id_map, s.tissue, s.bg_poly,
               s.dense_roi, s.inbounds, s.class_idx_of_id)
    # After the photometric stage on purpose: blurring an erased rectangle would soften its
    # edge into a learnable artefact.
    if rng.random() < cfg.get("cutout_p", 0.0):
        s = random_cutout(s, rng, cfg)
    return s
