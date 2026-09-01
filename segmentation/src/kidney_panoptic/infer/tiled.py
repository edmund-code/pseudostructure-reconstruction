"""Tiled inference: blended probability maps, windowed decode, stitched instances.

Two nested tilings, and they do different jobs:

* **Prediction tiles** (``infer.tile`` / ``infer.overlap``) produce the probability maps.
  Overlapping tiles are combined with a Gaussian window so a tile edge contributes least
  where it is least reliable. Only genuinely LOCAL fields are blended — semantic, boundary
  and center probabilities. Averaging a non-local field (a centroid-offset vector, say)
  across two tiles that disagree about which object a pixel belongs to would invent values
  neither tile predicted; that failure mode is designed out rather than mitigated.

* **Decode windows** (``infer.decode_window`` / ``decode_overlap``) run the instance decode
  on a large area at once, then hand the result to the stitcher.

Exactness of the blend
----------------------
Prediction tiles sit on a single slide-global grid. To decode a window, the tiles that
intersect it *plus a halo of one full tile* are run, so every pixel in the window receives
the exact same set of contributions it would have received from one slide-wide pass. The
halo costs ~2x the forward passes on a 2048 window and buys an exactly seam-free map,
which at this slide size is a few minutes of GPU time.
"""
from __future__ import annotations

from typing import Callable, Iterator, Optional

import numpy as np
import torch

from ..data.dataset import normalize_image
from ..data.classes import ClassSpec
from ..postprocess.decode import InstanceResult, decode_panoptic

FetchFn = Callable[[int, int, int, int], np.ndarray]


def gaussian_window(size: int, sigma_scale: float = 0.25) -> np.ndarray:
    """2D Gaussian tile weight, strictly positive so no pixel can end up unweighted."""
    sigma = max(size * float(sigma_scale), 1e-3)
    c = (size - 1) / 2.0
    ax = np.arange(size, dtype=np.float32) - c
    g = np.exp(-(ax ** 2) / (2 * sigma ** 2))
    w = np.outer(g, g).astype(np.float32)
    return np.maximum(w, 1e-4)


def dihedral_transforms(mode: str) -> list[tuple[int, bool]]:
    """(k, flip) pairs for test-time augmentation. ``k`` = quarter turns, then a horizontal flip.

    ``none`` -> identity only | ``rot4`` -> the 4 rotations | ``d4`` -> the full 8-element group.

    This is only sound because every predicted field is a spatially-equivariant SCALAR per
    channel: semantic probability, boundary/interior probability, center heatmaps. Rotating
    the input rotates each map and nothing else. A centroid-offset (HV) field would NOT
    survive this — its vector components would need rotating too, and averaging two
    disagreeing offset fields invents values neither pass predicted. Choosing the
    boundary+center representation (ADR 0001) is what makes TTA free here.
    """
    if mode in (None, "", "none"):
        return [(0, False)]
    if mode == "rot4":
        return [(k, False) for k in range(4)]
    if mode == "d4":
        return [(k, f) for f in (False, True) for k in range(4)]
    raise ValueError(f"Unknown infer.tta {mode!r}; expected none | rot4 | d4")


def _tta_forward(model, x: torch.Tensor, amp: bool, transforms) -> dict[str, torch.Tensor]:
    """Average ``model.predict_maps`` over the given dihedral transforms, in input space."""
    if len(transforms) == 1 and transforms[0] == (0, False):
        return model.predict_maps(x, amp=amp)
    acc: Optional[dict[str, torch.Tensor]] = None
    for k, flip in transforms:
        xt = torch.flip(x, dims=(-1,)) if flip else x
        xt = torch.rot90(xt, k, dims=(-2, -1)) if k else xt
        out = model.predict_maps(xt, amp=amp)
        for name, v in out.items():
            v = torch.rot90(v, -k, dims=(-2, -1)) if k else v
            v = torch.flip(v, dims=(-1,)) if flip else v
            if acc is None:
                acc = {}
            acc[name] = v if name not in acc else acc[name] + v
    n = float(len(transforms))
    return {name: v / n for name, v in acc.items()}


def tile_origins(start: int, stop: int, size: int, step: int, limit: int) -> list[int]:
    """Globally-aligned tile origins (multiples of ``step``) whose tile meets [start, stop)."""
    first = max(0, ((start - size + 1) + step - 1) // step)
    origins = []
    x = first * step
    while x < stop:
        o = min(x, max(0, limit - size))
        if o + size > start:
            origins.append(o)
        if x >= limit - size:
            break
        x += step
    return sorted(set(origins))


class TiledPredictor:
    """Runs the model over a region and returns blended probability maps."""

    def __init__(self, model, cfg: dict, spec: ClassSpec, device: str = "cuda",
                 amp: bool = True):
        self.model = model.to(device).eval()
        self.cfg = cfg
        self.spec = spec
        self.device = device
        self.amp = amp
        icfg = dict(cfg.get("infer", {}) or {})
        self.tile = int(icfg.get("tile", 512))
        self.overlap = int(icfg.get("overlap", 128))
        self.step = max(1, self.tile - self.overlap)
        self.batch_size = int(icfg.get("batch_size", 8))
        self.window = gaussian_window(self.tile, float(icfg.get("blend_sigma_scale", 0.25)))
        self.n_classes = spec.n_classes
        self.n_compact = len(spec.compact_indices)
        # Read the boundary width from the MODEL, never assume it. The head is 3 or 4
        # channels depending on model.boundary_classes, and hardcoding 3 here silently
        # mis-shapes the accumulator the moment that config changes.
        self.n_boundary = int(getattr(model, "boundary_classes", 3))
        # Test-time augmentation. Costs one forward pass per transform and buys sharper,
        # less noisy boundary probability -- which is what decides whether two touching
        # tubules separate. Off by default; `infer.tta: d4` is the 8x setting.
        self.tta = dihedral_transforms(icfg.get("tta", "none"))

    # ------------------------------------------------------------------ core
    @torch.no_grad()
    def predict_region(
        self,
        fetch: FetchFn,
        x0: int,
        y0: int,
        w: int,
        h: int,
        slide_w: int,
        slide_h: int,
        halo: Optional[int] = None,
        tissue_fetch: Optional[FetchFn] = None,
    ) -> dict[str, np.ndarray]:
        """Blended maps over [x0, x0+w) x [y0, y0+h) in slide (work) coordinates.

        ``halo`` defaults to one full tile, which makes the result identical to a slide-wide
        blended pass. ``tissue_fetch`` (optional) lets whole tiles of glass be skipped.
        """
        halo = self.tile if halo is None else int(halo)
        xs = tile_origins(x0 - halo, min(x0 + w + halo, slide_w), self.tile, self.step, slide_w)
        ys = tile_origins(y0 - halo, min(y0 + h + halo, slide_h), self.tile, self.step, slide_h)

        C, B, K = self.n_classes, self.n_boundary, self.n_compact
        acc = np.zeros((C + B + K, h, w), np.float32)
        wsum = np.zeros((1, h, w), np.float32)

        batch_imgs: list[np.ndarray] = []
        batch_pos: list[tuple[int, int]] = []

        def flush():
            if not batch_imgs:
                return
            x = torch.from_numpy(np.stack(batch_imgs)).to(self.device, non_blocking=True)
            out = _tta_forward(self.model, x, self.amp, self.tta)
            maps = torch.cat([out["semantic_prob"], out["boundary_prob"],
                              out["center_heatmaps"]], dim=1).cpu().numpy()
            for (tx, ty), m in zip(batch_pos, maps):
                self._accumulate(acc, wsum, m, tx - x0, ty - y0, w, h)
            batch_imgs.clear()
            batch_pos.clear()

        for ty in ys:
            for tx in xs:
                if tissue_fetch is not None:
                    if not tissue_fetch(tx, ty, self.tile, self.tile).any():
                        continue
                rgb = fetch(tx, ty, self.tile, self.tile)
                batch_imgs.append(normalize_image(rgb))
                batch_pos.append((tx, ty))
                if len(batch_imgs) >= self.batch_size:
                    flush()
        flush()

        np.maximum(wsum, 1e-6, out=wsum)
        acc /= wsum
        return {
            "semantic_prob": acc[:C],
            "boundary_prob": acc[C:C + B],
            "center_heatmaps": acc[C + B:],
            "weight": wsum[0],
        }

    def _accumulate(self, acc, wsum, maps, ox, oy, w, h):
        """Add one tile's weighted maps into the region accumulator, clipped to the region."""
        t = self.tile
        sx0, sy0 = max(0, -ox), max(0, -oy)
        dx0, dy0 = max(0, ox), max(0, oy)
        cw = min(t - sx0, w - dx0)
        ch = min(t - sy0, h - dy0)
        if cw <= 0 or ch <= 0:
            return
        win = self.window[sy0:sy0 + ch, sx0:sx0 + cw]
        acc[:, dy0:dy0 + ch, dx0:dx0 + cw] += maps[:, sy0:sy0 + ch, sx0:sx0 + cw] * win
        wsum[0, dy0:dy0 + ch, dx0:dx0 + cw] += win

    # ---------------------------------------------------------------- decode
    def predict_and_decode(
        self,
        fetch: FetchFn,
        x0: int,
        y0: int,
        w: int,
        h: int,
        slide_w: int,
        slide_h: int,
        tissue: Optional[np.ndarray] = None,
        tissue_fetch: Optional[FetchFn] = None,
    ) -> tuple[InstanceResult, dict]:
        maps = self.predict_region(fetch, x0, y0, w, h, slide_w, slide_h,
                                   tissue_fetch=tissue_fetch)
        res = decode_panoptic(
            maps["semantic_prob"], maps["boundary_prob"], maps["center_heatmaps"],
            self.cfg, self.spec, tissue_mask=tissue,
        )
        return res, maps


def array_fetcher(arr: np.ndarray) -> FetchFn:
    """Fetch function over an in-memory RGB array, zero-padded past the edges."""
    H, W = arr.shape[:2]

    def fetch(x: int, y: int, w: int, h: int) -> np.ndarray:
        out = np.zeros((h, w) + arr.shape[2:], dtype=arr.dtype)
        xs, ys = max(0, x), max(0, y)
        xe, ye = min(W, x + w), min(H, y + h)
        if xe > xs and ye > ys:
            out[ys - y:ye - y, xs - x:xe - x] = arr[ys:ye, xs:xe]
        return out

    return fetch


def reader_fetcher(reader) -> FetchFn:
    """Fetch function over a WSIReader at its working level."""
    def fetch(x: int, y: int, w: int, h: int) -> np.ndarray:
        return reader.read_region_work(int(x), int(y), int(w), int(h))
    return fetch


def decode_windows(width: int, height: int, window: int, overlap: int) -> Iterator[tuple[int, int, int, int]]:
    """Yield (x, y, w, h) decode windows covering the slide, overlapping by ``overlap``."""
    step = max(1, window - overlap)

    def axis(extent):
        last = max(0, extent - window)
        vals = list(range(0, last + 1, step))
        if not vals or vals[-1] != last:
            vals.append(last)
        return sorted(set(vals))

    for y in axis(height):
        for x in axis(width):
            yield x, y, min(window, width - x), min(window, height - y)
