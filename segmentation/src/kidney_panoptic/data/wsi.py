"""Whole-slide reader with NATIVE pyramid-level selection by microns-per-pixel.

Opens NDPI/SVS (and pyramidal TIFF) via openslide when possible, falling back to
tifffile, and plain images (PNG/JPG) as a single-level slide. ``target_mpp`` is used
ONLY to PICK the level whose native mpp is most appropriate (``prefer`` breaks ties);
the chosen level is then read at its NATIVE resolution -- there is NO resampling. The
working resolution therefore equals the chosen level's mpp (``working_mpp``), and
``scale_to_full`` is that level's integer downsample factor, giving exact integer
coordinate round-trips to full resolution.

Rationale: when all slides share the same chosen-level mpp (e.g. Level 1 @ 0.44068),
native reading avoids resampling artifacts, keeps exact integer coordinate scaling, and
is faster. Train and infer always read the same native level, so working mpp is
identical by construction. (If a slide's chosen-level mpp deviates, that needs explicit
resampling -- a separate path we intentionally do NOT implement here; callers should run
the cross-slide consistency check, e.g. scripts/check_mpp.py or extract_patches.)

Coordinate spaces
-----------------
* full-res  : level-0 pixel space; QuPath GeoJSON coords live here.
* work      : chosen-level pixel space. ``work = full / scale_to_full`` (integer factor).

Notes
-----
* Anisotropic mpp (mpp_x != mpp_y beyond tolerance): warn and use the mean.
* Non-integer chosen-level downsample: warn and round (coordinate round-trip approximate).

(Vendored unchanged from the v1 `kidneyseg` package, which validated it across the
NDPI cohort; kept self-contained here so kidney_panoptic has no cross-tree import.)
"""
from __future__ import annotations

import warnings
from pathlib import Path
from typing import Optional

import numpy as np

try:
    import openslide
except Exception:  # pragma: no cover - openslide optional at import time
    openslide = None

import tifffile

_OPENSLIDE_EXTS = {".ndpi", ".svs", ".mrxs", ".scn", ".vms", ".vmu", ".bif"}
_PLAIN_EXTS = {".png", ".jpg", ".jpeg", ".bmp"}


def _mpp_from_tiff(path: Path) -> Optional[tuple[float, float]]:
    """Read (mpp_x, mpp_y) from TIFF resolution tags. Returns None if unavailable."""
    try:
        with tifffile.TiffFile(path) as tif:
            page = tif.pages[0]
            tags = page.tags
            if "XResolution" not in tags or "YResolution" not in tags:
                return None
            xr = tags["XResolution"].value
            yr = tags["YResolution"].value
            unit = tags["ResolutionUnit"].value if "ResolutionUnit" in tags else 2
            unit = int(getattr(unit, "value", unit))

            def _to_mpp(res):
                num, den = (res if isinstance(res, tuple) else (res, 1))
                if num == 0:
                    return None
                ppu = num / den  # pixels per unit
                if unit == 3:      # centimeter
                    return 1.0e4 / ppu
                if unit == 2:      # inch
                    return 25400.0 / ppu
                return None        # unitless -> cannot derive microns

            mx, my = _to_mpp(xr), _to_mpp(yr)
            if mx is None or my is None:
                return None
            return float(mx), float(my)
    except Exception:
        return None


class WSIReader:
    """Reads a slide at a NATIVE pyramid level (no resampling).

    By default ``level=1`` is read directly (these slides always carry the working
    image at level 1 with a shared mpp). ``working_mpp`` equals that level's native mpp
    and ``scale_to_full`` is its integer downsample factor, used to round-trip
    coordinates to full-res. Set ``level=None`` to instead pick by ``target_mpp`` +
    ``prefer`` (scanner-agnostic fallback).
    """

    def __init__(
        self,
        path: str | Path,
        level: Optional[int] = 1,
        target_mpp: float = 0.5,
        manual_mpp: Optional[float] = None,
        anisotropy_tol: float = 0.02,
        prefer: str = "finer",
    ):
        self.path = Path(path)
        if not self.path.exists():
            raise FileNotFoundError(f"Slide not found: {self.path}")
        self.target_mpp = float(target_mpp)
        self.anisotropy_tol = anisotropy_tol
        if prefer not in {"finer", "coarser"}:
            raise ValueError(f"prefer must be 'finer' or 'coarser', got '{prefer}'")
        self.prefer = prefer

        self._osr = None  # openslide handle
        self._level_mpps: list[float] = []
        self._level_dims: list[tuple[int, int]] = []  # (w, h) per level
        self._level_downsamples: list[float] = []

        self._open(manual_mpp)

        # Read the chosen level NATIVELY -- no resampling. A fixed `level` (default 1) is
        # used directly; otherwise the level is picked by mpp. working_mpp = level mpp.
        if level is not None:
            n = len(self._level_mpps)
            if not (0 <= level < n):
                raise ValueError(
                    f"{self.path.name}: requested level {level} but slide has {n} level(s) "
                    f"(mpps={[round(m, 5) for m in self._level_mpps]}). "
                    f"Set wsi.level: null to pick by target_mpp, or check the slide pyramid."
                )
            self.level = level
        else:
            self.level = self.level_for_mpp(self.target_mpp)
        self.working_mpp = self._level_mpps[self.level]

        # Integer downsample of the chosen level relative to full-res (level 0); used to
        # round-trip coordinates exactly (full = work * scale_to_full).
        ds = self._level_downsamples[self.level]
        self.scale_to_full = int(round(ds))
        if abs(ds - self.scale_to_full) > 1e-3:
            warnings.warn(
                f"{self.path.name}: chosen level downsample {ds:.4f} is not integer; "
                f"rounding to {self.scale_to_full}. Coordinate round-trip may be approximate.",
                stacklevel=2,
            )
        # Multiplicative full-res -> work factor (work = full / scale_to_full).
        self.scale_full_to_work = 1.0 / self.scale_to_full

    # ------------------------------------------------------------------ open
    def _open(self, manual_mpp: Optional[float]) -> None:
        ext = self.path.suffix.lower()
        # openslide for true vendor WSI formats (NDPI/SVS/...); generic & pyramidal TIFF
        # go through tifffile, which correctly exposes subifd pyramid levels (openslide's
        # generic-tiff reader often collapses them to a single level).
        use_openslide = openslide is not None and ext in _OPENSLIDE_EXTS

        if use_openslide:
            if self._osr is None:
                self._osr = openslide.OpenSlide(str(self.path))
            self._level_dims = [tuple(d) for d in self._osr.level_dimensions]
            self._level_downsamples = list(self._osr.level_downsamples)
            base_mpp = self._read_openslide_mpp(manual_mpp)
            self.base_mpp = base_mpp
            self._level_mpps = [base_mpp * ds for ds in self._level_downsamples]
            self._backend = "openslide"
            return

        # tifffile / plain image path
        self._backend = "tifffile"
        self._tiff_levels = self._load_tiff_levels(ext)
        self._level_dims = [(arr.shape[1], arr.shape[0]) for arr in self._tiff_levels]
        self._level_downsamples = [
            self._level_dims[0][0] / dims[0] for dims in self._level_dims
        ]
        base_mpp = self._read_tiff_mpp(manual_mpp, ext)
        self.base_mpp = base_mpp
        self._level_mpps = [base_mpp * ds for ds in self._level_downsamples]

    def _load_tiff_levels(self, ext: str) -> list[np.ndarray]:
        if ext in _PLAIN_EXTS:
            import cv2

            img = cv2.imread(str(self.path), cv2.IMREAD_COLOR)
            if img is None:
                raise RuntimeError(f"Failed to read image: {self.path}")
            return [cv2.cvtColor(img, cv2.COLOR_BGR2RGB)]
        # Pyramidal TIFF: collect series[0] levels if present.
        with tifffile.TiffFile(self.path) as tif:
            series = tif.series[0]
            levels = []
            if getattr(series, "levels", None):
                for lvl in series.levels:
                    levels.append(self._rgb(lvl.asarray()))
            else:
                levels.append(self._rgb(series.asarray()))
        return levels

    @staticmethod
    def _rgb(arr: np.ndarray) -> np.ndarray:
        arr = np.asarray(arr)
        if arr.ndim == 2:
            arr = np.stack([arr] * 3, axis=-1)
        elif arr.ndim == 3 and arr.shape[2] >= 4:
            arr = arr[:, :, :3]
        elif arr.ndim == 3 and arr.shape[0] in (1, 3) and arr.shape[2] not in (1, 3):
            arr = np.moveaxis(arr, 0, -1)[:, :, :3]
        if arr.dtype != np.uint8:
            mx = arr.max() if arr.size else 1
            arr = (arr / mx * 255).astype(np.uint8) if mx > 255 else arr.astype(np.uint8)
        return np.ascontiguousarray(arr)

    # ----------------------------------------------------------------- mpp
    def _check_anisotropy(self, mx: float, my: float) -> float:
        if abs(mx - my) / max(mx, my, 1e-9) > self.anisotropy_tol:
            warnings.warn(
                f"Anisotropic mpp (x={mx:.4f}, y={my:.4f}) on {self.path.name}; "
                f"using mean.", stacklevel=3,
            )
        return 0.5 * (mx + my)

    def _read_openslide_mpp(self, manual_mpp: Optional[float]) -> float:
        props = self._osr.properties
        mx = props.get(openslide.PROPERTY_NAME_MPP_X)
        my = props.get(openslide.PROPERTY_NAME_MPP_Y)
        if mx and my:
            return self._check_anisotropy(float(mx), float(my))
        if manual_mpp is not None:
            warnings.warn(f"No mpp in {self.path.name}; using manual_mpp={manual_mpp}.", stacklevel=3)
            return float(manual_mpp)
        raise ValueError(
            f"Could not read mpp from {self.path.name} and no manual_mpp provided. "
            f"Set wsi.manual_mpp in config."
        )

    # Plausible range for a brightfield microscope pixel at 5x-100x. A value outside this
    # did not come from a microscope, so it must not be used as one.
    _PLAUSIBLE_MPP = (0.05, 5.0)

    def _read_tiff_mpp(self, manual_mpp: Optional[float], ext: str) -> float:
        # An explicit manual_mpp is caller intent and wins outright -- many TIFFs (esp.
        # ones exported by generic tools) carry a print-resolution tag (e.g. 300 DPI)
        # that has nothing to do with the microscope's actual pixel size and would
        # otherwise silently override a correct manual_mpp.
        if manual_mpp is not None:
            return float(manual_mpp)
        if ext not in _PLAIN_EXTS:
            mpp = _mpp_from_tiff(self.path)
            if mpp is not None:
                value = self._check_anisotropy(*mpp)
                lo, hi = self._PLAUSIBLE_MPP
                if not (lo <= value <= hi):
                    # A 300 DPI print tag yields ~84.7 um/px. Accepting it would silently
                    # scale every reported area by ~37000x, and nothing downstream would
                    # notice -- so refuse rather than guess.
                    raise ValueError(
                        f"{self.path.name}: TIFF resolution tags imply {value:.2f} um/px, "
                        f"outside the plausible microscope range {lo}-{hi}. This is almost "
                        f"certainly a print-resolution tag (e.g. 300 DPI), not a pixel size. "
                        f"Set wsi.manual_mpp explicitly, e.g. --set wsi.manual_mpp=0.44068"
                    )
                return value
        raise ValueError(
            f"Could not read mpp from {self.path.name} and no manual_mpp provided. "
            f"Set wsi.manual_mpp in config."
        )

    # ------------------------------------------------------------- geometry
    def level_for_mpp(self, target_mpp: float) -> int:
        """Pick the level whose native mpp is closest to target_mpp (no resampling).

        On a near-tie between two levels, ``prefer`` decides: 'finer' picks the
        smaller-mpp level, 'coarser' the larger. Selection is by mpp, never by index,
        so it stays scanner-agnostic.
        """
        diffs = [abs(m - target_mpp) for m in self._level_mpps]
        best = min(diffs)
        # Levels within a small relative epsilon of the best are treated as ties.
        eps = 1e-6 + 0.02 * target_mpp
        tied = [i for i, d in enumerate(diffs) if d <= best + eps]
        if self.prefer == "finer":
            return min(tied, key=lambda i: self._level_mpps[i])
        return max(tied, key=lambda i: self._level_mpps[i])

    @property
    def level_dimensions_work(self) -> tuple[int, int]:
        """(width, height) of the chosen level (== working resolution, native)."""
        return self._level_dims[self.level]

    def full_to_work(self, coords: np.ndarray) -> np.ndarray:
        """Full-res (level-0) coords -> work coords via the integer downsample."""
        return np.asarray(coords, dtype=np.float64) / self.scale_to_full

    def work_to_full(self, coords: np.ndarray) -> np.ndarray:
        return np.asarray(coords, dtype=np.float64) * self.scale_to_full

    # --------------------------------------------------------------- reading
    def read_region_work(self, x: int, y: int, w: int, h: int) -> np.ndarray:
        """Read a (h, w, 3) uint8 region from the chosen level, NATIVELY (no resample).

        (x, y, w, h) are in WORK pixels (== chosen-level pixels).
        """
        if self._backend == "openslide":
            # openslide read_region location is in LEVEL-0 coords.
            lx, ly = x * self.scale_to_full, y * self.scale_to_full
            region = self._osr.read_region((lx, ly), self.level, (w, h)).convert("RGB")
            return np.ascontiguousarray(np.asarray(region))

        arr = self._tiff_levels[self.level]
        tile = arr[y:y + h, x:x + w]
        if tile.shape[0] != h or tile.shape[1] != w:
            pad = np.zeros((h, w, 3), np.uint8)
            pad[: tile.shape[0], : tile.shape[1]] = tile
            tile = pad
        return np.ascontiguousarray(tile)

    def read_full_work(self) -> np.ndarray:
        w, h = self.level_dimensions_work
        return self.read_region_work(0, 0, w, h)

    def close(self) -> None:
        if self._osr is not None:
            self._osr.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
