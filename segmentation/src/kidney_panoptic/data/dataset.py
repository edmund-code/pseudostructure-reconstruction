"""Torch dataset producing panoptic targets from stored patches.

Each item::

    image     (3, H, W)  float32, ImageNet-normalized
    semantic  (H, W)     int64, 255 = ignore
    boundary  (H, W)     int64, 255 = ignore  (0 bg / 1 interior / 2 boundary)
    centers   (K, H, W)  float32
    valid     (H, W)     float32 {0, 1}
    weight    (H, W)     float32  per-pixel loss weight (touching pixels upweighted)
    id_map    (H, W)     int32    kept for instance-level evaluation
    meta      dict

Order of operations matters and is deliberate:
    load -> copy-paste -> geometry+photometric augmentation -> RELABEL id_map ->
    build valid mask -> generate targets.
Targets are generated LAST, from the warped instance map, so no target can disagree with
the image it supervises.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import torch
from torch.utils.data import Dataset

from .augment import Sample, augment_sample, copy_paste
from .classes import ClassSpec
from .ignore import apply_ignore_to_semantic, build_valid_mask
from .patching import load_patch
from .targets import targets_from_config

IGNORE_INDEX = 255
_MEAN = np.array([0.485, 0.456, 0.406], np.float32)
_STD = np.array([0.229, 0.224, 0.225], np.float32)


def normalize_image(rgb: np.ndarray) -> np.ndarray:
    """(H, W, 3) uint8 -> (3, H, W) float32, ImageNet-normalized."""
    x = rgb.astype(np.float32) / 255.0
    x = (x - _MEAN) / _STD
    return np.ascontiguousarray(x.transpose(2, 0, 1))


def denormalize_image(x: np.ndarray) -> np.ndarray:
    arr = np.asarray(x).transpose(1, 2, 0) * _STD + _MEAN
    return np.clip(arr * 255.0, 0, 255).astype(np.uint8)


def relabel_sequential(id_map: np.ndarray, class_idx_of_id: dict) -> tuple[np.ndarray, dict]:
    """Compact ids to 1..N after augmentation may have removed some entirely."""
    present = [i for i in np.unique(id_map) if i != 0]
    out = np.zeros_like(id_map, dtype=np.int32)
    classes: dict[int, int] = {}
    for new_id, old_id in enumerate(sorted(present), start=1):
        out[id_map == old_id] = new_id
        classes[new_id] = int(class_idx_of_id.get(int(old_id), 0))
    return out, classes


class PanopticPatchDataset(Dataset):
    def __init__(
        self,
        meta_paths: Sequence[str | Path],
        cfg: dict,
        class_spec: ClassSpec,
        *,
        train: bool = True,
        seed: int = 0,
    ):
        self.metas = [Path(p) for p in meta_paths]
        if not self.metas:
            raise ValueError("PanopticPatchDataset got an empty patch list.")
        self.cfg = cfg
        self.spec = class_spec
        self.train = train
        self.seed = seed
        self.aug_cfg = dict(cfg.get("augment", {}) or {})
        self.augment = bool(train and cfg.get("train", {}).get("augment", True))
        self.strategy = cfg.get("ignore_strategy", "dilated_gt_two_radius_plus_tissue")
        self.ann = dict(cfg.get("annotation", {}) or {})
        self.compact_indices = list(class_spec.compact_indices)
        self.touch_weight = float(
            cfg.get("loss", {}).get("boundary", {}).get("touching_weight", 1.0)
        )
        self.index = self._build_index()
        self.paste_bank = self._build_paste_bank() if self.augment else []

    # ------------------------------------------------------------------ index
    def _build_index(self) -> list[int]:
        """Sampling plan: cap background-only patches, then oversample rare classes.

        Extraction keeps glassy background patches generously because they are what
        supervises trusted background, but the resulting 3:1 background:foreground ratio
        spends most of each epoch on tiles with almost nothing to learn from.
        ``train.bg_patch_ratio`` caps them relative to the foreground patches; the subset is
        chosen deterministically so epochs are reproducible.
        """
        if not self.train:
            return list(range(len(self.metas)))

        counts = [_read_counts(mp) for mp in self.metas]
        fg = [i for i, c in enumerate(counts) if sum(c.values()) > 0]
        bg = [i for i, c in enumerate(counts) if sum(c.values()) == 0]

        ratio = self.cfg.get("train", {}).get("bg_patch_ratio")
        if ratio is not None and bg:
            keep = int(round(float(ratio) * len(fg)))
            if keep < len(bg):
                rng = np.random.default_rng(self.seed)
                bg = sorted(rng.permutation(bg)[:keep].tolist())
        base = sorted(fg + bg)

        over = dict(self.cfg.get("patch", {}).get("oversample", {}) or {})
        if not over:
            return base
        idx: list[int] = []
        for i in base:
            factor = max([1] + [int(over[c]) for c in over if counts[i].get(c, 0) > 0])
            idx.extend([i] * factor)
        return idx

    def _build_paste_bank(self) -> list[tuple[int, str]]:
        """(patch index, class name) pairs holding at least one rare-class instance."""
        cp = dict(self.aug_cfg.get("copy_paste", {}) or {})
        wanted = set(cp.get("classes", []) or [])
        if not wanted or cp.get("p", 0) <= 0:
            return []
        bank: list[tuple[int, str]] = []
        for i, mp in enumerate(self.metas):
            counts = _read_counts(mp)
            for c in wanted:
                if counts.get(c, 0) > 0:
                    bank.append((i, c))
        return bank

    def __len__(self) -> int:
        return len(self.index)

    # ------------------------------------------------------------------- item
    def __getitem__(self, i: int) -> dict:
        pi = self.index[i]
        rng = np.random.default_rng((self.seed * 1_000_003 + i * 7919) % (2 ** 32))
        patch = load_patch(self.metas[pi])

        s = Sample(
            image=patch["image"],
            id_map=patch["id_map"].astype(np.int32),
            tissue=patch["tissue"],
            bg_poly=patch["bg_poly"],
            dense_roi=patch["dense_roi"],
            inbounds=np.ones(patch["id_map"].shape, dtype=bool),
            class_idx_of_id=patch["class_idx_of_id"],
        )

        if self.augment:
            s = self._maybe_paste(s, rng)
            s = augment_sample(s, rng, self.aug_cfg)

        id_map, class_idx_of_id = relabel_sequential(s.id_map, s.class_idx_of_id)

        valid = build_valid_mask(
            id_map,
            strategy=self.strategy,
            tissue_mask=s.tissue,
            background_mask=s.bg_poly,
            dense_roi_mask=s.dense_roi,
            dilate_px_inner=int(self.ann.get("dilate_px_inner", 5)),
            dilate_px=int(self.ann.get("dilate_px", 30)),
        )
        # Content warped in from outside the extracted patch is not labelled: ignore it.
        valid &= s.inbounds

        t = targets_from_config(id_map, class_idx_of_id, self.compact_indices, self.cfg)

        weight = np.ones(id_map.shape, dtype=np.float32)
        if self.touch_weight != 1.0:
            weight[t.touching] = self.touch_weight

        return {
            "image": torch.from_numpy(normalize_image(s.image)),
            "semantic": torch.from_numpy(apply_ignore_to_semantic(t.semantic, valid, IGNORE_INDEX)),
            "boundary": torch.from_numpy(apply_ignore_to_semantic(t.boundary, valid, IGNORE_INDEX)),
            "centers": torch.from_numpy(t.centers),
            "valid": torch.from_numpy(valid.astype(np.float32)),
            "weight": torch.from_numpy(weight),
            "id_map": torch.from_numpy(id_map.astype(np.int32)),
            "patch_id": patch["meta"]["patch_id"],
        }

    def _maybe_paste(self, s: Sample, rng: np.random.Generator) -> Sample:
        cp = dict(self.aug_cfg.get("copy_paste", {}) or {})
        if not self.paste_bank or rng.random() >= float(cp.get("p", 0.0)):
            return s
        n = int(rng.integers(1, int(cp.get("max_objects", 1)) + 1))
        for _ in range(n):
            di, cname = self.paste_bank[int(rng.integers(0, len(self.paste_bank)))]
            donor = load_patch(self.metas[di])
            cls_idx = self.spec.index(cname)
            candidates = [k for k, v in donor["class_idx_of_id"].items() if v == cls_idx]
            if not candidates:
                continue
            inst = int(candidates[int(rng.integers(0, len(candidates)))])
            s = copy_paste(s, donor["image"], donor["id_map"], inst, cls_idx, rng,
                           blend_px=int(cp.get("blend_px", 3)))
        return s


def _read_counts(meta_path: Path) -> dict:
    with open(meta_path) as f:
        return json.load(f).get("class_counts", {}) or {}


def split_metas(metas: Sequence[Path], manifest_rows: Sequence[dict]) -> dict[str, list[Path]]:
    """Group patch metas into train/val by their slide's manifest split."""
    split_of = {r["slide_id"]: r.get("split", "train") for r in manifest_rows}
    out: dict[str, list[Path]] = {}
    for m in metas:
        slide = m.parent.name
        out.setdefault(split_of.get(slide, "train"), []).append(m)
    return out
