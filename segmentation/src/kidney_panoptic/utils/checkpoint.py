"""Checkpoint save/load.

The config and the class vocabulary travel WITH the weights. A checkpoint whose class order
came from a config that has since been edited would silently relabel every prediction, so
the vocabulary is read back from the checkpoint at inference time rather than from whatever
config happens to be on disk.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import torch

from ..data.classes import ClassSpec

FORMAT_VERSION = 2


def save_checkpoint(
    path: str | Path,
    model,
    cfg: dict,
    spec: ClassSpec,
    *,
    epoch: int = 0,
    metrics: Optional[dict] = None,
    optimizer=None,
    extra: Optional[dict] = None,
) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # The frozen encoder is 4.5 GB and fully reproducible from its name — don't store it.
    state = {k: v for k, v in model.state_dict().items() if not k.startswith("encoder.")}
    payload = {
        "format_version": FORMAT_VERSION,
        "model_state": state,
        "config": dict(cfg),
        "classes": list(spec.names),
        "compact_classes": list(spec.compact),
        "elongated_classes": list(spec.elongated),
        "class_map": dict(spec.mapping),
        "epoch": int(epoch),
        "metrics": metrics or {},
        **(extra or {}),
    }
    if optimizer is not None:
        payload["optimizer_state"] = optimizer.state_dict()
    torch.save(payload, path)
    return path


def load_checkpoint(path: str | Path, map_location: str = "cpu") -> dict:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"Checkpoint not found: {path}\n"
            f"Train one with scripts/train.py, or pass --checkpoint explicitly."
        )
    ckpt = torch.load(path, map_location=map_location, weights_only=False)
    for key in ("model_state", "config", "classes"):
        if key not in ckpt:
            raise ValueError(
                f"{path} is not a kidney_panoptic checkpoint (missing '{key}'). "
                f"Keys present: {sorted(ckpt)[:10]}"
            )
    return ckpt


def build_from_checkpoint(path: str | Path, overrides: Optional[dict] = None,
                          device: str = "cpu"):
    """Rebuild (model, cfg, spec) from a checkpoint, restoring its class vocabulary."""
    from ..models.panoptic import build_model

    ckpt = load_checkpoint(path, map_location="cpu")
    cfg = dict(ckpt["config"])
    if overrides:
        cfg.update(overrides)
    cfg["classes"] = list(ckpt["classes"])
    cfg["compact_classes"] = list(ckpt.get("compact_classes", []))
    cfg["elongated_classes"] = list(ckpt.get("elongated_classes", []))
    cfg["class_map"] = dict(ckpt.get("class_map", {}))

    spec = ClassSpec.from_config(cfg)
    model = build_model(cfg, spec)
    missing, unexpected = model.load_state_dict(ckpt["model_state"], strict=False)
    non_encoder_missing = [k for k in missing if not k.startswith("encoder.")]
    if non_encoder_missing:
        raise RuntimeError(
            f"Checkpoint {path} is missing {len(non_encoder_missing)} trainable "
            f"parameters (first few: {non_encoder_missing[:5]}). The config's decoder/head "
            f"shape does not match the one that was trained."
        )
    if unexpected:
        raise RuntimeError(f"Checkpoint {path} has unexpected keys: {unexpected[:5]}")
    return model.to(device), cfg, spec
