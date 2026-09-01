"""Frozen pathology foundation encoder producing a multi-layer token pyramid.

``openmidnight_vitg14_reg`` (default) is the OpenMidnight teacher — a DINOv2 ViT-g/14
**with registers** trained on pathology — loaded from the Hugging Face hub into the
DINOv2 architecture. ``dinov2_vitg14_reg`` is the natural-image DINOv2 for comparison, and
``random`` is a tiny CPU stub with the same interface so tests and the decode path can run
without downloading 4.5 GB of weights.

Multi-layer extraction
----------------------
A single final-layer token map is a poor pyramid: the last block of a ViT is heavily
specialised toward the pretraining objective, while mid-depth blocks retain more
localisable texture. ``encoder.layers`` pulls several blocks out and the decoder learns
how to weight them.

Everything here is ``requires_grad=False`` and runs under ``no_grad``. ``no_grad`` rather
than ``inference_mode`` matters: the output has to be usable as a constant leaf in the
trainable decoder's autograd graph.
"""
from __future__ import annotations

from pathlib import Path

import torch
import torch.nn as nn
import torch.nn.functional as F

_KNOWN = ("openmidnight_vitg14_reg", "dinov2_vitg14_reg", "random")


class FrozenEncoder(nn.Module):
    def __init__(self, cfg: dict):
        super().__init__()
        self.name = cfg.get("name", "openmidnight_vitg14_reg")
        if self.name not in _KNOWN:
            raise ValueError(f"Unknown encoder '{self.name}'. Known: {list(_KNOWN)}")
        self.input_size = int(cfg.get("input_size", 518))
        self.embed_dim = int(cfg.get("embed_dim", 1536))
        self.layers = list(cfg.get("layers", [9, 19, 29, 39]))
        self.patch_size = 14

        if self.input_size % self.patch_size:
            raise ValueError(
                f"encoder.input_size ({self.input_size}) must be a multiple of "
                f"{self.patch_size}; e.g. 518 = 14 x 37."
            )
        self.token_grid = self.input_size // self.patch_size

        if self.name == "random":
            self.backbone = _StubEncoder(self.embed_dim, len(self.layers))
        else:
            self.backbone = _load_vitg14_reg(self.name)

        for p in self.parameters():
            p.requires_grad = False
        self.eval()

    @property
    def out_dims(self) -> list[int]:
        return [self.embed_dim] * len(self.layers)

    def train(self, mode: bool = True):  # stay frozen even inside model.train()
        return super().train(False)

    @torch.no_grad()
    def forward(self, images: torch.Tensor) -> list[torch.Tensor]:
        """(B, 3, H, W) ImageNet-normalized -> list of (B, embed_dim, g, g) maps."""
        x = images
        if x.shape[-2:] != (self.input_size, self.input_size):
            x = F.interpolate(x, size=(self.input_size, self.input_size),
                              mode="bilinear", align_corners=False)
        if self.name == "random":
            return self.backbone(x, self.token_grid)
        feats = self.backbone.get_intermediate_layers(
            x, n=self.layers, reshape=True, norm=True
        )
        return [f.float() for f in feats]


class _StubEncoder(nn.Module):
    """Deterministic, untrained stand-in with the real output shapes (tests / CI)."""

    def __init__(self, embed_dim: int, n_layers: int):
        super().__init__()
        self.proj = nn.Conv2d(3, embed_dim, kernel_size=1, bias=False)
        self.n_layers = n_layers

    def forward(self, x: torch.Tensor, grid: int) -> list[torch.Tensor]:
        y = F.adaptive_avg_pool2d(x, (grid, grid))
        f = self.proj(y)
        return [f * (i + 1) / self.n_layers for i in range(self.n_layers)]


def _load_vitg14_reg(name: str) -> nn.Module:
    """Build DINOv2 ViT-g/14-with-registers and load the requested weights."""
    if name == "dinov2_vitg14_reg":
        return torch.hub.load("facebookresearch/dinov2", "dinov2_vitg14_reg",
                              pretrained=True).eval()

    from huggingface_hub import hf_hub_download

    try:
        ckpt_path = hf_hub_download(repo_id="SophontAI/OpenMidnight",
                                    filename="teacher_checkpoint_load.pt")
    except Exception as e:  # pragma: no cover - network / auth
        raise RuntimeError(
            "Could not fetch the OpenMidnight encoder from the Hugging Face hub "
            "(repo SophontAI/OpenMidnight, file teacher_checkpoint_load.pt). "
            "Set HF_HOME to a cache that already has it, or use "
            "--set encoder.name=dinov2_vitg14_reg. Original error: " + repr(e)
        ) from e

    model = torch.hub.load("facebookresearch/dinov2", "dinov2_vitg14_reg", pretrained=False)
    state = torch.load(ckpt_path, map_location="cpu", weights_only=False)
    # The released checkpoint carries its own positional embedding shape.
    model.pos_embed = nn.Parameter(state["pos_embed"])
    missing, unexpected = model.load_state_dict(state, strict=False)
    if missing:
        raise RuntimeError(
            f"OpenMidnight checkpoint is missing {len(missing)} expected parameters "
            f"(first few: {missing[:5]}). Refusing to train on a partially loaded encoder."
        )
    return model.eval()
