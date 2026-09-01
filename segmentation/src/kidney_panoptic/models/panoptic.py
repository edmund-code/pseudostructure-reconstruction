"""The assembled panoptic model.

    model(rgb) -> {
        semantic_logits : (B, C, H, W)
        boundary_logits : (B, 3, H, W)
        center_logits   : (B, K, H, W)      # K = number of compact classes
        center_heatmaps : (B, K, H, W)      # sigmoid of the above, for convenience
        hv              : (B, 2, H, W)      # only when model.legacy_hv is true
    }

All maps are returned at the INPUT patch resolution: when ``decoder.out_stride == 2`` the
logits are bilinearly upsampled by 2 as the last step, so downstream code never has to
think about stride. The encoder is frozen; only the stem, decoder and heads train.
"""
from __future__ import annotations

from typing import Optional

import torch
import torch.nn as nn
import torch.nn.functional as F

from ..data.classes import ClassSpec
from .decoder import FusionDecoder
from .encoder import FrozenEncoder
from .heads import BoundaryHead, CenterHead, HVHead, SemanticHead


class KidneyPanoptic(nn.Module):
    def __init__(self, cfg: dict, class_spec: ClassSpec):
        super().__init__()
        self.cfg = cfg
        self.spec = class_spec
        self.n_classes = class_spec.n_classes
        self.compact_indices = list(class_spec.compact_indices)

        self.encoder = FrozenEncoder(dict(cfg.get("encoder", {}) or {}))
        dec_cfg = dict(cfg.get("decoder", {}) or {})
        self.decoder = FusionDecoder(self.encoder.out_dims, dec_cfg)
        d = self.decoder.out_channels
        norm = dec_cfg.get("norm", "gn")

        self.boundary_classes = int(cfg.get("model", {}).get("boundary_classes", 3))
        self.semantic_head = SemanticHead(d, self.n_classes, norm=norm)
        self.boundary_head = BoundaryHead(d, self.boundary_classes, norm=norm)
        self.center_head = CenterHead(d, len(self.compact_indices), norm=norm)

        self.legacy_hv = bool(cfg.get("model", {}).get("legacy_hv", False))
        self.hv_head = HVHead(d, norm=norm) if self.legacy_hv else None

    # ------------------------------------------------------------------ params
    def trainable_parameters(self):
        return [p for p in self.parameters() if p.requires_grad]

    def n_trainable(self) -> int:
        return sum(p.numel() for p in self.trainable_parameters())

    # ----------------------------------------------------------------- forward
    def forward(self, rgb: torch.Tensor) -> dict[str, torch.Tensor]:
        tokens = self.encoder(rgb)
        feat = self.decoder(rgb, tokens)

        out = {
            "semantic_logits": self.semantic_head(feat),
            "boundary_logits": self.boundary_head(feat),
            "center_logits": self.center_head(feat),
        }
        if self.hv_head is not None:
            out["hv"] = torch.tanh(self.hv_head(feat))

        size = rgb.shape[-2:]
        for k, v in list(out.items()):
            if v.shape[-2:] != size:
                out[k] = F.interpolate(v, size=size, mode="bilinear", align_corners=False)
        out["center_heatmaps"] = torch.sigmoid(out["center_logits"])
        return out

    # -------------------------------------------------------------- inference
    @torch.no_grad()
    def predict_maps(self, rgb: torch.Tensor, amp: bool = True) -> dict[str, torch.Tensor]:
        """Probability maps for decoding: semantic_prob, boundary_prob, center_heatmaps."""
        was_training = self.training
        self.eval()
        dev_type = rgb.device.type
        with torch.autocast(device_type=dev_type, dtype=torch.float16,
                            enabled=bool(amp) and dev_type == "cuda"):
            out = self(rgb)
        res = {
            "semantic_prob": torch.softmax(out["semantic_logits"].float(), dim=1),
            "boundary_prob": torch.softmax(out["boundary_logits"].float(), dim=1),
            "center_heatmaps": torch.sigmoid(out["center_logits"].float()),
        }
        if was_training:
            self.train()
        return res


def build_model(cfg: dict, class_spec: Optional[ClassSpec] = None) -> KidneyPanoptic:
    spec = class_spec or ClassSpec.from_config(cfg)
    return KidneyPanoptic(cfg, spec)
