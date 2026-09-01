"""Trainable RGB stem + U-Net/FPN fusion decoder.

Why a stem at all
-----------------
The encoder's finest spatial unit is a 14 px token; at 0.44 µm/px that is ~6 µm, which is
wider than the wall between two touching tubules. Any decode that reads only from token
features has already thrown away the evidence needed to place the separating boundary. The
stem runs on the native patch pixels and carries stride-1/2/4 detail into the decoder via
skips, so the boundary head can localise a contour to a pixel or two while still being told
*what* it is looking at by the frozen encoder.

Shapes (patch 512, encoder input 518)
-------------------------------------
stem   : s1 512² x32, s2 256² x64, s4 128² x128, s8 64² x192, s16 32² x256
tokens : L x (1536, 37, 37) -> concat -> 1x1 -> (fuse_dim, 37, 37) -> resize -> 32²
fuse   : stride-16 tokens + stride-16 stem -> then up 8 -> 4 -> 2 -> 1
output : (decoder_channels[-1], 512, 512) when out_stride == 1
"""
from __future__ import annotations

from typing import Sequence

import torch
import torch.nn as nn
import torch.nn.functional as F


def _norm(kind: str, ch: int) -> nn.Module:
    if kind == "gn":
        return nn.GroupNorm(num_groups=min(32, max(1, ch // 8)), num_channels=ch)
    if kind == "bn":
        return nn.BatchNorm2d(ch)
    return nn.Identity()


class ConvBlock(nn.Module):
    def __init__(self, cin: int, cout: int, norm: str = "gn", stride: int = 1):
        super().__init__()
        self.block = nn.Sequential(
            nn.Conv2d(cin, cout, 3, stride=stride, padding=1, bias=False),
            _norm(norm, cout),
            nn.ReLU(inplace=True),
            nn.Conv2d(cout, cout, 3, padding=1, bias=False),
            _norm(norm, cout),
            nn.ReLU(inplace=True),
        )

    def forward(self, x):
        return self.block(x)


class RGBStem(nn.Module):
    """Pyramid of features computed from the raw patch at strides 1, 2, 4, 8, 16."""

    def __init__(self, channels: Sequence[int] = (32, 64, 128, 192, 256), norm: str = "gn"):
        super().__init__()
        if len(channels) != 5:
            raise ValueError(f"stem_channels must have 5 entries (strides 1,2,4,8,16), "
                             f"got {list(channels)}")
        c = list(channels)
        self.level0 = ConvBlock(3, c[0], norm, stride=1)
        self.level1 = ConvBlock(c[0], c[1], norm, stride=2)
        self.level2 = ConvBlock(c[1], c[2], norm, stride=2)
        self.level3 = ConvBlock(c[2], c[3], norm, stride=2)
        self.level4 = ConvBlock(c[3], c[4], norm, stride=2)
        self.channels = c

    def forward(self, x: torch.Tensor) -> list[torch.Tensor]:
        f0 = self.level0(x)
        f1 = self.level1(f0)
        f2 = self.level2(f1)
        f3 = self.level3(f2)
        f4 = self.level4(f3)
        return [f0, f1, f2, f3, f4]


class FusionDecoder(nn.Module):
    """Fuse the frozen token pyramid with the stem pyramid and upsample to out_stride."""

    def __init__(self, encoder_dims: Sequence[int], cfg: dict):
        super().__init__()
        self.out_stride = int(cfg.get("out_stride", 1))
        if self.out_stride not in (1, 2):
            raise ValueError(f"decoder.out_stride must be 1 or 2, got {self.out_stride}")
        norm = cfg.get("norm", "gn")
        fuse_dim = int(cfg.get("fuse_dim", 256))
        stem_ch = list(cfg.get("stem_channels", [32, 64, 128, 192, 256]))
        dec_ch = list(cfg.get("decoder_channels", [192, 128, 96, 64]))
        if len(dec_ch) != 4:
            raise ValueError(f"decoder.decoder_channels must have 4 entries "
                             f"(strides 8,4,2,1), got {dec_ch}")

        self.stem = RGBStem(stem_ch, norm)
        # Concat-then-project lets the decoder learn how much each encoder block matters,
        # rather than fixing a uniform sum.
        self.token_proj = nn.Sequential(
            nn.Conv2d(int(sum(encoder_dims)), fuse_dim, 1, bias=False),
            _norm(norm, fuse_dim),
            nn.ReLU(inplace=True),
        )
        self.fuse16 = ConvBlock(fuse_dim + stem_ch[4], fuse_dim, norm)
        self.up8 = ConvBlock(fuse_dim + stem_ch[3], dec_ch[0], norm)
        self.up4 = ConvBlock(dec_ch[0] + stem_ch[2], dec_ch[1], norm)
        self.up2 = ConvBlock(dec_ch[1] + stem_ch[1], dec_ch[2], norm)
        self.up1 = ConvBlock(dec_ch[2] + stem_ch[0], dec_ch[3], norm)
        self.dropout = nn.Dropout2d(float(cfg.get("dropout", 0.0)))
        self.out_channels = dec_ch[3] if self.out_stride == 1 else dec_ch[2]

    @staticmethod
    def _up_to(x: torch.Tensor, ref: torch.Tensor) -> torch.Tensor:
        return F.interpolate(x, size=ref.shape[-2:], mode="bilinear", align_corners=False)

    def forward(self, rgb: torch.Tensor, tokens: list[torch.Tensor]) -> torch.Tensor:
        s0, s1, s2, s3, s4 = self.stem(rgb)

        t = torch.cat(tokens, dim=1) if len(tokens) > 1 else tokens[0]
        t = self.token_proj(t)
        t = F.interpolate(t, size=s4.shape[-2:], mode="bilinear", align_corners=False)

        x = self.fuse16(torch.cat([t, s4], dim=1))
        x = self.up8(torch.cat([self._up_to(x, s3), s3], dim=1))
        x = self.up4(torch.cat([self._up_to(x, s2), s2], dim=1))
        x = self.up2(torch.cat([self._up_to(x, s1), s1], dim=1))
        if self.out_stride == 2:
            return self.dropout(x)
        x = self.up1(torch.cat([self._up_to(x, s0), s0], dim=1))
        return self.dropout(x)
