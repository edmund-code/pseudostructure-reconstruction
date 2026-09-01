"""Prediction heads: semantic, instance boundary/interior, compact-class centers."""
from __future__ import annotations

import torch
import torch.nn as nn

from .decoder import _norm


class ConvHead(nn.Module):
    """3x3 -> norm -> ReLU -> 1x1 logits."""

    def __init__(self, cin: int, cout: int, hidden: int | None = None, norm: str = "gn"):
        super().__init__()
        h = hidden or cin
        self.head = nn.Sequential(
            nn.Conv2d(cin, h, 3, padding=1, bias=False),
            _norm(norm, h),
            nn.ReLU(inplace=True),
            nn.Conv2d(h, cout, 1),
        )

    def forward(self, x):
        return self.head(x)


class SemanticHead(ConvHead):
    """Per-pixel class logits (index 0 = background)."""


class BoundaryHead(ConvHead):
    """Instance-structure logits.

    3 classes: background / interior / boundary.
    4 classes: background / interior / boundary-to-background / boundary-to-INSTANCE.

    The 4-class form gives the wall between two touching objects its own output. That wall
    is what must fire to keep two tubules apart, and it is a small minority of an already
    small boundary class — folded into one channel it is dominated by the far more common
    (and far easier) object-to-background contour.
    """

    def __init__(self, cin: int, n_classes: int = 3, hidden: int | None = None,
                 norm: str = "gn"):
        if n_classes not in (3, 4):
            raise ValueError(f"BoundaryHead n_classes must be 3 or 4, got {n_classes}")
        super().__init__(cin, n_classes, hidden, norm)
        self.n_classes = n_classes


class CenterHead(nn.Module):
    """Gaussian center heatmaps for compact classes, one channel per class.

    The final bias is initialised strongly negative so the head starts predicting ~0
    everywhere. The targets are almost entirely zero, and without this the first
    iterations are dominated by shrinking a uniform 0.5 response.
    """

    def __init__(self, cin: int, n_compact: int, hidden: int | None = None, norm: str = "gn"):
        super().__init__()
        self.n_compact = n_compact
        self.head = ConvHead(cin, max(n_compact, 1), hidden, norm)
        with torch.no_grad():
            self.head.head[-1].bias.fill_(-4.0)

    def forward(self, x):
        return self.head(x)


class HVHead(ConvHead):
    """EXPERIMENTAL centroid-offset branch (horizontal/vertical).

    Gated behind ``model.legacy_hv`` and deliberately NOT part of the default decode path
    (docs/adr/0001). It exists only so the comparison can be run without a code fork.
    """

    def __init__(self, cin: int, hidden: int | None = None, norm: str = "gn"):
        super().__init__(cin, 2, hidden, norm)
