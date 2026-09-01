import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from kidney_panoptic.data.classes import ClassSpec  # noqa: E402


# Mirrors configs/default.yaml (v4). Index meanings used throughout the tests:
#   0 background | 1 tubule_proximal | 2 tubule_distal | 3 tubule_collecting | 4 glomerulus
# Vessel / Blood Cells / Nerve are IGNORED, not a class — see the default.yaml rationale.
CFG = {
    "classes": ["background", "tubule_proximal", "tubule_distal",
                "tubule_collecting", "glomerulus"],
    "compact_classes": ["glomerulus"],
    "elongated_classes": ["tubule_proximal", "tubule_distal", "tubule_collecting"],
    "class_map": {"Tubules": "tubule_proximal", "Tubules2": "tubule_distal",
                  "Tubules3": "tubule_collecting", "Glomeruli": "glomerulus"},
    "background_names": ["Background"],
    "dense_roi_names": ["ROI", "DenseROI", "Region*"],
    "ignore_names": ["ignore", "Vessel", "Blood Cells", "Nerve"],
}

PROXIMAL, DISTAL, COLLECTING, GLOMERULUS = 1, 2, 3, 4
# Back-compat aliases so the geometry-focused tests keep reading naturally.
TUBULE_1, TUBULE_2, TUBULE_3 = PROXIMAL, DISTAL, COLLECTING


@pytest.fixture
def spec():
    return ClassSpec.from_config(CFG)


@pytest.fixture
def cfg():
    return dict(CFG)


@pytest.fixture
def touching_pair():
    """Two 40x40 squares sharing an edge at x=40 — zero gap, the tubule failure case."""
    lbl = np.zeros((64, 96), np.int32)
    lbl[12:52, 8:48] = 1
    lbl[12:52, 48:88] = 2
    return lbl
