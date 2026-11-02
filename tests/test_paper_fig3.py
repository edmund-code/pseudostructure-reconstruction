"""Synthetic checks for the Fig. 2-3 paper analyses (analysis/paper/)."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'analysis/paper'))


def test_capsule_depth_excludes_the_hilum_notch():
    from fig3_depth import capsule_and_depth

    yy, xx = np.mgrid[:201, :201]
    mask = (yy - 100) ** 2 + (xx - 100) ** 2 <= 90 ** 2
    mask &= ~((np.abs(yy - 100) < 8) & (xx > 140))            # a deep notch, like a hilum
    capsule, depth = capsule_and_depth(mask, px_um=1.0, capsule_um=5.0)
    assert not capsule[100, 150:185].any()                    # notch walls are not capsule
    assert depth[100, 100] > 80     # nearest capsule is the rim or the notch lips, not the notch floor (~41 px away)
    assert depth[100, 135] > 40                               # beside the notch: still far from the capsule


def test_tissue_mask_keeps_the_largest_saturated_object():
    from fig3_depth import tissue_mask

    rgb = np.full((120, 120, 3), 245, np.uint8)
    rgb[20:100, 20:100] = (230, 120, 200)
    rgb[50:60, 50:60] = (245, 245, 245)                        # an unstained lumen: filled
    rgb[5:8, 110:113] = (230, 120, 200)                        # a speck: dropped
    m = tissue_mask(rgb)
    assert m[55, 55] and m[30, 30] and not m[6, 111]


def test_densest_window_finds_the_cluster():
    from fig2_gallery import densest_window

    rng = np.random.default_rng(0)
    xy = np.vstack([rng.uniform(0, 1000, (200, 2)), rng.uniform(600, 650, (60, 2))])
    x, y = densest_window(xy, 100, 10)
    assert 550 <= x <= 610 and 550 <= y <= 610
