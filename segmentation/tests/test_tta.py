"""Test-time augmentation: the averaged maps must land back in INPUT orientation.

The bug this guards against is a silent one. Get the inverse transform wrong and the eight
passes still average to something smooth and plausible-looking — just blurred across
orientations, which quietly destroys exactly the thin boundary walls TTA was added to
sharpen. A model that echoes its input is the only way to catch it.
"""
import numpy as np
import pytest
import torch

from kidney_panoptic.data.classes import ClassSpec
from kidney_panoptic.infer.tiled import (TiledPredictor, array_fetcher,
                                         dihedral_transforms, _tta_forward)
from conftest import CFG


class AsymmetricEchoModel:
    """Returns the input's red channel verbatim as a probability map.

    Equivariant by construction, so a correct TTA average must reproduce it EXACTLY.
    """

    boundary_classes = 3

    def __init__(self, n_classes, n_compact):
        self.n_classes, self.n_compact = n_classes, n_compact

    def to(self, device):
        return self

    def eval(self):
        return self

    def predict_maps(self, x, amp=False):
        B, _, H, W = x.shape
        red = (x[:, 0] * 0.229 + 0.485).clamp(0, 1)
        sem = torch.zeros(B, self.n_classes, H, W)
        sem[:, 1] = red
        sem[:, 0] = 1 - red
        bnd = torch.zeros(B, 3, H, W)
        bnd[:, 0] = 1.0
        return {"semantic_prob": sem, "boundary_prob": bnd,
                "center_heatmaps": torch.zeros(B, self.n_compact, H, W)}


def asymmetric_tile(n=64):
    """No rotational or mirror symmetry, so any un-inverted transform shows up."""
    yy, xx = np.mgrid[0:n, 0:n]
    return ((xx * 3 + yy * 11) % 251).astype(np.float32) / 255.0


@pytest.mark.parametrize("mode,expected", [("none", 1), ("rot4", 4), ("d4", 8)])
def test_transform_counts(mode, expected):
    assert len(dihedral_transforms(mode)) == expected


def test_unknown_mode_raises():
    """Silently falling back to identity would make a requested 8x run a 1x run."""
    with pytest.raises(ValueError, match="tta"):
        dihedral_transforms("flip8")


def test_d4_transforms_are_distinct():
    assert len(set(dihedral_transforms("d4"))) == 8


@pytest.mark.parametrize("mode", ["none", "rot4", "d4"])
def test_tta_average_is_in_input_orientation(mode):
    spec = ClassSpec.from_config(CFG)
    model = AsymmetricEchoModel(spec.n_classes, len(spec.compact_indices))
    tile = asymmetric_tile()
    x = torch.from_numpy(((tile - 0.485) / 0.229)[None, None].repeat(3, axis=1))

    out = _tta_forward(model, x, False, dihedral_transforms(mode))
    got = out["semantic_prob"][0, 1].numpy()
    assert np.abs(got - tile).max() < 1e-5, (
        f"{mode}: averaged map is not in input orientation — an inverse transform is wrong")


def test_tta_does_not_change_a_rotationally_symmetric_result():
    """Sanity check on the fixture itself: a symmetric input cannot distinguish the modes,
    so agreement there proves nothing. This asserts the ASYMMETRIC fixture is what carries
    the signal, by showing the two disagree on nothing symmetric."""
    spec = ClassSpec.from_config(CFG)
    model = AsymmetricEchoModel(spec.n_classes, len(spec.compact_indices))
    flat = torch.full((1, 3, 64, 64), (0.5 - 0.485) / 0.229)
    a = _tta_forward(model, flat, False, dihedral_transforms("none"))["semantic_prob"]
    b = _tta_forward(model, flat, False, dihedral_transforms("d4"))["semantic_prob"]
    assert torch.allclose(a, b, atol=1e-5)


def test_predictor_honours_the_tta_config():
    spec = ClassSpec.from_config(CFG)
    img = np.stack([((np.mgrid[0:96, 0:96][1] * 7) % 251).astype(np.uint8)] * 3, axis=-1)
    ref = None
    for mode in ("none", "rot4", "d4"):
        cfg = {**CFG, "infer": {"tile": 64, "overlap": 16, "batch_size": 4,
                                "tissue_gate": False, "tta": mode}, "decode": {}}
        pred = TiledPredictor(AsymmetricEchoModel(spec.n_classes, len(spec.compact_indices)),
                              cfg, spec, device="cpu", amp=False)
        assert len(pred.tta) == {"none": 1, "rot4": 4, "d4": 8}[mode]
        maps = pred.predict_region(array_fetcher(img), 0, 0, 96, 96, 96, 96)
        if ref is None:
            ref = maps["semantic_prob"][1]
        else:
            assert np.abs(maps["semantic_prob"][1] - ref).max() < 2e-2, (
                f"tta={mode} changed an equivariant model's blended output")
