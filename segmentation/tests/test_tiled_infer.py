"""Blending correctness, tested with a fake model — no weights, no GPU.

The fake model reports back what it was shown, so a correct blend must reconstruct the
source exactly. That pins the three things that actually go wrong in tiled inference:
weights that do not normalise to 1, tiles placed off the global grid, and a halo too small
for a region's edge pixels to receive every contribution they would get from one big pass.
"""
import numpy as np
import pytest
import torch

from kidney_panoptic.data.classes import ClassSpec
from kidney_panoptic.infer.tiled import TiledPredictor, array_fetcher, gaussian_window
from conftest import CFG


class EchoModel:
    """Returns the tile's own red channel as semantic class-1 probability."""

    def __init__(self, n_classes, n_compact):
        self.n_classes, self.n_compact = n_classes, n_compact

    def to(self, device):
        return self

    def eval(self):
        return self

    def predict_maps(self, x, amp=False):
        B, _, H, W = x.shape
        # undo the ImageNet normalisation of the red channel -> [0, 1]
        red = (x[:, 0] * 0.229 + 0.485).clamp(0, 1)
        sem = torch.zeros(B, self.n_classes, H, W)
        sem[:, 1] = red
        sem[:, 0] = 1 - red
        bnd = torch.zeros(B, 3, H, W)
        bnd[:, 0] = 1.0
        return {"semantic_prob": sem, "boundary_prob": bnd,
                "center_heatmaps": torch.zeros(B, self.n_compact, H, W)}


def make_predictor(cfg_extra=None):
    spec = ClassSpec.from_config(CFG)
    cfg = {**CFG, "infer": {"tile": 64, "overlap": 16, "batch_size": 4,
                            "blend_sigma_scale": 0.25, "tissue_gate": False},
           "decode": {}}
    if cfg_extra:
        cfg["infer"].update(cfg_extra)
    model = EchoModel(spec.n_classes, len(spec.compact_indices))
    return TiledPredictor(model, cfg, spec, device="cpu", amp=False), spec


def gradient_image(h, w):
    """A pattern with no symmetry, so any misalignment shows up immediately."""
    yy, xx = np.mgrid[0:h, 0:w]
    r = ((xx * 7 + yy * 13) % 251).astype(np.uint8)
    return np.stack([r, r, r], axis=-1)


def test_blend_reconstructs_the_source_exactly():
    pred, _ = make_predictor()
    img = gradient_image(200, 240)
    maps = pred.predict_region(array_fetcher(img), 0, 0, 240, 200, 240, 200)
    got = maps["semantic_prob"][1]
    want = img[:, :, 0].astype(np.float32) / 255.0
    assert np.abs(got - want).max() < 2e-2, "blended maps must reproduce the tile outputs"


def test_weights_sum_to_a_positive_value_everywhere():
    pred, _ = make_predictor()
    img = gradient_image(150, 150)
    maps = pred.predict_region(array_fetcher(img), 0, 0, 150, 150, 150, 150)
    assert (maps["weight"] > 0).all(), "an unweighted pixel would divide by ~zero"


def test_probabilities_stay_normalized_after_blending():
    pred, _ = make_predictor()
    img = gradient_image(128, 128)
    maps = pred.predict_region(array_fetcher(img), 0, 0, 128, 128, 128, 128)
    total = maps["semantic_prob"].sum(axis=0)
    assert np.allclose(total, 1.0, atol=1e-3), "a convex blend of distributions is a distribution"


def test_a_sub_region_matches_the_same_area_of_the_whole():
    """The halo is what makes a windowed decode equal to one slide-wide pass."""
    pred, _ = make_predictor()
    img = gradient_image(256, 256)
    whole = pred.predict_region(array_fetcher(img), 0, 0, 256, 256, 256, 256)
    part = pred.predict_region(array_fetcher(img), 64, 64, 128, 128, 256, 256)
    a = whole["semantic_prob"][1][64:192, 64:192]
    b = part["semantic_prob"][1]
    assert np.abs(a - b).max() < 1e-4, "windowed inference must match the monolithic pass"


def test_region_at_the_slide_edge_is_handled():
    pred, _ = make_predictor()
    img = gradient_image(100, 100)
    maps = pred.predict_region(array_fetcher(img), 40, 40, 60, 60, 100, 100)
    assert maps["semantic_prob"].shape[-2:] == (60, 60)
    assert (maps["weight"] > 0).all()


def test_gaussian_window_downweights_tile_edges():
    """The point of weighting: a tile's own edge is where its prediction is least reliable."""
    w = gaussian_window(64)
    centre = w[32, 32]
    edge = max(w[0, 32], w[32, 0], w[63, 32], w[32, 63])
    # sigma = 0.25 * tile, so a mid-edge pixel sits 2 sigma out: exp(-2) ~ 0.135
    assert edge < 0.2 * centre
    assert w[0, 0] < edge, "corners, two sigma out on both axes, are weighted least of all"
    assert w.min() > 0, "still strictly positive so every pixel gets normalised"


def test_blend_reconstructs_a_constant_map_exactly():
    """A convex blend of identical predictions must return that prediction untouched."""
    pred, _ = make_predictor()
    img = np.full((150, 170, 3), 128, np.uint8)
    maps = pred.predict_region(array_fetcher(img), 0, 0, 170, 150, 170, 150)
    assert np.allclose(maps["semantic_prob"][1], 128 / 255.0, atol=1e-3)


def test_predictor_channel_layout_follows_the_models_boundary_width():
    """Regression: the accumulator width was hardcoded to 3 boundary channels.

    Switching model.boundary_classes to 4 then produced a shape mismatch deep inside the
    blend accumulator, which only surfaced at region-eval time on a real slide. The
    predictor must take the width from the MODEL, never assume it.
    """
    spec = ClassSpec.from_config(CFG)
    C, K = spec.n_classes, len(spec.compact_indices)

    for n_boundary in (3, 4):
        class Model(EchoModel):
            boundary_classes = n_boundary

            def predict_maps(self, x, amp=False):
                B, _, H, W = x.shape
                red = (x[:, 0] * 0.229 + 0.485).clamp(0, 1)
                sem = torch.zeros(B, C, H, W)
                sem[:, 1] = red
                sem[:, 0] = 1 - red
                bnd = torch.zeros(B, n_boundary, H, W)
                bnd[:, 0] = 1.0
                return {"semantic_prob": sem, "boundary_prob": bnd,
                        "center_heatmaps": torch.zeros(B, K, H, W)}

        cfg = {**CFG, "infer": {"tile": 64, "overlap": 16, "batch_size": 4,
                                "tissue_gate": False}, "decode": {}}
        pred = TiledPredictor(Model(C, K), cfg, spec, device="cpu", amp=False)
        assert pred.n_boundary == n_boundary
        maps = pred.predict_region(array_fetcher(gradient_image(100, 100)),
                                   0, 0, 100, 100, 100, 100)
        assert maps["boundary_prob"].shape[0] == n_boundary
        assert maps["semantic_prob"].shape[0] == C
        assert maps["center_heatmaps"].shape[0] == K
