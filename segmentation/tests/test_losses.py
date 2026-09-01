import numpy as np
import pytest
import torch

from kidney_panoptic.losses.panoptic_loss import PanopticLoss
from kidney_panoptic.losses.seg_losses import (
    class_weights_from_counts, masked_cross_entropy, masked_dice, masked_mse,
    masked_tversky,
)

IGNORE = 255


def test_masked_ce_is_zero_for_a_perfect_confident_prediction():
    tgt = torch.zeros(1, 8, 8, dtype=torch.long)
    tgt[:, 2:6, 2:6] = 1
    logits = torch.full((1, 3, 8, 8), -20.0)
    logits.scatter_(1, tgt.unsqueeze(1), 20.0)
    assert masked_cross_entropy(logits, tgt).item() < 1e-4


def test_masked_ce_ignores_the_ignore_region_entirely():
    tgt = torch.zeros(1, 8, 8, dtype=torch.long)
    tgt[:, 4:, :] = IGNORE
    logits = torch.full((1, 3, 8, 8), -20.0)
    logits[:, 0] = 20.0                      # correct on the valid half
    a = masked_cross_entropy(logits, tgt).item()
    logits[:, :, 4:, :] = torch.randn(1, 3, 4, 8) * 50   # garbage in the ignored half
    b = masked_cross_entropy(logits, tgt).item()
    assert a == pytest.approx(b, abs=1e-6)


def test_masked_ce_normalizes_over_valid_pixels_not_all_pixels():
    """Same errors, different annotation density -> same loss."""
    def build(n_ignored):
        tgt = torch.zeros(1, 8, 8, dtype=torch.long)
        tgt[0, 0, 0] = 1                                     # one wrong pixel
        if n_ignored:
            tgt[0, -n_ignored:, :] = IGNORE
        logits = torch.full((1, 3, 8, 8), -5.0)
        logits[:, 0] = 5.0
        return tgt, logits

    t0, l0 = build(0)
    t1, l1 = build(4)
    # loss = (sum over valid) / (count valid); the single error is in both valid regions
    assert masked_cross_entropy(l1, t1).item() > masked_cross_entropy(l0, t0).item()


def test_pixel_weight_upweights_touching_pixels():
    tgt = torch.zeros(1, 8, 8, dtype=torch.long)
    logits = torch.zeros(1, 3, 8, 8)
    w = torch.ones(1, 8, 8)
    base = masked_cross_entropy(logits, tgt, pixel_weight=w).item()
    w[:, 0, 0] = 10.0
    assert masked_cross_entropy(logits, tgt, pixel_weight=w).item() == pytest.approx(base, abs=1e-5)
    # uniform logits -> per-pixel loss identical, so weighting cannot change the mean;
    # make one pixel wrong and it must.
    logits[0, 1, 0, 0] = 10.0
    assert masked_cross_entropy(logits, tgt, pixel_weight=w).item() > \
           masked_cross_entropy(logits, tgt).item()


def test_masked_dice_perfect_is_zero_and_inverted_is_one():
    tgt = torch.zeros(1, 8, 8, dtype=torch.long); tgt[:, 2:6, 2:6] = 1
    valid = torch.ones(1, 8, 8)
    good = torch.full((1, 3, 8, 8), -20.0); good.scatter_(1, tgt.unsqueeze(1), 20.0)
    assert masked_dice(good, tgt, valid).item() < 0.02
    bad = torch.full((1, 3, 8, 8), -20.0); bad[:, 0] = 20.0
    assert masked_dice(bad, tgt, valid).item() > 0.9


def test_masked_dice_skips_classes_absent_from_the_target():
    tgt = torch.zeros(1, 8, 8, dtype=torch.long); tgt[:, 2:6, 2:6] = 1
    valid = torch.ones(1, 8, 8)
    logits = torch.full((1, 4, 8, 8), -20.0); logits.scatter_(1, tgt.unsqueeze(1), 20.0)
    # Class 3 never appears; correctly predicting none of it must not be penalised.
    assert masked_dice(logits, tgt, valid).item() < 0.02


def test_tversky_beta_penalizes_false_negatives_more():
    tgt = torch.zeros(1, 8, 8, dtype=torch.long); tgt[:, 2:6, 2:6] = 1
    valid = torch.ones(1, 8, 8)
    miss = torch.full((1, 2, 8, 8), -20.0); miss[:, 0] = 20.0          # all FN
    over = torch.full((1, 2, 8, 8), -20.0); over[:, 1] = 20.0          # all FP
    fn_heavy = masked_tversky(miss, tgt, valid, alpha=0.1, beta=0.9).item()
    fp_heavy = masked_tversky(over, tgt, valid, alpha=0.1, beta=0.9).item()
    assert fn_heavy > fp_heavy


def test_masked_mse_respects_valid_and_handles_no_compact_classes():
    pred = torch.ones(1, 2, 4, 4)
    tgt = torch.zeros(1, 2, 4, 4)
    valid = torch.zeros(1, 4, 4); valid[:, :2, :] = 1
    assert masked_mse(pred, tgt, valid).item() == pytest.approx(1.0)
    assert masked_mse(torch.zeros(1, 0, 4, 4), torch.zeros(1, 0, 4, 4), valid).item() == 0.0


def test_class_weights_are_normalized_and_favour_rare_classes():
    counts = torch.tensor([1e6, 1e5, 1e4, 1e3, 1e3, 1e2])
    w = class_weights_from_counts(counts, "inverse_sqrt_freq", clip=12.0)
    assert w.mean().item() == pytest.approx(1.0, abs=1e-5)
    assert w[-1] > w[0], "the rarest class must get the largest weight"
    assert w.max() / w.min() <= 12.0 + 1e-3, "clip must bound the ratio"
    assert torch.allclose(class_weights_from_counts(counts, "none"), torch.ones(6))


def test_panoptic_loss_runs_and_reports_every_head(spec):
    cfg = {"loss": {"semantic": {"ce": 1.0, "dice": 1.0},
                    "boundary": {"ce": 1.0, "dice": 1.0, "boundary_class_weight": 3.0},
                    "center": {"mse": 25.0}, "tversky": {"enabled": False}}}
    C, K, H, W = spec.n_classes, len(spec.compact_indices), 16, 16
    out = {"semantic_logits": torch.randn(2, C, H, W, requires_grad=True),
           "boundary_logits": torch.randn(2, 3, H, W, requires_grad=True),
           "center_heatmaps": torch.rand(2, K, H, W, requires_grad=True)}
    batch = {"semantic": torch.zeros(2, H, W, dtype=torch.long),
             "boundary": torch.zeros(2, H, W, dtype=torch.long),
             "centers": torch.zeros(2, K, H, W),
             "valid": torch.ones(2, H, W),
             "weight": torch.ones(2, H, W)}
    loss, parts = PanopticLoss(cfg, C)(out, batch)
    loss.backward()
    assert torch.isfinite(loss)
    for k in ("sem_ce", "sem_dice", "bnd_ce", "bnd_region", "center_mse",
              "loss_semantic", "loss_boundary"):
        assert k in parts
    assert out["semantic_logits"].grad is not None


def test_fully_ignored_batch_gives_a_finite_loss(spec):
    """A patch with no trusted pixels must not produce NaN and stall training."""
    cfg = {"loss": {}}
    C, K, H, W = spec.n_classes, len(spec.compact_indices), 8, 8
    out = {"semantic_logits": torch.randn(1, C, H, W),
           "boundary_logits": torch.randn(1, 3, H, W),
           "center_heatmaps": torch.rand(1, K, H, W)}
    batch = {"semantic": torch.full((1, H, W), 255, dtype=torch.long),
             "boundary": torch.full((1, H, W), 255, dtype=torch.long),
             "centers": torch.zeros(1, K, H, W),
             "valid": torch.zeros(1, H, W),
             "weight": torch.ones(1, H, W)}
    loss, _ = PanopticLoss(cfg, C)(out, batch)
    assert torch.isfinite(loss)
