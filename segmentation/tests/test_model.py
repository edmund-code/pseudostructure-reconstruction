"""Model shape/wiring tests. All use encoder.name='random' so nothing is downloaded."""
import pytest
import torch

from kidney_panoptic.models.panoptic import build_model


def base_cfg(**kw):
    cfg = {
        "classes": ["background", "tubule", "glomerulus", "vessel", "blood_cell", "nerve"],
        "compact_classes": ["glomerulus", "blood_cell"],
        "elongated_classes": ["tubule", "vessel"],
        "class_map": {},
        "encoder": {"name": "random", "input_size": 518, "embed_dim": 64,
                    "layers": [0, 1], "freeze": True},
        "decoder": {"out_stride": 1, "stem_channels": [8, 16, 24, 32, 40],
                    "fuse_dim": 32, "decoder_channels": [24, 16, 12, 8], "norm": "gn"},
        "model": {"legacy_hv": False},
    }
    for k, v in kw.items():
        cfg[k] = {**cfg.get(k, {}), **v} if isinstance(v, dict) else v
    return cfg


@pytest.mark.parametrize("size", [128, 256])
def test_output_maps_are_at_input_resolution(size):
    m = build_model(base_cfg())
    out = m(torch.randn(2, 3, size, size))
    assert out["semantic_logits"].shape == (2, 6, size, size)
    assert out["boundary_logits"].shape == (2, 3, size, size)
    assert out["center_logits"].shape == (2, 2, size, size)
    assert out["center_heatmaps"].shape == (2, 2, size, size)


def test_out_stride_2_still_returns_full_resolution():
    m = build_model(base_cfg(decoder={"out_stride": 2}))
    out = m(torch.randn(1, 3, 128, 128))
    assert out["semantic_logits"].shape[-2:] == (128, 128)


def test_encoder_is_frozen_and_stays_frozen_in_train_mode():
    m = build_model(base_cfg())
    assert all(not p.requires_grad for p in m.encoder.parameters())
    m.train()
    assert not m.encoder.training, "the frozen encoder must not enter train mode"
    enc = {id(p) for p in m.encoder.parameters()}
    assert all(id(p) not in enc for p in m.trainable_parameters())
    assert m.n_trainable() > 0


def test_gradients_reach_the_stem_but_not_the_encoder():
    m = build_model(base_cfg())
    out = m(torch.randn(1, 3, 128, 128))
    out["boundary_logits"].sum().backward()
    stem_w = m.decoder.stem.level0.block[0].weight
    assert stem_w.grad is not None and stem_w.grad.abs().sum() > 0, \
        "the full-res stem must receive gradient — it is what localises boundaries"
    assert all(p.grad is None for p in m.encoder.parameters())


def test_legacy_hv_is_off_by_default_and_gated_on_by_config():
    assert "hv" not in build_model(base_cfg())(torch.randn(1, 3, 64, 64))
    m = build_model(base_cfg(model={"legacy_hv": True}))
    out = m(torch.randn(1, 3, 64, 64))
    assert out["hv"].shape == (1, 2, 64, 64)
    assert out["hv"].abs().max() <= 1.0, "hv is tanh-bounded"


def test_predict_maps_returns_normalized_probabilities():
    m = build_model(base_cfg())
    p = m.predict_maps(torch.randn(1, 3, 128, 128), amp=False)
    assert torch.allclose(p["semantic_prob"].sum(1), torch.ones(1, 128, 128), atol=1e-4)
    assert torch.allclose(p["boundary_prob"].sum(1), torch.ones(1, 128, 128), atol=1e-4)
    assert (p["center_heatmaps"] >= 0).all() and (p["center_heatmaps"] <= 1).all()


def test_predict_maps_restores_training_mode():
    m = build_model(base_cfg())
    m.train()
    m.predict_maps(torch.randn(1, 3, 64, 64), amp=False)
    assert m.training and not m.encoder.training


def test_center_head_starts_near_zero():
    """The heatmap target is ~all zeros; the head must not start at 0.5 everywhere."""
    m = build_model(base_cfg())
    out = m(torch.zeros(1, 3, 64, 64))
    assert out["center_heatmaps"].mean().item() < 0.15


def test_bad_out_stride_fails_loudly():
    with pytest.raises(ValueError, match="out_stride"):
        build_model(base_cfg(decoder={"out_stride": 4}))


def test_encoder_input_size_must_be_a_multiple_of_14():
    with pytest.raises(ValueError, match="multiple of"):
        build_model(base_cfg(encoder={"input_size": 512}))


def test_unknown_encoder_name_fails_loudly():
    with pytest.raises(ValueError, match="Unknown encoder"):
        build_model(base_cfg(encoder={"name": "resnet50"}))
