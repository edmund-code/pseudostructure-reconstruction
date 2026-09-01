# ADR 0002 — Full-resolution RGB stem with skips into the decoder

**Status:** accepted
**Date:** 2026-07-29

## Context

The encoder is a frozen DINOv2 ViT-g/14 (OpenMidnight pathology weights). Its finest
spatial unit is one 14×14 patch token. At the working resolution of 0.44 µm/px that token
covers **6.2 µm**.

The wall between two touching tubules is thinner than that. So is the capsule between
adjacent glomeruli. A decoder reading only token features is being asked to place a
boundary inside a cell it cannot see into — the evidence was averaged away before the
decoder ever ran.

The predecessor decoded from a 64×64 feature grid for a 1024 px patch (stride 16) and
upsampled. Its boundary placement error was bounded below by that stride regardless of how
well it was trained.

## Decision

A trainable **convolutional stem runs on the native patch pixels** and produces a pyramid at
strides 1, 2, 4, 8, 16. The frozen token map is projected, resized to the stride-16 level,
and fused there; the decoder then upsamples through the stem pyramid with a **skip at every
level, including stride 1**. Output stride is 1 by default (`decoder.out_stride`, 1 or 2).

## Why

The two paths carry different information and both are needed:

- **Frozen tokens: semantics.** Pathology-pretrained features are what make a 5-class
  distinction learnable from a few hundred annotated instances. They are not localisable
  below 14 px.
- **Stem: geometry.** Where exactly the wall is. Trivially learnable from pixels, and it
  needs no pretraining — edges and texture gradients are what the first two conv layers of
  anything learn.

Fusing them means the semantic head can say *what* while the boundary head says *where*,
each from the source that actually carries that information.

Concat-then-project (rather than summing) is used for the multi-layer token fusion so the
decoder learns how much each encoder block contributes; the last block of a ViT is
specialised toward the pretraining objective and mid-depth blocks retain more localisable
texture.

## Consequences

- Memory and compute rise: stride-1 activations at 512² dominate the decoder cost. It fits
  comfortably in 48 GB at batch 8, and `out_stride: 2` halves it if needed.
- The stem is trained from scratch on a small dataset. It is deliberately small
  (~1.3 M parameters) and heavily augmented; a large stem here would overfit before the
  frozen features contributed anything.
- `tests/test_model.py::test_gradients_reach_the_stem_but_not_the_encoder` pins the wiring:
  gradient must reach the stride-1 stem convolution and must not reach the encoder.
- The encoder still resamples 512 → 518 to hit a multiple of 14. That 1.2 % rescale affects
  only the token path; the stem sees native pixels.
