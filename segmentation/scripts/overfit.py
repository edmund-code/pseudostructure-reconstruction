#!/usr/bin/env python
"""Sanity check: overfit a handful of dense patches and decode them.

    python scripts/overfit.py --config configs/data.yaml --steps 300 --gpu 1

This is the gate between "the code runs" and "the code learns". It trains on N patches with
augmentation OFF and no ignore mask, which makes the task memorisation. If the loss does not
collapse and the decode does not recover the instances, something is wired wrong — and it is
far cheaper to find that here than four hours into a real run.

Passing means: loss falls by >10x, and PQ against the training patches' own ground truth
exceeds --min-pq (default 0.80).
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import _bootstrap  # noqa: F401
import numpy as np
import torch

from kidney_panoptic.data.classes import ClassSpec
from kidney_panoptic.data.dataset import PanopticPatchDataset, denormalize_image
from kidney_panoptic.data.patching import find_patch_metas
from kidney_panoptic.eval.instance_metrics import evaluate_instances
from kidney_panoptic.losses.panoptic_loss import PanopticLoss
from kidney_panoptic.models.panoptic import build_model
from kidney_panoptic.postprocess.decode import decode_panoptic
from kidney_panoptic.utils.config import load_config, resolve_path
from kidney_panoptic.utils.seed import set_seed


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", required=True)
    p.add_argument("--patches", default=None)
    p.add_argument("--n-patches", type=int, default=4)
    p.add_argument("--steps", type=int, default=300)
    p.add_argument("--lr", type=float, default=1e-3)
    p.add_argument("--gpu", type=int, default=None)
    p.add_argument("--min-pq", type=float, default=0.80)
    p.add_argument("--out", default="results/overfit")
    p.add_argument("--set", nargs="*", dest="overrides", default=[])
    args = p.parse_args()

    overrides = list(args.overrides) + [
        "train.augment=false",       # memorisation, not generalisation
        "ignore_strategy=none",      # dense supervision on these few patches
    ]
    cfg = load_config(args.config, overrides)
    set_seed(int(cfg.get("seed", 1337)))
    device = f"cuda:{args.gpu}" if args.gpu is not None else (
        "cuda" if torch.cuda.is_available() else "cpu")
    spec = ClassSpec.from_config(cfg)

    patches_dir = Path(args.patches or resolve_path(cfg, cfg["paths"]["patches_dir"]))
    metas = find_patch_metas(patches_dir)
    dense = pick_dense(metas, args.n_patches)
    print(f"overfitting {len(dense)} dense patches: {[m.stem for m in dense]}")

    ds = PanopticPatchDataset(dense, cfg, spec, train=False)
    batch = collate([ds[i] for i in range(len(ds))], device)

    model = build_model(cfg, spec).to(device)
    crit = PanopticLoss(cfg, spec.n_classes).to(device)
    opt = torch.optim.AdamW(model.trainable_parameters(), lr=args.lr, weight_decay=0.0)
    print(f"trainable params: {model.n_trainable() / 1e6:.2f} M")

    model.train(); model.encoder.eval()
    first = None
    t0 = time.time()
    for step in range(args.steps):
        out = model(batch["image"])
        loss, parts = crit(out, batch)
        opt.zero_grad(set_to_none=True)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.trainable_parameters(), 5.0)
        opt.step()
        if first is None:
            first = parts["loss"]
        if step % 25 == 0 or step == args.steps - 1:
            print(f"  step {step:4d}  loss {parts['loss']:.4f}  "
                  f"sem {parts['loss_semantic']:.4f}  bnd {parts['loss_boundary']:.4f}  "
                  f"ctr {parts['center_mse']:.4f}", flush=True)
    last = parts["loss"]
    print(f"loss {first:.4f} -> {last:.4f}  ({first / max(last, 1e-9):.1f}x) "
          f"in {time.time() - t0:.0f}s")

    # ---- decode the memorised patches and score them ----
    out_dir = Path(args.out); out_dir.mkdir(parents=True, exist_ok=True)
    maps = model.predict_maps(batch["image"], amp=False)
    pqs = []
    for i in range(batch["image"].shape[0]):
        res = decode_panoptic(maps["semantic_prob"][i].cpu().numpy(),
                              maps["boundary_prob"][i].cpu().numpy(),
                              maps["center_heatmaps"][i].cpu().numpy(), cfg, spec)
        gt = batch["id_map"][i].cpu().numpy()
        m = evaluate_instances(gt, res.id_map)
        pqs.append(m["PQ"])
        print(f"  patch {i}: n_gt {m['n_gt']}  n_pred {m['n_pred']}  "
              f"PQ {m['PQ']:.3f}  AJI {m['AJI']:.3f}  split {m['split_rate']:.2f}")
        save_panel(out_dir / f"patch_{i}.png", batch, maps, res, i, spec)

    mean_pq = float(np.mean(pqs))
    ok = last < first / 10 and mean_pq >= args.min_pq
    with open(out_dir / "overfit.json", "w") as f:
        json.dump({"first_loss": first, "last_loss": last, "mean_pq": mean_pq,
                   "per_patch_pq": pqs, "steps": args.steps, "passed": bool(ok)}, f, indent=2)
    print(f"\nmean PQ {mean_pq:.3f}  ->  {'PASS' if ok else 'FAIL'}   (panels in {out_dir})")
    return 0 if ok else 1


def pick_dense(metas, n):
    """The n patches with the most annotated instances."""
    scored = []
    for m in metas:
        with open(m) as f:
            j = json.load(f)
        scored.append((j.get("n_instances", 0), j.get("fg_fraction", 0.0), m))
    scored.sort(key=lambda t: (-t[0], -t[1]))
    return [m for _, _, m in scored[:n]]


def collate(items, device):
    out = {}
    for k in ("image", "semantic", "boundary", "centers", "valid", "weight", "id_map"):
        out[k] = torch.stack([it[k] for it in items]).to(device)
    return out


def save_panel(path, batch, maps, res, i, spec):
    """RGB | GT instances | predicted instances | boundary probability."""
    import cv2

    rgb = denormalize_image(batch["image"][i].cpu().numpy())
    gt = batch["id_map"][i].cpu().numpy()
    pred = res.id_map
    pbnd = maps["boundary_prob"][i, 2].cpu().numpy()

    def colorize(lbl):
        rng = np.random.default_rng(0)
        lut = rng.integers(40, 255, size=(int(lbl.max()) + 2, 3)).astype(np.uint8)
        lut[0] = 0
        return lut[np.clip(lbl, 0, len(lut) - 1)]

    panel = np.concatenate([
        rgb, colorize(gt), colorize(pred),
        np.repeat((pbnd * 255).astype(np.uint8)[..., None], 3, axis=2),
    ], axis=1)
    cv2.imwrite(str(path), cv2.cvtColor(panel, cv2.COLOR_RGB2BGR))


if __name__ == "__main__":
    sys.exit(main())
