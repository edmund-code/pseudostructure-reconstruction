#!/usr/bin/env python
"""Train the panoptic model.

    python scripts/train.py --config configs/data.yaml --out runs/v2_panoptic --gpu 1

Checkpoint selection is by validation **PQ**, not validation loss. Each ignore strategy
normalises its loss over its own valid region, so loss is not even comparable between
configurations, let alone a proxy for instance quality — a model can lower its loss by
smoothing boundaries and simultaneously fuse every pair of touching tubules
(docs/adr/0004).
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from collections import Counter
from pathlib import Path

import _bootstrap  # noqa: F401
import numpy as np
import torch
from torch.utils.data import DataLoader

from kidney_panoptic.data.classes import ClassSpec
from kidney_panoptic.data.dataset import PanopticPatchDataset
from kidney_panoptic.data.manifest import load_manifest
from kidney_panoptic.data.patching import find_patch_metas, load_patch
from kidney_panoptic.eval.instance_metrics import evaluate_instances
from kidney_panoptic.losses.panoptic_loss import PanopticLoss
from kidney_panoptic.losses.seg_losses import class_weights_from_counts
from kidney_panoptic.models.panoptic import build_model
from kidney_panoptic.postprocess.decode import decode_panoptic
from kidney_panoptic.utils.checkpoint import save_checkpoint
from kidney_panoptic.utils.config import load_config, resolve_path
from kidney_panoptic.utils.seed import set_seed, worker_init_fn


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", required=True)
    p.add_argument("--out", required=True, help="run directory for checkpoints and logs")
    p.add_argument("--patches", default=None, help="overrides paths.patches_dir")
    p.add_argument("--manifest", default=None)
    p.add_argument("--gpu", type=int, default=None)
    p.add_argument("--resume", default=None)
    p.add_argument("--max-val-patches", type=int, default=64,
                   help="cap on val patches decoded per evaluation (decode is the slow part)")
    p.add_argument("--set", nargs="*", dest="overrides", default=[])
    return p.parse_args()


def main():
    args = parse_args()
    overrides = list(args.overrides)
    if args.gpu is not None:
        overrides.append(f"device=cuda:{args.gpu}")
    cfg = load_config(args.config, overrides)
    set_seed(int(cfg.get("seed", 1337)))

    device = cfg.get("device", "cuda")
    if device.startswith("cuda") and not torch.cuda.is_available():
        print("CUDA unavailable — falling back to CPU.")
        device = "cpu"

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    spec = ClassSpec.from_config(cfg)

    # ------------------------------------------------------------------ data
    patches_dir = Path(args.patches or resolve_path(cfg, cfg["paths"]["patches_dir"]))
    metas = find_patch_metas(patches_dir)
    if not metas:
        raise SystemExit(f"No patches under {patches_dir}. Run scripts/extract_patches.py.")

    manifest_path = args.manifest or cfg.get("manifest")
    rows = load_manifest(resolve_path(cfg, manifest_path), check_paths=False)
    split_of = {r["slide_id"]: r["split"] for r in rows}
    train_metas = [m for m in metas if split_of.get(m.parent.name, "train") == "train"]
    val_metas = [m for m in metas if split_of.get(m.parent.name, "train") == "val"]
    if not train_metas:
        raise SystemExit("No training patches (check the manifest `split` column).")
    if not val_metas:
        print("WARNING: no val patches — PQ selection will fall back to the training set.")
        val_metas = train_metas
    val_metas = build_val_set(val_metas, args.max_val_patches)
    print(f"train patches: {len(train_metas)}   val patches: {len(val_metas)} "
          f"({sum(1 for m in val_metas if _n_instances(m) > 0)} with annotations)")

    tcfg = cfg["train"]
    train_ds = PanopticPatchDataset(train_metas, cfg, spec, train=True,
                                    seed=int(cfg.get("seed", 1337)))
    val_ds = PanopticPatchDataset(val_metas, cfg, spec, train=False)
    print(f"train items (after oversampling): {len(train_ds)}")

    loader = DataLoader(
        train_ds, batch_size=int(tcfg["batch_size"]), shuffle=True,
        num_workers=int(tcfg.get("num_workers", 8)), pin_memory=(device != "cpu"),
        drop_last=len(train_ds) > int(tcfg["batch_size"]), worker_init_fn=worker_init_fn,
        persistent_workers=int(tcfg.get("num_workers", 8)) > 0,
    )

    # ----------------------------------------------------------------- model
    model = build_model(cfg, spec).to(device)
    print(f"trainable parameters: {model.n_trainable() / 1e6:.2f} M "
          f"(encoder frozen: {cfg['encoder']['name']})")

    weights = compute_class_weights(train_metas, spec, cfg).to(device)
    print("semantic class weights: " +
          ", ".join(f"{n}={w:.2f}" for n, w in zip(spec.names, weights.tolist())))
    criterion = PanopticLoss(cfg, spec.n_classes, weights).to(device)

    decay, no_decay = [], []
    for n, p in model.named_parameters():
        if not p.requires_grad:
            continue
        (no_decay if p.ndim <= 1 else decay).append(p)
    opt = torch.optim.AdamW(
        [{"params": decay, "weight_decay": float(tcfg.get("weight_decay", 0.01))},
         {"params": no_decay, "weight_decay": 0.0}],
        lr=float(tcfg["lr"]),
    )
    scaler = torch.amp.GradScaler("cuda", enabled=bool(tcfg.get("amp", True)) and device != "cpu")

    epochs = int(tcfg["epochs"])
    steps_per_epoch = max(1, len(loader))
    start_epoch = 0
    best_pq, best_epoch = -1.0, -1
    tie_eps = float(tcfg.get("select_tie_eps", 0.0))
    tie_break = str(tcfg.get("select_tie_break", "merge_rate"))
    best_tie = float("inf")
    history = []

    if args.resume:
        ck = torch.load(args.resume, map_location=device, weights_only=False)
        # A vocabulary change that keeps n_classes the same loads SILENTLY and relabels
        # every prediction; strict=False would also quietly skip a head whose shape moved
        # (e.g. the center head when compact_classes shrinks), leaving it randomly
        # initialised with no warning. Both are invisible at runtime, so refuse instead.
        if list(ck.get("classes", [])) != list(spec.names):
            raise SystemExit(
                f"--resume checkpoint was trained on classes {ck.get('classes')}\n"
                f"but this config declares                   {list(spec.names)}.\n"
                f"Resuming would silently relabel every prediction. Train from scratch."
            )
        model.load_state_dict(ck["model_state"], strict=False)
        if "optimizer_state" in ck:
            opt.load_state_dict(ck["optimizer_state"])
        start_epoch = int(ck.get("epoch", 0)) + 1
        best_pq = float(ck.get("metrics", {}).get("PQ", -1.0))
        print(f"resumed from {args.resume} at epoch {start_epoch} (best PQ {best_pq:.4f})")

    # -------------------------------------------------------------- training
    for epoch in range(start_epoch, epochs):
        model.train()
        model.encoder.eval()          # belt and braces: the encoder never trains
        t0 = time.time()
        agg = Counter()
        n_steps = 0

        for step, batch in enumerate(loader):
            lr = lr_at(epoch + step / steps_per_epoch, tcfg)
            for g in opt.param_groups:
                g["lr"] = lr

            batch = {k: (v.to(device, non_blocking=True) if torch.is_tensor(v) else v)
                     for k, v in batch.items()}
            with torch.autocast("cuda", dtype=torch.float16,
                                enabled=scaler.is_enabled()):
                out = model(batch["image"])
                loss, parts = criterion(out, batch)

            opt.zero_grad(set_to_none=True)
            scaler.scale(loss).backward()
            if tcfg.get("grad_clip"):
                scaler.unscale_(opt)
                torch.nn.utils.clip_grad_norm_(model.trainable_parameters(),
                                               float(tcfg["grad_clip"]))
            scaler.step(opt)
            scaler.update()

            for k, v in parts.items():
                agg[k] += v
            n_steps += 1
            if step % 20 == 0:
                print(f"  e{epoch:03d} s{step:04d}/{steps_per_epoch}  "
                      f"loss {parts['loss']:.4f}  sem {parts['loss_semantic']:.4f}  "
                      f"bnd {parts['loss_boundary']:.4f}  ctr {parts['center_mse']:.4f}  "
                      f"lr {lr:.2e}", flush=True)

        means = {k: v / max(n_steps, 1) for k, v in agg.items()}
        line = {"epoch": epoch, "lr": lr, "secs": round(time.time() - t0, 1), **means}

        if (epoch + 1) % int(tcfg.get("eval_every", 2)) == 0 or epoch == epochs - 1:
            val = validate(model, val_ds, cfg, spec, device, args.max_val_patches)
            line.update({f"val_{k}": v for k, v in val.items()})
            print(f"  [epoch {epoch}] val PQ {val['PQ']:.4f}  SQ {val['SQ']:.4f}  "
                  f"DQ {val['DQ']:.4f}  AJI {val['AJI']:.4f}  "
                  f"split {val['split_rate']:.3f}  "
                  f"merge {val['merge_rate']:.3f} ({val['n_merged_gt']}/{val['n_gt']})  "
                  f"n_pred {val['n_pred']}", flush=True)
            # PQ remains the selection metric (ADR 0004). But PQ differences inside the
            # measurement noise of a few hundred validation objects are not real, and the
            # previous run sat in a 0.745-0.748 band for ten epochs. Within that band,
            # break the tie on merge_rate rather than on whichever epoch happened to
            # sample higher. This changes how ties are resolved, not what is optimised —
            # a merge PENALTY in the metric would be wrong, since DQ already charges a
            # merge twice (2 FN + 1 FP).
            better = val["PQ"] > best_pq + tie_eps or (
                val["PQ"] > best_pq - tie_eps and val.get(tie_break, 1.0) < best_tie - 1e-9
            )
            if better:
                best_pq = max(best_pq, val["PQ"])
                best_tie, best_epoch = val.get(tie_break, 1.0), epoch
                save_checkpoint(out_dir / "best_pq.pt", model, cfg, spec,
                                epoch=epoch, metrics=val, optimizer=opt)
                print(f"  ** new best: PQ {val['PQ']:.4f} {tie_break} "
                      f"{val.get(tie_break, float('nan')):.4f} (epoch {epoch}) -> best_pq.pt")

        history.append(line)
        save_checkpoint(out_dir / "last.pt", model, cfg, spec, epoch=epoch,
                        metrics={"PQ": best_pq}, optimizer=opt)
        if int(tcfg.get("save_every", 0)) and (epoch + 1) % int(tcfg["save_every"]) == 0:
            save_checkpoint(out_dir / f"epoch_{epoch:03d}.pt", model, cfg, spec, epoch=epoch)
        with open(out_dir / "history.json", "w") as f:
            json.dump(history, f, indent=2)
        print(f"[epoch {epoch}] loss {means.get('loss', 0):.4f}  "
              f"({line['secs']}s)", flush=True)

    print(f"\nDone. Best val PQ {best_pq:.4f} at epoch {best_epoch} -> "
          f"{out_dir / 'best_pq.pt'}")
    if best_epoch < 0:
        print("WARNING: no evaluation ever ran; best_pq.pt was not written.")


# --------------------------------------------------------------------- helpers
def _n_instances(meta_path: Path) -> int:
    with open(meta_path) as f:
        return int(json.load(f).get("n_instances", 0))


def build_val_set(metas, cap: int, annotated_frac: float = 0.75) -> list[Path]:
    """A validation subset that can actually resolve PQ.

    Only ~16% of the extracted patches on a val slide carry any annotation, so taking the
    first N in filename order would compute PQ from a handful of objects. Annotated patches
    are taken first, then background-only ones are added to make up the rest — the
    background ones are not filler, they are where a hallucinated instance shows up as a
    false positive and is charged against DQ.
    """
    annotated = sorted((m for m in metas if _n_instances(m) > 0), key=str)
    empty = sorted((m for m in metas if _n_instances(m) == 0), key=str)
    n_ann = min(len(annotated), max(1, int(round(cap * annotated_frac))))
    picked = _spread(annotated, n_ann) + _spread(empty, max(0, cap - n_ann))
    return picked or list(metas)[:cap]


def _spread(items: list, n: int) -> list:
    """Evenly spaced deterministic subsample — keeps spatial coverage of the slide."""
    if n <= 0 or not items:
        return []
    if n >= len(items):
        return list(items)
    idx = np.linspace(0, len(items) - 1, n).round().astype(int)
    return [items[i] for i in dict.fromkeys(idx.tolist())]


def lr_at(epoch_f: float, tcfg: dict) -> float:
    base = float(tcfg["lr"])
    min_lr = float(tcfg.get("min_lr", 0.0))
    warm = float(tcfg.get("warmup_epochs", 0))
    total = float(tcfg["epochs"])
    if warm > 0 and epoch_f < warm:
        return base * (epoch_f + 1e-8) / warm
    p = min(max((epoch_f - warm) / max(total - warm, 1e-6), 0.0), 1.0)
    if tcfg.get("schedule", "cosine") == "poly":
        return min_lr + (base - min_lr) * (1 - p) ** 0.9
    return min_lr + 0.5 * (base - min_lr) * (1 + np.cos(np.pi * p))


def compute_class_weights(metas, spec: ClassSpec, cfg: dict, sample: int = 300):
    """Inverse-frequency-style weights from instance PIXEL counts over a patch sample."""
    lc = cfg.get("loss", {}).get("semantic", {})
    mode = lc.get("class_weight_mode", "inverse_sqrt_freq")
    counts = np.zeros(spec.n_classes, dtype=np.float64)
    rng = np.random.default_rng(0)
    idx = rng.permutation(len(metas))[:sample]
    for i in idx:
        p = load_patch(metas[i])
        id_map, cmap = p["id_map"], p["class_idx_of_id"]
        counts[0] += float((id_map == 0).sum())
        for iid, ci in cmap.items():
            counts[ci] += float((id_map == iid).sum())
    return class_weights_from_counts(torch.from_numpy(counts), mode,
                                     float(lc.get("class_weight_clip", 12.0)))


@torch.no_grad()
def validate(model, val_ds, cfg, spec, device, max_patches: int) -> dict:
    """Decode val patches and pool instance metrics.

    Metrics are POOLED (summed TP/FP/FN, area-weighted) rather than averaged per patch: a
    512 px patch may hold three objects, and a mean over per-patch PQ lets a patch with one
    object swing the number as hard as a patch with thirty.
    """
    model.eval()
    # `decode.min_area` is calibrated on WHOLE objects, for a WSI decode where nothing is
    # truncated. Validation decodes 512 px patches, where a large share of objects is cut
    # by the patch border: with a median tubule ~100 px across, every object whose centre
    # sits within ~50 px of an edge comes out undersized. Applying the slide-scale floor
    # here would delete those predictions while their GT counterparts still count as
    # misses, inflating FN in a way that has nothing to do with model quality.
    #
    # So validation decodes with its own, smaller floors. This does NOT relax what ships —
    # scripts/eval_regions.py and the WSI decode both use the calibrated values. It keeps
    # the epoch-ranking signal clean, which is all validate() is for.
    vcfg = dict(cfg)
    scale = float(cfg["train"].get("val_min_area_scale", 0.25))
    if scale != 1.0:
        dec = dict(cfg.get("decode", {}) or {})
        dec["min_area"] = {k: max(1, int(v * scale))
                           for k, v in (dec.get("min_area", {}) or {}).items()}
        vcfg = {**cfg, "decode": dec}
    cfg = vcfg

    n = len(val_ds)
    tp = fp = fn = 0
    sq_sum, aji_sum, split_sum, gt_sum, n_pred = 0.0, 0.0, 0, 0, 0
    merge_sum = merging_sum = 0

    for i in range(n):
        item = val_ds[i]
        img = item["image"].unsqueeze(0).to(device)
        maps = model.predict_maps(img, amp=bool(cfg["train"].get("amp", True)))
        res = decode_panoptic(
            maps["semantic_prob"][0].cpu().numpy(),
            maps["boundary_prob"][0].cpu().numpy(),
            maps["center_heatmaps"][0].cpu().numpy(),
            cfg, spec, tissue_mask=None,
        )
        gt = item["id_map"].numpy()
        valid = item["valid"].numpy() > 0.5
        m = evaluate_instances(gt, res.id_map, valid=valid)
        tp += m["TP"]; fp += m["FP"]; fn += m["FN"]
        sq_sum += m["SQ"] * m["TP"]
        aji_sum += m["AJI"] * max(m["n_gt"], 1)
        split_sum += m["n_split"]
        merge_sum += m["n_merged_gt"]
        merging_sum += m["n_merging_preds"]
        gt_sum += m["n_gt"]
        n_pred += m["n_pred"]

    sq = sq_sum / tp if tp else 0.0
    dq = tp / (tp + 0.5 * fp + 0.5 * fn) if (tp + fp + fn) else 0.0
    return {
        "PQ": sq * dq, "SQ": sq, "DQ": dq,
        "AJI": aji_sum / max(gt_sum, 1),
        "split_rate": split_sum / max(gt_sum, 1),
        # Merging is the dominant failure mode; raw counts travel with the rate because
        # on ~96 val patches the denominator is small enough that the rate alone hides
        # its own noise.
        "merge_rate": merge_sum / max(gt_sum, 1),
        "n_split": split_sum, "n_merged_gt": merge_sum, "n_merging_preds": merging_sum,
        "TP": tp, "FP": fp, "FN": fn, "n_gt": gt_sum, "n_pred": n_pred,
        "n_patches": n,
    }


if __name__ == "__main__":
    sys.exit(main())
