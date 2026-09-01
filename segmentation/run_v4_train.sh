#!/usr/bin/env bash
# v4 retrain: +Ctrl_1A2 (955 objects), other_structure dropped to ignore, anatomical names.
#
# Deliberately ONE experiment. The training-side variables are the new slide and the smaller
# label space; everything else that changed this round (min_area, interior_threshold,
# semantic_split off, TTA) is decode-time and applies to any checkpoint. The touching-channel
# loss rebalance is held back so the next run has a single attributable variable.
set -euo pipefail
cd "$(dirname "$0")"

conda run --no-capture-output -n tubule_segmentation python -u scripts/train.py \
  --config configs/data.yaml \
  --out runs/v4_panoptic \
  --gpu 1 \
  --max-val-patches 160
