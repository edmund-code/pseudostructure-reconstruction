# Architecture overview

## Data flow

```text
Private Visium HD outputs + v4 tubule GeoJSON
  -> analysis/notebooks/01_segmentation_to_gene_matrix.ipynb
  -> private tubule-by-gene H5AD matrices
  -> analysis/mouse_only_pseudospace.py
  -> private Harmony/DPT/result outputs + curated tables in docs/results/
```

The active analysis uses only the four mouse samples (`Ctrl1A2`, `Ctrl1A4`, `IR2A2`, `IR2A4`).
The independent `segmentation/` project trains the boundary-based panoptic model that creates the
v4 GeoJSON inputs. Its public package remains `kidney_panoptic`.

## Design boundaries

- `pseudospace/` is reusable, data-location agnostic analysis logic.
- `analysis/` contains workflow orchestration, notebook entry points, and utilities that require
  the study-specific data contract.
- `segmentation/` owns its own Python environment because its Torch/CUDA/whole-slide stack is
  intentionally independent of Scanpy and R/Harmony.
- `legacy/` is preserved solely for provenance. It is unsupported and must not be mixed with the
  current mouse-only results.

## Reproducibility model

Synthetic tests verify core behavior in GitHub CI. Full scientific reproduction requires
authorized source data and is intentionally external to Git. This separation prevents both
accidental publication of restricted material and uncloneable repository history.
