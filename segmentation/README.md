# Kidney panoptic segmentation

This subproject segments mouse-kidney H&E slides into proximal, distal, collecting-duct, and
glomerular instances. It combines a frozen OpenMidnight encoder with a trainable native-resolution
stem/decoder and local interior/boundary watershed decoding.

## Setup and tests

```bash
conda env create -f environment.yml
conda activate kidney-panoptic
pytest tests
```

The repository excludes slides, annotation GeoJSONs, extracted patches, and checkpoints. Use the
private manifest template in `data/` and the workflow guide at
[`../docs/workflows/segmentation.md`](../docs/workflows/segmentation.md).

Model selection is pooled class-agnostic PQ, not validation loss. D4 TTA is required for reported
inference. Current limitations include residual same-subtype merges and underestimated object
areas; see the architecture decisions in `../docs/architecture/`.
