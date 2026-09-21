# Private data directory

This directory is a mount point for authorized data; do not commit its contents. The only tracked
files are this guide and `catalog.csv`.

The active analysis expects the following layout beneath the path supplied as `--data-root`:

```text
tubule_by_gene/<sample>_tubule_by_gene_caleb.h5ad
human_mouse_hcop_fifteen_column.txt.gz
mouse_vs_human/pathway_gene_sets/<library>.json
Ctrl1A2_kept_tubules_labeled_fine.geojson, Ctrl1A4_kept_tubules_labeled_fine.geojson,
IR2A2_kept_tubules_labeled_fine.geojson, IR2A4_kept_tubules_labeled_fine.geojson
HUK1_COR1_v2.geojson, HUK1_MED1_v2.geojson
```

Use only the `*_kept_tubules_labeled_fine.geojson` mouse segmentations. These are already
quality controlled upstream: each feature carries `kept_tubule`, `doublet_tier1_flag`,
`coarse_class` and `segment_class`, and low-support structures have been removed. Notebook 1
consumes their geometry only. The retired `*_v4.geojson` and older processed segmentations are
not index-compatible even when their feature indices appear valid.

The mouse Visium HD outputs are not under this data root: they live in their own project
directory, passed through `PSEUDOSPACE_VISIUM_ROOT`. That directory must contain one subdirectory
per mouse sample (`Ctrl1A2`, `Ctrl1A4`, `IR2A2`, `IR2A4`), each holding the 2 um `square_002um`
outputs -- directly, nested under `outs/binned_outputs/`, or at the sample root. Notebook 1 checks
the barcode prefix (`s_002um_`) rather than trusting the path, because a sample directory can also
carry a coarser 8/16 um bin.

Human segmentations are the `*_v2.geojson` re-segmentations. `obs['feature_index']` is a positional
index into the exact GeoJSON file Part 1 enumerated, so never pair a segmentation with an H5AD built
from a different file — re-run notebook 1 after replacing one.

The human-versus-healthy-mouse workflow consumes `Ctrl1A2` and `Ctrl1A4` as healthy mouse
controls and `HUK1_COR1`/`HUK1_MED1` as human samples. The Visium matrix is named
`HUK1_MED`; the corresponding processed segmentation uses `HUK1_MED1`; keep that naming
translation in local manifests. The HCOP table is read-only input. Optional reviewed
overrides can be supplied as a private two-column CSV via
`PSEUDOSPACE_ORTHOLOG_OVERRIDES=...`.
