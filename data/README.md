# Private data directory

This directory is a mount point for authorized data; do not commit its contents. The only tracked
files are this guide and `catalog.csv`.

The active analysis expects the following layout beneath the path supplied as `--data-root`:

```text
tubule_by_gene/<sample>_tubule_by_gene_caleb.h5ad
human_mouse_hcop_fifteen_column.txt.gz
mouse_vs_human/pathway_gene_sets/<library>.json
Ctrl_1A2_v4.geojson, Ctrl_1A4_v4.geojson, IR_2A2_v4.geojson, IR_2A4_v4.geojson
HUK1_COR1_tubules_processed_caleb.geojson, HUK1_MED1_tubules_processed_caleb.geojson
```

Use only the `*_v4.geojson` mouse segmentations. The older processed segmentations are not
index-compatible even when their feature indices appear valid.

The human-versus-healthy-mouse workflow consumes `Ctrl1A2` and `Ctrl1A4` as healthy mouse
controls and `HUK1_COR1`/`HUK1_MED1` as human samples. The Visium matrix is named
`HUK1_MED`; the corresponding processed segmentation uses `HUK1_MED1`; keep that naming
translation in local manifests. The HCOP table is read-only input. Optional reviewed
overrides can be supplied as a private two-column CSV via
`PSEUDOSPACE_ORTHOLOG_OVERRIDES=...`.
