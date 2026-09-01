# Private data directory

This directory is a mount point for authorized data; do not commit its contents. The only tracked
files are this guide and `catalog.csv`.

The active analysis expects the following layout beneath the path supplied as `--data-root`:

```text
tubule_by_gene/<sample>_tubule_by_gene_caleb.h5ad
mouse_vs_human/pathway_gene_sets/<library>.json
Ctrl_1A2_v4.geojson, Ctrl_1A4_v4.geojson, IR_2A2_v4.geojson, IR_2A4_v4.geojson
```

Use only the `*_v4.geojson` mouse segmentations. The older processed segmentations are not
index-compatible even when their feature indices appear valid.
