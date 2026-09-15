# Human versus healthy-mouse workflow

Run `analysis/human_vs_healthy_mouse.py` or the matching
`analysis/notebooks/03_human_vs_healthy_mouse.ipynb` from the repository root:

```bash
python analysis/human_vs_healthy_mouse.py \
  --data-root /path/to/private-data \
  --results-root /path/to/private-results
```

The workflow loads `Ctrl1A2` and `Ctrl1A4` as healthy mouse controls and
`HUK1_COR1`/`HUK1_MED1` as human samples. It performs species-aware QC, maps human symbols
to one-to-one mouse orthologs using HCOP, integrates the shared expression matrix, scores
nephron segment markers, reconstructs one shared DPT axis, and generates marker, gene-curve,
pathway-curve, and sample-support outputs under `human_vs_healthy_mouse/`.

## Review order

1. `ortholog_map_used.csv` and `human_gene_mapping_report.csv` — verify the mapping and any
   private override file.
2. `cohort_summary.csv` and `diagnostics/sample_support.csv` — check sample sizes, capture,
   and overlapping pseudospace support.
3. `cross_species_integration` metadata and `segment_cluster_marker_scores.csv` — inspect
   integration and marker resolution before interpreting DPT.
4. `heatmaps/`, `curves/gene_trajectory_comparison.csv`, and
   `curves/pathway_trajectory_comparison.csv` — these are descriptive comparisons.

## Interpretation limits

There is one human donor versus two mouse specimens. Sampling is confounded with species because
the human inputs are two slices of the *same* healthy cortex while the mouse controls are
mixed-kidney samples. `HUK1_MED1` is only named medulla at source: the processed segmentation is
cortex tissue, so it is not a medullary sample and must not be analysed or described as one (the
workflow assigns `region = 'cortex'` to both human slices). `HUK1_MED1` has weaker capture in the
current data. Do not report cellwise GAM F statistics or curve differences as confirmatory species
significance. Repeat the analysis after excluding `HUK1_MED1`, and use anatomically matched
mouse annotations if they become available.

The archived `legacy/notebooks/4_mouse_vs_human.ipynb` documents the earlier analysis but has
stale paths and a deprecated permutation engine. It is provenance only; use the current script
and notebook above.
