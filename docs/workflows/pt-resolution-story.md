# PT resolution story

`analysis/notebooks/14_pt_paper_figures_resolution_story.ipynb` is a code-driven exploratory notebook using saved results from notebooks 12 and 13. It validates the integrated pathway table, saved gene fits, common gene universe, and reviewed PT metadata, then ranks and freezes two conventionally supported and three spatial examples. It does not refit expression models, recompute coordinates, rerun DESeq2, or rerun matched-null enrichment inference.

Run the generated mirror in the analysis environment after notebook 12 has written its artifacts:

```bash
python analysis/pt_paper_figures_resolution_story.py --results-root results
```

The `PSEUDOSPACE_RESULTS_ROOT` fallback is supported. Tables and exploratory PDF, SVG, and PNG plots are written under `pt_paper_figures_resolution_story/` in that results root. Tidy plot inputs remain in `source_data/`; these private outputs must stay outside Git. The working notebook keeps its rendered cell outputs. For repository staging, use the output-free staging-copy workflow described in `AGENTS.md`.

Candidate selection retains its saved-evidence gates: original and common spatial enrichment, residual-correlation support, and at least 75% retention across all planned sensitivity runs. Spatial examples without native whole-PT or S1/S2/S3 signed-GSEA evidence also exclude native segment ORA and unsigned-magnitude calls. Selection stops if evidence changes or examples become redundant. Candidate ordering uses the observed common-universe spatial AUC effect as an exploratory evidence ordering; this adds no test.

The pathway comparison presents native signed-GSEA NES beside the common-universe gene-level `T_spatial` percentile-rank ECDF for pathway members and nonmembers. The saved matched-null q-value belongs to the rank-AUC enrichment test and is not computed from the plotted ECDF. NES is signed and directional; spatial rank-AUC is unsigned and tests enrichment of spatial heterogeneity against a covariate-matched null. Their effect sizes and nulls differ, so the plots do not imply that one result explains the other. Before plotting, the notebook recomputes each selected pathway's observed rank AUC with the same helper used by notebook 12 and asserts exact agreement (within floating-point tolerance) with the saved common-universe AUC.

Driver curves show three pathway members selected by their observed common-universe `T_spatial` values. They are illustrative fitted gene trajectories, not pathway enrichment evidence and not proof that aggregation caused a nonsignificant native GSEA result. The full selected pathway membership and per-gene rank tables are saved alongside conventional and framework summaries. All conclusions are exploratory with one human donor.
