# PT resolution story

`analysis/notebooks/14_pt_paper_figures_resolution_story.ipynb` is a code-driven exploratory notebook using saved results from notebooks 12 and 13. It validates the integrated pathway table, saved gene fits, common gene universe, and reviewed PT metadata, then ranks and freezes two conventional-supported and three spatial examples. It does not refit expression models, recompute coordinates, DESeq2, or enrichment.

Run the generated mirror in the analysis environment after notebook 12 has written its artifacts:

```bash
python analysis/pt_paper_figures_resolution_story.py --results-root results
```

The `PSEUDOSPACE_RESULTS_ROOT` fallback is supported. Tables and exploratory PDF, SVG, and PNG plots are written under `pt_paper_figures_resolution_story/` in that results root. Tidy plot inputs remain in `source_data/`; these private outputs must stay outside Git. The working notebook keeps its rendered cell outputs. For repository staging, use the output-free staging-copy workflow described in `AGENTS.md`.

Candidate selection retains its saved-evidence gates: original and common spatial enrichment, residual-correlation support, and at least 75% retention across all planned sensitivity runs. Spatial examples without native whole-PT or S1/S2/S3 signed-GSEA evidence also exclude native segment ORA and unsigned-magnitude calls. Selection stops if evidence changes or examples become redundant.

The pathway mean curves average saved member-gene fits. Pointwise HC3 bands average residuals before variance estimation, retaining cross-gene residual covariance; they are conditional on observed structures and the coordinate, not independent-donor confidence intervals. Segment guides are descriptive midpoints between specimen-balanced label medians; distributions overlap. Same-curve compression summaries illustrate attenuation and cancellation without attributing native q-value differences to aggregation alone.

Whole-PT, coarse segment, and pseudostructure evidence target different questions, so their significance sets need not nest. Notebook 12's level, spatial, total, signed-GSEA, and conventional FDR conventions are retained. The controlled discrete-versus-continuous benchmark is supporting robustness evidence, not the central conventional baseline. All conclusions are exploratory with one human donor.
