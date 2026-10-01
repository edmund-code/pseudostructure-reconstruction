# PT manuscript resolution figures

`analysis/notebooks/14_pt_paper_figures_resolution_story.ipynb` generates five
figures from notebook 12's saved results and notebook 13's reviewed PT metadata.
It does not recompute coordinates, expression fits, DESeq2, or enrichment.
The generated mirror is `analysis/pt_paper_figures_resolution_story.py`.

Run in the analysis environment, after notebook 12 has written the integrated
pathway classification, atlas, sensitivity results, signed level comparison,
and gene trajectory archive:

```bash
python analysis/pt_paper_figures_resolution_story.py --results-root results
```

The usual `PSEUDOSPACE_RESULTS_ROOT` fallback is supported. Outputs stay under
`pt_paper_figures_resolution_story/` in that results root: PDF, editable SVG,
600-dpi PNG, legends, a figure index, input fingerprints, candidate rankings,
selected examples, and tidy plotted data in `source_data/`. These are private
results and must remain outside Git.

The notebook ranks all common-universe pathways before freezing two shared and
three conventionally undetected examples. Both classes require original and
common spatial enrichment, residual-correlation support, and retention in at
least 75% of all planned sensitivity runs. The second class also excludes native
whole-PT and segment signed-GSEA calls, native species-within-segment ORA calls,
and unsigned magnitude calls. Selection stops if the saved evidence changes or
examples duplicate a program or substantially overlap in genes.

Continuous pathway curves average saved member-gene fits. Pointwise HC3 bands
retain cross-gene residual covariance by averaging residuals before forming
variance; they remain conditional on observed structures and the coordinate,
not independent-donor confidence intervals. Segment transition guides are
descriptive midpoints between specimen-balanced label medians. The actual
segment distributions overlap and are plotted separately. Same-curve compression
measurements illustrate attenuation and cancellation without claiming that
aggregation alone explains differences between native statistical pipelines.

The information-profile figure keeps level, spatial and total tests visible;
screen unions are descriptive, not jointly FDR controlled. The cluster-only
diagnostic and matched discrete benchmark are confined to the supplement.

When the Nature figure auditor scripts are installed, every export runs panel
alignment, PDF glyph, and rendered collision checks. Their directory can be set
with `PSEUDOSPACE_FIGURE_QA_SCRIPTS`; an unavailable auditor is explicitly recorded.
Final visual inspection remains necessary. Keep working notebook outputs;
stage output-free copies using the repository staging tool.
