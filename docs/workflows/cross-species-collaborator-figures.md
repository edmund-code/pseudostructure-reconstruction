# Cross-species collaborator figures

The companion exporter reads the completed notebook 03 results and makes a small, annotated
figure pack. It does not change the analysis notebook, its cached fits, or its result tables.
Generated figures, source data, and the illustrated methods guide remain under the private
results directory and must not be committed.

## Figure contract

Audience: collaborators assessing descriptive human–mouse PT expression patterns. The biological
units are two mouse specimens and one human donor with two cortex slices. Tubules, genes, and
pseudospace grid points are not additional biological replicates.

The four quantitative grids answer separate questions:

1. **Metabolic expression:** do selected metabolic pathway scores differ, and do multiple member
   genes support the direction? Show fitted pathway curves, individual specimens, and the full
   distribution of member-gene level effects. Selection is a stated thematic shortlist, not an
   exhaustive or independent set of discoveries.
2. **Gene gradients:** can an overall offset coexist with an opposite spatial gradient? Show the
   same two illustrative, previously identified genes on their expression scale and after
   mean-centering only. Keep individual specimens visible and retain amplitude.
3. **Response modules:** what shape is shared by genes in the modules associated with inflammatory
   and metabolic gene sets? Show median standardized, centered difference curves with the
   interquartile spread across genes, then observed versus background-expected pathway membership.
   This spread is not a confidence interval. Module enrichment is exploratory.
4. **Robustness:** which findings survive weighting, normalization, and specimen omission? Separate
   stable level/rank summaries from the more sensitive module partition. Do not use the duplicated
   reference rows of the saved module-stability table as a pooled-versus-balanced comparison.

Backend: Python/matplotlib, matching the notebook. Exports: slide-sized PDF and SVG with editable
text, 300-dpi PNG previews, a combined PDF brief, an illustrated HTML guide, and per-figure CSV
source data. No microscopic images or AI-generated data are used. Each multi-panel plot must pass
the rendered panel-alignment gate; PDFs receive font-size and collision audits and visual review.

## Run

From the analysis environment, after notebook 03 has completed:

```bash
python analysis/scripts/export_cross_species_collaborator_figures.py \
  --results-root <results-root> --qa-tools <nature-figure-skill/scripts>
```

`--results-root` defaults to `PSEUDOSPACE_RESULTS_ROOT`, then the repository's `results` directory.
The default output is `human_vs_healthy_mouse/collaborator_figures` beneath that root.
The optional `--qa-tools` supplies the figure skill's rendered alignment auditor; its PDF audit
tools should also be run on the exports. The exporter copies its own source into the output
directory so the private deliverable remains reproducible.

The exporter reconstructs specimen curves with the existing GAM implementation and verifies them
against saved balanced effects and saved pathway curves before plotting. It selects a cached fit
only when its complete gene-level effects match the current result table. It never selects a cache
by modification time alone. If a prerequisite is absent or inconsistent, it stops.

## What the pathway ordinate means

For structure `i` and gene `g`, the notebook starts with

`x[i,g] = ln(1 + 10000 * count[i,g] / retained_panel_total[i])`.

Each gene is standardized using its mean and population standard deviation across **all pooled
mouse and human PT structures**. A pathway score is the arithmetic mean of those gene z-scores
over its available, tested members. The displayed curve is the fitted GAM expectation of that
score along pseudospace.

Use the label **Relative pathway expression (mean gene z-score)**. The old label, “fitted module
z-score”, is ambiguous: these are predefined pathway gene sets, and the pathway score itself has
not been standardized to unit variance. Zero is the pooled member-gene reference, not no
expression. A score of one is not a pathway standard deviation, a fold change, or a p-value.

## Testing boundaries

- Gene and pathway GAM contrasts are descriptive. Cellwise F statistics do not provide calibrated
  species p-values for this cohort.
- Response-module annotations use a hypergeometric overlap test against discovery-eligible genes,
  with BH correction across all tested module–pathway pairs. Data-driven discovery and dependent
  genes limit confirmatory interpretation.
- Signed-ranking enrichment uses the custom competitive normal approximation implemented in
  `pseudospace.enrichment`: a background-centered set mean, finite-population correction, and a
  variance factor based on squared residual correlations. BH correction covers all pooled and
  balanced questions jointly. This is not the validated CAMERA implementation.
- Specimen-omission checks assess sensitivity, not population uncertainty. Recomputed DPT uses the
  existing Harmony representation and an anchored root; it is not independent anatomical validation.
- Neither total-count normalization nor pooled z-scoring identifies absolute RNA abundance or
  pathway activity. Both human slices are cortex, and their common donor remains confounded with
  species and assay characteristics.
