# Which lines of evidence predict recognised biology? (notebook 113)

**The question.** Notebook 112 counted how many independent lines call each gene and showed that the count predicts
whether the gene is in a named literature, disease, transporter or hormone set: the anchored share rises from 0.21
at background to 0.62 at five lines and 0.78 at six. That is a validity check on the *count*. It does not say which
lines carry the information, and it treats a call by the androgen instrument and a call by the level screen as
worth the same.

This notebook replaces the count with a fitted model:

1. **Which lines predict recognised biology, and by how much?** A regression of "this gene is in a named set" on the
   line calls gives a weight per line. That is a measurement of each method's information content about real
   biology, on the same scale for every method — something neither the scoreboard's calibration yardsticks nor
   notebook 112's yield profile provides.
2. **Does the pattern of calls predict recognition out of sample?** The annotation sources are curated separately,
   so a model trained on four of them and tested on the fifth is not testing itself.
3. **Which unanchored genes does the model rank highest?** That is the prioritised list the paper can carry: genes
   that look like recognised biology by their evidence pattern but are in no named set.

**Why the outcome is legitimate.** "In a named set" is not truth; it is recognition, which is biased towards
abundant, old and human-disease-linked genes. Expression is therefore carried as a covariate in every model and as
the baseline every model must beat. A line that only reproduces the expression advantage will show no coefficient.

**The leakage rule.** A line built from the same source as the held-out annotation is dropped for that fold: the
androgen line with the androgen sets, the handoff line with the handoff annotation, the positional-marker
annotation with no line (none is derived from it). Without this the model would score a source against itself.

Frozen on 2026-10-19, before any model was fitted. The inputs are the two atlases committed in de4f0d1; no new
statistic is computed on the expression data.

## 1 · Inputs

- Notebook `analysis/notebooks/113_pt_evidence_model.ipynb`; logic `113.evidence_model.1`; seed 113.
- Results `results/pt_evidence_model/`.
- Inputs: `results/pt_gene_atlas/tables/aki_atlas.csv.gz` and `species_atlas.csv.gz` (notebook 112), plus the
  relabeled line definitions rebuilt as in 112 § C1 for the control.
- New helper `evidence_model` in `pseudospace/gene_atlas.py`, tested: an L2 logistic fit with the line calls and
  expression-decile indicators, returning coefficients and held-out scores.

## 2 · Design (fixed now)

- **Outcome:** membership of a named annotation set, one source at a time.
  - AKI sources: literature programs (79–82), disease genes (108), transporters (101), androgen sets (111),
    positional markers (105).
  - Species sources: literature programs, disease genes, transporters, handoffs (101), positional markers.
- **Predictors:** every line of that atlas, as 0/1, plus nine expression-decile indicators (the tenth is the
  reference). No interactions.
- **Model:** logistic regression, L2 penalty C = 1, `liblinear`, on the genes testable on at least one line.
  A gene not testable on a line enters that line as 0; the testable mask is carried as its own covariate per line
  so "not looked at" and "looked at and not called" are distinguished.
- **Folds:** leave-one-source-out. Train on the union of the other sources as the positive class, test on the
  held-out source's membership among genes not in any training source (so a gene cannot be positive in both).
- **Baselines, fitted on the same folds:** expression deciles alone, and the plain line count of notebook 112.

## 3 · Tests (fixed now)

**Q1 · Do the lines add information about recognised biology beyond expression?**
- Statistic: held-out AUC of the full model minus the expression-only model, per source.
- **Holds** when the full model is ≥ 0.05 AUC above expression-only in at least three of the five AKI sources, and
  the mean difference is positive on both atlases.

**Q2 · Does the fitted weighting beat the plain count?**
- Statistic: held-out AUC of the full model minus the count model (count plus expression deciles).
- **Holds** when the full model is above the count model on average across sources. A null or negative result means
  notebook 112's count is already sufficient, which is a useful answer.

**Q3 · Which lines carry the information?**
- Statistic: the coefficient of each line in the model fitted on all sources, with a 1,000-draw bootstrap over
  genes for its interval; reported beside each line's call count.
- **Prediction fixed now:** the external replication lines (Kirita, KPMP, atlas level) will carry the largest
  coefficients, because external replication is the closest available proxy for a real effect. If a screen line
  carries more than they do, that screen is finding recognised biology the external data do not.

**Q4 · The relabeled control.**
- Refit with the level, location, variability and androgen lines rebuilt from the balanced relabelings
  (notebook 112 § C1, which reproduces notebook 91's own counts).
- **Holds** when those lines' coefficients fall towards zero and the held-out AUC drops.

**Q5 · The prioritised unanchored list.**
- The model fitted on all sources, scored on genes in no named set, ranked. Report the top 50 per atlas with every
  line call, the injury or species statistics, and whether the gene passes the 104 and 109 filters.
- **Descriptive.** This list is the deliverable, not a test. It is read against the literature after the run and
  labelled as a post-hoc reading.

## 4 · Outputs and figures

- **Tables:** `fold_auc.csv`, `coefficients.csv`, `coefficients_relabeled.csv`, `ranked_unanchored_aki.csv`,
  `ranked_unanchored_species.csv`, `how_the_protocol_was_read.csv`, `verdicts.csv`, `run_record.json`.
- **Figures (plots only):**
  - 113a: held-out AUC per source for the three models, both atlases;
  - 113b: line coefficients with bootstrap intervals, condition against relabeled;
  - 113c: the recognition rate against the model score, in score deciles, with the count model beside it.

## 5 · Caveats

- Recognition is not truth. A gene can be real and uncurated; the model will rank it high only if its evidence
  pattern resembles curated genes'.
- The sources are not independent of each other: a disease gene is often also a transporter. Leave-one-source-out
  reduces but does not remove this.
- The species atlas has only seven lines, two of which are small (rewired 25, handoff 58), so its coefficients will
  be imprecise.
- Everything inherits the 2-vs-2 design and the single human donor.
