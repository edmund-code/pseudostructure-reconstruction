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

## Notebook 113 results

`analysis/notebooks/113_pt_evidence_model.ipynb` (logic `113.evidence_model.1`) writes to
`results/pt_evidence_model/`. It followed the protocol frozen in b629793 and ran in 56 s. Most of its frozen tests
fail, two of them because the protocol itself was wrong; what survives is a sharp ranking of the lines.

**Two flaws in the frozen design, found on the first run.**

1. **Three of the five AKI annotation sources cannot be held out.** The fold rule tests on genes that are positive
   in no training source, and transporters, positional markers and named positions are almost entirely nested
   inside the literature and disease sources. Held-out positives: transporters 0, positional markers 4, named
   positions 0. Only the literature programs (877), disease genes (395) and androgen sets (341) are testable, and
   the species atlas leaves only two.
2. **The leakage rule was too narrow.** It dropped a line only in the fold where its own source is held out, so the
   all-sources coefficient fit still scored the androgen line against the androgen sets. That line's coefficient was
   4.94 — by far the largest — and falls to **0.32** once the androgen sets are removed from the outcome. All
   coefficients below are from the corrected fit; the leaking values are kept in `coefficients.csv` beside them.
   The species atlas is unaffected: every handoff gene is also a transporter, so removing the handoff source
   changes nothing.

### Q1 · The lines do add to expression, where the question can be asked

| source | held-out positives | expression only | count model | full model | full − expression |
|---|---|---|---|---|---|
| AKI literature programs | 877 | 0.630 | 0.763 | 0.745 | **+0.115** |
| AKI disease genes | 395 | 0.517 | 0.553 | 0.549 | +0.032 |
| AKI androgen sets | 341 | 0.507 | 0.681 | 0.612 | **+0.106** |
| species literature programs | 912 | 0.685 | 0.743 | 0.695 | +0.010 |
| species disease genes | 420 | 0.630 | 0.572 | 0.599 | −0.031 |

- **Q1 fails as written** (it asked for three of five AKI sources at ≥ 0.05 AUC, and only three sources exist).
  **Two of the three testable sources clear the bar**, the AKI mean gain is +0.084, and the species mean is −0.011.
- The honest reading: on the AKI axis the evidence pattern predicts recognition well beyond expression; on the
  species axis it does not, where expression alone already reaches 0.63–0.69.

### Q2 · The plain count beats the fitted weighting — a useful negative

- **Q2 fails.** The count model (notebook 112's line count plus the same expression deciles) is above the fitted
  model on **four of the five** testable folds: AKI literature +0.018, AKI disease +0.004, AKI androgen +0.068,
  species literature +0.049; only species disease genes go the other way (−0.027). *(Corrected 2026-10-20: this
  paragraph said "three of the four" while the table above it lists all five folds. `fold_auc.csv` has five rows
  with a non-null AUC; the AKI disease-gene fold is testable, with 395 held-out positives, and the count model wins
  it too.)*
- Fitting nine weights on one set of curated sources does not transfer to another. **Notebook 112's unweighted
  count is the statistic to use**, which is the simpler recommendation and the one already in the scoreboard.

### Q3 · Which lines carry information about recognised biology

Coefficients from the corrected fit, with 1,000-gene-bootstrap intervals:

| line | calls | coefficient | 95% interval |
|---|---|---|---|
| A7 injury phase (Kirita as a ruler, 99) | 1,010 | **0.79** | 0.63 to 0.94 |
| A2 location screen (91) | 362 | **0.76** | 0.51 to 1.02 |
| A4 position-split selection (93) | 2,719 | 0.60 | 0.46 to 0.75 |
| A6 human AKI, KPMP (98) | 2,461 | 0.40 | 0.29 to 0.52 |
| A5 mouse IRI, Kirita 12 h (94) | 3,899 | 0.40 | 0.29 to 0.51 |
| A1 level screen (91) | 3,159 | 0.38 | 0.27 to 0.50 |
| A9 androgen instrument (111) | 195 | 0.32 | −0.11 to 0.69 |
| A8 variability (107) | 1,217 | 0.25 | 0.08 to 0.43 |
| **A3 positional F-test (71)** | **4,244** | **0.05** | **−0.09 to 0.19** |

- **The frozen prediction was wrong.** I predicted the external replication lines would carry the most. They sit in
  the middle (0.40). The phase assignment and the location screen carry twice as much.
- **The positional F-test carries no measurable information about whether a gene is recognised biology.** It makes
  the largest number of calls (4,244, the largest of any line) and its interval covers zero. Beside notebook 112,
  where it had the largest yield, this is the sharp result: **the F-test finds the most genes and says the least
  about which of them are real.**
- The two lines at the top are the two that place a change in *time* (phase) or in *space* (location) rather than
  measuring its size. Specificity about where and when beats volume.
- Species: F-test 0.88, zonation-class disagreement 0.64, atlas level 0.54, handoff 0.40, level screen 0.25,
  rewired 0.19 (interval crosses zero), reversed gradient −0.20 (crosses zero). The species F-test is the opposite
  of the AKI one: there it is the most informative line.

### Q4 · The relabeled control is weak by construction

- **Q4 fails as framed.** Of the four lines with stored relabelings, location collapses cleanly (362 calls → 0 and
  1, coefficient 0.76 → 0.00 and −0.05) and variability collapses (0.25 → −0.02 and −0.17), but the level line
  keeps 40 and 3 calls and its coefficient *rises* to 0.49 on R1. A logistic coefficient is not comparable across
  list sizes: a 40-gene list that happens to be well curated gets a large weight.
- **The scale-free version of the control holds but only just:** the relabeled model's mean held-out AUC is 0.627
  against the condition model's 0.635.
- **Why it is weak:** only four of the nine lines have stored relabelings, and the model's information comes mostly
  from the other five (phase, position-split, the two external lines, the F-test). A full control would need every
  notebook to store its relabeled statistic, which most do not.

### Q5 · The prioritised unanchored list

The corrected model scored on clean genes in no named set. Top of the AKI list, all falling after injury:

- **Peroxisomal and lipid**: Pipox, Far1, Hsd17b4, Cideb, Lpgat1, Chpt1, Acss2, Idh2.
- **Sulfur and one-carbon amino acid**: Cdo1, Cth, Gss, Aspa, Thnsl2, Aldh8a1.
- **Other**: Hps5, Bphl, Susd3, Tmem205, Tmem106b, Fmo2, Cdhr5, Fuca2, Csrp2.
- **Regulators**: Rorc, Zfp36, Lrrk2.

Top of the species list: Zbtb20, Acaa1a, Acss2, Igfbp4, Nfib, Rgl1, Selenop, App, Rhobtb1, Ctsb, Tmem64, Rtn4,
Bhmt2, Bicdl1 — eight of the fourteen carry a mouse gradient the human does not share.

The AKI list is the same catabolic program as notebook 112's strict core, which is a consistency check rather than
a new finding: the model ranks highly the genes that sit beside the curated members of a program without being
curated themselves. Post hoc, that is what "unanchored" means here — under-curated company, not new biology.

### What this answers

1. **Count, do not weight.** The unweighted convergence count of notebook 112 generalises across annotation sources
   better than a fitted model. This closes the question the atlas opened.
2. **Volume and information are different axes, and the positional F-test sits at opposite ends of them.** It
   delivers the most supported genes (notebook 112) and carries no information about which genes the field
   recognises (here). It should be used to find candidates, never to rank them.
3. **The most informative lines are the ones that localise a change** — the injury phase and the location screen —
   not the ones that measure its size, and not external replication.
4. **The species and AKI axes rank their methods differently.** The positional F-test is the least informative line
   on AKI and the most informative on human vs mouse, which matches the scoreboard's existing split (the F-test is
   the only calibrated species test and fails calibration on AKI) and is now measured on one scale.
5. **For the paper:** use the k ≥ 5 convergence core from notebook 112 as the gene table; quote the line
   coefficients here as the justification for which screens the paper runs; do not present a fitted gene score.
