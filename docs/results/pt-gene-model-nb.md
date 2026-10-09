# Gene model: negative binomial on raw counts

Decided by the user on 2026-10-09. This note freezes the specification before notebooks 37 onward
are re-run. It supersedes the Gaussian gene models described in SM12 of the supplementary information.

## Why the model changed

Until now every gene model was a weighted Gaussian regression spline on log1p(count / library × 10⁴).
A scratch pilot compared it with a negative-binomial (NB) GLM on raw counts. The pilot ran on the
private data, so this switch is **post hoc** and is recorded as a deviation in Supplementary Note 8.
Notebook 65 reproduces the pilot from tracked code.

**Pilot results (SCF13, species split and the two balanced relabelings):**

| Method | Gaussian on log1p | NB on raw counts |
|---|---|---|
| T_spatial + joint test (A6·JNT) | 60 calls; relabeled 2, 0; calibration 0.060 | 31; 0, 0; 0.067 |
| S1/S2/S3 steps + joint test (A4·JNT) | 50; 0, 0; 0.064 | 15; 1, 1; 0.079 |
| Continuous adds over steps / steps add over continuous | 13 / 3 | 16 / 0 |

**Why NB is better.** On the log1p scale, a constant species ratio is not a constant difference.
It bends where counts are low and where a gene is nearly absent in one species. The Gaussian
T_spatial therefore scores some offsets as shape:
- **Simulation.** Genes given a pure constant ratio reached the Gaussian top 10% at 10.1%, and
  27% in the largest-ratio tertile. Under NB the rate was 0%. The simulation draws counts from the
  NB model, so it favours NB by construction.
- **Real data.** 930 genes are detected in under 1% of one species' structures. 138 of them are in
  the Gaussian top 10% and 6 in the NB top 10%.
- **The relabeling control cannot catch this.** The artefact follows species (depth, probe panels,
  absent genes), and the relabelings mix species.

The voom split-plot model works on specimen × bin pseudobulk and keeps its own mean–variance
weights. Its significant genes are not enriched for nearly-absent genes: 5.1%, against a base rate
of 8.4%; the Gaussian top 10% was 12.5%. It is therefore unchanged.

## Frozen specification

**Model.** Each gene gets an NB GLM:
- log link and offset log(library). The library is the total over the 15,567 orthologs measured in
  all specimens.
- The same designs as before: `contrast_designs`, `segment_designs` and the species trajectory model.
- The same prior weights, wᵢ = n / (2 · J_s · n_j).
- Code: `pseudospace/nb_gam.py`, `GENE_MODEL = 'nb.1'`.

**Dispersion.**
- One α per gene per coordinate run, estimated by weighted maximum likelihood under the
  species-split continuous M_full on all structures. There is no Cox–Reid adjustment, since
  n = 12,272 and p ≤ 22.
- The same α is reused for the nested models, both relabelings, the within-species comparisons,
  every refit (stability, warped coordinate, cortical removal, power draws), the step screen and the
  curves.
- Notebook 45's external pools have no coordinate. Each pool gets its own α under the saturated
  donor × segment model.
- No shrinkage. A DESeq2-style check shrank each α toward a trend α = 0.024 + 0.011/mean, with a
  normal prior on log α. On SCF13 it gave the same 31 joint-test calls (relabeled 0, 0; calibration
  0.068), and Spearman ρ = 0.999 between gene ranks. The rule fixed beforehand (both pass and
  Jaccard ≥ 0.8) keeps the unshrunk estimate.

**Statistic.**
- The likelihood-ratio statistic, LR = max(D_reduced − D_full, 0), replaces the partial F.
- The quasi-likelihood F is computed and reported, and used for nothing.
- Gene statistics remain **ranks**. Their residual degrees of freedom count structures, so they
  are not calibrated gene tests.

**Failures.**
- A gene is non-converged when a fit fails R's glm criterion, the LR is not finite, or α does not
  converge. Its statistics are replaced by the median of the converged genes.
- A gene with zero counts in one species is **separated**. Its fit lies on the boundary and its
  limiting deviances are valid: T_spatial ≈ 0, because there is no positional information in a
  species without counts, and T_level is large. Separated genes are flagged and never filled.

**VIF.** The joint test's VIF comes from the Pearson residuals, (y − μ)/√(μ + αμ²), of the species
M_full, through the unchanged `residual_pathway_correlations`. The Gaussian-residual VIF is reported
in notebook 65 as a sensitivity.

**Species-difference curves.**
- Δ(s) is a natural-log rate ratio, human over mouse.
- Its standard error is an HC3-type sandwich with an expected-information bread; the model-based
  form is also reported.
- Display and descriptor thresholds keep their values, now read on the natural-log rate-ratio scale.
  They are SM17's |Δ| ≥ 0.25 and notebook 58's amplitude, range and prominence cut-offs. For
  moderately expressed genes the two scales coincide.

**Primary method.**
- T_spatial (NB LR) with the joint test is the primary pathway method **by design** (user,
  2026-10-09).
- Notebook 61 compares the methods side by side and no longer selects one by rule. Before this
  decision, its "most calls" step could have picked the split-plot ORA (35 calls) over the pilot
  NB joint test (31).
- Notebook 37's E2 fallback ("keep the earlier list as primary" when fewer than 20 pathways are
  robust) is removed, because the earlier list comes from the superseded model.

**Stop rule.** If the NB joint test fails the relabeling rule or the calibration criterion on SCF13
in notebook 37, the cascade stops and the result is reported. Nothing is tuned.

## Order of work

1. Notebook 65 (new) reproduces the Gaussian results while they are still live, then the NB
   comparison, the gene-level diagnostics and the simulation.
2. Notebook 37 on SCF13, DPT13 and SCF13-ED.
3. Notebooks 39, 43 and 47. Notebook 47 becomes the tracked writer of `coordinate_robustness.csv`,
   its funnel and `d4_coordinate_runs.csv`, which previously came from scratch scripts.
4. Notebook 50 (primary list), then 45 (many-draw and external nulls).
5. Notebooks 56 and 57, then 61 and 63.
6. Notebooks 48, 38, 58–60 and 41; then the draft, the SI and the figure plan.

The Gaussian results, the working notebooks with their outputs and the scratch provenance scripts
are archived in `results/_archive_gaussian_2026-10-09/`, which is outside Git.
