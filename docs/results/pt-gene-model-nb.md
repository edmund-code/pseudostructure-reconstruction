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
- **Real data.** 930 genes are detected in under 1% of one species' structures. 139 of them are in
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

## Notebook 65 results

`analysis/notebooks/65_pt_gene_model_comparison.ipynb` (logic version `65.gene_model.1`) writes to
`results/pt_gene_model_comparison/`. It recomputes the Gaussian side from code.
- On its first run every Gaussian statistic, VIF and call matched the live notebooks 37, 45 and 50.
  Notebook 45's pool was rebuilt and matched too.
- The Gaussian reference is saved in `reference/`, and later runs must reproduce it.
- The live-notebook guards switch off once notebook 37's gene table carries `gene_model = nb.1`, so
  notebook 37's NB re-run must write that column.

**Pathway screens on SCF13** (species calls; relabeled calls; share of relabeled p ≤ 0.05, limit 0.075):

| Screen | Gaussian | NB, Gaussian VIF | NB, NB Pearson VIF (frozen) |
|---|---|---|---|
| A6·JNT | 60; 2, 0; 0.060 | 31; 0, 0; 0.067 | 74; 13, 1; **0.085, fails calibration** |
| A4·JNT | 50; 0, 0; 0.064 | 15; 1, 1; 0.079, fails | 62; 9, 9; 0.105, fails |
| Continuous adds / steps add | 13 / 3 | 16 / 0 | 21 / 9 |

- **The pilot is reproduced exactly.** It used the Gaussian VIF, so the pilot table above
  describes the middle column.
- **The frozen primary fails criterion 2.** The NB LR with the Pearson-residual VIF, computed here
  for the first time, passes the relabeling rule (NCDR 0.09) but not the calibration limit.
- **Why.** The median VIF falls from 1.30 to 1.06, and the median within-pathway residual
  correlation from 0.010 to 0.002.
  - Gaussian log1p residuals share a depth component: within each specimen, a structure's mean
    residual correlates with its log library at r = 0.75–0.94.
  - The NB offset absorbs depth, giving r = −0.05 to −0.29.
  - Much of the Gaussian residual correlation was therefore shared depth. It happened to deflate the
    joint test enough to pass the relabelings.
- **Stop rule.** Notebook 37 would compute the same screen. Under the stop rule above, the cascade
  would therefore stop there with this specification. Nothing was tuned, and the choice is open.

**Gene level, as in the pilot:**
- T_spatial Spearman between the models is 0.36.
- The Gaussian-only top-10% genes have median |log ratio| 2.11 (NB-only: 0.83). 25% of them are
  nearly absent in one species (NB-only: 0.2%).
- Of the 930 nearly absent genes, 139 are in the Gaussian top 10% (12.5%) and 6 in the NB top 10%.
- The 4,588 split-plot genes are 5.1% nearly absent, against the 8.4% base rate.
- In the simulation, pure offsets reach the top-10% cut at 10.1% under Gaussian (27% in the
  largest-offset tertile) and at 0% under NB.

**Shrinkage.**
- With the Gaussian VIF it reproduces the pilot check: the same 31 calls (0, 0; 0.068; Jaccard 1.0),
  so the unshrunk α stays primary.
- With the NB Pearson VIF both versions fail calibration (0.085 and 0.087), so the rule cannot be met.

**What size the VIF should be.**
- Notebook 37 validated the Gaussian VIF by fitting E[z²] = a(1 + (m − 1)ρ) to the matched z of all
  1,513 pathways under the two relabelings, using `pathway_decoys.decoy_inflation`.
- The same fit applied to the NB statistic gives a = 1.04, ρ = 0.0096, and an implied median VIF of
  1.35. For the Gaussian statistic it gives 1.03, 0.0086 and 1.30.
- The NB statistic's own null correlation therefore matches the Gaussian-residual estimate (median
  ρ 0.0096). The Pearson-residual estimate (0.0018) is about 5 times too small.
- This uses only the relabelings, never the species split. The choice of VIF source is the user's
  (2026-10-09).

## Notebook 66: packaged set tests and clustering-first on the NB statistic

`analysis/notebooks/66_pt_pathway_pipelines_compared.ipynb` (logic version `66.pipelines.1`) writes to
`results/pt_pathway_pipelines_compared/`. It runs every pathway method on one shared NB gene model
(`pathway_pipelines.standard_gene_model`), for human vs mouse on SCF13 and for AKI vs control mouse
(notebook 45's inputs). Each row gives condition calls; relabeled calls; share of relabeled p ≤ 0.05:

| Method | Human vs mouse | AKI vs control |
|---|---|---|
| limma cameraPR on rank-normal T_spatial | 98; 46, 29; 0.121 (fails both) | 144; 87, 86; 0.170 (fails both) |
| gseapy prerank GSEA on rank-normal T_spatial | 146; 75, 80; 0.147 (fails both) | 203; 120, 146; 0.212 (fails both) |
| NB joint test, VIF fitted on the relabelings | 32; 3, 0; 0.067 | 12; 0, 0; 0.072 |
| Clustering-first ORA (Hallinan's package, single-split null) | 162 and 159; 0, 0 | 469 and 482; 2, 0 |

- **The packaged rank tests are not specific here.** Their relabeled z-scores are inflated by a
  constant factor (1.51 and 2.13), so no setting of `inter.gene.cor` repairs them. Ranking within
  expression × detection × coverage strata (post hoc) nearly passes human vs mouse (share 0.082) but
  not AKI (0.135).
- **The joint test's VIF is fitted on the relabelings that judge it.** Fitted on one and judged on
  the other, its share is 0.067 and 0.075; the AKI value sits at the limit.
- **Clustering-first is specific under a fair null, but returns no p-value per pathway.** It calls
  pathways whose genes share one direction, including constant offsets (oxidative phosphorylation)
  that T_spatial removes by design. It also finds AKI proliferation (G2-M, M phase), which T_spatial
  ranks low. Overlap with cameraPR is Jaccard 0.15 over redundancy clusters on both datasets.

## Notebook 67 protocol: T_spatial normalized by specimen shape variation

Frozen on 2026-10-09, before any fit. The user asked whether `T_spatial` should be normalized.
Notebook 66 suggested two reasons why packaged set tests on `T_spatial` are not specific:
- it grows with counts;
- it treats structures as independent, so shape differences between individual specimens count as
  group differences.

**Statistic.** The NB fits stay the same: the shared gene model, the species-split α, and the
weights and offset of `contrast_designs`.
- **Numerator.** `T_spatial`, unchanged: LR of M_full against M_full without group × basis, with
  df₁ = 6.
- **Denominator.** LR_spec = D(M_full) − D(M_spec), where M_spec adds a specimen-specific shape to
  M_full: specimen-within-group × basis. For the condition split df₂ = 12; for a relabeling, which
  keeps the condition shape as nuisance, df₂ = 6.
- s² = LR_spec / df₂ is moderated by `limma::squeezeVar`, with log mean count per structure as
  covariate. This is the primary. Sensitivities: no covariate, and `robust = TRUE`.
- The moderated F is (T_spatial / df₁) / s²_post, with df (6, df₂ + d₀). It converts to
  z = Φ⁻¹(1 − p).

**Pathway tests on z**, with gene sets as in notebook 66 (sizes 10–300):
- `limma::cameraPR`, with `inter.gene.cor` 0.01 and `use.ranks` FALSE. A call is Up with FDR ≤ 0.05.
  This is the primary test.
- gseapy prerank on z.
- ORA with `gseapy.enrich` on genes with BH ≤ 0.05; the background is all tested genes.
- Sensitivity: `cameraPR` on rank-normal scores of the F.

**Datasets.** Human vs mouse PT on SCF13 and AKI vs control mouse PT, each with the condition split
and both balanced relabelings, exactly as in notebook 66.

**Pass rule.** Notebook 66's rule, unchanged:
- **Relabeling.** NCDR < 0.25, and each relabeling calls < 5% of tested pathways.
- **Calibration.** The pooled share of relabeled one-sided p ≤ 0.05 is ≤ 0.075.

**Gene-level checks**, reported for information and not used for selection:
- under the relabelings, the share of genes with p ≤ 0.01 and p ≤ 0.05, and the median-based λ;
- the Spearman correlation of relabeled z with mean expression.

Genes with zero counts in any specimen are flagged. They stay in the primary and are dropped in a
sensitivity.

**Decision.** If the primary passes on both datasets, it becomes the user's candidate replacement
for the joint test as primary; the user decides. Otherwise the result is reported as is. Nothing is
tuned after the fits, and any later check is labelled post hoc.
