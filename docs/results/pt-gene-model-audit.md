# Gene model audit: tradeSeq, the NB assumptions and the parameters (notebook 87)

Frozen on 2026-10-13, before any fit of notebook 87 on our data. The user asked: "can you check our gams follow
tradeseq exactly? So it's defensible our gams are peer reviewed accepted already. Check our tubules fit the nb
assumptions. Check if parameters at the best. Show me the figures".

**What I had seen before freezing.**
- Our gene model's code (`pseudospace/nb_gam.py`, `pathway_pipelines.standard_gene_model`, `specimen_shape_lr`,
  `moderated_f`, `pathway_calibration.contrast_designs`) and its frozen specification (`pt-gene-model-nb.md`).
- tradeSeq 1.20.0's source, read in an isolated conda environment (`tradeseq-audit`: R 4.4.3, tradeSeq 1.20.0,
  mgcv 1.9.4, edgeR 4.4.0; nothing was installed into `kidney-pseudospace`). One fit on simulated counts showed the
  fitted model's settings (method REML, family `Negative Binomial(θ)`, 6 unconstrained coefficients per condition
  smooth). A timing run on simulated counts of our size: about 0.4 s per gene.
- 71's run record (compute times) and its stored tables. No new fit on our data.

## Labels

- Notebook: `analysis/notebooks/87_pt_gene_model_audit.ipynb`
- Logic: `87.gene_model_audit.1`
- Module: `pseudospace/gene_model_audit.py`, tested in `tests/test_gene_model_audit.py`
- R script: `analysis/scripts/tradeseq_audit.R`, run with `Rscript` from the `tradeseq-audit` environment
  (`analysis/scripts/tradeseq_audit_environment.yml`)
- Results: `results/pt_gene_model_audit/` (`tables/`, `figures/`, `tradeseq/`, `stage_cache/`, `run_record.json`)

## Inputs

- `hvm_inputs` exactly as 71 (mt-free library, support-scaled elastic-curve position, 11,383 units, 71's universe).
  Guard: the universe, units and gene sets equal 71's.
- 71's stored species and relabeling models, denominators and gene statistics (read only).
- AKI is out of scope: the same code runs both datasets, and 71's AKI screen already fails calibration.

## A · Our model against tradeSeq

### A1 · The modelling choices, side by side (from the code)

| Choice | tradeSeq 1.20.0 (`fitGAM`, defaults) | Ours (71: `standard_gene_model`) |
|---|---|---|
| Count model | NB, log link, `mgcv::gam(family = "nb")` | NB, log link, own IRLS (`nb_gam`, R glm's convergence rule) |
| Dispersion | one θ per gene, estimated inside mgcv's REML outer iteration, no shrinkage | one α per gene, ML under the species M_full, then fixed for every nested model; no shrinkage |
| Smoother | `s(t, bs = "cr", k = 6)`: cubic regression spline, 6 knots at pseudotime quantiles (0, 0.2, …, 1) | cubic B-spline, df 6: 3 interior knots at the position quartiles, boundaries 0 and 1 |
| Penalty | second-derivative penalty; smoothing parameter by REML, shared across conditions (`id = 1`) | none (fixed df) |
| Conditions | one unconstrained smooth per condition (level and shape) | group main effect + group × basis |
| Offset | log(library × TMM factor) (`edgeR::calcNormFactors`) unless given | log(library over the measured orthologs, mt-free) |
| Weights | none (cell weights only assign lineages) | each group half the weight, each specimen an equal share |
| Replicates (specimens) | none by default; `U` can add parametric covariates | specimen intercepts within group (sum-to-zero) |
| Condition test | `conditionTest`: Wald χ² on the differences between the conditions' smooth coefficients at the knots (level + shape), Bayesian covariance `Vp`, eigenvalue-thresholded | LR at fixed α of shape only (T_spatial, 6 df), divided by the specimen-shape LR (the positional F-test, limma-moderated) |
| Association test | `associationTest`: Wald χ² that the smooth is constant (differences from the start at 2 × nknots points) | within species: LR of the spline against specimen intercepts only |
| Unit of replication | the cell (here the tubule unit) | the specimen (the F denominator) |

**Read:** our model is not tradeSeq. Both are NB regression splines along a pseudotime coordinate with a log
library offset, a per-gene dispersion without shrinkage and about 6 basis dimensions. They differ in the penalty,
the test (Wald against LR), the condition test's target (level + shape against shape only), the weights and, most,
the unit of replication. A1 is printed in the notebook from the code, and section A2 measures what the differences
do.

### A2 · tradeSeq run on our units

- **Units and genes.** 71's 11,383 units and its full universe. Counts, position, specimen, species and library
  are exported once (`tradeseq/`) and read by `tradeseq_audit.R`.
- **Fits** (`fitGAM`, nknots 6, single lineage, `cellWeights` 1, one BLAS thread per worker):
  - **TS-default:** conditions = species; offset = tradeSeq's default (TMM over the exported universe genes).
  - **TS-relabel 1, 2:** the same with conditions = 71's two balanced relabelings.
  - **TS-U:** conditions = species with `U` = intercept + within-species specimen contrasts (tradeSeq's way of adding
    replicates); and its two relabelings.
  - **TS-mouse:** the two mice only, no conditions, offset as TS-default (recomputed on mouse units).
- **Tests:** `conditionTest` (global, l2fc 0) for every two-condition fit; `associationTest` (global) for TS-mouse.
- **Curves and smoothness** (gene subset G below): the same TS-default model refitted with `sce = FALSE`, to read
  each condition smooth's effective degrees of freedom, θ, AIC and the two condition curves on 71's grid.

**Comparisons.**
- **A2a · Ranking.** Spearman of −log10 p between TS-default's `conditionTest` and our statistics: the raw T_total LR
  (level + shape, χ² with 7 df; the closest analogue of `conditionTest`), the raw T_spatial LR (χ², 6 df) and the
  positional F-test (71's z). Overlap of the top 10% genes.
- **A2b · Calibration (primary).** The share of genes with p ≤ 0.05 under each relabeling, for TS-default, TS-U, our
  raw LRs and the positional F-test. 71's gene-level rule: ≤ 0.075 is calibrated. The relabelings group one human
  and one mouse specimen per side, so species cancels and only between-specimen differences remain.
- **A2c · Calls.** Genes with BH q ≤ 0.05 under the species split, and the ratio of mean relabeled calls to species
  calls (the NCDR) for each method.
- **A2d · Within mouse.** Spearman between TS-mouse's `associationTest` and our within-mouse LR (the spline against
  specimen intercepts, α under the spline model), and the top-10% overlap. Descriptive.
- **A2e · Penalty.** The edf of each condition smooth (TS-default, gene subset): the share of genes with edf ≥ 4 of
  the 5 possible beyond the level. Curve agreement: Pearson r and RMS difference between TS-default's species
  curves and 71's group curves on 71's grid (log scale, both centred).

### A3 · Precedents (PubMed-verified)

- tradeSeq: Van den Berge et al. 2020 ([PMID 32139671](https://pubmed.ncbi.nlm.nih.gov/32139671/)), penalised NB-GAM
  and Wald tests.
- Monocle: Trapnell et al. 2014 ([PMID 24658644](https://pubmed.ncbi.nlm.nih.gov/24658644/)); Qiu et al. 2017
  ([PMID 28825705](https://pubmed.ncbi.nlm.nih.gov/28825705/), [PMID 28114287](https://pubmed.ncbi.nlm.nih.gov/28114287/)).
  `differentialGeneTest`'s defaults are `~sm.ns(Pseudotime, df=3)` against `~1`, a VGAM NB at a fixed size, and
  `VGAM::lrtest`: an unpenalised regression spline with a likelihood-ratio test, the same construction as ours.
- PseudotimeDE: Song & Li 2021 ([PMID 33926517](https://pubmed.ncbi.nlm.nih.gov/33926517/)) shows that tests which
  ignore pseudotime uncertainty or treat cells as independent can be miscalibrated.
- Lamian: Hou et al. 2023 ([PMID 37949861](https://pubmed.ncbi.nlm.nih.gov/37949861/)) models sample-level
  variation for multi-sample pseudotime comparisons, which is the role of our specimen-shape denominator.
- Pseudoreplication: Squair et al. 2021 ([PMID 34584091](https://pubmed.ncbi.nlm.nih.gov/34584091/)); Zimmerman et
  al. 2021 ([PMID 33531494](https://pubmed.ncbi.nlm.nih.gov/33531494/)).
- NB for UMI counts: Svensson 2020 ([PMID 31937974](https://pubmed.ncbi.nlm.nih.gov/31937974/)); Hafemeister &
  Satija 2019 ([PMID 31870423](https://pubmed.ncbi.nlm.nih.gov/31870423/)); Lause et al. 2021
  ([PMID 34488842](https://pubmed.ncbi.nlm.nih.gov/34488842/)).

## B · Do the tubule units fit the NB assumptions?

Our units are polygon aggregates of many Visium HD bins, so their library varies widely. All checks use 71's species
M_full (the condition model's design, weights, knots and α), refitted for the gene subset G, at unit level, without
weights in the residuals.

- **Gene subset G.** From 71's universe: 150 genes per decile of log mean count (seed 87), plus the markers in the
  universe among Slc34a1, Slc5a2, Slc5a12, Slc22a6, Slc22a8, Slc22a7, Slc7a13, Lrp2, Cubn, Gpx3, Fxyd2, Aqp1,
  Slc13a3, Cyp4b1, Hnf4a, Spp1, Umod.
- **B1 · Zeros.** Per gene: observed zero fraction against the NB expectation, mean of (1 + αμ)^(−1/α), and the
  Poisson expectation, mean of e^(−μ). Excess zeros: z = (observed − expected zeros) / √Σp₀(1 − p₀), one-sided, BH
  over G. **Misfit:** q < 0.05 and observed − expected ≥ 0.02.
- **B2 · Variance.** The Pearson dispersion Σ(y − μ)² / (μ + αμ²) / (n − p) per gene (1 under NB), and the same with
  the Poisson variance. Pooled over G: the empirical variance against the mean in 20 bins of fitted μ, with the
  Poisson line and the NB curve at the median α.
- **B3 · Randomized quantile residuals** (Dunn & Smyth 1996): r = Φ⁻¹(u), u uniform on (F(y − 1), F(y)) under the
  fitted NB (seed 87). Per gene, the Kolmogorov–Smirnov D against N(0, 1), with its p (BH over G). **Misfit:** q <
  0.05 and D ≥ 0.05. Pooled QQ plot over G.
- **B4 · Offset slope.** The M_full with log(library) added as a covariate (the offset kept): the slope is 1 + its
  coefficient. LR test (1 df) at the gene's α, BH over G. **Misfit:** q < 0.05 and |slope − 1| ≥ 0.2.
- **B5 · Dispersion along the curve.** Given μ, one ML α per position third (equal-count thirds per specimen) against
  one α: LR (2 df), BH over G. **Misfit:** q < 0.05 and the largest pairwise ratio of the thirds' α ≥ 2.
- **B6 · Dispersion and unit size.** Given μ, one ML α per tertile of unit size (`n_spots`, within specimen). If units
  are sums of k bins each with dispersion α_b, the unit's α is about α_b / k, so α should fall with size. The slope of
  log α on log median `n_spots` across tertiles, per gene; LR (2 df) as in B5. Reported with the median slope.
- **B7 · Dispersion by specimen and species.** Given μ, one ML α per specimen; the human / mouse ratio of α.
- **B8 · Dispersion against mean.** α against the mean count, with a trend (B·P3).

## C · Are the parameters the best?

### C1 · Model fit by the number of basis functions (gene subset G)

- basis_df ∈ {3, 4, 5, 6, 7, 8, 10, 12}, knots at quantiles (df 3 is a cubic polynomial).
- **AIC and BIC** of the species M_full fitted without weights (α re-estimated by ML per df), as tradeSeq's
  `evaluateK` uses AIC. Per gene, the df with the smallest value; the median ΔAIC and ΔBIC against df 6.
- **Held-out deviance.** 5-fold cross-validation over units, folds drawn within specimen (seed 87). Per fold, α and β
  are fitted on the training units (71's weights), and the NB deviance of the held-out units is summed. Per gene, the
  df with the smallest total; the median relative change against df 6.
- Knot placement at df 6: quantile (ours, tradeSeq) against equally spaced interior knots (0.25, 0.5, 0.75), by
  held-out deviance.

### C2 · Do the choices change 71's results? (full universe)

Each variant reruns 71's positional F-test on all genes and partitions with one change: the T_spatial and LR_spec
LRs (the module's `positional_lrs`, which reproduces 71's statistics with 71's settings: guard G-C, rtol 1e-6),
`moderated_f` (covariate log mean count) and `camera_pr` on z with 71's gene sets.

| Variant | Change |
|---|---|
| df 4, 5, 8, 10 | basis_df |
| equal knots | equally spaced interior knots at df 6 |
| unit weights | no specimen-balancing weights (tradeSeq has none) |
| trended α | α from a fitted trend α(mean) = a + b / mean (least squares on log α) |
| shrunken α | DESeq2-style: log α shrunk toward the trend with a normal prior; the sampling variance of log α is trigamma((n − p) / 2) |
| TMM offset | log(library × TMM factor), tradeSeq's default |
| poscounts offset | DESeq2's poscounts size factors over the universe genes |

**Read per variant:** Spearman of the condition z with 71's; Jaccard of genes with q ≤ 0.05; the pathway calls,
their Jaccard with 71's, the relabeled calls, the NCDR and the calibration share (71's pass rule: NCDR < 0.25, each
relabeling < 5% of sets, calibration ≤ 0.075); the gene-level calibration share.

### C3 · Penalised or not

A2e's edf and curve agreement answer it: if REML leaves the smooths at about our df, the penalty does not change the
fit.

## What would make us change the model

- **tradeSeq.** If TS-default's relabeled share is ≤ 0.075 and it reproduces our ranking (Spearman ≥ 0.8 against
  the positional F-test), tradeSeq could replace our test outright. If its relabeled share is above 0.075 and ours is
  not, our test stays and the difference is the reason.
- **NB.** If more than 10% of G misfits in B1 (zeros) or B3 (residuals), the NB is not adequate for these units. If
  B4's median slope differs from 1 by ≥ 0.2, the offset is wrong. If B5 or B6 misfit more than 25% of G, a constant
  dispersion is wrong; the fix would be an α that depends on position or unit size.
- **Parameters.** A variant is better when it improves held-out deviance for most genes **and** keeps 71's pass rule.
  A variant that only moves calls (Jaccard < 0.8) with no gain in fit is reported, not adopted. Any change to the
  primary model would be a new protocol decided by the user.

## Bounds

- 71's gene statistics are ranks. Their residual df counts units, not animals (`pt-gene-model-nb.md`). Every NB
  diagnostic is at unit level.
- Cohort: two mice and one human donor (two sections). The calibration checks use the only two balanced relabelings.
- Everything here is descriptive apart from the stated per-gene tests.

## Figures (plots only)

- 87a: tradeSeq against ours: −log10 p scatter, relabeled shares, edf histogram, curve agreement examples.
- 87b: NB diagnostics: zeros observed against expected, variance against mean, pooled randomized-quantile QQ,
  offset slopes, α by position third and by unit size.
- 87c: parameters: ΔAIC, ΔBIC and held-out deviance by df; the variants' agreement with 71 and their calibration.
