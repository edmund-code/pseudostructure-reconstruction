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

## Notebook 87 results

`analysis/notebooks/87_pt_gene_model_audit.ipynb` (logic `87.gene_model_audit.1`) writes to
`results/pt_gene_model_audit/`. It followed the protocol above (0c8d28a), with the deviations listed at the end. The
tradeSeq run took 30.6 min (11,055 genes × 7 fits plus the subset). The cached rerun takes 100 s.

**Guards.**
- G1: 71's 11,055 genes, 11,383 units and 1,513 sets; 71's knots (0.240, 0.410, 0.649) and grid recomputed exactly.
- G-C (deviation 1): our code under 71's settings reproduces 71's T_spatial and LR_spec for ≥ 99.91% of genes in
  every partition at 1e-4. 23 gene statistics differ more; the condition z keeps Spearman 0.9991. The 71 calls are
  identical (83, with 1 and 0 relabeled).
  - The cause is α. Re-estimating it changes it for 310 near-Poisson genes on a flat likelihood (median α of those
    1e-8).
  - This is the run-to-run floor that every variant in C2 is read against.

### A · It is not tradeSeq, and tradeSeq is not calibrated on these units

**A1.** Table A1 of the protocol was read from both codes and stands. The two models share:
- the NB family, log link and log-library offset;
- a per-gene dispersion without shrinkage;
- cubic splines with about 6 basis functions and knots at position quantiles.

They differ in:
- **the penalty:** tradeSeq uses REML; ours has none;
- **the test:** tradeSeq uses a Wald test; ours uses LR;
- **conditionTest's target:** level + shape, against shape only in ours;
- **the weights;**
- **the unit of replication:** the unit in tradeSeq, the specimen in our F denominator.

**A2b and A2c · Calibration and calls (human vs mouse, the two balanced relabelings).**

| Method | Species genes (BH ≤ 0.05) | Relabeled genes | Relabeled share p ≤ 0.05 | Calibrated (≤ 0.075) |
|---|---|---|---|---|
| tradeSeq `conditionTest` (defaults) | 10,819 of 11,053 | 3,504, 5,339 | 0.476 | no |
| tradeSeq + specimen `U` | 10,761 of 11,000 | 1,923, 4,024 | 0.360 | no |
| ours: raw T_total LR (χ², 7 df) | 10,215 of 10,738 | 1,702, 1,717 | 0.274 | no |
| ours: raw T_spatial LR (χ², 6 df) | 5,144 | 405, 372 | 0.139 | no |
| **ours: positional F-test** | **3,526** | **0, 0** | **0.065** | **yes** |

- **Any test that treats the tubule unit as the replicate is miscalibrated here, tradeSeq's and our own raw LRs
  alike.** Under the relabelings, which put one human and one mouse specimen on each side, tradeSeq gives p ≤ 0.05
  for 48% of genes (36% with the specimens in `U`).
- Only the specimen-shape denominator (the positional F-test) brings the share to 0.065.
- This is the pseudoreplication problem (Squair 2021; Zimmerman 2021). It is why the positional F-test exists, and why
  it does not reduce to tradeSeq.

**A2a · Ranking** (species split, Spearman of −log10 p with tradeSeq's `conditionTest`).

| Ours | Spearman | Top-10% Jaccard |
|---|---|---|
| raw T_total LR (level + shape, its closest analogue) | 0.573 | 0.371 |
| raw T_spatial LR | 0.294 | 0.163 |
| positional F-test | 0.238 | 0.141 |

`conditionTest` mostly ranks level differences. Most genes differ in level between the species and the probe panels,
which is why it calls 98% of genes.

**A2d · Within mouse,** where both methods ask the same question (is the gene positional?), they agree:
- Spearman 0.821 between tradeSeq's `associationTest` and our spline-against-specimen LR, over 6,270 genes;
- top-10% Jaccard 0.676;
- BH ≤ 0.05: ours 4,178, tradeSeq 3,653.

**A2e · The penalty** (subset G, 1,516 genes).
- tradeSeq's REML leaves a median edf of 3.79 per condition smooth (IQR 2.28–4.84, maximum 6).
  - 22% of smooths keep ≥ 5.
  - edf rises with expression: from a median of 2.4 in the lowest decile to 5.0 in the highest.
- **Mouse curves agree:**
  - tradeSeq against 71: median r 0.88;
  - RMS difference 0.063 ln, against a median range of 0.52.
- **Human curves do not:**
  - median r 0.15;
  - RMS difference 0.26, against a range of 0.43.
  - The difference is the offset, not the penalty (PH2 below).

**Precedents** (protocol A3, PubMed-verified). Our construction is not new:
- an NB GLM on an unpenalised regression spline in pseudotime, tested by LR at a fixed dispersion, is the construction
  of Monocle 2's `differentialGeneTest` (`~sm.ns(Pseudotime, df=3)` against `~1`, VGAM NB, `lrtest`);
- tradeSeq is the penalised NB-GAM version with Wald tests;
- a sample-level denominator is the multi-sample principle of Lamian (Hou 2023).

"Defensible" therefore means: a published model class, plus a replicate-aware test that tradeSeq lacks, shown here to
be needed on these data. It does not mean "tradeSeq exactly".

### B · The tubule units fit the NB well (subset G; 1,465 usable genes)

| Check | Result | Misfit share (protocol rule) |
|---|---|---|
| B1 zeros: observed − NB expected (median) | −0.00003 (Poisson: +0.0011) | 0.3% |
| B2 Pearson dispersion (median; 1 under the model) | NB 1.015; Poisson 1.030 | — |
| B2 variance about the fitted mean ÷ NB variance, 20 bins of μ | 1.00–1.12 in 18 of 20 bins; 1.32 in one; 0.82 in the top bin (Poisson: 7.2 there) | — |
| B3 randomized quantile residuals (median KS D) | 0.007 | 0.0% |
| B4 library slope (median; 1 = offset correct) | 0.987 (IQR 0.94–1.02) | 3.4% |
| B5 α by position third | tests mostly null | 3.1% |
| B6 α by unit size tertile | falls with size: median slope of log α on log spots −0.31 | 1.5% |
| B7 α human ÷ mouse (median) | 1.14 | — |

- **Zero inflation:** none. The NB, and even the Poisson, predict the zeros.
- **Overdispersion is small.**
  - The median α is 0.037, and 10% of genes sit at the Poisson bound.
  - The Poisson fails only for the most expressed genes (7.2-fold excess variance in the top bin); the NB does not.
- **The log-library offset has the right slope.**
- **Dispersion is effectively constant** along the curve and across unit sizes. The per-third α ratios are large for
  weakly expressed genes (median max/min 8.0) because α is barely identified there, but the tests reject constancy for
  only 3.1% of genes.
- The direction of B6 (larger units, smaller α) is what aggregation of bins predicts, but the effect is weak.
- **B8:** α against the mean follows α = 0.0037 + 0.0001 / mean (Spearman −0.39). With 11,383 units the DESeq2
  shrinkage weight on each gene's own estimate is 1.0000, so shrinkage cannot matter here.

None of the protocol's "change the model" thresholds is reached (zeros or residuals > 10%; slope off by ≥ 0.2;
dispersion misfit > 25%).

### C · The parameters

**C1 · Fit by basis df** (species M_full; medians over G against df 6).

| df | ΔAIC | ΔBIC | Δ held-out deviance | Genes better than df 6 (held-out) |
|---|---|---|---|---|
| 3 | −3.2 | −47.3 | −0.16% | 92% |
| 4 | −2.7 | −32.1 | −0.11% | 93% |
| 5 | −1.7 | −16.4 | −0.05% | 87% |
| 6 | 0 | 0 | 0 | — |
| 7 | +2.1 | +16.8 | +0.06% | 11% |
| 8 | +4.2 | +33.5 | +0.12% | 3.5% |
| 10 | +7.1 | +65.8 | +0.24% | 0.9% |
| 12 | +10.6 | +98.7 | +0.35% | 0.5% |
| 6, equal knots | +0.25 | +0.25 | +0.00002% | 50% |

- Held-out, df 3 is best for 79% of genes; AIC picks it for 48% and BIC for 91%.
- **PH1:** this holds for the genes the positional F-test calls too (median best df 3; 8% best at ≥ 6).
- Gene curves along this curve are simple. The cost of df 6 is ≤ 0.16% of held-out deviance.
- Knot placement does not matter.

**C2 · Does the choice change 71's results?** Full universe, positional F-test and cameraPR. 71 has 83 calls with
relabeled 1 and 0, calibration share 0.072.

| Variant | z Spearman | Gene Jaccard | Calls (relabeled) | Calibration share | Call Jaccard | 71's rule |
|---|---|---|---|---|---|---|
| 71 as stored (run-to-run floor) | 0.999 | 0.997 | 83 (1, 0) | 0.071 | 1.00 | pass |
| df 4 | 0.930 | 0.778 | 99 (2, 1) | 0.089 | 0.80 | **fail** |
| df 5 | 0.962 | 0.868 | 87 (2, 1) | 0.069 | 0.93 | pass |
| df 8 | 0.947 | 0.843 | 88 (1, 0) | 0.052 | 0.84 | pass |
| df 10 | 0.909 | 0.787 | 95 (2, 3) | 0.043 | 0.80 | pass |
| equal knots | 0.980 | 0.910 | 86 (1, 0) | 0.064 | 0.92 | pass |
| unit weights (tradeSeq's) | 0.995 | 0.939 | 86 (1, 0) | 0.053 | 0.97 | pass |
| trended α | 0.998 | 0.981 | 80 (1, 1) | 0.080 | 0.96 | **fail** |
| shrunken α (DESeq2-style) | 0.999 | 0.994 | 83 (1, 0) | 0.072 | 1.00 | pass |
| **TMM offset (tradeSeq's)** | **0.446** | **0.369** | **166 (2, 0)** | 0.058 | **0.25** | pass |
| **poscounts offset (DESeq2's)** | **0.346** | **0.345** | **160 (5, 0)** | 0.077 | **0.15** | **fail** |

- **Basis df, knots, weights and dispersion move the calls by 3–20%.** This is the size of a change of df by one or
  two, and each choice keeps most of 71's list. df 5 to 8 all pass.
- **The normaliser is the one choice that matters.** TMM or poscounts size factors double the calls and keep only 15–25%
  of them.
- **PH2 shows why.** tradeSeq's TMM factor trends along the curve in human (Spearman with position 0.525; the mean
  log factor rises from −0.03 to +0.58–0.75 ln across position deciles) and not in mouse (0.037).
  - A factor that moves with position in one species enters every gene's species × position term.
  - It turns a normalisation difference into "shape" for most genes: 8,112 genes against 3,591.
  - It is also why tradeSeq's human curves disagree with ours (A2e).
- **Our library offset is supported at the unit level** (B4 slope 0.987). It is also the choice under which most genes
  are not called, as expected if most genes keep their shape.
  - The composition question is real: the human composition changes along the curve.
  - 68, 76 and 81 already handle it for programs, by centring on the median gene's Δ.

### Readout against the protocol's rules

- **tradeSeq:** its relabeled share is 0.476 (> 0.075), and ours is 0.065. **Our test stays,** and the calibration is
  the reason.
- **NB:** no threshold is reached. **The NB with a constant per-gene α and a log-library offset is adequate** for these
  units.
- **Parameters.** By the frozen rule, a variant is better when it improves held-out deviance for most genes and keeps
  71's rule. **df 5 meets it:**
  - 87% of genes improve, by a median 0.05%;
  - it passes the rule (calibration 0.069);
  - it keeps 93% of the calls.

  df 3 and 4 fit better still, but df 4 fails calibration and df 3 was not run downstream. The gain is negligible and
  the call overlap is high, so df 6 stays unless the user decides otherwise. That would be a new protocol.
  - Equal knots, unit weights and shrunken α change nothing material.
  - The normaliser is the open question. The current choice is supported (B4). It is recorded as an issue for the
    composition-aware analyses, not changed here.

### Recorded deviations

1. **G-C tolerance.**
   - The protocol asked rtol 1e-6.
   - Reproduction holds for ≥ 99.91% of genes at 1e-4. 23 gene statistics differ more, because α is re-estimated for
     310 near-Poisson genes.
   - The run stops when > 0.5% of genes differ. The floor is reported as "71 as stored".
2. **tradeSeq p-values.** tradeSeq reports 1 − pchisq(W, df), which is exactly 0 below about 1e-16. The upper tail of
   the same W is recomputed with `chi2.sf`. Shares and calls at 0.05 are unchanged; rankings are no longer tied.
3. **B2's display.** The protocol's "variance against mean in 20 bins of μ" mixes genes within a bin, which inflates the
   top bin. The figure and table show the mean squared deviation from each unit's own fitted mean, divided by the NB
   (and Poisson) expectation.
4. **The R run.** The first run lost a forked worker, probably to memory. The script now fits genes in chunks of 1,500
   and retries a failed chunk once with half the workers. On synthetic data the chunked and whole-matrix outputs are
   identical.
5. **Post hoc PH1 and PH2** (above), added after the first complete run.

### Outputs

- `tables/`:
  - `a1_side_by_side.csv`;
  - `a2a_ranking.csv`, `a2bc_calibration_and_calls.csv`, `a2d_within_mouse.csv`, `a2e_penalty_and_curves.csv`;
  - `b_diagnostics_subset.csv`, `b_summary.csv`, `b2_mean_variance_pooled.csv`;
  - `c1_fit_by_df.csv`, `c1_best_df_shares.csv`, `c1_best_df_by_decile.csv`;
  - `c2_variants.csv`;
  - `guard_c_reproduction.csv`;
  - `ph1_best_df_by_call.csv`, `ph1_heldout_change_by_call.csv`, `ph2_tmm_factor_along_curve.csv`.
- `tradeseq/`: tradeSeq's per-gene tests, the subset's edf, θ, AIC and curves, the TMM offsets, the session.
- `figures/`: `fig87a_tradeseq`, `fig87b_nb_diagnostics`, `fig87c_parameters` (png and pdf).
