# Pseudostructure reconstruction: literature survey, problem definition and method hypotheses

Status: design document, written 2026-10-03 **before any real-data outputs** for the proposed
method. It extends the [research ledger](pseudostructure-research-ledger.md) (notebooks 15–19) and
should be read with it. Nothing here is a result.

Paper goal: describe a pseudostructure reconstruction method suited to ordering segmented tubule
cross-sections along the nephron, compare it with Scanpy DPT (notebooks 02/03) and nonbranching
scFates (notebook 13), then repeat the control-vs-AKI and human-vs-mouse program comparisons
(notebooks 05/06/07/12) on the new coordinate. The contribution should be modest and fitted to the
problem: a problem-specific estimator plus a validation protocol that works without ground-truth
position. It is not new transport mathematics.

## 1. What makes our problem different

Generic trajectory inference (TI) assumes that cells are snapshots of a dynamic process and that
sampling density reflects dwell time. None of that holds here. Our observations have these
properties:

| # | Property | Consequence for the estimator |
|---|---|---|
| P1 | The unit is a **segmented tubule cross-section**: an aggregate of 2 µm Visium HD bins inside an H&E polygon. | Library size scales with cross-section area and capture. Exposure must be modelled, not normalised away (ledger: NB16 coverage ordering, NB17 thinning). |
| P2 | Many **physical realizations** (nephrons) of one ordered object, with unknown nephron identity. | The target is a shared canonical coordinate. Profile count is not replication; specimen is. |
| P3 | **Known topology and orientation**: a single nonbranching line, glomerulus → S1 → S2 → S3, with coarse ordinal labels from upstream segmentation. | No topology search is needed. Branching TI machinery (PAGA, trees) is unnecessary freedom. |
| P4 | **Non-uniform sampling along the axis**: convoluted S1/S2 coils are cut many times, while straight S3 runs through the outer stripe. Slice planes differ by specimen; the human slices are cortex only. | Occupancy is a nuisance, not information about position. Density-driven coordinates (DPT, principal curves) can be warped by it, and partial coverage must not be stretched to [0, 1]. |
| P5 | **Physical tissue context** is measured: polygon centroids, adjacency between neighbouring tubules, glomerular centroids, morphology. | This is the only information independent of expression. Depth is a poor longitudinal proxy (winding nephrons), but *local adjacency* is a different and weaker assumption (§4, H3). |
| P6 | **Domains** (healthy mouse, AKI mouse, healthy human) change expression at a conserved anatomy. | Injury and species programs can move sections along an expression-only axis. The notebook-05 note that `Krt20`/`Krt8` shape effects are "a property of the DPT coordinate itself" is a symptom. Position and domain deviation must be separated by restriction, not by Harmony erasure. |
| P7 | Two specimens per mouse condition; one human donor. | Validation must use withheld information (genes, specimens, reads, sections), not significance. |

## 2. Literature map: what each family assumes and what we can borrow

### 2.1 Graph/diffusion trajectory inference (the current baseline)

- **DPT** (Haghverdi et al. 2016, Nat Methods): random-walk distance on a diffusion map from a root.
  It is robust to noise and outliers, but the diffusion operator depends on sampling density unless
  the density is explicitly normalised. Density-corrected diffusion variants exist (e.g. target-measure
  diffusion maps, [Banisch et al.](https://arxiv.org/pdf/1710.03484)). Pseudotimes from random
  subsets of cells can vary substantially
  ([Campbell & Yau 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC5117567)).
- **Benchmark**: [Saelens et al. 2019](https://www.biorxiv.org/content/10.1101/276907v1) compared
  29 methods. Slingshot, TSCAN and Monocle DDRTree did well on simple topologies, PAGA on complex ones.
  Performance depended on topology, which our problem fixes (P3).
- *Mismatch with our problem:* density sensitivity (P4), no exposure model (P1), no domain
  restriction (P6).
- *Borrow:* DPT as baseline and as one possible initialiser. A density-invariance stress test
  (§5, E2) that DPT should fail.

### 2.2 Principal curves and graphs (scFates)

- [ElPiGraph](https://arxiv.org/pdf/1804.07580) elastic principal graphs, wrapped by
  [scFates](https://scfates.readthedocs.io/en/latest/Basic_Curved_trajectory_analysis.html)
  `tl.curve`. It fits a smooth curve through the point cloud and projects points onto it.
- *Mismatch:* the curve follows the dominant variance and density. It has no likelihood, exposure
  term or anchors. A reproducible nuisance program is fitted as readily as anatomy (our NB19 PC1
  failure is the same mechanism).
- *Borrow:* the "curve in expression space plus projection" structure is what our mean atlas already is.

### 2.3 Probabilistic latent pseudotime

- GPLVM-based pseudotime and capture-time priors ([GrandPrix](https://pmc.ncbi.nlm.nih.gov/articles/PMC6298059)).
- [Ouija](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6298060/): marker-gene switch/transient
  parametric curves with Bayesian latent position.
- [PhenoPath](https://www.nature.com/articles/s41467-018-04696-6): a shared latent axis plus
  covariate-by-position interactions.
- *Borrow:* an explicit generative likelihood, informative priors from known markers (Ouija),
  and covariate interaction terms for domains (PhenoPath). These show that latent position plus
  covariates is not itself novel (as the ledger already states).

### 2.4 Landmark-based zonation reconstruction (the closest biological analogue)

- **Liver lobule**: [Halpern et al. 2017](https://pubmed.ncbi.nlm.nih.gov/28166538/) inferred
  lobule layer from smFISH-calibrated landmark genes. About 50% of genes were zonated, many non-monotone.
- **Intestinal villus**: [Moor et al. 2018](https://www.biorxiv.org/content/10.1101/261529v1)
  calibrated landmark genes by laser-capture microdissection of five villus compartments.
- These are the closest problem class: repeated physical units, unknown unit identity, and a 1D
  axis inferred from expression using external anatomical calibration.
- *Borrow:* **external calibration**. For PT, the public microdissected segment transcriptomes of
  [Chen et al. 2021 (JASN)](https://pubmed.ncbi.nlm.nih.gov/33769951/)
  ([data](https://esbl.nhlbi.nih.gov/MRECA/Nephron)) give independent S1/S2/S3 reference profiles.
  These are an ordered external anchor that we have not yet used. Scope check:
  [Park et al. 2018](https://pmc.ncbi.nlm.nih.gov/articles/6188645) found that microdissected
  PT-S2 is a transcriptomic mixture (about 18% S1-like, 58% S2-like, 19% S3-like), which supports
  a continuum rather than three states.

### 2.5 Optimal transport (OT) families

| Family | Representative | Assumption | Relevance |
|---|---|---|---|
| Static OT between time points | [Waddington-OT](https://pmc.ncbi.nlm.nih.gov/articles/PMC6402800), [moscot](https://pubmed.ncbi.nlm.nih.gov/39843746/) | Snapshots at *known* times, with mass growth | We have no time labels; the "time" is the unknown. Not directly applicable. |
| Gromov–Wasserstein (GW) / fused GW (FGW) | [SCOT](https://www.biorxiv.org/content/10.1101/2020.04.28.066787v1), [PASTE](https://www.biorxiv.org/content/10.1101/2021.03.16.435604.full.pdf), [PASTE2 (partial FGW)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9881963/) | Matches two spaces through their *internal* geometry; FGW adds a feature cost. PASTE fuses expression and spatial distance. | GW to a 1D target is a soft seriation. FGW lets a feature cost and a structural cost share one coupling. |
| Reconstruction onto a target geometry | [novoSpaRc](https://www.nature.com/articles/s41586-019-1773-3) | Cells mapped to a known target shape (including a 1D villus axis) by GW plus optional landmarks; balanced marginals | The **direct template for our formulation**. Its balanced marginal forces uniform occupancy, which is wrong for P4. |
| Unbalanced / partial OT | [FUGW](https://arxiv.org/pdf/2206.09398) (aligning individual brains), PASTE2, unbalanced Monge maps | KL-relaxed marginals absorb size differences and missing regions | Handles P4 (occupancy) and human cortex-only partial coverage. |
| Dynamic OT / flows | [TrajectoryNet](https://arxiv.org/abs/2002.04461v1), MIOFlow | Continuous paths between known time snapshots | Our NB15 probability-flow test did not earn its complexity. Deprioritised. |
| Differentiable sorting | [Cuturi, Teboul & Vert 2019](https://research.google/pubs/differentiable-ranking-and-sorting-using-optimal-transport/) | Ranking as entropic OT to a sorted template | Shows that "assign N sections to K ordered grid points by entropic OT" is a soft sort. It gives a clean computational form. |
| Spatial OT trajectory | [SpaTrack](https://www.biorxiv.org/content/10.1101/2023.09.04.556175.full.pdf) | Transition cost mixes expression and physical distance | Precedent for a tissue-distance term; still a generic TI. |

Seriation theory ([spectral seriation and latent orderings](https://arxiv.org/pdf/1807.07122),
[line-embeddings of graphons](https://arxiv.org/pdf/2007.06444)): recovering an order from noisy
similarities is a 1D manifold-learning problem with known recovery conditions. This is useful for
stating identifiability honestly. A nuisance axis with a larger signal than position will be
recovered instead (our NB19 failure).

### 2.6 Spatially aware trajectory and gradient methods

- [stLearn PSTS](https://www.biorxiv.org/content/10.1101/2020.05.31.125658.full.pdf) and
  [SpaceFlow](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9283532/) mix spatial and expression
  distances in pseudotime.
- GASTON ([Chitra et al. 2025](https://www.biorxiv.org/content/10.1101/2023.10.10.561757.full.pdf))
  learns a tissue "isodepth" and piecewise-linear gene functions of it. GASTON-Mix gives one
  isodepth per domain.
- *Mismatch:* these model *tissue-plane* coordinates. Our coordinate runs along a winding tube, so
  tissue position is not the latent variable (the user's depth caveat). Only *local adjacency* is
  plausibly informative.

### 2.7 Comparing conditions along a trajectory

- [Lamian](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10638410/) accounts for multi-sample
  variability; condiments and tradeSeq do not.
- [Genes2Genes](https://www.biorxiv.org/content/10.1101/2023.03.08.531713v1.full.pdf) aligns
  trajectories gene by gene with match and mismatch states (already used in notebook 09).
- *Borrow:* the downstream comparison layer should stay specimen-aware (Lamian-style), and
  domain deviations should be estimated **at a frozen coordinate**, not by re-fitting position per domain.

## 3. What our own evidence already establishes (ledger, NB15–19)

1. Count-exposure observation models make positions robust to read thinning. Under 25% read
   retention, the NB likelihood moves sections 0.02 versus 0.07 under the Gaussian surrogate.
   This does not improve gene prediction.
2. Within-segment molecular programs transfer across the two controls beyond a coverage oracle on
   training-selected panels: about 4–18% gains versus both linear and quadratic coverage oracles.
   Whole-transcriptome gains are tiny.
3. Strong coarse baselines (soft segment mixture) explain much of the signal.
4. PC1 initialisation locks onto a dominant nuisance program and fails to converge (NB19).
   Reversed nuisance coupling silently flips order. Uncertainty, model agreement and thinning
   stability do not detect the flip.

So the method needs **a source of positional information that a molecular nuisance program cannot
mimic** (anchors, physical adjacency, external calibration). It also needs occupancy invariance.
More expression-model flexibility alone is not enough.

## 4. Proposed method and hypotheses

### 4.1 Formulation: anchored unbalanced transport of sections onto a canonical line

Let sections $i$ have counts $x_i$, exposure $\ell_i$, specimen $s_i$, domain $d_i$, coarse label
$a_i\in\{S1,S2,S3\}$, and tissue-adjacency weights $w_{ij}$. Use a grid $z_1<\dots<z_K$ on
$[0,1]$ with a canonical count-rate atlas $\mu(z)$; this is the existing
`count_rate_atlas` / `pseudostructure_atlas` stack. Estimate a coupling $P\in\mathbb{R}_+^{N\times K}$
(section × grid) by minimising

$$
\underbrace{\sum_{ik}P_{ik}\,[-\log \mathrm{NB}(x_i\mid \ell_i,\mu(z_k)+\alpha_{d_i}+b_{s_i})]}_{\text{exposure-aware feature cost (P1)}}
+\lambda_a\underbrace{\sum_{ik}P_{ik}\,c_{\mathrm{ord}}(a_i,z_k)}_{\text{ordinal anchors (P3)}}
+\lambda_s\underbrace{\sum_{ij}w_{ij}\sum_{kl}P_{ik}P_{jl}(z_k-z_l)^2}_{\text{tissue-adjacency GW-type term (P5)}}
+\tau\,\mathrm{KL}(P^\top\mathbf 1\,\|\,\nu)-\varepsilon H(P)
$$

- **Row constraint** $P\mathbf 1 = \frac{1}{N}\mathbf 1$: every section is placed.
- **Unbalanced grid marginal**: KL to a broad reference $\nu$, with $\tau$ small. Occupancy is not
  forced uniform (P4, partial coverage). This differs from novoSpaRc.
- **Domain/specimen terms** $\alpha_d, b_s$ are restricted to *position-constant* gene offsets in
  the primary model, with healthy mouse $\alpha=0$. Position-dependent domain effects are estimated
  *afterwards*, at the frozen coordinate, by the existing GAM/level-shape machinery (P6).
- **Optimisation**: alternate between the coupling and the atlas. The coupling step is entropic
  (mirror descent / Sinkhorn on the linearised quadratic term, as in FGW solvers). The atlas step
  is a specimen-balanced weighted NB rate update. Initialisation uses the training-only
  S3−S1 endpoint contrast already specified in the ledger, not PC1.

Each term targets one property and is a separate ablation. With $\lambda_s=0$ and balanced
marginals, the method reduces almost exactly to our existing NB atlas projection. That keeps the
paper's incremental claims checkable.

### 4.2 Hypotheses (each with a falsifying outcome)

- **H1 — exposure.** DPT and scFates coordinates on log-normalised Harmony data carry
  within-segment association with library size and detected genes; the count-exposure estimator
  reduces it.
  *Falsified if* DPT/scFates show no larger within-segment coverage association than the new coordinate.
- **H2 — occupancy invariance.** Removing or duplicating sections from one region (for example,
  randomly dropping 50% of S1 sections) moves the coordinates of *untouched* sections under
  DPT/scFates (density-dependent operators and curve refits). Under the anchored unbalanced
  estimator with a frozen atlas, it moves them substantially less. This test needs no ground truth:
  a physical section's position does not change when other sections are discarded.
  *Falsified if* DPT/scFates displacement is comparable to ours.
- **H3 — tissue adjacency.** Physically adjacent cross-sections carry positional information that
  expression alone does not. Adding $\lambda_s$ improves held-out-gene prediction and within-segment
  program transfer across specimens. This holds against two controls:
  1. degree-preserving shuffled adjacency;
  2. a spatially smooth nuisance check, because neighbours also share local tissue quality.

  *Falsified if* the shuffled-adjacency control gives the same gain, or if the gain tracks
  spatially autocorrelated coverage.
- **H4 — identifiability under nuisance.** Anchors plus endpoint-contrast initialisation recover
  order when a strong nuisance program is present, where DPT, scFates and PC1-initialised atlases
  fail. Tested **semi-synthetically on real counts**: spike a reproducible non-spatial program into
  real healthy-control counts, either independent of the coarse labels or coupled to them in one
  specimen only. Measure agreement with the *un-spiked* real coordinate and with the coarse order.
  *Falsified if* the anchored estimator also tracks the spike.
- **H5 — domain separation.** On AKI and human sections, a frozen healthy template with
  position-constant domain offsets keeps the anatomical order: agreement with the upstream
  segment-label order and with the external Chen 2021 S1/S2/S3 calibration. It does so better than
  DPT and scFates refitted on Harmony-pooled data. Human cortex-only sections are not stretched to
  the full axis.
  *Falsified if* DPT/scFates order AKI/human sections at least as well by these criteria, or if
  injury-program genes (e.g. `Havcr1`, `Krt20`, `Krt8`) predict the new coordinate within segment
  as strongly as they predict DPT.
- **H6 — downstream stability.** Program/pathway comparisons (control vs AKI, human vs mouse)
  computed on the new coordinate are more stable under leave-one-specimen-out (LOSO) and gene-panel
  perturbation than the same comparisons on DPT/scFates. Discoveries that change are reported, not
  hidden.

## 5. Real-data experimental plan

All experiments use the private data via `--data-root/--results-root`. Inputs are the validated
PT structures already loaded by `pseudospace.repeated_inputs`, plus polygon geometry for adjacency.
Baselines are recomputed **on the identical structure set**:

- DPT as in notebook 03 (pass-2 Harmony, PT root, marker orientation);
- nonbranching scFates as in notebook 13;
- the existing count-Gaussian and NB atlases;
- the soft-segment and coverage oracles.

| Exp | Data | Question | Primary metrics |
|---|---|---|---|
| E0 | All | Reproduce DPT/scFates on the PT set; build tissue adjacency (shared polygon boundary or centroid kNN within specimen) | Agreement with saved 03/13 coordinates; adjacency degree, distribution of label agreement across adjacent pairs |
| E1 | 2 controls | H1 coverage association | Within-segment Spearman of coordinate vs log library, detected genes, `n_spots` |
| E2 | 2 controls | H2 occupancy invariance | Displacement of untouched sections after region-selective subsampling and duplication (5 seeds × 3 regions) |
| E3 | 2 controls, both transfer directions × 3 gene folds | Core accuracy and H3 | Held-out-gene error; within-segment transfer vs coverage oracle (NB18 protocol); λ_s ablation with shuffled-adjacency control |
| E4 | 2 controls + spiked nuisance | H4 | Rank agreement with un-spiked coordinate and coarse order; convergence |
| E5 | External: Chen 2021 S1/S2/S3 | Independent calibration | Ordering of the projected bulk segment profiles; correlation of inferred position curves with segment-ordered bulk means for non-marker genes |
| E6 | AKI mice, human cortex | H5 | Coarse-order agreement; injury-gene association within segment; occupied axis range; Gaussian/NB disagreement |
| E7 | All | H6 + the paper's biology | Re-run 05/12-style comparisons on each coordinate; LOSO and gene-panel stability of ranked programs |

Prespecified rules, matching ledger practice:

- No hyperparameter ($\lambda_a,\lambda_s,\tau,\varepsilon$) is chosen by physical depth or by query
  gene error. Choose them on training-only cross-fits, then freeze them before transfer evaluation.
- Failed fits are retained and reported.
- Depth and glomerular distance remain descriptive only.

## 6. Paper framing (if the hypotheses survive)

"An anchored, exposure-aware transport estimator for reconstructing a repeated anatomical axis
from segmented spatial cross-sections." The contributions would be:

1. the problem formulation (P1–P7);
2. an estimator whose components each have an ablation;
3. a ground-truth-free validation protocol: specimen transfer, within-segment programs, occupancy
   invariance, read thinning, semi-synthetic nuisance and external microdissection calibration;
4. the AKI and cross-species program comparisons on the frozen coordinate, against DPT and scFates.

Claims that the cohort cannot support (population-level human inference, confirmatory significance
with 2 vs 2 mice) stay out.
