# Pseudostructure reconstruction research ledger

## Target and evidence boundary

Recover a canonical ordered anatomical coordinate from unordered molecular sections of many realizations of the same physical PT structure. Specimens are independent tissue realizations; physical nephron identities are unknown. Multiple sections may come from the same nephron, so profile count does not establish biological replication. Healthy mouse is the initial estimator development cohort. Disease and species may change molecular distributions at a conserved coordinate; that conservation remains a hypothesis to validate, not a constraint that proves itself.

The intended contribution is a problem-specific reconstruction and validation procedure, not a claim that latent curves, hierarchical offsets, ordinal likelihoods or probability flow are new mathematics. Design choices must earn their complexity through withheld information, especially held-out specimen transfer and within-segment prediction.

## 2026-10-03 — From reference atlas to joint reconstruction

Measured notebook 15 results favor a canonical mean progression over the tested state-dependent flows. The reference mean atlas beats the segment-only gene predictor only slightly (approximately 0.19% relative error improvement), and it still inherits training DPT. Gaussian covariance evolution has direction-dependent transport benefits and worsens coordinate-based held-out gene prediction. These findings do not establish a latent estimator independent of DPT.

Notebook 16 tests the next hypothesis: equal-specimen marginal fitting of a shared smooth mean, weak learned ordinal anatomy likelihood and constant shrunken specimen offsets can reconstruct a useful coordinate without DPT. Each unknown nephron realization contributes to position-specific residual variation. Uniform finite-grid coordinate probability is a computational gauge and simplifying prior; neither uniform sampling nor physical longitudinal length is claimed.

Prespecified first screen:

- Known marker genes excluded; three deterministic gene folds.
- Joint fitting uses 80% of profiles from each control and expression-only projection of the remaining profiles. This is interpolation, not independent specimen validation.
- Separately fit one control and project the other with frozen representation and no query anatomy, DPT, physical depth or offset calibration; swap directions.
- Predict held-out genes using shared training-reference curves and compare with an explicitly label-informed segment oracle using equal specimen contributions.
- Report physical depth after fitting, including within-segment correlations; depth never selects a model.
- Test overlapping weak-ordinal random initializations, specimen-offset and anatomy-likelihood ablations. Overall agreement is insufficient if fine order is unstable within segments.
- Monitor the complete penalized observed Gaussian objective, convergence and failures. Conditional latent posteriors are not calibrated anatomical confidence intervals.

All runs and intermediate data remain in the private results directory. Notebook 15 retains its existing outputs; notebook 16 adds a separate experiment rather than overwriting the previous evidence.

## Relevant methodological precedents

[PhenoPath, Campbell and Yau, Nature Communications (2018)](https://www.nature.com/articles/s41467-018-04696-6) jointly models latent progression and measured covariates, including static and progression-dependent molecular effects. It establishes that a shared latent axis with covariate terms is not by itself a novel formulation. Our proposed use instead targets repeated physical anatomy, restricted registration and position-versus-domain-deviation validation.

[MEFISTO, Velten et al., Nature Methods (2022)](https://www.nature.com/articles/s41592-021-01343-9) is a precedent for smooth shared factors, group variation and alignment across continuous covariates. Its treatment of continuous covariates does not remove our need to infer unknown anatomical position and test it independently.

[novoSpaRc, Nitzan et al., Nature (2019)](https://www.nature.com/articles/s41586-019-1773-3) is a precedent for probabilistic spatial reconstruction from expression. Spatial reconstruction itself is not a new claim; the question here is whether a shared ordered object can be recovered across repeated realizations while retaining position-specific molecular differences.

These are a small design-oriented source check, not an exhaustive novelty review. Conclusions about uniqueness require broader comparisons once the estimator demonstrates utility.

### First joint-latent run: measured result

All 13 fits converged by the specified observed-objective criterion, with no objective increases or recorded fitting failures. Joint profile-holdout gene error improved over the segment oracle by 0.56–0.61%; independent specimen transfer improved by 0.43–0.60%, across the three gene folds and both directions. These small gains do not establish spatial specificity.

On matched fold-0 projected profiles, mean-three-glomeruli depth Spearman correlations were 0.667/0.637 for transfer into Ctrl1A2/Ctrl1A4, compared with supplied DPT 0.725/0.688 on the same eligible profiles. Joint profile-holdout correlations were 0.696/0.662, compared with matched DPT 0.758/0.714. Within-segment S2/S3 depth associations were generally weaker than DPT.

A random start (seed 29) agreed almost perfectly in rank with PCA initialization, but seed 15 changed S1 ordering substantially (within-S1 Spearman 0.39/0.56). Its training objective was worse (43.074 versus approximately 42.75 for the other two fold-0 starts), so training-objective multistart selection can identify this local solution without selecting against validation data. This is not proof that all alternative anatomical orders are excluded.

The nominal offset penalty of 4 shrank offsets to norms around 0.05 in the modeling PC units; removing offsets produced almost identical scores. That ablation does not establish that specimen offsets are generally unnecessary: it shows this particular strong penalty effectively suppressed them. Removing the weak anatomy likelihood had small score effects; its PCA orientation still inherited training endpoint labels.

The decisive diagnostic is coverage association. Within-segment ordering correlated strongly with detected-gene count and aggregate library size, including detected-gene Spearman approximately −0.94 in held-out Ctrl1A4 S2 and +0.94 in held-out Ctrl1A2 S3. Overall associations were small because directions reversed between segments. Coverage could include real morphology as well as measurement variation, but these findings make a fine anatomical interpretation unsafe. A generic smooth mean manifold can predict genes while ordering a nonspatial amplitude/detection program.

### Next hypothesis: observation model before a more flexible atlas

Test a frozen negative-binomial Pearson-residual representation of raw aggregate counts with library exposure; retain the shared latent spatial estimator and validation boundaries. The reference gene probabilities, clipping and PCA are fitted on modeling genes and training profiles only; query data do not refit them. Total library exposure uses all measured genes, so gene holdouts retain the existing compositional dependence. Residual PCA does not apply a second per-gene variance standardization. The pilot is an exploratory response to observed failure, not a new independent validation cohort or evidence that coverage is entirely technical.

[Analytic Pearson residuals, Lause, Berens and Kobak, Genome Biology (2021)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8419999/) derives a closed-form count normalization under a negative-binomial offset model and examines it for single-cell UMI data. Its use on aggregated tubule counts here is a methodological test, not an established transfer of the paper's biological conclusions. A fixed dispersion of 100 and training-only clipping are prespecified for this first pilot; neither depth nor gene-prediction outcomes choose those parameters.

### Count-exposure atlas: rotated-fold result and interpretation boundary

All 12 count-atlas fits converged: three joint-fold fits, six independent specimen-transfer fits, two additional random starts and an endpoint-orientation-only ablation. Joint count-atlas gene prediction beats the segment oracle in all folds and both specimens, but gains are small; independent transfer gains are also consistently positive and smaller than the log-normalized latent model. The latter comparison is compatible with losing coverage-predictive structure, but does not prove it.

The count atlas removes most S2 count/detection association and stabilizes all three tested joint initializations within every segment (rank correlations above 0.999). One transfer direction retains a substantial S3 count/detection association, so the observation nuisance is not solved in general. Removing the ordinal likelihood, while retaining training S1/S3 PCA orientation, has little effect on fold-0 held-out gene errors; S2 supervision is unnecessary for that tested fit.

**User steering, 2026-10-03:** cortical depth is an especially imperfect supporting measurement. A nephron winds through tissue, and independent cross-sections can sample different nephrons; nearest/mean-three glomerular distance is not true longitudinal position. Do not rank methods or interpret higher depth correlation as improved anatomical reconstruction. Keep the measurements as descriptive checks, subordinate to cross-specimen molecular reproducibility, withheld within-segment programs, stability and observation-coverage robustness. Neither the representation nor its parameters were chosen by depth optimization.

Next experiment: binomially thin counts of the exact same held-out physical sections, freeze each reference representation/atlas, and quantify coordinate shifts, within-segment rank stability and conditional posterior behavior. The physical section does not move when its reads are thinned. Compare log-normalized and count-exposure atlases on identical query sections. This is a measurement robustness check, not proof of correct anatomical fine order.

### Same-section thinning: a direct observation-robustness result

At 50% retained reads, count-atlas mean absolute position shifts are approximately 0.034 in each control, versus 0.100/0.079 for the log-normalized atlas. At 25% retained reads, count-atlas shifts are 0.073/0.069 versus log-normalized 0.235/0.183. Conditional entropy increases for the count atlas and decreases for the log-normalized atlas, despite its larger movements. These are descriptive perturbations of the same query sections, not additional biological replicates or calibrated uncertainty tests.

The count atlas remains inadequate under severe thinning: 21–27% of sections move by more than 0.1 at 25% coverage; S3 within-segment rank stability averages only 0.64. Point stability can differ from uncertainty changes, and the representation plus isotropic Gaussian likelihood is not a full exposure-conditioned count observation model. A frozen canonical gene-rate atlas with count likelihood is a possible next simplification/test, rather than adding a more flexible molecular vector field.

Before further complexity, compare the same count PC1 plus identical gene-curve predictor and an expression-only soft three-segment mixture. The latter can model softened label boundaries without discovering fine within-segment anatomy. Fit only training labels and model genes; do not feed query labels into these predictors. This will distinguish continuous within-segment contribution from a normalization improvement or a better coarse classifier.

### Stronger baselines narrow the claim

The count atlas beats the frozen count PC1 plus gene-spline predictor in every tested gene fold and protocol, with relative aggregate-error improvements approximately 0.11–0.58% for independent specimen transfer. It also beats the soft segment mixture for joint profile interpolation and for Ctrl1A2→Ctrl1A4 transfer. In Ctrl1A4→Ctrl1A2 transfer it nearly ties that stronger coarse predictor and loses one fold by approximately 0.006%.

Within-segment comparisons are heterogeneous: S3 prediction favors the continuous atlas, while the soft segment mixture is better for some S1/S2 transfer comparisons. This does not establish a uniformly superior continuous anatomical method. The reconstruction is a promising candidate with small incremental molecular prediction benefits; softened coarse labels explain a substantial part of the signal. The count representation and its measurement robustness currently have clearer support than a large additional advantage from fine latent position.

Next bounded prototype: freeze a shared gene-rate path estimated from training count-atlas positions, and project raw counts using a library-conditioned NB likelihood. Estimate local rates separately per original specimen, then average equally at each position to avoid confounding canonical molecular means with specimen-specific position occupancy. Unknown nephron variation is only approximated by dispersion. This is a decoder/projection pilot before any full latent-count fitting; do not claim its conditional posteriors are calibrated, or that density improvement identifies real longitudinal position. The same-section thinning and cross-specimen/strong-baseline comparisons decide whether this observation-model upgrade earns its complexity.

## 2026-10-03 — Conditional-count projection pilot (notebook 17)

The canonical rate path uses the training-only count-Gaussian axis as initialization; DPT is absent. Local rates are estimated separately per original specimen and averaged equally at every grid point. Raw query counts are scored conditional on all-gene library exposure with fixed NB dispersion 100. This is a frozen decoder/projection pilot, not a fully joint latent-count estimator or a calibrated posterior for unknown physical nephron positions.

On the identical joint-held-out sections, retaining 25% of reads gives NB mean absolute coordinate movements 0.020/0.015, versus count-Gaussian 0.072/0.069. The fraction moving more than 0.1 falls from approximately 22–25% to 0.2–1.8%. Within-S3 rank stability improves from approximately 0.62 to 0.96. S2 rank stability is slightly weaker under NB in this perturbation. Conditional entropy increases with thinning, but the lower baseline NB entropy may reflect the unverified independent-gene likelihood rather than justified precision.

Held-out gene MSE is slightly worse than count-Gaussian in all four fold-0 protocol/specimen combinations (relative differences approximately 0.03–0.07%). Both retain small gains over count PC1 and the segment oracle; soft segment prediction remains nearly tied in one independent transfer direction. The observation likelihood has earned further testing through measurement robustness, not through superior molecular prediction or anatomical proof.

All three count-rate references fit without recorded errors, and minimum per-specimen kernel effective support is well above the prespecified threshold. That is a numerical support diagnostic, not evidence of known true positional coverage. The repeated-structure prior, imperfect anatomy labels, unknown nephron identities, fixed coordinate gauge and two-control replication limit remain unchanged.

Next research choices: test rotated count-rate gene folds and reference perturbations; assess count-likelihood calibration and dispersion misspecification; investigate groupwise latent-count refinement only with objective/stability safeguards. Preserve strong coarse/PC1 baselines. Disease/species projection must remain on a frozen healthy structural coordinate and report abnormality/ambiguity separately; none of these pilot results guarantees correct location for injured or cross-species sections.

### Next test specified before dispersion-sensitivity outputs

Rotate all three deterministic gene folds for the joint profile-holdout reference and hold each learned canonical rate path fixed. Compare NB dispersion 10, 100 and 1000, plus the Poisson limit, without selecting a parameter from query gene scores or physical depth. Repeat quarter-read thinning on the same query sections and report within-segment ordering, point movement, entropy and cross-dispersion agreement. This tests the observation approximation; it is not a new independently validated cohort or a full latent-count optimization.

A stable mean coordinate with dispersion-sensitive posterior width would support reporting reference-based ordering separately from conditional uncertainty. Disagreement in fine order would instead expose observation-model dependence. Neither outcome establishes true longitudinal location. Full latent-count fitting should follow these checks, rather than repeatedly reshaping the atlas against the same evaluation sections.

### Dispersion and gene-panel sensitivity: measured results

The initial sensitivity attempt stopped with strict coordinate-domain failures: a sharply peaked posterior could sum slightly above one after exponentiating log-normalized likelihoods, yielding an expected position just above 1. The shared projector now renormalizes probability rows explicitly, clips only the convex mean's floating-point endpoint overshoot, and computes entropy with zero-safe `xlogy`. A deterministic synthetic count-likelihood regression reproduces the former overshoot. All 202 synthetic tests pass, and the complete notebook rerun records no failed sensitivity fits (12 settings and 72 per-segment dispersion comparisons).

Across all three gene folds and both joint-held-out specimens, quarter-read mean coordinate movements are 0.016–0.019 for the four count observation settings, versus 0.069 for the Gaussian surrogate. The fraction moving over 0.1 is approximately 0.7–1.2%, versus 22.7% for Gaussian. Within-S3 thinning rank agreement is 0.95–0.97, versus 0.64. These averages describe repeated perturbations of the existing sections, not independent biological replication.

Changing NB shape from 100 to 10, 1000 or the Poisson limit moves coordinates by approximately 0.006–0.009 on average. Within-segment mean rank agreement across these choices is approximately 0.98–1.00; conditional entropy changes appreciably and remains uncalibrated. Count-likelihood gene MSE remains slightly worse than Gaussian in every fold/specimen combination, though better than frozen PC1 and soft segment mixtures in this joint interpolation experiment. This robustness check does not erase the near-tie with soft segments in independent specimen transfer. No setting was selected using these outcomes.

Gene-panel perturbation is more consequential than dispersion changes: NB shape-100 coordinate agreement across pairs of folds is high overall but, within Ctrl1A2 S3, mean rank agreement is 0.904 and the weakest pair is 0.857. Ctrl1A2 S2 averages 0.919. Ctrl1A4 S2/S3 average approximately 0.984/0.982. Gene panels overlap substantially and share all-gene exposure; these are useful sensitivity checks, not independent validations. Robust observation handling does not imply an identified fine anatomical coordinate.

Next priorities are training-reference perturbation and a more focused test of transferable within-segment molecular programs. Retain the repeated-structure prior and strong coarse baselines; do not optimize depth or advance AKI/human/TF claims on these results. A full latent-count estimator remains a candidate, not an automatic upgrade: the current frozen path inherits the count-Gaussian training axis, and moving training positions to improve its own count fit could reinforce a nonspatial program.

Visual inspection of all gene-fold pairs shows that large shifts concentrate in a small set of profiles: 14/548 joint-held-out Ctrl1A2 sections and 4/556 Ctrl1A4 sections have a three-panel coordinate range above 0.1; 8 Ctrl1A2 sections exceed 0.3, with a maximum range of 0.913. Some of these profiles map almost exactly to an endpoint with entropy below 0.1 in one panel, yet map near the opposite endpoint under another panel. The independent-gene likelihood can therefore be sharply conditional on an unreliable reconstruction. Notebook 17 now reports per-section gene-panel range alongside conditional entropy, without excluding any section or claiming calibrated interval coverage. Unknown-nephron/anatomical specificity is still unresolved. Reference-support and perturbation-based ambiguity must accompany final location estimates.

## 2026-10-03 — Within-segment transfer test specified before outputs

Test the count-Gaussian canonical axis across the two independent controls in both directions and all three gene folds. This is the simpler atlas with slightly stronger molecular prediction, not a new flexible flow or a full latent-count upgrade. Query inference remains expression-only. Permute training coordinates independently within original specimen×segment for five fixed seeds, refit the identical withheld-gene spline decoder, and keep query coordinates fixed. This retains each segment's coordinate distribution while removing its fine response ordering. Report descriptive prediction differences, not biological permutation p-values.

Include count PC1, soft segment prediction, and a deliberately privileged segment-plus-coverage oracle. The latter uses query segment labels plus log total library, log1p spot count and log1p detected modeling-gene count, with training-only, within-segment covariate normalization and equal original-specimen weighting. It diagnoses whether a count/structure-amplitude predictor explains the apparent fine signal; spot count and coverage may include real morphology and cannot automatically be declared technical artifacts. Total all-gene library preserves the disclosed compositional dependence of held-out responses.

Primary evaluation retains every withheld gene. Supplementary candidate programs are selected separately per training segment using a fixed two-way profile cross-fit of response decoders along the frozen modeling-gene axis: rank positive predictive gains over the coverage oracle and take at most 50 genes. Freeze selection before the independent query evaluation. Training profile splits are not independent biological replicates, nor are selected genes independent repetitions. This is a focused exploratory test of transferable within-segment molecular organization, not proof of longitudinal anatomical location. Previous gene-panel ambiguity remains a limitation; no unstable section is silently dropped.

### Within-segment transfer: measured result

All six count-Gaussian references converge; no transfer failures are recorded. Actual coordinates beat the mean of the five within-segment training-order shuffles in every segment, fold and transfer direction. Whole-gene error gains are small: approximately 0.57% for Ctrl1A2→Ctrl1A4 and 0.21% for the reverse direction. Against the stronger coverage oracle, whole-gene atlas error is approximately 0.16%/0.38% worse. This rules out presenting the atlas as a uniformly better expression decoder.

The supplementary training-selected programs transfer more clearly. Exactly 50 positive-CV-gain genes are selected per segment per fold in these runs, without query-based selection. Across folds, atlas error gains versus the coverage oracle are approximately 14%, 18% and 16% for S1/S2/S3 in Ctrl1A2→Ctrl1A4, and 9%, 10% and 4% in the reverse direction. Every segment/fold/direction shows positive selected-panel gains versus both the coverage oracle and the mean coordinate-shuffle null. These are descriptive effects in two controls, not confirmatory significance or anatomical ground truth.

Individual selected-gene improvement frequencies versus the coverage oracle are approximately 85–98% in the first transfer direction and 77–92% in the reverse direction. Among the 150 selected genes across disjoint response folds per segment, the separately trained lists overlap by 71 (S1), 107 (S2) and 54 (S3). That overlap is a consistency description; neither response list is used to reselect the independent query benchmark. The axis is built without those response genes, but global all-gene library exposure retains compositional dependence.

Interpretation: there is transferable continuous molecular organization within the supplied coarse anatomy beyond these coverage covariates and coarse response means. This is stronger support for the repeated-structure atlas hypothesis than tiny global gene-error gains alone. It does not prove that every fine ordering corresponds to longitudinal anatomy; shared nonspatial programs, unknown nephron identity, restricted linear nuisance modeling and conditional-posterior ambiguity remain. The original probability-flow extension still has no validated advantage.

Next method-development step: consolidate the selected architecture into a reproducible repeated-structure estimator with explicit count exposure, equal-original-specimen canonical fitting, weak anatomical orientation, frozen expression-only reference projection, and per-profile reference/gene-panel sensitivity. Keep canonical position and domain deviation separate, preserve partial coverage, and validate the estimator itself on synthetic repeated structures plus the existing control-transfer benchmarks. A fully joint count-position objective is optional if it adds measured value; do not invent it simply to replace a working two-stage estimator.
