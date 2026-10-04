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

### Next check specified before nonlinear-oracle outputs

The current coverage oracle is linear in log library, log1p spot count and log1p detected modeling genes within each supplied segment. A latent curve may predict nonlinear coverage relationships even when its apparent fine signal is nonspatial. Strengthen that diagnostic with training-only centered/scaled quadratic terms and interactions, retaining the same segment label privilege and frozen query transforms. Keep the notebook-18 training-selected panels fixed for the first comparison; do not choose new genes after inspecting query errors. Compare primary all-gene and these original panels on the same controls/response folds. If nonlinear coverage explains the gains, treat that as a model-development result, not a reason to silently change the benchmark. Consolidation of the estimator should preserve this diagnostic limit.

### Nonlinear coverage check on the original panels

The prespecified quadratic oracle adds squares and pairwise interactions to the same three coverage covariates, with training-only specimen-balanced scaling and ridge 0.01. The original gene selections and benchmark scores reproduce exactly; neither atlas positions nor query-based gene selections were changed. All six oracle fits complete without failures.

Across folds, original selected-panel atlas gains versus the quadratic oracle are approximately 13.2%, 18.1% and 9.8% for S1/S2/S3 in Ctrl1A2→Ctrl1A4, and 8.4%, 9.9% and 4.0% in the reverse direction. The S3 gain in the first direction falls from roughly 16% with the linear oracle to 10% with the quadratic oracle. Nonlinear coverage therefore explains part of that result, while the tested original panels retain transferable molecular structure beyond this diagnostic. This does not exclude other nonlinear, morphological or nonspatial explanations.

The quadratic oracle is slightly worse than the linear oracle in whole-gene transfer; adding features is not evidence of a universally stronger predictor. Whole-gene atlas gains against quadratic coverage remain tiny (approximately 0.10%/0.07%). Around 3–8% of query profiles lie outside the training 1st–99th percentile range of at least one coverage covariate; none were removed. These comparisons remain descriptive in two mice, not anatomical validation. Physical-depth correlations are secondary context only: winding nephrons and unknown section membership prevent interpreting them as longitudinal coordinate accuracy.

### Estimator consolidation and next identification test

The prototype public interface combines training-only count residual PCs, a specimen-balanced shared Gaussian mean atlas with restricted constant offsets and weak ordinal anatomy, and a frozen specimen-balanced local count-rate decoder. Unknown nephron identities remain integrated into the residual approximation. Query counts and exposure are projected without query labels, coordinates, depth, PCA refitting or range rescaling. The count-rate stage uses the Gaussian training axis; it is not a joint latent-count fit. Both conditional posteriors are retained, with model disagreement and aligned gene-panel ranges as descriptive sensitivity rather than calibrated intervals.

Before domain projection, the next bounded identification experiment should freeze this architecture and challenge **partial coverage** and competing repeated molecular programs on synthetic physical structures. A truncated query population must not be stretched to occupy [0,1]; count/exposure shifts should not silently redefine position. Include a reproducible nonspatial program independent of true position, and a deliberately ambiguous case where it is correlated with position in training but decoupled in query. Report recovery against known simulated position and failures separately from real-control molecular prediction. These tests investigate what the design can identify, not whether the simulation proves the biology. Keep AKI, human and regulator analysis behind this identification work.

Notebook 18 verifies the public interface on Ctrl1A2→Ctrl1A4, gene fold 0: all four Gaussian withheld-gene benchmark scores reproduce within relative tolerance 1e-8 across 2,782 query structures. Independent synthetic query draws recover simulated broad ordering without querying training labels, and sparse-input/batch projection reproduces full-batch results. The 216-test synthetic suite passes. These checks establish implementation consistency, not additional biological replication.

In that real transfer, Gaussian-versus-count mean absolute coordinate disagreement is 0.0268, with a maximum 0.845. Conditional count entropy is lower on average, but the extreme disagreement reinforces that a sharp model-conditional posterior is not sufficient evidence of reliable location. Preserve both estimates and sensitivity information rather than selecting the more confident decoder by its entropy.

## 2026-10-03 — Synthetic partial-coverage and competing-program test specified before outputs

Freeze the public estimator and its defaults before simulation results. Generate independent, unordered cross-sections from two training specimens and a third query specimen, retaining true position only for evaluation. Gene counts follow a multinomial exposure model with a shared spatial logit direction and an orthogonal nuisance direction. Independent per-section gene perturbations approximate unlinked realization variation; this does not simulate tracked nephron functions or prove unknown sections are independent in tissue.

Use three fixed seeds (19, 41, 73), 500 training sections per specimen and 400 query sections. Compare three training scenarios: spatial signal alone; an equally parameterized independent nuisance program; and nuisance perfectly coupled to position in training. For the last case, query both preserved and reversed coupling. Training anatomy is obtained by thirds of simulated position solely as coarse ordinal supervision; true fine position never enters fitting, initialization or projection. The artificially exact coarse labels make this an optimistic identification test.

For each reference, query full and restricted [0.35,0.65] coverage, and paired 25%-read count thinning with matching exposure. Report true-position rank agreement, endpoint saturation, conditional entropy, Gaussian/count disagreement and same-section thinning movement. Compare frozen count PC1 and an expression-only soft segment predictor, with fixed coarse ordinal response codes rather than true fine positions. Query labels remain evaluation-only. Verify that projecting a subset from the full query batch reproduces the identical section coordinates exactly; a separate truncated population cannot be rescaled to fill the axis.

Do not tune parameters against simulated recovery or discard failed fits. A repeated nonspatial program can be reproducible across specimens, so cross-specimen consistency alone is insufficient to identify spatial anatomy. Reversed coupling deliberately violates the frozen shared molecular realization and tests the ambiguity that domain-specific changes could cause. Failures should constrain claims and motivate explicit identification safeguards, not be hidden by a more favorable simulator.

### Synthetic identification stress: measured results (notebook 19)

Six of the nine prespecified references converge. All three independent-nuisance references fail the fixed 60-iteration convergence guard and are not accepted for atlas projection. Their PC1 initialization correlates approximately 0.997 with nuisance and only 0.02–0.05 with true position. Final objective changes remain outside tolerance. This demonstrates a failure of the current initialization/solver settings, not nonexistence of a spatial atlas or proof that running longer would fix it.

In spatial-only data, Gaussian/count full-axis position ranks average approximately 0.995/0.997, with partial-population ranks approximately 0.972/0.972. Partial Gaussian positions occupy roughly 0.35–0.64 rather than the entire axis. Every same-profile subset projection is exactly unchanged (maximum change zero). These checks support frozen population-invariant projection under this simulator; they do not identify physical length or establish correct mapping under domain change.

Coupled training/query nuisance yields apparently excellent full-axis Gaussian/count ranks (0.997/0.999). Reversing the nuisance coupling in query makes them approximately −0.991/−0.997; partial-query ranks are approximately −0.978/−0.974. Gaussian/count disagreement remains small (partial average 0.006) despite reversed ordering. Conditional uncertainty and inter-model agreement therefore cannot identify all failures. Quarter-read count coordinate movements remain small even in the reversed case: measurement robustness is distinct from correctness.

The underlying failure diagnostic retains the same objective, tolerance and iteration limit. On the independent-nuisance queries, frozen PC1 has essentially zero position rank, whereas the soft segment baseline reaches approximately 0.992 overall and 0.957 in the partial population. Within-segment full-population ranks are approximately 0.93/0.96/0.93 for that baseline. The representation contains spatial information that the PC1 initialization does not select. Its coarse labels are artificially exact and the simulator is deliberately simple; this is not evidence that soft segmentation solves real fine anatomical reconstruction.

All 227 synthetic tests pass. Notebook 19 retains failed fits, independent-nuisance objective histories, reference/baseline results, and a separate population-invariance guard. The cached notebook must also run successfully from existing cached JSON arrays; this guard checks reproducibility, not biology.

### Next initialization ablation specified before outputs

Keep the count representation, shared-mean objective, ordinal anatomy strength, constant specimen offsets, grid, tolerance and 60-iteration limit fixed. Compare the existing PC1 initialization with (i) a training-only S3-minus-S1 endpoint-mean contrast and (ii) the existing soft-segment predictor's expected ordinal code. Balance original specimens when forming endpoint means. Freeze all initializers and their numerical conventions using training data only; do not use true fine simulated position or query labels to choose an initialization.

First repeat all nine synthetic references and report convergence, solver objective and projected rank against the existing frozen scenarios, including reversed nuisance and partial coverage. Do not accept a nonconverged fit or choose the winning initialization by query recovery. This separates a PC1 initialization defect from limitations of the shared Gaussian likelihood. If an ordinal initialization helps the independent-nuisance case, validate it on the unchanged healthy-control gene-fold benchmark before adopting it. It cannot resolve the deliberately unidentifiable domain-coupling case by itself. More flexible flow or residual covariance is deferred until this simpler cause is tested.

### Ordinal initialization ablation: measured results (notebook 20)

All 27 prespecified solver attempts are retained. Endpoint contrast converges in all nine references; existing PC1 converges in six, and soft-segment initialization in three. Of the latter three, one is spatial-only and two are independent-nuisance. Its recovery summaries therefore condition on successful seeds and must not be compared as if all nine succeeded. The three coupled-nuisance soft starts remain nonconverged at 60 iterations. Neither iteration limits nor label strength were changed.

For independent nuisance, endpoint initialization's true-position rank is approximately 0.982, but accepted Gaussian query rank averages 0.778 on full coverage and 0.775 on partial coverage. Full-query ranks range 0.769–0.789 across the three seeds. Count-rate query ranks are approximately 0.763/0.682. The frozen soft-segment baseline on the same datasets is substantially stronger (approximately 0.992/0.957). Endpoint contrast fixes the tested convergence failure but does not establish a better reconstruction method. Partial endpoint-started Gaussian gauge error averages approximately 0.23, another sign that passing population invariance does not ensure reliable placement.

The two successful soft-start independent-nuisance references yield Gaussian full-query ranks 0.932 and 0.941, but partial ranks 0.746 and 0.760; the third seed does not converge. Their observed training objectives are lower than the endpoint-started fits. This is evidence of local-start dependence and residual/model limitations; it does not isolate a single cause or justify selecting an initializer from query truth. A lower objective can coexist with weaker fine-position recovery than a coarse-label baseline.

Endpoint and PC1 results are effectively unchanged in the spatial-only and coupled cases. Reversed coupling still reverses ordering (Gaussian full-query rank approximately −0.991, partial −0.978). Different starts cannot resolve the intentionally unmodeled domain ambiguity. All accepted subset projections remain unchanged by query population composition. Known simulated fine position is used only for evaluation, not fitting or initialization.

The new training-only initializer helper remains experimental. The public estimator's default architecture is unchanged; no initializer has been adopted or selected on query recovery. All 234 synthetic tests pass. The notebook records objectives, convergence failures, initialization ranks, full/partial query recovery and paired read perturbations. Before adoption, the unchanged healthy-control gene-fold benchmark and further initialization/residual-variation ablations would need to distinguish a numerical-start improvement from a better repeated-structure estimator. No new healthy-control, AKI, human or regulator experiment was started in this iteration.

## 2026-10-03 — Resume with actual-data trajectory benchmarks, specified before outputs

Following user steering, biological method-development experiments now use the actual healthy controls; synthetic fixtures remain code regression checks. A targeted primary-source survey is recorded in `docs/research/pseudostructure-methods-survey.md`, covering DPT/scFates, principal-curve initialization, covariate-aware latent models, groupwise smooth factor models, temporal/spatial OT, partial transport and trajectory alignment. No comprehensive novelty claim is made.

Notebook 21 tests all three fixed atlas starts against DPT and nonbranching scFates on the identical marker-excluded, training-only count-residual representation. Use both healthy transfer directions and all three existing gene folds, with no query coordinate/label/feature refitting. All methods receive training anatomical orientation; the repeated atlas additionally uses weak ordinal labels in its likelihood, so supervision is not identical. DPT uses 30 training neighbors and an actual S1 structure nearest its training S1 centroid as root. scFates uses 30 principal-curve nodes in all 15 supplied components, roots at the tip nearest the training S1 centroid, and uses one seeded pseudotime mapping. This is a matched-input algorithm benchmark, not reproduction of notebook 13's pooled Harmony/marker workflow. Save parameters, installed versions, roots, failures and all coordinate support; never discard disconnected DPT points without reporting a failed fit.

The primary comparison extends every learned training axis to query expression using the same inverse-distance-weighted 15-neighbor projector, then predicts every held-out gene with the same training spline/ridge and variance-normalized MSE. This isolates axis fitting from query projection. The atlas's native Gaussian projector is a separately named secondary comparison on the same fitted axis; DPT/scFates plus the shared neighbor extension are not advertised as native out-of-sample APIs. Retain count PC1, soft three-segment prediction and the privileged segment-plus-coverage oracle. Compare error overall and within each supplied segment; evaluate markers separately without selecting genes or hyperparameters on query results. Depth is excluded from scoring and selection.

Hypotheses: training-only ordinal starts may improve actual-data transfer over PC1; repeated-atlas axis fitting may improve held-out molecular programs over DPT/scFates; native distributional projection may improve on neighbor interpolation. Report failures and null results rather than adopt the best-looking query method. These are repeated development experiments on the same two mice, not confirmatory tests of true longitudinal position. Further OT/covariance changes must follow the measured comparison, rather than precede it.

### Actual-data DPT/scFates comparison: measured results (notebook 21)

All six DPT and scFates training fits complete on current healthy inputs, with no non-finite coordinates or silently excluded profiles. All six PC1 and endpoint-started atlas fits converge. Soft ordinal atlas starts fail in Ctrl1A2 folds 0 and 2, leaving only one successful fold in that direction; their reported means are incomplete and not equivalent to three-fold method comparisons. All 30 attempted algorithm records and failures are preserved.

Through the common neighbor projector, scFates beats DPT in all three Ctrl1A2→Ctrl1A4 folds (mean all-gene error reduction approximately 0.284%), and two of three reverse folds with a near-zero mean difference (0.003%). The original joint atlas does not consistently improve on scFates: native PC1 atlas error is approximately 0.047% worse in the first direction and 0.022% better in the reverse. Common-neighbor atlas error is approximately 0.107%/0.005% worse than scFates. Endpoint and PC1 results differ negligibly; the initialization improvement seen synthetically does not produce a useful actual-data improvement here.

Native Gaussian atlas projection reduces error versus the same atlas axis with neighbor projection by approximately 0.060%/0.027%. That is a small projection effect, not proof of superior anatomical axis inference. Segment results are heterogeneous: the atlas improves some S2/S3 comparisons but loses others. The privileged segment-plus-coverage oracle is better overall than scFates by approximately 0.11%/0.40%. Atlas marker prediction is stronger, but coarse labels inherited marker information; it is not independent validation of spatial reconstruction.

This result supports using an established scFates backbone as a candidate special case of the repeated-object estimator, without adopting it solely because of these query scores. The methodological design can reside in specimen-balanced canonical distributions, count exposure, restricted nuisance terms, frozen partial-coverage projection and explicit ambiguity, rather than a novel curve optimizer. No OT or biological force field has earned its complexity yet. The healthy data and query evaluation are reused for development and remain descriptive.

### Next real-data hypothesis specified before outputs

On the same two healthy controls and three gene folds, learn a fixed scFates backbone from training modeling genes alone. Compare its existing shared-neighbor projection with a frozen mean-reference projector and a conditional NB count-rate reference using the existing bandwidth 0.08, dispersion 100 and support threshold 12. Estimate reference rates separately per original specimen when jointly fitted; the outer transfer has one training specimen. Do not jointly refit positions or tune these parameters against query genes.

Evaluate all held-out response genes, supplied segment-specific errors, and paired quarter-read thinning on the same query sections. Retain the current DPT/soft-segment comparison and verify that the scFates-neighbor benchmark reproduces. Hypotheses are that distributional projection may improve transferable gene programs and/or observation robustness on a conventional backbone. A failed prediction gain but better read robustness is distinct evidence; neither proves fine anatomical location. Preserve conditional posterior entropy alongside model/gene-panel disagreement. Subsequent OT should address specimen distribution matching or partial support only if this simpler reference improves the relevant behavior.

### Fixed scFates reference: first measured results (notebook 22)

All six source-only references complete and reproduce notebook 21's scFates-neighbor response scores within the prespecified numerical tolerance. Every projector uses the same decoder at the original training scFates positions; training rows are not reprojected to refit that decoder. The mean-reference projector increases all-gene error by approximately 0.039%/0.105% in the two directions. The NB reference increases it by approximately 0.253%/0.446%; both alternatives lose overall in every fold. Segment-specific results remain mixed. These findings do not justify replacing the neighbor projector for prediction.

The conditional count reference does improve paired quarter-read stability: coordinate movement is approximately 0.014/0.024 versus 0.043/0.047 for neighbor projection and 0.053/0.072 for the mean reference. Fractions moving more than 0.1 are approximately 1.3%/3.5%, versus 9.7%/9.8% for neighbors. NB conditional entropy increases with reduced reads, but it remains an uncalibrated model summary. Stable predictions can still be wrong.

Gene-panel coordinate ranges are substantially larger when training on Ctrl1A4 (mean range approximately 0.114–0.115 across methods) than Ctrl1A2 (0.031–0.034). Overall pairwise panel ranks remain approximately 0.963–0.991 in the former direction and 0.991–0.996 in the latter. Numerical range therefore combines axis stretching and changed ordering; it is not anatomical uncertainty or a confidence interval. The NB maximum range reaches approximately 0.99 despite small conditional entropy. A sharp conditional posterior does not resolve cross-model ambiguity. The two specimens and overlapping gene folds are reused development data.

### Focused program check specified before its outputs

Reuse the already frozen notebook-18 training-selected positional gene panels, extracting only source, target, fold, segment, gene and training-CV gain, without reading their historical query-error columns into selection. Their selection used the earlier atlas and a linear coverage oracle, so they are development-informed and can favor that earlier backbone; they are not a neutral universal trajectory benchmark. Do not change panels in response to notebook-22 predictions.

Score the three notebook-22 projectors on these identical panels with the unchanged scFates training decoder. Add the existing privileged quadratic segment-plus-coverage oracle, using all-gene library, structure size and detected modeling genes. Fit all coefficients on the source only; query anatomy enters that named oracle and segment scoring, never molecular projection. Recompute the source scFates axis and guard reproduction of the original neighbor scores and saved query positions. This distinguishes broadly averaged error from positional-program transfer; it cannot establish physical longitudinal position. Preserve all missing-panel/fit failures and require six complete source-fold records. No flow, AKI or human branch is enabled by this check.

The overall panel ranks mask weaker fine ordering: in the Ctrl1A4-trained S1 query sections, pairwise gene-panel ranks range approximately 0.583–0.907 for neighbors, 0.601–0.953 for the mean reference and 0.626–0.877 for NB. Thus stretching is only part of the sensitivity; neither overall rank nor small read perturbation validates within-segment order. Quarter-read gene errors against fixed original response measurements change only modestly: NB is slightly better than neighbor projection in one direction and worse in the other. Read robustness does not become a consistent prediction gain.

### Next residual-metric ablation specified before outputs

A mean reference's isotropic residual likelihood treats every count-PC direction as equally informative about position. The repeated-realization model instead permits structured variation around the canonical path. Before adding position-dependent covariance or transport, test the smaller special case of one training-estimated residual covariance, holding the scFates axis, kernel mean and response-gene decoder fixed. Estimate residuals at the original supplied training positions by interpolating the fitted mean; fit a Ledoit–Wolf covariance with zero residual mean assumed and a numerical positive-definiteness floor. Use its Gaussian likelihood to project query profiles, without query fitting, labels, depth or genes from the response fold.

Compare against the existing isotropic mean, neighbor and NB references on the same six real-control transfers. Require reproduction of saved neighbor positions/scores, retain frozen positional-program panels and the quadratic coverage oracle, and evaluate paired quarter-read perturbation. Parameters remain training-only; no query-based shrinkage sweep or winning-arm adoption is allowed. This tests the metric for integrated repeated-structure variation, not covariance evolution, a flow field or identifiable physical nephron lengths. It is a candidate noise model inside the same repeated-object framework.

### Frozen program check: measured results (notebook 22 appendix)

All six fits and all 18 segment panels reproduce their reference positions and all-gene scores exactly in this run. Each panel retains its original 50 genes. scFates neighbor projection beats the privileged quadratic coverage oracle on these selected programs by approximately 14.6%/8.6%/4.8% (S1/S2/S3, Ctrl1A2 reference) and 6.4%/7.5%/3.0% in the reverse direction. This is evidence of transferable molecular organization beyond the tested coverage model; the panels were selected with the earlier atlas and are reused development evidence, not neutral algorithm rankings or true anatomical coordinates.

NB versus scFates is heterogeneous: approximately −1.4%/+1.9%/+4.2% in the first direction and −1.5%/−20.8%/+0.4% in the reverse. The reverse S2 comparison falls below the quadratic oracle by approximately 11.9%. Mean-reference projection also loses that reverse S2 comparison, by approximately 11.9% against scFates. Robust read perturbation alone would have hidden this weakness. No NB replacement or segment-specific winning-arm selection is made. The specified residual-metric ablation addresses a simpler reference likelihood before any transport or covariance-evolution model.

### Residual-metric first run and matched-scale control

The first six real-control covariance fits complete and reproduce the original references. Overall error is approximately 0.070%/0.103% worse than neighbors. In the Ctrl1A4 reference all three Ledoit–Wolf shrinkage estimates equal one, yielding a spherical fitted covariance. The earlier isotropic mean uses nearest-curve residual variance, whereas this estimator uses residuals at supplied training positions; thus their difference includes residual scale as well as geometry. Before interpreting covariance, add a matched spherical control with covariance `trace(fitted_covariance)/p * I`, keeping the same mean, supplied positions, decoder and all query data fixed. This is a diagnostic specified after the first covariance outputs, not a new preregistered discovery. No tuning or arm selection is introduced. Bump notebook logic for the added computed arm and retain the original null result.

### Residual metric and matched control: measured results (notebook 23)

All six updated fits complete with matching saved neighbor/mean coordinates and gene scores; no query fitting is introduced. The global covariance increases all-gene error versus the matched spherical residual control by approximately 0.029% with Ctrl1A2 as reference and has exactly zero mean difference with Ctrl1A4, where shrinkage is one. Both covariance and matched spherical references remain worse overall than neighbors. Geometry therefore earns no benefit in this tested likelihood. A full covariance cannot be credited for changes attributable to residual scale.

On frozen selected programs, covariance versus neighbors is approximately −0.35%/−0.13%/−1.06% (first direction, S1/S2/S3), and +0.32%/−11.46%/−1.03% (reverse). The matched spherical control is approximately −0.25%/+0.91%/−0.97% in the first direction and identical in reverse. The reverse S2 weakness persists below the quadratic oracle. Quarter-read movement remains approximately 0.051/0.072 for the covariance reference, versus 0.043/0.047 for neighbors. No covariance, local covariance or flow is adopted. All 248 regression tests pass; working notebooks retain outputs while staged copies are output-free.

### Next real-data hypothesis: restricted specimen registration

Keep the established source-only scFates backbone and neighbor projector. A shared repeated-object model also needs a restricted specimen intercept; the negative covariance result does not test that intercept. Next compare identity calibration with one constant PC-space offset, fitted only from a fixed structure-ID hash half of the healthy query specimen; evaluate the disjoint half. Query response genes, DPT/scFates coordinates, physical measurements and labels must not enter offset fitting.

For a conservative candidate, estimate the offset only in the orthogonal complement of the source mean-path variation subspace, with rank chosen from source curve singular values to retain 99% of its variation and shrinkage penalty 4 fixed from the existing offset model. Residuals use source-only nearest-curve assignments; perform five fixed update steps. Include an unrestricted constant-offset control with the same calibration rows, penalty and iterations to test whether the restriction matters. If the source curve spans all supplied components, report that the restricted offset is necessarily zero. These are geometric restrictions, not identification of true spatial versus injury directions.

Evaluate original/full calibration and an S3-only subset of those calibration rows as a partial-coverage stress test. Annotation selects that stress-test population only; neither estimator receives its labels. The evaluation half remains identical across calibration regimes. Score all response genes, the existing frozen segment programs and quadratic coverage oracle, and paired quarter-read movement with calibration parameters frozen. Do not force matched occupancy or use population ranks to stretch the query into the full axis. This is an explicit calibration-data allowance, distinct from the fully frozen zero-calibration transfer benchmark. Healthy calibration cannot automatically be used to erase injury or species differences; those branches remain deferred.

## 2026-10-03 — Restricted specimen registration: measured results (notebook 24)

All six source-fold fits reproduce notebook 21's DPT/scFates scores and per-profile coordinates. The evaluation halves contain 1,418/1,419 structures; calibration halves contain 1,364/1,397, with 427/433 used in the S3-only stress test. The same evaluation IDs appear in both regimes. All source mean paths retain four PC-space directions at the fixed 99% threshold; restricted offsets satisfy the orthogonality guard. There are 24 complete calibration records and no failures.

Restricted full-calibration all-gene error reductions versus identity are approximately 0.0006%/0.0015%, and S3-only reductions approximately −0.0004%/+0.0006%. Fine-program gains are small: full calibration approximately +0.013%/+0.040%/+0.049% (Ctrl1A2 reference, S1/S2/S3), and +0.089%/+0.031%/+0.199% in reverse. The unrestricted arm gives larger but heterogeneous effects, including an S2 loss under full first-direction calibration that becomes a gain under S3-only calibration. This does not justify adopting either offset.

Restricted coordinate changes average approximately 0.0005/0.0031 for full calibration. Their quarter-read movements remain approximately 0.0442/0.0464, almost the identity values 0.0440/0.0465; no material observation-robustness improvement is established. The repeated-object prior is preserved, but this particular offset restriction removes little variation relevant to placement. Current canonical neighbor geometry still transfers selected gene programs better than the tested coverage oracle; these remain development-informed molecular tests, not physical longitudinal accuracy. All 253 regression tests pass.

### Goal checkpoint before another experiment

The user explicitly asks whether each experiment advances a method designed for unordered cross-sections of repeated physical structures, rather than a generic trajectory or preprocessing improvement. The prepared reference-rate normalization experiment was never launched. Its unexecuted prototype was removed, and no biological result is claimed. Normalization remains a possible diagnostic only if a measured failure of the shared reconstruction requires it.

The object to estimate remains one canonical ordered PT distribution path from independent specimen point clouds, with unknown nephron memberships integrated out. A useful next experiment must test whether this shared-object inference contributes beyond an established pooled curve. Source-only projection and offset refinements have not yet established that contribution. Negative offset and covariance results must remain visible; they do not justify adding flow or flexible warps.

### Next core hypothesis: groupwise repeated-structure reconstruction

Before any outputs, compare a shared smooth latent mean atlas fitted jointly to source sections and a fixed, unlabeled calibration half of the second healthy specimen. Evaluate only the disjoint query half. Keep the existing source-only marker-excluded count-PC representation and three gene folds; calibration expression does not refit its transform. Supply S1/S2/S3 anchors from the source only. Query calibration labels, response-fold genes, existing query trajectories and depth must not enter reconstruction. Query-label information may select a separately named partial-coverage stress population, never provide its locations.

Compare three existing Gaussian-atlas special cases using the same initialization and hyperparameters: equal-specimen weighting with constant shrunk specimen offsets; equal-specimen weighting without offsets; and observation-weighted pooled fitting without offsets. Keep grid 101, roughness 0.001, offset penalty 4, anatomy strength 1, anatomy width 0.1, tolerance 1e-5 and iteration cap 60 from the existing model. These ablations ask whether identifying specimen point clouds matters beyond merely supplying more profiles. Original nephron identities remain unknown, so the model describes integrated variation, not individual reconstructed nephrons. Constant offsets are restricted in their functional form, but are not identified injury directions.

Fit pooled DPT and nonbranching scFates to the identical source plus calibration rows, root and orient them using source anchors only, and retain their original source-only counterparts as context. This controls the extra unlabeled data allowance. Source labels orient the graph baselines, whereas the atlas also uses their ordinal likelihood; report this supervision difference rather than describe the algorithms as identically supervised. Initialize the joint atlas with the pooled scFates coordinate, then allow the fitted latent positions to change; scFates is a backbone candidate, not anatomical truth. A graph anchor mask and opt-in missing anatomy in the joint model are prerequisites: unknown labels must never be silently replaced with S2 or contribute an ordinal likelihood term.

The primary comparison uses the same frozen neighbor extension for all reconstructed axes and the same source-only response-gene decoder. Calibrating the decoder on query response genes would defeat the transfer test. Report native expression-only atlas projection separately, with the learned calibration-specimen intercept frozen. Retain fixed source-selected positional programs, the privileged quadratic segment-plus-coverage oracle, solver failures, gene-panel disagreement and paired quarter-read perturbation. Depth is descriptive only and excluded from all selection. Conditional posterior concentration is not calibrated anatomical certainty.

The go/no-go question is whether the repeated-specimen construction improves reproducible within-segment molecular prediction or robustness over matched pooled DPT/scFates, in both transfer directions without erasing incomplete coverage. Global all-gene averages and attractive coordinates are insufficient; report heterogeneous and null outcomes. No query-selected parameter sweep, OT, AKI, human or regulator branch follows from a small gain. If the specimen-aware and pooled fits behave alike, simplify around an established curve and state that the repeated-structure contribution has not yet been demonstrated. This protocol is specified, not yet executed.

The groupwise API prerequisites are implemented and reviewed: opt-in unknown
anatomy is excluded from all ordinal likelihood calculations, and graph anchors
control rooting/orientation without excluding unlabeled expression observations.
All 266 regression tests pass, including both strict-default and masked DPT/scFates
smoke checks. The groupwise biological protocol remains unexecuted.

### Notebook 25 execution details fixed before results

The groupwise comparison will include full and S3-only unlabeled calibration on
the same evaluation half. Source anatomical labels alone provide anchors; all
calibration labels passed to fitters are unknown. The primary common-neighbor
extension uses the entire source-plus-calibration point cloud and each method's
fitted coordinates; each response decoder is fitted only on source response
genes at that method's source positions. Native Gaussian projection is reported
separately using the frozen learned query-specimen intercept.

Alongside the quadratic coverage oracle, include a literal privileged segment-only
baseline: fit three response-gene means on source S1/S2/S3 and assign them using
supplied evaluation labels. This is not a label-free projector. Unconverged joint
atlases retain their objective histories and failure records but provide no
accepted biological predictions; successful arms and graph comparators continue
independently. Aggregate only accepted fits and show denominators explicitly.
No solver or parameter changes will be selected from these outputs.

Before examining results, note one remaining model assumption: the Gaussian
atlas likelihood integrates over a uniform latent grid mixture. There is no
query rank stretching or hard occupancy matching, but the fitter also does not
estimate specimen-specific cut-position frequencies. Uniform numerical grid
spacing does not imply uniform physical sampling. Incomplete-calibration behavior
therefore tests this particular likelihood as well as curve and offset fitting;
a failure would not reject the existence of a shared anatomical object.

### Groupwise reconstruction: measured results (notebook 25)

All six source/fold references reproduce notebook 21, and all 12 full/S3-only
pooled graph comparisons complete. Of 36 intended joint-atlas attempts, 21
converge and 15 reach the fixed 60-iteration cap; there are no fitting exceptions.
Full calibration accepts 14/18 arms, versus 7/18 for S3-only calibration. Every
objective history decreases; nonconvergence means the specified stopping
criterion was not met, not numerical divergence. Failed fits remain excluded
from biological prediction summaries. Counts and exact matched-fold comparisons
are saved, so unequal convergence sets cannot masquerade as a weighting effect.

The balanced-offset atlas does not give a consistent improvement over the
source-only scFates reference. On matching accepted full-calibration folds, its
selected-program gains are approximately −3.37%/+8.46%/−9.78% (S1/S2/S3,
Ctrl1A2 reference; two folds) and +0.16%/−4.84%/+0.57% (reverse; three folds).
Its all-gene gains are approximately −0.129%/−0.034%. Paired quarter-read movement
is approximately 40%/68% greater on those same successful folds. Partial
calibration has only one accepted balanced fit per direction and cannot establish
robust transfer. Native Gaussian projection gives mixed program changes and
generally greater read sensitivity; no projector or winning segment arm is adopted.

The apparent reverse-direction gains versus pooled scFates require particular
care. Pooled scFates has weak gene-panel agreement: with Ctrl1A4 as source,
full-calibration S1/S2 pairwise ranks can fall to 0.071/0.167; S3-only calibration
S3 reaches −0.154. These are ordering disagreements, not just gauge stretching.
Source-only scFates is substantially more repeatable in several of those
comparisons, and source/pooled DPT remain more consistent across panels. There
is no true fine anatomical coordinate here, so stability alone cannot choose DPT.
The matched calibration-data allowance reveals that extra profiles do not
automatically improve the reconstruction.

Offsets contribute very little over the balanced no-offset arm (absolute mean
selected-program changes below 0.07% on matched cases). Equal-specimen versus
observation weighting has mixed small effects on the intersection of converged
cases: full first-direction effects approximately −0.13%/+0.19%/−0.75% (one
common fold), and reverse +0.55%/+0.12%/+0.23% (two). Those denominators are too
limited to credit the repeated-specimen prior with a reliable improvement.

There remains useful evidence for continuous molecular organization: source-only
scFates beats the literal three-mean segment oracle on the frozen selected panels
by approximately 16.3%/11.6%/16.3% in the first direction and 6.6%/10.5%/5.4%
in reverse, across all three folds. It also beats the stronger quadratic coverage
oracle on those panels. Selection favored an earlier atlas, shared count exposure
retains aggregate dependence, and two mice have been reused throughout
development. This validates neither true longitudinal distance nor a novel
reconstruction algorithm. Depth remains excluded.

The current unconstrained groupwise Gaussian prototype is not adopted. A useful
shared-object method must preserve repeatable local organization across
references and gene panels; pooling and a converged mean curve are insufficient.
No OT, flexible covariance, disease, human or regulator extension follows.

### Next goal-aligned model check: interpretable regularization units

One concrete limitation is visible in the implemented hierarchy. Both curve and
offset penalties are in absolute count-PC units. The balanced offset update
shrinks by `1 + 4*sigma2`, and observed residual variances span approximately
10.5–20.0, giving shrinkage factors around 43–81. Fitted offset norms are only
0.012–0.044. The null offset result therefore describes this strong prior; it
does not establish absence of specimen variation.

Analytically, rescaling the supplied PCs, mean and offsets by c and residual
variance by c squared leaves likelihood responsibilities unchanged, but multiplies
these fixed quadratic penalties by c squared. An arbitrary representation scale
therefore changes the effective reconstruction prior. Before another biological
experiment, check the optimizer/objective's behavior under a uniform unit change
and specify a training-defined reference-unit convention. This is a requirement
for a reproducible shared-object model across gene panels, not a new expression
normalizer or a query-tuned penalty sweep. Any corrected convention must return
to the same matched real-control benchmark and still earn improvement in fine
programs and reproducibility. Uniform latent mixing and absent fine anatomical
ground truth remain separate limitations.

### Reference-unit prerequisite: actual-data and regression evidence

A fixed probe uses Ctrl1A2, fold 0, full unlabeled Ctrl1A4 calibration and the
original pooled-scFates initialization. With all existing raw-unit parameters
unchanged, multiplying PCs by seven changes average fitted positions by 0.0360
(maximum 0.5266); accepted fits stop after 16 versus 7 iterations. This changes
the effective dimensional priors, not the underlying point geometry. It is a
numerical model check, not a biological prediction gain.

An explicit `prior_scale` now computes the model in declared reference units and
restores molecular parameters/variances to input units. History objectives and
stopping criteria remain in reference units. The previous default scale 1 is
preserved. Scaling both data and reference unit by seven reproduces the actual
probe's positions to maximum approximately 1.3e-9, objectives to approximately
2.7e-10, and the same 16 iterations. Code regression fixtures cover parameter
restoration, frozen query projection, strict defaults and invalid/unrepresentable
scales; they are not synthetic biological experiments. Changing reference scale
for fixed data intentionally changes prior strength.

### Next actual-data test specified before outputs: source-relative priors

Notebook 26 retains notebook 25's exact data, gene folds, source-only transform,
hidden calibration anatomy, source-plus-calibration graph inputs, initialization,
three model arms, full/S3-only regimes, source-only response decoders, frozen
programs, privileged segment/coverage comparators and paired read perturbation.
Set the prior reference scalar to the root mean square of centered source PC
values across rows and all 15 components. This uses only source modeling features
and uniformly converts units for the latent fitter; it neither changes the PCA
encoder nor reweights graph geometry. The priors now refer to total source
molecular variation, not to an estimated measurement-noise standard deviation.
No query-based scale, penalty, variance estimator or iteration sweep is allowed.

Reproduce both source-only and pooled graph positions/scores from existing
notebooks. Compare new joint arms against the old raw-unit arms on identical
accepted cases, showing missing fits and successful-fold counts. Require gains
against the source-only references as well as the extra-data pooled comparators;
report molecular program, read and gene-panel stability separately. The object
remains one ordered physical structure with integrated unknown-nephron variation
and restricted constant specimen offsets. This fixes interpretable prior units;
it does not identify anatomy, fix uniform sampling assumptions or automatically
validate the hierarchy. The estimator remains unadopted until the actual results
earn its use.

### Actual-control reference-unit result (notebook 26)

Executed `26.source_relative_priors.1`, cache key `dc5f39df312f7c54`.
All six source references and twelve pooled-graph regimes reproduce notebook 21
and 25 scores, evaluation coordinates and frozen quarter-read coordinates.
Source-only reference scalars range 6.329–7.898. All 36 joint fits return without
exceptions, but only **5/36 converge**, versus **21/36** under raw-unit priors.
All 24 specimen-balanced fits fail the fixed convergence criterion. Accepted
fits are observation-weighted without offsets: Ctrl1A2 full 2/3, S3-only 2/3;
Ctrl1A4 full 0/3, S3-only 1/3. Every objective history descends. Nonconvergence
at the fixed cap is not evidence of divergence or rejection of the physical
shared-object prior. No cap, scale or penalty tuning follows these outcomes.

Only two source/fold/regime cases are accepted under both prior conventions.
For their common-neighbor projection, selected-program gains versus raw-unit
fits are −3.84%/−8.91%/−2.93% (S1/S2/S3, Ctrl1A2 full) and
+1.23%/+3.41%/−2.90% (Ctrl1A4 S3-only), each one fold. Aggregate all-gene
gains are only +0.049% and +0.040%. Matched read movement improves overall
by 4.65% and 42.38%, respectively, while S2 movement worsens 80.77% in the
first case. These are separate descriptive metrics on sparse survivor sets.

Against source-only scFates, accepted observation-weighted common-neighbor
program gains are −5.72%/+1.13%/−1.26% (Ctrl1A2 full, two folds),
−0.90%/+1.90%/−0.73% (Ctrl1A2 S3-only, two) and
+1.63%/−11.81%/−0.09% (Ctrl1A4 S3-only, one). Some improvements over
DPT remain, but the stronger scFates reference and opposite direction do not
support adoption. Full Ctrl1A2 read movement is 32.02% greater than source-only
scFates on matching accepted folds. Native projection is reported separately.
Both plots were inspected; source-only context remains visible so gains over a
weak pooled comparator cannot stand in for a method-level improvement.

The reference-unit API passes code regression checks, including joint unit
rescaling and restoration; the old default reproduces HEAD exactly. The full
synthetic code suite passes 280 tests. The API remains explicit and default 1;
the source-RMS convention is experimental, not a silent public-default change.

### Goal decision after the unit check

The unit convention is now explicit, but relaxing the dimensional priors does
not rescue this unconstrained groupwise Gaussian estimator. Stop expanding this
branch with flow, covariance flexibility, specimen warps or outcome-selected
optimizer settings. The unresolved task is reconstructing reproducible local
organization from repeated cuts while permitting unequal coverage and molecular
variation. A next model experiment needs a concrete restriction or observation
model addressing that task, and must justify itself against the completed
source-only and pooled benchmarks before it is run. The project goal is still
unmet; a numerical property and a reusable negative benchmark are progress, not
a validated pseudostructure method. Depth remains unused.

## 2026-10-04 — Assignment-aware prediction: goal check before outputs

The next check is a statistical contract of repeated-cut reconstruction, not
another Gaussian hierarchy. Unknown section positions may have several plausible
assignments on the canonical object. Notebooks 22/23 keep Gaussian/NB assignment
weights but score responses at `decoder(E[Z|x])`; kNN similarly collapses neighbor
assignments to their mean coordinate. For nonlinear response curves this differs
from `E[decoder(Z)|x]`. The earlier negative scores concern point-coordinate
prediction, not the full probabilistic reference prediction. Uncertainty
propagation has precedent in Campbell and Yau's 2016 probabilistic pseudotime
work; no mathematical novelty is claimed here.

Notebook 27 will freeze notebook 22's source-only scFates axis, count transform,
mean/NB projectors, parameters and response decoder. Add DPT's existing source
reference as required context. Compare point-coordinate and assignment-averaged
decoding for **every** projector, including both graph baselines using their
actual inverse-distance neighbor weights. Repeat both control directions, three
gene folds, all-gene and frozen within-segment panels, privileged segment and
quadratic coverage predictors, and the same paired quarter-read realization.
Reproduce earlier coordinates and point-response scores before new comparisons;
no fitter, coordinate, query label, prior, parameter or sampling distribution is
changed. All-gene exposure and reused source-selected panels retain the earlier
independence limitations. kNN weights are assignments, not calibrated posteriors.

This addition is rejected as a reason to adopt probabilistic projection if it
fails to improve the relevant held-out programs beyond assignment-aware graph
comparators. Any prediction improvement is a response-reduction effect: it cannot
prove better spatial ordering, calibrated anatomical uncertainty or the final
joint estimator. Even a positive result must return to the unresolved shared
local-order problem. Do not revive the unconstrained hierarchy, flow, AKI or TF
branches from this diagnostic. Depth remains excluded.

### Assignment-aware prediction: measured result (notebook 27)

All six real-control transfers complete and pass coordinate, point-score,
selected-program and paired quarter-read reproduction guards. Cache key
`4ab201f39d67f090`. The source axes, fitted projectors and every spatial ordering
remain unchanged. The response reduction is now available for shared-grid
posteriors and per-query neighbor assignments without allocating a
query-by-support-by-gene tensor. Default neighbor projection remains unchanged;
293 synthetic code regression tests pass.

scFates assignment averaging improves all-gene errors over prediction at its mean
coordinate by +0.143%/+0.210% in the two directions. Selected-program gains are
+0.34%/+7.50%/+3.20% (S1/S2/S3, Ctrl1A2 source) and
+4.52%/+4.25%/+3.42% (reverse), each with all three folds. DPT's aggregate
gains are +0.043%/+0.059%, with smaller and heterogeneous program changes.
These are gains from averaging nonlinear response curves, not improved positions.

Gaussian/NB projection is not rescued. Versus assignment-aware scFates, the mean
reference changes all-gene errors by −0.116%/−0.157%, and NB by
−0.388%/−0.645%. NB selected-program gains are −1.75%/−5.97%/+0.91%
and −5.94%/−25.74%/−2.72%; mean-reference gains are
−1.05%/−6.37%/−1.89% and −1.19%/−13.57%/−2.52%. Each uses three
folds. Thinned-query S1/S3 comparisons sometimes favor NB, but reverse S2
remains worse. Strong conditional concentration and read stability do not
establish correct local assignments.

Assignment-aware scFates improves the frozen selected programs over the
privileged quadratic coverage oracle by +14.85%/+15.44%/+7.81% and
+10.61%/+11.37%/+6.31%. This strengthens a development comparator; the
source-selected panels and shared all-gene exposure still preclude independent
anatomical claims. The program plot was inspected, and cropped labels corrected.
Future comparisons should explicitly name point versus assignment-aware response
reduction for every arm. Retain the historical point results, and do not adopt
another projector or revive the unconstrained Gaussian hierarchy.

### Next method hypothesis: anatomy in canonical curve construction

The next method experiment must change reconstruction, not prediction again.
Our current scFates baseline fits a generic elastic curve and uses source anatomy
only for root and orientation. A concrete physical-object prior is to initialize
a nonbranching curve with a polyline through the source modeling-feature means
of S1, S2 and S3, in that order. No marker-expression profiles or fixed segment
widths are supplied. The installed ElPiGraph/scFates backend already accepts
initial node positions and edges; its author documentation and primary paper
give precedent for changing initialization. This is a tailored initialization,
not invented curve mathematics or an identification theorem.

Before running this experiment, fix a matched source-only and pooled-calibration
protocol using the existing actual controls, gene folds and partial-coverage
evaluation halves. Initialization must use source anchors only; pooled query
calibration anatomy remains hidden. Keep existing graph parameters, root,
orientation, response decoder and neighbor weights fixed. Require old-baseline
reproduction and report both response reductions, DPT/scFates comparisons,
within-segment program prediction, paired read perturbation and gene-panel rank
agreement. Reject the addition if it does not improve reproducible local order
and prediction beyond the established references. Distinguish initialization
from enforced anatomical order: subsequent elastic fitting can still move the
nodes. Unknown nephron variation remains integrated; no occupancy matching,
specimen warps, covariance or flow is introduced. A successful three-anchor
result must later face the weaker endpoint-only ablation. No new fit has been run
for this hypothesis yet.

### Notebook 28 protocol fixed before actual-data outputs

Use the two healthy controls, three SHA256 gene folds, 15 source-only count PCs
(dispersion 100), and the same hash calibration/evaluation halves as notebook 25.
Full and S3-only calibration retain identical evaluation rows. Source-only and
pooled default scFates are matched to source-only and pooled anatomy-initialized
curves. DPT source/pooled context remains. Source curve initial means use only
known source S1/S2/S3 rows; every pooled query calibration label is −1 and
excluded from initialization, rooting and orientation. Query S3 labels only
select the stress population and enter stratified scoring or privileged oracles.

Keep 30 curve nodes, elastic lambda 0.01 and mu 0.1, seed 15, root/orientation
rules, backend defaults and 15-neighbor extensions fixed. The ordered initial
graph has three source-centroid nodes and edges 0–1–2; the fitted graph is still
free to move. Require a connected nonbranching curve, finite coordinates, and
report source segment medians without filtering out anatomical inversions.
No specimen weights, ordinal likelihood, offset, covariance, flow or sampling
prior is added. This is a mean-atlas special case for repeated physical objects;
individual nephron identities and physical longitudinal distances are unknown.

Source fits reproduce notebook 27's full-query point/assignment predictions,
coordinates and paired thinning. Pooled/default half-query paths reproduce
notebook 25's coordinates, read coordinates and point scores. Then score every
arm with both response reductions, source-only response splines (grid 101,
ridge 0.001), frozen programs, literal segment and quadratic coverage oracles.
Retain failures and attempt counts for 18 new anatomical curve fits (six source,
twelve pooled) and their default/DPT comparators. Report matched program gains,
paired read changes, calibration-regime sensitivity and within-segment gene-panel
rank agreement. No query-based initialization, parameter or arm selection is
allowed. Repeated development panels and shared exposure remain limitations;
gene-panel rank agreement is stability, not true anatomical order. This test
can reject this initialization without rejecting the ordered-object prior.

### Notebook 28 actual-control readout (2026-10-04)

All 54 graph attempts completed: six source and twelve pooled fits each for DPT,
default scFates and anatomy-initialized curves. All eighteen new curves were
valid nonbranching graphs with ordered source segment medians. Uncached guards
reproduced notebook 27's full-query point/assignment scores, selected programs,
IDs and original/thinned coordinates, and notebook 25's evaluation-half default
scores, programs and coordinates. Cache key: `05066900892d77de`. The program
figure was inspected; no fits, inversions or unfavorable panels were omitted.

Source-only anatomy initialization changes assignment-aware all-gene prediction
error by −0.014%/+0.003% relative to generic source scFates (Ctrl1A2/Ctrl1A4
training directions). Its S1/S2/S3 selected-program gains are
+0.12%/+0.86%/−0.74% and +0.14%/+0.06%/−0.08%, each with three folds.
Source predictions are identical between calibration regimes by construction;
those duplicates are not extra replication.

Pooled anatomy versus source-only assignment-aware scFates has full-calibration
program gains +0.55%/+1.79%/−1.23% and +0.97%/−10.02%/−0.07%.
With S3-only calibration these are −0.76%/+0.19%/−0.02% and
−0.47%/−5.23%/+1.11%. Relative to its own matched pooled generic scFates,
the new initialization worsens ten of twelve program cells (two directions ×
two regimes × three segments); gains range −1.35% to +0.46%. Each cell has all
three folds. Improvements over a weaker pooled fit do not establish useful
canonical reconstruction relative to a frozen source reference.

Some stability changes are favorable, but not consistently. Source quarter-read
coordinate MAE changes 0.0440→0.0384 and 0.0465→0.0478. Source within-S3
gene-fold rank agreement changes 0.943→0.818 and remains approximately 0.794
in the reverse direction. Full-calibration reverse pooled S2 rank agreement
improves 0.447→0.924, despite worse held-out S2 program prediction. Partial
reverse pooled S3 agreement instead falls 0.235→0.089. These rank summaries
average three gene-fold pairs and measure stability, not true longitudinal order.
DPT's within-segment rank agreement is at least 0.915 across these cells;
this is another stability comparator, not proof of anatomical accuracy.

Pooled calibration-regime coordinate MAE improves 0.0631→0.0597 and
0.1319→0.1016 versus pooled generic scFates, but remains considerable. Query
read perturbation also does not uniformly favor the new curves. Shared all-gene
exposure, source-selected development programs, reused controls and only two
biological specimens remain limitations. Source anatomy derives from upstream
marker annotation; source-centroid initialization is not independent anatomical
validation. All labels used for query calibration are hidden from curve fitting.

**Decision:** do not adopt this initialization as the reconstruction method and
do not expand it into endpoint-only, parameter-search, covariance, flow, AKI or TF
branches. Ordered initial nodes and ordered segment medians are insufficient to
protect useful fine reconstruction against specimen and sampling effects. This
rejects an initialization addition, not the repeated-physical-object prior.
Before the next experiment, revisit how that prior enters the estimator during
fitting and how canonical position can be kept separate from sampling density.
A new hypothesis must identify a reconstruction failure it changes and state a
matched measurable retention criterion before fitting; improving generic curve
appearance or prediction reduction alone will not satisfy the project goal.
