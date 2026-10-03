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
