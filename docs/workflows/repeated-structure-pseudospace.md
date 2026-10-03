# Notebook 15: healthy canonical atlas first

The canonical notebook is
`analysis/notebooks/15_spatial_probability_flow_proof_of_concept.ipynb`.
It tests healthy cross-specimen molecular position mapping independently of
probability-flow success. The earlier OT-distilled flow result remains a historical
transport diagnostic. There is no notebook 15 mirror.

## Hypotheses and model hierarchy

The primary hypothesis is that one healthy PT realization teaches a canonical
mean molecular atlas that locates profiles from another realization using
expression alone. The secondary hypothesis is that spatial covariance evolution
adds held-out distributional predictive value beyond the same mean path.

Fit training-only PCA and a kernel-smoothed mean curve on training DPT. Project
held-out profiles by nearest molecular mean, without their DPT or S1/S2/S3 labels.
Compare with full-covariance Gaussian likelihood mapping on the same grid and mean
estimator. The Gaussian path also defines an affine probability-flow field through
the symmetric Lyapunov equation. Its flow and density atlas share the same
marginals: an exact marginal-preserving field does not independently change the
position likelihood. Therefore coordinate mapping compares mean versus Gaussian
atlas; transport compares mean motion versus mean-plus-covariance evolution.

The training axis is still upstream DPT. This is an inductive reference-atlas test,
not a wholly DPT-free joint reconstruction of training locations. Two controls
support descriptive comparisons only. No true longitudinal spatial metric or
biological dynamical mechanism is identified.

## Inputs and execution

Use `kidney-pseudospace`, run cells top to bottom, and provide roots through
`--data-root` / `--results-root` or `PSEUDOSPACE_DATA_ROOT` /
`PSEUDOSPACE_RESULTS_ROOT` (flags take precedence). Only Ctrl1A2 and Ctrl1A4 count
matrices under `tubule_by_gene/` and their matching
`*_kept_tubules_labeled_fine.geojson` segmentations are required. AKI is deferred.
The loader validates IDs and centroids, applies the pooled healthy-cohort mouse p5
spot-count structural floor, retains all common measured genes, normalizes total
counts to 10,000 + log1p before PT selection, and performs no tubule gene-count QC.

The default training-coordinate reference is
`<results-root>/human_vs_healthy_mouse/cross_species_pt_dpt.h5ad`, column
`total_scanpy_dpt`. Override via `PSEUDOSPACE_15_DPT_REFERENCE` and
`PSEUDOSPACE_15_DPT_COLUMN`. A strict unique within-specimen centroid join within
one pixel links current validated PT to the export; feature indices must agree
when available. DPT availability is recorded, but profiles absent from it remain
held-out projection targets. Only DPT-available training profiles fit the atlas.
Training DPT is scaled using that specimen alone, not the other control.

The historical fixed-axis transport diagnostic still uses a common supplied DPT
window grid; this is explicitly distinct from expression-only projection. Both
human-labeled sections in 03 are healthy cortex; neither contributes molecular
profiles to notebook 15's fits.

Physical validation uses the reviewed glomerular pass-1 export, default
`<results-root>/human_vs_healthy_mouse/cross_species_harmony_pass1.h5ad`, overridden
with `PSEUDOSPACE_15_GLOMERULI`. Anchor centroids/IDs are validated against current
segmentations. Nearest and mean-three-glomeruli distances enter only after
projection. They are cortical-depth proxies, not longitudinal nephron ground truth.

## Primary atlas comparison

- Exchange training and held-out specimens and rotate three deterministic gene
  folds. Canonical PT markers and evaluation genes are excluded before training
  scaler/PCA (15 PCs). All training observations enter kernel curve estimation;
  effective kernel support is checked at every point on a 101-point grid.
- Fit the mean atlas with fixed bandwidth 0.08. Use Euclidean nearest-mean MAP
  projection as primary; retain the isotropic model posterior, entropy and residual.
  Fit a Gaussian atlas with identical kernel means, full covariance, 25% diagonal
  shrinkage and a global variance floor, and retain its posterior separately.
- Estimate an optional single global specimen intercept from a deterministic half
  of held-out modeling-feature profiles. Penalize it by factor 4 and freeze it
  before projecting the disjoint validation half. Unadjusted projection is primary.
  This is transductive nuisance calibration, not zero-shot mapping or an evaluation
  gene offset. Sampling composition and intercept/location confounding remain.
- Fit held-out-gene mean curves on the training coordinate and interpolate at
  projected positions; scale MSE by training gene variance. Compare to literal
  training S1/S2/S3 gene means, using test labels only for this explicit oracle
  anatomy baseline. Also compare supplied DPT on a common available-profile subset.
  Report all-segment and within-segment errors, all test profiles, and the disjoint
  offset-validation subset separately. No test centering or hyperparameter tuning.
- After projection evaluate withheld segment rank/order, segment median positions,
  physical depth, excluded marker curves, and training-versus-held-out gene curve
  agreement. Label-derived anatomy and inherited training DPT retain indirect
  marker information; gene holdouts retain compositional normalization dependence.

## Gaussian covariance-flow ablation

Compare identity, raw training centroid increments, smooth mean-only increments,
and Gaussian affine flow on identical held-out source/target clouds. Use overlapping
and non-overlapping windows in both specimen directions. The field uses a cubic
mean interpolation and piecewise-linear positive-definite covariance interpolation;
the symmetric Lyapunov solution matches covariance derivatives within each interval.
RK4 integrates the field. Score energy distance, mean error, covariance Frobenius
error and centered-cloud energy. A smooth Gaussian atlas can reproduce an estimated
marginal path without establishing a unique particle mechanism.

Flow failure never gates mean-atlas fitting or validation. A descriptive flow screen
requires better energy distance than the same smooth mean path in both directions
and both window schemes. A separate gene-prediction screen asks whether the mean
atlas beats the strong segment oracle in every gene fold and direction, and whether
Gaussian likelihood beats the mean atlas. These screens do not substitute for
scientific review of depth, ordering, uncertainty and stability.

## Stability, outputs and reproducibility

Bootstrap training tubules, modeling genes, and both; refit PCA and the mean atlas
and reproject the untouched test profiles without test coordinates, labels or
intercept fitting. Sparse-support failures remain explicit. Posterior uncertainty
and bootstrap position variation are distinct conditional quantities, not calibrated
spatial credible intervals. Fully latent training-coordinate refinement is deferred
until this simpler atlas test establishes its value.

Outputs live under
`<results-root>/spatial_probability_flow_proof_of_concept/`: per-fold gene errors,
positions/full posteriors, nuisance offsets, failed fits, segment/depth validation,
marker/gene curves, distribution transport ablations, bootstrap positions, diagnostic
PDFs and decision JSONs. Expensive stages cache by actual inputs, parameters and
implementation under `stage_cache/`. `PSEUDOSPACE_STAGE_CACHE=0` rebuilds; validation
is never cached. `PSEUDOSPACE_15_BOOTSTRAPS` defaults to 10 and must be at least two.

Reusable implementations are `pseudospace/canonical_mean_atlas.py`,
`pseudospace/gaussian_spatial_flow.py`, and the historical
`pseudospace/spatial_probability_flow.py`; measured loading is unchanged in
`pseudospace/repeated_inputs.py`. Tests use only synthetic fixtures.

AKI projection and TF Jacobian interpretation are deferred in this version.
Do not infer superiority, correct injured location, regulator activity or causality
from a smooth curve. Stage output-free notebook copies without clearing live
outputs, run the notebook stage-order and staged hygiene gates, and keep private
matrices, segmentations and result artifacts out of Git.
