# Notebook 15: repeated-structure reconstruction

The canonical notebook is
`analysis/notebooks/15_spatial_probability_flow_proof_of_concept.ipynb`.
It tests healthy cross-specimen molecular position mapping independently of
probability-flow success. The earlier OT-distilled flow result remains a historical
transport diagnostic. There is no notebook 15 mirror.

## Repeated-structure prior and umbrella model

Every observed section is a partial observation of a physical realization of the
same ordered anatomical object. Let $r$ index physical nephron, $s$ specimen,
$d$ an observed species/condition domain, and $z$ canonical position:

$$
x_{d,s,i}=X_{r_i,s,d}(z_{d,s,i}), \qquad
X_{r,s,d}(z)=\mu(z)+\alpha_d(z)+b_s(z)+u_{r,s,d}(z)+\epsilon.
$$

The common program is $\mu$; domain biology is $\alpha_d$; restricted specimen
variation is $b_s$; realization-level variation is $u$. Unknown nephron identities
and absent linked longitudinal profiles prevent fitting individual nephron
functions or their across-position covariance. Integrating that variation out
produces the canonical distributional object:

$$
x_{d,s,i}\mid z_{d,s,i}\sim\rho_z^{d,s}, \qquad
\rho_z^{d,s}=(T_{d,s,z})_\#\rho_z.
$$

$T$ is a restricted domain/specimen deviation model. A translation by
$\alpha_d(z)+b_s$ is one special case; position-dependent covariance is another
model component. The prototype does not currently fit these hierarchical
transformations jointly or shrink specimen distributions toward a pooled latent
canonical distribution. A shared ordered coordinate is the prior; identical
molecular states, equal section sampling density, and complete axis coverage
are not required.

Healthy mouse, AKI mouse and healthy human are molecular realizations of homologous
structure under this framework. Domain-specific biology can differ at the same
position without defining a separate disease/species trajectory. Both human-labeled
slices are healthy cortex from one donor; they are not independent human donors.
The available domains do not support extrapolating human AKI or identifying an
unobserved species-by-injury interaction.

The hierarchy is **prior → inference → variation → discovery**: posit a shared
ordered physical object; reconstruct position from repeated cross-sections; allow
restricted specimen and domain deviations; then investigate conserved programs,
injury deviations, species differences and candidate regulators. This is a
problem formulation and estimator design, not an established novelty claim.

## What the current prototype estimates

Notebook 15 models the PT substructure; its numerical [0,1] gauge is not a
whole-nephron coordinate. It fixes domain to healthy mouse and fits a mean or Gaussian reference
from one specimen's supplied DPT coordinate. Unknown realization variation is
absorbed into residual density. The unadjusted reference sets specimen offset to
zero; calibration adds the special case $b_s(z)=b_s$. Gaussian covariance also
contains measurement and finite-window positional variation, so it cannot be
interpreted as an estimated nephron-specific random process.

The exchanged training/test controls are a validation experiment for shared
structure, not the final estimator. The intended full model jointly estimates a
specimen-balanced canonical object with restricted nuisance/deviation terms.
A joint fit on both controls would use both specimens for estimation; the same
specimens cannot then be advertised as an independent held-out assessment of
that joint estimator. Future joint latent fitting needs its own validation.

## Identification and replication limits

A healthy-mouse reference can define $\alpha_{healthy\ mouse}=0$, with an explicitly
chosen coordinate gauge and restricted/centered specimen deviations. Such
constraints distinguish the canonical object from nuisance terms by convention;
they do not prove that biological structural position is uniquely identified.
Orientation and monotone spatial reparameterization remain ambiguous without
anchors, and the inferred coordinate is not calibrated physical length.

Domain expression changes parallel to the canonical mean-path tangent can mimic
position shifts. Injury is not guaranteed to be geometrically orthogonal to the
axis. Freezing a reference or assigning a high residual density score does not
establish correct injured/human location. Restricted deviations, common measured
features, appropriate anatomical anchors, and independent validation must support
that separation; unrestricted domain functions or spatial warps could absorb
any alignment. Do not force partially sampled domains to occupy the entire axis.
Cross-species modeling must retain the existing measured-gene/ortholog-availability
contract and cannot treat structural-zero columns as expression observations.

Unknown nephron membership can induce dependence between sections. Treating
realizations as exchangeable is a modeling assumption, not proof that every
section is from a different nephron. Specimen/donor remains the biological
replication unit; thousands of sections do not create additional independent
mice or human donors. Current atlas transfer supports a limited reproducible
reference mapping, not estimation of all hierarchy levels or separation of
structural and disease/species variation.

## Hypotheses and model hierarchy

The primary hypothesis is that one healthy specimen containing many PT realizations
teaches a canonical mean molecular atlas that locates profiles from another realization using
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

## Notebook 16: jointly latent repeated-structure atlas

`analysis/notebooks/16_joint_repeated_structure_atlas.ipynb` develops the next
estimator independently of notebook 15's DPT reference. Its shared B-spline mean,
isotropic residual distribution, constant shrunken specimen offsets and learned
smooth ordinal probabilities are fitted using `pseudospace/joint_repeated_atlas.py`.
Equal specimen mass is used in both the modeling representation and latent
likelihood. Unknown physical nephron identities are marginalized; the objective
is a descriptive specimen-balanced marginal fit, not evidence of thousands of
independent biological replicates. The discrete uniform coordinate prior is a
computational gauge, not physical distance or a claim of uniform observed sampling.

Two protocols are kept separate: joint registration of 80% of profiles in both
controls followed by expression-only projection of the remaining sections, and
one-control training followed by frozen projection of the other. Three withheld
gene folds are rotated. The reconstruction never receives DPT, x/y or glomerular
depth. Training segment labels are weak supervision with learned transition
cutpoints; their marker-derived provenance remains explicit. Query labels are
consumed only by the reported segment-oracle predictor and evaluation.

Fixed profile/gene splits, same-reference segment means, random initialization
checks and anatomy/offset ablations test whether continuous positions earn their
complexity. Physical depth and coverage diagnostics are evaluated after inference;
DPT comparisons use matching eligible query profiles. Posterior assignments are
conditional model quantities, not calibrated anatomical confidence intervals.

Results are written to `<results-root>/joint_repeated_structure_atlas/`, including
input manifest, coordinates, held-out gene scores, objective histories, failures,
initialization stability, matched physical validation and coverage diagnostics.
Expensive fits are content-addressed caches; validation cells always run. The
notebook has no generated mirror. Follow the same output-free staging convention
as notebook 15. Decisions and sources accumulate in
`docs/results/pseudostructure-research-ledger.md` rather than replacing earlier
measured evidence.

The notebook now also retains an exposure-aware representation pilot:
`pseudospace/count_representation.py` learns reference-only, specimen-balanced gene
probabilities and frozen NB Pearson residual PCs from raw counts. Total all-gene
library exposure is supplied as observation normalization, preserving the disclosed
compositional dependence of gene holdouts. This is an approximate representation,
not a full count likelihood. Rotated gene folds and an endpoint-orientation-only
ablation preserve the original validation boundaries.

Same-section binomial thinning directly tests read-depth robustness with frozen
references; it changes observed counts, not physical section identity. Strong PC1
and soft segment-mixture comparators use `pseudospace/repeated_atlas_baselines.py`.
Keep depth associations subordinate to these molecular/stability/measurement checks:
winding anatomy and unknown nephron identity make cortical distance an especially
imperfect supporting proxy. Stage-specific logic versions retain unchanged caches;
new/changed stages have their own bumped notebook version and full input/code keys.

## Notebook 17: conditional-count projection pilot

`analysis/notebooks/17_exposure_conditioned_count_atlas.ipynb` tests a frozen gene-rate
atlas from the training-only count-Gaussian coordinate. It requires no DPT or physical
proxy exports. `pseudospace/count_rate_atlas.py` estimates exposure-normalized rates
separately in each original specimen and averages them equally at every position;
projection evaluates a full NB count likelihood conditional on query library size.
This is a projection/decoder upgrade, not a fully joint latent-count estimator.
Independent genes and fixed dispersion are approximations; posterior precision is
not calibrated anatomical uncertainty.

The notebook repeats fold-0 joint/profile and independent-specimen validation, strong
PC1/soft-segment/oracle comparisons, and same-section thinning. Outputs are under
`<results-root>/count_rate_projection_pilot/`, with the same content-addressed fit cache
and uncached input/guard validation. It has no generated mirror. Inference does not
consume private outputs of notebook 16: the reference is recomputed from validated
inputs and fixed settings. Keep the working notebook outputs and stage a cleaned copy.

A separate `17.count_rate_dispersion_sensitivity.1` stage rotates all three gene
folds for joint profile holdout and scores the frozen rate path at NB shape 10,
100 and 1000 plus the Poisson limit. It reuses each fold's identical all-gene
quarter-read thinning draw across observation settings, keeps PC1/soft-segment
baselines, and exports within-segment dispersion and gene-fold coordinate agreement.
Completeness guards run outside the cache; partial failures are not a completed
benchmark. These diagnostics do not choose a setting from query outcomes or depth.

The sensitivity report includes all gene-fold pair plots and per-section coordinate
ranges beside conditional entropy. Large panel-dependent shifts are retained in
benchmarks and reported as ambiguity; a sharp likelihood posterior alone is not
a trustworthy location confidence score. Private CSVs retain section IDs; committed
notebooks remain output-free.

## Notebook 18: continuous programs within coarse PT segments

`analysis/notebooks/18_within_segment_repeated_structure_validation.ipynb`
validates the count-Gaussian mean atlas on both independent-control transfer
directions and all three gene folds. Its `18.within_segment.1` cache is separate
from earlier workflows, and outputs reside in
`<results-root>/within_segment_repeated_structure_validation/`. No DPT, depth,
AKI or human exports are required, and no mirror is generated.

`pseudospace/within_segment_validation.py` provides an exact specimen×segment
coordinate-shuffle null and a privileged segment-plus-covariate oracle. The
null retains coarse coordinate distributions but removes fine training response
ordering. The oracle explicitly receives query segment labels and coverage
metadata; these labels never enter query atlas or soft-segment inference.

All withheld genes remain the primary prediction endpoint. Training-only
two-way profile cross-fitting selects at most 50 positive-gain genes per segment
and fold for supplementary transfer tests. Shared all-gene exposure preserves
compositional dependence. Profile splits, genes, null draws and folds do not
create independent biological replication. Inspect all-gene and selected-program
results together and retain the ambiguity limitations from notebook 17.

## Frozen pseudostructure atlas API

`pseudospace/pseudostructure_atlas.py` exposes a wrapper for fitting and
projecting a frozen repeated-structure reference. It encodes repeated physical
objects as a shared ordered-position prior while integrating unknown nephron
identity out; it does not recover linked nephron trajectories. The caller must
provide an explicit boolean `model_gene_mask`. Only those genes enter the count
residual representation and rate atlas. `library` remains the raw exposure summed
over **all measured genes**, including excluded genes, so this contract retains
the documented compositional dependence. Sparse count matrices are densified
inside the wrapper; use filtered structure aggregates whose dimensions fit memory.

Fit on training structures with `fit_pseudostructure_atlas(counts, library,
specimens, anatomy, gene_names, model_gene_mask, ...)`. The three ordered anatomy
codes are weak supervision during training only (S1=0, S2=1, S3=2); marker-derived
label provenance remains. The fitted Gaussian mean atlas and NB rate atlas are
frozen. Project new structures with `project_pseudostructure_atlas(reference,
counts, library, gene_names, structure_ids, specimens=None)`: query counts must
have exactly the same measured gene names in exactly the same order as training,
even when some genes were excluded from modeling. Projection does not accept
query segment labels or refit either model. Supplying query specimen names only
selects a learned constant training offset where available; an unseen specimen
uses zero offset.

Gaussian and conditional-count outputs include model-conditional position
posteriors and entropy, not calibrated confidence intervals. Agreement or
disagreement between the two frozen decoders is descriptive. For sensitivity to
feature choice, fit multiple references with aligned profile order and apply
`summarize_panel_sensitivity`; its across-panel coordinate range is a sensitivity
summary, not an uncertainty interval. Changes to an unmodeled domain, injury
state, or tissue composition remain unresolved by this wrapper and require
separate validation rather than being treated as known spatial displacement.

## Identification stress tests (notebook 19)

`19_repeated_structure_identification_stress.ipynb` generates unordered, independent
count observations with known simulated position, variable library exposure and
an orthogonal nuisance program. It freezes the current public estimator and
compares full/partial query coverage, paired read thinning and reversed nuisance
coupling. It preserves failed convergence cases and diagnoses their original
solver history; they are not accepted projections. An explicit same-profile
subset guard checks that changing query population composition cannot rescale
coordinates. The notebook is synthetic-only and needs no private input data.

The independent-nuisance scenario exposes the current PC1 initialization's
failure; the soft segment baseline still recovers simulated ordering. Reversed
nuisance coupling yields wrong ordering even when the two model estimates agree.
These are identification limits to address before disease/species mapping, not
additional biological validation. See the research ledger for the fixed design,
measured results and the next training-only initialization ablation.

## Training-only initialization ablation (notebook 20)

`20_repeated_structure_initialization_ablation.ipynb` holds the model objective
and convergence settings fixed while comparing PC1, specimen-balanced S3-minus-S1
endpoint contrast and a soft predicted ordinal code. The reusable
`initialize_atlas_positions` helper uses training expression, original specimen
IDs and coarse labels only; it does not accept true fine positions or query data.
The public estimator retains its original default. The notebook exposes failed
solver histories and never projects a nonconverged fit as an accepted atlas.

Endpoint contrast makes the tested synthetic references converge, but the
independent-nuisance reconstruction remains weaker than the soft segment
baseline. Soft-start summaries include only their explicitly reported converged
seeds. Better starting coordinates and solver convergence are therefore useful
numerical results, not sufficient evidence for an improved anatomical method.
The existing reversed-domain-coupling limitation remains.

## Actual-data trajectory comparison (notebook 21)

`21_healthy_reconstruction_trajectory_benchmark.ipynb` compares the three current
atlas initializations with DPT and nonbranching scFates in both healthy-control
transfer directions and all three gene folds. Each method uses the same
training-only, marker-excluded count representation. DPT/scFates receive coarse
endpoint orientation; the atlas additionally uses weak ordinal labels in fitting.
These refits are distinct from historical pooled Harmony/marker coordinates.

The primary comparison uses an identical frozen 15-neighbor query extension;
native Gaussian atlas projection is a named secondary comparison. Response
splines and variance-normalized errors are identical. Coarse-label and privileged
coverage baselines remain. Failed fits are visible; gene folds are not biological
replication. Physical depth is excluded. The real-control result does not justify
replacing scFates with the current latent optimizer; the next hypothesis tests
conditional count references on a conventional backbone. Literature rationale and
limits are recorded in `docs/research/pseudostructure-methods-survey.md`.

## Fixed-backbone reference ablation (notebook 22)

`22_scfates_repeated_count_reference.ipynb` holds the training scFates coordinate
and response-gene decoder fixed. Its three frozen query projectors are neighbor
interpolation, a kernel mean reference in count PCs and a conditional NB
count-rate reference. It checks reproduction of notebook 21, paired quarter-read
perturbations and gene-panel coordinate sensitivity. Raw-library exposure is
computed over all measured genes; positional markers are excluded from modeling.
This is a single-source validation of a repeated-object reference, not the final
joint multispecimen estimator or a demonstration of reconstructed nephrons.

The NB projector improves read robustness but worsens all-gene transfer error in
every fold. No alternative is adopted from this experiment. Cross-panel ranges
include differences in numerical parametrization; inspect ranks as well as
coordinate differences, and never interpret conditional entropy as calibrated
anatomical precision. Its separate focused-program check keeps notebook 18's
training-selected panels frozen and compares all projectors against the
privileged quadratic segment-plus-coverage oracle. Those panels were selected
using the earlier atlas, so this is development-informed evidence rather than a
neutral method ranking or independent validation. Depth does not enter any score
or selection decision.

## Training-residual metric (notebook 23)

`23_repeated_structure_residual_metric.ipynb` tests one global covariance around
the same frozen scFates mean path using actual control transfers. The reusable
`residual_metric_atlas` helper interpolates that mean at the supplied training
positions, fits a zero-mean Ledoit–Wolf residual covariance and floors singular
eigenvalues numerically. It returns frozen Gaussian grid probabilities and
conditional/marginal log densities. This integrates observed variation; without
nephron membership it does not estimate a nephron-specific random trajectory or
separate technical and biological variance.

A matched spherical control uses the covariance trace to separate directional
geometry from residual scale. It was specified after the first covariance run,
with no query parameter tuning; notebook logic was bumped for this added arm.
The notebook reproduces saved positions and response scores, retains the
previous NB/coverage comparators and frozen gene panels, and repeats paired read
thinning. The first covariance result does not improve the neighbor reference;
neither more flexible covariance nor flow follows automatically from this test.

For these count-reference benchmarks, held-out genes and markers are excluded as
individual modeling features, but raw-library exposure uses all measured genes,
including those panels. Response normalization shares that exposure. The tests
therefore retain aggregate compositional dependence and are not fully independent
of held-out measurements. Query labels also inherit upstream molecular annotation;
source anatomical orientation and frozen-panel selection limit claims of
independent anatomical validation. These limits apply to both notebooks 22 and 23.

## Restricted specimen calibration (notebook 24)

`24_repeated_structure_specimen_registration.ipynb` tests the next hierarchical
special case, a constant specimen offset around a frozen source-only canonical
PT curve. It learns the count representation and scFates backbone from one
control, fits calibration offsets using a deterministic hash half of the other,
and evaluates the disjoint half. The restricted arm removes components
orthogonal to the source mean-path subspace retaining 99% variation; the
unrestricted arm uses the same penalty 4 and five update steps. These are
geometric restrictions, not identified biological injury directions. A full-rank
source path permits an exactly zero restricted offset.

The S3-only regime constructs a deliberately incomplete calibration population
using annotation, then hides those labels from the estimator. Evaluation rows
stay identical between regimes. The notebook retains zero-calibration DPT and
scFates references and the privileged quadratic segment-plus-coverage oracle;
the calibrated arms explicitly receive extra unlabeled query data. Frozen
positional-program panels come from earlier training-only selection. Paired
quarter-read inference freezes calibration parameters and uses original response
measurements as targets. Unknown nephron membership and shared count exposure
continue to limit anatomical and independence claims. No joint final estimator,
disease transfer, human alignment or regulator analysis is adopted from this
single-source calibration experiment.

## Goal check before each biological experiment

State which part of repeated-structure reconstruction is being tested: the shared
ordered object, inference from unordered partial sections, restricted specimen
variation, or the reliability of frozen projection. Name the matched DPT/scFates
and segment/coverage comparators and specify a result that would reject the
proposed addition. A normalization or transport benchmark alone does not show a
new reconstruction method. Record the hypothesis before reading outputs.

The next core protocol in the research ledger concerns groupwise atlas fitting
with source-only anatomical anchors and disjoint unlabeled calibration/evaluation
rows. Its pooled trajectory comparators must receive the same calibration data.
Unknown calibration anatomy contributes no ordinal likelihood or root/orientation
information. Existing single-source null results remain evidence against adopting
unnecessary components. Depth remains weak descriptive context because a winding
nephron and unknown cross-section membership break its interpretation as
longitudinal ground truth.

The reusable joint-atlas fitter now permits `-1` anatomy only with explicit
`allow_unlabeled_anatomy=True`; those rows receive no ordinal term in its
posterior, cutpoint fitting or recorded objective. Trajectory baselines accept an
explicit boolean `anchor_mask`, validate S1/S2/S3 only among anchors, and use them
for rooting and orientation while retaining every point as a graph/curve
observation. Both APIs preserve their previous strict defaults. These changes
prepare the groupwise comparison; they are not results from that experiment.

## Groupwise reconstruction (notebook 25)

`25_groupwise_repeated_structure_reconstruction.ipynb` directly tests whether
separate specimen point clouds inform a shared latent PT object beyond a pooled
trajectory. The source-only count transform is frozen. Each fit uses all source
sections plus a deterministic calibration half of the other control, with the
remaining half reserved for evaluation. Full and S3-only calibration regimes keep
identical evaluation IDs. S3 annotation selects the stress population but is then
hidden from the fitter.

Pooled DPT/scFates receive the same calibration rows. Source labels root and
orient the graph baselines and provide a weak ordinal likelihood to the atlas;
these are different uses of the same source information. The three latent-atlas
variants separate equal-specimen weighting from constant specimen offsets. The
common-neighbor comparison isolates fitted axes; native expression-only atlas
projection is separately named. Every response decoder uses source responses
only. The privileged segment-only and quadratic segment-plus-coverage baselines
receive evaluation labels explicitly; molecular projection does not.

Unknown nephron identity is integrated variation, not estimated membership.
Uniform grid spacing fixes a numerical coordinate and does not infer physical
length or sampling uniformity. Constant offsets do not identify disease effects.
Frozen source-selected programs, gene-fold disagreement and paired read thinning
are descriptive development diagnostics on two repeatedly reused mice. Raw
all-gene exposure retains aggregate dependence on withheld responses. No depth
measurement contributes to selection, and no disease, human or regulator branch
is enabled by this experiment. Solver failures remain visible and are excluded
from accepted prediction summaries rather than silently treated as convergence.

The Gaussian likelihood still uses uniform grid mixing and does not estimate
specimen-specific sampling frequencies along the latent object. This is distinct
from hard occupancy matching and from physically uniform cuts, but may affect
partial-coverage fitting. Interpret the stress test as evidence about the tested
likelihood, not an identification guarantee or a rejection of the shared-object
prior.

The first notebook-25 run accepts 21/36 joint fits and does not establish a
consistent benefit over the source-only reference. Accepted-fit summaries show
denominators, and a post-output reporting appendix compares arms on identical
successful folds without changing any fit. Source-only comparators are visible
in the program plot because gains over an unstable pooled scFates baseline are
insufficient. The groupwise prototype remains experimental and unadopted;
full measured conclusions and the next unit-convention check are in the ledger.

## Reference units for the hierarchy (notebook 26)

`26_repeated_structure_reference_units.ipynb` keeps notebook 25's actual-control
design and changes only the joint fitter's prior reference scalar. For each
source/fold, it uses the RMS of centered source PCs across rows and components.
This is total source molecular variation, not a residual-noise estimate. PCA,
graph inputs, initialization, ordinal supervision and evaluation remain fixed.
The source-only and pooled graph coordinates, quarter-read coordinates and
scores must reproduce the preceding benchmarks before new joint fits proceed.

The reusable fitter's `prior_scale=1` default preserves previous behavior. With
an explicit scalar it fits in reference units and returns molecular parameters
and variances in input units; objectives and stopping criteria use reference
units. Scaling both inputs and the scalar preserves the fitted positions.
Changing the scalar for fixed data changes prior strength. Objective values
across conventions are therefore not model-selection evidence.

Raw-unit versus source-relative predictions are compared on identical accepted
cases, with unmatched cases and fold counts reported. The same source-only,
pooled, segment and coverage comparators remain required. Numerical invariance
does not identify the anatomical coordinate or resolve partial coverage.

The actual-control run accepts 5/36 source-relative fits, including none of the
24 specimen-balanced attempts, versus 21/36 raw-unit fits. All graph references
reproduce. Sparse accepted-fit gains do not consistently beat source-only
scFates, and the hierarchy remains unadopted. Numerical unit consistency is
retained as an explicit API property; see the ledger for matched denominators,
program and read results, and the decision to stop expanding this branch.

## Assignment-aware response prediction (notebook 27)

`27_repeated_structure_assignment_prediction.ipynb` freezes notebook 22's
scFates/count references, adds the existing source-only DPT comparator, and
distinguishes decoding at the mean position from averaging decoded responses over
assignments. `SplineAtlas.predict_distribution` averages spline basis values
before multiplying gene coefficients. `project_neighbor_axis` optionally exposes
its original supports/weights; its default coordinate output is unchanged.
Neighbor assignments are not calibrated probabilities.

Every arm uses both response reductions and identical frozen training decoders.
Actual-control guards reproduce earlier point positions/errors, selected programs
and paired thinning before interpretation. The six transfers complete; averaging
strengthens scFates' prediction while Gaussian/NB projection remains inferior in
important programs, especially reverse S2. Coordinates and ordering are unchanged.
Retain this explicit reduction distinction in subsequent comparisons; do not use
its prediction gains to claim a reconstructed physical coordinate. Both the
probabilistic projectors and unconstrained Gaussian hierarchy remain unadopted.
The next curve-construction hypothesis and all measured denominators are in the
research ledger. Depth remains excluded.
