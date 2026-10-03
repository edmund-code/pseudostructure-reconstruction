# Notebook 15: latent spatial probability-flow proof of concept

The canonical notebook is
`analysis/notebooks/15_spatial_probability_flow_proof_of_concept.ipynb`.
It replaces the previous repeated-structure spline exploration. There is no
notebook 15 mirror. Reusable numerical logic lives in
`pseudospace/spatial_probability_flow.py`; the measured-input loader is
`pseudospace/repeated_inputs.py`.

The decisive experiment is whether a smooth molecular velocity fitted on one
healthy PT specimen predicts adjacent molecular distributions in the other.
Ctrl1A2 and Ctrl1A4 exchange training/test roles. Two controls provide descriptive
evidence, not confirmatory inference or calibrated longitudinal position.

## Inputs and execution

Use `kidney-pseudospace` and run cells top to bottom. Root flags `--data-root` and
`--results-root` override `PSEUDOSPACE_DATA_ROOT` and `PSEUDOSPACE_RESULTS_ROOT`.
Required measured inputs are the two control matrices under
`<data-root>/tubule_by_gene/{sample}_tubule_by_gene_caleb.h5ad` and matching
`<data-root>/{sample}_kept_tubules_labeled_fine.geojson` files. AKI is optional.
The loader verifies segmentation joins and centroids, retains all measured common
genes, applies a pooled healthy-cohort p5 spot-count floor, normalizes to 10,000
counts and log1p, and selects upstream PT. There is no tubule gene-count QC.
Its structural-floor cohort differs from the default four-mouse loader and is
recorded in the input manifest.

The default temporary coordinate is
`<results-root>/human_vs_healthy_mouse/cross_species_pt_dpt.h5ad`, column
`total_scanpy_dpt`. Override it with `PSEUDOSPACE_15_DPT_REFERENCE` and
`PSEUDOSPACE_15_DPT_COLUMN` (for example a mouse-only export's
`pt_subset_scanpy_dpt`). This requires the existing upstream reviewed coordinate;
15 does not recompute DPT or bypass label checkpoints. A strict unique centroid
join within one pixel links the saved axis to current validated PT. Feature IDs
must also agree if available. Structures absent from that export are reported in
`scaffold_ineligible_structures.csv`; the fixed-axis experiment uses only the
explicit scaffold-eligible intersection. Healthy DPT is scaled jointly, never
independently by specimen. The default 03 scaffold has cross-species origins;
only control molecular profiles enter this experiment.

Physical evaluation requires the reviewed pass-1 glomerular export, default
`<results-root>/human_vs_healthy_mouse/cross_species_harmony_pass1.h5ad`.
Override with `PSEUDOSPACE_15_GLOMERULI`. Glomerular feature IDs and centroids are
verified against current segmentation. Nearest and mean-three-glomeruli distances
are computed after molecular fitting, never supplied to reconstruction. They are
2D depth proxies, not longitudinal ground truth. Both human-labeled slices are
healthy cortex; neither enters the molecular fits.

## Experiments and gates

1. Exclude canonical positional markers; rotate deterministic reconstruction and
   evaluation gene folds. Fit scaler/PCA only on the training control (15 PCs).
2. Inspect 12 overlapping DPT windows and a non-overlapping sensitivity. Unsupported
   windows stop the fixed-axis experiment. Bounded clouds use at most 128 structures.
3. Fit entropic-OT barycentric displacements between training windows, then distill
   them into a ridge-regularized affine velocity with smooth position dependence.
   Entropy is relative to each pair's median positive squared distance. Integrate
   with RK4 and score held-out distributions with multivariate energy distance.
4. Compare identity, training centroid displacement, and training DPT-local kNN
   mean displacement. Flow must strictly beat all three baselines, averaged across
   transitions, in both directions and both window schemes. The gate is descriptive.
5. Only after that gate, alternate density-atlas estimation and posterior position
   reassignment. Compare empirical Gaussian atlas moments against moments blended
   with one-step flow predictions. Ordered anatomical transition cutpoints are
   learned; no fixed segment thirds or expected marker curves. Compare ordinal,
   endpoint-only, and unsupervised label settings, and DPT/PCA1/random ordinal starts.
   Test-specimen projection uses neither labels nor its DPT.
6. Score held-out gene curves across three gene folds, initialization agreement,
   tubule/gene/both bootstraps, half-sampling, posterior uncertainty, excluded
   markers, and specimen-specific physical depth. Failures remain explicit.

The latent screen requires the prespecified DPT-initialized ordinal flow atlas
to beat atlas-only and supplied DPT on held-out genes in every fold and direction.
Endpoint-only PCA/random initializations omit S2 information; no-anatomy
PCA/random starts omit all labels. DPT initialization retains upstream information
in all ablations. Gene holdouts retain compositional dependence through the common
total-count normalization denominator.
This does not replace scientific review of support, stability or physical depth.
No specimen offsets are learned from held-out controls. A distributional source
window is observed in the test specimen in Experiment 1; its successful transport
is not itself de novo coordinate reconstruction.

## Optional extensions and artifacts

Outputs go to `<results-root>/spatial_probability_flow_proof_of_concept/`.
Manifests, per-transition scores, window support, gates, coordinate/posterior
probabilities, failed fits, bootstrap positions, marker curves, physical scores
and diagnostic PDFs remain private, untracked artifacts. Expensive stages cache
by actual inputs, parameters and code under `stage_cache/`.
`PSEUDOSPACE_STAGE_CACHE=0` rebuilds them. Validation is never cached.
`PSEUDOSPACE_15_BOOTSTRAPS` defaults to 10 and must be at least two; reducing it
is a smoke-run choice, not stronger evidence.

`PSEUDOSPACE_15_AKI=1` enables a gated frozen-healthy-reference projection of
IR2A2/IR2A4, using the healthy measured-gene normalization denominator.
The pooled healthy atlas uses exactly balanced specimen clouds at every window
and refuses sparse per-specimen support. Injury
never changes the axis. Position, posterior entropy and negative healthy mixture
log-density are separate exploratory outputs; the density score is not calibrated
injury severity. Injury's own cohort structure floor is recorded.

`PSEUDOSPACE_15_TF=1`, `PSEUDOSPACE_15_TF_LIST` (curated mouse CSV with unique
`gene`) and `PSEUDOSPACE_15_VALIDATION_REVIEWED=1` enable a gated TF feasibility
panel after reviewing validation. A chain-rule Jacobian measures low-dimensional
velocity sensitivity to modeling-gene log-expression, not causality or TF activity.
For this affine field it is position-dependent but state-independent. Geneformer
may later prioritize an external TF panel; it never creates the coordinate.

## Interpretation ceiling

This is finite-step OT distillation and a flow-smoothed Gaussian-density pilot,
not a reproduction of PFI, exact EM, or enforcement of the continuous continuity
PDE. Affine fields cannot represent arbitrary multimodal trajectories and
barycentric transport can contract variance. No drift/diffusion decomposition is
attempted. Metric coordinate, sampling density, and biological dynamics remain
unidentified. DPT and upstream annotations retain indirect marker information;
marker exclusion does not make those diagnostics independent. Overlapping windows
share observations, which motivates the non-overlap sensitivity. Synthetic tests
verify implementation only. No superiority or novelty follows from writing or
successfully executing the notebook.

Before commits run synthetic tests, the notebook stage-order gate and repository
hygiene. Stage output-free notebook copies without clearing a live notebook.
