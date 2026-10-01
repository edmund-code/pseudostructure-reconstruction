# Notebook 15: measured-data exploration

Notebook 15 now uses only the existing project measurements. The simulation
benchmark, generator, simulation-specific tests and simulation-derived verdicts
have been removed. Method conclusions must come from held-out specimen/gene
prediction and robustness on the measured tubule profiles.

The old requirement for a `mouse_only_v6` DPT export has been removed. Notebook 15
loads measured counts and segmentation annotations directly. The primary estimator
initializes independently of DPT. Newly fitted PCA/DPT methods are comparators;
a saved DPT coordinate is only an optional, input-matched diagnostic.

Actual input verification succeeded: 35,581 mouse structures matched their current
segmentation feature indices and centroids (maximum deviation 0 pixels). The pooled
5th-percentile spot floor was 74, leaving 33,834 structures before PT selection and
14,861 PT profiles with 19,059 common measured genes. All methods within 15 share
this newly normalized universe. No method superiority or novelty follows from
this input verification; the full reconstruction comparisons remain to be run.

Two healthy mice cannot support replicate-informed training plus an independent
healthy test specimen. Four-mouse LOSO remains an AKI invariance stress test.
Gene/marker holdout, anatomy ablation, reflected initialization, shuffled controls
and bootstrap resampling operate on the measured data. These checks do not supply
exact fine spatial ground truth or establish physical coordinate calibration.

The estimator combines standard spline regression, replicate-based feature
scoring, local projection and frozen-atlas mapping. The focused literature review
does not prove absence of an earlier equivalent objective. Geneformer representation
inference remains deferred until defensible cell-level inputs are available; its
optional external gene-prior branch stays separate.

New artifacts are written beneath `novel_repeated_structure_pseudospace/real_data/`
to keep them separate from superseded experiment artifacts. Previous working
notebook contents were preserved in an ignored local backup before replacement
of obsolete experiment/verdict cells. Unaffected live outputs were retained.

See the [workflow guide](../workflows/repeated-structure-pseudospace.md) for inputs,
execution and interpretation limits. Unit tests use small numerical fixtures to
verify implementation; they are not biological experiments or method evidence.
