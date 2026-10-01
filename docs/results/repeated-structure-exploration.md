# Notebook 15 initial exploration

The required current `mouse_only_v6` export was unavailable during implementation.
No v5 results were substituted. The real-data diagnostic, outer mouse comparison,
and biological consequences are pending. No method superiority or novelty is
claimed. Two healthy mice also cannot provide replicate-aware training with an
independent healthy test specimen.

The initial synthetic benchmark used six specimens, 80 sampled profiles per
specimen before missing-region filtering, 120 genes, three seeds, nuisance
strengths 0/2/6, and six scenarios. Each fit held out the final specimen, excluded
half the genes from coordinate construction, and selected a 20-gene panel inside
training. The seven comparisons produced 378 method/scenario/seed records;
disconnected DPT fits were explicitly recorded rather than repaired by invented
edges. Synthetic anatomy is supplied to the prototype and provides substantial
gross-order information; within-segment metrics are therefore essential.

The fixed exploratory criterion did not find consistent fine-order improvement
over the conventional baselines. This is a failure to establish superiority,
not evidence that the principle can never work. Initial nuisance stress tests
sometimes favored the atlas over DPT, while PCA remained competitive. Three seeds
and one held-out specimen per simulation are insufficient for a robust general
performance claim. The default settings were not optimized to rescue the story.

The current estimator is a heuristic combination of standard spline regression,
replicate-based feature scoring, constrained local projection and frozen-atlas
mapping. No new objective is established as mathematically identifiable. Any
monotone coordinate warp leaves the observation model equivalent, and orientation
requires an anchor. The literature table documents close spatial reconstruction,
latent trajectory and cross-sample alignment precedents; its focused search does
not prove absence of an exact earlier formulation.

Geneformer V2 representation inference remains deferred because aggregate mouse
tubule profiles do not satisfy its documented human single-cell training/input
setting. A separately sourced, documented external gene prior can be compared
without conflating feature-prior effects with pretrained-representation effects.
The implementation was exercised using mock prior and TF files only as code
contract tests; those runs provide no Geneformer or biological evidence.

Validation included the full synthetic notebook run, a separate synthetic AnnData
stand-in exercising the real-data code paths and optional prior/TF routes, all
analysis tests, and notebook stage-order checks. Existing working notebooks and
outputs were preserved. See the [workflow guide](../workflows/repeated-structure-pseudospace.md)
for execution and artifact contracts. Synthetic outputs and private results stay
outside Git.
