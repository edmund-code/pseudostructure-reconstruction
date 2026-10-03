# Notebook 15: healthy canonical atlas first

Notebook 15 now separates healthy molecular atlas mapping from the probability-flow
hypothesis. Its canonical source remains
`analysis/notebooks/15_spatial_probability_flow_proof_of_concept.ipynb`.
The October 2 OT-distilled affine flow lost to centroid displacement in both
specimen directions and both window schemes. That result no longer gates the
atlas experiment. The superseded spline-atlas verification numbers are not results
for this protocol.

## Measured atlas-first run, 2026-10-03

A measured source-cell run completed the revised comparisons using Ctrl1A2 and
Ctrl1A4, exchanging training and test roles across three marker-excluded gene
folds. All 5,598 validated healthy PT profiles were retained as held-out projection
targets, including 67 without saved DPT. Each training fit used only its own
DPT-available profiles, training-only scaler/PCA and training-only coordinate
normalization. Held-out projection received molecular features alone.

The canonical mean atlas had lower variance-scaled held-out-gene MSE than the
strong S1/S2/S3 mean oracle in every gene fold and direction. The relative reduction
averaged about **0.19%**, ranging approximately 0.17–0.20% across the six primary
comparisons. Within-segment reductions were also small. This is a consistent
but modest descriptive increment, not evidence of substantial method superiority.
The oracle baseline uses test segment labels; the atlas does not.

Expression-only mean-atlas positions preserved the physical-depth relationship.
Across gene folds, mean-three-glomeruli Spearman correlations on the DPT-available
comparison cohort were approximately **0.732 for Ctrl1A2** and **0.692 for Ctrl1A4**,
versus supplied DPT's 0.725 and 0.688. These are 2D cortical-depth proxies, not
longitudinal PT ground truth; the small differences are not confirmatory evidence.
Withheld segment order and excluded marker curves remained coherent. Inherited
training DPT and upstream annotation information still limit independence claims.

Gaussian likelihood mapping had higher held-out-gene error than the mean atlas in
all six primary comparisons. For the distributional transport ablation:

| Windows | Training control | Smooth mean energy | Gaussian-flow energy |
|---|---|---:|---:|
| Non-overlapping | Ctrl1A2 | 0.678 | 0.651 |
| Non-overlapping | Ctrl1A4 | 0.611 | 0.683 |
| Overlapping | Ctrl1A2 | 0.402 | 0.352 |
| Overlapping | Ctrl1A4 | 0.569 | 0.639 |

Lower energy distance is better. Covariance flow improved one training direction
but worsened the other in both schemes, so it did not pass the incremental screen
against its identical smooth mean path. Centered shape metrics are saved separately;
an exact Gaussian marginal-preserving flow shares its atlas's likelihood and is
not a separate coordinate estimator.

A fixed-penalty global-intercept sensitivity, calibrated on a deterministic half
of test modeling-feature profiles and applied to the disjoint validation half,
changed predictive scores very little. It did not learn evaluation-gene offsets
or position-dependent effects. Sampling composition and intercept/position
confounding remain possible; unadjusted projection is primary.

Ten bootstrap fits per scheme and direction (training tubules, genes, and both)
completed without sparse-support failures. Median per-profile position SD was
approximately **0.008–0.018**, but about **2.3–5.5%** of profiles had SD above 0.1,
and maximum SD reached approximately 0.44. The stable majority does not eliminate
large positional ambiguity in the tail. Per-profile SD/range and full conditional
posteriors are retained; neither is calibrated spatial confidence.

## Decision and artifacts

The simple mean atlas warrants further healthy-specimen validation. The small
continuous-versus-segment predictive increment and ambiguous minority prevent a
strong biological or methodological claim. Gaussian density/flow has not yet
earned a consistent added predictive benefit. AKI projection and TF Jacobian
interpretation remain deferred.

Two biological controls support descriptive comparisons only. Private matrices,
segmentations, per-profile results, curve tables and diagnostic plots remain under
the configured results root. Live notebook outputs are preserved locally; only an
output-free copy is committed. Synthetic tests verify numerical behavior and
input separation, not biological validity. See the
[workflow guide](../workflows/repeated-structure-pseudospace.md) for contracts,
comparisons and interpretation limits.
