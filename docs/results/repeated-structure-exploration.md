# Notebook 15: spatial probability-flow experiment status

Notebook 15 has been replaced from scratch by
`analysis/notebooks/15_spatial_probability_flow_proof_of_concept.ipynb`.
The previous spline-atlas design and its input-verification numbers are not
results for this new experiment.

The current protocol first tests cross-specimen distribution prediction on the
existing PT DPT scaffold using healthy controls Ctrl1A2 and Ctrl1A4. Latent-position
inference is gated on beating identity, centroid and DPT-local baselines in both
specimen directions, including non-overlapping windows. Atlas-only versus
flow-smoothed density, rotated held-out genes, uncertainty and stability are
separate tests. Physical glomerular-depth proxies enter only after inference.

A measured source-cell run on 2026-10-02 completed the fixed-axis experiment and
physical-depth diagnostics. The affine flow did not pass the prespecified gate:
training-centroid displacement had lower held-out energy distance in both
specimen directions and both window schemes. Latent fitting, bootstrap fitting,
AKI projection and TF analysis were therefore skipped. This rejects advancement
of this particular approximation; it does not rule out every probability-flow
model. Detailed scores and manifests remain private results artifacts.

A successful numerical or input check does not establish biological validity,
method superiority, causality or novelty. With two controls, all comparisons are
descriptive. AKI frozen-reference projection and TF sensitivity are optional
extensions requiring additional validation; neither defines the healthy axis.

The committed notebook carries no executed outputs. Measured results and run
manifests belong under the configured results root, outside Git. See the
[workflow guide](../workflows/repeated-structure-pseudospace.md) for the complete
input contract, gates and interpretation limits.
