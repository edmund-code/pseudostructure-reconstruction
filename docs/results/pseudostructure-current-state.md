# Pseudostructure research: current checkpoint

Updated 2026-10-05 after notebooks 52–53 (paper workstream A). This short entrypoint avoids reopening
closed branches; the [research ledger](pseudostructure-research-ledger.md) holds
protocols fixed before outcomes, complete denominators and measured results.

## Goal and assumptions

Recover one canonical ordered PT distribution path from unordered cross-sections
of repeated physical realizations. Individual nephron memberships and longitudinal
lengths are unknown. Healthy, AKI and human retain the same structural coordinate
while molecular realizations may differ. Human cortex/medulla-labelled specimens
are both healthy cortex from one donor. The final method is **not achieved**.

Keep physical-object inference distinct from generic molecular progression.
Before every experiment: identify the reconstruction failure it changes, explain
how the repeated-object prior enters inference, freeze its protocol, and specify
what measured result would justify retaining it. Existing mathematics is fine;
adding complexity or reframing an established curve is not enough by itself.

## Strongest current components and evidence limits

Actual healthy controls, marker-excluded source-only count PCs, independent gene
folds, frozen source-scFates/DPT axes, and source neighbor assignment distributions
form the strongest current molecular-reference benchmark. Average response
predictions over assignments explicitly; prediction at their mean is different.
Assignments are not calibrated anatomical posteriors. DPT generally has stronger
within-segment gene-panel rank stability; that is not anatomical ground truth.

Frozen development programs were source-selected in notebook 18. They are useful
comparators, not independent fine-position truth. Shared all-gene count exposure
couples modeling and response genes. Source segment labels inherit marker
information. Two repeatedly reused controls support descriptive development.
Depth is excluded because winding nephrons and unknown memberships weaken it as
a longitudinal target. Injury, species and regulator interpretation remain gated.

## Closed hypotheses

| Experiment | Measured decision |
|---|---|
| 15 OT probability flow | Centroid shifts predict better; no state-dependent field adopted. |
| 16–20 Gaussian/count atlas and starts | Early pilots encouraging; stronger real-data references later defeat broad claims. |
| 21–22 actual DPT/scFates/count comparison | Source scFates neighbor reference strongest for selected programs; NB projection loses important fine programs. |
| 23 residual covariance | No consistent prediction/read benefit. |
| 24 restricted/unrestricted constant offset | Small or heterogeneous effects; neither adopted. |
| 25–26 joint Gaussian hierarchy/prior units | Sparse convergence and no consistent reference advantage; retain explicit unit API, not the hierarchy. |
| 27 assignment-aware response reduction | Improves nonlinear response prediction, without changing positions; retain fairer benchmark. |
| 28 ordered anatomical initialization | All 18 new fits valid but no consistent benefit; coarse median order insufficient. |
| 29 independent query path/partial ordered registration | 11/12 eligible; full prediction worsens, partial ranks collapse and reverse S2 prediction fails severely. Do not adopt. |
| 51 landmark count-position model (LCP; no x/y) | Within-segment stability fixed (gene-fold 0.96 vs scFates 0.60; LOSO 0.99) and P1w up (0.34 vs 0.17), but abandonment rule (iii) triggered: held-out-gene gain over the segment + coverage oracle in only 3/6 mouse cells. No coordinate (scFates, DPT, LRS, LCP) beats that oracle (P4c ≤ 0 in 44/48 cells). A zero-cost landmark ratio score gets most of the gain. Not adopted. |
| 52 landmark model ablation ladder | Count likelihood, zone shapes, specimen balancing and joint species fitting each failed their pre-registered margins over the free landmark ratio; calibration abandoned (bootstrap 2–4× conservative). Not adopted as a model. |
| 53 zone-derived landmark ratio (LRS-Z; landmarks from the data's own zones, read-split, no atlas, no x/y) | Passes the amended adoption criteria (P4c as claim gate, decided before fits): within-segment gene-fold agreement 0.93/0.95 (scFates 0.60/0.21), human S1 0.96 (0.07), external within-segment concordance 0.38/0.19 (0.17/0.10), positional noise 0.21–0.31 of within-segment spread (0.83–1.49). P4c still ≤ 0 for every arm. Caveats: between-zone order equals the labels; registration gap tautological in zone units; uncertainty conservative. Adoption as the paper coordinate awaits the user. |

Notebook 29 uses full/S3-only calibration with identical evaluation halves,
free reference endpoints, and common expression-based support/frozen-reference
fallback. Exact source/oracle reproduction and fallback guards run after cache
hits too. Its support heuristic admits some S2 profiles during S3-only calibration,
more after thinning; it does not identify physical overlap. One low-support fit
is rejected without tuning. Local cache key: `bd25e5defdfbfd54`; 333 tests pass.
Working notebooks/results stay local; committed notebook copies are output-free.

## Next scientific checkpoint

Do not repeat initialization, mean-path registration, covariance, offset or flow
variants without a new reconstruction-specific hypothesis. The evidence favors
retaining local empirical reference information over forcing query means onto a
new independently ordered curve. The conditional measurement model remains a
candidate failure source: count exposure changes the Pearson-residual geometry,
while thinning cannot move the physical section. Its role has not been isolated.
Review that observation contract before another private fit; require the same
strong source/segment/coverage comparators and explicit partial-support behavior.

A future method still needs a coherent estimator, validated healthy transfer and
fine-order stability, repeated-specimen contribution beyond a generic pooled
curve, honest conditional uncertainty, and frozen position-versus-domain-deviation
behavior. None can be claimed solely from a smooth curve, a green regression
suite, preserved segment medians or a lower fitting cost.

## Stopped checkpoint

Work stopped at the user’s request. Notebook 30 contains only the prespecified
measurement-invariance protocol; it is not implemented or executed. Notebook 29
remains the latest measured result. No research processes are running.

## Paper re-centring (2026-10-05)

The paper now introduces a reconstruction method and asks what continuous analysis shows that
discrete segments do not. The user's constraint: the method must stay general for continuous repeated
structures (for example intestinal crypts), so tissue x/y is never used in inference. Notebook 51
(LCP) is the first experiment under this goal; notebooks 56–60 (pathway workstream) measure
continuous-only information per coordinate. See `docs/results/pt-reconstruction-v3-phase1.md` and
`docs/results/pt-pathway-continuous-vs-discrete.md`.
