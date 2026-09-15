# Active pseudospace workflow

## Inputs and roots

`analysis/mouse_only_pseudospace.py` accepts `--data-root` and `--results-root`; the equivalent
environment variables are `PSEUDOSPACE_DATA_ROOT` and `PSEUDOSPACE_RESULTS_ROOT`. Command-line
arguments take precedence. The repository's `data/` and `results/` directories are only safe
defaults for authorized local mounts.

The data root must follow [data/README.md](../../data/README.md). The result root receives a
`mouse_only_v5/` directory and is never committed.

## Run sequence

1. Regenerate tubule matrices with `analysis/notebooks/01_segmentation_to_gene_matrix.ipynb` only
   when the private Visium and v4 segmentation inputs change.
2. Create the pinned `kidney-pseudospace` environment from `environment.yml`, then run the
   mouse-only script or notebook from a clean kernel. The workflow requires R `harmony` 2.0.5;
   it deliberately refuses older user-level R installations.
3. At the coarse-label checkpoint, inspect the dotplot and explicitly confirm every Leiden
   cluster. The reference labels are pinned to the 13-cluster reference fingerprint; if it
   changes, the workflow stops before any cells are retained for DPT. Re-read the dotplot and
   update both the labels and fingerprint together—do not bypass the guard.
4. Review Harmony integration diagnostics, PAGA connectivity, DPT-by-segment ordering, and the
   physical-axis sensitivity analysis before interpreting condition effects.
5. Use the QuPath scripts with the same roots to export labels or perform spatial validation.

## Scientific guardrails

- Only `*_v4.geojson` mouse segmentations are valid; centroid verification is mandatory.
- The canonical coordinate is Scanpy DPT on the pass-2 Harmony embedding, rooted in PT and
  oriented by the early-to-late marker axis.
- The cohort is 2 versus 2 biological specimens. Shape and pathway ranks are effect-size
  descriptions, not significant confirmatory findings.
- The PT arm is currently the defensible trajectory. Treat the non-PT limitations detailed in
  [the results interpretation](../results/current-mouse-run.md) as active constraints.

## How results are reported (level, amplitude, pattern)

`levelshape.summarize_curve_effects` splits every fitted curve pair into a **level** offset, an
**amplitude** change and a **shape/pattern** change, with `level_fraction + shape_fraction = 1` from
the exact orthogonal split of the difference. Report all three: on the current mouse-human fits the
top of the ranking is 98-99 % vertical offset, and a pure amplitude loss keeps a non-zero
`shape_rms`. Every figure that compares two conditions is written twice — a common-scale view
(fitted lognorm units, one colour scale) and a shape-only view (each curve standardised) — because
one pooled z-score answers neither question.

## Pathways

Membership is resolved through the accepted ortholog map (not by case-collision), each pathway
carries its coverage stages (`n_requested`, `n_with_ortholog`, `n_assayed`, `n_tested`), excluded
pathways are kept with a reason, and no upper size bound is applied to membership — the number of
pathways such a bound would remove is reported instead. Pathways whose member sets largely coincide
are grouped in `pathway_redundancy_*.csv`, and each prioritized pathway ships member-gene evidence
(agreement fraction, strongest contributor, sign flip without it). Enrichment is competitive, uses
the genes eligible for the analysis as background, and reports both a label-permutation p and a
correlation-aware p; correction is applied across all tested pairs.

## Curve modules

`pseudospace.modules` discovers **positional** modules from the reference curves and **response**
modules from the mean-centred difference curves by correlation distance on standardised curves.
Peaks are annotations, never the grouping key, and dynamic time warping is not used. Stability is
reported under specimen omission and under a different data mixture.

## Specimen weighting and magnitude

Specimen-balanced curves (equal weight per specimen) are the primary descriptive summary; the pooled
fit and individual specimen curves stay visible. Pseudobulk profiles are written per specimen and
pseudospace bin for count-based modelling — bins from one specimen are repeats along one coordinate,
so the independent units remain the specimens. Before any condition or species claim, check
detection/abundance per side, ratio-versus-abundance, within-side gradients, the matched-position
per-specimen contrast, and the capture summary.

## Caveats that must travel with the numbers

- The 2-vs-2 permutation p-floor is **2/6**, not 1/6: six relabelings form three mirror-image pairs
  whose statistic is identical, and ties count as "at least as extreme".
- Peak positions are coordinates on one notebook's own DPT construction (notebook 02: global-nephron
  PT coordinate; notebook 03: PT-specific recomputed DPT) and are **not comparable across the two**.
- Normalisation divides by the retained panel's total after the expression filter; the native
  denominator is recomputed as a sensitivity check (`normalisation_sensitivity.csv`).
