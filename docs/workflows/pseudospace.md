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
