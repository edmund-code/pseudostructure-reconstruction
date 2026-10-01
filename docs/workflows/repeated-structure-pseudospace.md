# Notebook 15: repeated-structure pseudospace exploration

`analysis/notebooks/15_novel_repeated_structure_pseudospace.ipynb` asks whether
cross-specimen molecular prediction improves a canonical order of independent
segmented profiles. It does not interpret nephron position as time or assert a
novel method. Existing notebooks are unchanged.

The supplied probability-flow reference is Maddu, Chardès and Shelley,
[PNAS 2025](https://doi.org/10.1073/pnas.2420621122). The notebook also discusses
[Pheno-GS](https://arxiv.org/abs/2609.27633) as a preprint and includes a focused
precedence table. The review is not an exhaustive absence-of-precedent proof.

## Inputs and execution

Use the `kidney-pseudospace` environment and run the notebook from the repository
or one of its subdirectories. CLI root flags override `PSEUDOSPACE_DATA_ROOT`
and `PSEUDOSPACE_RESULTS_ROOT`, as in the existing workflows. Notebook 15 reads:

```
<results-root>/mouse_only_v6/all_mouse_tubules_scanpy_dpt.h5ad
```

This must be the reviewed notebook 02 export from current QC segmentations,
with its centroid verification and coarse-label checkpoint completed upstream.
Notebook 15 consumes `layers['lognorm']`, selects the existing PT cohort, and
retains all profiles. It does not rerun integration, reapply count QC, infer
new labels, or fall back to v5. Unknown PT intervals or incomplete specimen
contracts fail explicitly. If the export is absent, real analyses are marked
pending and synthetic experiments still run.

Outputs go to `<results-root>/novel_repeated_structure_pseudospace/`; none are
committed. The canonical notebook is output-free. There is no generated mirror
for notebook 15. Do not clear live notebook outputs when staging.

The full run evaluates every available reconstruction gene plus panels of
10/20/50/100 genes. It can require several GB of RAM and substantial CPU time.
It caches fold fits by expression, labels, settings, source metadata and code.
Bootstrap fits are deliberately separate. Smoke-run controls are:

- `PSEUDOSPACE_15_BOOTSTRAPS` (default 20 per resampling scheme per outer fold).
- `PSEUDOSPACE_15_SIM_SEEDS` (default 3; increase for robust synthetic inference).
- `PSEUDOSPACE_15_PANEL_DPT=0` to omit optional Scanpy DPT baselines.
- `PSEUDOSPACE_STAGE_CACHE=0` to rebuild cached fits.

A smoke run with fewer repetitions or disabled DPT is not the complete benchmark.

## Estimator and validation boundaries

Reusable helpers live in `pseudospace/repeated_structure.py` and
`pseudospace/repeated_baselines.py`. The prototype combines conventional
specimen-balanced P-spline regression, local grid projection, training-only
feature selection, and anatomy interval constraints. The optional graph is a
standard harmonic coordinate, not a novelty claim. Disconnected DPT graphs
are recorded as failed baselines; no edges are invented to make them run.

The saved Harmony/DPT coordinate is a transductive reference. The initial
healthy DPT perturbation diagnostic uses both controls in selection/refinement;
its leave-one-control-out expression fit is not independent coordinate validation.
The strict outer experiment initializes from reconstruction-gene PCA and anatomy,
refines only training profiles, and maps untouched test profiles to a frozen atlas.
The training graph baseline uses PCA/DPT and nearest-neighbor projection, without
claiming to implement a frozen out-of-sample Harmony transform.

Two healthy specimens cannot support replicate-informed training plus an
independent healthy specimen holdout. Four-mouse LOSO is an AKI invariance stress
test. Agreement across injured and healthy mice can suppress real biology or
confound injury with position. No confirmatory significance is reported.

Gene and marker exclusion protect direct use of evaluation expression, but
anatomy and the saved DPT inherit marker information. Random gene holdout also
leaves correlated-gene dependence. Independent anatomy and module/block gene
holdout are required before a strong independence claim. Iterative training
updates feed back across specimens; only the untouched outer holdout supplies
prospective separation.

Replication cannot distinguish a coordinate from any monotone reparameterization.
Anatomy thirds are a gauge, not physical lengths. Hard intervals are a prototype
ceiling; the unconstrained test mapping is an anatomy ablation. Shared nuisance,
symmetric programs, missing overlap and biological warping are explicit failures.
Bootstrap ranges describe conditional stability, not calibrated spatial uncertainty.

## Geneformer branch

Direct V2-104M inference on aggregated mouse profiles is deferred: the documented
input is raw human single-cell expression with version-specific rank encoding,
medians, token dictionary and Ensembl IDs. Ortholog conversion does not validate
aggregated profiles. Underlying 2µm bins are not automatically single cells.

A documented external single-cell-derived gene prior can be supplied with:

- `PSEUDOSPACE_15_GF_PRIOR`: CSV with unique `gene` and finite `score` columns.
- `PSEUDOSPACE_15_GF_MANIFEST`: JSON with `scope="external"`, `model_revision`,
  `dataset_accession`, `tokenizer`, `importance_procedure`, `ortholog_version`.
- `PSEUDOSPACE_15_TF_LIST`: optional curated, versioned mouse TF CSV with `gene`.

Dataset-dependent priors are refused by this external-prior route; they require
fold-specific computation. A gene-prior panel is compared separately with
Geneformer-ranked genes plus DPT. Cell embedding plus DPT remains an explicitly
deferred representation baseline until defensible cell-level input is available.
Perturbation sensitivities and TF prediction do not establish causal regulation.

## Artifacts and interpretation

Inspect `go_no_go_transductive.csv`, `heldout_gene_metrics.csv`,
`specimen_metrics.csv`, `selected_gene_panels.csv`, `heldout_positions.csv`,
`baseline_failures.csv`, `synthetic_metrics.csv` and
`synthetic_paired_ordering_gains.csv`. Conditional bootstrap artifacts cover
positions, fixed-pair ranks, fixed evaluation-gene curves and transition/peak
locations. The final verdict answers all twelve methodological questions and
keeps missing biological evidence explicitly pending.

Biological GAM examples are gated on excluded-gene/within-segment prediction
and synthetic fine-order improvement. No pathway or cross-species discovery is
attempted before that gate. The fixed 5% exploratory improvement threshold is
not statistical significance. A failed gate means retaining conventional DPT,
not tuning the simulation or adding a foundation model to rescue the story.

Figures A/B/E run without private inputs; C/D/F/G require the current real-data
export. All are exploratory panels with source tables for quantitative outputs;
publication assembly and final journal QA remain separate.
