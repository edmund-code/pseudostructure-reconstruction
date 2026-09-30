# AGENTS.md

Private, reproducible kidney-spatial research workspace. `analysis/` reconstructs nephron
"pseudospace" from Visium HD tubule aggregates; `segmentation/` provides H&E panoptic segmentation
and its GeoJSON inputs.

## Project

- Analysis uses Python 3.11, Scanpy/AnnData, and R `harmony` **2.0.5** via rpy2; older user-level R
  installs are refused. Segmentation has a separate Python 3.10/Torch environment.
- **Never put private data in Git** (slides, Visium matrices, H5AD, GeoJSON, checkpoints, or executed
  notebooks). Data stays outside the repo and is passed by path; see [data/README.md](data/README.md)
  and [docs/data/access.md](docs/data/access.md). `tools/check_repository_hygiene.py` enforces this.
- Canonical workflows are notebooks in `analysis/notebooks/`; generated mirrors include mouse-only
  (`analysis/mouse_only_pseudospace.py`), cross-species 03, GAM 05, spatial rewiring 06, gene modules
  07, scFates 13, and pathway remodeling 12. Notebook 13 compares nonbranching scFates and DPT on
  nephron and PT scopes; its PT scFates coordinate feeds 12. GAM 05 reads 03's PT pseudospace and
  recomputes no coordinate. Notebook 12 uses gene-first nested models on 13's PT scFates coordinate.
- **Cohort:** 2 control mouse, 2 AKI mouse, and 2 human kidney slices. The human slices are named
  cortex and medulla, but **both are healthy human cortex**; never treat the medulla-labelled slice
  as medullary tissue.
- The 2-vs-2 mouse comparison is descriptive, **not confirmatory**. Read
  [docs/results/current-mouse-run.md](docs/results/current-mouse-run.md) before quoting ranked
  genes or pathways.

## Commands and architecture

Create the analysis environment with `conda env create -f environment.yml`; activate
`kidney-pseudospace`. Run synthetic tests with `pytest tests`. Before commits, run the repository
hygiene check and stage output-free notebook copies with `python tools/stage_notebooks_without_outputs.py`;
after notebook edits, run `python tools/check_notebook_stage_order.py analysis/notebooks/*.ipynb`.
Workflow scripts accept `--data-root` and `--results-root` (CLI flags override
`PSEUDOSPACE_DATA_ROOT` and `PSEUDOSPACE_RESULTS_ROOT`). Segmentation uses
`conda env create -f segmentation/environment.yml`, then
`PYTHONPATH=segmentation/src pytest segmentation/tests -m "not slow"`. CI runs hygiene, synthetic
pseudospace, and CPU segmentation-unit jobs.

`pseudospace/` contains reusable marker, Harmony, trajectory, GAM/statistics, module, enrichment,
pathway, specimen, cross-species, heatmap, QC, cache, and scPrisma logic. `analysis/` holds workflow
notebooks/mirrors and QuPath/spatial-validation scripts. `segmentation/` is the self-contained
`kidney_panoptic` package. `docs/` holds architecture decisions, workflow guides, results, and
publication readiness. `legacy/` is provenance only; never mix it with current results.

## Conventions and scientific guardrails

- Never hardcode machine-specific paths. Scripts use the root flags/env vars above; hygiene rejects
  tracked machine-local absolute home paths.
- Notebooks are canonical and `.py` files are generated, one-way mirrors. Make durable code changes
  in `.ipynb`, then regenerate the mirror; `.py`-only changes are discarded. `comment_magics = true`
  in `jupytext.toml` keeps mirrors parseable while restoring magics in notebooks.
- Run the notebook stage-order gate before mirror regeneration. It catches use-before-definition and
  names assigned only in a conditional branch then read later (a latent `NameError`).
- Committed notebooks must be output-free, but never clear a working notebook to achieve that.
  `stage_notebooks_without_outputs.py` stages a cleaned copy without changing the working file; do
  not `git add` those notebook paths afterward. Check the index with
  `python tools/check_repository_hygiene.py --staged`. Clearing a live notebook in place destroys
  outputs that cannot be recovered.
- Tests must not need private data; use synthetic fixtures and `pytest.importorskip(...)` for heavy
  imports.
- Cross-species shared space contains only genes measured in every input. The accepted ortholog map
  creates a target column for every pair, so an unfed column is a structural zero, not an observation
  of zero expression; a measured zero stays. `combine_cross_species` records
  `measured_in_both_inputs` and `uns['cross_species_availability']` (03 also writes
  `diagnostics/cross_species_*.csv`). Notebook 03 keeps every gene in the object
  (`require_measured_in_both=False`) and masks the analysis gene set instead: removing genes changes
  Scanpy seurat-flavour HVG binning by mean expression, which changes PCA/Harmony, moves the Leiden
  partition, and invalidates the hand-reviewed cluster labels.
- Set-level enrichment is exploratory, never confirmatory. `camera_like_enrichment` contrasts a
  set's mean statistic with the background mean, corrects across all tested pairs (pooled and
  specimen-balanced rankings together), and estimates variance inflation from residual correlations
  given the fitted design. With two mice and one human donor, specimens are the replication units.
- Only `*_kept_tubules_labeled_fine.geojson` mouse segmentations are valid. They are quality
  controlled upstream, so mouse workflows do no tubule-level **gene-count** QC. They apply a
  structure filter: a per-species `n_spots` quantile floor (`MIN_SPOTS_PER_SPECIES_QUANTILE = 0.05`
  in 03, `MIN_SPOTS_QUANTILE` in 02). Use the quantile floor, not an absolute threshold, to preserve
  species balance and nephron ordering. Centroid verification is mandatory; see
  `docs/workflows/pseudospace.md`.
- Coarse-label checkpoint requires explicit per-cluster confirmation; update the reference
  13-cluster fingerprint with the labels. Never bypass the stop before DPT if it drifts.
- Original 02/03 coordinate is Scanpy DPT on pass-2 Harmony, rooted in PT and oriented by the
  early→late marker axis. Notebook 13 repeats two-pass integration and compares nonbranching scFates
  with DPT on nephron and PT subsets; its separate PT scFates coordinate feeds 12.
- Notebook 06 is the **frozen discovery notebook** for human vs healthy-mouse PT spatial analysis. It
  uses 03's shared PT DPT exactly as 03 oriented it (no cross-species re-registration); displacement
  is a supporting measurement, not a phenotype class. Validate leading pathway candidates in a
  separate analysis. It caches its own per-specimen fits and cross-checks against 05's saved curves
  and atlas without consuming them. Notebook 07 discovers gene-level normal-positional and human-vs-mouse response
  modules from 03's PT object and ortholog map; load pathway annotations only after cluster assignments
  and cluster-number decisions are frozen.
- In `segmentation/`, select models by pooled class-agnostic PQ (not validation loss); D4 TTA is
  required for reported inference.
- Public `obs` label contract: `segment_class` (fine), `coarse_class` (rollup),
  `broad_tubule_marker_call` (downstream alias). Keep scientific caveats beside result summaries.
- Prefer reusable logic in `pseudospace/` with thin orchestration in `analysis/`.
- Expensive stages cache by parameters, inputs, and code in `results/<workflow>/stage_cache/`.
  Unchanged notebooks print `[stage cache] hit`; `PSEUDOSPACE_STAGE_CACHE=0` forces rebuild, and
  `pseudospace.stage_cache.cache_status(...)` / `purge_stage_cache(...)` inspect or clear cache.
  Bump `NOTEBOOK_LOGIC_VERSION` in the config cell when a cached cell changes. Never cache validation
  or guard cells (Harmony version, input existence, fingerprint checks).
- **Commit and push every verified change to `origin/main`**; never push private data or executed
  notebooks. Follow the output-free notebook convention above.

## Mirror regeneration

Pairs are manual: for notebook IDs 02–09, 12, and 13, remove the `NN_` prefix and change `.ipynb`
to `.py` under `analysis/`. `jupytext --sync` cannot infer these partners. Regenerate explicitly:

```bash
jupytext --to py:percent --output <mirror> analysis/notebooks/<notebook>
python -m py_compile <mirror>
```

Set `JUPYTER_DATA_DIR=/tmp/jupyter-data` if nbformat cannot write its signature secret file.
Notebook 01 (`01_segmentation_to_gene_matrix.ipynb`) has no mirror. Notebook 02 owns
`_workflow_roots` (`argparse` plus `PSEUDOSPACE_*` fallback); restore it in the notebook, never the
generated `.py`, or regenerated mouse script root flags disappear.

<!-- Keep the repository root free of scratch files, temp scripts, and vendored archives. -->

## Implementation discipline

Trace the real flow and callers before editing; fix shared root causes. Prefer no new code,
existing helpers, stdlib, native features, installed dependencies, then the smallest working change.
Avoid unrequested abstractions, boilerplate, and dependencies. Question unnecessary complexity;
prefer deletion and edge-case-correct solutions. Do not cut corners on required behavior, trust-boundary
validation, data-loss prevention, security, accessibility, or hardware calibration. Comment deliberate
simplifications with their ceiling and upgrade path (a "ponytail"). Nontrivial logic needs one small
runnable regression check; trivial one-liners need none.
