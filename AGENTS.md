# AGENTS.md

Private, reproducible kidney-spatial research workspace: reconstructs a nephron "pseudospace"
(position along the nephron) from Visium HD tubule aggregates, plus the H&E panoptic segmentation
subproject that produces its GeoJSON inputs.

## Project

- Two connected workflows: `analysis/` (pseudospace) and `segmentation/` (boundary-based tubule /
  glomerulus panoptic segmentation, package `kidney_panoptic`).
- Python 3.11 + Scanpy/AnnData for analysis; R `harmony` **2.0.5** via rpy2 (the workflow refuses
  older user-level R installs). The `segmentation/` project has its own, deliberately separate
  Python 3.10 / Torch environment.
- **No private data in Git.** No slides, Visium matrices, H5AD, GeoJSON, checkpoints, or executed
  notebooks — enforced by `tools/check_repository_hygiene.py`. Real data lives outside the repo and
  is passed in by path; see [data/README.md](data/README.md) and [docs/data/access.md](docs/data/access.md).
- Entry points: the notebooks in `analysis/notebooks/`, with generated `.py` mirrors at
  `analysis/mouse_only_pseudospace.py` (active mouse-only workflow) and
  `analysis/human_vs_healthy_mouse.py` (cross-species).
- **Cohort**: 2 control mouse specimens, 2 AKI mouse specimens, and 2 human kidney slices. The two
  human slices are *named* cortex and medulla, but **both are in reality healthy human cortex** —
  do not treat the medulla-labelled slice as medullary tissue in any analysis or write-up.
- Scientific status: the 2-vs-2 mouse comparison is descriptive, **not** confirmatory. Read
  [docs/results/current-mouse-run.md](docs/results/current-mouse-run.md) before quoting any ranked
  gene/pathway list.

## Commands

```bash
# Analysis environment (repo root)
conda env create -f environment.yml && conda activate kidney-pseudospace

# Synthetic regression tests (no private data; needs anndata/scanpy/matplotlib)
pytest tests

# Repository hygiene gate (requires nbformat) — run before every commit
python tools/check_repository_hygiene.py

# Active workflows (flags win over env vars)
python analysis/mouse_only_pseudospace.py   --data-root <dir> --results-root <dir>
python analysis/human_vs_healthy_mouse.py   --data-root <dir> --results-root <dir>
# env-var equivalents: PSEUDOSPACE_DATA_ROOT, PSEUDOSPACE_RESULTS_ROOT

# QuPath export / spatial validation
python analysis/scripts/export_kept_tubules_to_qupath_geojson.py --data-root <dir> --results-root <dir> [--dry-run]
python analysis/scripts/spatial_validation_of_pseudospace.py      --data-root <dir> --results-root <dir> [--dry-run]

# Segmentation subproject (independent env)
conda env create -f segmentation/environment.yml && conda activate kidney-panoptic
PYTHONPATH=segmentation/src pytest segmentation/tests -m "not slow"
```

CI (`.github/workflows/ci.yml`) runs three jobs: `hygiene`, `pseudospace-synthetic`
(`pip install -e . && pytest tests`), and `segmentation-unit` (CPU-only, `-m "not slow"`).

## Architecture

| Path | Role |
| --- | --- |
| `pseudospace/` | Reusable, data-location-agnostic analysis logic: `markers` (alias resolution, marker axis, DE cluster annotation), `harmony` (R integration + version guard), `trajectory` (DPT root/orientation), `levelshape` + `stats_gam` (level/shape GAM decomposition, LOSO stability), `cross_species`, `heatmaps`, `io_qc`, `scprisma_pseudospace`. |
| `analysis/` | Canonical `notebooks/*.ipynb` workflows plus their generated `.py` mirrors and the study-specific data contract. |
| `analysis/scripts/` | QuPath label export and spatial validation of the pseudospace. |
| `segmentation/` | Self-contained `kidney_panoptic` project (`src/kidney_panoptic/{data,models,losses,postprocess,eval,utils,infer}`, `scripts/`, `configs/`). Frozen encoder + native-resolution decoder + watershed decode. |
| `docs/` | `architecture/000N-*.md` decision records, `workflows/` run guides, `results/` curated tables, `publication-readiness.md`. |
| `tools/` | `check_repository_hygiene.py` — the pre-commit gate. |
| `legacy/` | Superseded notebooks, provenance only. Unsupported; never mix with current results. |

Data flow: private Visium HD + v4 GeoJSON → `notebooks/01_segmentation_to_gene_matrix.ipynb` →
private tubule-by-gene H5AD → the mouse-only workflow → private Harmony/DPT outputs + curated
tables under `docs/results/`.

## Conventions

- **Never hardcode user-specific paths.** Scripts take `--data-root` / `--results-root` and fall back
  to the `PSEUDOSPACE_*` env vars; CLI arguments win. The hygiene checker rejects any tracked file
  containing a machine-local absolute home directory (never spell one out in tracked text — that is
  what makes the checker fail).
- **Notebooks are the canonical artifact; the `.py` files are generated mirrors.** You (the human)
  execute `analysis/notebooks/*.ipynb` — that is where the workflow lives, and running cells is
  always fine. The `.py` exists only as a token-efficient plain-text surface for agents. The **agent
  edits only the `.py`**, never the `.ipynb`. Regeneration is one-way, notebook → script (see ## Notes),
  so a `.py` edit survives only until the next regeneration: anything durable must reach the notebook.
- **`comment_magics = true`** (`jupytext.toml`): IPython magics appear as comments
  (`# %load_ext autoreload`) in the mirror, keeping the `.py` valid, parseable Python that agents can
  syntax-check. jupytext restores them as live magics in the notebook, so the notebook is unaffected.
- **Committed notebooks must be output-free.** Strip before committing with
  `python -m nbconvert --clear-output --inplace <notebook>`; tracked cell outputs or execution
  counts fail the hygiene gate.
- **Tests must not need private data.** Use synthetic fixtures; gate heavy imports with
  `pytest.importorskip(...)` (see `tests/test_pseudospace_synthetic.py`).
- **Scientific guardrails** (do not bypass):
  - Only `*_v4.geojson` mouse segmentations are valid; centroid verification is mandatory
    (`analysis/`, `docs/workflows/pseudospace.md`).
  - The coarse-label checkpoint requires explicit per-cluster confirmation, and the reference
    13-cluster fingerprint must be updated together with the labels. The workflow stops before DPT
    if it drifts — never bypass the guard.
  - The canonical pseudospace coordinate is Scanpy DPT on the pass-2 Harmony embedding, rooted in PT
    and oriented by the early→late marker axis.
  - Model selection in `segmentation/` is pooled class-agnostic PQ (not validation loss); D4 TTA is
    required for reported inference.
- **The public `obs` label contract** consumed by the package, QuPath export, and spatial
  validation: `segment_class` (fine), `coarse_class` (rollup), `broad_tubule_marker_call` (alias
  expected downstream).
- Keep scientific caveats beside result summaries, not in a separate file.
- `pseudospace/` modules are extracted from the notebooks and keep notebook-equivalent code;
  prefer adding reusable logic there and thin orchestration in `analysis/`.
- **Commit and push after every change.** Version control is the safety net here: once an edit is
  verified, commit it and push to `origin/main`
  (`https://github.com/edmund-code/pseudostructure-reconstruction.git`) rather than leaving work
  uncommitted — and never push private data or executed notebooks.

## Notes

- **Regenerating the mirrors is a manual, one-way step — the pairs are not auto-wired.** `jupytext.toml`
  declares a default `ipynb,py:percent` pair, but `jupytext --paired-paths` resolves to files that do
  not exist (the notebooks sit in `analysis/notebooks/` behind a `NN_` prefix; the scripts are
  `analysis/<name>.py`), so jupytext cannot infer the pairs and `jupytext --sync` has no partner.
  Pairing:
  `02_mouse_only_pseudospace.ipynb → analysis/mouse_only_pseudospace.py` and
  `03_human_vs_healthy_mouse.ipynb → analysis/human_vs_healthy_mouse.py`.

  ```bash
  jupytext --to py:percent --output analysis/mouse_only_pseudospace.py \
    analysis/notebooks/02_mouse_only_pseudospace.ipynb
  jupytext --to py:percent --output analysis/human_vs_healthy_mouse.py \
    analysis/notebooks/03_human_vs_healthy_mouse.ipynb
  ```

  Add `JUPYTER_DATA_DIR=/tmp/jupyter-data` when nbformat cannot write its signature secret file.
  Then confirm both compile: `python -m py_compile <script>`.
- **Notebook 02 owns the `_workflow_roots` helper** (`argparse` + `PSEUDOSPACE_*` fallback) that keeps
  `--data-root` / `--results-root` working in the regenerated mouse script. Restore it in the
  notebook, never in the generated `.py`, or those flags silently disappear from the script.

<!-- Leave a clean house: no stray scratch files, temp scripts, or vendored archives at the repo
     root. Anything that isn't project source belongs outside the workspace. -->
