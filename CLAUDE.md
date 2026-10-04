# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

@AGENTS.md

The imported `AGENTS.md` holds the authoritative project rules. This file adds what it leaves out.
Like every tracked file, it must contain no machine-local absolute paths: the hygiene gate rejects them.

## Commands not covered above

```bash
pytest tests -q                                   # full synthetic suite, about 20 s
pytest tests/test_ordered_registration.py -q      # one module
pytest "tests/test_ordered_registration.py::test_partial_reverse_path_does_not_stretch_to_reference_endpoints"
python tools/agent_gates.py --fast                # staged hygiene + stage order + 06 column audit + source-divergence warning
python tools/agent_gates.py                       # adds the synthetic dry-run of 06 and pytest
python tools/notebook_dryrun.py --notebook analysis/notebooks/<nb>.ipynb   # run cells on synthetic data (default: 06)
python tools/stage_notebooks_without_outputs.py --list                     # report notebooks with outputs, change nothing
```

`agent_gates.py` rewrites the tracked file `.agent/gates.json` on every run. Do not commit that change
unless you meant to. Run `check_repository_hygiene.py` without flags and it inspects the working tree. That check fails
while notebooks hold live outputs, which is normal. The commit gate is `--staged`.

## Working-tree safety

The working notebooks usually hold executed outputs and uncommitted runs. That is why `git status`
lists every notebook as modified. Stage only the paths you changed, by name. Never run `git add -A`,
`git add analysis/notebooks`, `git checkout -- <notebook>`, `git stash`, or `nbconvert --inplace`
on them. Any of these either commits outputs or destroys results that cannot be recovered. Before
regenerating or overwriting a notebook, check that its source does not diverge from HEAD; the
`agent_gates.py` warning reports this.

## Two analysis lines

**Established pseudospace pipeline (notebooks 01–14).**
- 01 builds `tubule_by_gene/*.h5ad` from Visium HD and segmentations.
- 02 is mouse-only DPT.
- 03 builds the cross-species shared space and the PT DPT that 05, 06 and 07 consume.
- 13 produces the PT scFates coordinate. 12 builds pathway remodeling on it. 14 draws figures from
  the saved outputs of 12 and 13 and refits nothing.
- 04 and 08–11 are comparisons and validations.

Generated `.py` mirrors exist for 02–14. The mirror list in AGENTS.md omits 10, 11 and 14, but they
follow the same naming rule. Notebooks 01 and 15+ have no mirror.

**Repeated-structure ("pseudostructure") reconstruction (notebooks 15–30).** This is the active
research line. It aims to recover one canonical ordered PT distribution path from unordered
cross-sections of repeated physical nephrons. The canonical notebook 15 is
`15_spatial_probability_flow_proof_of_concept.ipynb`. Each notebook is a self-contained experiment
with the same skeleton:
- It locates `REPO` by walking up to `pseudospace/` and adds it to `sys.path`.
- It takes `--data-root`/`--results-root` via `parse_known_args` with the `PSEUDOSPACE_*` env
  fallback.
- It calls `pseudospace.repeated_inputs.load_repeated_mouse_inputs(DATA_ROOT, samples=CONTROLS)` on
  the healthy controls `Ctrl1A2`/`Ctrl1A4`.
- Explicit input guards follow: S1/S2/S3 labels only, finite expression, `n_spots`.
- It fixes sha256 gene folds and excludes markers.
- It sets `NOTEBOOK_LOGIC_VERSION = 'NN.name.k'` and uses `stage_cache.cached_payload`.
- It writes outputs under `RESULTS_ROOT/<experiment_name>/`.

New reusable logic goes in a `pseudospace/<module>.py` with a matching `tests/test_<module>.py`.
An experiment commit (see `git show --stat HEAD`) normally touches the notebook, its module and
tests, and four documents:
- `docs/results/pseudostructure-current-state.md` is the short entry point: goal, strongest
  components, a **closed-hypotheses table**, and the next checkpoint. Read it first.
- `docs/results/pseudostructure-research-ledger.md` records each **protocol before any private-data
  fit**, then the measured results with their logic version and cache key.
- `docs/workflows/repeated-structure-pseudospace.md` covers the model and its identification limits.
- `docs/research/pseudostructure-methods-survey.md` holds methodological precedents.

Rules for this line (from the current-state doc):
- Before each experiment, name the reconstruction failure it addresses and how the repeated-object
  prior enters inference. Freeze the protocol in the ledger and state the result that would justify
  keeping it.
- Do not rerun closed variants without a new reconstruction-specific hypothesis. Closed so far:
  initialization, mean-path registration, covariance, offset, and flow variants.
- Compare against the strong source scFates/DPT neighbor references plus segment and coverage
  oracles. DPT/scFates are comparators, not anatomical truth. Depth never selects or validates a
  model.
- Record null and negative results in the ledger instead of tuning them away. Keep injury, species
  and regulator interpretation gated until healthy transfer is validated.
