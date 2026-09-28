# Kidney pseudospace reconstruction

Private, reproducible research workspace for two connected kidney-spatial workflows:

1. **Pseudospace analysis** reconstructs nephron position from Visium HD tubule aggregates and
   describes control-versus-AKI and human-versus-healthy-mouse changes along that coordinate.
2. **Panoptic segmentation** produces the v4 tubule/glomerulus GeoJSON segmentations consumed by
   the analysis workflow.

The repository deliberately contains no raw slides, Visium HD matrices, H5AD files, GeoJSON
annotations, or model checkpoints. See [data/README.md](data/README.md) for the private-data
contract and [docs/data/access.md](docs/data/access.md) for setup.

## Repository map

| Directory | Purpose |
| --- | --- |
| `analysis/` | Current preprocessing, mouse-only pseudospace workflow, and QuPath utilities. |
| `pseudospace/` | Reusable analysis package. |
| `segmentation/` | Boundary-based `kidney_panoptic` training and inference project. |
| `docs/` | Architecture, workflow guides, data access, current results, and scientific caveats. |
| `examples/`, `tests/` | Synthetic fixtures and no-private-data regression tests. |
| `legacy/` | Superseded notebooks retained for provenance only. |

## Quick start

Create the analysis environment, point it at authorized data, and run the active workflow:

```bash
conda env create -f environment.yml
conda activate kidney-pseudospace
python analysis/mouse_only_pseudospace.py \
  --data-root /path/to/private/pseudospace-data \
  --results-root /path/to/private/pseudospace-results
```

The workflow intentionally pauses at the documented cluster-label review checkpoint. It is not a
one-command black box. Read [docs/workflows/pseudospace.md](docs/workflows/pseudospace.md) first.

For the cross-species comparison, use:

```bash
python analysis/human_vs_healthy_mouse.py \
  --data-root /path/to/private/pseudospace-data \
  --results-root /path/to/private/pseudospace-results
```

Read [the human-versus-healthy-mouse workflow guide](docs/workflows/human_vs_healthy_mouse.md)
before interpreting its descriptive curves.

For the PT pathway workflow, run the two-pass nephron/scFates producer and then notebook 12:

```bash
python analysis/minimal_pt_scfates.py --data-root /path/to/private/pseudospace-data --results-root /path/to/private/pseudospace-results
python analysis/pt_pathway_remodeling.py --data-root /path/to/private/pseudospace-data --results-root /path/to/private/pseudospace-results
```

The first step repeats notebook 03's reviewed clustering, nephron filtering, and second Harmony
pass. It compares nonbranching scFates with DPT on the full nephron and PT subset, then writes
notebook 12's PT coordinate and ortholog-map inputs under `minimal_pt_scfates/`.

For synthetic regression checks that need no research data:

```bash
pytest tests
```

The segmentation project has an independent environment and workflow; see
[segmentation/README.md](segmentation/README.md).

## Scientific status

The current mouse run is descriptive (two control and two AKI specimens), not confirmatory. Its
proximal-tubule ordering is well supported, whereas the broader nephron axis has known weak
regions. Start with [the current results interpretation](docs/results/current-mouse-run.md)
before using any ranked gene or pathway list.

## Citation and access

Use [CITATION.cff](CITATION.cff) for internal citation metadata. This private repository has no
open-source license yet; do not redistribute code or data without maintainer approval.
