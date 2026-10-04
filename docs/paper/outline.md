# Paper plan: where human and mouse proximal tubule differ along the tubule axis

Working document owned by the coordinating session. Every analysis task must map to a claim,
figure or method below. Work that does not serve this plan is out of scope until the plan
changes.

## Central claim

Reconstructing a continuous proximal-tubule (PT) coordinate from segmented Visium HD tubules
lets a human–mouse comparison ask *where* along the tubule a program differs. With two
specimens per species, conventional pathway comparisons (whole-PT or S1/S2/S3 pseudobulk)
cannot be distinguished from specimen variation. Position-dependent differences can, and they
replicate in independent segment-resolved data.

## Results sections, evidence and status

| # | Section | Claim | Evidence | Status / owner |
|---|---|---|---|---|
| R1 | Reconstructing the PT axis | A continuous, cross-species PT coordinate recovers directly measured S1→S3 zonation in both species and does not absorb species differences | Notebook 13 (scFates). Notebook 34 validates segment labels only; coordinate-based evidence is pending (workstream 1, notebook 35) | **Open: evidence and method choice.** Workstream 1: validation layer for scFates first, then a conserved-anchor alternative |
| R2 | A specificity problem in cross-species pathway analysis | Average-level pathway screens report as many pathways for relabeled specimen groups as for species; only smooth position-dependent screens are specific | Notebook 31 (Figure 1, within-species controls) | Done; workstream 2 audits and finalizes |
| R3 | Position-dependent pathway differences | 82 confident → 66 robust → 42 core pathways; human-high programs peak early, mouse-high late | Notebooks 31–33 (atlas, funnel, landscape) | Done; workstream 2 finalizes the reporting set and figure |
| R4 | Robustness and independent replication | Results survive registration, region, detection, sex and trajectory-method checks and replicate in external segment data | Notebooks 33, 34 | Done; extend if workstream 1 changes the coordinate |
| R5 | Biology | Lead stories: creatine-synthesis zonation (human side new), mitochondrial β-oxidation (Acadm/Acaa2 opposite in mouse, flat in human), drug-handling placement (human front-loaded conjugation/oxidation versus mouse late glutathione and S1–S2 Cyp2e1), Dcxr reversal. Supporting: Rbp4, sterol synthesis, human late-PT retention of SLC6A19/SLC9A3/SLC4A4 | Notebooks 32–34, 38; Lake/KPMP human S1/S2/S3; `docs/results/pt-pathway-literature-novelty.md`; draft `docs/paper/r5_draft.md` | Evidence done (notebook 38). Text conditional on workstream 1 attenuation check |

## Figures (draft)

1. Study design, segmentation, reconstructed coordinate, external zonation concordance (R1).
2. Specificity of pathway strategies under specimen relabeling (R2).
3. Positional atlas of species pathway differences and the early/late landscape (R3).
4. Robustness and external replication (R4).
5. Biology panels: gene curves with specimen means and conventional fold changes (R5).

## Methods (to be written from the notebooks)

Segmentation (`segmentation/`), tubule-by-gene matrices (01), cross-species integration (03,
13), coordinate reconstruction (13 or workstream 1), nested GAMs and `T_spatial` (12), matched
rank-AUC tests (12), specimen relabeling control (31), robustness checks (33), external
segment datasets (34), literature-check procedure (literature-novelty doc).

## Standing caveats (must appear beside results)

Two control mice and two cortex sections from one male human donor; all specimens male;
human sections are cortex only, so mouse outer-stripe S3 has no human counterpart; Visium HD
probe panels differ between species; findings are descriptive.

## Definition of done

- Every figure is produced by a committed, output-free notebook from private inputs.
- R1 method decision made, with pre-specified criteria and literature/first-principles
  motivation.
- R5 novelty claims each carry a literature status (verified sources) and an
  independent-data check.
- A complete manuscript draft (`docs/paper/draft.md`) with Methods, Results, Discussion and
  figure legends, consistent with the results docs.

## Workstreams (current)

| Workstream | Goal | Notebook numbers | Results folder |
|---|---|---|---|
| 1 Reconstruction | Problem-specific PT reconstruction, or a justified decision to keep scFates | 35–36 | `results/pt_reconstruction_v2/` |
| 2 Pathway analysis | Final pathway method, reporting set and figure; audit 31–34 | 37 | `results/pt_pathway_final/` |
| 3 Literature and novelty | Verified novelty claims, deeper evidence for R5 stories | 38 | `results/pt_literature_deep/` |
