# Paper plan: where human and mouse proximal tubule differ along the tubule axis

Working document owned by the coordinating session. Every analysis task must map to a claim,
figure or method below. Work that does not serve this plan is out of scope until the plan
changes.

## Central claim

Reconstructing a continuous proximal-tubule (PT) coordinate from segmented Visium HD tubules lets a
human–mouse comparison ask *where* along the tubule a program differs. With two specimens per
species, average-level pathway comparisons (whole-PT or S1/S2/S3 pseudobulk) cannot be separated
from specimen variation; position-dependent comparisons can (relabeling control, split-plot
rationale; many-draw null pending, workstream 5). Applied to human and mouse PT, they show that
**axial zonation is largely species-specific beyond a conserved core of strong markers**: against
same-species reliability ceilings, the human–mouse correlation of gene gradients is only 0.17
(S1→S2) and 0.28 (S3) of the ceiling (notebook 43). Confidently species-only genes are few and not
skewed toward mouse (41 mouse-only, 76 human-only); human amplitudes are also lower overall. The
pathway screen mainly summarizes mouse zonation that human PT does not share (26 of 36 robust
pathways). The earlier "mouse zonates more programs" reading was a threshold artefact and is
withdrawn; rat data fail the reliability rule and are not used.

> **Status after internal referee review (v0.1):** major revision. The central claim is under test
> (`docs/paper/referee_report_v01.md`, issues M1–M8): human segment labels need anatomical validation,
> the cross-species correlation needs a same-species ceiling, one-species "flat" calls need external
> confirmation, and state confounders must be bounded. If conservation proves moderate rather than low,
> the claim becomes "zonation is partly conserved, with transport conserved and metabolic programs
> differing", which the paper can still support.

## Results sections, evidence and status

| # | Section | Claim | Evidence | Status / owner |
|---|---|---|---|---|
| R1 | Reconstructing the PT axis | The scFates coordinate recovers directly measured zonation (coordinate-based P1: mouse 0.745, human 0.25), absorbs no injected species-by-position effect and survives an independent count split (59/66 robust kept). Its fine within-segment order is weakly supported, it predicts held-out genes no better than segment labels, it is depth-dependent in human S1 and its registration is fragile. Its value is replication across positions within specimens and continuous localization, not finer-than-segment resolution | Notebooks 13, 35, 36 | **Decided: scFates primary**; DPT and conserved-anchor coordinates as sensitivity (DPT fails the downstream specificity gate, NCDR 0.40; anchor fails its adoption rule). Coordinator applied the pre-registered downstream gate to all candidates after seeing DPT's result; disclose. R1 text: workstream 1 |
| R2 | A specificity problem in cross-species pathway analysis | Average-level competitive screens report as many pathways for relabeled specimen groups as for species (NCDR ≥ 1.4; decoy FDP ≥ 0.75); position-dependent screens are specific (joint `T_spatial`: 60 versus 2 and 0). Split-plot argument: offsets have 2 df of specimen error, position effects 2(k−1) | Notebooks 31, 37 (Figure 2) | Done (workstream 2) |
| R3 | Position-dependent pathway differences | Joint test: 60 → 51 confident → **36 robust (17 programs)** → 26 core; human-high programs peak early (median 0.16), mouse-high late (0.63); differences mostly graded and one-signed | Notebook 37 (Figure 3); 31–33 superseded list as sensitivity | Done (workstream 2) |
| R4 | Robustness and independent replication | Robust pathways keep ≥ 70% of effect under all checks; specimen-level split-plot test calibrated and concordant (complement). External: measurement-level agreement; robust ≫ uncalled pathways in donor-level replication, but pathway-by-pathway replication not shown (17/36 Census, 6/36 Lake), and little beyond a shared global flattening of human gradients | Notebooks 33, 34, 37 (Figure 4) | Done; rerun notebook 37 if workstream 1 changes the coordinate |
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
