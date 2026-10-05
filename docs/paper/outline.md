# Paper plan: where human and mouse proximal tubule differ along the tubule axis

Working document owned by the coordinating session. Every analysis task must map to a claim,
figure or method below. Work that does not serve this plan is out of scope until the plan
changes.

## Central claim

Reconstructing a continuous proximal-tubule (PT) coordinate from segmented Visium HD tubules lets a
human–mouse comparison ask *where* along the tubule a program differs. With two specimens per
species, average-level pathway comparisons (whole-PT or S1/S2/S3 pseudobulk) cannot be separated
from specimen variation; position-dependent comparisons can when specimens of a group share their positional biology
(relabeling control; many-draw external null, notebook 45; the split-plot argument is motivation, not
proof). Applied to human and mouse PT, they show that
**axial zonation is largely species-specific beyond a conserved core of strong markers**: against
same-species reliability ceilings, the human–mouse correlation of gene gradients is only 0.17
(S1→S2) and 0.28 (S3) of the ceiling (notebook 43). Confidently species-only genes are few and not
skewed toward mouse (41 mouse-only, 76 human-only); human amplitudes are also lower overall. The
pathway screen mainly summarizes mouse zonation that human PT does not share (26 of 36 robust
pathways). The earlier "mouse zonates more programs" reading was a threshold artefact and is
withdrawn; rat data fail the reliability rule and are not used.

> **Status after the revision round (notebooks 42, 44, 45, 46).** The central claim survived its tests.
> Human segment labels have only weak reference support (transfer agreement 0.71), but with transferred
> labels and held-out genes the conservation index stays ≤ 0.31 (upper bound 0.38; notebook 46). Our
> donor sits at the median of healthy KPMP donors on a cross-platform injury score; time of day does not
> create the lead mouse gradients, but fasting is not excluded for sterol and PPARα programs (notebook
> 44). R2 holds in external many-draw nulls at segment resolution but is conditional: positional
> screens lose specificity when specimens differ in positional biology (age, injury severity; notebook 45).
>
> **Decisions after the independent review of v0.3 (user, 2026-10-04;
> `docs/paper/referee_report_v03.md`).**
> - **Framing: biology-led, with the pathway-specificity method as a full second contribution.** Lead with
>   the ceiling-calibrated human–mouse comparison and agreement by gradient strength; the few-specimen
>   pathway-specificity result keeps its own section and figure; the continuous coordinate is a tool
>   (one Figure 1 panel plus a supplementary note) that adds within-specimen replication, not finer
>   resolution.
> - **Primary pathway list: the 22 pathways robust under both scFates and DPT13**; peaks reported in
>   segment units; curves still drawn on scFates; each coordinate reported separately as a sensitivity.
>   This supersedes the earlier "keep scFates (36)" decision. Disclose that the pre-registered rule
>   selected DPT13 and that the override reason did not hold on the joint test.
> - Restructure toward ~4,500 words of Results and six main figures (referee report section 4).

## Results sections, evidence and status

| # | Section | Claim | Evidence | Status / owner |
|---|---|---|---|---|
| R1 | Reconstructing the PT axis | The scFates coordinate recovers directly measured zonation (coordinate-based P1: mouse 0.745, human 0.25), absorbs no injected species-by-position effect and survives an independent count split (59/66 robust kept). Its fine within-segment order is weakly supported, it predicts held-out genes no better than segment labels, it is depth-dependent in human S1 and its registration is fragile. Its value is replication across positions within specimens and continuous localization, not finer-than-segment resolution | Notebooks 13, 35, 36 | **Decided (user): scFates primary.** The pre-registered rule selected DPT13 (P5 0.80 vs 0.72; registration gap 0.14 vs 0.29); the override cited NCDR 0.40 from the matched-only test, and on the joint test DPT13 is specific too. Disclose both. DPT13 and the anchor coordinate are sensitivity analyses; 22/36 robust pathways are also DPT13-robust. Human segment labels: weak reference support (notebook 42), claim robust to transferred labels (notebook 46) |
| R2 | A specificity problem in cross-species pathway analysis | Average-level competitive screens report as many pathways for relabeled specimen groups as for species (NCDR ≥ 1.4; decoy FDP ≥ 0.75); position-dependent screens are specific (joint `T_spatial`: 60 versus 2 and 0). External many-draw nulls: whole-PT GSEA NCDR 1.34, segment-step joint screen 0.05; specificity is empirical and fails when specimens differ in positional biology (age, AKI severity). Split-plot argument is motivation only. NCDR measures vulnerability to specimen noise, not FDP | Notebooks 31, 37, 45 (Figure 2); `docs/paper/r2_revision_results.md` | Done (workstreams 2, 5) |
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
