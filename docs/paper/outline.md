# Paper plan: where human and mouse proximal tubule differ along the tubule axis

Working document owned by the coordinating session. Every analysis task must map to a claim,
figure or method below. Work that does not serve this plan is out of scope until the plan
changes.

## Central claim

Measured against same-species reproducibility, **human and mouse proximal tubule agree mainly on
their most strongly zonated genes** (notebook 49). With gradient strength cross-fitted on external
donors and ceilings estimated per strength bin, the human–mouse conservation index is about 0.55 in
the top 1% of genes, 0.1–0.3 among moderately zonated genes whose gradients still reproduce within
each species, and near zero below that. Over all genes the index is 0.1–0.35 depending on
winsorisation (0.17/0.28 at the published setting; module-block intervals reach 0.54 for the S3
contrasts), and 0.25–0.31 in public atlases alone. Among genes zonated in both species, direction
agreement rises from about 0.5 to 0.85–0.89 with strength. Whether mouse-only genes outnumber
human-only ones depends on absolute versus amplitude-relative thresholds. Second contribution: with
two specimens per species, average-level pathway screens cannot be shown specific by relabeling,
whereas position-dependent screens can when specimens within a group share positional biology
(notebooks 31, 37, 45). The primary list of 22 position-dependent pathways (robust under both
coordinates; 20 also called coordinate-free) mostly summarises mouse zonation that the human
sections do not share.

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

## Results sections, evidence and status (draft v0.5 structure)

| # | Section | Evidence | Figure | Status |
|---|---|---|---|---|
| R1 | Segmented PT structures and their segment labels: cohort, joint labels, glomerulus distance, reference transfer (weak human support), spillover; the coordinate as a tool | Notebooks 01, 03, 13, 42; Supplementary Note 1 (notebooks 35, 36) | 1 | Done; lab records pending |
| R2 | Human and mouse agree mainly on their most strongly zonated genes: ceilings, conservation index and sensitivity, agreement by gradient strength, cross-fitted labels, integration genes excluded | Notebooks 43, 46, 47, 49 | 2 | Done (v0.6 wording fixes pending) |
| R3 | Which genes are shared: symmetric classes, direction agreement, transporters over-represented (33 of 196 conserved), absolute versus relative asymmetry, tissue and physiological state | Notebooks 43, 44, 47, 48 | 3 | Done |
| R4 | With two specimens per species, average-level pathway screens cannot be shown specific by relabeling; positional screens pass when specimens share positional biology; external many-draw nulls; AKI read under both rules | Notebooks 31, 37, 45, 48 | 4 | Done; depends on unrecorded mouse ages |
| R5 | Position-dependent pathways mostly summarise mouse zonation: 22 robust under both coordinates (13/3/0/1/5; 12 core), 20 also called coordinate-free; replication beyond random sets, beyond relabel-called in one of two atlases | Notebooks 37, 47, 48, 50 | 5 | Done |
| R6 | Lead genes: Gatm (amplitude), Acadm/Acaa2, Dcxr and Ugt3a1 reversals, RBP4, Slc7a13 reference; rodent protein support | Notebooks 38, 43, 44 | 6 | Done |

Figures are built by notebook 41 from saved results (`results/paper_figures/`). Supplementary Notes 1–10 and
the protocol-deviation table are in `docs/paper/supplementary_information.md`.

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

## Notebooks behind the paper

| Notebooks | Role | Results folder |
|---|---|---|
| 13, 35, 36 | Coordinate and its validation (Supplementary Note 1) | `minimal_pt_scfates/`, `pt_reconstruction_v2/` |
| 31, 37 | Pathway method selection and final joint test | `pt_pathway_method_selection/`, `pt_pathway_final/` |
| 38, 44 | Literature, external genes, state and rodent protein | `pt_literature_deep/`, `pt_revision_state/` |
| 39, 40, 43 | Amplitude, v1 classes (superseded), ceilings and symmetric classes | `pt_zonation_amplitude/`, `pt_zonation_classes/`, `pt_revision_classes/` |
| 41 | All paper figures | `paper_figures/` |
| 42, 46, 47 | Label validation, cross-fitted labels, addenda | `pt_revision_labels/`, `pt_revision_classes/crossfit/`, `pt_revision_addenda/` |
| 45, 48, 50 | Many-draw nulls, overlap-preserving replication, primary list | `pt_revision_method/`, `pt_revision_addenda_pathways/`, `pt_primary_pathways/` |
| 49 | Agreement by gradient strength and index sensitivity | `pt_conservation_strength/` |
