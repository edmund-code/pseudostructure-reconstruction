# Figure plan: pseudostructure reconstruction paper

This plan follows the authors' manuscript outline (2026-10-05). It replaces `figure_plan.md`, which
followed the biology-led draft v0.6.

The PT ordering step is **scFates**. A new ordering model replaces it only if it passes its
pre-registered adoption criteria (notebook 52, running).

**Status codes**
- **ready**: the panel can be built from saved results.
- **rebuild**: the content exists, but the panel must be redrawn for this plan.
- **new**: needs a new analysis or render.
- **blocked**: waits on a decision or on data.

## Design rules for every figure

- One question per figure. The title states the answer.
- Show continuous and discrete views side by side whenever the claim is "beyond clustering":
  whole-PT bar, then S1/S2/S3 bars, then the continuous curve.
- Show genes, not pathway z-scores. Use gene heatmaps ordered by peak position, z-scored within
  species, with one strip per specimen so replication is visible.
- Give positions in segment units (S1 0–1, S2 1–2, S3 2–3), or on the coordinate with the
  per-specimen S1→S2 and S2→S3 transitions marked.
- Fixed colours: mouse blue, human orange; segments S1/S2/S3 on a grey ramp; a diverging map with a
  grey midpoint for z. Status colours are never reused for series.
- Every number in a panel is read from a saved table and asserted (notebook 41 convention).

---

## Figure 1: From 2D sections to a continuous pseudospace (outline: Overview)

| Panel | Content | Source | Status |
|---|---|---|---|
| a | Schematic: 3D nephron → section of disconnected cross-sections → segmentation → tubule-by-gene matrix → ordering on one pseudospace axis | Illustration | new (authors) |
| b | H&E tile with the instance-segmentation overlay (PT, DT, CD, glomerulus) and a zoomed inset | `segmentation/` inference outputs | new (render) |
| c | Segmentation performance: held-out PQ 0.837 (0.880 with D4 test-time augmentation; 98 objects), per class | Segmentation repository | rebuild |
| d | Visium HD 2 µm bins assigned to one tubule polygon (zoom), with distributions of bins and UMI per tubule by species | Notebooks 01, 42 | new (render) |
| e | Integrated tubule embedding coloured by class and specimen; PT subset with the scFates curve | Notebook 13 | rebuild (notebook 41 Fig 1b) |
| f | **Spatial map of pseudospace:** PT cross-sections on the tissue coloured by their coordinate, for one mouse and one human section. Physically scattered tubules land on one axis. x/y is used for display only | Notebook 13 coordinate + polygons | new (render, cheap) |

## Figure 2: Pseudospace recovers canonical PT organisation (outline: canonical segments, "higher resolution")

| Panel | Content | Source | Status |
|---|---|---|---|
| a | Canonical markers along the coordinate (for example Slc5a2, Slc5a12, Slc34a1, Slc22a6, Slc7a13, Slc22a7): species curves with specimen bin means | Notebook 37 fits | ready |
| b | Segment-label composition along the coordinate per specimen (stacked S1/S2/S3), with the transitions consistent across specimens | `segment_transitions_by_specimen.csv` | rebuild |
| c | External agreement: our gradients against microdissection and snRNA S1/S2/S3 references (gene scatter; mouse 0.74, human 0.25) | Notebook 35 | rebuild |
| d | **Within-segment order** (the "higher resolution" claim): within-segment external concordance, gene-fold stability and leave-one-specimen-out, for scFates, DPT, the landmark score and the new model, against segment labels | Notebooks 35, 36, 51, 52 | blocked (notebook 52) |
| e | Robustness: independent count split, and the injected species effect not absorbed by the coordinate | Notebook 35 | ready |

> **Wording risk.** Within-segment support for scFates is modest (external concordance 0.17 mouse,
> 0.10 human). "Higher resolution" should be stated as what panel d measures. If notebook 52 adopts a
> better ordering, this panel carries the method's improvement.

## Figure 3: Continuous analysis resolves gradients that clustering collapses (outline: beyond discrete clustering)

| Panel | Content | Source | Status |
|---|---|---|---|
| a | One gene or program shown three ways: whole-PT bar, S1/S2/S3 bars, continuous curve. One case where segment means look flat but the curve shows a replicated within-segment gradient or a sign change | Notebook 61 worked examples | rebuild |
| b | **Solute-handling micro-gradients:** heatmap of transporter genes ordered by peak, mouse and human blocks, with specimen strips | Notebook 61 heatmap code + new transporter selection | new |
| c | Within-segment gradients that replicate between specimens beyond step-only nulls, per species × segment | Notebook 58(i) | ready |
| d | Resolution budget: positional noise against within-segment spread, per coordinate, against the number of cells with replicated within-segment information | Notebook 60 | ready (rerun if notebook 52 adopts) |
| e | What continuous adds, by category: location of a difference, sign change along the PT, interior peak; counts with validation status. Also how much it adds over the S1/S2/S3 step screen: 13 more pathways, all near misses for the steps | Notebooks 57, 58, 61 | rebuild |

## Figure 4: Pathway analysis along pseudospace, methods compared

This figure exists because of the authors' interest. It could become Supplementary if space is
tight.

| Panel | Content | Source | Status |
|---|---|---|---|
| a | The specimen-relabeling control: balanced contrasts and what each relabeling tests | Notebooks 31, 37 | ready |
| b | Heat-table of methods × metrics, grouped by family (pathway-first, score-first, gene-first): calls, relabeled calls, pass/fail, programs, information returned | Notebook 63 (full grid) | blocked (running) |
| c | External many-draw null: average-level screens against positional screens (NCDR distributions) | Notebook 45 | ready |
| d | Convergence: positional methods collapse to the same ~25 programs; overlap (UpSet) by family | Notebook 63 | ready |
| e | Pathway-first then genes: one program heatmap (driver genes, member-pathway matrix), next to the best gene-first module | Notebook 61 | rebuild |

## Figure 5: Comparative analysis across species (outline: across species or samples)

| Panel | Content | Source | Status |
|---|---|---|---|
| a | Human and mouse on the shared coordinate: segment transitions per species, and canonical markers aligned | Notebooks 13, 37 | rebuild |
| b | Agreement rises with gradient strength: conservation index by strength bin, with per-bin same-species ceilings | Notebook 49 | ready (agent figure; restyle) |
| c | Gene classes: per-gene S1→S2 human against mouse, coloured conserved / reversal / species-only | Notebook 43 | ready |
| d | Where the species differ: the 12 primary programs as a peak-position (segment units) by direction summary, with one program heatmap | Notebooks 50, 61 | rebuild |
| e | Lead genes: Dcxr reversal, Gatm amplitude, Acaa2/Acadm; curves with specimen means and external references | Notebooks 38, 43 | ready |

## Figure 6: Segment-dependent transport loss and repair activation after AKI (outline: AKI)

**Blocked (on hold, by the authors' decision).** Needs notebook 02 rerun on the current segmentations. It stops at the
coarse-label checkpoint if the clusters drift from the reviewed labels.

| Panel | Content | Source | Status |
|---|---|---|---|
| a | Control against AKI on one coordinate: tubule density along the PT. Deep S3 is depleted in AKI sections (notebook 45: 248 and 91 against about 830) | Notebook 02 rerun; notebook 45 | blocked |
| b | Transport genes (heatmap, control against AKI) along the PT: where transport is lost | New | blocked |
| c | Cell-cycle and repair genes (for example Mki67, Top2a, Havcr1, Lcn2) along the PT | New | blocked |
| d | **Decoupling:** transport loss against cell-cycle activation by position | New | blocked |
| e | Specificity: control/AKI relabelings. The two AKI mice differ in severity (positional coherence 0.42; notebook 45) | New | blocked |

---

## Supplementary figures (candidates)

- **S1** Segmentation model: architecture, training data, per-class metrics, D4 test-time augmentation.
- **S2** Tubule QC: marker coherence, filters, bins and UMI per tubule.
- **S3** Integration diagnostics, the reviewed cluster labels and the fingerprint.
- **S4** Coordinate alternatives: DPT, the equal-depth refit, the landmark score and the new model; the coordinate decision.
- **S5** Label validation: glomerulus distance and reference transfer (notebook 42).
- **S6** Full pathway-method grid (notebook 63) and the program galleries (12 primary programs).
- **S7** Conservation-index sensitivity (notebooks 46, 47, 49).
- **S8** Physiological and tissue state (notebook 44).
- **S9** External datasets and lead-gene evidence (notebook 38).

## Gaps, in priority order

1. **Notebook 52 outcome** (new ordering model). Decides Figure 2d and whether Figures 3d and 5 are
   rebuilt on a new coordinate.
2. **Solute-transporter micro-gradient analysis** for Figure 3b. The abstract names "solute
   handling". Select transporters by a pre-registered rule and test within-segment replication.
3. **Renders:** the H&E overlay, the bins-in-polygon zoom and the spatial pseudospace map (Figure 1b, d, f).
4. **AKI section** (Figure 6). On hold.
5. **Title wording.** "Comparative 3D analysis" may draw objections. The method orders 2D cross-sections
   of a 3D structure on a 1D axis; it does not reconstruct 3D. Consider "along the nephron axis".
