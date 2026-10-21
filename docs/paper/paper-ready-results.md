# What is paper-ready, and what belongs in the reconstruction package

Written 2026-10-20 in answer to the user's question: "summarize the current best results we have that should go in
the paper and what the best analysis to package into the pseudospace package should be".

This page is a selection, not a catalogue. The catalogue is `docs/results/pt-discovery-methods.md` (every method
tried, with its verdict) and `docs/results/pt-package-pseudospace.md` (every notebook, with its document). This page
names the subset that is strong enough to publish, ranks it, and recommends what to ship as code.

## 0 · The framing problem to settle first

`docs/paper/outline.md` is the human-vs-mouse paper (draft v0.6, Results R1–R6, notebooks 13–50). Notebooks 70–114
rebuilt everything on the frozen reconstruction package and added an AKI line that the outline does not mention.
The outline's own re-centring note (2026-10-05) says the paper is about **the reconstruction method and what
continuous analysis shows that S1/S2/S3 does not**, with the biology as application. Nothing has rewritten the
Results table since.

**Recommendation: keep the methods-first frame and make AKI the lead application, with human vs mouse second.**

The reason is arbitration, not novelty. The AKI line has three independent external arbiters and the species line
has one weak one:

| | external arbiters | what they can judge |
|---|---|---|
| AKI along the PT | Kirita 2020 snRNA, 18 mice, 5 times, segment-resolved (94, 114); KPMP human AKI biopsies (98); Visium IRI sections from two labs (105) | level, shape, position, and **time** |
| Human vs mouse | cell atlases (Lake/KPMP, Census), Cheval 2012 microdissection | gradients only; the probe panels differ, so **levels can never be read** |

A methods paper needs the application where the method's claims can be checked by someone else's data. That is AKI.
Human vs mouse keeps its own section because the conservation-by-strength result is the strongest single biological
claim in the repository, but it cannot carry the validation argument.

This is a recommendation, not a decision. It is recorded here rather than put as a question, per the standing rule.

## 1 · Methods results that should go in the paper

Ranked by how much of the paper's argument each carries.

### M1 · A continuous axis reconstructed from unordered cross-sections, and what it costs
The frozen package (bins → units → QC → two-pass integration → elastic principal curve) with its validation
(notebook 72) and its honest limits: the curve is **noise-dominated within segments** (10 of 12 cells; human
gene-fold agreement 0.49), and stage 03 is chaotic in the thread count. Both must be stated. The within-segment
improvement from the zone-landmark ratio (notebook 53; human S1 gene-fold 0.96 against 0.07) is the benchmark panel
the user already assigned to Figure 2.
*Why it leads:* it is the contribution. Everything below is a consequence.

### M2 · Specificity by balanced specimen relabeling
With two specimens per group, average-level screens cannot be shown specific; position-dependent screens can, when
specimens within a group share positional biology (31, 37, 45). Per-call bounds make this concrete: notebook 91's
level screen calls 3,159 genes against 40 and 3 under the two relabelings, and its location screen 362 against 0
and 1.
*Why it matters:* it is the generalisable statistical result, it applies to any small-n spatial design, and it is
the reason the rest of the paper's calls are defensible. With the AKI caveat (95): the A2/A4 relabeling correlates
with the injury response itself (ρ 0.23–0.66), so for AKI external replication outranks internal calibration.

### M3 · Where continuous beats discrete, and where it does not
- Gene shape classes: continuous favoured for **28%** of non-flat genes against **3%** for segment steps (90 A).
- Positional information: **26** single pathways localise a unit better than its own segment label (90 C).
- Change density concentrates **beside, not at**, the label transitions (90 E) — the discrete model's assumption of
  within-segment constancy is what fails.
- Within-segment pathway positions (84) and transition sharpness (85) are only askable on a continuous axis.
- **The honest limit, which must be in the paper:** at set level the segment F-test matches the positional F-test
  (74, 75); the spline adds 12 sets and the steps add 5. Continuous wins at gene and position resolution, not at
  pathway-count resolution.
*Why it matters:* this is the paper's second title claim, and reporting the limit is what makes the rest credible.

### M4 · A planted-truth benchmark of the calling methods
Power and false calls with a known answer (100, 102, 103): the positional F-test finds shape, S3-only and S2-bump
changes at 0.5 log2 (≥ 0.96) and is blind to pure level (≤ 0.06); clustering finds level at 0.5 log2 and shape only
at 1 log2. One planted program makes 100–650 sets significant, so **raw set counts are never reportable**; the
network grouping returns 8–10 programs.
*Why it matters:* it is the only yardstick with ground truth, and it justifies running two screens instead of one.

### M5 · An external time course used as a ruler, and the position × time grid
Notebook 99 dates a single-time-point spatial dataset by projecting phase programs defined in independent data:
our IR mice sit at **12 h to 2 d**. Notebook 114 extends this to a position × time grid and shows the geography is
acute (strongest at 4 h, gone by 14 d) while the magnitude persists (median signed location ratio 0.16 against 0.54
for level), with the cell cycle the only late-peaking geography (2 d).
*Why it matters:* it is a reusable method for the common case of one spatial time point, and it sets the boundary
condition for every positional AKI claim in the paper (nothing beyond 2 days).

### M6 · Convergence across independent screens: count, do not weight
Notebook 112 joins every per-gene result into one atlas and nulls convergence inside expression deciles:
enrichment over the null rises 1.10× (3 lines) → 2.79× (5) → 5.08× (6) on AKI, collapses under relabeling
(795 → 150 at k ≥ 5), and is validated by annotation the lines never saw (anchored share 0.21 background → 0.62 at
k ≥ 5 → 0.78 at k ≥ 6). Notebook 113 shows a fitted weighting does **not** transfer across annotation sources (the
plain count wins on 3 of 4 testable folds), and ranks the screens on one scale: injury phase 0.79, location screen
0.76, external replication 0.40, **positional F-test 0.05 with an interval crossing zero** although it makes the
most calls.
*Why it matters:* "the method that finds the most genes says the least about which of them are real" is a result
about screens, not about our kidney. It also supplies the paper's gene table (the k ≥ 5 core).

### M7 · Transferring the axis to data without segmentation
Notebook 105 places 55 µm Visium spots from two other labs on our axis using injury-stable positional genes, with
no segmentation and no tissue x/y. This is the mechanism by which anyone else uses the coordinate.

## 2 · Biology that should go in the paper

### B1 · Acute PT injury is S3-weighted, measured continuously, and it replicates three times
- On our axis, in both directions and both mice, beyond the controls' own spread (99).
- In two independent labs' spatial IRI sections (105: Cheung 12 h ρ 0.41–0.66; Gharaie 24 h 0.45–0.46; controls
  −0.22 to 0.17).
- In 18 segment-resolved mice at 4 h (114: injury markers z +2.0, the 93 transporters −3.5, our own convergence
  core −7.2 with no external positional input).

This is the headline application result and the only claim with spatial, segment-resolved and human external
support. **Caveat to carry:** it is acute. The geography is gone by 14 days.

### B2 · Injury erases the S3 end of the gradient, most at the S2→S3 stretch
The IR PT keeps 0.52 of the healthy gradient overall and 0.31 from S2 to S3, against 0.9–1.2 under relabeling (91).
Replicated in Gharaie: S2→S3 survives at 0.18–0.20 against 0.36–0.46 for S1→S2 (105).

### B3 · Late S1 is the least responsive stretch of the PT
A replicated minimum at 0.85 units, in every injury phase and in both directions (91, 99). It does not appear at
55 µm (105), which is itself the resolution argument. Notebook 114 adds that our own down-called genes separate
S1 from S2 more than S3 from S2 (z −10.6 at 4 h), so the paper should say "S1 is spared" alongside "S3 is worst".

### B4 · The gene-level core: a catabolic loss, not a generic injury signature
The k ≥ 5 convergence core is **699 AKI genes** (435 anchored in named sets, 264 not). The strict core is a
coherent peroxisomal and mitochondrial catabolic program falling after injury (Csad, Pank1, Fmo5, Amacr, Dao,
Crot, Agps, Scp2, Acadm, Abat, Hgd, Haao, Kyat3, Mccc1) against a small induced set (Anxa2, Anxa3, Krt8, Tgm2).
This is the gene table for the paper; the ranked unanchored list (Pipox, Cdo1, Cth, Aspa, Aldh8a1, Hps5, Bphl,
Rorc, Zfp36, Lrrk2 …) is a candidate list, explicitly post hoc.

### B5 · The S3 loss splits into a sex-independent part and a male androgen program
AR-maintained S3 genes fall 0.33 log2 beyond matched genes; the knockout effect predicts the injury change beyond
gradient and transport; the fall is male-specific in an independent dataset; and the extra loss is early (4–12 h)
and gone by 2 d (111). All our specimens are male, so this is a necessary control as well as a finding.

### B6 · Proliferation is a late, patchy, outer-stripe wave — and there is no repair niche
Cell cycle is the only late-peaking geography (2 d; 114), it is patchy across tubules rather than graded (5 of 8
bins; 107), its spatial clustering is the late-S3 wave rather than a niche, and neighbours' injury and
interstitial lineages do not predict it beyond position and depth (110). Raw Moran's I misleads here.

### B7 · Human and mouse PT agree mainly on their most strongly zonated genes
Conservation index ≈ 0.55 in the top 1% of genes, 0.1–0.3 among moderately zonated genes, near zero below; direction
agreement rises from ≈ 0.5 to 0.85–0.89 with strength; transporters over-represented among conserved genes (33 of
43 zonated transporters — a count made by hand, still open item 18). The 22 position-dependent pathways mostly
summarise mouse zonation the human sections do not share.
*Keep as its own section; do not let it carry the validation argument (see §0).* The species axis is treated in
depth in `human-mouse-paper-ready.md`, which supersedes this entry: after the package switch the species result is
**gene-level, not pathway-level**, the reversed-gradient list is the strongest claim on that axis, and the
conservation-by-strength number above is the last one still computed on the retired coordinate.

### B8 · Two hypothesis-driven results that read as translational
- **Transporter handoffs are not a family-wide principle** (Q1 fails), but 11 conserved handoffs are confirmed in
  microdissection and atlases — amino acids, organic anions, SGLT2→SGLT5 — and they sit overlapping, inside S2 (101).
- **Mendelian proximal-tubulopathy genes sit toward S1 in both species** (mean E −1.47 / −0.56, beyond matched
  nulls); 25 of 27 literature positions reproduced, including OCT2 rising toward S3 in mouse but flat in human (108).

### B9 · Failed repair is transcriptionally dominant and spatially diffuse
Its magnitude more than doubles from 12 h to 14 d (Gerhardt z 14.3 → 25.9) while its S3 weighting disappears
(+3.3 → −3.7) (114). This settles an open question in favour of Kirita's "scattered throughout" reading and tells
the paper not to claim failed repair as an S3 phenomenon.

## 3 · Negatives worth publishing

They are short, they are cheap, and they are what distinguishes a careful paper.

- No repair niche at tubule resolution between 12 h and 2 d (110).
- No sequence-evolution prior on reversed gradients: reversed and conserved gradients are equally constrained
  (AUC 0.49–0.51), and an artefact forest cannot separate them (0.54) (106).
- Regulator activity is **worse** than the member sets it summarises (relabeled 0.62 per call against 0.34); treat
  regulator names as hypotheses (95).
- No pathway-level species program beyond random genes; species differences are gene-level (90 D, 91 S).
- A fitted gene score does not transfer; the unweighted count does (113).
- Three of five annotation sources cannot be held out at all, because curated sets are nested (113). This is a
  general warning about leave-one-source-out validation on curated annotation.

## 4 · What to ship in the reconstruction package

The package stops at a coordinate: stages 01–04 build it and 05–07 compare its clustering. Everything that makes
the coordinate *useful* lives in this repository and is not shippable as written. Recommendation, in order, as one
new `positional` layer rather than renumbered stages (05–07 are taken):

### P1 · The positional test with its relabeling calibration — ship these two together
- **What:** the nested NB model on raw counts with the library offset, the moderated partial F for the shape term
  (`T_spatial`), and set-level `cameraPR`; plus the balanced specimen relabeling that bounds every call.
- **Where it is now:** `pseudospace/nb_gam.py` (the NB/IRLS engine, dispersion, HC3 contrasts),
  `pseudospace/pathway_calibration.py` (`balanced_partitions`, `contrast_designs`, `nested_partial_f`),
  `pseudospace/enrichment.py`, `pseudospace/specimen_null.py`.
- **Why first:** it is the only calibrated positional test we have, it is the one the paper recommends, and it is
  useless without the relabeling — shipping the F-test alone would hand users a screen that fails calibration on
  one of our two datasets and nobody would know. Ship them as one call that returns calls **and** their relabeled
  bound, so the bound is not optional.
- **Generality:** nothing in it is kidney-specific. It needs counts, a coordinate, a group label and a specimen
  label. It works for crypts or lobules unchanged, which is the stated requirement for the method.

### P2 · The level / location / gradient-survival decomposition
- **What:** per gene, a whole-structure level change, a location change beyond level, and the survival of the
  healthy gradient, each with a matched null and a relabeled bound (notebook 91).
- **Where it is now:** `pseudospace/positional_decomposition.py` (`level_location`, `survival_parts`,
  `cross_fitted_survival`, `set_mean_z`, `interior_dip`).
- **Why second:** it is the best description of a condition effect along an axis that we have, it produced B2, B3
  and most of B4, and 113 scores its location term as the second most informative line of all nine (0.76). It is
  also the smallest, most readable output: three numbers per gene.

### P3 · Axis transfer: place external units or spots on a fitted coordinate
- **What:** choose condition-stable positional genes, score foreign data against them, return a rank position
  (notebook 105's machinery) with the PT-purity and depth guards.
- **Where it is now:** `pseudospace/external_spatial.py` (`type_markers`, `nnls_weights`, `rank_position`).
- **Why third:** it is the adoption driver. It lets someone with 55 µm Visium and no segmentation use a coordinate
  built from segmented data, and it is how the coordinate gets externally validated at all.

### P4 · The convergence atlas (optional, ship last)
`pseudospace/gene_atlas.py` (`join_lines`, `stratified_line_null`, `convergence_counts`, `yield_profile`) is
general and small, but it is a meta-analysis layer: it is only worth anything once a user has run several screens.
Ship it as an example rather than a stage.

### Do not ship
- Any **fitted gene score** or weighted line combination (113: it does not transfer; the count does).
- **Regulon and footprint activity** as a discovery screen (95: worse than member sets). PROGENy stays as a
  signalling add-on in the analysis repository only.
- **GSEA/ssGSEA on raw positional statistics** and the constant-offset pathway screens (66, 75: rejected; not
  specific in a few-specimen design).
- **Raw set counts** of any screen (100, 102).

### Two things the package must carry as warnings, not features
1. Stage 03 is chaotic: the BLAS/OpenMP thread count or an unsorted sparse layer changes the clusters and moves
   units. Pin threads and record them, or refits do not reproduce.
2. The elastic curve is weak within segments. Any within-segment claim needs the zone-landmark ratio or an
   equivalent, and the package should say so at the point where it returns the coordinate.

## 5 · What is still missing before any of this can be submitted

Unchanged from `docs/paper/open_items.md` section A, and none of it is recoverable from the repository: mouse strain
and age, human donor provenance and consent, ethics approvals, tissue processing and sequencing details, the human
pixel scale, and segmentation provenance for the human sections. The science is further along than the records.
