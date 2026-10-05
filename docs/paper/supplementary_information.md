# Supplementary Information, draft v0.6

Supplementary Notes 1–10 are followed by the Supplementary Methods. The standing caveat applies
throughout: two control mice and two cortex sections from one male human donor; all specimens male;
human cortex only; different probe panels; descriptive.

Markers: **[VERIFY]** (lab records), **[CITE]**.

---

## Supplementary Note 1 · The continuous PT coordinate

**Construction** (Supplementary Methods SM9).
- *Fit.* A nonbranching principal curve (scFates; Faure et al. 2023; ElPiGraph, Albergante et al. 2020)
  with 30 nodes, on the first five pass-2 Harmony dimensions of the 12,866 PT structures.
- *Orientation.* The curve was rooted at the tip with the higher early-marker score, mapped once, and
  scaled to [0, 1]. The marker axis is the mean within-species z-score of late genes (Slc22a7,
  Slc7a13, Cyp7b1) minus that of early genes (Slc5a2, Slc5a12, Gatm). The axis sets only the root.
- *Dimensionality.* The curve lives in five dimensions. Its trace in a two-dimensional projection
  (Supplementary Figure S1a) can loop back on itself.
- *Transitions.* Reviewed S1→S2 transitions fall at 0.26–0.28 in mouse and 0.35–0.38 in human; S2→S3
  transitions at 0.63–0.64 in all four specimens. In every specimen, the S3 median exceeds the S1 and
  S2 medians.

**What it recovers** (Supplementary Figure S1b, c). Each gene was scored on a coordinate rebuilt
without it (three gene folds). Its within-specimen rank correlation with the coordinate was then
correlated, across genes, with external S3-versus-early contrasts.

| Species | External reference | Coordinate | Segment labels |
|---|---|---:|---:|
| Mouse | Male snRNA | 0.74 | 0.73 |
| Mouse | Microdissection | 0.53 | 0.51 |
| Human | Cortex snRNA (convoluted PT against S3 only) | 0.25 | 0.27 |

- *Fold protection.* Without it, agreement was 0.75 in mouse and 0.31 in human, so part of the human
  agreement comes from the genes that built the coordinate.
- *Injected gradient.* A human-only gradient injected into a tenth of the construction genes survived
  refitting (absorption ratio 1.01).
- *Externally species-specific genes.* For 244 genes whose zonation differs between species in
  external data, the coordinate kept the external direction for 91%, against 95% for segment labels.
- *Count splitting* (Neufeld et al. 2023). Building the coordinate on half of the reads and testing on
  the other half kept 59 of 66 robust and 39 of 42 core pathways from the list of the method-selection
  analysis. Reusing the same reads inflated gene statistics by 34%, and the circular arm called more
  pathways than the independent arm (145 against 130).

**What it does not resolve** (Supplementary Figure S1d–f).
- *Within-segment order.* Slopes along the coordinate within segments agreed only weakly with adjacent
  external contrasts: 0.17 in mouse and 0.10 in human. Ordering by read depth gave about 0, and −0.16
  in human early PT.
- *Gene-fold stability.* Coordinates built from disjoint thirds of the genes ranked structures within
  segments only moderately alike (0.34–0.70).
- *Prediction.* The coordinate predicted held-out genes in the other specimen no better than segment
  means.
- *Depth.* After thinning every structure to equal depth, within-segment agreement with the original
  order was 0.36–0.46 in human S1. The correlation of position with library size fell from −0.31 and
  −0.38 to −0.04 and −0.13.
- *Registration.* Removing a third of the genes or half the reads moved human segments by up to 0.29
  relative to mouse. In one gene-fold refit the human S1 < S2 < S3 order broke.
- *Numerical sensitivity.* Input changes of order 10⁻¹⁶ moved individual structures by up to 0.44.
  Twenty resampled mappings moved positions by SD 0.007 and changed T_spatial only slightly (ρ 0.994).

**Choice of coordinate** (Supplementary Note 8, deviation 1).
- *Candidates.* The equal-depth refit of the curve (SCF13-ED), diffusion pseudotime on the same
  embedding (DPT13), a PT-only integration, and a conserved-anchor coordinate.
- *The pre-registered rule.* It required:
  - equal-depth agreement ≥ 0.80;
  - external agreement within 0.05 of the principal curve in each species;
  - a species-registration gap ≤ 0.15.
- *What it selected.* DPT13 met all three criteria (equal-depth agreement 0.80 against 0.72 for the
  principal curve; gap 0.14 against 0.29). It also had better gene-fold agreement (0.59 against 0.49)
  and better human external agreement (0.29 against 0.25).
- *Failures.* SCF13-ED failed the rule. The conserved-anchor coordinate did not register the species.
- *The specificity gate.* The rule did not contain the pre-registered downstream gate (NCDR < 0.25 and
  < 5% of pathways called under relabeling). After DPT13 was seen to fail it on the matched-only test
  (NCDR 0.40), the gate was applied to every candidate; none qualified, and the principal curve was
  kept.
- *The joint test.* On the joint test, which is the primary pathway test, DPT13 passes: 45 species calls
  against 3 and 0 relabeled (NCDR 0.03).
- *Decision after review.* The primary pathway list is the set robust under both the principal curve
  and DPT13, and each coordinate is reported separately.

**Pathway calls under three coordinates** (Supplementary Figure S2; Table S2). The final pathway
analysis was run unchanged on each coordinate.

| | Principal curve | DPT13 | SCF13-ED |
|---|---:|---:|---:|
| Joint T_spatial: species calls (relabeled) | 60 (2, 0) | 45 (3, 0) | 66 (1, 0) |
| Joint NCDR | 0.02 | 0.03 | 0.01 |
| Matched T_spatial: species (relabeled), NCDR | 116 (17, 13), 0.13 | 99 (41, 36), 0.39 | 123 (19, 2), 0.09 |
| Joint step screen: species (relabeled) | 16 (0, 0) | 14 (0, 0) | 31 (5, 6) |
| Within-species positional calls (mouse, human) | 16, 0 | 37, 5 | 12, 3 |
| Candidates / not specimen-sensitive / stable | 60 / 58 / 51 | 45 / 42 / 33 | 66 / 65 / 60 |
| Robust / core | 36 / 26 | 26 / 15 | 51 / 36 |
| Robust human-high / mouse-high at peak | 10 / 26 | 10 / 16 | 17 / 34 |
| Principal-curve robust kept as robust | – | 22 / 36 | 32 / 36 |
| Principal-curve core kept as core | – | 12 / 26 | 21 / 26 |
| Robust called beyond steps / by the step screen | 35 / 9 | 14 / 10 | 47 / 23 |
| Median peak, human-high / mouse-high (segment units; each coordinate's own robust list: 36, 26, 51) | 0.50 / 1.99 | 0.95 / 1.90 | 0.54 / 1.88 |

- *Within-mouse calls, 37 against 36.* The DPT13 within-mouse count is 37 in the full rerun and 36 in
  the separate matched-test diagnostic of the coordinate comparison. These are different runs of the
  same coordinate.
- *Primary list.* 22 pathways, of which 13 are summaries of mouse zonation, 3 shared, 0 human-led,
  1 mixed and 5 unresolved; 12 are core under both coordinates.
- *Direction.* At the peak it agrees under all three coordinates for all 40 pathways robust under the
  principal curve or DPT13.
- *Peak positions.* They correlate at Spearman 0.71 between the principal curve and DPT13 (22
  pathways; 7 change segment) and at 0.54 between the principal curve and SCF13-ED (32 pathways; 10
  change segment, 6 move by more than 0.25). Three pathways move from about 0.6–0.75 to near 0 under
  SCF13-ED: Vitamin digestion and absorption, Metabolism Of Fat-Soluble Vitamins, and Nitrogen
  metabolism. Peaks are therefore reported in segment units only.

**Unequal positional noise does not create the early/late split.** On the 36 principal-curve robust
pathways, human-high pathways peak early and mouse-high pathways late. Adding the estimated excess
human positional noise to mouse positions (σ 0.089) accounts for about 15% of the gap between the
median peaks and 7% of the mouse-high count (σ 0.10: 14% and 8%).

**What the coordinate adds.** Not finer-than-segment resolution:
- it replicates measurements across positions within each specimen;
- over the step screen it adds power, not specificity (60 against 16 joint-test calls, both with 0–2
  under relabeling).

---

## Supplementary Note 2 · Specificity details

**All strategies** are in Table S1. Figure 4b shows five; Supplementary Figure S3a shows all.

**Target–decoy estimates** (Elias & Gygi 2007; Supplementary Figure S3b).
- *Average-level screens.* No competitive average-level screen reached a decoy-estimated FDP of 10% at
  any threshold. The lowest attainable FDPs were 1.00, 0.75 and 0.95.
- *Pathway scores.* The self-contained pathway-score offset model reached it with 1,087 calls, because
  nearly every pathway score differs more between species than between sections.

**The pathway layer.**
- *Genes.* Whole-PT DESeq2 (Love et al. 2014) called 9,010 genes for species and only 2 and 8 for the
  relabelings. Specimen-level variation moves whole pathways together, which gene-permutation tests
  assume away (Goeman & Bühlmann 2007; Wu & Smyth 2012).
- *Decoy inflation* (Supplementary Figure S3d). Mean decoy z² rose with pathway size more steeply for
  the offset statistic (inflation 1.40; implied correlation 0.019) than for the position statistic
  (1.03; 0.009).

**Probe panels.** Per-gene probe efficiency differs between the human and mouse panels. It cancels in
within-species gradients but not in species offsets, so average-level species differences are hard to
interpret in any case.

**Within-species controls** (Supplementary Figure S3c, titled "Within-species controls"). The two mice
gave 48 offset and 16 positional pathway calls, and the two human sections gave 94 and 0. The human
sections come from one donor, so their 0 positional calls say nothing about donor-to-donor positional
variation.

**Matched and joint tests.**
- *Matched test.* Under relabeling it was anti-conservative alone: on average 59.5 pathways reached
  p ≤ 0.01, where about 15 are expected.
- *Joint test.* It divides the matched z by a variance-inflation factor estimated from model
  residuals. Its residual correlation (median ρ ≈ 0.010) matched the inflation estimated from decoys
  (ρ ≈ 0.009).

**Specimen-level split-plot test** (Supplementary Figure S3e). This test treats specimen × position
pseudobulks as replicates.
- *Genes.* It was calibrated under relabeling (0.6% and 0.4% of genes at p ≤ 0.01; λ ≤ 0.87). Its gene
  ranking agreed only weakly with T_spatial (Spearman 0.21).
- *Pre-specified pathway test.* It called nothing. Decoys showed that its pseudobulk correlation
  correction overcorrects (decoy inflation 1.10, against a median correction of 2.07).
- *Matched test alone.* It called 41 pathways, with 0 and 0 under relabeling. Of the 22 primary
  pathways, 7 were among them; of the 36 principal-curve robust pathways, 15.

**Many-draw nulls: our cohort as an anchor.** For the external comparison, notebook 45 re-implements
the strategies at unit level for our cohort.
- *Whole-PT GSEA.* The re-implementation's relabeled calls are 89 and 266 (Table S1: 193 and 255).
- *Two step-screen designs, 50 against 16 species calls.* Both are joint tests of species ×
  S1/S2/S3 steps beyond a species offset, on identical structures, genes and pathways.
  - The coordinate-free design (notebook 45; 50 calls; 0 and 0 relabeled) uses segment intercepts as
    the shared shape, with exact matched moments and a segment-centred variance inflation. The external
    nulls calibrate this design.
  - The spline-based design (notebook 37; 16 calls; 0 and 0 relabeled) adds the steps to a shared
    6-df spline along the coordinate, so it is not coordinate-free.
  - 15 of the 16 are among the 50.
  - The 50 contain 28 of the 36 principal-curve robust pathways (against 9 for the 16), 20 of the 22
    primary pathways (against 8) and 47 of the 60 T_spatial candidates.
- *Percentiles.* Those of our cohort in Figure 4c use the re-implementation.

**AKI positive control under both rules** (Figure 4e; 1,429 tested pathways).
- *R4 rule:* NCDR < 0.25 and mean relabeled calls < 5% of tested pathways.
- *Positive-control rule:* pre-specified; relabeled calls ≤ 5% of condition calls.

| Strategy | AKI calls | Relabeled calls | NCDR | Mean relabeled ÷ tested | R4 rule | Positive-control rule |
|---|---:|---|---:|---:|:-:|:-:|
| Whole-PT pseudobulk, GSEA | 649 | 208, 20 | 0.18 | 8.0% | fail | fail |
| S1/S2/S3 pseudobulk, GSEA | 790 | 292, 207 | 0.32 | 17.5% | fail | fail |
| Offset, joint test | 56 | 10, 0 | 0.09 | 0.3% | pass | fail |
| Step screen (coordinate-free), joint test | 26 | 12, 6 | 0.35 | 0.6% | fail | fail |
| Offset, matched test | 236 | 152, 11 | 0.35 | 5.7% | fail | fail |
| Step screen, matched test | 169 | 131, 120 | 0.74 | 8.8% | fail | fail |
| Whole-PT pseudobulk, joint rank test | 394 | 70, 0 | 0.09 | 2.4% | pass | fail |
| S1/S2/S3 pseudobulk, joint rank test | 447 | 81, 11 | 0.10 | 3.2% | pass | fail |
| Coordinate T_spatial, joint test | 50 | 10, 7 | 0.17 | 0.6% | pass | fail |

The AKI kidneys' within-group specimen × segment coherence is 0.42, against 0.013 in our cohort and at
least 0.14 in the external human pools.

---

## Supplementary Note 3 · Alternative zonation-class rules

The primary setting scales the human effect floor and flat margin to the lower human amplitude, and it
requires species-only calls to be flat-compatible in the other species' reference atlas. The table
gives class counts (probe-balanced genes; union of S2 − S1 and S3c − early) under every setting.

| Setting | Conserved | Reversal | Mouse-only | Human-only | Indeterminate | Neither |
|---|---:|---:|---:|---:|---:|---:|
| Unscaled thresholds, no flat confirmation (v1) | 155 | 52 | 247 | 66 | 5,590 | 1,297 |
| v1 + female mouse snRNA confirmation | 155 | 48 | 222 | 67 | 5,605 | 1,310 |
| v1 + female microdissection confirmation | 137 | 34 | 167 | 69 | 5,657 | 1,343 |
| Scaled thresholds only | 196 | 75 | 59 | 141 | 6,224 | 712 |
| Flat confirmation only | 155 | 52 | 187 | 31 | 5,685 | 1,297 |
| **Symmetric (scaled + flat confirmation; primary)** | **196** | **75** | **41** | **76** | **6,307** | **712** |
| Symmetric + female mouse snRNA | 187 | 62 | 38 | 71 | 6,334 | 715 |
| Symmetric + female microdissection | 161 | 47 | 31 | 78 | 6,371 | 719 |

**Per contrast** (symmetric setting):

| Contrast | Conserved | Reversal | Mouse-only | Human-only |
|---|---:|---:|---:|---:|
| S2 − S1 | 109 | 37 | 31 | 58 |
| S3c − early | 125 | 41 | 13 | 28 |
| S3 − early (mostly outer-stripe mouse S3) | 197 | 92 | 1 | 110 |

- *Strong mouse-only genes.* Under the symmetric rules, only 1 mouse-only gene has a strong gradient
  (|mouse| ≥ 1.5 log2 with human within ±0.3).
- *Enrichments.* Notebook 40's enrichment code was rerun on the symmetric classes. It reproduced all
  5,680 rows of its saved table when given notebook 40's classes.
  - *Naming rule.* An enrichment is named only if q ≤ 0.10, it is not expression-driven, and at least
    three class genes carry it. The rule was declared after the all-gene run had been seen, and before
    the probe-balanced run (Supplementary Note 8).
  - *Probe-balanced classes.*
    - Mouse-only: KEGG tyrosine metabolism (q = 2 × 10⁻⁶), carried by three genes: Comt and Maoa
      (catecholamine degradation) and Fah (tyrosine catabolism).
    - Human-only: Cardiac Conduction (Ahcyl1, Hipk2, Wwtr1) and Regulation Of TP53 Activity Thru
      Phosphorylation (Hipk2, Prkaa2, Taf4). These are generic labels on overlapping three-gene sets,
      so the main text names the genes.
    - Conserved: 12 counting pathways, among them amino-acid and SLC-mediated transport, transport of
      bile salts and organic acids, bile secretion, apical surface and renin secretion (Ace, Agt,
      Gnai1, Ppp3cc).
    - Reversal: none counting.
  - *Not named.* Sterol synthesis rests on Ebp alone (q 0.08), and peroxisomal β-oxidation on Hsd17b4
    alone. The v1 enrichments (sterol synthesis, peroxisomal β-oxidation, pyrimidine catabolism) do not
    hold under the symmetric classes and are not used.
- *Sex bias.* Class membership was cross-tabulated against mouse PT sex bias (Xiong et al. 2023) after
  the classes were computed.

---

## Supplementary Note 4 · Amplitude, physiological state and tissue state

**Amplitude** (Supplementary Figure S5). Ratios from regressing each species' gradients on a reference
gradient depend on the reference axis:
- 0.15 (0.11–0.19) on male mouse microdissection, the pre-specified primary;
- 0.15–0.23 on four mouse references;
- 1.96 (1.41–3.05) on the healthy human cortex atlas;
- external atlases alone: 0.14 on a mouse axis and 1.66 on a human axis.

Deming and orthogonal slopes on the human axis were 0.16–0.17. With a weak cross-species correlation,
total-least-squares slopes track variance ratios rather than covariance. The reference-free
noise-corrected spread (0.58, 0.47 and 0.84; main text R2) is therefore the amplitude measure.

**Measurement checks.**
- The top expression tertile gave a ratio of 0.18, and genes with equal probe counts 0.16, so log
  compression does not explain the difference.
- Low- and high-injury halves of our structures gave 0.18 and 0.19.

**Tissue state** (Supplementary Figure S6c).
- *Score.* A PT injury signature was defined from the cortex atlas: 30 up and 30 down genes of
  altered-state against healthy PT. The score is the mean within-sample percentile of the up genes
  minus that of the down genes.
- *Calibration.* A platform offset was estimated in mouse, from our control mice against healthy
  Census donors. As a check, our AKI kidneys scored 0.33 above our controls. The offset is assumed to
  transfer to human.
- *Our donor.* Both sections score −0.45 after calibration, at the median of the 14 healthy reference
  donors, with a predicted altered-state fraction of 0.50.
- *Injury and amplitude among reference donors.* Median human-axis amplitude was 0.69 in healthy,
  0.50 in chronic kidney disease and 0.35 in AKI donors. The Spearman correlation of injury with
  amplitude was −0.41 on the human axis (post hoc) and −0.24 on the pre-specified mouse axis.
- *Mouse AKI.* It kept 0.74 of control S2 − S1 amplitude (gradient correlation 0.83) and cut S3 − early
  amplitude to 0.41–0.44 (correlation 0.59).

**Physiological state** (Supplementary Figure S6a, b).
- *Programs* (fixed in advance): sterol/SREBP2, PPARα/fatty-acid oxidation, glutathione synthesis,
  tyrosine catabolism, the lead genes and structural markers.
- *Fasting.* Whole-kidney 24-h fasting against fed (GSE267280; Jiang et al. 2025) changes PPARα-target and sterol
  genes by amounts comparable to the species × position interaction. There is one pooled column per
  condition, so this is a point estimate.
- *Time of day, whole kidney.* ZT18 against ZT6 (GSE277302; Nguyen et al. 2025).
- *Time of day, segment level.* Mouse single-nucleus data at ZT4 and ZT16, two mice each (GSE331332;
  Wigger et al. 2026).
  - Zonation (S3 − S1) of fatty-acid oxidation, glutathione synthesis and tyrosine catabolism genes,
    the lead genes and the structural markers kept its sign in all four mice for every zonated gene,
    with amplitude ratios of 0.99–1.35.
  - Sterol-synthesis zonation was not stable: 73% of zonated genes were sign-stable, and the change was
    at least half the species gap.
- *A state marker.* PDK4, a fasting- and stress-responsive gene, is detected in 94% of human PT
  structures against 7% of mouse. The ketogenic gene HMGCS2 is low in human. These level differences
  suggest a difference in physiological state but do not establish one.
- *Literature.* Mouse PT PPARα is fasting-activated and drives renal fatty-acid oxidation (Aomura et
  al. 2026). Early- and late-PT mitochondrial fatty-acid oxidation capacities differ and change with
  fasting (Feola et al. 2026). The renal tubular clock controls β-oxidation (Bignon et al. 2023).

---

## Supplementary Note 5 · External agreement

**Gene level.** Spearman correlation of our species × segment interaction with the external
interaction (human cortex atlas − male mouse snRNA):

| Reference | Our top-10% genes by T_spatial | Expression-matched random, median (97.5th percentile) | Highest-expressed genes |
|---|---:|---|---:|
| Census human cortex | 0.67 | 0.57 (0.62) | 0.64 |
| Lake/KPMP cortex | 0.73 | 0.60 (0.66) | 0.72 |

**Pathway level** (BH ≤ 0.10 within the robust list):

| Reference | 36 principal-curve robust | 22 primary |
|---|---:|---:|
| Census cortex − male mouse | 17 | 12 |
| Census cortex − female mouse | 8 | 5 |
| Lake/KPMP cortex − male mouse | 6 | 5 |
| Lake/KPMP cortex − female mouse | 4 | 3 |

- *Replication with overlap kept* (notebooks 48 and 50). The score is each pathway's donor-level joint
  z. Each null draw relabels genes within the matching strata and rescores every set at once, so set
  sizes, strata and pairwise overlaps are kept (9,999 draws per reference). Male-mouse references:

  | | Census cortex | Lake/KPMP cortex |
  |---|---:|---:|
  | 22 primary, mean z (p against random) | 1.62 (1 × 10⁻⁴) | 1.44 (1 × 10⁻⁴) |
  | 12 primary core, mean z | 1.78 (1 × 10⁻⁴) | 1.42 (1 × 10⁻⁴) |
  | 21 relabel-called, not robust, mean z (p against random) | 0.59 (0.026) | 0.44 (0.070) |
  | 22 − relabel-called, pathways (p) | 1.03 (0.004) | 0.99 (0.007) |
  | Programs, 13 against 13 (label permutation p) | 0.039 | 0.065 |
  | Selection-matched null: 22, mean z (p against random) | 0.82 (0.001) | 0.69 (0.006) |
  | Selection-matched null: programs against relabel-called (label permutation p) | 0.27 | 0.37 |
  | Slope-free score, 22 − relabel-called (p) | 0.26 | 0.13 |

  - *Reading.* The 22 replicate beyond random sets in both atlases. They replicate beyond relabel-called
    pathways at the pathway level in both, but at the program level, the pre-specified unit, in one of
    two atlases. They do not once the shared flattening of human gradients is removed.
  - *Selection.* The 22 were defined after the coordinate results.
  - *The 36.* The same test on the 36 principal-curve pathways gave mean z 1.53 and 1.24 (p = 10⁻⁴),
    differences from relabel-called pathways with p = 0.005 and 0.015, and program-level p = 0.059 and
    0.17: better than random only.
  - *Earlier comparisons.* This test replaces the Mann–Whitney comparisons, which treated overlapping,
    selected pathways as independent.
- *Global component.* All three human datasets share it: the interaction is −0.81 to −0.90 times the
  mouse gradient (R² 0.60–0.82). Once it is removed, 1–3 of the 36 robust pathways replicate.
- *Conventional-only pathways.* Of 184 pathways called only by conventional screens, 6 replicated their
  average difference at nominal p ≤ 0.05 and none at BH ≤ 0.10. An overlap-naive comparison suggested
  that they replicate better as a group, but against overlap-preserving random collections they did not
  (p = 0.41 and 0.24). No group claim is made.

**Genome-wide against the cortex atlas** (Supplementary Figure S7).
- *Human S3 − S1.* Our human contrast correlated with the atlas at Spearman 0.43 (6,787 genes), with
  93.1% direction agreement for 335 strongly zonated genes.
- *Interaction.* For the 10% most position-dependent genes, our interaction correlated with atlas −
  male mouse snRNA at 0.70, with 96.7% direction agreement for 276 genes strong in both. With female
  mice the figures were 0.67 and 93.7%.

---

## Supplementary Note 6 · Additional human Visium donors and rat segment protein

**Visium donors** (Abedini et al. 2024; GSE211785).
- *Samples.* 10 samples with ≥ 50 PT_S1 and ≥ 50 PT_S3 spots: 2 controls and 8 with disease.
- *Label check, pre-specified.* SLC5A2 falls, and SLC5A1, SLC7A13 and SLC22A7 rise, in ≥ 70% of samples.
  It failed: the S3 markers did not rise in the authors' PT_S3 spots. Spot-level segment labels in
  these samples therefore do not separate S1 from S3 as the markers require.
- *Caveat, for every row below.* The labels failed the check, so no row is evidence for or against a
  positional claim. The table is reported for completeness and is not cited in the main text.
- *Columns.*
  - "Pre-specified": samples with |S3 − S1| ≥ 0.1, and the share agreeing with the cortex-atlas sign.
  - "Post hoc": the same after centring each sample on its own median gene (decided after the label
    check failed).

| Gene | Cortex atlas S3 − S1 | Visium mean S3 − S1 (SE) | Informative samples | Agree with atlas sign (pre-specified) | Centred mean (post hoc) | Agree with atlas sign (post hoc) | Agree with mouse sign (post hoc) |
|---|---:|---:|---:|---:|---:|---:|---:|
| GATM | +0.51 | −0.95 (0.09) | 10 | 0 of 10 | −1.02 | 0 of 10 | 10 of 10 |
| GAMT | +0.09 | −0.39 (0.13) | 7 | 0 of 7 | −0.46 | 0 of 8 | 0 of 8 |
| ACADM | +0.46 | −0.22 (0.10) | 10 | 2 of 10 | −0.29 | 2 of 10 | 2 of 10 |
| ACAA2 | +0.58 | −0.26 (0.09) | 7 | 0 of 7 | −0.33 | 0 of 7 | 7 of 7 |
| DCXR | +1.53 | −0.30 (0.15) | 7 | 2 of 7 | −0.36 | 2 of 8 | 6 of 8 |
| UGT3A1 | −0.98 | −0.67 (0.19) | 8 | 8 of 8 | −0.74 | 8 of 8 | 0 of 8 |
| ACOX2 | −0.46 | −0.48 (0.11) | 8 | 8 of 8 | −0.55 | 8 of 8 | 0 of 8 |
| GCLC | +0.61 | −0.09 (0.25) | 7 | 2 of 7 | −0.16 | 1 of 6 | 1 of 6 |
| GCLM | −0.76 | −0.04 (0.06) | 5 | 3 of 5 | −0.11 | 5 of 7 | 2 of 7 |
| GSS | +0.67 | −0.04 (0.11) | 8 | 5 of 8 | −0.11 | 3 of 6 | 3 of 6 |
| UGT1A9 | −0.67 | −1.03 (0.22) | 9 | 9 of 9 | −1.10 | 9 of 9 | 0 of 9 |
| RBP4 | +3.46 | +1.07 (0.28) | 7 | 7 of 7 | +1.00 | 7 of 8 | 1 of 8 |
| AOX1 | +1.13 | −0.30 (0.10) | 8 | 1 of 8 | −0.37 | 1 of 8 | 1 of 8 |
| SLC6A19 | −0.30 | −1.07 (0.30) | 10 | 10 of 10 | −1.14 | 10 of 10 | 10 of 10 |
| SLC22A6 | +0.32 | −1.32 (0.28) | 10 | 0 of 10 | −1.39 | 0 of 10 | 10 of 10 |
| SLC13A3 | −0.23 | −0.90 (0.15) | 10 | 10 of 10 | −0.97 | 10 of 10 | 10 of 10 |

Most genes fall from S1 to S3 in these samples, including genes that rise in the cortex atlas (GATM,
ACADM, ACAA2, DCXR, AOX1, SLC22A6). That pattern is what the failed label check predicts.

**Rat segment protein** (Limbutara et al. 2020).
- *Data.* Microdissected S1, S2 and S3, 3 samples per segment, in copies per cell. Log2 ratios are
  centred on the median ratio of proteins quantified in all nine PT samples.
- *Why S3/S2.* Rat S1 samples carry less protein, so S3/S2 is the primary ratio. Ratios involving a
  segment where the protein was not detected are unbounded and are given for completeness only.

| Gene | Detected in S1 / S2 / S3 (of 3) | log2 S3/S2 | log2 S3/S1 |
|---|---|---:|---:|
| Gatm | 3/3/3 | −3.7 | −2.9 |
| Gamt | 3/3/3 | +1.6 | +6.6 |
| Gclc | 0/3/3 | +2.8 | (not detected in S1) |
| Gclm | 3/3/3 | +1.9 | +10.0 |
| Gss | 3/3/3 | +2.4 | +4.0 |
| Dcxr | 3/3/3 | −1.1 | −1.8 |
| Acadm | 3/3/3 | +0.7 | +2.4 |
| Acaa2 | 3/3/3 | −2.7 | −2.8 |
| Slc22a6 (OAT1) | 3/3/1 | −6.5 | −3.1 |
| Slc22a7 (OAT2) | 0/2/3 | +1.0 | (not detected in S1) |
| Slc22a8 (OAT3) | 3/3/0 | (not detected in S3) | (not detected in S3) |
| Cyp2e1 | 0/3/1 | −5.5 | (not detected in S1) |
| Rbp4 | 3/0/0 | (S1 only) | (S1 only) |
| Slc6a19 | 0/3/0 | (S2 only) | (S2 only) |
| Slc5a2 (SGLT2) | 3/3/2 | −3.8 | −9.4 |
| Slc5a1 (SGLT1) | 0/3/3 | +0.7 | (not detected in S1) |
| Aqp1 | 3/3/3 | +3.7 | +5.2 |
| Ugt3a1 | not quantified | | |

---

## Supplementary Note 7 · Segmentation model

The model, training, selection and inference are described in Supplementary Methods SM3. In summary:
- *Encoder and decoder.* A frozen pathology foundation encoder (a DINOv2 ViT-g/14 with registers;
  Oquab et al. 2023; Kaplan et al. 2025) feeds a trainable decoder at output stride 1.
- *Decoding.* Instances come from a marker-controlled watershed.
- *Training.* Partially annotated mouse slides were used, with an ignore mask for unannotated tissue.
  Part of `Ctrl1A2` was in the training data.
- *Selection.* The checkpoint was chosen by pooled class-agnostic panoptic quality, never by validation
  loss.
- *Held-out performance.* On a held-out nephrectomy slide (98 objects), panoptic quality was 0.837
  without and 0.880 with D4 test-time augmentation. This is one slide, from mouse, and the human
  segmentation has no comparable evaluation.
- *Not recorded.* The human segmentation's provenance and pixel scale [VERIFY], and the rules of the
  upstream mouse quality control [VERIFY].

---

## Supplementary Note 8 · Deviations from pre-specified protocols

| # | Analysis | Pre-specified | What was done | Why | Effect on the results |
|---|---|---|---|---|---|
| 1 | Coordinate selection | A three-criterion rule; the principal curve kept if no candidate passed | The rule selected DPT13. A downstream specificity gate, pre-registered elsewhere but missing from the rule, was applied after DPT13 failed it on the matched-only test (NCDR 0.40); the principal curve was kept. On the joint test DPT13 passes (NCDR 0.03). After review, the primary pathway list became the 22 pathways robust under both coordinates | The gate was omitted from the rule; the override relied on the matched-only test | Pathway lists differ by coordinate (36 against 26; 22 shared). Directions agree for all 40. Zonation results do not use the coordinate |
| 2 | Specificity rule | NCDR < 0.25 | The condition "each relabeling calls < 5% of tested pathways" was added after the first run | A pathway-score model passed the NCDR condition alone | Average-level screens fail either way |
| 3 | Primary specificity metric | Decoy-estimated FDP | Replaced by NCDR | Pathway-score models reach decoy FDP 0 because nearly every score differs between species | Decoy FDP is reported (Supplementary Note 2) |
| 4 | Pathway test | Matched rank test | The joint test replaced it after the matched test's relabeling behaviour was seen | Matched test anti-conservative under relabeling (59.5 against 15 at p ≤ 0.01) | Both reported (Table S1) |
| 5 | Zonation classes | Unscaled thresholds, zonated calls confirmed externally (v1) | Symmetric classes: scaled thresholds and external flat confirmation, pre-registered before computing them | Reviewer: flat calls were unconfirmed and thresholds asymmetric | The mouse-only excess depends on the scale (main text R3); all settings in Supplementary Note 3 |
| 6 | Cross-fitted labels | Species-only counts "unskewed" (ratio < 2 in every setting, half and contrast) | Criterion failed in 1 of 12 checks (17 against 4; the other half 7 against 22) | – | Reported as failed (main text R3) |
| 7 | Amplitude | Projection ratios on mouse references | Noise-corrected spread and disattenuated correlations added post hoc | Projection ratios depend on the axis when correlation is low | The spread is the amplitude measure |
| 8 | Sex-bias check of classes | Not planned | Added after the classes were computed | Reviewer request | Descriptive |
| 9 | Many-draw nulls | Not in the original plan | Added after review, with a decision rule fixed before they ran. The age analysis of null calls was post hoc | Reviewer: two null draws | Main text R4 |
| 10 | AKI positive control | Specificity if relabeled ≤ 5% of condition calls | The main R4 rule is also reported; under it the coordinate screen passes | Reviewer: the two rules disagree | The R4 rule cannot detect positional nonspecificity of this size |
| 11 | Visium donors | Label check, then lead-gene replication | The label check failed; a centred table was added post hoc | – | Not used in the main text (Supplementary Note 6) |
| 12 | Lead-gene stories | Stories promoted by pre-specified verdicts | Drug handling was promoted against the rule; under symmetric classes it is a supplementary observation | – | Supplementary Note 9 |
| 13 | Primary pathway list | Robust list on the selected coordinate | Coordinate-robust intersection (22) | Item 1 | Main text R5 |
| 14 | Class enrichments | Enrichment of each class (q ≤ 0.10, robust to expression) | A naming rule (also ≥ 3 class genes carrying the enrichment) was declared after the all-gene run had been seen, and before the probe-balanced run | Single-gene enrichments | Only tyrosine catabolism is named for mouse-only genes |
| 15 | External replication | Mann–Whitney comparisons of pathway z | Replaced by an overlap-preserving collection null | Overlapping, selected pathways are not independent | Better than random only |
| 16 | Agreement by gradient strength | Rule (ii): the weakest 80% of genes are reproducible within species but not shared. If it fails, the verdict is "agreement is confined to genes that are reliably zonated, and few genes are" | Rule (ii) held for S3 − early and S3c − early and could not be evaluated for S2 − S1 (human ceiling not estimable from half the atlas donors), so the verdict applied. v0.5 quoted it without "and few genes are"; v0.6 quotes it in full and states that it holds only for S2 − S1 at half-atlas resolution. The title states the concentration result of rule (i) instead of the verdict | Human S2 reference too small when halved; the verdict does not describe the S3 contrasts | Title and R2 wording |
| 17 | Agreement by gradient strength | Strength cross-fitted on both references | A post hoc scheme defined strength on mouse data alone, so the full atlas could give the human ceiling (weakest 80%, S2 − S1: ceiling 0.68, index −0.11) | Rule (ii) was not estimable for S2 − S1 | Reported as post hoc (R2) |
| 18 | Index sensitivity | "Largely different genes" holds if every upper bound for S2 − S1 and S3 − early is < 0.5 under every setting | Failed: module-block upper bounds reach 0.50–0.54 for the S3 contrasts | – | The wording "largely different genes" is not used; agreement is described as low and concentrated |
| 19 | Primary-list replication | Overlap-preserving test on the robust list | Repeated on the 22 primary pathways, which were defined after the coordinate results | Primary list changed after review | Replication beyond random sets, and beyond relabel-called pathways in one of two atlases (program level) |
| 20 | Conventional-only replication | – | An earlier overlap-naive group claim was tested against overlap-preserving collections and dropped (p = 0.41 and 0.24) | Overlap | No group claim |

---

## Supplementary Note 9 · Supplementary gene panels

Classes are the symmetric union classes (main text R3). External values are S3 − S1 (log2) unless
stated. Panels are in Supplementary Figure S11.

**Glutathione synthesis** (Gclc and Gclm indeterminate, Gclc with unequal probes, 3 mouse against 6
human; Gss conserved).
- *Mouse.* Gclc, Gclm and Gss rise toward late PT: snRNA +1.9 to +3.3 and microdissection +2.0 to +3.7,
  in both sexes.
- *Rat protein* agrees: S3/S2 +2.8, +1.9 and +2.4.
- *Human.* The genes differ in the cortex atlas: GCLC +0.61, GCLM −0.76 (−1.05 in male donors) and GSS
  +0.67. Human GCLM therefore falls where mouse rises, while GSS is conserved.
- *Pathway.* Glutathione metabolism is a robust pathway, but it is shared rather than mouse-led, and it
  loses significance with equal-probe genes.
- *Literature.* In rabbit PT, synthesis is highest in S1 and cellular GSH highest in S3 (Parks et al.
  1998). In rat, glutathione is highest in convoluted and early straight PT (Brehe et al. 1976).

**Cyp2e1** (not classifiable).
- *Mouse.* Cyp2e1 is confined to S1–S2 in our data and in microdissection. Rat protein peaks in S2
  (S3/S2 −5.5).
- *Literature.* Mouse renal CYP2E1 is a known proximal-tubular, male and testosterone-regulated enzyme
  (Hu et al. 1990; Speerschneider et al. 1995).
- *Human.* CYP2E1 is barely expressed: detected in 1.6% of our human structures, with at most 3 CPM in
  the cortex atlas. The difference is one of expression, not position. The literature on human
  protein conflicts (Speerschneider et al. 1995; Akakpo et al. 2023), so we state low mRNA, not absent
  protein.

**UGT1A9** (expressed in human only). It is highest in early PT (atlas −0.67). UGT1A9 is among the most
abundant renal UGTs in human (Margaillan et al. 2015). UGT1A isoforms share exons 2–5, and we did not
check whether the panel's three UGT1A9 probes lie in the UGT1A9-specific exon 1. The signal may
therefore include other UGT1A transcripts.

**AOX1** (expressed in human only). It rises toward late human PT (atlas +1.13). The evidence is
supporting only, because the probe counts differ (2 human against 3 mouse). Humans carry one AOX gene
where mice express four (Terao et al. 2016).

**Sterol synthesis.** Hmgcr is conserved, Ebp mouse-only, and Hmgcs1, Cyp51 and Lss indeterminate.
Mouse sterol zonation varies with time of day (Supplementary Note 4).

**Mitochondrial β-oxidation, literature.** In rat, 3-hydroxyacyl-CoA dehydrogenase activity is
roughly uniform along cortical segments (Le Hir et al. 1982; Bastin et al. 1990), so activity assays
would not see gene-specific zonation such as that of Acadm and Acaa2 (main text R6).

**Peroxisomal and other fatty-acid oxidation genes.** Acsm3, Crot and Nudt19 rise in mouse but are
indeterminate. Acsm3's late rise is a male-mouse program. Hsd17b4 is mouse-only with reviewed labels
and indeterminate with transferred labels. An intravital study found that late mouse PT mobilises lipid
droplets with lipase-rich lysosomes (Kaminska et al. 2026).

**Organic anion uptake** (Slc22a6 indeterminate; Slc13a3 conserved).
- *OAT1 (Slc22a6)* peaks in S2 in mouse (970 TPM in the cortical medullary-ray segment) and in rat
  protein (S3/S2 −6.5). Human OAT1 mRNA shows no clear gradient, whereas human OAT1 protein is reported
  stronger in S1/S2 than in medullary-ray S3 (Breljak et al. 2016).
- *NaDC3 (Slc13a3)* falls toward late PT in both species in S3c − early. In S2 − S1 it rises in mouse
  only.
- *Whole-PT averages* (Supplementary Figure S11g, h). Slc22a6, Slc13a3 and Cyp24a1 have non-significant
  whole-PT DESeq2 differences but significant segment differences of opposite sign. The S3 bars compare
  human cortex S3 with mouse S3 that includes the outer stripe.

**Early-PT transporters in late PT.** Slc9a3 is human-only in S2 − S1 and changes class with transferred
labels. Slc4a4 is conserved and Slc6a19 indeterminate. Mouse B0AT1's S1–S2 restriction is known (Navarro
Garrido et al. 2022).

**Genes flagged by interaction sign alone.**

| Gene | Symmetric class |
|---|---|
| Pah, Cyp24a1 | conserved |
| Igfbp4 | mouse-only |
| Acox2, Glyat, Nt5e, Cyp7b1 | indeterminate |
| Vnn1 | reversal (S2 − S1) |

**Claims that did not hold.**
- *Serine synthesis is not human-specific.* Psat1 is conserved, with unequal probes. Mouse Phgdh and
  Psat1 each have one probe against three in human, and mouse microdissection shows both
  S1-restricted, as in human.
- *Not replicated by the pre-specified checks:* Gamt (in our data), Hadh, Ephx1 and Me1.

---

## Supplementary Note 10 · Agreement by gradient strength and sensitivity of the conservation index

Notebook 49 (Supplementary Methods SM34) first reproduced the published conservation index exactly.
Then:

**Index sensitivity** (Tables SN10a–c).
- *Winsorisation.* The published estimator winsorises gradients at the 1st and 99th percentiles, which
  clips the strongest gradients. Our index ranged from 0.09 to 0.35 across winsorisation settings and
  reliability conventions, and the public-atlas index from 0.25 to 0.31.
- *Gene bootstrap.* It treats genes as independent, and every upper bound stayed below 0.5.
- *Module-block bootstrap.* It resamples co-expression modules of our PT structures. Its upper bounds
  reached 0.50–0.54 for the S3 contrasts and 0.44 for S2 − S1. The public-atlas index stayed below
  0.46.
- *Donor bootstrap.* Resampling 7 atlas donors with replacement duplicates donors. That overstates
  reliability and shifts the intervals relative to the point estimate, so the leave-one-donor-out range
  is the more useful check.
- *Pre-specified decision.* "Largely different genes" requires every upper bound below 0.5 under every
  setting, so it does not hold (Supplementary Note 8, deviation 18).

**Agreement by strength** (Tables SN10d–i).
- *Design.* Strength was cross-fitted: defined on one half of the external donors and evaluated on the
  other.
- *Within-bin values.* Reliabilities were computed on the winsorised values, so the all-gene values
  differ slightly from the published index (0.16 against 0.17 for S2 − S1).
- *Estimability.* A dash means a correlation could not be estimated: the bin's reliability was below
  0.2 in at least one fold, or more than 10% of draws were not estimable.
- *Module-block intervals.* They are wide in the top bins (Table SN10d): −0.04 to 0.82, −0.02 to 0.74
  and 0.06 to 0.79 for the top 1%.
- *Precision of the strength references.* The human half-reference has 3–4 atlas donors. Its
  reliability is 0.10–0.78 below the 95th percentile, against 0.67–0.98 for the mouse half-reference,
  so the strength axis is set mainly by mouse gradients.
- *Bin classes.* "Reproducible, not shared" means both ceilings' lower bounds ≥ 0.3 and the index's
  upper bound < 0.5. "Not reliably zonated" means a ceiling was not estimable or its upper bound was
  below 0.3. "Uncertain" covers the remaining bins, where the index's upper bound is ≥ 0.5.

**Table SN10a · Index under each winsorisation and reliability setting** (gene bootstrap point estimates).

| Contrast | Version | Reliability | q = 0 | q = 0.01 | q = 0.05 |
|---|---|---|---:|---:|---:|
| S2 − S1 | ours | unwinsorised (published) | 0.25 | 0.17 | 0.10 |
| S2 − S1 | ours | matched | 0.25 | 0.17 | 0.09 |
| S2 − S1 | public atlases | unwinsorised (published) | 0.30 | 0.27 | 0.25 |
| S2 − S1 | public atlases | matched | 0.30 | 0.28 | 0.28 |
| S3 − early | ours | unwinsorised (published) | 0.35 | 0.28 | 0.19 |
| S3 − early | ours | matched | 0.35 | 0.28 | 0.19 |
| S3 − early | public atlases | unwinsorised (published) | 0.31 | 0.28 | 0.25 |
| S3 − early | public atlases | matched | 0.31 | 0.28 | 0.26 |
| S3c − early | ours | unwinsorised (published) | 0.34 | 0.29 | 0.21 |
| S3c − early | ours | matched | 0.34 | 0.29 | 0.21 |
| S3c − early | public atlases | unwinsorised (published) | 0.29 | 0.27 | 0.26 |
| S3c − early | public atlases | matched | 0.29 | 0.27 | 0.25 |

**Table SN10b · 95% intervals at the published setting** (q = 0.01, unwinsorised reliability).

| Contrast | Version | Gene | Module k = 50 | Module k = 100 | Module k = 200 | Donor | Donor + gene |
|---|---|---|---|---|---|---|---|
| S2 − S1 | ours (0.17) | 0.12 to 0.22 | −0.10 to 0.41 | −0.07 to 0.40 | 0.00 to 0.34 | 0.17 to 0.23 | 0.13 to 0.26 |
| S2 − S1 | public atlases (0.27) | 0.22 to 0.33 | 0.10 to 0.40 | 0.11 to 0.40 | 0.16 to 0.38 | 0.18 to 0.29 | 0.15 to 0.32 |
| S3 − early | ours (0.28) | 0.22 to 0.33 | −0.02 to 0.53 | 0.01 to 0.51 | 0.08 to 0.49 | 0.28 to 0.33 | 0.23 to 0.37 |
| S3 − early | public atlases (0.28) | 0.23 to 0.33 | 0.11 to 0.42 | 0.12 to 0.42 | 0.15 to 0.42 | 0.22 to 0.31 | 0.20 to 0.32 |
| S3c − early | ours (0.29) | 0.23 to 0.35 | −0.03 to 0.54 | 0.01 to 0.51 | 0.08 to 0.52 | 0.29 to 0.35 | 0.23 to 0.38 |
| S3c − early | public atlases (0.27) | 0.21 to 0.32 | 0.08 to 0.44 | 0.10 to 0.43 | 0.13 to 0.44 | 0.21 to 0.31 | 0.18 to 0.33 |

Across all settings, the largest upper bound of our index was 0.44 (S2 − S1), 0.54 (S3 − early) and 0.54 (S3c − early); for public atlases alone, 0.45.

**Table SN10c · Dependence on the human ceiling.**

| Contrast | Human ceiling | Leave-one-donor-out | Index, leave-one-donor-out | Ceiling at which the index's upper bound reaches 0.5 | Human datasets giving a ceiling |
|---|---:|---|---|---:|---:|
| S2 − S1 | 0.65 | 0.63–0.68 | 0.166–0.173 | 0.13 | 1 |
| S3 − early | 0.67 | 0.65–0.68 | 0.275–0.281 | 0.30 | 2 |
| S3c − early | 0.66 | 0.65–0.68 | 0.285–0.292 | 0.33 | 2 |

**Table SN10d · Agreement by strength** (primary scheme: larger of the human and mouse ranks, cross-fitted; q = 0.01). The index is given with gene-bootstrap and module-block 95% intervals; a dash means not estimable (reliability below 0.2 in at least one fold, or more than 10% of draws not estimable).

| Contrast | Strength bin | Genes | Index (gene interval) | Module-block interval | Cross r* | Human ceiling | Mouse ceiling | Index, public atlases |
|---|---|---:|---|---|---:|---:|---:|---:|
| S2 − S1 | 0–50% | 3443 | – | – | −0.13 | – | 0.65 | – |
| S2 − S1 | 50–80% | 2066 | – | – | 0.01 | – | 0.84 | – |
| S2 − S1 | 80–90% | 688 | 0.24 (0.14 to 0.33) | 0.03 to 0.41 | 0.20 | 0.73 | 0.91 | 0.23 |
| S2 − S1 | 90–95% | 344 | 0.22 (0.10 to 0.32) | −0.01 to 0.45 | 0.19 | 0.83 | 0.92 | 0.32 |
| S2 − S1 | 95–99% | 276 | 0.33 (0.18 to 0.44) | −0.01 to 0.53 | 0.29 | 0.85 | 0.93 | 0.27 |
| S2 − S1 | 99–100% | 69 | 0.57 (0.34 to 0.74) | −0.04 to 0.82 | 0.53 | 0.90 | 0.95 | 0.39 |
| S2 − S1 | weakest 80% | 5509 | – | – | −0.05 | – | 0.77 | – |
| S2 − S1 | top 5% | 345 | 0.43 (0.32 to 0.54) | 0.03 to 0.63 | 0.39 | 0.89 | 0.94 | 0.34 |
| S2 − S1 | all genes | 6886 | 0.16 (0.12 to 0.21) | – | 0.13 | 0.75 | 0.89 | 0.29 |
| S3 − early | 0–50% | 3358 | – | – | −0.07 | – | 0.65 | – |
| S3 − early | 50–80% | 2014 | 0.09 (0.02 to 0.17) | −0.13 to 0.29 | 0.07 | 0.64 | 0.83 | 0.21 |
| S3 − early | 80–90% | 671 | 0.26 (0.17 to 0.34) | 0.06 to 0.47 | 0.21 | 0.73 | 0.89 | 0.21 |
| S3 − early | 90–95% | 336 | 0.27 (0.15 to 0.39) | −0.01 to 0.54 | 0.23 | 0.75 | 0.91 | 0.20 |
| S3 − early | 95–99% | 268 | 0.40 (0.27 to 0.51) | 0.04 to 0.61 | 0.35 | 0.84 | 0.90 | 0.26 |
| S3 − early | 99–100% | 68 | 0.58 (0.41 to 0.73) | −0.02 to 0.74 | 0.54 | 0.93 | 0.93 | 0.56 |
| S3 − early | weakest 80% | 5372 | 0.03 (−0.03 to 0.08) | −0.18 to 0.22 | 0.02 | 0.61 | 0.78 | 0.18 |
| S3 − early | top 5% | 336 | 0.49 (0.36 to 0.59) | 0.10 to 0.66 | 0.44 | 0.89 | 0.91 | 0.40 |
| S3 − early | all genes | 6715 | 0.27 (0.22 to 0.33) | 0.01 to 0.50 | 0.22 | 0.74 | 0.88 | 0.28 |
| S3c − early | 0–50% | 2391 | – | – | −0.08 | – | 0.65 | – |
| S3c − early | 50–80% | 1434 | 0.10 (0.02 to 0.20) | −0.14 to 0.29 | 0.07 | 0.62 | 0.86 | 0.23 |
| S3c − early | 80–90% | 478 | 0.21 (0.10 to 0.32) | −0.01 to 0.44 | 0.18 | 0.74 | 0.91 | 0.18 |
| S3c − early | 90–95% | 239 | 0.30 (0.14 to 0.43) | 0.00 to 0.62 | 0.25 | 0.73 | 0.92 | 0.21 |
| S3c − early | 95–99% | 192 | 0.37 (0.19 to 0.49) | −0.06 to 0.61 | 0.32 | 0.81 | 0.93 | 0.28 |
| S3c − early | 99–100% | 48 | 0.55 (0.39 to 0.71) | 0.06 to 0.79 | 0.53 | 0.94 | 0.96 | 0.47 |
| S3c − early | weakest 80% | 3825 | 0.03 (−0.04 to 0.10) | −0.20 to 0.22 | 0.02 | 0.60 | 0.80 | 0.21 |
| S3c − early | top 5% | 240 | 0.46 (0.34 to 0.55) | 0.09 to 0.64 | 0.42 | 0.88 | 0.94 | 0.36 |
| S3c − early | all genes | 4782 | 0.29 (0.22 to 0.35) | 0.01 to 0.51 | 0.23 | 0.73 | 0.90 | 0.27 |

Reliability of the half-references (fold-averaged), per bin:

| Contrast | 0–50% | 50–80% | 80–90% | 90–95% | 95–99% | 99–100% |
|---|---|---|---|---|---|---|
| S2 − S1 (human / mouse) | 0.10 / 0.67 | 0.24 / 0.88 | 0.41 / 0.94 | 0.44 / 0.96 | 0.68 / 0.97 | 0.88 / 0.98 |
| S3 − early (human / mouse) | 0.23 / 0.76 | 0.47 / 0.92 | 0.68 / 0.96 | 0.75 / 0.97 | 0.87 / 0.98 | 0.96 / 0.99 |
| S3c − early (human / mouse) | 0.24 / 0.78 | 0.46 / 0.92 | 0.69 / 0.96 | 0.78 / 0.98 | 0.86 / 0.98 | 0.96 / 0.99 |

**Table SN10e · Bin classification** (rule (ii) applied to each bin).

| Contrast | Bin | Genes | Gene bootstrap | Module-block bootstrap |
|---|---|---:|---|---|
| S2 − S1 | 0–50% | 3,443 | not reliably zonated | not reliably zonated |
| S2 − S1 | 50–80% | 2,066 | not reliably zonated | not reliably zonated |
| S2 − S1 | 80–90% | 688 | reproducible, not shared | reproducible, not shared |
| S2 − S1 | 90–95% | 344 | reproducible, not shared | reproducible, not shared |
| S2 − S1 | 95–99% | 276 | reproducible, not shared | uncertain |
| S2 − S1 | 99–100% | 69 | uncertain | uncertain |
| S3 − early | 0–50% | 3,358 | not reliably zonated | not reliably zonated |
| S3 − early | 50–80% | 2,014 | reproducible, not shared | reproducible, not shared |
| S3 − early | 80–90% | 671 | reproducible, not shared | reproducible, not shared |
| S3 − early | 90–95% | 336 | reproducible, not shared | uncertain |
| S3 − early | 95–99% | 268 | uncertain | uncertain |
| S3 − early | 99–100% | 68 | uncertain | uncertain |
| S3c − early | 0–50% | 2,391 | not reliably zonated | not reliably zonated |
| S3c − early | 50–80% | 1,434 | reproducible, not shared | reproducible, not shared |
| S3c − early | 80–90% | 478 | reproducible, not shared | reproducible, not shared |
| S3c − early | 90–95% | 239 | reproducible, not shared | uncertain |
| S3c − early | 95–99% | 192 | reproducible, not shared | uncertain |
| S3c − early | 99–100% | 48 | uncertain | uncertain |

**Table SN10f · Direction agreement by strength** (genes zonated in both species; share conserved, Wilson 95% interval; raw sign agreement over all genes).

| Contrast | Bin | Zonated in both | Share conserved (95% interval) | Raw sign agreement, all genes |
|---|---|---:|---|---:|
| S2 − S1 | 0–50% | 0 | – | 0.48 |
| S2 − S1 | 50–80% | 14 | 0.57 (0.33 to 0.79) | 0.50 |
| S2 − S1 | 80–90% | 34 | 0.59 (0.42 to 0.74) | 0.55 |
| S2 − S1 | 90–95% | 27 | 0.70 (0.52 to 0.84) | 0.57 |
| S2 − S1 | 95–99% | 44 | 0.86 (0.73 to 0.94) | 0.58 |
| S2 − S1 | 99–100% | 27 | 0.89 (0.72 to 0.96) | 0.67 |
| S2 − S1 | all genes | 146 | 0.75 (0.67 to 0.81) | 0.51 |
| S3 − early | 0–50% | 1 | 1.00 (0.21 to 1.00) | 0.61 |
| S3 − early | 50–80% | 52 | 0.48 (0.35 to 0.61) | 0.57 |
| S3 − early | 80–90% | 66 | 0.68 (0.56 to 0.78) | 0.60 |
| S3 − early | 90–95% | 54 | 0.67 (0.53 to 0.78) | 0.62 |
| S3 − early | 95–99% | 78 | 0.74 (0.64 to 0.83) | 0.64 |
| S3 − early | 99–100% | 38 | 0.84 (0.70 to 0.93) | 0.72 |
| S3 − early | all genes | 289 | 0.68 (0.63 to 0.73) | 0.60 |
| S3c − early | 0–50% | 0 | – | 0.63 |
| S3c − early | 50–80% | 19 | 0.58 (0.36 to 0.77) | 0.58 |
| S3c − early | 80–90% | 38 | 0.74 (0.58 to 0.85) | 0.60 |
| S3c − early | 90–95% | 32 | 0.75 (0.58 to 0.87) | 0.67 |
| S3c − early | 95–99% | 50 | 0.76 (0.63 to 0.86) | 0.67 |
| S3c − early | 99–100% | 27 | 0.89 (0.72 to 0.96) | 0.79 |
| S3c − early | all genes | 166 | 0.75 (0.68 to 0.81) | 0.62 |

**Table SN10g · Strength from one species only** (weakest 80% and top 5%; index with gene-bootstrap 95% interval). These schemes were pre-specified as secondary. The contrast between the human-ranked and mouse-ranked weakest 80% (0.18 against 0.00 for S3 − early) is not interpretable as a species asymmetry: in a reviewer check, unequal precision of the two strength references reproduced it without any biological difference.

| Contrast | Scheme | Weakest 80%: index | Human ceiling | Mouse ceiling | Top 5%: index |
|---|---|---|---:|---:|---|
| S2 − S1 | larger of both ranks (primary) | – | – | 0.77 | 0.43 (0.32 to 0.54) |
| S2 − S1 | human rank only | – | – | 0.88 | 0.53 (0.40 to 0.65) |
| S2 − S1 | mouse rank only | – | – | 0.72 | 0.47 (0.34 to 0.57) |
| S2 − S1 | mean of both ranks | – | – | 0.82 | 0.53 (0.42 to 0.64) |
| S2 − S1 | our mouse gradient (absolute), same data | −0.06 (−0.12 to 0.00) | 0.69 | 0.71 | 0.46 (0.35 to 0.56) |
| S3 − early | larger of both ranks (primary) | 0.03 (−0.03 to 0.08) | 0.61 | 0.78 | 0.49 (0.36 to 0.59) |
| S3 − early | human rank only | 0.18 (0.12 to 0.25) | 0.59 | 0.87 | 0.50 (0.33 to 0.62) |
| S3 − early | mouse rank only | 0.00 (−0.06 to 0.06) | 0.70 | 0.71 | 0.62 (0.51 to 0.70) |
| S3 − early | mean of both ranks | 0.10 (0.04 to 0.16) | 0.62 | 0.82 | 0.56 (0.43 to 0.65) |
| S3 − early | our mouse gradient (absolute), same data | 0.03 (−0.02 to 0.09) | 0.67 | 0.70 | 0.53 (0.39 to 0.63) |
| S3c − early | larger of both ranks (primary) | 0.03 (−0.04 to 0.10) | 0.60 | 0.80 | 0.46 (0.34 to 0.55) |
| S3c − early | human rank only | 0.19 (0.12 to 0.26) | 0.59 | 0.88 | 0.46 (0.28 to 0.58) |
| S3c − early | mouse rank only | 0.02 (−0.05 to 0.10) | 0.68 | 0.72 | 0.57 (0.44 to 0.68) |
| S3c − early | mean of both ranks | 0.09 (0.03 to 0.16) | 0.63 | 0.83 | 0.50 (0.36 to 0.59) |
| S3c − early | our mouse gradient (absolute), same data | 0.05 (−0.03 to 0.13) | 0.66 | 0.77 | 0.48 (0.35 to 0.58) |

**Table SN10h · Post hoc: strength from mouse snRNA alone, full cortex atlas as human reference.**

| Contrast | Bin | Index (95% interval) | Human ceiling | Mouse ceiling |
|---|---|---|---:|---:|
| S2 − S1 | 0–50% | −0.23 (−0.31 to −0.15) | 0.66 | 0.56 |
| S2 − S1 | 50–80% | −0.02 (−0.09 to 0.05) | 0.70 | 0.81 |
| S2 − S1 | weakest 80% | −0.11 (−0.16 to −0.05) | 0.68 | 0.72 |
| S2 − S1 | top 5% | 0.46 (0.33 to 0.57) | 0.89 | 0.94 |
| S3 − early | 0–50% | −0.06 (−0.14 to 0.02) | 0.65 | 0.55 |
| S3 − early | 50–80% | 0.04 (−0.04 to 0.11) | 0.70 | 0.79 |
| S3 − early | weakest 80% | 0.00 (−0.06 to 0.06) | 0.67 | 0.71 |
| S3 − early | top 5% | 0.62 (0.51 to 0.70) | 0.83 | 0.92 |
| S3c − early | 0–50% | −0.10 (−0.22 to 0.03) | 0.66 | 0.53 |
| S3c − early | 50–80% | 0.09 (0.00 to 0.18) | 0.66 | 0.81 |
| S3c − early | weakest 80% | 0.02 (−0.05 to 0.10) | 0.66 | 0.72 |
| S3c − early | top 5% | 0.57 (0.44 to 0.68) | 0.84 | 0.94 |

**Table SN10i · Public atlases alone, on our strength bins** (noise-corrected correlation of the cortex atlas with male mouse snRNA).

| Contrast | Weakest 80% | Top 5% | Top 1% | All genes |
|---|---:|---:|---:|---:|
| S2 − S1 | 0.08 | 0.39 | 0.56 | 0.26 |
| S3 − early | 0.15 | 0.44 | 0.55 | 0.23 |
| S3c − early | 0.13 | 0.36 | 0.40 | 0.22 |

Cortex-atlas healthy donors had 59–780 S2 nuclei each (S1 95–1648; S3 54–608).

---

## Supplementary figure legends

Panel sources are in `figure_plan_v04.md`. The standing caveat applies to every figure.

**Supplementary Figure S1 | The continuous coordinate: construction and validation.**
- **a,** The pass-2 Harmony embedding (dimensions 1 and 2) of PT structures, coloured by segment label,
  with a coordinate trace (median position in 30 equal-count bins). The curve is fitted in five
  dimensions, so its trace can loop in this projection.
- **b,** Fold-protected agreement of the coordinate (filled) and of segment labels (open) with external
  S3-versus-early contrasts.
- **c,** Identifiability:
  - absorption of an injected human-only gradient;
  - sign agreement for 244 externally species-specific genes (coordinate 91%, labels 95%);
  - count splitting.
- **d,** Equal-depth thinning: agreement with the original order and the correlation of position with
  library size.
- **e,** Within-segment agreement between gene-fold coordinates, and the numerical floor.
- **f,** Held-out gene prediction against segment and coverage oracles.

**Supplementary Figure S2 | Pathway calls under three coordinates.**
- **a,** Robust and core counts under the principal curve, diffusion pseudotime and the equal-depth
  refit, and their overlaps (36, 26 and 51 robust; 22 shared by the first two).
- **b,** Peak positions in segment units under each coordinate, for human-high and mouse-high
  pathways.
- **c,** Within-species positional calls per coordinate.
- **d,** For each robust list, the share also called by the coordinate-free step screen, by the
  spline-based step screen, and by position beyond the steps.

**Supplementary Figure S3 | Specificity details in our cohort.**
- **a,** All strategies (Table S1): species and relabeled calls, with NCDR.
- **b,** Decoy-estimated FDP as the threshold is lowered; the line marks 10%.
- **c,** Within-species controls: offset and positional calls for the species split, mouse against
  mouse, and human section against section.
- **d,** Mean decoy z² against pathway size, with fitted inflation for the offset and position
  statistics.
- **e,** Specimen-level split-plot gene test: Q–Q plot for the species split and both relabelings.

**Supplementary Figure S4 | The 36 pathways robust on the principal curve.** Median member Δ(s) along
the coordinate in 17 programs. ★ marks core pathways and ● the 22 primary pathways. Peaks are marked
in segment units.

**Supplementary Figure S5 | Zonation amplitude.**
- **a,** Projection ratios on each reference axis (OLS, Deming, orthogonal).
- **b, c,** Gradients against male mouse microdissection (b) and the human cortex atlas (c).
- **d,** Measurement checks.
- **e,** Tissue-state ratios.
- **f,** Noise-corrected spread and disattenuated correlations.

**Supplementary Figure S6 | Physiological and tissue state.**
- **a,** Whole-kidney fasting and time-of-day effects against the species × position interaction, per
  program.
- **b,** Mouse S3 − S1 at ZT4 against ZT16.
- **c,** Our donor on the cortex-atlas injury axis.
- **d,** Injury and amplitude among atlas donors.

**Supplementary Figure S7 | Genome-wide agreement with the cortex atlas.**
- **a,** Our human S3 − S1 against the atlas.
- **b, c,** Our interaction against atlas − male (b) and atlas − female (c) mouse snRNA.

**Supplementary Figure S8 | Segment-label checks.**
- **a,** Glomerular distance by label, with morphological and expression-defined glomeruli. Human
  distances are shown at 0.5 µm per pixel (uncorrected) and at 0.44 µm per pixel [VERIFY].
- **b,** Glomerulus-attached PT at 2 µm and 10 µm.
- **c, d,** Reference-transfer confusion matrices for human (cortex atlas) and mouse (microdissection
  and snRNA).
- **e,** Spillover index per species and segment.

**Supplementary Figure S9 | Conservation index under transferred labels and without integration
genes.** The index with its 95% interval per setting, gene half and contrast. Settings: reviewed labels;
human labels from the atlas; human labels from the atlas with mouse labels from microdissection; each
with and without the 611 integration genes. The dashed line marks 0.5.

**Supplementary Figure S10 | Lead-gene verification.** One row per gene, with one column per
pre-specified test, the symmetric zonation class, probe counts and the rat protein S3/S2 ratio.

**Supplementary Figure S11 | Supplementary gene panels.**
- **a–f,** Cards as in Figure 6 for Gamt, Gclc, Gclm, Gss, Cyp2e1 and UGT1A9.
- **g, h,** Slc22a6, Slc13a3 and Cyp24a1: fitted curves (g) and DESeq2 log2 fold changes for whole PT
  and each segment (h). The S3 bars compare human cortex S3 with mouse S3 that includes the outer
  stripe.

**Supplementary Figure S12 | Agreement by strength under alternative strength schemes** (the
single-species schemes are secondary; their difference reflects the unequal precision of the strength
references, Supplementary Note 10). Index, human
and mouse ceilings, and cross-species correlation per strength bin, for strength defined from:
- the larger of both ranks (primary);
- the human rank only;
- the mouse rank only;
- the mean rank;
- our own mouse gradient (same data, showing the selection effect).

Post hoc full-atlas estimates are overlaid for the mouse-rank scheme (Table SN10h).

---

## Supplementary Methods

These Supplementary Methods give the full parameters behind the main Methods. Each subsection names
the notebook that implements it, and Table SM1 (end) maps every analysis to its notebook and results
folder. Parameters were read from the code.

### SM1 · Study design and specimens

**Cohort.**
- Mouse: two control mouse kidney sections, `Ctrl1A2` and `Ctrl1A4`.
- Human: two sections of renal cortex from one donor, `HUK1_COR1` and `HUK1_MED1`. The second is
  named "MED" at source, but its segmentation contains cortex only, and it is analysed as cortex.

All four specimens express Y-linked transcripts and are treated as male [VERIFY against records].
The mouse sections are whole-kidney sections, so they include cortex and outer medulla. The two
ischaemia–reperfusion kidneys of the same mouse experiment (`IR2A2`, `IR2A4`) are used only in the
tissue-state analyses (notebooks 39 and 44) and as the positive control of notebook 45 (SM14).
Specimens are the units of replication.

**Not yet recorded:**
- Mouse: strain (the paper plan says C57BL/6-type; no file records it), age, supplier, housing and
  euthanasia [VERIFY].
- Human: tissue source, donor age, kidney function, and whether the two sections are serial or come
  from separate blocks [VERIFY].
- Tissue processing: fixation and embedding, section thickness, H&E protocol and imaging system
  [VERIFY]. A CytAssist image exists for the human sections.
- Ethics: animal-protocol and human-tissue approvals [VERIFY].

### SM2 · Visium HD processing

Sections were profiled with Visium HD (Oliveira et al. 2025). Library preparation, sequencer and
depth are not recorded here [VERIFY]. Run metadata:

| | Mouse | Human |
|---|---|---|
| Space Ranger | 3.1.1 | 4.0.1 |
| Probe set | Visium Mouse Transcriptome v2.0 | Visium Human Transcriptome v2.1.0 |
| Reference | mm10 | GRCh38-2024-A |
| Chemistry | Visium HD v1 | Visium HD H1 slide, probe-based v1 |
| Features | 19,059 | 18,132 |

We used the 2-µm bin output (`square_002um`). The bin size was confirmed by the `s_002um_` barcode
prefix, and only bins under tissue were kept.

### SM3 · Tubule segmentation

**Model.** Tubules and glomeruli were segmented on H&E whole-slide images with `kidney_panoptic`
(repository `segmentation/`, checkpoint v4), at 0.44068 µm per pixel.

- *Encoder.* A frozen OpenMidnight encoder: a pathology-tuned DINOv2 ViT-g/14 with registers
  (Oquab et al. 2023; Kaplan et al. 2025). Tokens were taken from blocks 9, 19, 29 and 39, with
  518-px input.
- *Decoder.* A trainable convolutional stem on native pixels feeds a decoder at output stride 1:
  7.31 million parameters, GroupNorm, dropout 0.10.
- *Heads.* A 5-class semantic head (background, proximal, distal and collecting-duct tubule,
  glomerulus). A 4-class boundary head (background, interior, boundary to background, boundary to
  another instance). A centre heat map for glomeruli.

**Decoding.** Instances come from a marker-controlled watershed:

1. Foreground is 1 − p(background) ≥ 0.50.
2. Seeds are interior components with p(interior) ≥ 0.40 and area ≥ 24 px².
3. The watershed floods over the boundary probability.
4. Each instance takes its majority semantic class, with the boundary band excluded.
5. Per-class minimum areas (px²) are applied: proximal 1,950; distal 1,000; collecting duct 950;
   glomerulus 2,300.

Glomeruli with two or more centre peaks were split. Fragment merging and semantic splitting were
disabled.

**Whole-slide inference.**
- Tiles of 512 px with 128-px overlap were blended with a Gaussian window (σ = 0.25 × tile).
- Decode windows were 2,048 px with 256-px overlap, stitched at IoU ≥ 0.5.
- A tissue gate dropped instances lying more than half off tissue.
- D4 test-time augmentation averaged the probability maps over 8 orientations.
- Polygons were simplified at 1 px, with at most 400 vertices.

**Training and selection.**
- *Data.* 6,475 patches of 512 px (stride 256) from 7 partially annotated mouse slide regions
  (nephrectomy and control kidneys). Part of `Ctrl1A2` was in the training data.
- *Ignore mask.* Annotated instances were trusted foreground. A 5-px ring around each instance and
  off-tissue glass were trusted background. Unannotated tissue was ignored.
- *Loss.* Cross-entropy plus Dice on the semantic and boundary heads. Boundary pixels were weighted
  ×3 and instance-to-instance boundaries ×6. The centre head used mean squared error.
- *Optimiser.* AdamW, learning rate 3 × 10⁻⁴ with cosine decay, weight decay 0.03, batch 8, mixed
  precision. Augmentation was heavy: geometric, stain, blur, noise, JPEG, cutout, elastic and
  copy-paste.
- *Checkpoint selection.* The checkpoint was chosen by pooled class-agnostic panoptic quality (PQ)
  (Kirillov et al. 2019), never by validation loss. The selected checkpoint is epoch 35 of 54.
- *Held-out PQ.* On a fully held-out nephrectomy slide (98 objects), PQ was 0.837 without and 0.880
  with D4 test-time augmentation (one mouse slide; the human segmentation has no comparable
  evaluation).
- *Known errors.* Predicted objects are smaller than annotated ones, and merges of touching
  same-type tubules are the main error.

**Human segmentation and mouse quality control.**
- The human `*_v2.geojson` re-segmentations have no recorded provenance [VERIFY]. Their µm fields
  assume 0.5 µm per pixel, whereas the segmentation model works at 0.44068 µm per pixel (SM25)
  [VERIFY against scanner metadata].
- Mouse analyses use only the quality-controlled `*_kept_tubules_labeled_fine.geojson`. The rules
  and code of that upstream QC are not in the repository [VERIFY].
- Only polygon geometry is used downstream.

**Centroid check.** Every join between expression and segmentation recomputes the indexed polygon's
centroid and stops at any deviation above 1 px.

### SM4 · Tubule-by-gene matrices

Each 2-µm bin was assigned to the polygon covering its centre, in full-resolution image coordinates
(scale factor 1, shapely STRtree, at most one polygon per bin). Raw counts were summed per polygon.
The segmentation and the Visium HD image must share one coordinate frame. How they were registered
is not recorded [VERIFY]. Notebook 42 checked the result (SM25):
- every polygon centroid matches the centroid stored with the matrix (maximum offset 0 px);
- bins cover a median 97.4% of mouse polygon area and 76.5% of human polygon area.

The human coverage is close to the 75.7% expected if the human µm fields overstate the pixel size
(0.5 instead of 0.44068 µm), so it does not indicate a misregistration.

| Specimen | Polygons | Polygons with ≥ 1 bin |
|---|---:|---:|
| `Ctrl1A2` | 7,624 | 7,624 |
| `Ctrl1A4` | 7,002 | 7,002 |
| `HUK1_COR1` | 8,747 | 6,416 |
| `HUK1_MED1` | 14,186 | 7,223 |

Human polygons without bins (2,331 in `HUK1_COR1`, 6,963 in `HUK1_MED1`) lie almost entirely outside
the region covered by bins. Only 16 and 90 of them fall inside the convex hull of bin-covered
polygons, where they make up 0.2% and 1.2% of polygons. This is consistent with segmented tissue
beyond the Visium HD capture area. Mouse sections
have no polygons without bins (notebook 42, `registration_qc.csv`).

### SM5 · Cross-species ortholog space

**Ortholog map.** We used the HCOP human–mouse table (Yates et al. 2021; local copy dated 31 August 2026,
sha256 prefix `0cfb78e4eb273751`).

- Pairs supported by fewer than three databases were removed.
- Each pair was scored as 10 × (number of supporting databases) + 1 if the symbols are identical.
- Only mutual best pairs were kept. This gives a bijective map of 17,449 pairs.

**Mapping and measurement flags.** Human counts were moved to mouse-symbol columns by a sparse
transformation. Each gene is flagged as measured in every input or not. A column that no human
feature feeds is a structural zero, not an observation. Of the 17,449 pairs, 15,567 are measured in
all four specimens.

Gene-level analyses use measured genes only. The integration object keeps all columns, because
removing columns would change Scanpy's Seurat-flavour HVG binning and, through it, the reviewed
clustering.

### SM6 · Structure filters and normalisation for integration

**Structure filters** (reproduced from raw counts by notebook 13):
- Mouse structures (QC'd upstream) have no gene-count threshold.
- Human structures need ≥ 100 detected ortholog genes.
- The lowest 5% of bins per structure (`n_spots`) is removed within each species.

This retained 26,839 of 37,559 structures: `Ctrl1A2` 7,231; `Ctrl1A4` 6,676; `HUK1_COR1` 6,098;
`HUK1_MED1` 6,834.

**Gene filter and normalisation.**
- Genes needed detection in ≥ 5% of structures and ≥ 20 total counts, giving 10,076 genes (9,943
  measured in all inputs).
- Counts were scaled to 10⁴ over these genes and transformed with log1p.
- Mitochondrial and ribosomal genes (`mt-`, `Rpl`, `Rps`, `Mrpl`, `Mrps`) were excluded from
  integration features.

### SM7 · Integration, clustering and reviewed labels (pass 1)

**Integration.**
- *Features.* HVGs were selected within each species (Seurat flavour, batch-aware by specimen; mean
  0.0125–3; dispersion ≥ 0.5). Only genes selected in both species were kept: 732.
- *PCA.* 50 components on log-normalised values (Scanpy; Wolf et al. 2018), using all 732 genes.
- *Harmony* (Korsunsky et al. 2019) via rpy2, R harmony 2.0.5 enforced. Batch = specimen; θ = 6,
  λ = 1, at most 30 iterations, τ = 0, seed 0.

Specimen and species are confounded in this design. We therefore use the integrated embedding only
to group and order structures, and every gene comparison uses raw counts.

**Clustering.** Leiden clustering (Traag et al. 2019) on a 30-neighbour graph of the 50 Harmony
dimensions: igraph flavour, resolution 0.7, 2 iterations, seed 0. This gives 11 clusters. A
fingerprint (structures, parameters, number of clusters, SHA-1 of memberships) pins the partition,
and a changed fingerprint stops the run.

**Labels.** Clusters were labelled by manual review of marker detection and differential expression,
at the granularity the markers support. Reference panels:

| Segment | Markers |
|---|---|
| PT-S1 | Slc5a2, Slc5a12, Gatm, Lrp2, Cubn, Slc34a1 |
| PT-S2 | Slc22a6, Slc13a3, Cyp2e1 |
| PT-S3 | Slc22a7, Slc7a13, Cyp7b1, Slc6a18, Acsm3 |

The map was PT-S1, PT-S2, PT-S3, two ascending-limb clusters, DCT, two collecting-system clusters,
glomerulus, smooth muscle and one unresolved cluster. No cluster had thin-limb markers.

**What the labels are.**
- The S1/S2/S3 labels are memberships of joint cross-species Leiden clusters on the integrated
  embedding, named by marker review. They are not anatomical calls made separately in each species.
- How joint clustering merges or splits cell types across species depends on the integration
  strategy and its strength (Song et al. 2023).
- Notebook 42 tests the labels against distance to glomeruli and against expression-only reference
  transfer (SM25). Notebook 46 repeats the conservation index and class counts with transferred
  labels (SM28).
- θ = 6 was fixed in the pipeline (in notebook 03 since its first full version) before any analysis
  reported here; its original rationale is [VERIFY: authors]. An older diagnostic table that compared
  θ values is no longer produced and is not used. The labels were checked by reference
  transfer instead (SM25, SM28).

Notebook 13 repeats notebook 03's pass-1 clustering on the same 26,839 structures. The two
partitions differ in membership hash but agree closely:
- adjusted Rand index (ARI) 0.979 for the Leiden clusters and 0.981 for the reviewed labels;
- 99.1% of structures carry the same label.

Notebook 13 calls 12,866 structures PT and notebook 03 calls 12,871; 12,864 are PT in both. Of these,
98.5% have the same S1/S2/S3 label (ARI 0.955). Most disagreements are notebook-03 S1 structures that
notebook 13 labels S2 (125).

### SM8 · Nephron re-integration (pass 2) and PT subset

Glomerulus, smooth-muscle and unresolved clusters (2,504 structures) were removed. That left 24,335
tubular structures: PT 12,866; ascending limb 6,731; collecting system 3,363; DCT 1,375. Integration
was repeated on these structures:

- 752 species-intersected HVGs were selected.
- 50 principal components were computed on 611 of them: those also among pass 1's 2,000 Seurat HVGs.
  Scanpy 1.11 restricts PCA to the inherited `highly_variable` column.
- Harmony was run as in pass 1, followed by a 30-neighbour graph.

The PT subset (12,866 structures) was taken after pass 2:

| Specimen | S1 | S2 | S3 |
|---|---:|---:|---:|
| `Ctrl1A2` | 938 | 976 | 938 |
| `Ctrl1A4` | 873 | 1,017 | 918 |
| `HUK1_COR1` | 1,890 | 908 | 640 |
| `HUK1_MED1` | 2,076 | 905 | 787 |

### SM9 · PT coordinate

**Fit.** The PT coordinate is a nonbranching principal curve (scFates 1.2.5; Faure et al. 2023;
ElPiGraph, Albergante et al. 2020).
- Input: the first five pass-2 Harmony dimensions (computed from 611 PCA genes).
- Curve: 30 nodes, seed 0.
- The fit was required to be a single path with two tips and no forks.

**Orientation.**
- A marker axis was computed as the mean within-species z-score of late genes (Slc22a7, Slc7a13,
  Cyp7b1) minus that of early genes (Slc5a2, Slc5a12, Gatm).
- The curve was rooted at the tip with the higher early score, mapped once (`n_map = 1`) and scaled
  to [0, 1].
- The axis sets only the root. It is not a cross-species check: Gatm and Cyp7b1 are not zonated in
  human external data (human snRNA S3 − early −0.10 and +0.39).
- scFates' automatic tip rule sets interior nodes to 0. The per-species z-scored score avoids
  rooting at an interior node, and notebook 36 guards against it.

**Uncertainty.** Twenty resampled mappings moved positions by SD ≈ 0.007 and changed gene T_spatial
only slightly (ρ = 0.994). This is far less than the disagreement between coordinates built from
disjoint gene thirds, so `n_map` understates positional uncertainty.

**Required order.** In every specimen, the S3 median had to exceed the S1 and S2 medians.

**DPT comparator.** Scanpy diffusion pseudotime (Haghverdi et al. 2016) was run on the same five
dimensions: 30 neighbours, rooted among the structures most strongly assigned to the scFates root
node.

### SM10 · Coordinate validation and alternatives (notebooks 35 and 36)

All rules below were fixed before outputs were seen, except where marked.

**External concordance (P1, P1w).**
- Each gene's within-specimen Spearman correlation with the coordinate, and its correlation within
  reviewed segments, were compared across genes with external S3-versus-early or adjacent-segment
  contrasts.
- Each gene was scored on a coordinate refitted without its sha256 gene fold.
- External donors were split into an anchor-selection set (mouse microdissection, female mouse
  snRNA, human donor half A) and an evaluation set (male mouse snRNA, human donor half B).

**Further checks.**
- *Absorption.* A human-only multiplicative injection exp(±log 2·(2s − 1)) into 10% of the
  construction genes, followed by refitting.
- *Equal depth.* Binomial thinning of every structure to the pooled PT 20th-percentile library size.
- *Count splitting.* Poisson splitting (ε = 0.5) between building the coordinate and testing
  (Neufeld et al. 2023).
- *Mapping.* `n_map = 20` resampling.
- *Registration and seed stability.* These were added after notebook 35's outputs had been seen.

**Alternatives** (notebook 36). Four alternatives were compared under pre-specified rules:
- an equal-depth refit;
- DPT;
- a PT-only integration;
- a conserved-anchor coordinate (160 genes with the same external zonation direction in all
  selection datasets).

**Pre-registered selection rule.** A candidate was adopted only if it met all three:
- (i) equal-depth agreement P5 ≥ 0.80;
- (ii) fold-protected P1 ≥ the scFates value − 0.05 in each species;
- (iii) registration gap ≤ 0.15.

Among passing candidates, the one with the highest gene-fold within-segment agreement was chosen. If
none passed, scFates was kept.

**Outcome** (Supplementary Note 8, deviation 1).
- The rule selected DPT13 (P5 0.80; P1 0.74 mouse, 0.29 human; gap 0.14).
- The scFates equal-depth refit failed (P5 0.73; gap 0.36).
- scFates itself would also have failed (P5 0.72; gap 0.29).
- The rule omitted the pre-registered downstream specificity gate: criterion 5 of rule C6, which
  requires the T_spatial screen to keep NCDR < 0.25 and fewer than 5% of pathways called under
  relabeling.
- After DPT13 was seen to fail that gate on notebook 36's matched-only test (NCDR 0.40), we applied
  the gate to every candidate. No candidate then qualified, and scFates was kept as primary.
- Notebook 37's full pipeline was then run on DPT13 and on the scFates equal-depth refit (SCF13-ED;
  SM29). Under notebook 37's primary joint test, DPT13 is specific: 45 species calls against 3 and 0
  relabeled (NCDR 0.03). The matched-only test in the same run gives 99 calls against 41 and 36
  (NCDR 0.39). DPT13 gives 26 robust and 15 core pathways.
- Whether DPT13 passes the gate therefore depends on the test. After review, the primary pathway
  list was defined as the pathways robust under both scFates and DPT13 (22), with each coordinate
  reported separately as a sensitivity analysis.

**Outcomes.**
- *DPT* on notebook 13's embedding (notebook 36): 96 calls; 56 of 66 robust and 38 of 42 core
  pathways of the notebook 31–33 list retained; NCDR 0.40 (matched test).
- *DPT13 and SCF13-ED through notebook 37* (`d4_coordinate_runs.csv`):

  | | scFates | DPT13 | SCF13-ED |
  |---|---:|---:|---:|
  | Joint T_spatial: species calls (relabeled) | 60 (2, 0) | 45 (3, 0) | 66 (1, 0) |
  | Joint NCDR | 0.02 | 0.03 | 0.01 |
  | Matched T_spatial: species calls (relabeled), NCDR | 116 (17, 13), 0.13 | 99 (41, 36), 0.39 | 123 (19, 2), 0.09 |
  | Joint step screen: species (relabeled) | 16 (0, 0) | 14 (0, 0) | 31 (5, 6) |
  | Within-species positional calls (mouse, human) | 16, 0 | 37, 5 | 12, 3 |
  | Candidates / not specimen-sensitive / stable | 60 / 58 / 51 | 45 / 42 / 33 | 66 / 65 / 60 |
  | Robust / core | 36 / 26 | 26 / 15 | 51 / 36 |
  | Robust human-high / mouse-high at peak | 10 / 26 | 10 / 16 | 17 / 34 |
  | scFates robust kept as robust | – | 22 / 36 | 32 / 36 |
  | scFates core kept as core | – | 12 / 26 | 21 / 26 |

  In the 40-pathway union of scFates and DPT13 robust pathways, the direction at the divergence peak
  agrees across coordinates wherever the pathway is a candidate. For the 22 pathways robust under
  both, peak positions correlate with Spearman ρ = 0.71 (SM29).
- *Conserved-anchor coordinate:* 42 of 66 and 30 of 42 retained, on 7,167 structures, because it
  did not register the species.
- Refits whose curve could not be rooted were reported as failed, not retuned.

### SM11 · PT analysis set and gene universe

**Common support.** Gene models use the intersection of each specimen's 1st–99th percentile
coordinate interval, [0.018, 0.864]. This keeps 12,272 PT structures:

| Specimen | Structures | S1 | S2 | S3 |
|---|---:|---:|---:|---:|
| `Ctrl1A2` | 2,552 | 823 | 958 | 771 |
| `Ctrl1A4` | 2,685 | 821 | 1,014 | 850 |
| `HUK1_COR1` | 3,355 | 1,845 | 904 | 606 |
| `HUK1_MED1` | 3,680 | 2,041 | 902 | 737 |

The evaluation grid has 61 points, and every specimen needs ≥ 15 structures within ±8% of the
interval around each point.

**Expression.** Expression was rebuilt from raw counts for all 17,449 pairs:
- Library size is the total over the 15,567 orthologs measured in all specimens.
- Values are log1p(count / library size × 10⁴).
- A gene is eligible if measured in all specimens with mean equal-specimen detection ≥ 2%. This gives
  11,071 genes.

### SM12 · Nested gene models

**Models.** For each eligible gene, three nested Gaussian regression-spline models were fitted by
weighted least squares on coordinate s:

| Model | Terms |
|---|---|
| M₀ | f(s) + specimen contrasts |
| M_level | adds a human indicator |
| M_full | adds human × g(s) |

**Specification.**
- f and g share a fixed, unpenalised cubic B-spline basis: 6 columns plus intercept, internal knots
  at the 25th, 50th and 75th percentiles, boundaries 0 and 1.
- Specimen intercepts are contrasts that sum to zero within species.
- Weights give each species half the total weight, shared equally between its specimens:
  wᵢ = n / (2 · J_s · n_j).

**Statistics.**
- Partial F statistics: T_level (M₀ vs M_level, 1 df), T_spatial (M_level vs M_full, 6 df) and
  T_total (7 df). They rank genes. They are not tests, because their residual df count structures.
- Count-based trajectory tests such as tradeSeq and condiments (Van den Berge et al. 2020; Roux de
  Bézieux et al. 2024) test conditions along a trajectory with cells as units. With two specimens
  per species, that would treat structures as replicates. We therefore use the gene statistics only
  as ranks, and place inference at the pathway level under specimen relabeling (SM14) and at the
  specimen level (SM18).
- The species difference Δ(s) = μ̂_human − μ̂_mouse has HC3 standard errors (MacKinnon & White 1985),
  and Z(s) = Δ/SE.

### SM13 · Pathway libraries and the joint test

**Libraries.**
- Reactome 2022 (Gillespie et al. 2022), MSigDB Hallmark 2020 (Liberzon et al. 2015) and KEGG 2019
  Mouse (Kanehisa et al. 2019), as Enrichr library files (Kuleshov et al. 2016; local copies dated
  31 August 2026).
- Members were mapped through the ortholog table, and sets with 10–300 eligible members were tested:
  1,513 pathways (Reactome 1,181, Hallmark 50, KEGG 282).

**Matched test.**
- Rank-AUC of members' T_spatial against all other eligible genes, compared with 9,999 random sets of
  equal size. The random sets were drawn without replacement within tertile strata of equal-specimen
  mean expression, detection and positional coverage (the fraction of 10 coordinate bins with
  detection).
- Seed 12. The matched z is (AUC − null mean) / null SD.

**Correlation correction.** Following CAMERA (Wu & Smyth 2012), each set's mean pairwise residual
correlation ρ was estimated from M_full residuals:
- computed within each specimen;
- Fisher-averaged with equal weights;
- floored at 0.

This gives a variance inflation factor VIF(m, ρ) for the rank sum.

**Joint test.**
- Joint z = matched z / √VIF, with a one-sided normal p and BH over all 1,513 sets (Benjamini &
  Hochberg 1995).
- A candidate has a positive effect and joint q ≤ 0.05.
- The same inflation is applied to every contrast and sensitivity run. Gene-subset checks use the
  inflation estimated for their own reduced sets.

**Local divergence.** D(s) = AUC_P(|Z(s)|) − 0.5. The peak position and the median member Z at the
peak describe each pathway.

### SM14 · Specificity: relabeling, NCDR and target–decoy calibration

**Relabeling.** The four specimens have three balanced 2 + 2 splits: the species split and two
relabelings, {`HUK1_COR1` + `Ctrl1A2`} and {`HUK1_COR1` + `Ctrl1A4`} against the rest.
- Every strategy was rerun unchanged on both relabelings, as in Lamian's random-partition null
  (Hou et al. 2023).
- In relabeled gene models, species × spline enters as a nuisance term, and within-group specimen
  contrasts absorb the species offset.
- Relabeled DESeq2 used `~species + group`.
- Two within-species comparisons (mouse vs mouse, human section vs section) were added as controls.

**Motivation for relabeling (not a guarantee of specificity).**
- Write a specimen's deviations as an offset a_ij and a positional deviation b_ij(s), independent
  across specimens.
- The species split and the two relabelings are the three orthogonal ±½ contrasts of the four
  specimens. Under independence they have equal noise variance, so the relabelings are exchangeable
  draws of the specimen-level null.
- The design is a split plot (Altman & Krzywinski 2015): offsets have 2 df of error, and positional
  terms 2(k − 1) df.
- A term shared by both human sections cancels in the relabelings and is not covered.
- The argument does not predict that position-dependent screens are specific. Whether they are
  depends on how large and how coherent within pathways the specimen × position deviations are. This
  was measured, in our cohort and in external many-draw nulls (below).

**Summaries.**
- NCDR is the mean number of relabeled calls divided by the number of species calls.
- A strategy is called specific when NCDR < 0.25 and relabelings call < 5% of tested pathways. The
  second condition was added after the first run (Supplementary Note 8).
- Decoy-estimated FDP (Elias & Gygi 2007) is FDP(t) = (mean relabeled count ≥ t) / (species count
  ≥ t), made monotone as a q-value. We report calls at FDP ≤ 10%, with 300-draw pathway-resampling
  intervals, and the minimum attainable FDP.
- Decoy FDP was the pre-specified primary R2 metric. It was replaced by NCDR because pathway-score
  models reach FDP 0 (Supplementary Note 8).
- Decoy inflation a and implied ρ were estimated from the slope of mean decoy z² on pathway size.
- NCDR is a ratio of call counts, not a false-discovery proportion. A strong true effect dominates a
  competitive ranking, so relabeled counts can overstate how many species calls are specimen noise
  (see the precision check below). NCDR measures vulnerability to specimen noise.

**Many-draw specimen nulls (notebook 45; protocol recorded before it ran).**
- *Pools.* PT nuclei with author segment labels:
  - CELLxGENE Census 2025-11-08 mouse dataset `25818bf7…`, with embryonic and newborn stages
    removed (male and female pools);
  - Census human cortex dataset `09b518f9…`;
  - Lake et al. (2023) KPMP v1.5 healthy cortex.

  Donors needed ≥ 50 nuclei in every segment. Counts over the accepted ortholog panel were normalised
  as for our structures. The extraction script is `analysis/scripts/fetch_census_pt_cells.py`.
- *Fixed once per pool:* tested genes (≥ 2% detection, equal-donor mean), matching strata, notebook
  37's 1,513 pathway sets (10–300 tested members) and a CAMERA variance inflation from within-cell
  residual correlations.
- *Strategies.* For each 2 + 2 partition:
  - unit-level weighted partial F statistics for a group offset and for group × segment steps.
    These use notebook 37's designs with segment intercepts as the shared shape, computed exactly
    from cell means and within-cell sums of squares, and are tested with the covariate-matched
    rank-AUC (exact null moments) and the joint test;
  - whole-PT and per-segment pseudobulk moderated t statistics (voom-style weights, empirical-Bayes
    variances), tested with signed preranked GSEA.

  Calls are BH q ≤ 0.05 within partition and strategy.
- *Same-species splits* are nulls: 300 splits each for male and female mice, 105 for Lake cortex and
  45 for Census human cortex.
- *Cross-species draws* (100): the species split is compared with its two balanced relabelings
  (species segment shape as nuisance) and with an all-donor reference. Precision is the share of
  species calls also called with all donors.
- *Mechanism.* The median within-pathway correlation of member genes was computed for specimen
  offsets, specimen × segment deviations and within-unit residuals, in each pool and in our cohort.
- *AKI positive control.* The mouse AKI design (`mouse_only_v5`: two AKI against two control kidneys,
  on the mouse-only PT DPT) was run with the same strategies plus notebook 37's coordinate
  T_spatial. Detection passes if injury programs are called and change most in S3. Specificity was
  pre-specified as relabeled calls ≤ 5% of condition calls; the main R4 rule (NCDR < 0.25 and < 5% of
  tested pathways per relabeling; 1,429 tested) is reported beside it.
- The atlases have no coordinate, so the external nulls compare average-level screens with
  segment-step position-dependent screens. They do not test the smooth coordinate itself.

### SM15 · Conventional comparators

**Pseudobulk.**
- Raw counts were summed per specimen × reviewed segment; whole-PT pseudobulk is the sum over
  segments.
- Genes needed ≥ 10 counts in ≥ 2 of 4 pseudobulks within each segment, giving 10,933 genes and
  1,508 pathways.
- DESeq2 models (Love et al. 2014) were fitted with PyDESeq2 0.5.4 (Muzellec et al. 2023): design
  `~species`, Wald tests, no Cook's refit.
- Signed Wald statistics were tested by preranked GSEA (Subramanian et al. 2005) in GSEApy 1.3.1
  (Fang et al. 2023): weight 1, multilevel, seed 12, 10–300 genes. BH was applied over pathways for
  whole-PT and over pathway × segment for segments.

**Other screens.**
- Pathway-score models: per-structure mean of weighted-standardised member expression, tested with
  structure-level F.
- Gene-level step models: species × S1/S2/S3, and species × 6 equal coordinate bins.
- T_level and T_total.
- A size-only null and unweighted preranked GSEA on T_spatial.

### SM16 · Building the reporting set

The reporting set is built in five steps:

1. **Candidates.** Joint q ≤ 0.05 and a positive effect.
2. **Not specimen-sensitive.** Not called by the joint test in either relabeling or either
   within-species comparison.
3. **Stable.** Retained (positive effect, joint q ≤ 0.05) in ≥ 5 of 7 planned refits: each specimen
   left out, spline basis 4 and 8, and 5% detection. Strata were recomputed in each refit.
4. **Robust.** Two checks, both required:
   - *Registration.* Each human section is warped piecewise-linearly so that its S1→S2 and S2→S3
     transitions fall on the mouse transitions. Transitions are logistic fits per specimen
     (C = 10⁴). After refitting, the effect ratio must be ≥ 0.7 and joint q ≤ 0.10.
   - *Region.* Depth is the mean distance to the three nearest glomerulus-labelled structures, with
     centroids verified. Glomeruli were taken from notebook 03's labels. Every notebook-03
     glomerulus is also a glomerulus in notebook 13, which adds three. Using notebook 13's glomeruli
     would change the deep mouse S3 exclusion by 22 of 1,581 structures (Jaccard 0.97 in Ctrl1A2;
     identical in Ctrl1A4). Mouse S3
     deeper than the 95th percentile of same-specimen S1 depth is removed. The same number of mouse
     S3 structures is also removed at random 20 times. A pathway is region-sensitive if its effect
     turns non-positive or falls below the 5th percentile of the power-matched ratios.
5. **Core.** Robust, and also passes two gene-subset checks, each with effect ratio ≥ 0.7 and joint
   q ≤ 0.10:
   - genes detected in ≥ 5% of structures in each species;
   - removal of mouse PT sex-biased genes, defined from Xiong et al. (2023) per-segment tables as
     |log2FC| ≥ 1 and padj ≤ 0.05 in any segment (511 male-biased, 700 female-biased).

**Annotations, not filters:**
- removal of flagged ambient, neighbouring-segment and shared-exon genes;
- equal probe counts;
- 5-bin matching;
- an analysis on notebook 03's DPT coordinate (notebook 12 logic v5; 11,106 genes, 1,516
  pathways). This is a historical sensitivity only and is not cited as robustness support.
  Robustness to the coordinate is assessed instead by rerunning notebook 37 on DPT13 and SCF13-ED
  (SM29; Table S2 coordinate columns);
- the specimen-level split-plot test.

### SM17 · Programs and descriptors

**Programs** are built after the list is frozen.
- Contributing members are a pathway's members in the top 10% of T_spatial; a pathway with fewer than
  three uses all its members.
- Pathways are grouped by average linkage on 1 − overlap coefficient of contributing members, cut at
  0.5.
- Each program is named after its lowest-q pathway.

**Descriptors.**
- *Shape* (reversing, localised, graded) of the median member Δ(s).
- *Averaging loss*: the median over contributing members of 1 − |mean_s Δ| / mean_s |Δ|.

**Direction reversal.**
- *Display rule:* ≥ 3 contiguous grid points with |Z| ≥ 3 and |Δ| ≥ 0.25 on each side.
- *Segment DESeq2 confirmation:* opposite significant signs.

### SM18 · Specimen-level split-plot test

**Model.**
- Counts were summed per specimen × 8 equal-width coordinate bins, keeping cells with ≥ 20
  structures.
- log2-CPM values received voom-style weights (Law et al. 2014).
- Each gene was fitted with specimen + cubic spline(position, 3 df) + species × spline.
- The 3-df interaction was tested with a moderated F (Smyth 2004) against specimen × position
  deviations (22 df).
- The test is a Python implementation (`pseudospace/split_plot.py`), checked against simulations.

**Calibration** was judged under relabeling by median-F inflation λ and the fraction of genes at
p ≤ 0.01.

**Pathways.** A matched rank-AUC on the moderated F, with the inflation estimated from pseudobulk
residual correlations.

### SM19 · Zonation amplitude (notebook 39; protocol pre-registered)

**Gradients.**
- Pseudobulk log2 CPM ratios of summed counts per specimen and reviewed segment, computed over
  orthologs measured in all specimens.
- Contrasts: S2 − S1 (primary, cortical in both species); S3c − early, with mouse S3 restricted to
  cortical-like structures; and S3 − early.
- Genes needed ≥ 20 counts in every specimen × segment pseudobulk of both species, and a reference
  mean ≥ 3 log2 in one compared group.

**Reference-axis (projection) ratio.**
- Both of our species gradients were regressed on an independent reference gradient. The ratio is
  cov(human, ref) / cov(mouse, ref).
- This ratio measures how much of each species' gradient lies along the reference axis. When the
  cross-species correlation is low (main text R2), it is a projection ratio and understates the human
  amplitude off that axis. The noise-corrected SD ratio (below) is the amplitude measure used in R3.
- References were male mouse microdissection (Chen et al. 2021), microdissection by sex (Chen et al.
  2023), mouse snRNA by sex (Census) and Lake healthy human cortex.
- Fits used OLS (primary), Deming (error-variance ratio from replicate variances) and orthogonal
  regression.
- A two-stage bootstrap resampled donors, then genes (2,000 draws). The four human-section ×
  mouse-specimen pairings are also reported.

**Measurement checks.** Expression tertiles, the per-structure log1p-mean contrast, and genes with
equal probe counts.

**Tissue state.**
1. Mouse AKI against controls, within the mouse-only pipeline's labels and counts.
2. Lake donor injury scores against donor amplitude. Each marker's whole-PT log2 CPM was z-scored
   across donors, and the score is the mean over markers present in the gene space. Injury markers:
   Havcr1, Lcn2, Vcam1, Cd44, Krt8, Krt18, Krt20, Sox9, Vim, Spp1, Clu, Cdkn1a, Timp1, Lgals3, Ccn2
   (Ctgf). A stress score used Fos, Fosb, Jun, Junb, Egr1, Atf3, Ier2, Hspa1a, Hspa1b and Dnajb1.
   Donor amplitude is the slope on male microdissection.
3. Low- and high-injury halves of our structures.

**Post-hoc diagnostics (not pre-registered).**
- A symmetric external ratio on a human axis.
- The noise-corrected SD ratio of gradients: signal variance is the variance of the mean gradient
  minus the mean per-gene sampling variance.
- The disattenuated cross-species correlation: r divided by √(reliability_human ×
  reliability_mouse).
- Lake donor amplitude on a human axis.

**Flattening classes.** A gene is beyond flattening when its human contrast departs from
r × mouse contrast, with r = 0.152. For pathways, the gene statistic is the largest |human − r·mouse|
over the two contrasts, tested by the joint matched test within the 36 robust pathways.

### SM20 · Zonation classes, v1 setting (notebook 40)

Notebook 40's rules form the v1 alternative setting (Supplementary Note 3; Table S7). The primary
classes are notebook 43's symmetric classes (SM26), which reuse the calls below with amplitude-matched
human thresholds and external flat confirmation.

**Contrasts and calls.**
- Per specimen: S2 − S1, and S3c − early with cortical-like mouse S3. Contrasts are pseudobulk log2
  CPM ratios, and genes need ≥ 20 counts in every specimen × segment pseudobulk of both species.
- Within each species, the specimen-mean contrast was tested with a moderated one-sample t. Variances
  were shrunk within ten expression bins, in the manner of limma's `squeezeVar`. BH was applied across
  genes.

**Confirmation and status.**
- *Zonated:* q ≤ 0.05 and |contrast| ≥ 0.5 log2, and an independent same-species reference agrees in
  direction (one-sided moderated donor-level p ≤ 0.05). References were CELLxGENE Census mouse snRNA
  (12 male donors) and Lake et al. (2023) healthy human cortex (7 donors).
- *Flat:* the 95% interval lies within ±0.3 log2.
- *Indeterminate:* otherwise.

**Classes.**
- *Conserved:* both species zonated in the same direction.
- *Reversal:* zonated in opposite directions.
- *Mouse-only / human-only:* zonated in one species and flat in the other.
- A gene is assigned over the two contrasts combined.

The primary analysis uses 7,407 probe-balanced genes. Probe-imbalanced genes and genes expressed in
one species only are reported separately.

**Alternatives.**
- A1: no confirmation.
- A2: no effect floor.
- A3: human floor scaled by the amplitude ratio.
- A4: flat margin 0.5.
- A5: mouse confirmation by microdissection.

**Enrichment.** For each class, the class indicator was tested with a competitive rank-AUC against
expression- and detection-matched random sets (9,999 draws).
- Variance was inflated by a CAMERA factor estimated from residual correlations after removing
  specimen × segment means.
- BH was applied within each class.
- Robustness was checked with 5-bin matching and with the lowest expression tertile removed.
- Enriched pathways were grouped by member overlap.

**Sex check.** Class membership was cross-tabulated against the Xiong et al. (2023) sex-bias
annotation. This check was added after the classes were computed (Supplementary Note 8).

### SM21 · External segment-resolved datasets and verification of lead genes

**Datasets.**

| Dataset | Content |
|---|---|
| GSE150338 (Chen et al. 2021) | Microdissected male mouse PTS1–3; TPM replicates |
| GSE212213 (Chen et al. 2023) | Microdissected mouse PTS1–3, both sexes; mean TPM |
| GSE56743 (Lee et al. 2015) | Microdissected rat S1–S3; RPKM; excluded by notebook 43's reliability rule (SM26) |
| CELLxGENE Census 2025-11-08 (CZI Cell Science Program et al. 2025), mouse dataset `25818bf7…` (Chen S et al. 2025) | 12 male and 12 female donors with ≥ 50 nuclei per segment |
| Census human renal-cortex dataset `09b518f9…` (Acera-Mateos et al. 2026) | 6 donors; convoluted PT vs S3 |
| Lake et al. (2023) snRNA v1.5 (CELLxGENE `a12ccb9b…`) | Cortex nuclei (`region C`, cortex tissue), healthy donors, author labels PT-S1/S2/S3. 7 donors with ≥ 50 nuclei per segment (4 female, 3 male); 12 at ≥ 20 nuclei as sensitivity. AKI and CKD donors were used for tissue state. |

Datasets used only in the revision analyses are described with those analyses: the many-draw nulls
in SM14 and the state, Visium and rat protein data in SM27.

**Processing.**
- *Census and Lake.* Counts were summed per donor and segment, normalised by all-gene totals and
  log2(CPM + 1); embryonic donors were excluded. Lake genes were matched by Ensembl id.
- *Microdissection.* log2(TPM + 1).
- *Rat.* Symbols were matched case-insensitively (approximate).

**Agreement with our data** (notebooks 34 and 37).
- Gene-level Spearman of our species × segment interaction against the external interaction, on
  informative genes (detected in ≥ 10% of our structures in each species; external log2 CPM ≥ 1).
  Comparators were expression-matched random genes and top genes by expression.
- Pathway-level: a Welch t per gene (human against male mouse donors), signed by our direction, then
  a matched rank-AUC per pathway with donor-level inflation and BH within the robust list.
- The global component is the median contrast and the slope of the interaction on the mouse
  gradient. It was removed per donor and per specimen.
- Comparators: uncalled and relabel-called pathways; female mice; additional removal of the slope
  component.

**Lead-gene verification** (notebook 38; rules fixed before the external data were opened).
- Each claimed gene was tested for:
  - our arms in each specimen (or all four pairings);
  - robustness to attenuation;
  - the Lake human arm;
  - interaction sign against male and female mouse snRNA;
  - mouse arms in snRNA and microdissection of both sexes;
  - for one-species genes, absence (mouse microdissection maximum < 5 TPM).
- Probes per gene came from the public 10x mouse v2.0 panel and from our own human
  `molecule_info.h5`.
- Human Protein Atlas kidney IHC records presence only.
- Verdicts were headline, supporting or not replicated.
- Attenuation *a* was estimated on 70 control genes chosen only from external data.

### SM22 · Literature check

Each confident pathway (with up to eight driver genes), each direction-reversing gene, and each
pathway called only by conventional screens was searched.

- *Search.* Titles and abstracts in Europe PMC, using only gene, pathway and generic terms. Key
  papers were read in open-access full text where available.
- *Status.*
  - KNOWN: species difference and PT position reported.
  - PARTLY KNOWN: one part reported.
  - NOT FOUND: nothing after ≥ 2 targeted queries. This is a candidate, not proof of novelty.
  - CONTRADICTS: the literature reports the opposite.
- *Who searched.* The searches were run by AI agents (Claude, Anthropic) through the Europe PMC
  interface, with queries containing only gene, pathway, method and generic terms.
- *Verification.* Each cited source was then checked by PubMed identifier, and its abstract or open
  full text was read.
- *Records.* The literature-check table in the repository records the outcome and number of queries
  per claim. The query strings themselves are not part of the repository.

### SM23 · Statistics and multiple testing

**What the statistics mean.**
- Two mice and two sections from one human donor give no estimate of variation among human donors.
  Results are conditional on these specimens and descriptive.
- Structure-level statistics rank and describe. Pathway tests are competitive tests of a gene
  statistic against matched genes, not population tests of a species effect.

**Multiple testing.** BH at q ≤ 0.05 within each family:
- each statistic across all pathways;
- each partition × statistic for specificity;
- within each zonation class;
- within the robust list for external replication;
- genes within each species for zonation calls.

**Display rules and limits.**
- With 9,999 draws, the smallest empirical p is 10⁻⁴.
- Themes, peaks, reversal flags, programs and external concordance are display rules, not tests.
- Two relabelings screen only for gross nonspecificity. Notebook 45's many-draw nulls place them in
  an external distribution (SM14), but our cohort's own NCDR still rests on two draws.
- NCDR is not a false-discovery proportion (SM14).
- Noise-corrected correlations and the conservation index have gene-bootstrap intervals (1,000
  draws). They do not resample donors of our cohort, which has two per species.
- The symmetric classes, cross-fitted label checks and coordinate reruns are sensitivity analyses
  of fixed rules. None selects a model or a threshold.

**Coordinate limits.**
- Human S1 ordering depends on read depth: equal-depth agreement is 0.36–0.46, and the
  library-size correlation falls from −0.31 and −0.38 to −0.04 and −0.13 after thinning.
- Human placement moves by up to 0.29 when a third of the genes or half the reads are removed.

### SM24 · Software and reproducibility

**Analysis environment.**
- Python 3.11.16.
- Scanpy 1.11.5, AnnData 0.12.19, NumPy 1.26.4, SciPy 1.17.1, pandas 2.3.3, statsmodels 0.15.0,
  patsy 1.0.3, scikit-learn 1.9.0, python-igraph 1.0.0, leidenalg 0.12.0.
- scFates 1.2.5, PyDESeq2 0.5.4, GSEApy 1.3.1, shapely 2.1.2, rpy2 3.6.7.
- R 4.5.3 with harmony 2.0.5.
- `cellxgene_census` 1.18.0 for the Census extraction (committed as
  `analysis/scripts/fetch_census_pt_segments.py` for segment pseudobulks and
  `analysis/scripts/fetch_census_pt_cells.py` for the per-nucleus counts of notebook 45).

**Segmentation environment.** Python 3.10, PyTorch 2.5.1 (CUDA 12.1) and scikit-image (van der Walt
et al. 2014) [VERIFY the versions used in the inference run].

**Numerical reproducibility of the coordinate.** The scFates curve is numerically sensitive.
- Multiplying by 10⁴/library instead of dividing changes values by ≤ 9 × 10⁻¹⁶, yet lowers median
  within-segment agreement to 0.985 and shifts individual structures by up to 0.44.
- Unsorted CSR indices alone give an overall Spearman of 0.989 and a maximum shift of 0.48.
- Exact reproduction requires scanpy's normalisation and sorted indices. Notebooks 35–36 verify this
  before each use. DPT is unaffected.

**Workflow.**
- Expensive stages are cached by parameters, inputs and code.
- Every result comes from a committed, output-free notebook taking `--data-root` and
  `--results-root`.
- Notebooks 35–46 are committed.

### SM25 · Segment-label validation (notebook 42)

Notebook 42 runs in `results/pt_revision_labels/` (logic `42.label_validation.1`). Its protocol was
recorded before any result. Mouse uses the fine GeoJSONs and human the `v2` GeoJSONs, with
`feature_index` as a positional index.

**Geometry and pixel scale.**
- *Centroid check.* Every polygon centroid had to lie within 2 px of the centroid recorded in the
  matrix. The observed maximum offset was 0 px in all four specimens.
- *Pixel size.* It was read from each feature's `Area px` and `Area µm²` (median ratio per
  specimen): 0.441 µm (mouse) and 0.500 µm (human).
  - The segmentation configuration works at 0.44068 µm per pixel. The human µm fields are therefore
    probably mis-scaled: (0.44068/0.5)² × 0.974 = 0.757, matching the observed human bin coverage of
    0.765.
  - Human µm distances are converted at 0.44068 µm per pixel, i.e. multiplied by 0.881 (−11.9%);
    the uncorrected values are 13.5% too large [VERIFY against scanner metadata].
  - The bin-to-polygon join uses pixel coordinates, and the structure filters use counts, so neither
    is affected. Human areas in µm² (overstated by 29%) and notebook 11's extent analysis are
    affected.

**Glomeruli.**
- *Primary:* the morphological glomerulus polygons of the H&E segmentation. These are the human
  `v2` files, and for mouse the retired `v4` files, used for geometry only. The v4 files share the
  fine files' pixel frame: the median offset to matched polygons is 0.2–0.3 px.
- *Sensitivity:* notebook 13's expression-defined glomerulus cluster.
- *Size rule (pre-registered).* If glomerulus polygons were not glomerulus-sized (median area ≥ 3×
  the PT median), distances were to be reported descriptively only.
  - Human glomeruli are 5.6–6.1× the PT median.
  - Mouse glomeruli are 1.3–1.4×, so mouse results are descriptive.

**Label checks.**
- *(a) Distance.* Distance from each PT polygon centroid to the nearest glomerulus boundary in the
  same specimen, with the Spearman correlation of label order (S1 = 0, S2 = 1, S3 = 2) with distance.
- *(b) Attachment.* PT polygons with a boundary within 2 µm of a glomerulus (primary) or within
  10 µm (sensitivity), and the odds ratio for S1.
- *(c) Reference transfer.* Expression only, without the joint embedding.
  - *Genes:* measured in both species, detected in ≥ 10% of that species' PT structures, maximum
    reference log2 ≥ 1, and |difference| ≥ 1 log2 between some pair of reference segments.
  - *Primary score:* per structure, the Spearman correlation between log-normalised expression and
    each reference segment profile; the label is the argmax.
  - *Sensitivity score:* Pearson correlation between the specimen-centred profile and the
    reference profile centred over its three segments.
  - *References:* Lake/KPMP cortex S1/S2/S3 (human); microdissection (primary) and male Census
    snRNA (mouse).
  - *Pre-specified reading:* "supports the current human labels" requires human agreement ≥ 0.70
    and within 0.10 of mouse agreement.

**Integration genes.** The 611 PCA genes and 752 Harmony HVGs were cross-tabulated against the
notebook 40 classes.

**Spillover index.** Counts of fixed non-PT markers, divided by total counts per structure. Markers
were kept if detected in ≥ 1% of non-PT structures:
- ascending limb: Umod, Slc12a1;
- DCT: Slc12a3, Pvalb;
- CNT/CD: Aqp2, Calb1, Scnn1g;
- podocyte: Nphs1, Nphs2, Podxl;
- endothelium: Pecam1, Emcn, Kdr;
- immune: Ptprc;
- stroma: Pdgfrb, Col1a1, Col3a1;
- smooth muscle: Acta2, Myh11.

**Registration QC.** For each specimen, the analysis counted polygons with zero bins and the share of
these inside the convex hull of bin-covered polygons. Bin coverage was computed as n_spots × 4 µm² ÷
polygon area.

### SM26 · Same-species ceilings, conservation index and symmetric classes (notebook 43)

Notebook 43 runs in `results/pt_revision_classes/<labels>/` (logic `43.revision_classes.1`). Its
protocol was pre-registered in the notebook header. Labels are a parameter: the reviewed labels are
primary, and transferred labels are sensitivities (SM28).

**Gradients.**
- Contrasts: S2 − S1 and S3 − early, where early = S1 + S2. Our mouse data also get the
  region-matched S3c − early.
- Datasets:
  - human: ours (2 sections), Lake/KPMP healthy cortex (7 donors) and Census human cortex (6 donors;
    S3 − early only);
  - mouse: ours (2 mice), Census snRNA by sex (12 donors each), and microdissection (male GSE150338;
    sex-specific GSE212213, means only);
  - rat: microdissection (GSE56743).

**Correlation.**
- Pearson correlation of gene gradients, each winsorised at the 1st and 99th percentiles.
- Genes had to be expressed in both datasets: ≥ 20 counts in every relevant pseudobulk for our data,
  and mean log2 ≥ 3 in one segment group externally.
- Correlations were disattenuated by each dataset's replicate reliability (signal ÷ total variance
  across genes). Datasets without replicates are reported raw.

**Conservation index.**
- Index = r*(our human, our mouse) ÷ √(r*(our human, Lake) × r*(our mouse, mouse snRNA male)). It
  is computed on genes expressed in all four datasets, with 95% intervals from 1,000 gene bootstraps.
- The external-only version replaces the numerator with r*(Lake, mouse snRNA male).
- Decision: an upper bound < 0.5 in both primary contrasts keeps "largely different genes"; 0.5–0.8
  means "partly conserved"; ≥ 0.8 means "largely conserved".

**Rat rule.** Rat was used only with replicate reliability ≥ 0.5 and ≥ 70% symbol matching. It
failed: reliability was 0.03 (S2 − S1) and 0.45 (S3 − early), and 80% of genes matched.

**Symmetric classes.** Notebook 40's moderated per-species calls and reference confirmation of
zonated calls (SM20), plus two changes:
- *A3b:* the human effect floor and flat margin are multiplied by the contrast's noise-corrected
  human ÷ mouse amplitude (0.58, 0.84 and 0.70 for S2 − S1, S3c − early and S3 − early).
- *A8:* a species-only class also needs the other species' atlas to be flat-compatible: one-sided
  p > 0.05 in the same direction and |mean| < 0.5. Otherwise the gene is indeterminate.

The symmetric setting (A3b + A8) is primary. Female confirmations (A6, A7) are reported. Counts are
given for probe-balanced genes per contrast and for the union of S2 − S1 and S3c − early.

Asymmetry was judged by the mouse-only ÷ human-only ratio: ≥ 2 in both S2 − S1 and S3c − early means
"asymmetry supported"; in neither means "no consistent asymmetry".

**Pathway categories.** Each robust pathway's contributing members (top 10% by notebook 12's
T_spatial) were classed per gene, over S2 − S1 and S3c − early:
- *mouse-led:* zonated in mouse in some contrast where human is not zonated, and never zonated in
  both;
- *human-led:* the reverse;
- *shared:* conserved or reversal in some contrast.

Among members zonated in either species (at least 3 needed), the pathway takes the category of ≥ 50%
of members ("summary of mouse zonation", "shared" or "human-led"); otherwise it is "mixed". Pathways
with fewer than 3 zonated members are "unresolved". Pathways failing notebook 37's equal-probe check
are flagged as probe-sensitive.

**Transporter shares.** Shares of zonated Slc, Abc and Aqp genes per class were counted from
`gene_classes_symmetric.csv` (union class, probe-balanced genes). This count is descriptive and was
reproduced by notebook 47 (`transporter_class_shares.csv`, `transporter_direction_agreement.csv`),
which also reports direction agreement among genes zonated in both species.

### SM27 · Physiological state, donor tissue state, additional human donors and rodent protein (notebook 44)

Notebook 44 runs in `results/pt_revision_state/` (logic `44.revision_state.1`). Its protocol was
recorded before it ran. Every input is public except the tissue-state placement, which uses our
structures.

**Physiological state.**
- *Programs* (fixed in advance): sterol/SREBP2 (12 genes), PPARα/fatty-acid oxidation (13),
  glutathione synthesis (4), tyrosine catabolism (6), the lead genes (Gatm, Gamt, Acaa2, Dcxr, Ugt3a1,
  Acox2), and structural markers (Slc5a2, Slc5a12, Slc7a13, Slc22a7, Slc34a1, Lrp2).
- *Whole-kidney state effects:*
  - 24-h fasting against fed, GSE267280 (one pooled column per condition, so a point estimate)
    (Jiang et al. 2025);
  - ZT18 against ZT6, GSE277302 (Nguyen et al. 2025).

  Each was compared with the gene's external species × position interaction (Lake against male
  mouse snRNA). A program's state effect was called comparable if its median largest state effect was
  ≥ 0.5 × its median |interaction|. GSE54650 was dropped because its array annotation could not be
  retrieved.
- *Segment zonation at two times of day:* control mouse snRNA at ZT4 and ZT16 (GSE331332; 2 mice
  each; Wigger et al. 2026).
  - PT nuclei were called by a canonical PT score above the Otsu threshold.
  - S1 and S3 were called from Census-derived marker sets that exclude every tested gene.
  - Pseudobulk S3 − S1 was computed per mouse.
  - A program was state-stable if ≥ 80% of its zonated genes kept their sign in all four mice and
    the amplitude ratio ZT16 ÷ ZT4 lay in 0.67–1.5.

**Donor tissue state.**
- *Signature.* 30 up and 30 down genes, from the aPT-against-healthy PT log2 CPM differences of Lake
  cortex donors (≥ 30 aPT and ≥ 100 healthy PT nuclei).
- *Score.* Mean within-sample percentile of the up genes minus the down genes.
- *Validation.* Against the altered-state fraction across Lake donors.
- *Calibration.* A platform offset was estimated in mouse: our control mice against healthy Census
  donors. As a check, our AKI mice must score above our controls. The offset is assumed to transfer to
  human.

**Additional human donors with spatial data.**
- *Data.* Abedini et al. (2024) Visium (GSE211785), using samples with ≥ 50 PT_S1 and ≥ 50 PT_S3 spots.
- *Pre-specified label check.* SLC5A2 falls and SLC5A1, SLC7A13 and SLC22A7 rise in ≥ 70% of samples.
- *Lead genes.* The human arm replicates if its sign holds in ≥ 70% of informative samples (at least
  5 samples with |S3 − S1| ≥ 0.1).

**Rodent protein.** Rat Kidney Tubule Expression Atlas, microdissected S1/S2/S3 (Limbutara et al.
2020; 3 samples per segment), in copies per cell. Log2 ratios are centred on the median ratio of
proteins quantified in all nine PT samples. We quote S3/S2 because rat S1 samples carry less
protein.

### SM28 · Cross-fitted label sensitivity (notebook 46)

Notebook 46 runs in `results/pt_revision_classes/crossfit/`. Its protocol was recorded before any
result.

**Label settings.**
- (i) Human labels transferred from Lake (centred score), mouse reviewed.
- (ii) Human from Lake, mouse transferred from microdissection.
- (iii) Reviewed labels in both species, as the baseline for halving the genes.

**Cross-fitting.** Genes were split into two halves by sha256(gene) mod 2. Labels were transferred
using the reference-zonated genes of one half and evaluated on the other half, then the halves were
swapped. Notebook 46 re-implements notebook 43's primary conservation index and symmetric class
counts from the same modules. On the reviewed labels with all genes it reproduces notebook 43
exactly.

**Decision rule.**
- The claim stands if the primary index is ≤ 0.5 with an upper bound < 0.6 in every setting, half and
  primary contrast.
- Species-only counts are "unskewed" if mouse-only ÷ human-only < 2 in S2 − S1 and S3c − early.

Notebook 43 was also run unchanged with full-gene transferred labels. Those runs are circular and are
reported for completeness only.

### SM29 · Coordinate robustness of pathway calls

Notebook 37 was run unchanged through nbclient on two further coordinates, DPT13 and SCF13-ED.
Outputs are in `results/pt_pathway_final/<label>/`, and the summary is `d4_coordinate_runs.csv`.

For the union of the robust pathways under scFates and DPT13 (40 pathways),
`coordinate_robustness.csv` records:
- robust and core status under each coordinate;
- the joint z and peak position under each;
- sign agreement of the direction at the peak, wherever the pathway is a candidate;
- the notebook 43 category.

Notebook 47 joins these columns to every tested pathway (`table_S2_pathways_with_coordinate.csv`),
which is Table S2.

### SM30 · Figure assembly (notebook 41)

Notebook 41 draws every main and supplementary figure from saved results. It fits no model and runs no
test.
- Descriptive aggregates that no notebook saves (segment pseudobulks, specimen bin means, reference
  gradients) are recomputed with the source notebooks' functions. Each must reproduce saved
  statistics before it is drawn.
- Key counts are asserted against their tables.
- A layout check reports overlapping text.

Outputs go to `results/paper_figures/`, with a panel-source manifest.

### SM31 · Revision addenda (notebooks 47 and 48)

**Label and class addenda (notebook 47).** Notebook 47 runs in `results/pt_revision_addenda/` (logic
`47.addenda.1`). Its protocol was recorded before any result.
- *Index without integration genes.* Notebook 46's evaluation function, unchanged, with the 611
  pass-2 PCA genes removed from the evaluation genes. The moderation prior still uses all genes. The
  settings are reviewed labels and notebook 46's cross-fit label files (i) and (ii). The claim stands
  if the index is ≤ 0.5 with an upper bound < 0.6 for S2 − S1 and S3 − early in every evaluation.
- *Transporters.* `^(Slc|Abc|Aqp)` on mouse symbols, notebook 43's union classes, probe-balanced
  genes. It reports shares of zonated genes, the rate of zonation in both species, and direction
  agreement among genes zonated in both (conserved ÷ (conserved + reversal)).
- *Integration genes by class.* The 611 PCA genes and 752 Harmony HVGs counted per symmetric class.
- *Robust-list decomposition.* Each robust pathway is checked against the beyond-steps and step-screen
  joint tests of the same run (q_joint ≤ 0.05, species partition).
- *Peaks in segment units.* Each run's mean segment transitions are used, and the value is the segment
  index plus the fraction through the segment.
- *Reversals.* Reversal counts per contrast under reviewed, transferred, cross-fitted and
  no-integration-gene settings.
- *Table S2.* Notebook 37's table joined with robust, core and peak under DPT13 and SCF13-ED.

**Pathway addenda (notebook 48).** Notebook 48 runs in `results/pt_revision_addenda_pathways/` (logic
`48.pathway_addenda.1`). Its protocol was recorded before any statistic was computed.
- *Overlap-preserving replication.* The score is notebook 37's per-pathway donor-level joint z,
  recomputed with exact matched-null moments (median |Δ| 0.009 from the saved z). A null draw relabels
  genes within the matching strata and rescores every set at once, keeping set sizes, strata and
  pairwise overlaps. Each draw recomputes the matched z and the variance inflation (9,999 draws per
  reference and statistic). Collections are the 36 robust, 26 core, 21 and 28 relabel-called pathways,
  and the 17 robust and 13 relabel-called programs. The decision rule required robust > random and
  robust > relabel-called at both the pathway and the program level.
- *Symmetric-class enrichment.* Notebook 40's competitive enrichment code was run on notebook 43's
  symmetric classes. A guard reproduced notebook 40's saved table from its own classes. Naming rule:
  q ≤ 0.10, not expression-driven, and ≥ 3 class genes.
- *AKI under both rules,* from notebook 45's saved counts.
- *Step-screen reconciliation.* Notebook 45's coordinate-free step-screen calls for our cohort were
  mapped to pathway identifiers and compared with notebook 37's spline-based step screen and robust
  lists.

### SM32 · Analysis-to-notebook map (Table SM1)

| Analysis | Notebook | Results folder |
|---|---|---|
| Tubule-by-gene matrices | 01 | `tubule_by_gene/` |
| Cross-species integration (pass 1) | 03 | `human_vs_healthy_mouse/` |
| Integration pass 2, PT coordinate | 13 | `minimal_pt_scfates/` |
| Pathway method selection, relabeling | 31 | `pt_pathway_method_selection/` |
| Coordinate validation and selection | 35, 36 | `pt_reconstruction_v2/` |
| Final pathway analysis (per coordinate) | 37 | `pt_pathway_final/<coordinate>/` |
| Lead-gene verification | 38 | `pt_literature_deep/` |
| Zonation amplitude | 39 | `pt_zonation_amplitude/` |
| Zonation classes, v1 setting | 40 | `pt_zonation_classes/` |
| Figures | 41 | `paper_figures/` |
| Segment-label validation, coordinate robustness | 42 | `pt_revision_labels/` |
| Ceilings, index, symmetric classes | 43 | `pt_revision_classes/<labels>/` |
| Physiological and tissue state, Visium donors, rat protein | 44 | `pt_revision_state/` |
| Many-draw nulls, AKI control | 45 | `pt_revision_method/` |
| Cross-fitted labels | 46 | `pt_revision_classes/crossfit/` |
| Revision addenda | 47 | `pt_revision_addenda/` |
| Agreement by gradient strength, index sensitivity | 49 | `pt_conservation_strength/` |
| Primary 22-pathway list, its replication, coordinate-free calls | 50 | `pt_primary_pathways/` |
| Overlap-preserving replication, symmetric-class enrichment, AKI under both rules, step-screen reconciliation | 48 | `pt_revision_addenda_pathways/` |

### SM33 · Data and code availability

[VERIFY: deposition of the Visium HD data and segmentation polygons (controlled access for human
tissue), code archive DOI and licence.] The public datasets are listed in SM21.

### SM34 · Agreement by gradient strength and index sensitivity (notebook 49)

Notebook 49 runs in `results/pt_conservation_strength/` (logic `49.conservation_strength.1`). Its
protocol was written before any result. It rebuilds notebook 43's inputs with the reviewed labels and
first reproduces notebook 43's published conservation index exactly (all eight rows, to 10⁻⁹).

**Strength bins without selection bias.**
- *Why cross-fit.* Binning genes on one dataset's |gradient| and correlating that same dataset selects
  its noise extremes (regression to the mean).
- *Folds.* The healthy cortex-atlas donors (7) and the male mouse snRNA donors (12) were each split into
  two folds: sorted by sha256 of the donor identifier, then dealt alternately (4/3 and 6/6).
- *Strength.* In fold A, each species' |gradient| was converted to a percentile rank, and a gene's
  strength was the larger of its human and mouse ranks.
- *Bins.* Genes were binned at the 50th, 80th, 90th, 95th and 99th percentiles of strength.
- *Evaluation.* All correlations used fold B and our data. The folds were then swapped and the two
  directions averaged.
- *Secondary schemes.* Bins on the human rank only, the mouse rank only, the mean rank, and |our mouse
  gradient| (the same data) are also reported.

**Per-bin index.**
- Notebook 43's conservation index was recomputed within each bin, with gradients winsorised at the 1st
  and 99th percentiles within the bin and reliabilities computed on the same winsorised values.
- A noise-corrected correlation was not estimable if a dataset's reliability in the bin was below 0.2.
- Intervals came from 1,000 gene bootstraps and 1,000 module-block bootstraps. A bound was not
  reported when more than 10% of draws were not estimable.

**Decision rules (fixed in advance).**
- (i) *Concentration.* In S2 − S1 and S3 − early, the index in the top 5% of strength exceeds that in
  the bottom 80% by ≥ 0.25, with an interval excluding 0, under both bootstraps. If the bottom index is
  not estimable, the rule applies to the noise-corrected cross-species correlation.
- (ii) *Reproducible but not shared.* In the bottom 80%, both ceilings have lower bounds ≥ 0.3 and the
  index an upper bound < 0.5. If the bottom-80% ceilings are below 0.3 or not estimable, the weaker
  genes are "not reliably zonated".

**Module blocks.** Co-expression modules of our PT structures:
- log-normalised expression of the index genes, centred per specimen and scaled per gene within species
  (species weighted equally);
- embedded by a 50-component truncated SVD;
- clustered by k-means (k = 100; 50 and 200 as sensitivities).

A draw resamples modules with replacement.

**Index sensitivity.**
- Winsorisation q ∈ {0, 0.01, 0.05}, with reliability either unwinsorised (notebook 43) or computed on
  the winsorised values.
- Gene, module-block, external-donor (cortex-atlas and male mouse donors resampled; our specimens fixed)
  and combined donor + gene bootstraps.
- Leave-one-atlas-donor-out human ceilings.

**Direction agreement.** Among probe-balanced genes zonated in both species under notebook 43's
symmetric classes, the share conserved per strength bin (Wilson intervals), with strength from the full
references.

**Post hoc diagnostic** (Supplementary Note 8, deviation 17). The S2 − S1 human ceiling was not
estimable below the 80th percentile: atlas folds of 3–4 donors had reliability < 0.2. Strength was
therefore also defined from the mouse snRNA fold alone, so that the full atlas could give the human
ceiling without selection bias.

**Limitation of the donor bootstrap.** Resampling 7 atlas donors with replacement duplicates donors,
which understates within-reference noise and overstates reliability. Donor-bootstrap percentile
intervals are therefore shifted relative to the point estimate, and the leave-one-donor-out range is
reported alongside.

### SM35 · Primary pathway list and its replication (notebook 50)

Notebook 50 runs in `results/pt_primary_pathways/` (logic `50.primary_pathways.1`).
- *Primary list.* It saves the 22 pathways robust under both the principal curve and DPT13
  (`primary_pathway_list.csv`), with:
  - category;
  - segment-unit peaks under three coordinates;
  - probe flags;
  - the call of notebook 45's coordinate-free step screen on our cohort, recomputed per pathway
    (`cohort_step_screen_per_pathway.csv`).
- *Replication.* Notebook 48's overlap-preserving collection test was repeated on the 22, on their 12
  core pathways and on their 13 programs, against the 21 relabel-called pathways (9,999 draws; male
  mouse references; matched, selection-matched and slope-free scores).
- *Conventional-only pathways.* The group comparison of the 184 pathways called only by conventional
  screens was repeated against overlap-preserving collections.
