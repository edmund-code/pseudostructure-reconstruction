# R1 text: the PT coordinate (workstream 1, drafts for the coordinator)

Sources: notebooks 35 (`results/pt_reconstruction_v2/nb35/`) and 36 (`.../nb36/`), the DPT diagnostic
(`.../dpt_diagnostic/`), and `SCRATCH/ws1/plan.md`, which holds the pre-registered rules,
timestamps and the coordinator's decision.

---

## Results: Reconstructing the PT axis (draft)

We ordered the 12,866 PT structures along a nonbranching principal curve fitted to the
cross-species integrated embedding (Methods). We asked whether it recovers measured zonation,
whether it absorbs the species differences we test, and what it adds beyond S1/S2/S3 labels.

**Recovered zonation.** Each gene was scored only on a coordinate rebuilt without it, using three
gene folds. Its within-specimen rank correlation with the coordinate agreed across genes with
external S3-versus-early contrasts:
- mouse: 0.75 against male snRNA and 0.53 against microdissection;
- human: 0.25 against cortex snRNA, which resolves only convoluted PT versus S3.

Segment labels give similar values (0.73, 0.51, 0.27). Without fold protection, human agreement is
0.31, so part of it comes from the genes that built the coordinate.

**No absorption of species differences.**
- A human-only gradient injected into a tenth of the construction genes survived refitting
  (absorption ratio 1.01).
- Genes whose zonation differs between species externally kept that direction in 91% of 244 cases.
- When the coordinate was built and tested on independent halves of the reads, 59 of 66 robust and
  39 of 42 core pathways were still called. Circular use inflated gene statistics by about 30%
  without changing which pathways were called.

**Fine order within segments is weakly supported.**
- Within reviewed segments, slopes along the coordinate agreed only weakly with adjacent external
  contrasts: 0.17 in mouse, 0.10 in human, against about 0 for ordering by read depth.
- Coordinates from disjoint gene thirds ranked structures within segments only moderately alike
  (0.34–0.70).
- The coordinate predicted held-out genes in the other specimen no better than segment means.
- Thinning to equal depth exposed a depth dependence in human S1: within-segment agreement was
  0.36–0.46, and the correlation of position with library size fell from about −0.35 to about −0.1.
- Removing a third of the genes or half the reads moved human segments by up to 0.29 relative to
  mouse.
- scFates is numerically sensitive to 10⁻¹⁶ input changes.

**Alternatives.**
- DPT was more stable but made the pathway screen non-specific (NCDR 0.40). Its relabeled calls
  were metabolic programs, and it found 36 positional pathway differences between the two mice,
  against 16 for scFates.
- A coordinate built from 160 externally conserved zonation genes failed to register the species.
- Under these alternatives, 56 and 42 of the 66 robust pathways were retained.

**What the coordinate adds** is not resolution finer than segments.
- It replicates measurements across positions within each specimen, so species-by-position
  differences are judged against specimen-by-position deviations (the split-plot argument, R2).
- It localises differences continuously: 19 human-high pathways peak early and 47 mouse-high
  pathways peak late. Unequal positional noise between species explains at most 9–17% of this
  split.

---

## Methods: PT coordinate and its validation (draft paragraph)

**Construction.**
- The PT coordinate is notebook 13's nonbranching scFates 1.2.5 curve: 30 nodes, seed 0, fitted on
  the first five pass-2 Harmony dimensions of the PT structures.
- The pass-2 PCA used 611 genes. These are the 752 species-intersected HVGs that were also among
  pass 1's 2,000 Seurat HVGs; scanpy's PCA restricts to the inherited `highly_variable` column.
- It was rooted at the tip with the higher marker score and scaled to [0, 1].

**Validation** (notebook 35; rules fixed before outputs):
- Gene-wise Spearman correlations with the coordinate, within specimen and within reviewed segment,
  compared with external S3-versus-early or adjacent-segment contrasts. Each gene is scored on a
  coordinate refitted without its sha256 gene fold.
- External data were split into anchor-selection sets (mouse microdissection, female mouse snRNA,
  human donor half A) and evaluation sets (male mouse snRNA, human donor half B).
- A human-only multiplicative injection, exp(±log 2·(2s − 1)), into 10% of the construction genes.
- Binomial thinning of every structure to the pooled PT 20th-percentile library size.
- Poisson count splitting (ε = 0.5) between building the coordinate and testing.
- scFates `n_map = 20` resampling.

**Arm comparison** (notebook 36):
- Alternatives:
  - equal-depth refit;
  - DPT;
  - PT-only integration;
  - a conserved-anchor coordinate (160 genes with the same external zonation direction in all
    selection datasets).
- These were compared under pre-specified rules. Every candidate also had to keep the `T_spatial`
  screen specific (NCDR < 0.25 and fewer than 5% of pathways called under relabeling). This gate
  was applied to all candidates after DPT's result; it is stated as such.
- Registration stability and seed stability were added after notebook 35's outputs.
- Refits whose curve could not be rooted were reported as failed, not retuned.

---

## Figure plan for R1

**Figure 1 (main).**

| Panel | Content | Source |
|---|---|---|
| A | Study design, segmentation and the PT subset | Existing (notebooks 01, 13) |
| B | The pass-2 PT embedding with the scFates curve, coloured by reviewed segment and split by species; specimen segment medians on the coordinate | Notebook 13 outputs; `nb35/species_registration_across_refits.csv` (saved row) |
| C | P1 coordinate versus external zonation per species (fold-protected), with segment labels as reference; P1w within-segment bars with read-depth ordering as a null-like reference | `nb35/figures/p1_p1w_external_concordance.*`, `nb35/p1_*.csv` |
| D | Identifiability and circularity: absorption ratio (1.01); positive-control sign agreement (0.91); count-split arms (robust 59/58/57, core 39/40/39) | `nb35/p2_absorption.csv`, `positive_control_recovery.csv`, `g6_count_split.csv` |

**Figure S1 (supplementary, limitations).**

| Panel | Content | Source |
|---|---|---|
| A | Equal-depth agreement by specimen × segment, plus the library-size correlation before and after thinning (human S1 highlighted) | `nb35/p5_exposure_invariance.csv` |
| B | Species registration gap across refits (gene folds, half reads), SCF13 versus DPT13 | `nb35/species_registration_across_refits.csv` |
| C | Gene-fold within-segment agreement and the numerical floor | `nb35/p3a_*.csv`, `numerical_floor.csv` |
| D | Held-out gene prediction gains versus segment and coverage oracles (all near 0) | `nb36/p4_heldout_prediction.csv` |
| E | Attenuation decomposition (coordinate / labels / external) and the noise-matching dose-response | `nb35/figures/a3_attenuation.*`, `a3b_*.csv` |

**Table S1.** Arm comparison and the decision rules: every metric for SCF13, SCF13-ED, DPT13,
SCF-PT and CAC; C6; the selection rule with the gate; downstream calls, NCDR and retention.
Source: `nb36/arm_metrics.csv`, `c6_decision.csv`, `selection_rule.csv`, `downstream_summary.csv`.

**Figure S2 (optional).** DPT specificity diagnostic: relabeled-call categories; within-species
calls (DPT: mouse 36, human 5; scFates: 16, 0). Source: `dpt_diagnostic/*.csv`.

---

## DPT specificity diagnostic (task 1; descriptive)

- **(a) What the relabeled calls are.** DPT13's two relabelings called 41 and 36 pathways
  (scFates: 17 and 13). 72 of DPT's 77 relabeled calls are metabolic, transport or signalling
  programs, for example fatty-acid, amino-acid and steroid metabolism and SLC transport. Only 5 are
  ribosome, OXPHOS or housekeeping. They are zonation biology, not technical programs. 38 of DPT's
  96 species calls are also called under a relabeling (scFates: 20 of 116).
- **(b) The two mice disagree on DPT.** Within-species models gave 36 mouse-versus-mouse positional
  pathway calls on DPT, against 16 on scFates. Human-versus-human gave 5 against 0. 26 of DPT's 57
  relabeled calls are also within-species calls. Each relabeled group holds one mouse, so a
  specimen-specific positional shape in mouse appears as relabeled signal.
- **(c) It is not a gross shift.** Within segments, the specimens' median positions differ by at
  most 0.02 on both coordinates, and the KS statistics are similar (DPT 0.03–0.14; scFates
  0.05–0.13). However, DPT compresses mouse S3 (pooled IQR 0.055, against 0.10 for scFates) and
  stretches S1 (0.18 against 0.12). Small specimen differences in late-PT composition, such as
  how much outer stripe each mouse section samples, therefore become steep shape differences on
  DPT.
- This interpretation is consistent with the data but was not tested directly. The diagnostic is
  descriptive and was not used for selection.
