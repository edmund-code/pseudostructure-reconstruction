# Choosing the PT pathway analysis for the paper

Notebook `analysis/notebooks/31_pt_pathway_method_selection.ipynb`. Run 2026-10-04,
logic version `31.method_selection.1`. Inputs: notebook 13's PT scFates coordinate and
notebook 12's gene universe, pathway families and sensitivity runs. Nothing upstream is refit.

Cohort: 2 control mice, 2 cortex sections from 1 human donor, 12,272 PT structures,
11,071 genes, 1,513 pathways. All four specimens are male (Y-linked genes detected in every
one), so sex does not separate the species. Everything below is descriptive.

## Final analysis (notebook 37): read this first

`analysis/notebooks/37_pt_pathway_final.ipynb` (logic version `37.pathway_final.1`) recomputes the
whole pathway analysis from one coordinate in a single run. By default the coordinate is notebook
13's scFates; `--coordinate` takes another. Results are in `results/pt_pathway_final/scfates/`.
Its protocol and decision rules were fixed before it ran (workstream 2 plan). **It supersedes the
82 → 66 → 42 list below**, which stays as a sensitivity.

**Why the list changed.** Notebook 12 required a pathway to pass a covariate-matched test **and,
separately,** a CAMERA-style correlation test. Under relabeling, the matched test alone is
anti-conservative:
- 50–70 relabeled pathways reach p ≤ 0.01, where about 15 are expected.
- The decoy z² grows with set size, which is the signature of inter-gene correlation.

Notebook 37 uses one **joint test** instead: the matched z divided by √VIF, with the VIF from
notebook 12's residual correlations. The decoy-estimated inflation for `T_spatial` (ρ ≈ 0.009)
matches the residual-correlation estimate (median ρ ≈ 0.010). This independently supports the
size of the adjustment.

**R2: specificity.** With 2 + 2 specimens, the species split and the two relabelings are the
three orthogonal ±½ contrasts. Under independent specimen deviations, each relabeling is one draw
of the species contrast's **section-level** null.

The design is a split plot:
- a species offset is judged against specimens within species, which leaves 2 df;
- a species-by-position difference is judged against specimen-by-position deviations, which leave
  2(k − 1) df.

A term shared by both human sections, such as donor, probe panel or region, cancels in the
relabelings. Relabeling therefore cannot detect it. Per-pathway specimen-level inference is
impossible: a 2-df test calls nothing.

| Strategy | Species calls | Relabeled calls | NCDR | Calls at decoy FDP ≤ 10% |
|---|---:|---:|---:|---:|
| Whole-PT pseudobulk + signed GSEA | 158 | 193, 255 | 1.42 | 0 (min FDP 1.00) |
| S1/S2/S3 pseudobulk + signed GSEA | 204 | 343, 255 | 1.47 | 0 (min FDP 0.75) |
| Gene GAM offset `T_level`, joint test | 0 | 46, 2 | – | 0 (min FDP 0.96) |
| Gene GAM × S1/S2/S3 steps, joint test | 16 | 0, 0 | 0 | 6 |
| **`T_spatial`, joint test** | **60** | **2, 0** | **0.02** | **57 (interval 23–78)** |

The pre-specified rule would have made decoy FDP the primary R2 metric. It failed because of the
pathway-score models: their decoy FDP is 0, yet they call 1,444 and 1,509 of 1,513 pathways for
species and 67–1,317 under relabeling. These scores are self-contained, so nearly every pathway's
mean differs between species by more than any section-to-section difference. Their decoy FDP
measures only that. R2 therefore keeps NCDR with the null-rate condition as its primary metric.
Decoy FDP is a secondary metric for competitive strategies.

**R3: reporting set.**

| Stage | Pathways |
|---|---:|
| Tested | 1,513 |
| Joint-test candidates | 60 |
| Not specimen-sensitive | 58 |
| Stable (≥ 5 of 7 runs) | 51 |
| **Robust** (registration effect ratio ≥ 0.7 and q ≤ 0.10; not region-sensitive against 20 power-matched removals) | **36** |
| **Core** (also detected in both species and mouse sex-unbiased) | **26** |

- **Programs.** The robust pathways form **17 gene-overlap programs**. The largest are steroid
  hormone and xenobiotic metabolism (Cyp2e1, Hsd11b1, Adh1, Ugt genes), mitochondrial and
  peroxisomal fatty-acid oxidation (Acsm3, Ehhadh, Acadm, Acaa2), vitamin and cofactor metabolism
  (Slc22a13, Slc5a8, Pank1), bile acid and peroxisome (Crot, Nudt19, Slc27a2), fat-soluble
  vitamins (Rbp4, Lpl, Apob), glutathione (Gclc, Gclm, Gss) and BCAA degradation.
- **Early/late split.** The 10 human-high pathways peak early (median 0.16); the 26 mouse-high
  pathways peak late (median 0.63).
- **Shape.** 26 are graded, 4 localized and 6 reversing. The median averaging loss is 0: most
  contributing genes differ in one direction along the PT. A whole-PT average **dilutes** these
  differences more than it cancels them.
- **"New with pseudospace" is dropped** as a class. It was defined by non-calls from screens whose
  decoy FDP is 0.75–1.00.
- **Overlap with the old list.** 27 of the 36 robust pathways were in the old 66, and 17 of the 26
  core in the old 42. The joint test applied to the relabelings flags fewer pathways as
  specimen-sensitive, so some broad programs return. Table S2 gives each pathway's specificity
  margin.

**R4: robustness.** Every robust pathway keeps ≥ 70% of its effect under:
- the anatomy-anchored coordinate;
- cortical-like mouse S3;
- removal of flagged ambient, neighbouring-segment and mapping genes;
- 5-bin matching;
- equal probe counts (workstream 3's table);
- the DPT coordinate.

Detected-in-both keeps 94% of robust pathways at that threshold, and sex-unbiased keeps 97%.

**R4: specimen-level split-plot test (E4).** This is a pseudobulk of specimen × 8 coordinate bins
with voom-style weights and moderated F. Species × position is tested against specimen × position
(22 df).
- **Gene level.** It is calibrated under relabeling: 0.6% and 0.4% of genes at p ≤ 0.01, λ ≤ 0.87.
  Its gene ranking agrees only weakly with `T_spatial` (Spearman 0.21).
- **Pathway level, pre-specified joint test.** With the pseudobulk VIF it calls nothing. Decoys
  show that this VIF overcorrects (decoy a = 1.10, ρ = 0.002, against a median VIF of 2.07).
- **Pathway level, matched only.** It calls 41 pathways, with 0 and 0 under relabeling. 15 of the
  36 robust pathways are among them, and the robust pathways' split-plot z far exceeds that of
  uncalled pathways (median 2.6 against −0.5; Mann–Whitney p = 7 × 10⁻²⁰).
- **Use.** A complement, not a replacement.

**R4: external replication (E3).** Donors are the replicates. References are Census human cortex
snRNA (6 donors, primary) and Lake/KPMP human cortex (7 donors), each against 12 male or 12 female
mice.

- **Gene level.** Our top-10% `T_spatial` genes correlate with the external species × position
  contrast at ρ 0.68 (Census) and 0.73 (Lake). Expression-matched random genes give medians of
  0.57 and 0.60, with 97.5th percentiles of 0.62 and 0.66. The agreement is therefore above
  expression-matched chance, but similar to selecting the most highly expressed genes (0.64,
  0.72).
- **Pathway level.** Robust pathways replicate far better than uncalled pathways: Mann–Whitney
  p = 2 × 10⁻¹³ (Census) and 1 × 10⁻⁹ (Lake). However, only 17 of 36 replicate within the list at
  BH ≤ 0.10 against the Census reference, and 6 against Lake. The pre-specified criterion (≥ 50%)
  is not met.
- **Claim for R4.** Measurement-level agreement, plus the comparison with uncalled pathways. It is
  not pathway-by-pathway replication.
- **Global component.** In ours and in both external human references, the human S3 − early
  contrast is a small fraction of the mouse contrast. Regressing the interaction on the mouse
  gradient gives slopes of −0.81 (ours, segment labels), −0.86 (Census) and −0.88 (Lake), with
  R² 0.60–0.76.
- **Beyond the global component.** When this uniform flattening is removed from both sides, only
  1–3 robust pathways replicate in any reference. The pathway screen therefore identifies the
  programs where mouse zonation is strongest and the human profile is flat. Pathway-specific
  positional differences beyond that global pattern are not replicated.
- **Conventional-only pathways.** Pathways called only by conventional screens (184) rarely
  replicate their average difference: 6 at nominal p ≤ 0.05, none at BH ≤ 0.10.

## The control this adds

Notebook 12 compares screens to each other. It has no measurement of how many pathways a
screen reports when there is **no** species difference to find. With two specimens per
species, the four specimens can be relabeled into two balanced mouse+human groups. Both
relabelings keep the structures, coordinate, weights, gene universe, pathway family and test
unchanged; only the species difference is removed. The species' own smooth shape enters each
relabeled model as a nuisance term, so the tested term is a difference between specimens of
the same species. This is the random-partition control of
[Lamian (Hou et al. 2023)](https://doi.org/10.1038/s41467-023-42841-y).

The **negative-control discovery ratio** (NCDR) is mean relabeled discoveries ÷ species
discoveries. A strategy is called specific when NCDR < 0.25 and the relabeled groupings call
fewer than 5% of tested pathways. Two relabelings screen for gross nonspecificity; this is not
a calibrated false-discovery rate.

## Result: only smooth position-dependent screens are specific

| Strategy | Species | Relabeled (mean) | NCDR | Specific |
|---|---:|---:|---:|:--:|
| 1 Whole-PT pseudobulk + signed GSEA | 158 | 224 | 1.42 | no |
| 2 S1/S2/S3 pseudobulk + signed GSEA | 204 | 299 | 1.47 | no |
| 3a Pathway score GAM, offset | 1444 | 1152 | 0.80 | no |
| 3b Pathway score GAM, shape | 1509 | 281 | 0.19 | no (19% of pathways) |
| 4 Gene GAM offset `T_level` | 27 | 104 | 3.85 | no |
| 5 Gene GAM total `T_total` | 76 | 131 | 1.72 | no |
| 6 Gene GAM × S1/S2/S3 steps | 67 | 42 | 0.62 | no |
| 7 Gene GAM × 6 equal bins | 100 | 42 | 0.42 | no |
| **8 Gene GAM × smooth position `T_spatial`** | **116** | **15** | **0.13** | **yes** |
| 8b `T_spatial`, size-only random sets | 344 | 57 | 0.17 | yes |
| 8c `T_spatial`, preranked GSEA | 311 | 29 | 0.09 | yes |
| Smooth position beyond S1/S2/S3 steps | 176 | 37 | 0.21 | yes |

The gene statistics behave correctly: whole-PT DESeq2 calls 9,010 genes at padj ≤ 0.05 for
species and 2–8 for the relabelings. The nonspecificity is in the pathway layer. A competitive
set test asks only whether member genes **rank** above other genes, and specimen-level
variation moves whole pathways together, which gene-permutation tests assume away
([Goeman & Bühlmann 2007](https://doi.org/10.1093/bioinformatics/btm051);
[Wu & Smyth 2012](https://doi.org/10.1093/nar/gks461)). With two specimens per species a
pathway-level average difference cannot be separated from specimen variation. Structure-level
tests of pathway scores additionally treat 12,272 structures as replicates
([Squair et al. 2021](https://doi.org/10.1038/s41467-021-25960-2)).

Within-species controls agree: comparing the two mice gives 48 `T_level` pathways and 16
`T_spatial`; the two human sections give 94 and 0. Specimen differences are offsets, not
positional remodeling.

**Decision.** Primary analysis = gene-level `T_spatial` from the nested GAM on the scFates
coordinate, tested with notebook 12's covariate-matched rank-AUC. It keeps notebook 12's list
unchanged, is the only strategy passing specificity, answers the positional question, and
reads positions back from D(s)/S(s). Constant-offset results (whole-PT, S1/S2/S3, `T_level`)
are reported as context only: they may be real biology, but this cohort cannot separate them
from specimen variation. Conventional analyses are not wrong; they ask a question this design
cannot answer specifically.

## Final list (notebook 31; superseded by notebook 37 above)

129 notebook-12 spatial candidates → 104 not specimen-sensitive → 84 stable (retained in ≥5 of
7 planned sensitivity runs) → **82 correlation-supported**. Versus conventional screens:
54 new with pseudospace (no conventional species call), 26 whose conventional call is also
produced by a relabeling, 1 conventional call that no relabeling produces, 1 not testable
conventionally.

Themes (keyword rules applied after the list was frozen): amino-acid metabolism (12),
xenobiotic & drug metabolism (12), sterol/steroid/bile acid (12), vitamins & cofactors (9),
fatty-acid oxidation & peroxisome (6), glucose & energy metabolism (6), solute transport (6),
other/signaling (19).

70 of 11,071 genes (0.6%) reverse species direction along the PT under the display rule
(|Z| ≥ 3 and |difference| ≥ 0.25 over ≥3 contiguous grid points on each side). A whole-PT
comparison reports one direction per gene and cannot represent these.

## The paper story (notebook 32)

`analysis/notebooks/32_pt_pathway_story.ipynb` tests nothing new. It picks readable examples from
the frozen 82-pathway list, shows their genes, and makes the conventional contrast concrete.

**Where the species differ, not whether.** Among the 82 confident pathways, the 22 that are
human-high at their divergence peak have a median peak at 0.13 on the coordinate (S1). The 60
mouse-high ones have a median peak at 0.65 (S2/S3 boundary). Pathways share genes, so this is a
descriptive split, not a test. A whole-PT or three-segment summary assumes every pathway differs
in the same place.

**Flagships.** Chosen from pathways that are new with pseudospace and retained in all 7
sensitivity runs, one per theme where possible: Drug ADME, mitochondrial fatty-acid
β-oxidation, valine/leucine/isoleucine degradation, glutathione metabolism, retinoid
metabolism & transport. Each figure shows D(s), member-gene differences, fitted curves with
specimen means, and conventional log2FC.

Gene-level examples worth naming in the text (both specimens of each species agree):

- **Slc22a6 (OAT1)**, the main renal drug-uptake transporter: whole-PT DESeq2 reports no
  species difference (log2FC −0.09, ns), while the segments differ significantly in opposite
  directions (S1 +0.32, S2 −1.15, S3 +2.84). Mouse peaks sharply in S2; human stays flat. This
  is confirmed without the coordinate (segment labels only) and survives a per-gene shift.
  Slc13a3 and Cyp24a1 show the same whole-PT-null, regionally-opposite pattern (notebook 33,
  Figure 9).
- **Ugt3a1**: human-high in S1, mouse-high in S2/S3; segment DESeq2 confirms both signs.
  Whole-PT DESeq2 reports mouse-high (log2FC −1.3, significant), which is wrong for S1.
- **Cyp2e1**: mouse-expressed through S1–S2, falling in S3; essentially absent in human PT.
  Mouse renal CYP2E1 is androgen-regulated and bioactivates acetaminophen in S1/S2, while human
  kidney shows little CYP2E1 activity. This is a translational point for rodent nephrotoxicity
  models.
- **Gatm** (creatine synthesis, first step): human high along the whole PT; mouse high only
  early, falling toward S3. Segment log2FC grows from +1.2 (S1) to +7.0 (S3), so the species
  difference is a zonation shift rather than a uniform offset.
- **Acsm3, Crot, Nudt19, Slc27a2**: fatty-acid oxidation and peroxisomal genes that rise along
  the PT in mouse and stay flat in human. Rodent-specific PPARα-driven peroxisomal programs
  are a known species difference; Nudt19 is a mouse-kidney-specific peroxisomal CoA
  diphosphohydrolase.
- **Reversing genes (70)**: e.g. Acaa2, Dcxr, Cndp2, Slc22a6, Slc22a7. 33 belong to a
  confident pathway. 54 are confirmed by segment DESeq2 with significant opposite signs, and
  24 of those have a non-significant whole-PT comparison. 21 also survive a per-gene shift of
  the human curve; the others are better described as a species shift of the gene's expression
  boundary than as a change of direction. A single log2FC per gene cannot represent either.

Outputs: `results/pt_pathway_story/` (`figures/fig5_example_*`, `fig6_direction_reversing_genes`,
`fig7_confident_pathway_landscape`, `example_gene_evidence.csv`,
`flagship_member_genes.csv`, `direction_reversing_genes.csv`,
`confident_pathways_for_paper.csv`).

## Robustness to coordinate registration and sampled region (notebook 33; the counts here refer to the superseded list)

`analysis/notebooks/33_pt_pathway_registration_sensitivity.ipynb` tests the two alternatives
the relabeling control cannot detect, because both are species-level effects.

**Misregistration.** The reviewed S2→S3 transition falls at 0.63–0.64 on the coordinate in all
four specimens. The S1→S2 transition falls at 0.27–0.28 in mouse but 0.35–0.38 in human, so human
S1 spans more of the coordinate. A global affine stretch/shift barely improves human–mouse curve
agreement (mean shape correlation 0.18 → 0.20), so the curves genuinely differ in shape. Warping
each human section so its segment transitions land on the mouse ones and refitting everything:
gene `T_spatial` Spearman 0.98, 129 → 131 spatial calls, **79 of 82 confident and all 54
new-with-pseudospace pathways retained**. The three lost (OXPHOS, TCA cycle, glyoxylate) are
broad energy pathways with non-specific conventional calls.

**Sampled region.** Both human sections are cortex; mouse sections include the outer stripe.
Using distance to the three nearest reviewed glomeruli (centroids verified against the fine
GeoJSON), only 15–17% of mouse S3 lies within the depth of 95% of the same specimen's S1. After
removing the deeper mouse S3 and refitting, 69 of 82 confident and 44 of 54 new pathways are
retained. Losses may reflect region or the power lost by removing most of mouse S3.

**Paper-quotable set.** 66 of 82 confident pathways survive both checks (44 new with
pseudospace, 21 with non-specific conventional calls, 1 not testable conventionally). 13 are
region-sensitive, mostly late sterol/bile-acid terms; 3 are registration-sensitive. All five
flagships survive both (q ≈ 0.002). The early/late split holds in the robust set: 19 human-high
pathways peak at median 0.145, 47 mouse-high at 0.624. Table:
`results/pt_pathway_registration_sensitivity/confident_pathways_with_robustness.csv`.

**Trajectory method.** Notebook 12's earlier run on notebook 03's DPT coordinate (different
trajectory method and PT labels; `results/pt_pathway_remodeling/`) recovers 61 of the 66 robust
pathways, 42 of the 44 robust new ones, and all five flagships. Spatial pathway effects
correlate at Spearman 0.92 between the DPT and scFates runs (126 and 129 calls, 101 shared).

**Detection in both species and mouse sex bias (notebook 33, sections 10–11).**
- 47 of the 66 robust pathways keep their call when only genes detected in at least 5% of
  structures in each species are used ("shared-gene zonation"). 14 depend on genes that one
  species barely detects, which could be true species-specific expression or probe gaps.
- Removing all 1,211 genes that are sex-biased in mouse PT (Xiong et al. 2023 public tables)
  keeps 51 of 61 testable robust pathways and all five flagships. The drug-metabolism sets
  (Aspirin and Paracetamol ADME, Drug metabolism) depend on sex-biased genes.
- **Core set: 42 pathways** pass every check (robust, shared-gene, sex-robust); 25 are new with
  pseudospace.

Literature status of these findings: [pt-pathway-literature-novelty.md](pt-pathway-literature-novelty.md).

### What this changes for notebook 14

Notebook 14's headline benchmark (66 discrete → 89 continuous `T_total` calls) uses a
constant-offset-dominated statistic. Under relabeling, `T_total` reports more pathways than
for species (NCDR 1.72), so that count should not carry the "uncovers more pathways" claim.
Use the `T_spatial` specificity result and the 44 robust new-with-pseudospace pathways instead.

## Caveats to keep beside any number here

- Two mice and two sections from one human donor. Specimens, not structures, are the
  replication units; no result is population-level human inference.
- The relabeling control uses two partitions only.
- Themes, drivers and reversal counts are display rules on working statistics, not tests.
- Both human slices are healthy cortex despite the cortex/medulla labels; mouse S3 is mostly
  outer stripe, so late-PT species differences carry a region caveat (see robustness section).

Outputs: `results/pt_pathway_method_selection/` (tables, `figures/`, `run_manifest.json`).
