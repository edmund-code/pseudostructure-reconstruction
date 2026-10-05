# Manuscript draft v0.4

> **Editorial note (remove before submission).** v0.4 restructures v0.3 after the independent review
> (`docs/paper/referee_report_v03.md`) and the decisions in `docs/paper/outline.md`. Supplementary
> Notes and Supplementary Methods are in `supplementary_notes_v04.md`, and the figure plan is in
> `figure_plan_v04.md`. Fixes and unverified numbers are listed in `draft_v04_changelog.md`.
>
> - **Markers.**
>   - **[PENDING WS6]**: slots for the index-sensitivity and agreement-by-strength analysis, still
>     running. No number in a slot is filled. The notebook 47 and 48 results are filled in.
>   - **[VERIFY]**: needs lab records, or a number no saved table supports.
>   - **[CITE]**: needs a reference.
> - **Standing caveat**, which sits in short form beside each result:
>   - two control mice and two cortex sections from one male human donor; all specimens male;
>   - the human section labelled "medulla" is cortex, so mouse outer-stripe S3 has no human
>     counterpart;
>   - the Visium HD probe panels differ between species;
>   - all findings are descriptive.

---

## 1. Title options

The final wording depends on [PENDING WS6]: whether agreement rises with gradient strength.

1. Human and mouse proximal tubule share the zonation of their most strongly zonated genes and little
   else
2. Measured against same-species reproducibility, human and mouse proximal tubule share little axial
   zonation beyond their strongest segment genes
3. A reproducibility-calibrated comparison of human and mouse proximal tubule zonation

## 2. Abstract (≈ 270 words, excluding slots)

Mouse proximal tubule (PT) is the main experimental model for human PT transport, metabolism and
injury. Whether the two species arrange their genes alike along the S1→S3 axis has not been measured
against how reproducible such gradients are within one species.

We profiled two control mouse kidney sections and two cortical sections from one human donor with
Visium HD, segmented tubule cross-sections on the paired histology, and compared segment gradients of
gene expression with same-species reproducibility across independent datasets. Relative to that
reproducibility, human–mouse agreement was low: a conservation index of 0.17 (95% interval 0.12–0.22)
for S1→S2 and 0.28 (0.22–0.33) for S3 against early PT [PENDING WS6: range under index
sensitivity], and 0.27–0.28 in public single-nucleus atlases alone. Agreement was concentrated in a
small set of strongly zonated genes [PENDING WS6: agreement by gradient strength]; among genes zonated
in both species, 72% agreed in direction. Confidently species-only genes were few, and whether
mouse-only genes outnumber human-only ones depends on whether zonation is judged on an absolute scale
or relative to each species' amplitude.

A second finding concerns method. With two specimens per species, pathway screens that compare
averages could not be shown to be specific: relabeling the specimens into mixed-species groups gave
as many calls as the species contrast, in our cohort and in 100 cross-species designs drawn from
public atlases. Screens for position-dependent differences passed this control when the specimens
within each group shared their positional biology. Most of their 22 robust pathways summarise mouse
zonation that the human sections do not share.

With one human donor, these results are descriptive. Outside a small set of strongly zonated genes,
mouse PT is an incomplete guide to how human PT zonates its genes.

---

## 3. Introduction

**Translation of PT biology depends on knowing where species differ.**
- *Function.* The PT reabsorbs about 60–70% of filtered water and NaCl, more of the filtered
  bicarbonate and nearly all filtered nutrients. It also secretes solutes, produces hormones and
  carries out most of the kidney's metabolism (Curthoys & Moe 2014).
- *Systemic roles.* During fasting and stress the kidney supplies up to 40% of systemic glucose
  production, and PT gluconeogenesis is impaired in acute kidney injury (AKI) (Legouis et al. 2020).
  Defective fatty-acid oxidation in tubular epithelium contributes to fibrosis (Kang et al. 2015). The
  SLC22 organic anion transporters handle common drugs, toxins, nutrients and putative uremic toxins
  (Nigam et al. 2015).
- *Injury.* The PT is regarded as the primary target of injury in AKI and chronic kidney disease
  (Chevalier 2016), and as the primary site of drug-induced kidney injury (Maass et al. 2019).
- *Translation from rodents.* Much of this biology is studied in rodents:
  - no intervention has yet prevented AKI in patients, which has raised the question of whether animal
    models predict human responses (de Caestecker et al. 2015);
  - rodent AKI models do not fully recapitulate the human disease (Fu et al. 2024), and experimental
    models differ in where along the nephron tubular injury falls (Heyman et al. 2010);
  - preclinical models often fail to predict drug nephrotoxicity (Maass et al. 2019);
  - species differences in renal drug transporters and metabolising enzymes contribute to poor
    translation (Basit et al. 2019; Thakur et al. 2024).

**The PT is zonated, and some of its zonation is known to be conserved.**
- *Distribution of function.* Along its length the PT is divided into S1, S2 and S3, and functions are
  distributed along this axis. In mouse, SGLT2 sits in the early PT (Vallon et al. 2011). In isolated
  human PT segments, S2 and S3 make more glucose from lactate than S1 (Conjard et al. 2001).
- *Segment-resolved data.* Segment-resolved transcriptomes exist for microdissected rat (Lee et al.
  2015) and mouse tubules (Chen et al. 2021), and protein for rat (Limbutara et al. 2020). Single-cell
  atlases describe PT diversity in mouse (Ransick et al. 2019) and human (Lake et al. 2023).
- *Human protein maps.*
  - They place SGLT2 in S1/S2 and SGLT1 in S3, much as in rodents (Vrhovac et al. 2015).
  - They place AQP1 highest in the proximal straight tubule (Maunsbach et al. 1997).
  - They place all three OATs more intensely in S1/S2 than in medullary-ray S3, with outer-stripe S3
    unstained (Breljak et al. 2016).
- *Other programs.* In mouse, sexually dimorphic gene activity maps predominantly to PT segments, but
  only a limited set of genes shows conserved sex-linked regulation in human (Xiong et al. 2023). The
  largest mouse PT sex differences are in S2/S3 (Chen et al. 2023).
- *The gap.* These observations come from different studies, assays and segment definitions.
  Computational models of the human nephron still take most transporter densities from rat (Layton et
  al. 2019).

**Comparisons of cell types do not say how similar zonation is.**
- *Atlases.* Cross-species single-cell atlases align conserved cell states between human and rodent
  kidney (Klötzer et al. 2025). A unified mouse–human kidney atlas found little overlap of single-cell
  expression changes but agreement at the pathway level (Zhou et al. 2023). These comparisons ask
  which cells correspond, not whether genes change alike along a segment.
- *The missing benchmark.* A low cross-species correlation of segment gradients means "different
  genes" only if two datasets of the same species agree much better. Segment annotation, dissociation
  and platform all differ between datasets, so that benchmark has to be measured.
- *Continuous zonation elsewhere.* Reconstructions in other organs show what segment bins can hide.
  About half of mouse liver genes are zonated along the lobule (Halpern et al. 2017), and enterocyte
  functions are broadly zonated along the intestinal villus (Moor et al. 2018).

**Spatial transcriptomics keeps tissue context.**
- *In kidney.* It has mapped cell types and injury niches together with single-nucleus data (Melo
  Ferreira et al. 2021), with Slide-seqV2 (Marshall et al. 2022), and in the human kidney atlas (Lake
  et al. 2023).
- *Visium HD.* It measures the whole transcriptome in 2-µm bins (Oliveira et al. 2025), which can be
  aggregated into cells or larger structures (Polański et al. 2024).
- *Position along the tubule.* A section cuts each tubule at an arbitrary point, so a structure's
  position along the nephron can be read from segment markers or reconstructed. Approaches include:
  - landmark genes (Halpern et al. 2017; Moor et al. 2018);
  - expression-only spatial reconstruction (Nitzan et al. 2019);
  - unrolling of spatial data (Guo et al. 2026);
  - trajectory methods (Saelens et al. 2019), such as diffusion pseudotime (Haghverdi et al. 2016) and
    principal curves (Faure et al. 2023).
- *Condition comparisons.* Condition-by-pseudotime tests exist (Van den Berge et al. 2020; Roux de
  Bézieux et al. 2024), but they treat cells, not specimens, as the units.

**Few specimens make species comparisons hard to separate from specimen differences.**
- *Replication.* Cells or structures from one specimen are not independent replicates. Methods that
  ignore variation between biological replicates can report hundreds of differentially expressed
  genes where no biological difference exists (Squair et al. 2021; Zimmerman et al. 2021).
- *Pseudobulk.* Aggregating to pseudobulk has been described both as conservative and underpowered
  (Zimmerman et al. 2021) and as the better balanced choice (Murphy & Skene 2022). Both approaches have
  been benchmarked for multi-sample designs (Crowell et al. 2020).
- *Random partitions.* For pseudotime, Lamian randomly partitioned samples into groups with no
  expected difference; methods that ignored sample variation then reported thousands of
  false-positive genes (Hou et al. 2023).
- *Correlated gene sets.* Competitive gene-set tests that treat genes as independent are
  anti-conservative when genes in a set are correlated (Goeman & Bühlmann 2007; Wu & Smyth 2012).
- *Consequence.* Human tissue is scarce, and spatial studies often have a few sections per species, so
  an analysis has to show that its calls are specific to the species contrast.

**This study.** We profiled two control mouse kidney sections and two cortex sections from one human
donor with Visium HD. We segmented tubule cross-sections from the paired H&E images and aggregated
2-µm bins into tubule-level profiles in a shared one-to-one ortholog space. Then:
- we validated the segment labels anatomically and by reference transfer;
- we compared human and mouse segment gradients against same-species reproducibility measured across
  independent datasets, in our data and in public atlases alone;
- we called each gene's zonation separately in each species, confirming every zonated call in an
  independent atlas of the same species;
- we asked when pathway screens are specific with two specimens per species, using a
  specimen-relabeling control in our cohort and many-draw nulls from public atlases;
- we checked lead genes against microdissected mouse tubules of both sexes, mouse single-nucleus data,
  a cortex-only human atlas and rat segment proteomics.

The two species share the zonation of a small set of strongly zonated genes, in which transporters are
over-represented, and agree little elsewhere.

### Box 1 · Glossary

- **S1, S2, S3; early PT.** Proximal-tubule segments in order from the glomerulus. "Early" is S1 + S2.
- **S3c.** Cortical-like mouse S3: mouse S3 structures no deeper than the 95th percentile of S1
  depth in the same section. Depth is the mean distance to the three nearest glomeruli. S3c is used
  to match human cortex.
- **Gradient.** A gene's log2 pseudobulk ratio between two segments (S2 − S1, S3 − early or
  S3c − early).
- **Ceiling.** The noise-corrected correlation of gradients between two independent datasets of the
  same species.
- **Conservation index.** The noise-corrected human–mouse correlation of gradients divided by the
  geometric mean of the human and mouse ceilings. A value of 1 means the species agree as well as two
  datasets of one species.
- **Zonation classes.**
  - *Conserved:* zonated in the same direction in both species.
  - *Reversal:* zonated in opposite directions.
  - *Mouse-only, human-only:* zonated in one species and flat in the other, with flatness confirmed in
    that species' reference atlas.
  - *Indeterminate:* neither call is possible with two specimens.
- **Relabeling.** Regrouping the four specimens into the two balanced mixed-species pairs. A
  relabeled contrast has no species difference.
- **NCDR (negative-control discovery ratio).** The mean number of pathways called under relabeling
  divided by the number called for species. It measures vulnerability to specimen noise. It is not a
  false-discovery proportion.
- **T_spatial.** A gene's statistic for a species difference that changes with position along the PT
  coordinate. It is used to rank genes and is not a gene-level test.
- **Step screen.** The same test with position replaced by S1/S2/S3 labels.
- **Joint test.** A pathway's matched rank test of its members' statistics, corrected for
  correlation among members.
- **Robust pathway.** A joint-test call that is not specimen-sensitive, is stable across refits, and
  survives anatomical registration and removal of deep mouse S3.

---

## 4. Results

### R1 · Segmented PT structures and their segment labels

**Cohort and structures** (Figure 1a).
- *Specimens.* Two control mouse kidney sections, which include cortex and outer medulla, and two
  sections of renal cortex from one male human donor. All four specimens are male. The human section
  named "medulla" at source contains cortex only.
- *Measurement.* Visium HD measured the whole transcriptome in 2-µm bins, with species-specific probe
  panels.
- *Segmentation and matrices.* Tubule cross-sections and glomeruli were segmented on the paired H&E
  images (Supplementary Note 7), and bins were summed per structure. After structure filters, 26,839
  structures remained.
- *PT structures.* Two passes of cross-species integration in a one-to-one ortholog space, each
  followed by clustering, identified 12,866 PT structures: mouse 2,852 and 2,808; human 3,438 and
  3,768.

**What the segment labels are.** Each PT structure carries an S1, S2 or S3 label. The labels are
memberships of joint cross-species clusters, named by marker review. They are not anatomical
annotations made separately in each species. We tested them without the joint embedding (Figure 1b–d;
Supplementary Figure S8).

1. **Distance to the nearest glomerulus** followed the label order in both species, but separated human
   S1 from S2 only weakly. Spearman correlation of label order with distance:
   - mouse 0.65, with median distances (S1 / S2 / S3) 56 / 93 / 297 µm;
   - human 0.23, with 131 / 148 / 228 µm.

   The human distances assume a pixel size of 0.44 µm, which is unconfirmed [VERIFY against scanner
   metadata]; at the recorded 0.5 µm they are 149 / 168 / 259 µm. Mouse distances are descriptive,
   because the mouse glomerulus polygons fail the pre-specified size rule.
2. **Glomerulus-attached PT.** Within 10 µm of a glomerulus, S1 was enriched less in human (odds ratio
   2.0, 237 polygons) than in mouse (4.5, 543 polygons). At the pre-specified 2 µm, too few human
   polygons were attached to judge (9).
3. **Reference transfer.** Labels transferred by expression alone from reference segment profiles
   agreed with ours in 71% of human structures (κ 0.54), against 87% (κ 0.80) and 92% (κ 0.87) in mouse.
   - Human S1 and S3 transferred mostly to themselves (81% and 99%).
   - 65% of human "S2" structures transferred to the reference's S3. The reference annotates few S2
     nuclei, so part of this may reflect the reference.

By the pre-specified reading, the human labels have weak reference support, and human S2 is the least
supported. We therefore repeat the main comparison with transferred labels (R2), and we flag gene
claims that rest on S2 − S1 alone (R3, R6).

**Contamination and capture.**
- *Spillover.* Counts from non-PT markers made up a small and similar share in both species: 0.07–0.45%
  in human and 0.13–0.29% in mouse, by segment.
- *Polygons without bins.* Human polygons without bins lie outside the area covered by bins, consistent
  with tissue beyond the capture area.

**A continuous coordinate as a supporting tool** (Figure 1e; Supplementary Note 1). We also ordered the
PT structures along a nonbranching principal curve fitted to the integrated embedding.
- *What it recovers.* Measured zonation, about as well as the segment labels: its agreement with
  external S3-versus-early contrasts was 0.74 in mouse and 0.25 in human, against 0.73 and 0.27 for the
  labels.
- *What it does not resolve.* It does not reliably order structures within segments, and in human S1
  its order depends on read depth.
- *How we use it.* Only to replicate measurements across positions within each specimen in the
  pathway screen (R4, R5).

The zonation comparisons (R2, R3) and the lead genes (R6) use segment labels only.

### R2 · Measured against same-species reproducibility, human–mouse agreement of PT gradients is low

**Same-species ceilings.** A low human–mouse correlation of gene gradients means "different genes" only
relative to what two datasets of the same species give. For S2 − S1 and S3 − early (early = S1 + S2),
we computed noise-corrected correlations of gene gradients for every pair of datasets with replicates
(Figure 2a; Table 1).

| Pair | S2 − S1 | S3 − early |
|---|---|---|
| Within human (ours, cortex snRNA atlases) | 0.70 | 0.70–0.93 |
| Within mouse (ours, snRNA by sex, microdissection) | 0.63–0.88 | 0.73–0.93 |
| Human against mouse | 0.09–0.20 | 0.13–0.24 |

**The conservation index** divides the human–mouse correlation by the geometric mean of a human ceiling
(our sections against the cortex atlas) and a mouse ceiling (our mice against male mouse snRNA).

| Contrast | Our data (95% interval) | Public atlases alone |
|---|---|---|
| S2 − S1 | 0.17 (0.12–0.22) | 0.27 |
| S3 − early | 0.28 (0.22–0.33) | 0.28 |
| S3c − early (region-matched) | 0.29 (0.23–0.35) | 0.27 |

Every upper bound lies below 0.5, so the pre-specified wording "largely different genes" applies. Three
properties of the index bear on how it is read:
- *The S2 − S1 human ceiling rests on one dataset pair* (ours against the cortex atlas), because the
  second human atlas has no S2.
- *Platforms differ.* The ceilings compare across platforms (Visium HD against single-nucleus RNA),
  while the index's numerator compares within one platform. This probably biases the index upward,
  which is conservative for the claim.
- *Intervals.* They come from gene bootstraps. These treat genes as independent, although gradients
  are correlated within modules, and they condition on the specimens.

[PENDING WS6: index sensitivity. Report the index under winsorisation at q = 0, 0.01 (published) and
0.05, with reliability computed on the same winsorised values; module- or pathway-block bootstrap
intervals; a donor bootstrap of the external index. State whether every upper bound stays below 0.5.]

**Agreement is concentrated in a small set of genes** (Figure 2b).
- *Direction.* Among the genes called zonated in both species (R3), 196 of 271 (72%) agree in
  direction, against 50% expected without association. A low index therefore does not mean that most
  zonated genes disagree.
- *Gradient strength.* [PENDING WS6: agreement by gradient strength. Correlation and sign agreement
  of human and mouse gradients in bins of |mouse gradient|, with per-bin same-species ceilings, for
  S2 − S1, S3 − early and S3c − early. Expected reading if it holds: agreement concentrated among the
  most strongly zonated genes and absent from the weakly zonated majority.]

**Robustness to the human labels** (Figure 2c; Supplementary Figure S9).
- *Design.* We transferred human labels from the cortex atlas, alone or with mouse labels transferred
  from microdissection, using one half of the genes, and evaluated the index on the other half.
- *Pre-specified contrasts.* The index stayed at or below 0.31 for S2 − S1 and S3 − early in every
  setting and half (largest upper bound 0.38).
- *Region-matched contrast.* For S3c − early it reached 0.33 (upper bound 0.43).
- *Human ceiling.* With transferred labels, the S2 − S1 human ceiling fell from 0.64–0.67 to
  0.53–0.55.

**Robustness to the genes that built the integration** (Figure 2c).
- *Reviewed labels.* Without the 611 genes used to build the integrated embedding, the index was
  0.15 (upper bound 0.20) for S2 − S1, 0.22 (0.27) for S3 − early and 0.25 (0.31) for S3c − early.
- *Cross-fitted labels.* It stayed at or below 0.24 for S2 − S1 and S3 − early (largest upper bound
  0.33), and at or below 0.29 for S3c − early (upper bound 0.38).

**Other species as a benchmark.** Rat microdissection data failed the pre-specified reliability rule
(replicate reliability 0.03 for S2 − S1 and 0.45 for S3 − early, against ≥ 0.5) and were not used. We
therefore cannot say whether agreement this low is unusual for two mammals.

**Human gradients are also smaller** (Figure 2d; Supplementary Figure S5).
- *Spread.* The noise-corrected spread of S2 − S1 gradients in human was 0.58 of the mouse value in our
  data and 0.47 in public single-nucleus data. For S3c − early it was 0.84. With transferred human
  labels, the S2 − S1 value was about 0.50.
- *Ratios on a reference axis depend on the axis.* Projected on mouse reference gradients, human
  amplitude was 0.15–0.23 of mouse. Projected on the human cortex atlas, it was 1.96 (1.41–3.05).
  With low cross-species correlation, these are projection ratios, so the reference-free spread is the
  amplitude measure.
- *Log compression.* It does not explain the difference: the ratio was 0.18 in the top expression
  tertile and 0.16 for genes with equal probe counts.

### R3 · Which genes are shared

**Calling zonation separately in each species** (Figure 3a).
- *Calls.* Each species' gradients were tested against that species' own specimen noise.
- *Confirmation.* Every zonated call had to be confirmed in an independent atlas of the same species:
  12 male mouse single-nucleus donors, and 7 healthy human cortex donors.
- *Species-only calls.* A gene was called zonated in one species only when the other species was flat
  in our data and flat-compatible in its own atlas.
- *Thresholds.* The effect floor and flat margin were scaled to each species' amplitude (R2).
- *Counts.* 7,526 genes have equal probe counts in both panels. Of these, 7,407 could be tested in
  S2 − S1 or S3c − early.

| Class (union of S2 − S1 and S3c − early) | Genes |
|---|---:|
| Conserved (same direction) | 196 |
| Reversal (opposite directions) | 75 |
| Mouse-only | 41 |
| Human-only | 76 |
| Flat in at least one species, zonated in neither | 712 |
| Not classifiable with two specimens per species | 6,307 |

**Agreement among genes zonated in both species.**
- *Direction.* 72% of the 271 genes zonated in both species agree in direction.
- *Reversals.* They make up 28% of the 271, so they are not rare.
  - *S2 − S1 reversals depend on the human labels.* There are 37 with reviewed labels; 24 and 15 with
    labels transferred on all genes; and 5–11 in cross-fitted gene halves, against 15 and 21 for
    reviewed labels on the same halves.
  - *S3c − early reversals are stable.* There are 41 with reviewed and with transferred labels, and
    11–26 in cross-fitted halves, against 16 and 24 for reviewed labels.

**Transporters are over-represented among conserved genes, but they are not most of them** (Figure 3b).
- *Share.* 33 of the 196 conserved genes (17%) are transporters (Slc, Abc and Aqp genes).
- *Power.* Transporters are highly expressed and strongly zonated, so 40 of the 43 zonated transporters
  (93%) are zonated in both species, against 271 of 388 zonated genes overall (70%).
- *Direction agreement.* Among genes zonated in both species, transporters agree in direction more
  often than other genes: 33 of 40 (83%) against 196 of 271 (72%).
- *Examples* of conserved zonation: SGLT2 and SGLT1 (Slc5a2, Slc5a1), Slc5a8, Slc7a13, Slc6a18,
  Slc22a7 (OAT2), Slc22a12, Slc34a3, and the intrarenal renin–angiotensin genes Agt, Ace and Enpep.
- *Human protein maps* agree with these calls for SGLT1/2 but not for two of the three OATs (Breljak et
  al. 2016):
  - OAT1 protein is stronger in S1/S2 than in medullary-ray S3, but human OAT1 (SLC22A6) mRNA is
    indeterminate in our data;
  - OAT2 protein is also S1/S2 > S3, but OAT2 mRNA rises toward S3 in our sections and in the cortex
    atlas.

  Possible explanations are mRNA–protein differences, antibody specificity (the antibodies were
  validated only on transfected cells), and the interindividual variability the authors report.
- *Conserved zonation is not conserved usage.* In adult human kidney, SLC34A3 carries about 40% of
  sodium–phosphate cotransport but is negligible in mouse (Fernandes et al. 2026), although its
  zonation is conserved here.

**The balance of species-only genes depends on the scale** (Figure 3c).

| Thresholds | Species-only call needs a flat reference atlas | Mouse-only : human-only |
|---|---|---|
| Absolute | no | 247 : 66 |
| Absolute | yes | 187 : 31 |
| Scaled to each species' amplitude | no | 59 : 141 |
| Scaled to each species' amplitude (primary) | yes | 41 : 76 |

- *Absolute scale.* Judged on the same log2 scale in both species, more genes are zonated only in
  mouse.
- *Relative scale.* Judged relative to each species' own amplitude, there is no mouse excess.
- *Which is right.* Scaling is appropriate if the lower human amplitude is technical (capture, bin
  coverage or label impurity). It is not appropriate if the lower amplitude is biological. We cannot
  distinguish the two.
- *Failed pre-specified check.* With labels transferred from reference atlases on held-out genes, the
  pre-specified criterion "species-only counts unskewed" (ratio < 2 in both contrasts, every setting
  and half) failed. One check of 12 gave 17 mouse-only against 4 human-only; the same setting in the
  other gene half gave 7 against 22.

Either way, confidently species-only genes are few.
- *Mouse-only examples:* Fah, Comt, Maoa, Igfbp4 and Ebp; Hsd17b4 is mouse-only with reviewed labels
  and indeterminate with transferred labels.
- *Human-only examples:* Sel1l3 and Ace2; Slc9a3 is human-only in S2 − S1 and changes class with
  transferred labels.
- *Pathways.* We named an enrichment only if it had q ≤ 0.10, was not driven by expression, and was
  carried by at least three class genes (Supplementary Note 3).
  - *Mouse-only:* tyrosine catabolism (Comt, Fah, Maoa; q = 2 × 10⁻⁶). Sterol synthesis rests on Ebp
    alone and peroxisomal β-oxidation on Hsd17b4 alone, so neither is a program-level result.
  - *Human-only:* the two enrichments that pass are generic labels on overlapping three-gene sets, so we
    name the genes instead: Ahcyl1, Hipk2, Wwtr1, Prkaa2 and Taf4.
  - *Conserved:* enriched for amino-acid and organic-anion transport, bile secretion, apical-surface
    and renin-secretion genes.
- *Integration genes.* The 611 genes used to build the integrated embedding make up 12% of conserved
  genes (24 of 196), 5% of reversals, 1% of human-only genes and none of the mouse-only genes, against
  6% of genes that could not be classified. Removing them lowers the index slightly (R2).
- *Caveat.* Mouse S3c calls are confirmed against reference atlases whose S3 lies mostly in the outer
  stripe.

**Log compression is ruled out; tissue state does not explain the pattern in the data available**
(Supplementary Note 4; Supplementary Figure S6).
- *Our donor.* On a PT injury score calibrated across platforms in mouse, the donor lies at the median
  of 14 healthy reference donors.
- *Injury and amplitude.* Amplitude falls with injury among reference donors, but healthy human donors
  remain far flatter than healthy mice.
- *Time of day.* Mouse zonation of the lead genes and of the fatty-acid oxidation, glutathione and
  tyrosine programs is the same at rest and in the active phase. Sterol-synthesis zonation is not.
- *Fasting.* It changes PPARα-target and sterol genes at whole-kidney scale by amounts comparable to
  the species differences. No segment-resolved fasted-kidney data exist.
- *A state marker.* Human PT expresses the fasting- and stress-responsive gene PDK4 strongly (94% of
  human structures, against 7% of mouse).

### R4 · With two specimens per species, average-level pathway screens cannot be shown specific by relabeling

**The control.**
- *Relabelings.* With two specimens per species, a pathway screen can only be judged against specimen
  variation. The four specimens allow exactly two balanced mixed-species relabelings, and each should
  call nothing.
- *What it cannot see.* A term shared by both human sections (donor, probe panel, sampled region)
  cancels in the relabelings and is not covered.
- *The rule.* We call a strategy specific when NCDR < 0.25 and each relabeling calls fewer than 5% of
  the tested pathways.

**Our cohort** (Figure 4a, b; Supplementary Note 2; Supplementary Figure S3).

| Strategy | Species calls | Relabeled calls | NCDR |
|---|---:|---:|---:|
| Whole-PT pseudobulk, signed GSEA | 158 | 193, 255 | 1.42 |
| S1/S2/S3 pseudobulk, signed GSEA | 204 | 343, 255 | 1.47 |
| Gene-level offset, joint test | 0 | 46, 2 | – |
| Step screen without a coordinate (segment intercepts), joint test | 50 | 0, 0 | 0 |
| Step screen added to the coordinate spline, joint test | 16 | 0, 0 | 0 |
| Position along the coordinate (T_spatial), joint test | 60 | 2, 0 | 0.02 |

Average-level screens called more pathways for the relabelings than for species. Screens for
position-dependent differences called few or none for the relabelings.
- *Two step-screen designs.* Both test species × S1/S2/S3 steps beyond a species offset, on the same
  structures, genes and pathways. The coordinate-free design uses segment intercepts as the shared
  shape; it is the design the external nulls calibrate. The spline-based design adds the steps to a
  shared spline along the coordinate. 15 of its 16 calls are among the 50.
- *The coordinate screen* adds calls (60), not specificity.

**Many-draw nulls from public atlases** (Figure 4c, d). We drew 750 same-species and 100 cross-species
2 + 2 designs from single-nucleus atlases with segment labels (24 mice; 13 human donors). The atlases
have no coordinate, so these nulls compare average-level screens with the step screen, not with the
coordinate.
- *Average-level screens* called pathways in almost every same-species split (whole-PT GSEA: median 13
  to 40 per pool). In cross-species draws they called as many under relabeling as for species (median
  NCDR 1.34).
- *The step screen* called none in adult human cortex (median 0) and few in mice (medians 5 and 11). In
  cross-species draws it called a median of 15 pathways for species against 0.5 under relabeling
  (NCDR 0.05). 80% of draws met the specificity rule, exactly the pre-specified 80%.
- *Calibration.* At p ≤ 0.01, the step screen called 1.0% and 0.6% of pathways in the human pools and
  3.8–4.8% in the mouse pools, against 1% expected. Average-level screens called 2.4–7.5%.
- *Our cohort's step-screen relabelings* (0 and 0) sit at the median of the external relabeling
  distribution.

**"Cannot be shown specific" is not "wrong."** NCDR measures vulnerability to specimen noise, not the
false-discovery proportion. In cross-species draws, 71% of whole-PT GSEA species calls were also
called with all donors (base rate 3%), against 81% for the step screen (base rate 2%). A strong true
effect dominates a competitive ranking, so relabeled counts overstate how many species calls are noise.

**Why the screens differ.** Specimen offsets were far more coherent within pathways than
specimen-by-segment deviations: the median correlation of member genes was 0.26–0.50 for offsets,
against 0.03–0.20 for positional deviations. Variance inflation estimated from within-unit residuals
cannot absorb offsets this coherent.

In our cohort, within-species positional deviations were small and incoherent (0.013). This is lower
than in the external human pools (0.14–0.20), because our human pair is two sections of one donor. Our
relabelings therefore under-represent donor-level positional variation, and the external nulls
compensate only for the step screen.

**When positional screens lose specificity.**
- *Age* (post hoc). Male mouse splits unbalanced for pre-pubertal (3-week-old) mice gave median 22 and
  62 step-screen calls with an imbalance of one or two such mice, against 1 when balanced. A 2 + 2 design cannot
  separate this real pre- versus post-pubertal positional difference from the group factor.
- *Injury severity* (positive control; Figure 4e). In two AKI against two control mouse kidneys:
  - *Detection passed.* 13 of 14 injury markers rose most in S3; the positional screens called the
    injury Hallmarks TNF-α/NF-κB signalling, hypoxia and epithelial–mesenchymal transition; and
    77–80% of their calls changed most in S3.
  - *Specificity under two rules* (Supplementary Note 2).
    - The two AKI kidneys differ in positional injury: within-group specimen × segment coherence is
      0.42.
    - Even so, the coordinate screen passed the R4 rule (NCDR 0.17; 10 and 7 relabeled calls of 1,429
      tested), as did the average-level joint tests (NCDR 0.09–0.10).
    - These screens failed only the stricter positive-control rule specified before the run, which
      allows relabeled calls of at most 5% of condition calls; every strategy failed that rule.
    - The R4 rule did catch the step screen (26 calls against 12 and 6; NCDR 0.35).
  - *Reading.* The R4 rule cannot detect a positional nuisance of this size in the coordinate screen.
    Our cohort's specificity rests on its margin (NCDR 0.02; within-species coherence 0.013), not on
    passing the rule.

**Conclusion.** With two specimens per species, average-level pathway screens cannot be shown specific
by relabeling, in our cohort or in external designs. Position-dependent screens can be, when the
specimens within each group share their positional biology. That holds for adult human cortex,
age-matched mice and, as far as relabeling can show, our two control mice. Whether the two control
mice are age-matched is not recorded [VERIFY].

### R5 · Position-dependent pathways mostly summarise mouse zonation

**The primary list: pathways robust under two coordinates** (Figure 5a; Supplementary Figure S2;
Table S2).
- *Coordinates.* The pathway screen was run on two coordinates: the principal curve and a
  diffusion-pseudotime coordinate built on the same embedding.
- *Funnel on the principal curve.* Of 1,513 tested pathways, 60 were candidates, 58 not
  specimen-sensitive, 51 stable across planned refits and 36 robust to anatomical registration and to
  removal of deep mouse S3.
- *Diffusion-pseudotime coordinate.* It gave 26 robust pathways. All 36 principal-curve pathways are
  shown in Supplementary Figure S4.
- *Primary list.* The 22 pathways robust under both coordinates. Each coordinate is reported separately
  as a sensitivity analysis (Supplementary Note 1). The direction of the species difference at the
  peak agreed between coordinates for every pathway robust under either (40 of 40).

**What the 22 pathways contain.** Each pathway's contributing members (top 10% by T_spatial) were
classed as mouse-led, human-led or shared by their zonation classes (R3).

| Category | Pathways |
|---|---:|
| Summary of mouse zonation (≥ 50% of zonated members mouse-led) | 13 |
| Shared (Nitrogen metabolism, Aspirin ADME, Glutathione metabolism) | 3 |
| Mixed (Drug ADME) | 1 |
| Unresolved (fewer than 3 zonated members) | 5 |
| Human-led | 0 |

- *The mouse summaries* include fatty-acid metabolism, peroxisome, bile acid, cholesterol, vitamin and
  retinoid metabolism, phase II conjugation and the estrogen-response sets.
- *Core.* 12 of the 22 are also robust to restricting to genes detected in both species and to removing
  mouse sex-biased genes, under both coordinates.
- *Direction.* At the peak, 10 pathways are higher in human and 12 higher in mouse.
- *Peaks* (Figure 5d; segment units: 0–1 is S1, 1–2 is S2, 2–3 is S3).
  - Human-high pathways peak at a median of 0.50 on the principal curve, 0.95 on diffusion pseudotime
    and 0.54 on an equal-depth refit of the curve; that is, in S1 under every coordinate.
  - Mouse-high pathways peak at a median of 1.77–1.95, in late S2 near the S2→S3 transition.
  - Individual peaks move between coordinates (Spearman 0.71 between the two primary coordinates; 7
    of 22 change segment), so we read them only at segment resolution.
- *Why mouse-led.* A species × position statistic cannot recover a conserved core. A conserved gene
  enters T_spatial only through an amplitude difference, and the interaction is about −0.8 to −0.9 times
  the mouse gradient. The list is best read as "programs whose mouse zonation the human sections do not
  share".

**Coordinate-free support.** The calls might rest on within-segment order, which the coordinate
resolves weakly (R1). Two step-screen designs test this.
- *Coordinate-free step screen.* It is the screen that the external nulls calibrate, and on our cohort
  it calls 50 pathways with 0 and 0 under relabeling (R4).
  - It contains 20 of the 22 primary pathways, including all 12 core pathways.
  - The two exceptions, Formation Of Cornified Envelope and Metabolism of xenobiotics by cytochrome
    P450, are both probe-sensitive (below).
  - It also contains 28 of the 36 pathways robust on the principal curve.
- *Spline-based step design.* It calls 8 of the 22 (9 of the 36).
- *Position beyond the steps.* All 22 are also called by this statistic.

The primary list is therefore supported by an externally calibrated screen that uses no coordinate.

**Robustness** (Figure 5b).
- *Four checks.* All 22 keep at least 70% of their effect, with q ≤ 0.10, under anatomical
  registration, removal of deep mouse S3, removal of ambient and shared-exon genes, and finer
  covariate matching.
- *Genes detected in both species:* 20 keep the effect and 14 keep significance.
- *Removal of mouse sex-biased genes:* 21 and 20.
- *Equal probe counts:* all 22 keep the effect, but 7 lose significance: Glutathione metabolism, Aspirin
  ADME, Formation of Cornified Envelope, Steroid hormone biosynthesis, Retinol metabolism, Metabolism of
  xenobiotics by cytochrome P450 and Chemical carcinogenesis.

**Agreement with independent data is at the level of measurements** (Figure 5c; Supplementary Note 5).
- *Genes.* Our 10% most position-dependent genes agreed with the species × segment contrasts of two
  human cortex atlases (Spearman 0.67 and 0.73). That is above expression-matched random genes (0.57
  and 0.60), but no better than the most highly expressed genes (0.64 and 0.72).
- *Pathways.* Of the 22, 12 and 5 replicated against the two atlases with male mice, and 5 and 3 with
  female mice (BH ≤ 0.10).
- *Overlap-preserving collections* (tested on the 36 principal-curve pathways). Against random
  collections that keep set sizes, strata and overlaps, the robust pathways replicated in both atlases
  (mean z 1.53 and 1.24; p = 10⁻⁴). They also beat pathways called only under relabeling (p = 0.005 and
  0.015). They did not beat them with programs as units (p = 0.06 and 0.17), nor once the shared
  flattening of human gradients was removed (p = 0.42).

By the pre-specified rule this is replication better than chance, not pathway-specific replication.
The pathways replicate mainly because they are where mouse zonation is strongest.

### R6 · Lead genes

Each gene is given with its zonation class (R3) and, separately, its T_spatial rank, which can be high
from an amplitude difference alone.
- *Checks.* They were specified before the external data were opened. They used segment labels against
  the cortex-only human atlas (Lake et al. 2023; 7 donors), mouse snRNA by sex (12 male and 12 female
  donors) and microdissected mouse segments of both sexes (Chen et al. 2021, 2023).
- *Rodent protein.* Rat segment protein (Limbutara et al. 2020) is quoted as S3/S2 ratios, because rat
  S1 samples carry less protein (Supplementary Note 6).
- *Where the rest is.* The verification matrix is in Supplementary Figure S10, and further genes are in
  Supplementary Note 9.

Public Visium data from additional human donors (Abedini et al. 2024) could not be used. Their
authors' segment labels failed a pre-specified marker check: the S3 markers SLC5A1, SLC7A13 and SLC22A7
did not rise in their S3 spots (Supplementary Note 6).

Standing caveat: one human donor, two male mice, human cortex only, different probe panels.

**Creatine synthesis: a large amplitude difference** (Gatm indeterminate; Figure 6).
- *Mouse arm.* Gatm (AGAT, the first step) falls steeply from S1 to S3 in mouse: by 6.8 and 6.3 log2
  units in male and female snRNA, and by 12.1 and 11.0 in male and female microdissection. The
  cortical medullary-ray segment has 23.5 TPM, against 4,270 in S1.
- *Rodent protein.* Rat AGAT protein falls at S3 (S3/S2 −3.7). Transamidinase activity is confined to
  S1 and S2 in rat (Takeda et al. 1992).
- *Human arm.* GATM changes little and not consistently: −0.78 in our sections, +0.51 in the cortex
  atlas and −0.10 in a second cortex atlas (convoluted PT against S3). It is neither confidently
  zonated nor confidently flat.
- *Ranking.* Gatm ranks in the top 0.1% of genes by T_spatial.
- *Claim.* Mouse confines the first step of creatine synthesis to early PT, whereas human GATM is
  expressed along the PT. This is an amplitude difference, not a reversal.
- *Second step.* In rodents the second step, GAMT, rises toward late PT: in mouse microdissection and in
  rat protein (S3/S2 +1.6). Our data do not resolve Gamt.
- *Literature.* GATM is a human PT enzyme (Reichold et al. 2018), and the human kidney releases
  guanidinoacetate in vivo (Edison et al. 2007).

**Mitochondrial β-oxidation** (Acadm conserved; Acaa2 reversal, label-dependent; Figure 6).
- *Acadm (MCAD)* rises toward late PT in both species, much more in mouse (snRNA +2.1 and +2.3) than in
  the human atlas (+0.46).
- *Acaa2* (mitochondrial 3-ketoacyl-CoA thiolase) is confined to early PT in mouse (snRNA −4.0 and
  −4.2) and rises modestly in the human atlas (+0.58; +0.95 in male donors).
  - With reviewed labels, the reversal rests on S2 − S1, the contrast with the weakest human labels;
    the region-matched S3c − early is indeterminate.
  - With human labels transferred from the atlas on all genes, Acaa2 is a reversal in S3c − early
    instead.
- *Rodent protein* matches the mouse directions: Acadm S3/S2 +0.7 and Acaa2 −2.7.
- *State.* Neither gene is sex-biased in mouse PT, and their mouse zonation does not change with time
  of day. Whether fasting flattens it in human cannot be tested with public data.
- *Ranking.* Both rank in the top 1% by T_spatial.

**Genes whose difference changes sign** (Dcxr, Ugt3a1 reversal; Figure 6).
- *Dcxr* (dicarbonyl/L-xylulose reductase) falls toward late PT in mouse and rises in human. It is a
  reversal in S3c − early and S3 − early.
  - Human: +0.97 in our sections and +1.53 in the cortex atlas (+1.59 in men, +1.48 in women).
  - Mouse: −0.46 and −0.54 in snRNA, −1.94 and −2.03 in microdissection. Rat protein agrees (S3/S2
    −1.1).
  - Probes are balanced. DCXR protein is known in mouse PT (Nakagawa et al. 2002), and DCXR loss tracks
    human CKD (Perco et al. 2019); its axial distribution is not reported in either species.
- *Ugt3a1* rises toward late PT in mouse (snRNA +1.0 and +1.5) and falls in human (ours −0.62, atlas
  −0.98). It is a reversal in S3c − early. Rat protein does not quantify it.
- *Ranking.* Both rank in the top 1% by T_spatial.

**Expressed and zonated only in human** (RBP4; Figure 6).
- *Human.* RBP4 rises steeply toward late human PT (ours +2.46, atlas +3.46). Mouse Rbp4 is absent:
  0.0 TPM in every microdissected segment of both sexes. Because it is expressed in one species only,
  RBP4 is outside the zonation classes.
- *Rat.* RBP protein is found only in S1, where it is reabsorbed ligand. Rat in situ hybridisation
  placed RBP mRNA in outer-stripe S3 and also in perinephric fat (Makover et al. 1989).
- *Alternatives.* This is an mRNA claim. Ambient RNA from perirenal fat, transcript diffusion between
  bins and, for protein, uptake of filtered RBP4 are untested alternatives.

**A conserved reference** (Slc7a13; Figure 6). Slc7a13 rises toward late PT in both species and is
shown for comparison.

**Supplementary gene panels** (Supplementary Note 9; Supplementary Figure S11).
- *Glutathione synthesis.* Gclc and Gclm are indeterminate and Gss is conserved. Human GCLM falls in the
  cortex atlas where mouse rises.
- *Cyp2e1.* It is expressed in mouse PT and essentially not in human: an expression difference, not a
  positional one.
- *Also covered there:* sterol synthesis, organic anion uptake, early-PT transporters in late PT,
  UGT1A9, AOX1, and genes whose claims did not hold.

---

## 5. Discussion

**Main finding.** Measured against how well two datasets of the same species agree, human and mouse PT
agree little in how their genes change along the S1→S3 axis. The conservation index was 0.17–0.29 in
our data and 0.27–0.28 in public single-nucleus atlases alone. It stayed at most 0.33 when the human
segment labels were transferred from a reference atlas on held-out genes, and it fell slightly when
the genes that built the integration were removed.
- *Where the agreement is.* It sits in a small set of genes. Among genes zonated in both species,
  72% agree in direction, and transporters are over-represented among them. [PENDING WS6: whether
  agreement rises with gradient strength, and the index range under winsorisation and block or donor
  bootstraps. If agreement is confined to the most strongly zonated genes, state it here as the main
  result.]
- *What it implies.* Outside that set, mouse PT zonation is a weak guide to human PT zonation.

**What our data add to public atlases.** The central comparison can be reproduced from public atlases
alone, which supports it. Our data add three things:
- a measurement of both species in intact tissue on one platform;
- region matching by physical depth, which restricts the mouse S3 comparison to cortex;
- replication of measurements across structures and positions within each specimen, which the
  position-dependent pathway screen uses.

**Absolute or relative zonation.**
- *Two scales.* Human gradients were smaller than mouse gradients (spread 0.58 of mouse for S1→S2).
  Whether more genes are zonated only in mouse therefore depends on the scale. On an absolute scale,
  mouse-only genes outnumber human-only ones. Relative to each species' amplitude, they do not.
- *What would decide it.* The choice turns on whether the lower human amplitude is technical (capture,
  bin coverage, label impurity) or biological. Our data cannot settle it. The human amplitude was
  similar in an independent single-nucleus atlas (0.47), which makes a cause specific to our platform
  less likely, but not causes shared with that atlas, such as procurement or segment definition.
- *What does not depend on it.* In either case, confidently species-only genes are few, and the low
  agreement is a genome-wide property rather than a long list of species-specific genes.

**An evolutionary benchmark is missing.** Rat microdissection data were too noisy to serve as a third
species, so we cannot say whether agreement this low is unusual for two mammalian species. PT zonation programs may be evolutionarily labile in general. A rodent pair measured
with replicates would answer this.

**Pathway specificity with few specimens.** The second finding is methodological and general.
- *Average-level screens.* With two specimens per species, pathway screens that compare averages
  cannot be shown specific by relabeling. They called as many pathways for mixed-species relabelings as
  for species, in our cohort and in external designs.
- *Not wrong.* Most of their species calls recurred when all donors were used (71%), so the problem is
  that specificity cannot be demonstrated, not that the calls are false.
- *Position-dependent screens* passed the relabeling control when the specimens within each group
  shared their positional biology. They failed when the specimens differed in age or injury severity.
- *Why.* Specimen offsets move whole pathways together, whereas positional deviations are small and
  incoherent in healthy, matched specimens.
- *Limits of the control.* Our specificity rule did not detect the positional nonspecificity in the
  AKI design. A relabeling control therefore needs the external calibration and the margin we report,
  not only a pass.
- *Applies widely.* The argument holds for any spatial or single-cell comparison with a few specimens
  per group in which position or state varies within specimens.

**Position-dependent pathways.**
- *What the 22 robust pathways summarise.* Mostly mouse zonation that the human sections do not share.
  None is led by human zonation.
- *Coordinate-free support.* Twenty of the 22 are also called by a step screen that uses no coordinate
  and is calibrated by the external nulls.
- *Expected from the statistic.* A species × position statistic cannot recover a conserved component,
  and the human–mouse interaction is dominated by a global reduction of the mouse gradient.
- *What the list is useful for.* Naming the mouse programs whose positional organisation does not
  carry over (fatty-acid and peroxisomal metabolism, bile acid and cholesterol handling, vitamin and
  retinoid metabolism). It does not show where human PT places these programs.

**Relation to prior work.**
- *Human atlases* label PT-S1/S2/S3 as discrete classes of dissociated nuclei (Lake et al. 2023).
  Mouse microdissection gives three anatomical points (Chen et al. 2021, 2023). Rat segment proteomics
  gives protein along the rodent tubule (Limbutara et al. 2020).
- *Cross-species atlases* align cell states (Klötzer et al. 2025). They find little overlap of
  single-cell expression changes but agreement at the pathway level (Zhou et al. 2023).
- *Cross-species integration* trades species mixing against conservation of biology (Song et al. 2023).
  That is one reason we validated the labels without the integrated embedding.
- *Prior comparisons.* A structured literature search found no position-by-position comparison of
  human and mouse PT gene programs. The PT annotation depth of Klötzer et al. (2025) could not be
  checked, because the article is not open access.

**Implications for mouse models of PT physiology and toxicity.**
- *Transport.* For transporters zonated in both species, mouse predicts the direction of human
  zonation more often than for other genes (83% against 72%). Conserved zonation is not conserved usage, however:
  SLC34A3 carries a large share of phosphate transport in human but not in mouse (Fernandes et al.
  2026). OAT1 and OAT2 transcript profiles also disagree with human protein maps (Breljak et al.
  2016).
- *Metabolism.* Several rodent metabolic placements have weaker counterparts in the human sections,
  including early-PT creatine synthesis and mitochondrial β-oxidation steps. An effect confined to one
  mouse segment may be weaker, or spread along the human PT.
- *Expression.* Some differences are differences of expression, not of position. Mouse PT expresses
  Cyp2e1, which human PT barely expresses, and human PT expresses RBP4, which mouse PT does not.
- *Physiological state.* Mouse PPARα and sterol programs respond strongly to fasting (Aomura et al.
  2026), and PPARα targets differ between mouse and human (Rakhshandehroo et al. 2009). The renal
  tubular clock controls β-oxidation (Bignon et al. 2023). Laboratory mice are fed ad libitum, so part
  of any species difference in these programs may be physiological state.
- *Injury models.* Experimental injury models already differ in where injury falls along the nephron
  (Heyman et al. 2010). Our results add that many genes are zonated differently.
- *Computational models.* Human nephron models borrow rodent axial parameters (Layton et al. 2019).
  Strong transporters support that choice more than other genes do.
- *All of these are mRNA-level observations.*

**Limitations.**
- *Replication.* Two male control mice and two cortex sections from one male human donor. Specimens,
  not structures, are the replicates. Agreement between the two human sections is not donor
  replication, and no statement here is a population estimate for humans. The ages of the two mice are
  not recorded [VERIFY]. If they differ, the specificity of the position-dependent screen in our cohort
  (R4) must be revisited.
- *Sex.* All specimens are male, so species and sex are confounded. Female mouse references gave
  similar class counts.
- *Region.* The human sections are cortex only, so mouse outer-stripe S3 has no counterpart. We matched
  the region physically for S3 (S3c); the reference atlases' S3 lies mostly in the outer stripe.
- *Labels.* The segment labels are cross-species cluster memberships with weak reference support in
  human, especially for S2. The conservation index is robust to transferred labels. Gene claims that
  rest on S2 − S1 alone, such as the Acaa2 reversal and several species-only calls, are
  label-dependent.
- *Classification.* Most genes (6,307 of 7,407) cannot be classified with two specimens per species.
  The human confirmation atlas has fewer donors than the mouse one.
- *Measurement.* The species were measured with different probe panels and Space Ranger versions.
  Genes with unequal probes are reported separately.
- *Coordinate.* Its within-segment order is weakly supported and depth-dependent in human S1. Pathway
  calls depend on the coordinate, so the primary list is the set robust under two coordinates.
- *External references.* External human positions are single-nucleus cluster annotations, not
  anatomy. Pathway-level replication beyond the shared global component is not shown.

**Shared confounders.** Our human tissue and the external human atlases share surgical or biopsy
procurement with warm and cold ischaemia, older donors, comorbidity and medication, and
pre-operative fasting. Mouse kidneys came from young, ad-libitum-fed animals. Fasting changes
PPARα-target and sterol-synthesis genes in mouse kidney by amounts comparable to the species
differences we report, and human PT expresses the fasting- and stress-responsive gene PDK4
strongly. Segment-resolved data from fasted animals or fed humans do not exist. A contribution of
physiological state to the species difference, especially for sterol synthesis and PPARα
targets, is therefore not excluded. Caloric restriction also feminizes male mouse kidney gene
expression (Xiong et al. 2023), and PPARα regulates largely different target genes in mouse and
human hepatocytes (Rakhshandehroo et al. 2009).

The source, age and procurement of our donor are not recorded here [VERIFY]. Deviations from the
pre-specified protocols are listed in Supplementary Note 8.

**What would move this forward.**
- Several human donors with spatial data at single-cell resolution and anatomical segment annotation.
- Female animals and female human donors.
- Fasted-mouse or fed-human references.
- A replicated second rodent species.
- In situ validation on human sections where S1 and S3 can be distinguished: GATM, ACAA2, DCXR,
  UGT3A1 and RBP4.

---

## 6. Methods

Parameters, software versions and the analysis-to-code map are in the Supplementary Methods.

### Specimens and ethics

**Specimens.**
- *Mouse.* Two control kidney sections from two mice (`Ctrl1A2`, `Ctrl1A4`). These are whole-kidney
  sections and include cortex and outer medulla. Strain, age, supplier, housing and euthanasia are not
  yet recorded [VERIFY]. Two ischaemia–reperfusion kidneys from the same experiment were used only as
  a positive control and for tissue-state calibration.
- *Human.* Two sections of renal cortex from one donor (`HUK1_COR1`, `HUK1_MED1`). The second is named
  "MED" at source, but its segmentation contains cortex only. Tissue source, donor age, kidney function,
  procurement, and whether the sections are serial are not yet recorded [VERIFY].
- *Sex.* All four specimens express Y-linked transcripts and are treated as male [VERIFY against
  records].

**Ethics.** Animal-protocol and human-tissue approvals [VERIFY].

### Visium HD and segmentation

**Visium HD.** Sections were profiled with Visium HD (Oliveira et al. 2025), and the 2-µm bin output
was used.

| | Mouse | Human |
|---|---|---|
| Space Ranger | 3.1.1 | 4.0.1 |
| Probe set | Mouse Transcriptome v2.0 | Human Transcriptome v2.1.0 |
| Reference | mm10 | GRCh38-2024-A |

Library preparation and sequencing depth are not yet recorded [VERIFY].

**Segmentation** (Supplementary Note 7).
- *Model.* Tubule cross-sections and glomeruli were segmented on H&E whole-slide images with a panoptic
  model at 0.44068 µm per pixel. It combines a frozen pathology foundation encoder (Oquab et al. 2023;
  Kaplan et al. 2025) with a trainable decoder, and decodes instances by marker-controlled watershed.
  Inference used D4 test-time augmentation.
- *Training and selection.* The model was trained on partially annotated mouse slides, and the
  checkpoint was selected by panoptic quality (Kirillov et al. 2019).
- *Mouse polygons* passed an upstream quality-control step whose rules are not yet recorded [VERIFY].
- *Human polygons* come from a re-segmentation whose provenance and pixel scale are not yet recorded.
  Their µm fields assume 0.5 µm per pixel [VERIFY].

**Tubule-by-gene matrices.**
- *Assignment.* Each bin was assigned to the polygon covering its centre, in full-resolution image
  coordinates, and raw counts were summed per polygon. Every join recomputes each polygon's centroid
  and stops at any deviation above 1 px.
- *Coverage.* Bins covered a median 97.4% of mouse polygon area and 76.5% of human polygon area. The
  human value is close to the 75.7% expected if the human µm fields overstate the pixel size.
- *Polygons without bins.* Human polygons without bins lie outside the area covered by bins: 0.2% and
  1.2% of the polygons inside the bin-covered hull lack bins. The method used to register the H&E
  images to the Visium HD frame is not yet recorded [VERIFY].
- *Structure filters.* Human structures needed ≥ 100 detected ortholog genes. Within each species,
  structures below the 5th percentile of bins per structure were removed. Mouse polygons had no
  gene-count threshold, because they had passed upstream quality control. 26,839 of 37,559 structures
  were retained.

### Ortholog space, integration and segment labels

**Ortholog space.** Human counts were mapped to mouse symbols through mutual-best HCOP pairs supported
by at least three databases (Yates et al. 2021; 17,449 pairs) [VERIFY download date].
- An ortholog column that no human feature feeds is a structural zero, not an observed zero.
- Gene-level analyses use the 15,567 orthologs measured in all four specimens.

**Integration and clustering.**
- *First pass.* Highly variable genes were selected within each species and intersected (732 genes).
  The first 50 principal components were integrated with Harmony (Korsunsky et al. 2019; θ = 6, λ = 1,
  batch = specimen; R harmony 2.0.5) [VERIFY θ rationale]. Leiden clusters (Traag et al. 2019;
  resolution 0.7) were labelled by marker review.
- *Second pass.* After removal of glomerulus, smooth-muscle and unresolved clusters, integration was
  repeated on the 24,335 tubular structures, and 12,866 PT structures were taken.
- *Caveat.* Specimen and species are confounded in the integration, so the integrated embedding is used
  only to group and order structures. Every gene comparison uses raw counts.

**Segment labels.** S1, S2 and S3 are memberships of the joint PT clusters. Marker panels:
- S1: Slc5a2, Slc5a12, Gatm, Lrp2, Cubn, Slc34a1;
- S2: Slc22a6, Slc13a3, Cyp2e1;
- S3: Slc22a7, Slc7a13, Cyp7b1, Slc6a18, Acsm3.

**Label validation.** The protocol was specified before any result.
- *Glomerular distance.* Distance from each PT centroid to the nearest glomerulus polygon in the same
  section, and the Spearman correlation of label order with distance.
  - A pre-registered size rule required glomerulus polygons to be at least 3× the median PT area. Mouse
    polygons failed it, so mouse distances are descriptive.
  - Human distances are converted at 0.44068 µm per pixel, a factor of 0.881 on the recorded 0.5 µm
    [VERIFY against scanner metadata].
- *Attachment.* The share of S1 among PT polygons within 2 µm (primary) or 10 µm of a glomerulus.
- *Reference transfer.* Each structure was assigned the reference segment whose profile correlated best
  with its own (Spearman), on genes that differ by ≥ 1 log2 between reference segments. References
  were healthy human cortex single-nucleus data (Lake et al. 2023), mouse microdissection (Chen et al.
  2021) and male mouse single-nucleus data. Support for the human labels required agreement ≥ 0.70 and
  within 0.10 of mouse.
- *Spillover.* The share of counts from fixed non-PT markers in each PT structure.

### Segment gradients, ceilings and the conservation index

**Gradients.**
- *Definition.* Pseudobulk log2 CPM ratios of summed counts, per specimen or donor, for S2 − S1,
  S3 − early (early = S1 + S2) and S3c − early.
- *S3c.* Mouse S3 structures no deeper than the 95th percentile of same-section S1 depth, where depth is
  the mean distance to the three nearest glomeruli.
- *Genes.* In our data, ≥ 20 counts in every relevant specimen × segment pseudobulk. Externally, a mean
  log2 ≥ 3 in one compared segment group.

**Datasets.**
- *Human:* our sections, healthy cortex single-nucleus data (Lake et al. 2023; 7 donors with ≥ 50
  nuclei per segment), and a second cortex atlas (Acera-Mateos et al. 2026; 6 donors; convoluted PT
  against S3 only).
- *Mouse:* our sections, single-nucleus data by sex (Chen S et al. 2025; 12 donors each, via CELLxGENE
  Census; CZI Cell Science Program et al. 2025), and microdissection (Chen et al. 2021, 2023).

**Correlations and the index.**
- *Correlations.* Gene gradients were winsorised at the 1st and 99th percentiles and correlated
  (Pearson). Correlations were disattenuated by each dataset's replicate reliability (signal ÷ total
  variance across genes).
- *Index.* r*(our human, our mouse) ÷ √(r*(our human, Lake) × r*(our mouse, male mouse snRNA)),
  computed on genes expressed in all four datasets, with 95% intervals from 1,000 gene bootstraps.
- *External-only version.* The numerator is r*(Lake, male mouse snRNA).
- *Decision rule (pre-specified).* An upper bound < 0.5 in both primary contrasts gives "largely
  different genes"; 0.5–0.8 gives "partly conserved"; ≥ 0.8 gives "largely conserved".
- *Other species.* A further species was used only with replicate reliability ≥ 0.5 and ≥ 70% symbol
  matching.
- *Integration genes.* The index and classes were also evaluated without the 611 genes that built the
  integrated embedding. The moderation prior still used all genes.

[PENDING WS6: winsorisation sensitivity with reliability on matched values, block and donor
bootstraps, and agreement by gradient strength with per-bin ceilings.]

**Cross-fitted labels.**
- *Splitting.* Genes were split into halves by a hash of the gene symbol.
- *Transfer and evaluation.* Human labels (and, in one setting, mouse labels) were transferred using
  reference-zonated genes of one half. The index and class counts were evaluated on the other half, and
  then the halves were swapped.
- *Criteria (pre-specified).* Index ≤ 0.5 with an upper bound < 0.6 in every setting, half and primary
  contrast; and mouse-only ÷ human-only < 2.

**Amplitude.** The noise-corrected SD ratio of human to mouse gradients. Signal variance is the
variance of the mean gradient minus the mean per-gene sampling variance. Regression of each species on
a reference gradient gives a projection ratio, reported in Supplementary Note 4.

### Zonation classes

**Per-species calls.**
- *Test.* The specimen-mean contrast was tested with a moderated one-sample t. Variances were shrunk
  within ten expression bins, as in limma's `squeezeVar` (Smyth 2004), and BH was applied across genes.
- *Zonated.* q ≤ 0.05 and |contrast| ≥ 0.5 log2, confirmed by an independent same-species reference
  (one-sided donor-level p ≤ 0.05). The references were 12 male mouse snRNA donors and 7 healthy human
  cortex donors.
- *Flat.* The 95% interval lies within ±0.3 log2.
- *Primary (symmetric) setting.*
  - The human effect floor and flat margin were multiplied by the noise-corrected human ÷ mouse
    amplitude ratio of the contrast (0.58, 0.84 and 0.70 for S2 − S1, S3c − early and S3 − early).
  - A species-only call also required the other species' reference to be flat-compatible: one-sided
    p > 0.05 in the same direction, and |mean| < 0.5. Otherwise the gene was indeterminate.

**Classes.** Conserved, reversal, mouse-only, human-only, indeterminate or neither, per contrast and as
a union over S2 − S1 and S3c − early. Class enrichment used a competitive rank test against
expression- and detection-matched random sets, with CAMERA-style variance inflation. An enrichment was
named only with q ≤ 0.10, robustness to expression matching, and at least three class genes carrying
it. Primary counts use genes with equal probe counts in both panels.
Alternative settings (unscaled thresholds; no reference confirmation of flat calls; female mouse
references) are in Supplementary Note 3.

### PT coordinate and pathway screens

**Coordinate** (Supplementary Note 1).
- *Principal curve.* A nonbranching principal curve (scFates; Faure et al. 2023) on the first five
  Harmony dimensions of the PT structures, with 30 nodes. It was rooted at the tip with the higher
  early-marker score and scaled to [0, 1].
- *Diffusion pseudotime.* A diffusion-pseudotime coordinate (Haghverdi et al. 2016) was built on the
  same dimensions.
- *Common support.* Gene models use the intersection of each specimen's 1st–99th percentile interval
  (12,272 structures).

**Gene models.**
- *Models.* For each of 11,071 eligible genes (log1p of counts per 10⁴, mean detection ≥ 2%), nested
  weighted regression-spline models were fitted:
  - M₀: shared shape plus specimen contrasts;
  - M_level: adds a species offset;
  - M_full: adds a species × position spline.
- *Statistic.* T_spatial is the partial F of M_full against M_level. It ranks genes and is not a test,
  because its degrees of freedom count structures.
- *Step screens.* Two designs test species × S1/S2/S3 steps beyond a species offset. The spline-based
  design adds the steps to the shared spline along the coordinate. The coordinate-free design uses
  segment intercepts as the shared shape, as in the external nulls, with exact matched moments.
- *Weights.* Each species carries half the total weight, shared equally between its specimens.

**Joint pathway test.**
- *Libraries.* Reactome 2022, MSigDB Hallmark 2020 and KEGG 2019 Mouse (Gillespie et al. 2022;
  Liberzon et al. 2015; Kanehisa et al. 2019) in Enrichr format (Kuleshov et al. 2016) [VERIFY
  download source and date]. Sets with 10–300 eligible members were tested: 1,513 pathways.
- *Matched test.* The rank-AUC of members' statistics was compared with 9,999 random sets drawn within
  strata of expression, detection and positional coverage.
- *Joint test.* The matched z was divided by a CAMERA-style variance inflation (Wu & Smyth 2012)
  estimated from within-specimen residual correlations. BH was applied over all pathways.
- *Peaks.* The peak is the position of the largest member divergence. It is reported in segment units:
  the segment index plus the fraction through the segment, using the run's mean segment transitions
  (0–1 is S1, 1–2 is S2, 2–3 is S3).
- *Beyond-steps and step screens.* Each pathway was also tested with two statistics from the same
  models: position beyond the S1/S2/S3 steps (the species × spline term given species × segment
  steps), and the step screen.

**Reporting set.**
- *Candidates.* Joint q ≤ 0.05 with a positive effect.
- *Not specimen-sensitive.* Not called under either relabeling or either within-species comparison.
- *Stable.* Retained in at least 5 of 7 planned refits.
- *Robust.* Effect ratio ≥ 0.7 and q ≤ 0.10 after two checks:
  - piecewise-linear registration of human segment transitions onto mouse;
  - removal of mouse S3 deeper than the S1 range, compared with power-matched random removal.
- *Primary list.* Pathways robust under both the principal curve and the diffusion-pseudotime
  coordinate. The pathway pipeline was run unchanged on each coordinate.
- *Categories.* Contributing members (top 10% by T_spatial) were classed as mouse-led, human-led or
  shared from the zonation classes. A pathway takes a category when ≥ 50% of its zonated members fall
  in it.

### Specificity of pathway screens

**Relabeling.**
- *Contrasts.* The four specimens have three balanced 2 + 2 contrasts: the species split and two
  mixed-species relabelings. Every strategy was rerun unchanged on both relabelings, as in Lamian's
  random-partition null (Hou et al. 2023). In relabeled models, species-specific shape enters as a
  nuisance term.
- *Motivation.* The design is a split plot (Altman & Krzywinski 2015): offsets have 2 df of error, and
  positional terms more. This motivates the comparison of offset and positional screens but does not
  guarantee specificity.
- *NCDR* is the mean relabeled count divided by the species count.
- *Rule.* A strategy is called specific when NCDR < 0.25 and each relabeling calls < 5% of the tested
  pathways. The second condition was added after the first run (Supplementary Note 8).

**Many-draw nulls.**
- *Pools.* PT nuclei with author segment labels from:
  - Census mouse male and female (12 donors each; embryonic and newborn stages removed);
  - Census human cortex (6 donors);
  - Lake healthy cortex (7 donors).

  Donors needed ≥ 50 nuclei in every segment. Tested genes, strata, pathways and variance inflation
  were fixed once per pool.
- *Strategies.* For each 2 + 2 partition:
  - unit-level offset and segment-step statistics, tested with the matched and joint tests;
  - whole-PT and per-segment pseudobulk moderated t statistics, tested with signed preranked GSEA
    (Subramanian et al. 2005; Fang et al. 2023).
- *Designs.* All same-species splits were used, sampled up to 300 per mouse pool, together with 100
  cross-species draws of 2 human + 2 mouse donors.
- *Precision.* The share of a draw's species calls also called with all donors.
- *Coherence.* The median within-pathway correlation of member genes, for specimen offsets and for
  specimen × segment deviations.

**AKI positive control.** Two ischaemia–reperfusion against two control mouse kidneys, on the
mouse-only PT pseudotime, with the same strategies and the coordinate screen.
- *Detection* passed if injury programs were called and changed most in S3.
- *Specificity* was judged by the R4 rule and by the pre-specified positive-control rule (relabeled
  calls ≤ 5% of condition calls).

### External agreement and lead genes

**External agreement.**
- *Genes.* The Spearman correlation of our species × segment interaction with the external
  interaction, over informative genes. Comparators: expression-matched random genes and the most highly
  expressed genes.
- *Pathways.* A donor-level Welch t per gene, followed by the matched rank-AUC per pathway, with BH
  within the robust list.
- *Replication with overlap kept.* Mean replication z of the robust pathways was compared with 9,999
  random collections. Each collection relabels genes within matching strata and keeps set sizes,
  strata and pairwise overlaps. The comparators were random sets and pathways called only under
  relabeling, each at the pathway and the program level.

**Lead genes.** The checks were specified before the external data were opened:
- our arms in each specimen;
- robustness to attenuation;
- the human atlas arm;
- the interaction sign against male and female mouse snRNA;
- mouse arms in snRNA and microdissection of both sexes;
- absence for one-species genes.

Probe counts per gene came from both panels. Human protein data record presence only. Rat segment
protein (Limbutara et al. 2020) is quoted as log2 S3/S2 ratios, centred on the median ratio of proteins
quantified in all PT samples.

### Statistics

- *What the statistics mean.* Two mice and two sections of one human donor give no estimate of
  variation among human donors. Results are conditional on these specimens and descriptive.
- *Pathway tests.* Competitive tests of a gene statistic against matched genes, not population tests
  of a species effect.
- *Multiple testing.* BH at q ≤ 0.05 within each family: each statistic across pathways, each partition
  and strategy, genes within each species for zonation calls, and the robust list for external
  replication.
- *Intervals.* Gene bootstraps treat genes as independent and condition on the specimens.
- *Display rules.* Peaks, categories and programs are descriptions, not tests.

### Data and code availability

[VERIFY: deposition of the Visium HD data and segmentation polygons (controlled access for human
tissue), code archive DOI and licence.]

Every result is produced by committed, output-free notebooks that take the data and results roots as
arguments. The Supplementary Methods map each analysis to its notebook. The public datasets are cited
above.

---

## 7. Figure legends and tables

The source table for each panel is in `figure_plan_v04.md`. Every figure carries the standing caveat:
two control mice (M1, M2) and two cortex sections (H1, H2) from one male human donor; H2 is labelled
medulla but is cortex; all findings are descriptive.

**Figure 1 | Data and segment labels.**
- **a,** Design. Visium HD 2-µm bins from two control mouse kidney sections and two human cortex
  sections were summed into segmented tubule cross-sections: 26,839 structures, of which 12,866 are PT.
- **b,** Distance of PT structures to the nearest glomerulus by segment label: median and interquartile
  range. The Spearman correlation of label order with distance is 0.65 in mouse and 0.23 in human.
  - Human distances use 0.44 µm per pixel [VERIFY].
  - Mouse distances are descriptive, because the mouse glomerulus polygons fail the size rule.
- **c,** Share of S1 among PT structures within 10 µm of a glomerulus against all PT structures, per
  species (odds ratios 4.5 in mouse and 2.0 in human).
- **d,** Labels transferred by expression alone against our labels, as row shares. Human labels come
  from the cortex atlas and mouse labels from microdissection. Agreement is 71% (κ 0.54) in human and
  87% (κ 0.80) in mouse.
- **e,** The supporting coordinate: per-specimen segment medians (points) and segment transitions
  (ticks). The curve is fitted in five integrated dimensions.

**Figure 2 | Measured against same-species reproducibility, human–mouse agreement of PT gradients is
low.**
- **a,** *Left:* noise-corrected correlations of gene gradients for dataset pairs within human, within
  mouse and between species, for S2 − S1 and S3 − early. Each point is one pair; the within-human
  S2 − S1 value is a single pair. *Right:* the conservation index with its 95% interval, for our data
  and for public atlases alone. Dashed lines mark the pre-specified thresholds of 0.5 and 0.8.
- **b,** [PENDING WS6] Correlation and sign agreement of human and mouse gradients in bins of mouse
  gradient strength, against per-bin same-species ceilings.
- **c,** The index under cross-fitted transferred labels and without the 611 integration genes, per
  setting and gene half.
- **d,** Noise-corrected human ÷ mouse spread of gradients, for our data and public atlases.

**Figure 3 | Which genes are shared.**
- **a,** Per-gene S2 − S1 in mouse against human, coloured by zonation class (7,391 genes tested in
  S2 − S1; conserved 109, reversal 37, mouse-only 31, human-only 58).
- **b,** Among genes zonated in both species, the share agreeing in direction: transporters 33 of 40
  against all genes 196 of 271. Bars at left give the share of zonated genes that are zonated in both
  species (93% against 70%).
- **c,** Mouse-only and human-only counts under absolute or amplitude-scaled thresholds, with or
  without flat confirmation in the other species' atlas, and with cross-fitted transferred labels.
- **d,** [PENDING WS6] Class composition by mouse gradient strength.

**Figure 4 | With two specimens per species, average-level pathway screens cannot be shown specific by
relabeling.**
- **a,** The species split and the two balanced relabelings of four specimens, with split-plot error
  terms. These motivate the design and do not guarantee specificity.
- **b,** Pathways called (BH q ≤ 0.05) for species (filled) and each relabeling (open), with NCDR. The
  six strategies are whole-PT and segment GSEA, the offset joint test, the coordinate-free and
  spline-based step screens, and the coordinate screen.
- **c,** Pathways called per same-species 2 + 2 split in four public pools (300, 300, 105 and 45 splits).
  Grey marks average-level screens and black the step screen. Stars mark our cohort's relabelings, as
  re-implemented for this comparison.
- **d,** Cross-species draws (100). *Left:* median NCDR per strategy; the dashed line marks 0.25. *Right:*
  precision against the all-donor reference: whole-PT GSEA 71% and step screen 81%, with base rates of
  3% and 2%.
- **e,** AKI positive control: calls for AKI against control and for each relabeling, per strategy,
  marked as passing or failing the R4 rule and the pre-specified positive-control rule. The coordinate
  screen passes the R4 rule (NCDR 0.17) and fails the positive-control rule.

**Figure 5 | Position-dependent pathways mostly summarise mouse zonation.**
- **a,** The 22 pathways robust under both coordinates. Bars split each pathway's contributing members
  into mouse-led, shared, human-led, mixed, not zonated and not classifiable. ★ marks core pathways,
  and ◆ marks the 20 also called by the coordinate-free step screen.
- **b,** Effect retained under seven checks. A separate marker shows robustness under each coordinate.
- **c,** Our species × segment interaction against the external interaction (cortex atlas − male mouse
  snRNA). The 10% most position-dependent genes are black.
- **d,** Peak positions in segment units under the two coordinates, for human-high and mouse-high
  pathways.

**Figure 6 | Lead genes.** Gatm, Acadm, Acaa2, Dcxr, Ugt3a1 and RBP4, with Slc7a13 as a conserved
reference.
- *Left of each card:* fitted human and mouse curves along the coordinate, with each specimen's
  12-bin means.
- *Right:* each external dataset's S1, S2 and S3 values centred on its own S1. The datasets are the
  human cortex atlas, mouse microdissection and mouse snRNA (both sexes), and rat protein as S3/S2.
- *Titles:* zonation class, probe counts (mouse/human) and T_spatial rank.

**Table 1 | Conservation index and zonation-class counts.**
- *Part A:* the index for each contrast and version (our data, public atlases, alternative ceilings,
  cross-fitted labels, without integration genes), with 95% intervals and both ceilings.
- *Part B:* class counts for each rule setting and contrast (probe-balanced genes).

**Supplementary tables.**

| Table | Content |
|---|---|
| S1 | Specificity of every pathway strategy: species and relabeled calls, NCDR, null rate, decoy calls, minimum FDP |
| S2 | All 1,513 tested pathways: joint test, specificity margin, robustness ratios, stability, program, shape, external replication, and robust, core and peak under each coordinate |
| S3 | External replication per pathway |
| S4 | Programs, with drivers and member pathways |
| S5 | External replication summary and global component |
| S6 | Coordinate comparison and selection rule |
| S7 | Probe counts; zonation classes of every gene under every setting; genes expressed in one species only |
| S8 | Segmentation performance [VERIFY whether to report] |
| S9 | Many-draw nulls: call distributions, cross-species summary, AKI calls, coherence |
| S10 | Label validation and cross-fitted index and classes |
| S11 | Physiological and tissue state; additional Visium donors; rat segment protein |

---

## 8. References

The list covers the main text, the Supplementary Notes and the Supplementary Methods. Every entry was
checked in Europe PMC, Crossref or on the publisher page. Datasets without a citable publication are
identified by accession (GSE267280, GSE277302) [CITE: publications, if any].

- Abedini A et al. Single-cell multi-omic and spatial profiling of human kidneys implicates the fibrotic microenvironment in kidney disease progression. *Nat Genet* 56:1712–1724 (2024). https://pubmed.ncbi.nlm.nih.gov/39048792/
- Acera-Mateos M et al. Systematic evaluation of single-cell multimodal data integration enhances cell type resolution and discovery of clinically relevant states in complex tissues. *Genome Biol* 27:64 (2026). https://pubmed.ncbi.nlm.nih.gov/41821037/ (Census collection cites the preprint https://doi.org/10.1101/2025.03.06.637075.)
- Akakpo JY et al. Lack of mitochondrial Cyp2E1 drives acetaminophen-induced ER stress-mediated apoptosis in mouse and human kidneys. *Toxicology* 500:153692 (2023). https://pubmed.ncbi.nlm.nih.gov/38042273/ †
- Albergante L et al. Robust and Scalable Learning of Complex Intrinsic Dataset Geometry via ElPiGraph. *Entropy* 22:296 (2020). https://pubmed.ncbi.nlm.nih.gov/33286070/
- Altman N, Krzywinski M. Split plot design. *Nat Methods* 12:165–166 (2015). https://pubmed.ncbi.nlm.nih.gov/25879095/
- Aomura D et al. Proximal Tubule-Specific Genetic Deficiency of PPARα Worsens Systemic Lipid and Glucose Metabolism During Fasting. *FASEB J* 40:e71333 (2026). https://pubmed.ncbi.nlm.nih.gov/41460648/
- Basit A et al. Kidney Cortical Transporter Expression across Species Using Quantitative Proteomics. *Drug Metab Dispos* 47:802–808 (2019). https://pubmed.ncbi.nlm.nih.gov/31123036/
- Bastin J et al. Postnatal development of oxidative enzymes in various rat nephron segments. *Am J Physiol* 259:F895–F901 (1990). https://pubmed.ncbi.nlm.nih.gov/2260682/ †
- Benjamini Y, Hochberg Y. Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. *J R Stat Soc B* 57:289–300 (1995). https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
- Bignon Y et al. Multiomics reveals multilevel control of renal and systemic metabolism by the renal tubular circadian clock. *J Clin Invest* 133:e167133 (2023). https://pubmed.ncbi.nlm.nih.gov/36862511/
- Brehe JE et al. Effect of methionine sulfoximine on glutathione and amino acid levels in the nephron. *Am J Physiol* 231:1536–1540 (1976). https://pubmed.ncbi.nlm.nih.gov/998800/ †
- Breljak D et al. Distribution of organic anion transporters NaDC3 and OAT1-3 along the human nephron. *Am J Physiol Renal Physiol* 311:F227–F238 (2016). https://pubmed.ncbi.nlm.nih.gov/27053689/
- Chen L et al. A Comprehensive Map of mRNAs and Their Isoforms across All 14 Renal Tubule Segments of Mouse. *J Am Soc Nephrol* 32:897–912 (2021). https://pubmed.ncbi.nlm.nih.gov/33769951/
- Chen L et al. Multiomics Analyses Reveal Sex Differences in Mouse Renal Proximal Subsegments. *J Am Soc Nephrol* 34:829–845 (2023). https://pubmed.ncbi.nlm.nih.gov/36758122/
- Chen S et al. Multi-omic and spatial analysis of mouse kidneys highlights sex-specific differences in gene regulation across the lifespan. *Nat Genet* 57:1213–1227 (2025). https://pubmed.ncbi.nlm.nih.gov/40259083/
- Chevalier RL. The proximal tubule is the primary target of injury and progression of kidney disease: role of the glomerulotubular junction. *Am J Physiol Renal Physiol* 311:F145–F161 (2016). https://pubmed.ncbi.nlm.nih.gov/27194714/
- Conjard A et al. Gluconeogenesis from glutamine and lactate in the isolated human renal proximal tubule: longitudinal heterogeneity and lack of response to adrenaline. *Biochem J* 360:371–377 (2001). https://pubmed.ncbi.nlm.nih.gov/11716765/
- Crowell HL et al. muscat detects subpopulation-specific state transitions from multi-sample multi-condition single-cell transcriptomics data. *Nat Commun* 11:6077 (2020). https://pubmed.ncbi.nlm.nih.gov/33257685/
- Curthoys NP, Moe OW. Proximal tubule function and response to acidosis. *Clin J Am Soc Nephrol* 9:1627–1638 (2014). https://pubmed.ncbi.nlm.nih.gov/23908456/
- CZI Cell Science Program et al. CZ CELLxGENE Discover. *Nucleic Acids Res* 53:D886–D900 (2025). https://pubmed.ncbi.nlm.nih.gov/39607691/
- de Caestecker M et al. Bridging Translation by Improving Preclinical Study Design in AKI. *J Am Soc Nephrol* 26:2905–2916 (2015). https://pubmed.ncbi.nlm.nih.gov/26538634/
- Edison EE et al. Creatine synthesis: production of guanidinoacetate by the rat and human kidney in vivo. *Am J Physiol Renal Physiol* 293:F1799–F1804 (2007). https://pubmed.ncbi.nlm.nih.gov/17928413/ †
- Elias JE, Gygi SP. Target-decoy search strategy for increased confidence in large-scale protein identifications by mass spectrometry. *Nat Methods* 4:207–214 (2007). https://pubmed.ncbi.nlm.nih.gov/17327847/
- Fang Z et al. GSEApy: a comprehensive package for performing gene set enrichment analysis in Python. *Bioinformatics* 39:btac757 (2023). https://pubmed.ncbi.nlm.nih.gov/36426870/
- Faure L et al. scFates: a scalable python package for advanced pseudotime and bifurcation analysis from single-cell data. *Bioinformatics* 39:btac746 (2023). https://pubmed.ncbi.nlm.nih.gov/36394263/
- Feola K et al. Dynamic and differential renal cortical cell-specific mitochondrial metabolism in response to fasting. *Am J Physiol Renal Physiol* 330:F210–F222 (2026). https://pubmed.ncbi.nlm.nih.gov/41428383/
- Fernandes AL et al. Renal Phosphate Reabsorption in Humans Depends on at Least Three Distinct Transporters Unlike in Mice. *Acta Physiol (Oxf)* 242:e70271 (2026). https://pubmed.ncbi.nlm.nih.gov/42494142/
- Fu Y et al. Rodent models of AKI and AKI-CKD transition: an update in 2024. *Am J Physiol Renal Physiol* 326:F563–F583 (2024). https://pubmed.ncbi.nlm.nih.gov/38299215/
- Gillespie M et al. The reactome pathway knowledgebase 2022. *Nucleic Acids Res* 50:D687–D692 (2022). https://pubmed.ncbi.nlm.nih.gov/34788843/
- Goeman JJ, Bühlmann P. Analyzing gene expression data in terms of gene sets: methodological issues. *Bioinformatics* 23:980–987 (2007). https://pubmed.ncbi.nlm.nih.gov/17303618/
- Guo J et al. SMURF: soft-segmentation for single-cell reconstruction and topological analysis of spatial transcriptomic data. *Nat Commun* 17:7990 (2026). https://pubmed.ncbi.nlm.nih.gov/42350392/
- Haghverdi L et al. Diffusion pseudotime robustly reconstructs lineage branching. *Nat Methods* 13:845–848 (2016). https://pubmed.ncbi.nlm.nih.gov/27571553/
- Halpern KB et al. Single-cell spatial reconstruction reveals global division of labour in the mammalian liver. *Nature* 542:352–356 (2017). https://pubmed.ncbi.nlm.nih.gov/28166538/
- Heyman SN et al. Experimental ischemia-reperfusion: biases and myths—the proximal vs. distal hypoxic tubular injury debate revisited. *Kidney Int* 77:9–16 (2010). https://pubmed.ncbi.nlm.nih.gov/19759527/
- Hou W et al. A statistical framework for differential pseudotime analysis with multiple single-cell RNA-seq samples. *Nat Commun* 14:7286 (2023). https://pubmed.ncbi.nlm.nih.gov/37949861/
- Hu JJ et al. Mouse renal cytochrome P450IIE1: immunocytochemical localization, sex-related difference and regulation by testosterone. *Biochem Pharmacol* 40:2597–2602 (1990). https://pubmed.ncbi.nlm.nih.gov/2260985/ †
- Kaminska M et al. Heterogeneity in lysosomal dynamics and metabolic functions along the kidney proximal tubule. *Nat Commun* 17:3677 (2026). https://pubmed.ncbi.nlm.nih.gov/41803103/ †
- Kanehisa M et al. New approach for understanding genome variations in KEGG. *Nucleic Acids Res* 47:D590–D595 (2019). https://pubmed.ncbi.nlm.nih.gov/30321428/
- Kang HM et al. Defective fatty acid oxidation in renal tubular epithelial cells has a key role in kidney fibrosis development. *Nat Med* 21:37–46 (2015). https://pubmed.ncbi.nlm.nih.gov/25419705/
- Kaplan et al. (Sophont). OpenMidnight model card and blog post, "How to Train a State-of-the-Art Pathology Foundation Model with $1.6k" (2025). https://huggingface.co/SophontAI/OpenMidnight [CITE: no peer-reviewed or arXiv record found]
- Kirillov A et al. Panoptic Segmentation. CVPR 2019; arXiv:1801.00868. https://arxiv.org/abs/1801.00868
- Klötzer KA et al. Analysis of individual patient pathway coordination in a cross-species single-cell kidney atlas. *Nat Genet* 57:1922–1934 (2025). https://pubmed.ncbi.nlm.nih.gov/40775269/
- Korsunsky I et al. Fast, sensitive and accurate integration of single-cell data with Harmony. *Nat Methods* 16:1289–1296 (2019). https://pubmed.ncbi.nlm.nih.gov/31740819/
- Kuleshov MV et al. Enrichr: a comprehensive gene set enrichment analysis web server 2016 update. *Nucleic Acids Res* 44:W90–W97 (2016). https://pubmed.ncbi.nlm.nih.gov/27141961/
- Lake BB et al. An atlas of healthy and injured cell states and niches in the human kidney. *Nature* 619:585–594 (2023). https://pubmed.ncbi.nlm.nih.gov/37468583/
- Law CW et al. voom: Precision weights unlock linear model analysis tools for RNA-seq read counts. *Genome Biol* 15:R29 (2014). https://pubmed.ncbi.nlm.nih.gov/24485249/
- Layton AT et al. A computational model of epithelial solute and water transport along a human nephron. *PLoS Comput Biol* 15:e1006108 (2019). https://pubmed.ncbi.nlm.nih.gov/30802242/
- Le Hir M et al. Peroxisomal and mitochondrial beta-oxidation in the rat kidney: distribution of fatty acyl-coenzyme A oxidase and 3-hydroxyacyl-coenzyme A dehydrogenase activities along the nephron. *J Histochem Cytochem* 30:441–444 (1982). https://pubmed.ncbi.nlm.nih.gov/7200500/ †
- Lee JW et al. Deep Sequencing in Microdissected Renal Tubules Identifies Nephron Segment-Specific Transcriptomes. *J Am Soc Nephrol* 26:2669–2677 (2015). https://pubmed.ncbi.nlm.nih.gov/25817355/
- Legouis D et al. Altered proximal tubular cell glucose metabolism during acute kidney injury is associated with mortality. *Nat Metab* 2:732–743 (2020). https://pubmed.ncbi.nlm.nih.gov/32694833/
- Liberzon A et al. The Molecular Signatures Database (MSigDB) hallmark gene set collection. *Cell Syst* 1:417–425 (2015). https://pubmed.ncbi.nlm.nih.gov/26771021/
- Limbutara K et al. Quantitative Proteomics of All 14 Renal Tubule Segments in Rat. *J Am Soc Nephrol* 31:1255–1266 (2020). https://pubmed.ncbi.nlm.nih.gov/32358040/
- Love MI et al. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol* 15:550 (2014). https://pubmed.ncbi.nlm.nih.gov/25516281/
- Maass C et al. Translational Assessment of Drug-Induced Proximal Tubule Injury Using a Kidney Microphysiological System. *CPT Pharmacometrics Syst Pharmacol* 8:316–325 (2019). https://pubmed.ncbi.nlm.nih.gov/30869201/
- MacKinnon JG, White H. Some heteroskedasticity-consistent covariance matrix estimators with improved finite sample properties. *J Econometrics* 29:305–325 (1985). https://doi.org/10.1016/0304-4076(85)90158-7
- Makover A et al. Localization of retinol-binding protein messenger RNA in the rat kidney and in perinephric fat tissue. *J Lipid Res* 30:171–180 (1989). https://pubmed.ncbi.nlm.nih.gov/2469758/ †
- Margaillan G et al. Quantitative profiling of human renal UDP-glucuronosyltransferases and glucuronidation activity. *Drug Metab Dispos* 43:611–619 (2015). https://pubmed.ncbi.nlm.nih.gov/25650382/ †
- Marshall JL et al. High-resolution Slide-seqV2 spatial transcriptomics enables discovery of disease-specific cell neighborhoods and pathways. *iScience* 25:104097 (2022). https://pubmed.ncbi.nlm.nih.gov/35372810/
- Maunsbach AB et al. Aquaporin-1 water channel expression in human kidney. *J Am Soc Nephrol* 8:1–14 (1997). https://pubmed.ncbi.nlm.nih.gov/9013443/
- Melo Ferreira R et al. Integration of spatial and single-cell transcriptomics localizes epithelial cell-immune cross-talk in kidney injury. *JCI Insight* 6:e147703 (2021). https://pubmed.ncbi.nlm.nih.gov/34003797/
- Moor AE et al. Spatial Reconstruction of Single Enterocytes Uncovers Broad Zonation along the Intestinal Villus Axis. *Cell* 175:1156–1167 (2018). https://pubmed.ncbi.nlm.nih.gov/30270040/
- Murphy AE, Skene NG. A balanced measure shows superior performance of pseudobulk methods in single-cell RNA-sequencing analysis. *Nat Commun* 13:7851 (2022). https://pubmed.ncbi.nlm.nih.gov/36550119/
- Muzellec B et al. PyDESeq2: a python package for bulk RNA-seq differential expression analysis. *Bioinformatics* 39:btad547 (2023). https://pubmed.ncbi.nlm.nih.gov/37669147/
- Nakagawa J et al. Molecular characterization of mammalian dicarbonyl/L-xylulose reductase and its localization in kidney. *J Biol Chem* 277:17883–17891 (2002). https://pubmed.ncbi.nlm.nih.gov/11882650/ †
- Navarro Garrido A et al. Aristolochic acid-induced nephropathy is attenuated in mice lacking the neutral amino acid transporter B0AT1 (Slc6a19). *Am J Physiol Renal Physiol* 323:F455–F467 (2022). https://pubmed.ncbi.nlm.nih.gov/35979966/ †
- Neufeld A et al. Inference after latent variable estimation for single-cell RNA sequencing data. *Biostatistics* 25:270–287 (2023). https://pubmed.ncbi.nlm.nih.gov/36511385/
- Nigam SK et al. The organic anion transporter (OAT) family: a systems biology perspective. *Physiol Rev* 95:83–123 (2015). https://pubmed.ncbi.nlm.nih.gov/25540139/
- Nitzan M et al. Gene expression cartography. *Nature* 576:132–137 (2019). https://pubmed.ncbi.nlm.nih.gov/31748748/
- Oliveira MF et al. High-definition spatial transcriptomic profiling of immune cell populations in colorectal cancer. *Nat Genet* 57:1512–1523 (2025). https://pubmed.ncbi.nlm.nih.gov/40473992/
- Oquab M et al. DINOv2: Learning Robust Visual Features without Supervision. arXiv:2304.07193 (2023). https://arxiv.org/abs/2304.07193
- Parks LD et al. Heterogeneity of glutathione synthesis and secretion in the proximal tubule of the rabbit. *Am J Physiol* 274:F924–F931 (1998). https://pubmed.ncbi.nlm.nih.gov/9612330/ †
- Perco P et al. Identification of dicarbonyl and L-xylulose reductase as a therapeutic target in human chronic kidney disease. *JCI Insight* 4:128120 (2019). https://pubmed.ncbi.nlm.nih.gov/31217356/ †
- Polański K et al. Bin2cell reconstructs cells from high resolution Visium HD data. *Bioinformatics* 40:btae546 (2024). https://pubmed.ncbi.nlm.nih.gov/39250728/
- Rakhshandehroo M et al. Comparative analysis of gene regulation by the transcription factor PPARalpha between mouse and human. *PLoS One* 4:e6796 (2009). https://pubmed.ncbi.nlm.nih.gov/19710929/
- Ransick A et al. Single-Cell Profiling Reveals Sex, Lineage, and Regional Diversity in the Mouse Kidney. *Dev Cell* 51:399–413 (2019). https://pubmed.ncbi.nlm.nih.gov/31689386/
- Reichold M et al. Glycine Amidinotransferase (GATM), Renal Fanconi Syndrome, and Kidney Failure. *J Am Soc Nephrol* 29:1849–1858 (2018). https://pubmed.ncbi.nlm.nih.gov/29654216/ †
- Roux de Bézieux H et al. Trajectory inference across multiple conditions with condiments. *Nat Commun* 15:833 (2024). https://pubmed.ncbi.nlm.nih.gov/38280860/
- Saelens W et al. A comparison of single-cell trajectory inference methods. *Nat Biotechnol* 37:547–554 (2019). https://pubmed.ncbi.nlm.nih.gov/30936559/
- Smyth GK. Linear models and empirical bayes methods for assessing differential expression in microarray experiments. *Stat Appl Genet Mol Biol* 3:Article3 (2004). https://pubmed.ncbi.nlm.nih.gov/16646809/
- Song Y et al. Benchmarking strategies for cross-species integration of single-cell RNA sequencing data. *Nat Commun* 14:6495 (2023). https://pubmed.ncbi.nlm.nih.gov/37838716/
- Speerschneider P et al. Renal tumorigenicity of 1,1-dichloroethene in mice: the role of male-specific expression of cytochrome P450 2E1. *Toxicol Appl Pharmacol* 130:48–56 (1995). https://pubmed.ncbi.nlm.nih.gov/7839370/ †
- Squair JW et al. Confronting false discoveries in single-cell differential expression. *Nat Commun* 12:5692 (2021). https://pubmed.ncbi.nlm.nih.gov/34584091/
- Subramanian A et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *PNAS* 102:15545–15550 (2005). https://pubmed.ncbi.nlm.nih.gov/16199517/
- Takeda M et al. Intranephron distribution of glycine-amidinotransferase activity in rats. *Ren Physiol Biochem* 15:113–118 (1992). https://pubmed.ncbi.nlm.nih.gov/1378964/ †
- Terao M et al. Structure and function of mammalian aldehyde oxidases. *Arch Toxicol* 90:753–780 (2016). https://pubmed.ncbi.nlm.nih.gov/26920149/ †
- Thakur A et al. Sex and the Kidney Drug-Metabolizing Enzymes and Transporters: Are Preclinical Drug Disposition Data Translatable to Humans? *Clin Pharmacol Ther* 116:235–246 (2024). https://pubmed.ncbi.nlm.nih.gov/38711199/
- Traag VA et al. From Louvain to Leiden: guaranteeing well-connected communities. *Sci Rep* 9:5233 (2019). https://pubmed.ncbi.nlm.nih.gov/30914743/
- Vallon V et al. SGLT2 mediates glucose reabsorption in the early proximal tubule. *J Am Soc Nephrol* 22:104–112 (2011). https://pubmed.ncbi.nlm.nih.gov/20616166/
- Van den Berge K et al. Trajectory-based differential expression analysis for single-cell sequencing data. *Nat Commun* 11:1201 (2020). https://pubmed.ncbi.nlm.nih.gov/32139671/
- van der Walt S et al. scikit-image: image processing in Python. *PeerJ* 2:e453 (2014). https://pubmed.ncbi.nlm.nih.gov/25024921/
- Vrhovac I et al. Localizations of Na+-D-glucose cotransporters SGLT1 and SGLT2 in human kidney and of SGLT1 in human small intestine, liver, lung, and heart. *Pflugers Arch* 467:1881–1898 (2015). https://pubmed.ncbi.nlm.nih.gov/25304002/
- Wigger L et al. Expression Landscape and Circadian Regulation of lncRNAs in the Kidney. *Acta Physiol (Oxf)* 242:e70273 (2026). https://pubmed.ncbi.nlm.nih.gov/42400070/
- Wolf FA et al. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol* 19:15 (2018). https://pubmed.ncbi.nlm.nih.gov/29409532/
- Wu D, Smyth GK. Camera: a competitive gene set test accounting for inter-gene correlation. *Nucleic Acids Res* 40:e133 (2012). https://pubmed.ncbi.nlm.nih.gov/22638577/
- Xiong L et al. Direct androgen receptor control of sexually dimorphic gene expression in the mammalian kidney. *Dev Cell* 58:2338–2358 (2023). https://pubmed.ncbi.nlm.nih.gov/37673062/
- Yates B et al. Updates to HCOP: the HGNC comparison of orthology predictions tool. *Brief Bioinform* 22:bbab155 (2021). https://pubmed.ncbi.nlm.nih.gov/33959747/
- Zhou J et al. Unified Mouse and Human Kidney Single-Cell Expression Atlas Reveal Commonalities and Differences in Disease States. *J Am Soc Nephrol* 34:1843–1862 (2023). https://pubmed.ncbi.nlm.nih.gov/37639336/
- Zimmerman KD et al. A practical solution to pseudoreplication bias in single-cell studies. *Nat Commun* 12:738 (2021). https://pubmed.ncbi.nlm.nih.gov/33531494/
