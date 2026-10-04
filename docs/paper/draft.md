# Manuscript draft: where human and mouse proximal tubule differ along the tubule axis

**Version 0.2 (partial).** Text fixes from the internal referee report (section C, M5, M4) and
queued items 1, 2 and 4 of `v02_changes.md`. R2, R3 and the abstract were not rewritten; they
change after the revision analyses. Sentences that wait on a revision analysis are marked
**[PENDING D4]** or **[PENDING D8]**. The change log is in `draft_v02_changelog.md`.

Assembled by workstream 4 on 2026-10-04 from the committed results documents (`docs/results/`,
`docs/paper/`) and the saved result tables they cite. Each number was checked against those tables.

How to read the markers:

- **[VERIFY]**: needs lab records, or a number that no committed document or table supports.
- **[CITE]**: needs a reference.
- **(disclosed)**: a decision taken after results were seen. It must stay in the paper.

Citations are author–year. Section 8 gives PubMed or DOI links; every entry was checked in Europe
PMC, Crossref or on the publisher page. Section 9 lists everything still open.

---

## 1. Title options

1. Beyond a conserved transport core, human and mouse proximal tubule zonate different genes
2. Human and mouse proximal tubule share their transport zonation but not their metabolic zonation
3. A positional comparison of human and mouse proximal tubule from segmented spatial
   transcriptomics
4. Mouse proximal tubule zonates metabolic programs that human proximal tubule keeps flat

## 2. Abstract (≈ 230 words)

The proximal tubule (PT) carries out most renal reabsorption, drug secretion and metabolism, and is
studied mainly in mice. Whether the two species arrange these functions alike along the S1→S3 axis
is unknown, because comparisons of cell types or segments collapse position.

We profiled two control mouse kidneys and two cortical sections from one human donor with Visium HD.
We segmented tubule cross-sections and ordered 12,866 PT structures along a continuous
cross-species coordinate. The coordinate recovered measured mouse zonation and did not absorb an
injected species difference, but it did not resolve order finer than segments.

With two specimens per species, average-level pathway screens called as many pathways for
mixed-species relabelings of the specimens as for species. A position-dependent test called 60
pathways for species and 2 and 0 for the relabelings, as the split-plot design predicts.

Human and mouse S1→S2 gene gradients correlated weakly (noise-corrected r = 0.14; 0.23 in
independent atlases), and human gradients were about half as large. Of 7,407 probe-balanced genes:

- 155 were zonated alike in both species (apical transporters, intrarenal renin–angiotensin genes);
- 247 were zonated only in mouse (sterol synthesis, peroxisomal β-oxidation, Gatm);
- 66 were zonated only in human;
- 52, including DCXR and UGT3A1, were zonated in opposite directions;
- 5,590 could not be classified.

Mouse-only zonation was not explained by male-biased genes. These differences agree with
independent segment-resolved data at the gene-measurement level, though not pathway by pathway.
Mouse PT predicts where human PT places transport, but not metabolism. All results are
descriptive.

---

## 3. Introduction

**Translation of PT biology depends on knowing where species differ.** The PT reabsorbs about
60–70% of filtered water and NaCl, more of the filtered bicarbonate and nearly all filtered
nutrients. It also secretes solutes, produces hormones and carries out most of the kidney's
metabolism (Curthoys & Moe 2014). During fasting and stress the kidney supplies up to 40% of
systemic glucose production, and PT gluconeogenesis is impaired in acute kidney injury (AKI)
(Legouis et al. 2020). Defective fatty-acid oxidation in tubular epithelium contributes to fibrosis
(Kang et al. 2015). The SLC22 organic anion transporters handle common drugs, toxins, nutrients
and putative uremic toxins (Nigam et al. 2015).

The PT is regarded as the primary target of injury in AKI and chronic kidney disease (Chevalier
2016), and as the primary site of drug-induced kidney injury (Maass et al. 2019). Much of this
biology is studied in rodents. No intervention has yet prevented AKI in patients, which has raised
the question of whether animal models predict human responses (de Caestecker et al. 2015). Rodent
AKI models do not fully recapitulate the human disease (Fu et al. 2024), and experimental models
differ in where along the nephron tubular injury falls (Heyman et al. 2010). Preclinical models
often fail to predict drug nephrotoxicity (Maass et al. 2019), and species differences in renal
drug transporters and metabolising enzymes contribute to poor translation (Basit et al. 2019; Thakur
et al. 2024).

**The PT is zonated, and some of its zonation is known to be conserved.** Along its length the PT
is divided into S1, S2 and S3, and functions are distributed along this axis:

- In mouse, SGLT2 sits in the early PT (Vallon et al. 2011).
- In isolated human PT segments, S2 and S3 make more glucose from lactate than S1 (Conjard et al.
  2001).
- After ischaemia in mice, damage is most severe in S3 of the outer stripe (Kim et al. 2011).

Segment-resolved transcriptomes exist for microdissected rat (Lee et al. 2015) and mouse tubules
(Chen et al. 2021). Single-cell atlases describe PT diversity in mouse (Ransick et al. 2019) and
human (Lake et al. 2023).

Some canonical transport markers sit at the same position in humans and rodents. SGLT2 lies in
S1/S2 and SGLT1 in S3, in locations that "largely resembled those in rats and mice" (Vrhovac et al.
2015). Human protein maps also place AQP1 highest in the proximal straight tubule (Maunsbach et al.
1997) and OAT1–3 higher in S1/S2 than in S3 (Breljak et al. 2016). These two maps have not been
compared position by position with rodents.

Other programs differ. In mouse, sexually dimorphic gene activity maps predominantly to PT segments,
but only a limited set of genes shows conserved sex-linked regulation in human (Xiong et al. 2023).
The largest mouse PT sex differences are in S2/S3 (Chen et al. 2023).

These observations come from different studies, assays and segment definitions. Computational
models of the human nephron still take most transporter densities from rat (Layton et al. 2019). To
our knowledge, no study has compared human and mouse PT gene programs position by position along
the S1→S3 axis (literature check, Methods).

**Comparing cell types or segments cannot show where along the PT the species differ.**
Cross-species single-cell atlases align conserved cell states between human and rodent kidney
(Klötzer et al. 2025). A unified mouse–human kidney atlas found little overlap in single-cell
expression changes but agreement at the pathway level (Zhou et al. 2023). These comparisons ask
which cells correspond, not where along a segment a program changes. Dissociation also discards the
spatial relationships between cells (Nitzan et al. 2019).

Assigning structures to S1, S2 and S3 turns the axis into three bins. From those labels alone one
cannot tell whether expression changes in steps or gradually, or whether a boundary sits at the
same place in both species. Continuous reconstructions in other organs show what binning can hide.
About half of mouse liver genes are zonated along the lobule, and many of them peak in the middle
layers (Halpern et al. 2017). Enterocyte functions are broadly zonated along the intestinal villus
(Moor et al. 2018).

**Spatial transcriptomics and pseudospace reconstruction.** Spatial transcriptomics keeps tissue
context. In kidney it has mapped cell types and injury niches together with single-nucleus data
(Melo Ferreira et al. 2021), with Slide-seqV2 (Marshall et al. 2022) and in the human kidney atlas
(Lake et al. 2023). Visium HD measures the whole transcriptome at single-cell-scale resolution
(Oliveira et al. 2025).

A section, however, cuts each tubule at an arbitrary point, so each structure's position along the
nephron has to be reconstructed. Several approaches exist:

- inferring lobule and villus coordinates from landmark genes (Halpern et al. 2017; Moor et al.
  2018);
- recovering spatial arrangements from expression alone (Nitzan et al. 2019);
- unrolling spatial data to expose zonated programs (Guo et al. 2026);
- trajectory methods, which order cells along continuous processes (Saelens et al. 2019). These
  include diffusion pseudotime (Haghverdi et al. 2016) and principal curves (Faure et al. 2023).

**Few specimens make species comparisons hard to separate from specimen differences.** Cells or
structures from one specimen are not independent replicates. Methods that ignore variation between
biological replicates are biased, and they can report hundreds of differentially expressed genes
where no biological difference exists (Squair et al. 2021; Zimmerman et al. 2021). Pseudobulk
aggregation is conservative but underpowered relative to mixed models (Zimmerman et al. 2021); both
approaches have been benchmarked for multi-sample designs (Crowell et al. 2020).

For pseudotime, Lamian showed what happens when samples are randomly partitioned into groups with
no expected difference: methods that ignore sample variation report thousands of false-positive
genes (Hou et al. 2023). Competitive gene-set tests that treat genes as independent are
anti-conservative when genes in a set are correlated (Goeman & Bühlmann 2007; Wu & Smyth 2012).

Human tissue is scarce, and spatial studies often have only a few sections per species. An analysis
therefore has to show that its calls are specific to the species contrast, not merely that they
are numerous.

**This study.** We profiled two control mouse kidney sections and two cortex sections from one human
donor with Visium HD. We segmented individual tubule cross-sections from the paired H&E images,
aggregated 2-µm bins into tubule-level profiles in a shared one-to-one ortholog space, and
reconstructed a continuous PT coordinate across species. We then tested what that coordinate does
and does not resolve.

To decide which pathway screens are specific with two specimens per species, we used a split-plot
argument and a specimen-relabeling control. We called zonation separately in each species and
confirmed every call in an independent atlas of the same species. Finally, we tested lead genes
against microdissected mouse tubules of both sexes, mouse single-nucleus data and a cortex-only human
atlas.

Human PT is not simply a flattened mouse PT. The two species share the zonation of apical solute
transport and of the intrarenal renin–angiotensin genes. Mouse PT additionally zonates a set of
metabolic programs that human PT keeps flat. Human PT has fewer zonated genes of its own, and a
handful of genes change direction.

---

## 4. Results

### R1 · A continuous cross-species PT coordinate

**Cohort.**

- *Segmentation and quality control.* We segmented tubule cross-sections in two control mouse kidney
  sections and two cortex sections from one human donor. After structure quality control, 26,839
  structures remained (Methods).
- *PT structures.* Cross-species integration and reviewed clustering identified 12,866 PT structures
  (mouse 2,852 and 2,808; human 3,438 and 3,768), each labelled S1, S2 or S3 by marker review.

**Coordinate.** We ordered the PT structures along a nonbranching principal curve fitted to the
integrated embedding (scFates; Methods; Figure 1a, b).

- In every specimen the median coordinate of S3 exceeded those of S1 and S2.
- S2→S3 transitions fell at 0.63–0.64 in all four specimens.
- S1→S2 transitions fell at 0.26–0.28 in mouse and 0.35–0.38 in human, so human S1 occupies more of
  the coordinate.

**Recovered zonation.** Each gene was scored only on a coordinate rebuilt without it (three gene
folds). Across genes, its within-specimen rank correlation with the coordinate was then compared
with external S3-versus-early contrasts (Figure 1c):

| Species | External reference | Coordinate | Segment labels |
|---|---|---|---|
| Mouse | Male snRNA | 0.74 | 0.73 |
| Mouse | Microdissection | 0.53 | 0.51 |
| Mouse | Rat microdissection (male) | 0.11 | 0.10 |
| Human | Cortex snRNA (resolves only convoluted PT versus S3) | 0.25 | 0.27 |

Without fold protection, agreement was 0.75 for mouse snRNA and 0.31 for human, so part of the human
agreement came from the genes that built the coordinate. Rat zonation agreed poorly with our mouse
data: in the segment-label comparison of notebook 34, Spearman was 0.13, with 55% direction
agreement for strongly zonated genes. Rat is therefore not used as a proxy for mouse.

**No absorption of species differences.**

- *Injected difference.* A human-only gradient injected into a tenth of the construction genes
  survived refitting (absorption ratio 1.01; Figure 1d).
- *External sign agreement.* Among 244 genes whose zonation differs between species in external
  data, the coordinate kept the external direction for 91%.
- *Independent reads.* We built the coordinate on one half of the reads and tested on the other.
  59 of 66 robust and 39 of 42 core pathways from an earlier list were still called.
- *Circular use.* Reusing the same reads inflated gene statistics by about 30%. The circular arm
  called more pathways than the independent arm (145 against 130), but kept a similar share of the
  earlier robust list (57 against 59 of 66).

**What the coordinate does not resolve** (Figure S1).

- Within reviewed segments, slopes along the coordinate agreed only weakly with adjacent external
  contrasts: 0.17 in mouse and 0.10 in human, against about 0 for ordering by read depth.
- Coordinates built from disjoint thirds of the genes ranked structures within segments only
  moderately alike (0.34–0.70).
- The coordinate predicted held-out genes in the other specimen no better than segment means.
- After thinning to equal depth, human S1 showed a depth dependence. Within-segment agreement was
  0.36–0.46, and the correlation of position with library size fell from −0.31 and −0.38 to −0.04
  and −0.13.
- Removing a third of the genes or half the reads moved human segments by up to 0.29 relative to
  mouse.
- The curve is numerically sensitive: input changes of order 10⁻¹⁶ moved structures by up to 0.44
  (Methods).

**Alternatives** (Table S1).

- *The pre-registered coordinate-selection rule selected DPT.* The rule required three things:
  - equal-depth agreement ≥ 0.80;
  - external agreement within 0.05 of scFates in each species;
  - a species-registration gap ≤ 0.15.

  Diffusion pseudotime (DPT) met all three (0.80; 0.74 mouse and 0.29 human; gap 0.14). It was
  better than scFates on every pre-specified coordinate-quality metric:

  | Metric | DPT | scFates |
  |---|---:|---:|
  | Human external agreement | 0.29 | 0.25 |
  | Gene-fold within-segment agreement | 0.59 | 0.49 |
  | Equal-depth agreement | 0.80 | 0.72 |
  | Registration gap | 0.14 | 0.29 |

  scFates itself would not have met the rule. The equal-depth refit of scFates also failed it
  (0.73; gap 0.36).
- *DPT failed downstream specificity.* On DPT the positional pathway screen was not specific
  (NCDR 0.40, defined in R2; matched test). Its relabeled calls were metabolic and transport
  programs, and it found 36 positional pathway differences between the two mice, against 16 on
  scFates.
- *A conserved-anchor coordinate,* built from 160 externally conserved zonation genes, failed to
  register the two species.
- *Choice (disclosed).* The selection rule did not include the pre-registered downstream
  specificity gate: NCDR < 0.25 and fewer than 5% of pathways called under relabeling (criterion 5
  of the conserved-anchor rule C6). After DPT was seen to fail that gate, the coordinator applied it
  to every candidate. No candidate then qualified, and the rule's fallback, the scFates coordinate,
  was kept as primary.
- *Consequence.* Specificity of the positional screen is a property of the coordinate as well as of
  the design. Notebook 37's joint-test pipeline is being run end to end on DPT [PENDING D4].

**What the coordinate adds.** Not resolution finer than segments.

- It replicates measurements across positions within each specimen. That lets species-by-position
  differences be judged against specimen-by-position deviations (R2).
- It localises differences continuously, but peak positions should be read at about segment
  resolution.
- The early/late split of the 36 robust pathways is not an artefact of unequal positional noise.
  Ten human-high pathways peak at a median of 0.16 and 26 mouse-high pathways at 0.63. Equalising
  positional noise between species changes this only slightly. We added the estimated excess human
  noise to mouse positions (σ = 0.089, from gene-fold refits). The split moved to 11.9 / 24.1,
  which accounts for about 15% of the gap between the median peaks and 7% of the mouse-high count.
  With σ = 0.10 the figures were 14% and 8%.

Zonation classes (R3) and the R5 stories use reviewed segment labels and do not depend on the
coordinate.

### R2 · Average-level pathway screens are not specific with two specimens per species

**The design is a split plot.** With two specimens per species, the four specimens form exactly
three balanced 2 + 2 contrasts (Figure 2a): the species split, and two relabelings that each pair one
mouse with one human section. If specimen deviations are independent, each relabeling is one draw of
the species contrast's section-level null.

In this split-plot design (Altman & Krzywinski 2015), each kind of species effect has its own error
term:

- a constant species offset is judged against specimens within species, which leaves 2 degrees of
  freedom;
- a species-by-position difference is judged against specimen-by-position deviations, which leave
  2(k − 1) degrees of freedom over k positions.

A term shared by both human sections, such as donor, probe panel or sampled region, enters the
species contrast but cancels in the relabelings. The relabeling control cannot detect it; R4
addresses these terms separately.

**Average-level screens called as many pathways for relabeled groups as for species.** We reran each
pathway strategy unchanged on both relabelings (Figure 2b; Table S2). The negative-control discovery
ratio (NCDR) is the mean relabeled count divided by the species count.

| Strategy | Species calls | Relabeled calls | NCDR |
|---|---:|---:|---:|
| Whole-PT pseudobulk with signed GSEA | 158 | 193, 255 | 1.42 |
| S1/S2/S3 pseudobulk with signed GSEA | 204 | 343, 255 | 1.47 |
| Gene-level constant-offset statistic (T_level), matched test | 27 | 127, 81 | 3.85 |

Target–decoy estimates (Elias & Gygi 2007) agreed (Figure 2c). No competitive average-level screen
reached a decoy-estimated false-discovery proportion (FDP) of 10% at any threshold. The
self-contained pathway-score offset model, also average-level, reached it with 1,087 calls; see
the disclosed changes below. The lowest attainable FDPs
were 1.00, 0.75 and 0.95; with the joint test, T_level reached 0.96.

**The nonspecificity is in the pathway layer, not in the gene statistics.**

- Whole-PT DESeq2 called 9,010 genes for species and only 2 and 8 for the relabelings.
- Specimen-level variation moves whole pathways together, which gene-permutation tests assume away
  (Goeman & Bühlmann 2007; Wu & Smyth 2012). Decoy z² rose with pathway size, more steeply for the
  offset statistic (inflation a = 1.40, implied correlation ρ = 0.019) than for the position
  statistic (a = 1.03, ρ = 0.009) (Figure 2e).
- Pathway scores tested at the structure level called 1,444 and 1,509 of 1,513 pathways, and 67–1,317
  under relabeling, consistent with treating 12,272 structures as replicates (Squair et al. 2021).

Within-species controls agreed (Figure 2d). Comparing the two mice gave 48 offset and 16 positional
pathway calls; comparing the two human sections gave 94 and 0. Specimen differences are offsets.

**The smooth position-dependent test was specific.** It combines two parts:

- *Statistic.* T_spatial measures how much the species difference changes along the coordinate,
  beyond a constant offset.
- *Joint test.* A covariate-matched rank test with a correction for inter-gene correlation
  (Methods).

The joint test called 60 pathways for species and 2 and 0 for the relabelings (NCDR 0.02), and 57
at a decoy FDP of 10% (interval 23–78). Under relabeling, the matched test alone was
anti-conservative: on average 59.5 pathways reached p ≤ 0.01 where about 15 are expected. The
correlation estimated from model residuals (median ρ ≈ 0.010) matched the inflation estimated from
decoys (ρ ≈ 0.009).

A specimen-level split-plot test on specimen × 8-bin pseudobulks was calibrated under relabeling. It
ranked genes only weakly like T_spatial (Spearman 0.21), and we use it as a complement (R4).

**Disclosed changes to the plan.**

- Decoy FDP was the pre-specified primary metric. It failed for pathway-score models: their decoy FDP
  is 0, because almost every self-contained pathway score differs more between species than between
  sections. NCDR, together with the condition that relabelings call fewer than 5% of tested
  pathways, stayed primary.
- That 5% condition was itself added after the first run (Methods).

Conventional analyses are not wrong. They ask whether a pathway differs on average, and this design
cannot answer that question specifically.

### R3 · Beyond a conserved transport core, PT zonation is largely species-specific

All analyses in this section use reviewed segment labels and pseudobulk counts, not the
coordinate.

**Human PT is not simply a flatter mouse PT** (notebook 39; Figure S4). We asked first whether human
gradients are scaled-down mouse gradients. Both of our gradients (S2 − S1, cortical in both species)
were regressed on independent reference gradients, which keeps reference noise out of the ratio
(Methods). The answer depended on which species supplied the gene axis:

| Reference axis | Human ÷ mouse amplitude |
|---|---|
| Male mouse microdissection | 0.15 (95% CI 0.11–0.19) |
| Male mouse snRNA | 0.15 |
| Female mouse snRNA | 0.16 |
| Female mouse microdissection | 0.23 |
| Healthy human cortex atlas (Lake et al. 2023) | 1.96 (1.41–3.05) |

Without our data, the same split appeared: the human atlas divided by mouse snRNA gave 0.14 on a
mouse axis and 1.66 on a human axis (S3 − early).

Such a split arises when the two species zonate different genes. Post-hoc diagnostics, which were
not pre-registered, confirmed it:

- The noise-corrected cross-species correlation of S2 − S1 gradients was 0.14 (0.10–0.19) in our
  data and 0.23 in independent snRNA, against 0.83 between mouse AKI and control mice. For
  S3c − early in our data it was 0.29 (0.23–0.34).
- The noise-corrected spread of S2 − S1 gradients was about half as large in human: 0.58 in our data
  and 0.47 independently. For S3 − early with cortical mouse S3 it was 0.84.

**Log compression and tissue state** (pre-specified checks).

- *Log compression is not the cause.* The top expression tertile gave 0.18 and probe-balanced genes
  0.16.
- *Injury does not reproduce the S2 − S1 pattern, but tissue state is not excluded* [PENDING D8].
  - *Mouse AKI.* AKI kept 0.74 of control S2 − S1 amplitude (gene-gradient correlation 0.83). It
    cut S3 − early amplitude to 0.41–0.44 of control on mouse axes (correlation 0.59). That is
    close to the human/mouse S3c − early ratio of 0.28–0.33.
  - *Lake donors.* Injury scores correlated weakly with mouse-axis amplitude (ρ = −0.24, the
    pre-specified test). They correlated more strongly with human-axis amplitude (ρ = −0.41, post
    hoc). Median human-axis amplitude was 0.69 in healthy, 0.50 in CKD and 0.35 in AKI donors.
  - *Workstream 3 check.* A separate check on another amplitude measure met its pre-specified rule
    for tissue-state flattening (ρ = −0.42, 95% CI −0.66 to −0.10), although the effect it implies
    is small.
  - *Our structures.* Low- and high-injury halves gave the same S2 − S1 ratio (0.18 and 0.19).
  - *Unmeasured differences.* Our donor's tissue state is unknown. Procurement, age, fasting and
    medication differ between human tissue and mice, and the external human atlases share these
    differences.

Injury lowers zonation amplitude in both species, and healthy human donors remain flatter than
healthy mice. In these data injury does not reproduce the human S2 − S1 pattern. Whether tissue
state contributes remains open [PENDING D8].

**Zonation classes** (notebook 40; Figure 3; Table 1). We then called zonation separately in each
species:

- *Contrasts.* S2 − S1, cortical in both species; and S3 − early, with mouse S3 restricted to
  cortical-like structures so that the region matches human cortex.
- *Noise scale.* Each call was tested against that species' own specimen noise.
- *Confirmation.* Each call was confirmed in an independent atlas of the same species: 12 male mouse
  snRNA donors, and 7 healthy human cortex donors from KPMP.
- *Species-only genes.* A gene was called zonated in only one species when the other species was
  confidently flat, not merely non-significant. This prevents the noisier human data from inflating
  the mouse-only class.

Of 7,407 probe-balanced genes measurable in both species:

| Class | Genes |
|---|---:|
| Zonated in the same direction in both species (conserved) | 155 |
| Zonated in opposite directions (reversal) | 52 |
| Zonated only in mouse | 247 |
| Zonated only in human | 66 |
| Flat in at least one species and zonated in neither | 1,297 |
| Could not be classified with two specimens per species | 5,590 |

*Mouse-only genes outnumbered human-only genes under every alternative rule (A1–A7) and in both
region-matched contrasts* (Figure 3b). They did not when all mouse S3 was included: in S3 − early,
which is mostly outer stripe in mouse, the counts were 77 mouse-only against 111 human-only.

- Scaling the human effect threshold to the lower human amplitude narrowed the gap from 247 against
  66 to 243 against 141.
- Most of the excess came from S1 → S2 (219 against 44 genes).
- For the region-matched S3 contrast, the two classes were similar (48 against 35).
- The excess held when mouse calls were confirmed in female rather than male mouse donors:
  - female snRNA confirmation kept 214 of 247 mouse-only, 151 of 155 conserved and 44 of 52
    reversal calls;
  - female microdissection kept 149, 127 and 30.

  Male confirmation remains primary because our mice are male.

*Not a male-mouse artefact.* Mouse-only genes were not enriched for male-biased mouse PT genes (Xiong
et al. 2023 tables): 5.3% were male-biased, against 15.5% of conserved genes. 199 of the 247 were not
sex-biased in either direction.

**What each class contains.**

- *Conserved core: apical solute transport* (Figure 3d), with enrichments that survived finer
  expression matching and removal of low-expression genes:
  - neutral and cationic amino-acid transport (Slc7a8, Slc6a18, Slc7a7; Slc3a1 is
    probe-imbalanced);
  - organic anion and urate transport (Slc22a7, Slc22a12);
  - phosphate (Slc34a3), Slc5a2 and Slc5a10;
  - the intrarenal renin–angiotensin genes (Agt, Ace, Enpep).

  This agrees with the conserved localisation of SGLT2 and SGLT1 (Vrhovac et al. 2015). It does not
  extend to every classical marker:
  - Aqp1 is indeterminate in our classes.
  - Slc22a6 (OAT1) is mouse-only.
  - Slc22a7 (OAT2) is conserved in our transcript data, rising toward S3 in both species, whereas
    human OAT2 protein stains more strongly in S1/S2 than in S3 (Breljak et al. 2016).
- *Mouse-only zonation: metabolism.*
  - sterol synthesis (Lss, Ebp; Cyp51, probe-imbalanced);
  - peroxisomal very-long-chain β-oxidation (Acox1, Acox2, Hsd17b4, Slc27a2);
  - tyrosine and pyrimidine catabolism.

  It also includes creatine synthesis (Gatm), glutathione synthesis (Gclm; Gclc, probe-imbalanced),
  and the organic-anion uptake genes Slc22a6 (OAT1) and Slc13a3 (NaDC3).
- *Human-only and reversal classes are small.* Their pathway enrichments were mostly not robust to
  expression, so the human side is best described gene by gene:
  - reversals: Dcxr, Ugt3a1, Vnn1;
  - zonated only in human: Sel1l3, Pkhd1;
  - expressed and zonated only in human: Rbp4 and Aox1, both rising toward late PT. Both are absent
    from mouse PT.

Most genes remained indeterminate (5,590 of 7,407). Indeterminate genes include several earlier
leads, such as Acadm, Acaa2, Acsm3, Crot and Nudt19. Probe-imbalanced genes (1,622, 86 of them in a
zonated class) and genes expressed in only one species are reported separately (Table S6).

**Position-dependent pathways recover both the conserved core and the mouse-specific programs**
(notebook 37; Figure S3; Tables S3, S4).

*Funnel.* Of 1,513 tested pathways:

| Stage | Pathways |
|---|---:|
| Joint-test candidates | 60 |
| Not specimen-sensitive | 58 |
| Stable in ≥ 5 of 7 planned refits | 51 |
| Robust to coordinate registration and sampled region | 36 |
| Core: robust, and also robust to genes detected in both species and to removing mouse sex-biased genes | 26 |

*Programs.* The 36 robust pathways form 17 gene-overlap programs. The largest are:

- steroid-hormone and xenobiotic metabolism (Cyp2e1, Hsd11b1, Adh1, Ugt genes);
- mitochondrial and peroxisomal fatty-acid oxidation (Acsm3, Ehhadh, Acadm, Acaa2);
- vitamin and cofactor metabolism (Slc22a13, Slc5a8, Pank1);
- bile acid and peroxisome (Crot, Nudt19, Slc27a2);
- fat-soluble vitamins (Rbp4, Lpl, Apob);
- glutathione (Gclc, Gclm, Gss);
- branched-chain amino-acid degradation.

*By zonation class.* Of the 36 robust pathways:
- 19 have mostly mouse-only zonated members;
- 15 have mostly conserved members;
- 1 has mostly reversal members (Vitamin B5 metabolism);
- 1 has no zonated member (Chemical carcinogenesis).

Two of them, Steroid hormone biosynthesis (9 classified members) and Nitrogen metabolism (8), fall
below the 10-member floor for class enrichment. They are still assigned a dominant class: mouse-only
and conserved respectively.

*Peaks and shape.* The 10 human-high pathways peak early (median 0.16) and the 26 mouse-high
pathways late (median 0.63). Most differences are graded and one-signed: 26 graded, 4 localised and
6 reversing, with a median averaging loss of 0. A whole-PT average therefore dilutes these
differences more than it cancels them.

*Relative to global flattening.* Taking human amplitude as 0.15 of mouse, 32 of the 36 are consistent
with "mouse zonation, flat human". Four go beyond it (BH ≤ 0.10): Estrogen Response Early and Late,
Glutathione metabolism, and Formation of Cornified Envelope.

The early/late split is therefore mainly the signature of mouse late-PT programs that human PT does
not zonate. It is not evidence that the species place one shared program at different positions.

*Earlier list.* An earlier list built with a separate correlation filter (82 → 66 → 42) is kept as a
sensitivity analysis. It contains 27 of the 36 robust pathways and 17 of the 26 core pathways.

*Flagged drivers.* Some robust programs carry drivers flagged as neighbouring-segment, ambient or
hepatocyte-type transcripts (for example Akr1c18 and Ttr in the fat-soluble vitamin program). These
pathways keep ≥ 70% of their effect when flagged genes are removed (R4).

### R4 · Robustness and agreement with independent data

**Robustness.** Every robust pathway keeps at least 70% of its effect under each of these checks
(Figure 4a):

- an anatomy-anchored coordinate, in which human sections are warped so that their segment
  transitions fall on the mouse transitions;
- cortical-like mouse S3 only;
- removal of flagged ambient, neighbouring-segment and shared-exon genes;
- 5-bin covariate matching;
- equal probe counts in the two Visium HD panels.

The first two checks define "robust", so they also hold significance (joint q ≤ 0.10). So do the
flagged-gene and 5-bin checks. With equal-probe genes only, 8 of the 36 lose significance
(q > 0.10) although their effect is retained.

The two gene-subset checks that define the core set were passed by fewer pathways (effect ratio
≥ 0.7 and q ≤ 0.10):
- genes detected in both species: 28 of 36;
- mouse sex-biased genes removed: 32 of 36.

Together they leave the 26 core pathways.

DPT is not cited as robustness support until notebook 37 has been run on it end to end [PENDING
D4].

**Specimen-level split-plot test** (Methods).

- *Gene level.* Calibrated under relabeling: 0.6% and 0.4% of genes at p ≤ 0.01, λ ≤ 0.87
  (Figure 4e).
- *Pathway level, pre-specified joint test.* It called nothing. Decoys showed that its pseudobulk
  correlation correction overcorrects (decoy a = 1.10 against a median correction of 2.07).
- *Pathway level, matched test.* It called 41 pathways, against 0 and 0 under relabeling. 15 of the
  36 robust pathways were among them.
- *Robust versus uncalled pathways.* The robust pathways' split-plot z far exceeded that of uncalled
  pathways (median 2.6 against −0.5; Mann–Whitney p = 7 × 10⁻²⁰).

**External agreement is at the measurement level** (Figure 4c, d; Table S5). Donors are the
replicates. Human cortex snRNA (Census, 6 donors, primary; Lake/KPMP, 7 donors) was compared with 12
male or 12 female mice.

*Genes.* Our top 10% of genes by T_spatial correlated with the external species × position contrast
above expression-matched chance, but about as well as the most highly expressed genes:

| Reference | Our top-10% genes | Matched random genes, median (97.5th pct) | Highest-expressed genes |
|---|---:|---|---:|
| Census | 0.67 | 0.57 (0.62) | 0.64 |
| Lake/KPMP | 0.73 | 0.60 (0.66) | 0.72 |

*Pathways.* Robust pathways replicated far better than uncalled ones (Mann–Whitney p = 2 × 10⁻¹³
for Census and 1 × 10⁻⁹ for Lake). They did not clearly outperform pathways called only under
relabeling (p = 0.016 and 0.19 against Census with male and female mice; 0.033 and 0.21 against
Lake). Further:

- Only 17 of 36 replicated individually against Census, and 6 against Lake (BH ≤ 0.10 within the
  list). The pre-specified criterion (≥ 50%) was not met.
- All three human datasets share a global component: the interaction is −0.81 to −0.90 times the
  mouse gradient (R² 0.60–0.82; Table S5). Once it is removed, 1–3 robust pathways replicate in any
  reference.
- Pathways called only by conventional screens (184) replicated their average difference 6 times
  at nominal p ≤ 0.05 and never at BH ≤ 0.10.

We therefore claim agreement at the level of gene measurements, and a contrast with uncalled
pathways, but not pathway-by-pathway replication.

**Genome-wide comparison with the cortex-only human atlas** (segment labels).

- Our human S3 − S1 contrasts correlated with Lake at Spearman 0.43 (6,787 genes; 93.1% direction
  agreement for 335 strongly zonated genes).
- For the 10% most position-dependent genes, our interaction correlated with Lake minus male mouse
  snRNA at 0.70, with 96.7% direction agreement for 276 genes strong in both. With female mice as
  the reference, the values were 0.67 and 93.7%.

### R5 · Gene programs placed differently along the PT

**How the claims were tested.** Every claim below uses reviewed segment labels. Each was tested
against:

- the cortex-only human atlas (Lake et al. 2023; 7 donors, 4 female and 3 male);
- mouse snRNA by sex (12 male and 12 female donors);
- microdissected mouse segments from both sexes (Chen et al. 2021, 2023).

Rules were fixed before these data were opened. Of 39 claimed genes, 30 met every pre-specified test,
5 were supporting and 4 did not replicate (Figure S5). Each gene also carries its R3 zonation class.

**Caveats that apply to every story.**

- *Probe panels.* Genes whose panels carry unequal probe numbers cannot be headline results.
- *Attenuation.* A segment-label attenuation check (a = 0.86, 95% interval 0.64–1.35) does not
  overturn any "human flat" claim.

**Creatine synthesis is zoned in mouse and PT-wide in human** (mouse-only; Figure 5).

- *Mouse.* Gatm falls steeply from S1 to S3: by 6.8 (male snRNA), 6.3 (female snRNA), 12.1 (male
  microdissection) and 11.0 (female microdissection) log2 units. Even the cortical medullary-ray
  segment of mouse has 23.5 TPM against 4,270 TPM in S1.
- *Human.* Gatm stays high along the PT; the atlas shows no fall in cortex (S3 − S1 = +0.51 log2).
- *Confounders.* Probes are balanced. Mouse Gatm is, if anything, female-biased, so male mice may
  understate the effect.
- *Literature.* The rodent half is known: rat transamidinase activity is confined to S1 and S2
  (Takeda et al. 1992). GATM is a human PT enzyme (Reichold et al. 2018), and the human kidney
  releases guanidinoacetate in vivo (Edison et al. 2007). The human axial profile and the species
  contrast were not found in the literature.

**Peroxisomal and mitochondrial fatty-acid oxidation.**

- *Peroxisomal very-long-chain β-oxidation is mouse-only zonation* (Acox1, Acox2, Hsd17b4,
  Slc27a2), and it is one of the enriched mouse-only programs.
- *Mitochondrial steps are mouse-zoned in the external references but indeterminate in our
  classification.*
  - In mouse references, Acadm (MCAD) rises from S1 to S3 (+2.1 and +2.3 in male and female snRNA)
    and Acaa2 is confined to S1 (−4.0 and −4.2); microdissection agrees in both sexes.
  - In the human atlas ACADM is nearly flat (S3 − S1 +0.46).
  - ACAA2 rises modestly in human (+0.58; +0.95 in male donors). That is above our own 0.5 log2
    zonation floor, so ACAA2 is not flat in human and may reverse direction. Notebook 39 classes it
    as a reversal.
  - Both genes hold against the cortical mouse segment, neither is sex-biased in mouse PT, probes are
    balanced, and both proteins are high in human tubules (Human Protein Atlas). Low human zonation
    therefore does not mean absence.
  - With two specimens per species, our own data leave both indeterminate. The claim rests on the
    external references.
- *Acsm3, Crot and Nudt19 are indeterminate.* Acsm3's late rise is moreover a male-mouse program
  (female microdissection 0.0 / 1.3 / 0.3 TPM). Hadh did not replicate.
- *Literature.* Nothing on these genes. Classical rat microchemistry found 3-hydroxyacyl-CoA
  dehydrogenase activity roughly uniform along cortical segments (Le Hir et al. 1982; Bastin et
  al. 1990), so an activity assay would not see gene-specific zonation. An intravital study
  found that late mouse PT mobilises lipid droplets with lipase-rich lysosomes (Kaminska et al. 2026).

**Drug handling sits in different places** (lead story by coordinator decision, disclosed).

*UGT1A9 is expressed only in human and falls toward late PT.*

- Human: S3c − early is "down" in our data; atlas S3 − S1 = −0.67.
- Mouse: no expression (microdissection ≤ 0.8 TPM in every segment, both sexes; probes balanced).
- Literature: UGT1A9 is among the most abundant renal UGTs in human (Margaillan et al. 2015); its
  axial profile is not reported.

*Mouse places glutathione synthesis late.*

- The glutamate–cysteine ligase subunits are mouse-only zonation: Gclm, and Gclc (probe-imbalanced;
  human GCLC atlas +0.61 against mouse +3.3 and +2.9). The third enzyme, Gss, is conserved.
- Glutathione metabolism is one of the four robust pathways beyond global flattening, and is
  enriched for mouse-only genes.
- Axial glutathione biology is species-dependent. In rabbit PT, synthesis is highest in S1 and
  cellular GSH highest in S3 (Parks et al. 1998). In rat, glutathione is highest in convoluted and
  early straight PT (Brehe et al. 1976).

*Mouse Cyp2e1 is confined to S1–S2.*

- Mouse renal CYP2E1 is a known proximal-tubular, male and testosterone-regulated enzyme (Hu et al.
  1990; Speerschneider et al. 1995).
- Human CYP2E1 mRNA is essentially absent from our donor and the atlas, so its zonation cannot be
  classified. The literature on human protein conflicts (Speerschneider et al. 1995; Akakpo et al.
  2023).

*Decision and its basis.* The pre-specified rule marked this story supporting because Gclc has 3
mouse against 6 human probes. The coordinator promoted it: Gclm, Ugt1a9 and Cyp2e1 are
probe-balanced, and the Gclc imbalance works against the observed mouse-high direction. EPHX1 did
not replicate and is dropped.

*Reading.* A rodent model reproduces neither where human PT conjugates (early) nor where mouse PT
makes glutathione (late) and activates Cyp2e1 substrates (S1–S2).

**Genes whose species difference changes sign** (reversal).

*Dcxr is the clearest case.* It falls from S1 to S3 in mouse and rises in human. All five external
references agree:

- human atlas +1.53 (+1.59 in men, +1.48 in women);
- mouse snRNA −0.46 and −0.54;
- mouse microdissection −1.94 and −2.03;
- rat −1.45.

Detection is high in both species and probes are balanced. DCXR protein is known in mouse PT
(Nakagawa et al. 2002), and DCXR loss tracks human CKD (Perco et al. 2019). Its axial distribution
is not reported in either species.

*Other reversals.* Ugt3a1 reverses in our data and in the atlas, and Vnn1 is a reversal in our
classification.

*Genes that earlier rules called reversals.* Earlier interaction-sign rules flagged more genes, but
the zonation classes reassign several of them:

| Gene | R3 class |
|---|---|
| Pah, Cyp24a1 | Conserved |
| Acox2, Igfbp4, Glyat | Mouse-only |
| Nt5e, Cyp7b1 | Indeterminate |

Only Dcxr and Ugt3a1 are named as reversals in the text.

**Expressed and zonated only in human.**

- *RBP4* rises steeply toward late human PT (atlas +3.46). It also rises in rat, where in situ
  hybridisation placed RBP mRNA in outer-stripe S3 (Makover et al. 1989). Mouse Rbp4 is 0.0 TPM in
  every microdissected segment of both sexes. The same rat study notes that perinephric fat expresses
  RBP mRNA. This is therefore an mRNA claim, and ambient signal from perirenal fat is untested.
- *AOX1* rises toward late human PT (atlas +1.13). It is supporting only, because the panels carry
  2 human against 3 mouse probes. Humans carry one AOX gene where mice express four (Terao et al.
  2016).

**The organic-anion module is mouse-only zonation.**

- OAT1 (Slc22a6) and NaDC3 (Slc13a3) rise from S1 to S2 in mouse and stay flat in human. In mouse,
  OAT1 peaks in the cortical medullary-ray segment (970 TPM).
- An earlier reading, "human keeps the module into late PT", came from comparing human S3 with
  outer-stripe mouse S3. It is withdrawn.
- The module remains the clearest example of a whole-PT null that hides opposite regional
  differences (Figure S6). Human OAT1–3 are known to be stronger in S1/S2 than in S3 (Breljak et al.
  2016).

**Supporting observations.**

- *Sterol synthesis* is mouse-only zonation (Lss, Ebp; Cyp51, probe-imbalanced). It holds in female
  mice and in the cortical segment. Hmgcr is conserved, and rat runs the opposite way, so this is a
  mouse program, not a rodent program.
- *Early-PT transporters in late PT.* Human late PT keeps SLC6A19, SLC9A3 and SLC4A4 more than mouse
  does in the external references. In our classes, Slc6a19 and Slc9a3 are indeterminate and Slc4a4
  is conserved. This stays a hypothesis. Mouse B0AT1's S1–S2 restriction is known (Navarro Garrido et
  al. 2022).

**What did not hold.**

- *Serine synthesis is not human-specific.* Psat1 is conserved, though probe-imbalanced. Mouse Phgdh
  and Psat1 each have one probe against three in human, and mouse microdissection shows both
  S1-restricted, as in human.
- *Not replicated:* Gamt (in our data), Hadh, Ephx1 and Me1.

---

## 5. Discussion

**What is new.** To our knowledge this is the first comparison of human and mouse PT gene programs as
a function of position along the tubule. It yields three results.

1. **A design argument.** With two specimens per species, pathway screens that compare averages
   cannot be separated from specimen variation. Position-dependent screens can, because each specimen
   supplies replicated positions and specimen offsets cancel (split plot; relabeling control).
2. **A biological result.** Beyond a conserved core of apical solute transport and intrarenal
   renin–angiotensin genes, PT zonation is largely species-specific.
   - S1→S2 gradients correlate at only 0.14–0.23 after noise correction.
   - Mouse-only zonated genes (247) outnumber human-only genes (66), with the excess concentrated in
     S1 → S2 metabolism.
   - Human S1→S2 gradients are about half as large overall (0.47–0.58). The region-matched S3
     contrast is reduced less (0.84).
3. **A short list of positional differences, at stated levels of evidence.**
   - *Headline under pre-specified gene-level rules,* tested against both mouse sexes, direct
     microdissection and a cortex-only human atlas:
     - Gatm: mouse-only zonation of creatine synthesis;
     - Gclm: glutamate–cysteine ligase subunit;
     - Dcxr and Ugt3a1: reversals;
     - UGT1A9 and RBP4: expressed and zonated only in human.
   - *Supporting under those rules:* Gclc and AOX1. The glutathione story was promoted to lead after
     the rules were applied (disclosed).
   - *Supported by class enrichment (R3) but not verified gene by gene:* mouse-only sterol synthesis
     (Lss, Ebp) and peroxisomal very-long-chain β-oxidation (Acox1, Hsd17b4).
   - *Headline genes that are indeterminate in our own classification:* Acadm, Acaa2, Acsm3 and
     Nudt19. These rest on the external references.

Rodent zonation of several of these programs was already known from microdissection or enzyme
activity (Takeda et al. 1992; Parks et al. 1998). The human arm and the species contrast are new.

**Relation to prior work.**

- *Human atlases* label PT-S1/S2/S3 as discrete classes of dissociated nuclei (Lake et al. 2023).
- *Mouse microdissection* gives three anatomical points (Chen et al. 2021, 2023).
- *Cross-species atlases* align cell states (Klötzer et al. 2025). They find little overlap in
  single-cell expression changes but agreement at the pathway level (Zhou et al. 2023).
- *The mouse sex atlas* compares whole-PT sex differences with human, not positions (Chen S et al.
  2025).

Our map adds a positional axis measured on intact tubules in both species, and uses the newest atlases
to test our claims rather than as competing maps. The 2025 cross-species atlas is not open access,
and its PT annotation depth could not be checked; this is the main residual novelty risk.

**Why "flattening" is mostly a mouse-centric view.** Measured on mouse-zonated genes, human amplitude
looks like 0.15 of mouse; measured on human-zonated genes, it looks like 2–3. The remainder is a
moderate global reduction. It is about half for S1→S2, and smaller (0.84) for the region-matched
S3 contrast. Injury lowers amplitude in both species and, in our checks, does not reproduce the
human pattern. A contribution of tissue state is not excluded [PENDING D8].

**Implications for mouse models of PT physiology and toxicity.**

- *Transport.* For the apical transporters in the conserved class, mouse predicts where human PT
  places them. This does not hold for every transporter: OAT1 and NaDC3 are mouse-only zonation.
- *Metabolism.* Mouse PT confines several metabolic programs to one region that human PT keeps flat:
  creatine and glutathione synthesis, sterol synthesis and peroxisomal β-oxidation. A metabolic or
  toxic effect confined to one mouse segment may therefore be spread along the human PT, or absent.
- *Drug handling.* Human early PT carries UGT1A9, which mouse PT lacks. Mouse PT activates Cyp2e1
  substrates in S1–S2 and synthesises glutathione late.
- *Injury models* already differ in where injury falls along the nephron (Heyman et al. 2010). Our
  map adds that the metabolic programs themselves sit elsewhere.
- *Sex.* The mouse-only excess is not driven by male-biased genes. Individual programs are, however
  (Acsm3; Xiong et al. 2023), so sex matters for specific genes.
- *Computational models.* Human nephron models borrow rodent axial parameters (Layton et al. 2019).
  The conserved transport core supports that choice for transport, but not for metabolism.
- *Caution.* All of these are mRNA-level observations. Protein evidence for the lead genes is
  presence or absence in human tubules, not zonation.

**Limitations.**

- *Replication.* Two male control mice and two cortex sections from one male human donor. Specimens,
  not structures, are the replicates. Agreement between the two human sections is not donor
  replication, and no statement here is a population estimate for humans.
- *Sex.* All specimens are male, so species and sex are confounded. Female mouse references confirm
  most late mouse programs but not Acsm3, and the mouse confirmation atlas for the zonation classes
  used male donors.
- *Region.* Human sections are cortex only, so mouse outer-stripe S3 has no counterpart. Region
  checks (cortical-like mouse S3; the cortical medullary-ray mouse segment) changed the conclusion
  for the OAT module.
- *Classification.* Most genes (5,590 of 7,407) remain indeterminate with two specimens per species.
  Class sizes depend on rules (Figure 3b), and the human reference has fewer donors than the mouse.
- *Probe panels and processing.* The two species were measured with different panels. Genes with
  unequal probes are reported separately; one apparent species difference, serine synthesis, was a
  panel artefact. Space Ranger versions also differ.
- *Coordinate.* Within-segment order is weakly supported. Human S1 ordering depends on read depth,
  human–mouse registration moves under perturbation, and the curve is numerically sensitive. Its
  value is replication across positions and continuous localisation.
- *External references.* External human positions are snRNA cluster annotations, not anatomy.
  Pathway-level replication is not shown beyond the shared global component.
- *Decisions taken after seeing results (disclosed).*
  1. The pre-registered coordinate-selection rule selected DPT. The pre-registered downstream
     specificity gate (C6 criterion 5) was then applied to every candidate after DPT was seen to
     fail it, and scFates was kept as primary.
  2. The 5% condition was added to the specificity rule after the first run.
  3. NCDR replaced the pre-specified decoy FDP as the primary R2 metric.
  4. The drug-handling story was promoted to lead.
  5. The notebook 39 diagnostics that turned "flattening" into "different genes" were post hoc.
- *Unknown donor context.* The human donor's source, age and tissue state are not recorded here
  [VERIFY].

**What would move this forward.**

- Female animals and female human donors on the same platform.
- Several human donors.
- Matched cortical sampling, or human medulla.
- In situ validation of position-defining genes on human sections where S1 and S3 can be
  distinguished: GATM, GCLM, DCXR, UGT3A1, RBP4 and UGT1A9.

---

## 6. Methods

This section merges `docs/paper/methods_draft.md` (workstream 4), `docs/paper/methods_pathway_draft.md`
(workstream 2) and the corrections in `docs/paper/methods_corrections_ws1.md` (workstream 1). It
adds the zonation amplitude (notebook 39), zonation class (notebook 40) and biology verification
(notebook 38) methods. Parameters were read from the code.

### 6.1 Study design and specimens

**Cohort.**
- Mouse: two control mouse kidney sections, `Ctrl1A2` and `Ctrl1A4`.
- Human: two sections of renal cortex from one donor, `HUK1_COR1` and `HUK1_MED1`. The second is
  named "MED" at source, but its segmentation contains cortex only, and it is analysed as cortex.

All four specimens express Y-linked transcripts and are treated as male [VERIFY against records].
The mouse sections are whole-kidney sections, so they include cortex and outer medulla. The two
ischaemia–reperfusion kidneys of the same mouse experiment (`IR2A2`, `IR2A4`) are used only in the
tissue-state analysis of notebook 39. Specimens are the units of replication.

**Not yet recorded:**
- Mouse: strain (the paper plan says C57BL/6-type; no file records it), age, supplier, housing and
  euthanasia [VERIFY].
- Human: tissue source, donor age, kidney function, and whether the two sections are serial or come
  from separate blocks [VERIFY].
- Tissue processing: fixation and embedding, section thickness, H&E protocol and imaging system
  [VERIFY]. A CytAssist image exists for the human sections.
- Ethics: animal-protocol and human-tissue approvals [VERIFY].

### 6.2 Visium HD processing

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

### 6.3 Tubule segmentation

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
  with D4 test-time augmentation [VERIFY whether to report].
- *Known errors.* Predicted objects are smaller than annotated ones, and merges of touching
  same-type tubules are the main error.

**Human segmentation and mouse quality control.**
- The human `*_v2.geojson` re-segmentations have no recorded provenance [VERIFY].
- Mouse analyses use only the quality-controlled `*_kept_tubules_labeled_fine.geojson`. The rules
  and code of that upstream QC are not in the repository [VERIFY].
- Only polygon geometry is used downstream.

**Centroid check.** Every join between expression and segmentation recomputes the indexed polygon's
centroid and stops at any deviation above 1 px.

### 6.4 Tubule-by-gene matrices

Each 2-µm bin was assigned to the polygon covering its centre, in full-resolution image coordinates
(scale factor 1, shapely STRtree, at most one polygon per bin). Raw counts were summed per polygon.
The segmentation and the Visium HD image must share one coordinate frame [VERIFY how they were
registered].

| Specimen | Polygons | Polygons with ≥ 1 bin |
|---|---:|---:|
| `Ctrl1A2` | 7,624 | 7,624 |
| `Ctrl1A4` | 7,002 | 7,002 |
| `HUK1_COR1` | 8,747 | 6,416 |
| `HUK1_MED1` | 14,186 | 7,223 |

The many human polygons without bins are unexplained [VERIFY].

### 6.5 Cross-species ortholog space

**Ortholog map.** We used the HCOP human–mouse table (Yates et al. 2021) [VERIFY download date;
local sha256 prefix `0cfb78e4eb273751`].

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

### 6.6 Structure filters and normalisation for integration

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

### 6.7 Integration, clustering and reviewed labels (pass 1)

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

Notebook 13 repeats notebook 03's pass-1 clustering on the same 26,839 structures. The two
partitions differ in membership hash but agree closely:
- adjusted Rand index (ARI) 0.979 for the Leiden clusters and 0.981 for the reviewed labels;
- 99.1% of structures carry the same label.

Notebook 13 calls 12,866 structures PT and notebook 03 calls 12,871; 12,864 are PT in both. Of these,
98.5% have the same S1/S2/S3 label (ARI 0.955). Most disagreements are notebook-03 S1 structures that
notebook 13 labels S2 (125).

### 6.8 Nephron re-integration (pass 2) and PT subset

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

### 6.9 PT coordinate

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

### 6.10 Coordinate validation and alternatives (notebooks 35 and 36)

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
- *Count splitting.* Poisson splitting (ε = 0.5) between building the coordinate and testing.
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

**Outcome (disclosed).**
- The rule selected DPT13 (P5 0.80; P1 0.74 mouse, 0.29 human; gap 0.14).
- The scFates equal-depth refit failed (P5 0.73; gap 0.36).
- scFates itself would also have failed (P5 0.72; gap 0.29).
- The rule omitted the pre-registered downstream specificity gate: criterion 5 of rule C6, which
  requires the T_spatial screen to keep NCDR < 0.25 and fewer than 5% of pathways called under
  relabeling.
- After DPT13 was seen to fail that gate (NCDR 0.40), the coordinator applied it to every
  candidate. No candidate then qualified, and scFates was kept as primary.

**Outcomes.**
- *DPT* on notebook 13's embedding: 96 calls; 56 of 66 robust and 38 of 42 core pathways from the
  earlier list retained; NCDR 0.40 (matched test), so not specific. Notebook 37's joint-test
  pipeline is being run end to end on DPT13 and the equal-depth refit [PENDING D4].
- *Conserved-anchor coordinate:* 42 of 66 and 30 of 42 retained, on 7,167 structures, because it
  did not register the species.
- Refits whose curve could not be rooted were reported as failed, not retuned.

### 6.11 PT analysis set and gene universe

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

### 6.12 Nested gene models

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
- The species difference Δ(s) = μ̂_human − μ̂_mouse has HC3 standard errors (MacKinnon & White 1985),
  and Z(s) = Δ/SE.

### 6.13 Pathway libraries and the joint test

**Libraries.**
- Reactome 2022 (Gillespie et al. 2022), MSigDB Hallmark 2020 (Liberzon et al. 2015) and KEGG 2019
  Mouse (Kanehisa et al. 2019), in Enrichr format (Kuleshov et al. 2016) [VERIFY download source and
  date].
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

### 6.14 Specificity: relabeling, NCDR and target–decoy calibration

**Relabeling.** The four specimens have three balanced 2 + 2 splits: the species split and two
relabelings, {`HUK1_COR1` + `Ctrl1A2`} and {`HUK1_COR1` + `Ctrl1A4`} against the rest.
- Every strategy was rerun unchanged on both relabelings, as in Lamian's random-partition null
  (Hou et al. 2023).
- In relabeled gene models, species × spline enters as a nuisance term, and within-group specimen
  contrasts absorb the species offset.
- Relabeled DESeq2 used `~species + group`.
- Two within-species comparisons (mouse vs mouse, human section vs section) were added as controls.

**Why relabeling is a valid null.**
- Write a specimen's deviations as an offset a_ij and a positional deviation b_ij(s), independent
  across specimens.
- The species split and the two relabelings are the three orthogonal ±½ contrasts of the four
  specimens. Under independence they have equal noise variance.
- The design is a split plot (Altman & Krzywinski 2015): offsets have 2 df of error, and positional
  terms 2(k − 1) df.
- A term shared by both human sections cancels in the relabelings and is not covered.

**Summaries.**
- NCDR is the mean number of relabeled calls divided by the number of species calls.
- A strategy is called specific when NCDR < 0.25 and relabelings call < 5% of tested pathways. The
  second condition was added after the first run (disclosed).
- Decoy-estimated FDP (Elias & Gygi 2007) is FDP(t) = (mean relabeled count ≥ t) / (species count
  ≥ t), made monotone as a q-value. We report calls at FDP ≤ 10%, with 300-draw pathway-resampling
  intervals, and the minimum attainable FDP.
- Decoy FDP was the pre-specified primary R2 metric. It was replaced by NCDR because pathway-score
  models reach FDP 0 (disclosed).
- Decoy inflation a and implied ρ were estimated from the slope of mean decoy z² on pathway size.

### 6.15 Conventional comparators

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

### 6.16 Building the reporting set

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
- the earlier analysis on notebook 03's DPT coordinate (notebook 12 logic v5; 11,106 genes, 1,516
  pathways). This is a historical sensitivity only and is not cited as robustness support
  [PENDING D4];
- the specimen-level split-plot test.

### 6.17 Programs and descriptors

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
- *Earlier display rule:* ≥ 3 contiguous grid points with |Z| ≥ 3 and |Δ| ≥ 0.25 on each side.
- *Segment DESeq2 confirmation:* opposite significant signs.

### 6.18 Specimen-level split-plot test

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

### 6.19 Zonation amplitude (notebook 39; protocol pre-registered)

**Gradients.**
- Pseudobulk log2 CPM ratios of summed counts per specimen and reviewed segment, computed over
  orthologs measured in all specimens.
- Contrasts: S2 − S1 (primary, cortical in both species); S3c − early, with mouse S3 restricted to
  cortical-like structures; and S3 − early.
- Genes needed ≥ 20 counts in every specimen × segment pseudobulk of both species, and a reference
  mean ≥ 3 log2 in one compared group.

**Dilution-free ratio.**
- Both of our species gradients were regressed on an independent reference gradient. The ratio is
  cov(human, ref) / cov(mouse, ref).
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
2. Lake donor injury scores against donor amplitude. Scores are within-dataset z-scores of
   injury-marker pseudobulks [VERIFY marker list from notebook 39], and donor amplitude is the slope
   on male microdissection.
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

### 6.20 Zonation classes (notebook 40)

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
annotation. This check was added by the coordinator after the classes were computed.

### 6.21 External segment-resolved datasets and verification of lead genes

**Datasets.**

| Dataset | Content |
|---|---|
| GSE150338 (Chen et al. 2021) | Microdissected male mouse PTS1–3; TPM replicates |
| GSE212213 (Chen et al. 2023) | Microdissected mouse PTS1–3, both sexes; mean TPM |
| GSE56743 (Lee et al. 2015) | Microdissected rat S1–S3; RPKM; reported only |
| CELLxGENE Census 2025-11-08 (CZI Cell Science Program et al. 2025), mouse dataset `25818bf7…` (Chen S et al. 2025) | 12 male and 12 female donors with ≥ 50 nuclei per segment |
| Census human renal-cortex dataset `09b518f9…` (Acera-Mateos et al. 2026) | 6 donors; convoluted PT vs S3 |
| Lake et al. (2023) snRNA v1.5 (CELLxGENE `a12ccb9b…`) | Cortex nuclei (`region C`, cortex tissue), healthy donors, author labels PT-S1/S2/S3. 7 donors with ≥ 50 nuclei per segment (4 female, 3 male); 12 at ≥ 20 nuclei as sensitivity. AKI and CKD donors were used for tissue state. |

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

### 6.22 Literature check

Each confident pathway (with up to eight driver genes), each direction-reversing gene, and each
pathway called only by conventional screens was searched.

- *Search.* Titles and abstracts in Europe PMC, using only gene, pathway and generic terms. Key
  papers were read in open-access full text where available.
- *Status.*
  - KNOWN: species difference and PT position reported.
  - PARTLY KNOWN: one part reported.
  - NOT FOUND: nothing after ≥ 2 targeted queries. This is a candidate, not proof of novelty.
  - CONTRADICTS: the literature reports the opposite.
- *Verification.* Workstream 3's phase 2 re-verified each citation used in R5. The search was done
  with AI agents [VERIFY wording for the paper].

### 6.23 Statistics and multiple testing

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
- Two relabelings screen only for gross nonspecificity.

**Coordinate limits.**
- Human S1 ordering depends on read depth: equal-depth agreement is 0.36–0.46, and the
  library-size correlation falls from −0.31 and −0.38 to −0.04 and −0.13 after thinning.
- Human placement moves by up to 0.29 when a third of the genes or half the reads are removed.

### 6.24 Software and reproducibility

**Analysis environment.**
- Python 3.11.16.
- Scanpy 1.11.5, AnnData 0.12.19, NumPy 1.26.4, SciPy 1.17.1, pandas 2.3.3, statsmodels 0.15.0,
  patsy 1.0.3, scikit-learn 1.9.0, python-igraph 1.0.0, leidenalg 0.12.0.
- scFates 1.2.5, PyDESeq2 0.5.4, GSEApy 1.3.1, shapely 2.1.2, rpy2 3.6.7.
- R 4.5.3 with harmony 2.0.5.
- `cellxgene_census` 1.18.0 for the Census extraction (committed as
  `analysis/scripts/fetch_census_pt_segments.py`).

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
- Notebook 40 is not yet committed [VERIFY].

### 6.25 Data and code availability

[VERIFY: deposition of the Visium HD data and segmentation polygons (controlled access for human
tissue), code archive DOI and licence.] The public datasets are listed in 6.21.

---

## 7. Figures and tables

Every figure must come from a committed, output-free notebook. Sources are given relative to
`results/`.

### Main figures

**Figure 1 | Study design and a continuous cross-species PT coordinate (R1).**
- **a,** Design. Visium HD on two control mouse kidney sections and two cortex sections from one
  human donor; panoptic segmentation of tubule cross-sections from H&E; aggregation of 2-µm bins
  into tubule-by-gene profiles; reviewed clustering; the PT subset.
- **b,** Pass-2 PT embedding with the nonbranching scFates curve, coloured by reviewed segment and
  split by species. Inset: specimen segment medians and transitions on the coordinate.
- **c,** Coordinate agreement with external zonation:
  - per-gene, fold-protected within-specimen correlation with the coordinate, compared across genes
    with external S3 − early contrasts (mouse male snRNA, mouse microdissection, human cortex snRNA
    donor half B);
  - segment labels shown as a reference;
  - within-segment agreement (P1w), with read-depth ordering as a null-like reference.
- **d,** Identifiability checks:
  - absorption ratio for an injected human-only gradient (1.01);
  - sign agreement for 244 externally species-specific genes (0.91);
  - count-split arms (robust and core pathways of the earlier list retained).

*Sources:*
- Notebooks 01 and 13: `minimal_pt_scfates/pt_scfates_dpt_embedding.pdf`,
  `pt_segment_order_scfates_dpt.pdf`.
- Notebook 35: `pt_reconstruction_v2/nb35/figures/p1_p1w_external_concordance.*`,
  `p2_absorption.csv`, `positive_control_recovery.csv`, `g6_count_split.csv`,
  `species_registration_across_refits.csv`.

[VERIFY: no committed notebook yet assembles Figure 1 as one figure; panel a needs a schematic.]

**Figure 2 | Only position-dependent pathway screens are specific with two specimens per species
(R2).**
- **a,** The three balanced contrasts of four specimens (species split and two relabelings) and the
  split-plot error terms.
- **b,** Pathways called (BH q ≤ 0.05) for the species split (filled; orange = specific) and both
  relabelings (open), per strategy. NCDR, the null rate and calls at decoy FDP ≤ 10% are shown at
  right.
- **c,** Decoy-estimated FDP as the threshold is lowered, for average-level and positional
  strategies.
- **d,** Within-species controls. Offset (T_level) and smooth positional (T_spatial) calls for the
  species split, mouse versus mouse, and human section versus section.
- **e,** Mean decoy z² against pathway size. Inflation a and implied correlation ρ for T_level and
  T_spatial; DESeq2 gene calls shown for reference.

*Sources:* notebook 37, `pt_pathway_final/scfates/figures/fig2_specificity.*`; tables
`table_S1_strategy_specificity.csv`, `within_species_controls.csv`,
`decoy_inflation_by_statistic.csv`.

**Figure 3 | Beyond a conserved transport core, PT zonation is largely species-specific (R3).**
- **a,** Per-gene S1 → S2 zonation in each species: mean of the two mice against mean of the two
  human sections (log2). The panel shows the 7,391 probe-balanced genes measurable in both
  species. Colour gives the S2 − S1 class: mouse-only 219, human-only 44, conserved 65, reversal 20
  (the same gene set as `class_sizes_by_contrast.csv`). Calls are noise-scaled and
  reference-confirmed.
- **b,** Class sizes under the primary rule and seven alternatives: no confirmation; no effect floor;
  amplitude-scaled floor; flat margin 0.5; male microdissection confirmation; female snRNA
  confirmation; female microdissection confirmation. Probe-balanced genes,
  S2 − S1 and S3c − early combined.
- **c,** Representative genes chosen by rule, in each specimen at S1, S2 and cortical S3 (log2 CPM;
  orange human sections, blue mice):
  - conserved: Slc5a10, Tnfsf10;
  - reversal: Vnn1, Ugt3a1;
  - mouse-only: Slc5a12, Slc8a1;
  - human-only: Sel1l3, Tiparp.
- **d,** Programs enriched in each class, collapsed by gene overlap. −log10 q from the joint matched
  test; ✓ marks enrichment that also holds with 5-bin matching and without low-expression genes.
- **e,** The robust position-dependent pathways by the class of their classified members (★ core).

*Sources:* notebook 40 (uncommitted at assembly [VERIFY]),
`pt_zonation_classes/figures/fig_R3_zonation_classes.*`; tables
`table_R3_zonation_classes.csv` (Table 1), `class_sizes_by_setting.csv`, `class_programs.csv`,
`robust_pathways_by_class.csv`.

**Figure 4 | Robustness and agreement with independent data (R4).**
- **a,** Effect retained (AUC effect ÷ original) for the 36 robust pathways under each check, with
  the share ≥ 0.7 at right. The DPT row is to be replaced by the notebook 37 run on DPT
  [PENDING D4]. Core pathways dark, robust-only light. The green bar marks power-matched
  random removal of mouse S3.
- **b,** Reviewed segment transitions on the coordinate, per specimen.
- **c,** Gene-level agreement. Our species × segment interaction against Census (human − male
  mouse)(S3 − early), centred. Top 10% T_spatial genes in black; ρ compared with expression-matched
  random genes.
- **d,** Pathway-level external replication z, with donors as replicates. Robust and uncalled
  pathways against Census and Lake/KPMP; relabel-called pathways against Census.
- **e,** Specimen-level split-plot gene test. Q–Q plot for the species split and both relabelings
  (relabeling λ ≤ 0.87).

*Sources:* notebook 37, `pt_pathway_final/scfates/figures/fig4_robustness_replication.*`; tables
`table_S2_pathways.csv`, `external_replication_summary.csv`, `split_plot_diagnostics.csv`.

**Figure 5 | Lead genes against direct segment measurements (R5).**
- *Left,* fitted human (orange) and mouse (blue) curves along the PT coordinate, with each
  specimen's binned means. Shading marks the reviewed S1/S2/S3 ranges.
- *Right,* each external dataset's S1, S2 and S3 values, centred on its own S1 (log2): Lake human
  cortex (7 donors), mouse snRNA by sex, mouse microdissection by sex, and rat microdissection.
- *Panel titles* give story, gene, pre-specified verdict and probes (mouse/human).

The current version shows Gatm, Gamt, Acadm, Acaa2, Gclc, Ugt1a9, Ephx1, Cyp2e1 and Dcxr.
[Coordinator/WS3: add each gene's R3 class to the titles. Consider replacing the two non-replicating
genes (Gamt, Ephx1) with Rbp4 and Ugt3a1 or Gclm, and moving Gamt and Ephx1 to Figure S5.]

*Sources:* notebook 38, `pt_literature_deep/figures/fig11_r5_story_genes.*`; tables
`story_gene_evidence.csv`; R3 classes from `pt_zonation_classes/tables/lead_gene_classes.csv`.

### Main table

**Table 1 | Zonation classes.** For each class:
- genes under the primary rule (probe-balanced);
- range across alternatives;
- counts by contrast;
- enriched pathways and programs;
- top programs and representative genes.

*Source:* `pt_zonation_classes/tables/table_R3_zonation_classes.csv`.

### Supplementary figures and tables

| Item | Content | Source |
|---|---|---|
| Figure S1 | Coordinate limitations: (a) equal-depth agreement and library-size correlation by specimen × segment; (b) species registration gap across refits; (c) gene-fold within-segment agreement and the numerical floor; (d) held-out gene prediction against segment and coverage oracles; (e) attenuation decomposition and noise matching | `pt_reconstruction_v2/nb35/p5_exposure_invariance.csv`, `species_registration_across_refits.csv`, `p3a_*.csv`, `numerical_floor.csv`, `figures/a3_attenuation.*`; `nb36/p4_heldout_prediction.csv` [VERIFY: no assembled figure yet] |
| Table S1 | Coordinate arm comparison and decision rules (SCF13, SCF13-ED, DPT13, SCF-PT, CAC), including the post-hoc specificity gate | `nb36/arm_metrics.csv`, `c6_decision.csv`, `selection_rule.csv`, `downstream_summary.csv` |
| Figure S2 | DPT specificity diagnostic: relabeled-call categories; within-species calls (DPT mouse 36, human 5; scFates 16, 0) | `pt_reconstruction_v2/dpt_diagnostic/*.csv` [VERIFY: no assembled figure yet] |
| Table S2 | Strategy specificity: species and relabeled calls, NCDR, null rate, decoy calls and interval, minimum FDP | `pt_pathway_final/scfates/tables/table_S1_strategy_specificity.csv` |
| Figure S3 | Position-dependent pathways: (a) joint-test funnel; (b) median member Δ(s) for the 36 robust pathways in 17 programs, with peak ticks; (c) peak position against median member Z (human-high early, mouse-high late); (d) the six largest programs with a leading driver gene; (e) averaging loss | `pt_pathway_final/scfates/figures/fig3_positional_programs.*`, `fig3e_averaging_loss.*` |
| Table S3 | All 1,513 pathways: joint test, specificity margin, every robustness ratio, stability, program, shape, external replication | `table_S2_pathways.csv` |
| Table S4 | The 17 programs with drivers and member pathways | `table_S4_programs.csv` |
| Table S5 | External replication summary and global component | `table_S3_external_replication.csv`, `table_S5_global_component.csv`, `external_replication_summary.csv` |
| Figure S4 | Zonation amplitude: (a) human ÷ mouse ratio against each reference (OLS, Deming, orthogonal); (b, c) both species on male microdissection and on Lake human cortex; (h) reference-free noise-corrected SD ratios; (i) disattenuated cross-species correlations (post hoc); flattening checks (tissue state, measurement) | `pt_zonation_amplitude/figures/fig_amplitude_ratio.*`, `fig_reference_free_amplitude.*`, `fig_flattening_checks.*`; `tables/*.csv` |
| Figure S5 | Lead-gene verification matrix: one row per claimed gene, one column per pre-specified test, with verdict, probe counts and annotations | `pt_literature_deep/figures/fig12_r5_replication_matrix.*`; `story_verdicts.csv` |
| Figure S6 | Opposite regional differences that a whole-PT comparison averages away (Slc22a6, Slc13a3, Cyp24a1, Slc22a7). Counts in the source notebook refer to the superseded list. | `pt_pathway_registration_sensitivity/figures/fig9_opposite_regional_differences.*` |
| Figure S7 | Segment-label external validation (notebook 34): zonation within each species and species × segment against external data | `pt_external_segment_validation/figures/fig10_external_segment_validation.*` |
| Table S6 | Probe counts per gene in both panels; zonation classes for all genes, including probe-imbalanced and one-species-expressed genes | `pt_literature_deep/probe_counts.csv`; `pt_zonation_classes/tables/gene_zonation_classes.csv`, `one_species_expressed_genes.csv` |
| Table S7 | Segmentation performance (PQ, SQ, DQ, merge and split rates) | `segmentation/model_summary_v2.txt` [VERIFY: not produced by a notebook] |

---

## 8. References

All references were checked in Europe PMC, Crossref or on the publisher page by workstream 4. Those
marked † were also verified by workstream 3 (`SCRATCH/ws3/literature_phase2.md`).

- Acera-Mateos M et al. Systematic evaluation of single-cell multimodal data integration enhances cell type resolution and discovery of clinically relevant states in complex tissues. *Genome Biol* 27:64 (2026). https://pubmed.ncbi.nlm.nih.gov/41821037/ (Census collection cites the preprint https://doi.org/10.1101/2025.03.06.637075.)
- Akakpo JY et al. Lack of mitochondrial Cyp2E1 drives acetaminophen-induced ER stress-mediated apoptosis in mouse and human kidneys. *Toxicology* 500:153692 (2023). https://pubmed.ncbi.nlm.nih.gov/38042273/ †
- Albergante L et al. Robust and Scalable Learning of Complex Intrinsic Dataset Geometry via ElPiGraph. *Entropy* 22:296 (2020). https://pubmed.ncbi.nlm.nih.gov/33286070/
- Altman N, Krzywinski M. Split plot design. *Nat Methods* 12:165–166 (2015). https://pubmed.ncbi.nlm.nih.gov/25879095/
- Basit A et al. Kidney Cortical Transporter Expression across Species Using Quantitative Proteomics. *Drug Metab Dispos* 47:802–808 (2019). https://pubmed.ncbi.nlm.nih.gov/31123036/
- Bastin J et al. Postnatal development of oxidative enzymes in various rat nephron segments. *Am J Physiol* 259:F895–F901 (1990). https://pubmed.ncbi.nlm.nih.gov/2260682/ †
- Benjamini Y, Hochberg Y. Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. *J R Stat Soc B* 57:289–300 (1995). https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
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
- Kim J et al. Intra-renal slow cell-cycle cells contribute to the restoration of kidney tubules injured by ischemia/reperfusion. *Anat Cell Biol* 44:186–193 (2011). https://pubmed.ncbi.nlm.nih.gov/22025970/
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
- Love MI et al. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol* 15:550 (2014). https://pubmed.ncbi.nlm.nih.gov/25516281/
- Maass C et al. Translational Assessment of Drug-Induced Proximal Tubule Injury Using a Kidney Microphysiological System. *CPT Pharmacometrics Syst Pharmacol* 8:316–325 (2019). https://pubmed.ncbi.nlm.nih.gov/30869201/
- MacKinnon JG, White H. Some heteroskedasticity-consistent covariance matrix estimators with improved finite sample properties. *J Econometrics* 29:305–325 (1985). https://doi.org/10.1016/0304-4076(85)90158-7
- Makover A et al. Localization of retinol-binding protein messenger RNA in the rat kidney and in perinephric fat tissue. *J Lipid Res* 30:171–180 (1989). https://pubmed.ncbi.nlm.nih.gov/2469758/ †
- Margaillan G et al. Quantitative profiling of human renal UDP-glucuronosyltransferases and glucuronidation activity. *Drug Metab Dispos* 43:611–619 (2015). https://pubmed.ncbi.nlm.nih.gov/25650382/ †
- Marshall JL et al. High-resolution Slide-seqV2 spatial transcriptomics enables discovery of disease-specific cell neighborhoods and pathways. *iScience* 25:104097 (2022). https://pubmed.ncbi.nlm.nih.gov/35372810/
- Maunsbach AB et al. Aquaporin-1 water channel expression in human kidney. *J Am Soc Nephrol* 8:1–14 (1997). https://pubmed.ncbi.nlm.nih.gov/9013443/
- Melo Ferreira R et al. Integration of spatial and single-cell transcriptomics localizes epithelial cell-immune cross-talk in kidney injury. *JCI Insight* 6:e147703 (2021). https://pubmed.ncbi.nlm.nih.gov/34003797/
- Moor AE et al. Spatial Reconstruction of Single Enterocytes Uncovers Broad Zonation along the Intestinal Villus Axis. *Cell* 175:1156–1167 (2018). https://pubmed.ncbi.nlm.nih.gov/30270040/
- Muzellec B et al. PyDESeq2: a python package for bulk RNA-seq differential expression analysis. *Bioinformatics* 39:btad547 (2023). https://pubmed.ncbi.nlm.nih.gov/37669147/
- Nakagawa J et al. Molecular characterization of mammalian dicarbonyl/L-xylulose reductase and its localization in kidney. *J Biol Chem* 277:17883–17891 (2002). https://pubmed.ncbi.nlm.nih.gov/11882650/ †
- Navarro Garrido A et al. Aristolochic acid-induced nephropathy is attenuated in mice lacking the neutral amino acid transporter B0AT1 (Slc6a19). *Am J Physiol Renal Physiol* 323:F455–F467 (2022). https://pubmed.ncbi.nlm.nih.gov/35979966/ †
- Nigam SK et al. The organic anion transporter (OAT) family: a systems biology perspective. *Physiol Rev* 95:83–123 (2015). https://pubmed.ncbi.nlm.nih.gov/25540139/
- Nitzan M et al. Gene expression cartography. *Nature* 576:132–137 (2019). https://pubmed.ncbi.nlm.nih.gov/31748748/
- Oliveira MF et al. High-definition spatial transcriptomic profiling of immune cell populations in colorectal cancer. *Nat Genet* 57:1512–1523 (2025). https://pubmed.ncbi.nlm.nih.gov/40473992/
- Oquab M et al. DINOv2: Learning Robust Visual Features without Supervision. arXiv:2304.07193 (2023). https://arxiv.org/abs/2304.07193
- Parks LD et al. Heterogeneity of glutathione synthesis and secretion in the proximal tubule of the rabbit. *Am J Physiol* 274:F924–F931 (1998). https://pubmed.ncbi.nlm.nih.gov/9612330/ †
- Perco P et al. Identification of dicarbonyl and L-xylulose reductase as a therapeutic target in human chronic kidney disease. *JCI Insight* 4:128120 (2019). https://pubmed.ncbi.nlm.nih.gov/31217356/ †
- Ransick A et al. Single-Cell Profiling Reveals Sex, Lineage, and Regional Diversity in the Mouse Kidney. *Dev Cell* 51:399–413 (2019). https://pubmed.ncbi.nlm.nih.gov/31689386/
- Reichold M et al. Glycine Amidinotransferase (GATM), Renal Fanconi Syndrome, and Kidney Failure. *J Am Soc Nephrol* 29:1849–1858 (2018). https://pubmed.ncbi.nlm.nih.gov/29654216/ †
- Saelens W et al. A comparison of single-cell trajectory inference methods. *Nat Biotechnol* 37:547–554 (2019). https://pubmed.ncbi.nlm.nih.gov/30936559/
- Smyth GK. Linear models and empirical bayes methods for assessing differential expression in microarray experiments. *Stat Appl Genet Mol Biol* 3:Article3 (2004). https://pubmed.ncbi.nlm.nih.gov/16646809/
- Speerschneider P et al. Renal tumorigenicity of 1,1-dichloroethene in mice: the role of male-specific expression of cytochrome P450 2E1. *Toxicol Appl Pharmacol* 130:48–56 (1995). https://pubmed.ncbi.nlm.nih.gov/7839370/ †
- Squair JW et al. Confronting false discoveries in single-cell differential expression. *Nat Commun* 12:5692 (2021). https://pubmed.ncbi.nlm.nih.gov/34584091/
- Subramanian A et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *PNAS* 102:15545–15550 (2005). https://pubmed.ncbi.nlm.nih.gov/16199517/
- Takeda M et al. Intranephron distribution of glycine-amidinotransferase activity in rats. *Ren Physiol Biochem* 15:113–118 (1992). https://pubmed.ncbi.nlm.nih.gov/1378964/ †
- Terao M et al. Structure and function of mammalian aldehyde oxidases. *Arch Toxicol* 90:753–780 (2016). https://pubmed.ncbi.nlm.nih.gov/26920149/ †
- Thakur A et al. Sex and the Kidney Drug-Metabolizing Enzymes and Transporters: Are Preclinical Drug Disposition Data Translatable to Humans? *Clin Pharmacol Ther* 116:235–246 (2024). https://pubmed.ncbi.nlm.nih.gov/38711199/
- Traag VA et al. From Louvain to Leiden: guaranteeing well-connected communities. *Sci Rep* 9:5233 (2019). https://pubmed.ncbi.nlm.nih.gov/30914743/
- Vallon V et al. SGLT2 mediates glucose reabsorption in the early proximal tubule. *J Am Soc Nephrol* 22:104–112 (2011). https://pubmed.ncbi.nlm.nih.gov/20616166/
- van der Walt S et al. scikit-image: image processing in Python. *PeerJ* 2:e453 (2014). https://pubmed.ncbi.nlm.nih.gov/25024921/
- Vrhovac I et al. Localizations of Na+-D-glucose cotransporters SGLT1 and SGLT2 in human kidney and of SGLT1 in human small intestine, liver, lung, and heart. *Pflugers Arch* 467:1881–1898 (2015). https://pubmed.ncbi.nlm.nih.gov/25304002/
- Wolf FA et al. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol* 19:15 (2018). https://pubmed.ncbi.nlm.nih.gov/29409532/
- Wu D, Smyth GK. Camera: a competitive gene set test accounting for inter-gene correlation. *Nucleic Acids Res* 40:e133 (2012). https://pubmed.ncbi.nlm.nih.gov/22638577/
- Xiong L et al. Direct androgen receptor control of sexually dimorphic gene expression in the mammalian kidney. *Dev Cell* 58:2338–2358 (2023). https://pubmed.ncbi.nlm.nih.gov/37673062/
- Yates B et al. Updates to HCOP: the HGNC comparison of orthology predictions tool. *Brief Bioinform* 22:bbab155 (2021). https://pubmed.ncbi.nlm.nih.gov/33959747/
- Zhou J et al. Unified Mouse and Human Kidney Single-Cell Expression Atlas Reveal Commonalities and Differences in Disease States. *J Am Soc Nephrol* 34:1843–1862 (2023). https://pubmed.ncbi.nlm.nih.gov/37639336/
- Zimmerman KD et al. A practical solution to pseudoreplication bias in single-cell studies. *Nat Commun* 12:738 (2021). https://pubmed.ncbi.nlm.nih.gov/33531494/

---

## 9. Remaining open items

### A. [VERIFY]: facts that need lab records

1. **Mouse specimens.** Strain, age, supplier, housing, euthanasia and animal approval. Sex was
   inferred from Y-linked genes only.
2. **Human donor.** Tissue source, age, kidney function, sex from records, consent and IRB. Whether
   the two sections are serial or come from separate blocks. Tissue state, which matters for R3's
   amplitude.
3. **Tissue and assay.** Fixation, section thickness, H&E protocol, scanner, the mouse Visium HD
   workflow, library preparation and sequencing depth.
4. **Segmentation provenance.** How the human `v2` segmentations were made. Rules and code of the
   upstream mouse QC (`kept_tubules_labeled_fine`). How the H&E images were registered to the
   Visium HD frame. Why many human polygons have no bins. Whether to report segmentation PQ, and
   versions of the segmentation inference environment.
5. **Reference inputs.** HCOP download date. Source and date of the pathway libraries. Data and code
   availability, licence and deposition.

### B. [VERIFY]: numbers or statements not backed by a committed document or table

6. **Resolved (v0.2).** The attenuation bound on the 36 robust pathways is about 15% of the peak gap
   and 7–8% of the mouse-high count (workstream 1 final checks).
7. **Resolved (notebook 40 v2).** The class table now has all 36 robust pathways: 19 mouse-only, 15
   conserved, 1 reversal, 1 none. The two pathways below the 10-member floor are kept with a note.
8. **Resolved (notebook 40 v2).** Figure 3a now shows the 7,391 probe-balanced genes, with counts
   219 / 44 / 65 / 20.
9. **The coordinator's R3 note on mouse sex is wrong.** `docs/paper/r3_draft.md` says that the mouse
   snRNA confirmation atlas pools male and female donors. Notebook 40's code uses the 12 male donors
   only (`male_sn`). I wrote "12 male donors"; please correct r3_draft.md.
10. **Resolved (v0.2).** Partition concordance (ARI 0.979 / 0.981; 98.5% of S1/S2/S3 labels agree)
    and the glomerulus-source check (22 of 1,581) are now in Methods.
11. **Lake injury-marker list** used in notebook 39. Notebook 39 gives ρ = −0.24 and −0.41, while
    workstream 3's separate check gives −0.42 on another amplitude measure. Report one, and say which
    metric it uses.
12. **Uncommitted notebook.** Notebook 40 is not committed. Figures 1, S1, S2 and Table S7 are not
    yet assembled by a committed notebook.

### C. [CITE]

13. **OpenMidnight** has only a model card and blog post (Kaplan et al. 2025); choose the citation
    form. The three [CITE] gaps in the earlier Introduction are now closed:
    - drug injury: Maass et al. 2019;
    - injury distribution across models: Heyman et al. 2010;
    - the continuous-zonation sentence was rephrased so that it no longer makes a factual claim.

### D. Coordinator decisions to keep disclosed in the paper

14. The pre-registered coordinate-selection rule selected DPT13. The coordinator applied the
    downstream specificity gate (C6 criterion 5) to every candidate after DPT13 was seen to fail it,
    and scFates was kept (R1, Methods 6.10, Limitations).
15. The 5% condition was added to the specificity rule after the first run.
16. NCDR replaced the pre-specified decoy FDP as the primary metric.
17. The drug-handling story was promoted to lead against the pre-specified rule.
18. The notebook 39 diagnostics are post hoc.
19. The R3 sex-bias check was added after the classes were computed.

### E. Consistency items for workstreams 2 and 3

20. **Earlier reversal list.** The draft now uses notebook 40 classes throughout (R5 table). Only
    Dcxr, Ugt3a1 and Vnn1 are named as reversals; Acaa2 is noted as a possible reversal (notebook 39).
    `pt-pathway-literature-novelty.md` (Tier 1) and `r5_draft.md` still carry the old labels for Pah,
    Acox2, Igfbp4, Nt5e, Glyat and Cyp24a1 (outside this file).
21. **OAT module.** The "human retains it into late PT" reading in the novelty doc is withdrawn: the
    module is mouse-only zonation in notebook 40.
22. **Figure 5 genes.** Add R3 classes to the titles. Consider replacing Gamt and Ephx1, which do not
    replicate, with Rbp4, Ugt3a1 or Gclm.
23. **Resolved in this draft.** R1 now uses 10 / 26 (median peaks 0.16 / 0.63). `docs/paper/r1_draft.md`
    still says 19 / 47.

### F. Pending revision analyses (sentences marked in the text)

24. **[PENDING D4].** Notebook 37 end to end on DPT13 and the equal-depth refit. Affects R1
    (Consequence), R4 (DPT robustness), Methods 6.10 and 6.16, and the Figure 4a DPT row.
25. **[PENDING D8].** Tissue state of the human sections, and its contribution to the amplitude
    pattern. Affects R3 (Log compression and tissue state) and the Discussion.
