# Manuscript draft v0.3: where human and mouse proximal tubule differ along the tubule axis

**Version 0.3.** This version integrates the revision analyses of notebooks 41–46 and the
coordinator's decisions in `v02_changes.md` (items 1–14). It answers the internal referee report
(`docs/paper/referee_report_v01.md`); the point-by-point response is in `referee_response_v01.md`
and the change log in `draft_v03_changelog.md`.

Every number was checked against the committed documents (`docs/results/`, `docs/paper/`) or the
saved result tables they cite. Correctness came before length; the draft has not yet been shortened.

**Markers.**
- **[VERIFY]**: needs lab records, or a number that no committed document or table supports.
- **[CITE]**: needs a reference.
- **(disclosed)**: a decision taken after results were seen. It must stay in the paper.

**Citations.** Author–year in the text. Section 8 gives PubMed or DOI links, and every entry was
checked in Europe PMC, Crossref or on the publisher page. Section 9 lists what remains open.

**Standing caveat.** Two control mice and two cortex sections from one male human donor; all
specimens male. The human section labelled "medulla" is healthy cortex, so mouse outer-stripe S3 has
no human counterpart. The Visium HD probe panels differ between species. All findings are
descriptive. A short form of this caveat sits beside each result.

---

## 1. Title options

1. Beyond a conserved core of strong markers, human and mouse proximal tubule zonate largely
   different genes
2. Human and mouse proximal tubule share the zonation of strong transport markers and little else
3. A positional comparison of human and mouse proximal tubule from segmented spatial transcriptomics

Option 4 of v0.2 ("Mouse … zonates metabolic programs that human … keeps flat") is withdrawn. It
overclaims, because most metabolic leads are indeterminate in our own data.

## 2. Abstract (≈ 270 words)

The proximal tubule (PT) performs most renal reabsorption, drug secretion and metabolism, and much
of its biology is studied in mice. Whether human and mouse PT arrange these functions alike along the
S1→S3 axis has not been measured directly. We profiled two control mouse kidney sections and two
cortical sections from one human donor with Visium HD, segmented tubule cross-sections from the
paired histology, and ordered 12,866 PT structures along a cross-species coordinate.

With two specimens per species, average-level pathway screens called as many pathways after
relabeling the specimens into mixed-species groups as for the species contrast, in our cohort and in
100 cross-species designs drawn from public atlases. Position-dependent screens were specific only
when the specimens within a group shared their positional biology.

We compared gene gradients against same-species, cross-dataset ceilings. Relative to these ceilings
(1 = as similar as two datasets of one species), human–mouse conservation was 0.17 (95% CI
0.12–0.22) for S1→S2 and 0.28 (0.22–0.33) for S3 against early PT. It was 0.27–0.28 in public atlases
alone, and at most 0.31 with human segment labels transferred from a reference atlas. Of 388 genes
classifiable as zonated in either species, 196 were zonated alike and 75 in opposite directions; only
41 were confidently mouse-only and 76 human-only. Human S1→S2 gradients were also smaller
(noise-corrected spread 0.58 of mouse). Of 36 robust position-dependent pathways, 26 summarised mouse
zonation that the human sections did not share.

With one human donor and weakly supported human S2 labels, these results are descriptive. They
suggest that, beyond a conserved core of strong markers, mouse PT is an incomplete guide to how human
PT zonates its genes.

---

## 3. Introduction

**Translation of PT biology depends on knowing where species differ.** The PT reabsorbs about
60–70% of filtered water and NaCl, more of the filtered bicarbonate and nearly all filtered
nutrients. It also secretes solutes, produces hormones and carries out most of the kidney's
metabolism (Curthoys & Moe 2014). During fasting and stress the kidney supplies up to 40% of
systemic glucose production, and PT gluconeogenesis is impaired in acute kidney injury (AKI)
(Legouis et al. 2020). Defective fatty-acid oxidation in tubular epithelium contributes to fibrosis
(Kang et al. 2015). The SLC22 organic anion transporters handle common drugs, toxins, nutrients and
putative uremic toxins (Nigam et al. 2015).

The PT is regarded as the primary target of injury in AKI and chronic kidney disease (Chevalier
2016), and as the primary site of drug-induced kidney injury (Maass et al. 2019). Much of this
biology is studied in rodents. No intervention has yet prevented AKI in patients, which has raised
the question of whether animal models predict human responses (de Caestecker et al. 2015). Rodent
AKI models do not fully recapitulate the human disease (Fu et al. 2024), and experimental models
differ in where along the nephron tubular injury falls (Heyman et al. 2010). Preclinical models also
often fail to predict drug nephrotoxicity (Maass et al. 2019). Species differences in renal drug
transporters and metabolising enzymes contribute to poor translation (Basit et al. 2019; Thakur et al.
2024).

**The PT is zonated, and some of its zonation is known to be conserved.** Along its length the PT
is divided into S1, S2 and S3, and functions are distributed along this axis:
- in mouse, SGLT2 sits in the early PT (Vallon et al. 2011);
- in isolated human PT segments, S2 and S3 make more glucose from lactate than S1 (Conjard et al.
  2001);
- after ischaemia in mice, damage is most severe in S3 of the outer stripe (Kim et al. 2011).

Segment-resolved transcriptomes exist for microdissected rat (Lee et al. 2015) and mouse tubules
(Chen et al. 2021). Single-cell atlases describe PT diversity in mouse (Ransick et al. 2019) and
human (Lake et al. 2023).

Human protein maps place SGLT2 in S1/S2 and SGLT1 in S3, much as in rodents (Vrhovac et al. 2015).
They place AQP1 highest in the proximal straight tubule (Maunsbach et al. 1997). They also place all
three OATs more intensely in S1/S2 than in medullary-ray S3, with outer-stripe S3 unstained (Breljak
et al. 2016). No study has compared these maps position by position with rodent data. Our transcript
calls agree with this map for SGLT1/2 but not for OAT1 and OAT2 (Results).

Other programs differ. In mouse, sexually dimorphic gene activity maps predominantly to PT segments,
but only a limited set of genes shows conserved sex-linked regulation in human (Xiong et al. 2023).
The largest mouse PT sex differences are in S2/S3 (Chen et al. 2023).

These observations come from different studies, assays and segment definitions. Computational
models of the human nephron still take most transporter densities from rat (Layton et al. 2019). To
our knowledge, no study has compared human and mouse PT gene programs position by position along the
S1→S3 axis (literature check, Methods).

**Comparing cell types or segments cannot show where along the PT the species differ.**
Cross-species single-cell atlases align conserved cell states between human and rodent kidney
(Klötzer et al. 2025). A unified mouse–human kidney atlas found little overlap of single-cell
expression changes but agreement at the pathway level (Zhou et al. 2023). These comparisons ask
which cells correspond, not where along a segment a program changes. Dissociation also discards the
spatial relationships between cells (Nitzan et al. 2019).

Assigning structures to S1, S2 and S3 turns the axis into three bins. From those labels alone one
cannot tell whether expression changes in steps or gradually, or whether a boundary sits at the same
place in both species. Continuous reconstructions in other organs show what binning can hide. About
half of mouse liver genes are zonated along the lobule, and many peak in the middle layers (Halpern
et al. 2017). Enterocyte functions are broadly zonated along the intestinal villus (Moor et al. 2018).

**Spatial transcriptomics and pseudospace reconstruction.** Spatial transcriptomics keeps tissue
context. In kidney it has mapped cell types and injury niches, together with single-nucleus data
(Melo Ferreira et al. 2021), with Slide-seqV2 (Marshall et al. 2022), and in the human kidney atlas
(Lake et al. 2023). Visium HD measures the whole transcriptome at single-cell-scale resolution
(Oliveira et al. 2025), and its 2-µm bins can be aggregated into cells or larger structures (Polański
et al. 2024).

A section cuts each tubule at an arbitrary point, so a structure's position along the nephron has to
be reconstructed. Several approaches exist:
- landmark genes have been used to infer lobule and villus coordinates (Halpern et al. 2017; Moor et
  al. 2018);
- spatial arrangements have been recovered from expression alone (Nitzan et al. 2019);
- spatial data have been unrolled to expose zonated programs (Guo et al. 2026);
- trajectory methods order cells along continuous processes (Saelens et al. 2019), including
  diffusion pseudotime (Haghverdi et al. 2016) and principal curves (Faure et al. 2023).

Condition-by-pseudotime tests exist for single-cell data (Van den Berge et al. 2020; Roux de Bézieux
et al. 2024), but they treat cells, not specimens, as the units.

**Few specimens make species comparisons hard to separate from specimen differences.** Cells or
structures from one specimen are not independent replicates. Methods that ignore variation between
biological replicates are biased, and can report hundreds of differentially expressed genes where no
biological difference exists (Squair et al. 2021; Zimmerman et al. 2021). Aggregating to pseudobulk
has been described both as conservative and underpowered (Zimmerman et al. 2021) and as the better
balanced choice (Murphy & Skene 2022). Both approaches have been benchmarked for multi-sample designs
(Crowell et al. 2020).

For pseudotime, Lamian randomly partitioned samples into groups with no expected difference. Methods
that ignored sample variation then reported thousands of false-positive genes (Hou et al. 2023).
Competitive gene-set tests that treat genes as independent are anti-conservative when genes in a set
are correlated (Goeman & Bühlmann 2007; Wu & Smyth 2012).

Human tissue is scarce, and spatial studies often have only a few sections per species. An analysis
therefore has to show that its calls are specific to the species contrast, not merely that they are
numerous.

**This study.** We profiled two control mouse kidney sections and two cortex sections from one human
donor with Visium HD. We segmented individual tubule cross-sections from the paired H&E images,
aggregated 2-µm bins into tubule-level profiles in a shared one-to-one ortholog space, and
reconstructed a continuous PT coordinate across species. We then tested what that coordinate does and
does not resolve.

To learn when pathway screens are specific with two specimens per species, we used a
specimen-relabeling control in our cohort and many-draw nulls from public atlases. To compare
zonation, we called it separately in each species against same-species ceilings from independent
datasets. Every zonated call was confirmed in an independent atlas of the same species. A call of
"zonated in one species only" also required the other species to be flat in its own atlas. We checked
the human segment labels anatomically and by reference transfer, and tested lead genes against
microdissected mouse tubules of both sexes, mouse single-nucleus data, a cortex-only human atlas and
rat segment proteomics.

Beyond a conserved core of strong markers, mostly apical transporters, the two species zonate largely
different genes. Confidently species-only genes are few, and they are not skewed toward mouse. The
pathway screen mainly reports mouse zonation that the human sections do not share.

---

## 4. Results

### R1 · A continuous cross-species PT coordinate

**Cohort and labels.** We segmented tubule cross-sections in two control mouse kidney sections and
two cortex sections from one human donor. After structure quality control, 26,839 structures
remained (Methods).

Cross-species integration and reviewed clustering identified 12,866 PT structures: mouse 2,852 and
2,808; human 3,438 and 3,768. Each carries an S1, S2 or S3 label. These labels are memberships of
joint Leiden clusters on the cross-species Harmony embedding, named by marker review (Methods 6.7).
They are not anatomical annotations. Their validation is described below, and the zonation analyses
(R3, R5) were repeated with labels transferred from reference atlases.

**Coordinate.** We ordered the PT structures along a nonbranching principal curve fitted to the
integrated embedding (scFates; Methods; Figure 1a, b).
- In every specimen, the median coordinate of S3 exceeded those of S1 and S2.
- S2→S3 transitions fell at 0.63–0.64 in all four specimens.
- S1→S2 transitions fell at 0.26–0.28 in mouse and 0.35–0.38 in human, so human S1 occupies more of
  the coordinate.

**Recovered zonation.** Each gene was scored only on a coordinate rebuilt without it (three gene
folds). Its within-specimen rank correlation with the coordinate was then correlated, across genes,
with external S3-versus-early contrasts (Figure 1c):

| Species | External reference | Coordinate | Segment labels |
|---|---|---|---|
| Mouse | Male snRNA | 0.74 | 0.73 |
| Mouse | Microdissection | 0.53 | 0.51 |
| Human | Cortex snRNA (resolves only convoluted PT versus S3) | 0.25 | 0.27 |

Without fold protection, agreement was 0.75 for mouse snRNA and 0.31 for human. Part of the human
agreement therefore came from the genes that built the coordinate.

Rat microdissection data were not used. Their replicate reliability is 0.03 for S2 − S1 and 0.45 for
S3 − early, which fails the pre-specified rule (≥ 0.5; Methods 6.26). Their concordance with our
mouse zonation is correspondingly low (0.10–0.13).

**No absorption of species differences.**
- A human-only gradient injected into a tenth of the construction genes survived refitting
  (absorption ratio 1.01; Figure 1d).
- For 244 genes whose zonation differs between species in external data, the coordinate kept the
  external direction for 91%.
- We built the coordinate on one half of the reads and tested on the other half. 59 of 66 robust and
  39 of 42 core pathways from an earlier list were still called.
- Reusing the same reads inflated gene statistics by about 30%. The circular arm called more pathways
  than the independent arm (145 against 130), but kept a similar share of the earlier robust list
  (57 against 59 of 66).

**What the coordinate does not resolve** (Figure S1).
- *Within-segment order.* Slopes along the coordinate within reviewed segments agreed only weakly
  with adjacent external contrasts: 0.17 in mouse and 0.10 in human, against about 0 for ordering by
  read depth.
- *Gene-fold stability.* Coordinates built from disjoint thirds of the genes ranked structures within
  segments only moderately alike (0.34–0.70).
- *Prediction.* The coordinate predicted held-out genes in the other specimen no better than segment
  means.
- *Depth.* After thinning every structure to equal depth, human S1 showed a depth dependence.
  Within-segment agreement was 0.36–0.46, and the correlation of position with library size fell
  from −0.31 and −0.38 to −0.04 and −0.13.
- *Registration.* Removing a third of the genes or half the reads moved human segments by up to 0.29
  relative to mouse. In one gene-fold refit the human S1 < S2 < S3 order broke.
- *Numerical sensitivity.* Input changes of order 10⁻¹⁶ moved individual structures by up to 0.44
  (Methods).

**Choice of coordinate (disclosed).**
- *What the rule selected.* The pre-registered selection rule picked the diffusion-pseudotime
  coordinate DPT13, not scFates (Methods 6.10). DPT13 met all three criteria: equal-depth agreement
  0.80 against 0.72 for scFates; external agreement within 0.05 of scFates in each species; and a
  species-registration gap of 0.14 against 0.29. It was also better than scFates on gene-fold
  agreement (0.59 against 0.49) and human external agreement (0.29 against 0.25). The equal-depth
  refit of scFates (SCF13-ED) failed the rule.
- *The override.* The rule did not contain the pre-registered downstream specificity gate: NCDR
  < 0.25, with fewer than 5% of pathways called under relabeling (NCDR is defined in R2). After DPT13
  was seen to fail that gate, we applied it to every candidate. The failure was NCDR 0.40 under the
  matched-only test of notebook 36. No candidate then qualified, and scFates was kept
  as the rule's fallback.
- *What the later rerun showed.* On the paper's primary pathway test (notebook 37's joint test),
  DPT13 is also specific: 45 species calls against 3 and 0 relabeled (NCDR 0.03). It gives 26 robust
  and 15 core pathways. Evaluated on the joint test, the rule plus gate would therefore have selected
  DPT13.
- *Decision.* The authors kept scFates, the coordinate on which the analysis was designed.
- *How the choice matters for pathways* (Table S2, coordinate columns):
  - Of the 36 robust pathways under scFates, 22 are also robust under DPT13 and 32 under SCF13-ED.
  - The direction at the divergence peak agrees for all 40 pathways in the scFates–DPT13 union.
  - Peak positions on the 22 shared pathways correlate at Spearman 0.71.
  - The 22 shared pathways split 13 / 3 / 0 / 1 / 5 across the R3 categories (summary of mouse
    zonation / shared / human-led / mixed / unresolved).
  - Within-species controls differ by coordinate. DPT13 gives 37 positional pathway calls between the
    two mice and 5 between the two human sections; scFates gives 16 and 0.
  - A conserved-anchor coordinate built from 160 externally conserved zonation genes failed to
    register the two species and was not used.

**Segment labels: anatomy and reference transfer** (notebook 42; Figure S10). The labels were tested
without the joint embedding in three ways.

1. *Distance to the nearest glomerulus* followed the label order in both species, but separated human
   S1 from S2 only weakly. Median distances (S1 / S2 / S3) and the Spearman correlation of label
   order with distance:

   | Species | S1 (µm) | S2 (µm) | S3 (µm) | Spearman |
   |---|---:|---:|---:|---:|
   | Mouse | 56 | 93 | 297 | 0.65 |
   | Human | 131 | 148 | 228 | 0.23 |

   - The human values are corrected for pixel scale. The human GeoJSON µm fields assume 0.50 µm per
     pixel, whereas the segmentation ran at about 0.44. The uncorrected values are 149 / 168 / 259 µm
     [VERIFY pixel scale against scanner metadata].
   - Mouse distances are descriptive only. The mouse morphological glomerulus polygons are not
     glomerulus-sized, but expression-defined glomeruli give nearly the same distances.

2. *Glomerulus-attached PT.* At the pre-specified 2 µm, too few human polygons touched a glomerulus
   to judge (9 polygons, odds ratio for S1 1.0; mouse 46 polygons, 3.0). Within 10 µm, S1 was enriched
   less in human (odds ratio 2.0, 237 polygons) than in mouse (4.5, 543 polygons).

3. *Reference transfer.* Labels transferred by expression alone from reference S1/S2/S3 profiles
   agreed with the reviewed labels:

   | Species | Reference | Agreement | κ |
   |---|---|---:|---:|
   | Human | Lake cortex atlas | 0.71 | 0.54 |
   | Mouse | Microdissection | 0.87 | 0.80 |
   | Mouse | snRNA | 0.92 | 0.87 |

   By the pre-specified reading, this is weak support for the human labels.
   - Human S1 and S3 transferred mostly to themselves (81% and 99%).
   - 65% of human "S2" transferred to Lake S3. Lake annotates few S2 nuclei, so part of this may
     reflect the reference.
   - Even with transferred labels, human S1 stays larger than mouse S1 (0.44–0.49 against 0.25–0.34).

The genes used to build the integration are concentrated in the conserved class. They make up 16% of
the notebook 40 conserved genes, against 2% of mouse-only genes. R3 is therefore repeated with
labels transferred on held-out genes.

**Spillover and capture.** Non-PT marker spillover into PT structures was small and similar between
species. Median shares of counts were 0.07% / 0.10% / 0.45% in human S1 / S2 / S3, and 0.14% / 0.13% /
0.29% in mouse. Late PT picked up some thick-ascending-limb signal in both. Human polygons without
bins lie outside the Visium HD capture area: only 0.2% and 1.2% of polygons inside the captured
region lack bins (Methods 6.4).

**What the coordinate adds.** Not resolution finer than segments.
- It replicates measurements across positions within each specimen.
- It localises differences, though peak positions should be read at about segment resolution.
- Over the segment-step screen, it adds power, not specificity: 60 against 16 joint-test pathway
  calls, both with 0–2 calls under relabeling (R2).
- The early/late split of the 36 robust pathways is not an artefact of unequal positional noise.
  - Ten human-high pathways peak at a median of 0.16 and 26 mouse-high pathways at 0.63.
  - Adding the estimated excess human noise to mouse positions (σ = 0.089) moved the split to
    11.9 / 24.1. That accounts for about 15% of the gap between the median peaks and 7% of the
    mouse-high count (σ = 0.10: 14% and 8%).

The zonation analyses (R3) and the gene stories (R5) use segment labels and do not depend on the
coordinate.

### R2 · Average-level pathway screens are not specific with two specimens per species; positional screens are specific only conditionally

**The question.** With two specimens per species, a pathway screen can only be judged against
specimen variation. We measured how often each strategy reports pathways when no group difference
exists in two settings:
- in our cohort, using its two balanced relabelings;
- in 750 same-species and 100 cross-species 2 + 2 designs drawn from public snRNA atlases with
  segment labels (24 mice; 13 human donors).

The atlases have no coordinate. The external test therefore compares average-level screens with a
segment-step positional screen, not with the smooth coordinate.

**The design as motivation, not proof.** With 2 + 2 specimens, the species split and the two
relabelings are the three orthogonal ±½ contrasts of the four specimens (Figure 2a). If specimen
deviations are independent, each relabeling is an exchangeable draw of the species contrast's
section-level null. The design is a split plot (Altman & Krzywinski 2015):
- a species offset is judged against specimens within species, leaving 2 degrees of freedom, so
  per-pathway offset inference is impossible;
- positional terms are replicated within specimens.

This argument motivates position-dependent screens, but it cannot predict their specificity. That
depends on two empirical quantities: the size of the specimen-by-position deviations and how
coherently they move within pathways. A term shared by both human sections (donor, probe panel,
sampled region) cancels in the relabelings. The relabeling control cannot detect such terms; R3 and
R4 address them separately.

**Our cohort.** We reran each pathway strategy unchanged on both relabelings (Figure 2b; Table S1).
The negative-control discovery ratio (NCDR) is the mean relabeled count divided by the species
count.

| Strategy | Species calls | Relabeled calls | NCDR |
|---|---:|---:|---:|
| Whole-PT pseudobulk with signed GSEA | 158 | 193, 255 | 1.42 |
| S1/S2/S3 pseudobulk with signed GSEA | 204 | 343, 255 | 1.47 |
| Gene-level offset (T_level), matched test | 27 | 127, 81 | 3.85 |
| Gene-level offset (T_level), joint test | 0 | 46, 2 | – |
| Segment steps (S1/S2/S3), joint test | 16 | 0, 0 | 0 |
| Smooth position (T_spatial), matched test | 116 | 17, 13 | 0.13 |
| Smooth position (T_spatial), joint test | 60 | 2, 0 | 0.02 |

*Target–decoy estimates* (Elias & Gygi 2007) agreed (Figure 2c). No competitive average-level screen
reached a decoy-estimated false-discovery proportion (FDP) of 10% at any threshold; the lowest
attainable FDPs were 1.00, 0.75 and 0.95. The self-contained pathway-score offset model reached it
with 1,087 calls, because nearly every pathway score differs more between species than between
sections.

*The pathway layer.* Whole-PT DESeq2 called 9,010 genes for species and only 2 and 8 for the
relabelings. Specimen-level variation moves whole pathways together, which gene-permutation tests
assume away (Goeman & Bühlmann 2007; Wu & Smyth 2012). Decoy z² rose with pathway size more steeply
for the offset statistic (inflation a = 1.40, implied correlation ρ = 0.019) than for the position
statistic (a = 1.03, ρ = 0.009) (Figure 2e).

*The probe panels.* These make average-level species comparisons hard to interpret in any case.
Per-gene probe efficiency differs between the human and mouse panels. It cancels in within-species
gradients but not in species offsets, so the 9,010 DESeq2 genes are specific only in the relabeling
sense.

*Within-species controls* (Figure 2d). The two mice gave 48 offset and 16 positional pathway calls;
the two human sections gave 94 and 0. The human sections come from one donor, so their 0 positional
calls say nothing about donor-to-donor positional variation.

*The joint test.* Under relabeling the matched test alone was anti-conservative: on average 59.5
pathways reached p ≤ 0.01 where about 15 are expected. The joint test divides the matched z by a
variance-inflation factor estimated from model residuals. Its residual correlation (median
ρ ≈ 0.010) matched the inflation estimated from decoys (ρ ≈ 0.009). On the joint test all three
coordinates were specific:
- scFates: 60 against 2 and 0;
- DPT13: 45 against 3 and 0;
- SCF13-ED: 66 against 1 and 0.

The coordinate-free segment-step screen was specific as well (16 against 0 and 0), so the coordinate
adds power rather than specificity.

**Many-draw external nulls** (notebook 45; Figure 2f–h). Average-level competitive screens reported
pathways in almost every null split. Whole-PT pseudobulk GSEA gave a median of 13–40 calls per split
across the four same-species pools. In cross-species draws they reported as many under relabeling as
for species (median NCDR 1.34).

The segment-step screen with the joint test behaved differently:
- *Adult human cortex:* no calls (median 0 in 150 splits).
- *Mice:* few calls (median 5 in females, 11 in males).
- *Cross-species draws:* a median of 15 species calls against 0.5 relabeled (NCDR 0.05). 80% of
  draws met the specificity rule. That is exactly the pre-specified 80%, with no margin.
- *Our cohort:* the step-screen relabelings (0 and 0) sit at the median of the external relabeling
  distribution. Our average-level relabelings sit at its 70th–100th percentile.
- *Calibration under the null.* The step screen put 1.0% and 0.6% of pathways at p ≤ 0.01 in the two
  human pools (1% expected) and 3.8–4.8% in the mouse pools. Average-level and offset screens put
  2.4–7.5% there everywhere.

NCDR measures vulnerability to specimen noise, not the false-discovery proportion. In cross-species
draws, 71% of whole-PT GSEA species calls recurred in the all-donor reference (7 human against 12
mouse donors), against a base rate of 3%.

**Mechanism.** Specimen offsets were far more coherent within pathways than specimen-by-segment
deviations. The median correlation of member genes was 0.26–0.50 for offsets against 0.03–0.20 for
positional deviations. Unit-level variance inflation cannot absorb offsets this coherent. In our
cohort, within-species positional deviations were small and incoherent (0.013).

**Where positional screens fail.** Positional screens lose specificity when specimens within a group
differ in positional biology.
- *Age* (post hoc). Male mouse splits unbalanced for pre-pubertal (3-week-old) mice gave median 22
  and 62 null calls with one or two such mice per pair, against 1 when balanced. The step screen is
  detecting a real pre- versus post-pubertal positional difference, which a 2 + 2 design cannot
  separate from the group factor. Sex-biased PT programs emerge after 3 weeks of age (Chen S et al.
  2025).
- *Injury severity* (positive control; Figure 2h).
  - *Detection passed.* In our AKI design (two control and two ischaemia–reperfusion mice), 13 of 14
    injury markers rose most in S3. The positional screens called the injury Hallmarks (TNF-α/NF-κB,
    hypoxia, EMT), and 77–80% of their called pathways changed most in S3.
  - *Specificity failed.* Relabeled calls were 35% (segment steps: 26 against 12 and 6) and 17%
    (T_spatial: 50 against 10 and 7) of the condition calls, against the 5% condition. The two AKI
    mice differed in injury severity, and that difference is positional (specimen × segment coherence
    0.42).

**Conclusion.** Position-dependent screens were specific in our cohort and in the external nulls when
the specimens within a group shared their positional biology: adult human cortex, age-matched mice and
our two control mice. They are not specific by design. In our cohort, specificity rests on small,
incoherent within-species positional deviations (coherence 0.013; 16 and 0 within-species positional
calls), which we report beside every call.

**Disclosed changes to the plan.**
- Decoy FDP was the pre-specified primary metric. It failed for pathway-score models (FDP 0 because
  almost every self-contained score differs between species). NCDR, with the condition that
  relabelings call fewer than 5% of tested pathways, was kept as primary.
- That 5% condition was added after the first run.
- The joint test replaced the matched test after the matched test's relabeling behaviour had been
  seen.
- The many-draw nulls were added in response to the referee. Their decision rule was fixed before
  they ran.

### R3 · Beyond a conserved core of strong markers, PT zonation is largely species-specific

All analyses in this section use segment labels and pseudobulk counts, not the coordinate. The labels
are the reviewed cross-species cluster labels (R1); transferred labels are a sensitivity analysis.
Standing caveat: one human donor (two sections), two male mice, human cortex only.

**Same-species ceilings and a conservation index** (notebook 43; Figure 3b). A low cross-species
correlation means "different genes" only relative to what two datasets of the same species give. We
therefore computed noise-corrected (disattenuated) correlations of gene gradients for every pair of
datasets with replicates:

| Pair | S2 − S1 | S3 − early |
|---|---|---|
| Within human (ours vs Lake; Lake vs Census human) | 0.70 | 0.70–0.93 |
| Within mouse (ours, Census snRNA by sex, microdissection) | 0.63–0.88 | 0.73–0.93 |
| Human vs mouse | 0.09–0.20 | 0.13–0.24 |

Relative to these ceilings, conservation was low. The conservation index is the cross-species
correlation divided by the geometric mean of the human ceiling (ours against Lake) and the mouse
ceiling (ours against Census male snRNA).

| Contrast | Index, our data (95% CI) | Index, external atlases alone |
|---|---|---|
| S2 − S1 | 0.17 (0.12–0.22) | 0.27 |
| S3 − early | 0.28 (0.22–0.33) | 0.28 |
| S3c − early, region-matched (co-primary) | 0.29 (0.23–0.35) | 0.27 |

Every upper bound lies below 0.5, so the pre-registered wording "largely different genes" stands. Rat
microdissection failed the reliability rule and was excluded; all rat-based arguments of earlier
drafts are withdrawn.

**Robustness to the human labels** (notebook 46; Figure S9). Labels transferred from Lake (human), alone or with
microdissection labels for mouse, were cross-fitted: genes in one half defined the labels and only the
other half was evaluated.
- In every setting and half, the conservation index stayed at or below 0.31 for S2 − S1 and
  S3 − early, the two pre-specified contrasts, with a largest upper 95% bound of 0.38. For the
  region-matched S3c − early it reached 0.33 (upper bound 0.43).
- With transferred human labels, the human ceiling for S2 − S1 fell from 0.64–0.67 to 0.53–0.55, and
  human amplitude fell to about 0.50 of mouse.

The claim does not depend on the reviewed human labels.

**Human gradients are smaller as well as different** (notebook 39; Figure S4).
- *Spread.* The noise-corrected spread of S2 − S1 gradients in human was 0.58 of the mouse value in
  our data and 0.47 in independent snRNA. For the region-matched S3c − early contrast it was 0.84.
- *Which genes define the axis matters.* The human/mouse amplitude ratio was about 0.15 when genes
  were placed on a mouse reference axis (0.15 to 0.23 across four mouse references). It was 1.96
  (1.41–3.05) on the healthy human cortex atlas axis. The external atlases alone gave 0.14 on a mouse
  axis and 1.66 on a human axis.
- *Errors-in-variables fits.* Deming and orthogonal slopes on the human axis were 0.16–0.17, not
  about 2. With a weak cross-species correlation, total-least-squares slopes track variance ratios
  rather than covariance.

These ratios project one species' gradients onto the other's axis, so they are not pure amplitude
ratios. The reference-free spread is the amplitude measure.

**Log compression is ruled out; tissue state does not explain the pattern in the data available**
(notebooks 39 and 44; Figures S4, S11).
- *Log compression.* The top expression tertile gave an amplitude ratio of 0.18 and probe-balanced
  genes 0.16.
- *Our donor's injury state.* Our donor's PT injury score, calibrated across platforms in mouse, lies
  at the median of the 14 healthy KPMP reference donors. The score is −0.45 in both sections and
  predicts an altered-state fraction of 0.50. The score separates our AKI mice from our control mice
  (+0.33). Its calibration offset is estimated in mouse and assumed to transfer to human.
- *Injury and amplitude.* Within the Lake donors, amplitude falls with injury: median human-axis
  amplitude is 0.69 in healthy, 0.50 in CKD and 0.35 in AKI donors (ρ = −0.41 post hoc; −0.24 on the
  pre-specified mouse axis). Healthy human donors nevertheless remain far flatter than healthy mice.
- *Mouse AKI.* AKI kept 0.74 of control S2 − S1 amplitude (gene-gradient correlation 0.83). It cut
  S3 − early amplitude to 0.41–0.44 (correlation 0.59).
- *Time of day.* Mouse zonation of fatty-acid oxidation, glutathione synthesis, tyrosine catabolism,
  the lead genes and the structural markers is unchanged between rest (ZT4) and active (ZT16) phases.
  The sign was stable in all four mice for every zonated gene, and the amplitude ratio was 0.99–1.35.
  Mouse sterol-synthesis zonation is not stable across these times: 73% sign-stable, with a change at
  least half the species gap.
- *Fasting.* At whole-kidney scale, 24-h fasting changes the PPARα-target and sterol genes by amounts
  comparable to the species differences. This is an unreplicated point estimate, and no
  segment-resolved fasted-kidney data exist. Human PT expresses the fasting- and stress-responsive
  gene PDK4 strongly: 94% of our human structures against 7% of mouse. The ketogenic gene HMGCS2 is
  low in human. These level differences are suggestive of physiological state, not decisive.

**Classification of individual genes** (notebook 43; Figure 3a, c; Table 1). Zonation was called
separately in each species against that species' own specimen noise. Each zonated call had to be
confirmed in an independent atlas of the same species: 12 male mouse snRNA donors, and 7 healthy
human cortex donors.

The primary, symmetric rules differ from earlier drafts in two ways:
- the human effect floor and flat margin are scaled to the lower human amplitude;
- a gene called zonated in one species only also needs the other species to be flat-compatible in its
  own atlas.

Of 7,407 probe-balanced genes measurable in both species, summed over S2 − S1 and S3c − early:

| Class | Genes |
|---|---:|
| Zonated in the same direction in both species (conserved) | 196 |
| Zonated in opposite directions (reversal) | 75 |
| Zonated only in mouse | 41 |
| Zonated only in human | 76 |
| Flat in at least one species and zonated in neither | 712 |
| Not classifiable with two specimens per species | 6,307 |

By contrast:

| Contrast | Mouse-only | Human-only |
|---|---:|---:|
| S2 − S1 | 31 | 58 |
| S3c − early | 13 | 28 |
| Full S3 − early (mostly outer-stripe mouse S3) | 1 | 110 |

Female confirmation gave 38 against 71 (snRNA) and 31 against 78 (microdissection). Mouse-only genes
with a strong gradient (|mouse| ≥ 1.5 log2, human within ±0.3) number 1 under these rules.

*The earlier mouse-only excess was a threshold artefact and is withdrawn.* Notebook 40 reported 247
mouse-only against 66 human-only genes. Scaling the human thresholds alone gives 59 against 141, and
the symmetric rule gives 41 against 76. With cross-fitted transferred labels, mouse-only genes did not
outnumber human-only genes by 2-fold in 11 of 12 setting × half × contrast checks. The single exception
(17 against 4) reversed in the other half (7 against 22).

*Species specificity is a genome-wide property.* It shows as low gradient correlations and lower human
amplitude, not as a large set of species-only genes. Most genes cannot be classified with two
specimens per species.

**What is shared and what is not.**
- *The conserved core is mostly strong transport and segment-identity genes.* Examples include SGLT2
  and SGLT1 (Slc5a2, Slc5a1), Slc5a8, Slc5a10, Slc7a13, Slc7a8, Slc6a18, Slc22a7 (OAT2), Slc22a12,
  Slc34a3, and the intrarenal renin–angiotensin genes Agt, Ace and Enpep. Slc7a13 is the conserved
  reference in Figure 5.
- *Transporters are conserved more often than other genes.* Under the symmetric rules, 33 of 43
  zonated transporters (Slc, Abc and Aqp genes; 77%) are conserved, against 51% of all zonated genes.
  Under the superseded notebook 40 rules the figures were 30 of 54 (56%) against 29%.
- *Agreement with human protein maps is limited to SGLT1/2.* Two of the three OATs disagree with the
  protein data (Breljak et al. 2016):
  - OAT1 protein is stronger in S1/S2 than in medullary-ray S3, but our human OAT1 (Slc22a6) mRNA
    shows no clear gradient (indeterminate);
  - OAT2 protein is likewise S1/S2 > S3, but OAT2 mRNA rises toward S3 in our data and in the KPMP
    atlas.

  Possible explanations are mRNA–protein differences, antibody specificity (the antibodies were
  validated only on transfected cells), and the interindividual variability the authors report. AQP1
  is indeterminate in our classes.
- *Conserved zonation is not conserved usage.* In adult human kidney SLC34A3 carries about 40% of
  sodium–phosphate cotransport and is negligible in mouse (Fernandes et al. 2026), although Slc34a3
  zonation is conserved here.
- *Species-only and reversal genes are few.* Examples:
  - mouse-only: Hsd17b4, Ebp, Fah, Comt, Maoa and Igfbp4;
  - human-only: Sel1l3, Slc9a3 and Ace2;
  - reversals: Dcxr, Ugt3a1, Acaa2 and Vnn1.

  No class enrichment has been computed under the symmetric rules; the notebook 40 enrichments
  (sterol synthesis, peroxisomal β-oxidation and tyrosine catabolism among mouse-only genes) belong
  to the superseded classes.

**Position-dependent pathways mostly summarise mouse zonation** (notebooks 37 and 43; Figure 3d;
Tables S2, S4).

*The reporting set.* Of 1,513 tested pathways, the joint test gave:
- 60 candidates;
- 58 not specimen-sensitive;
- 51 stable in at least 5 of 7 planned refits;
- 36 robust to coordinate registration and sampled region (17 gene-overlap programs);
- 26 core, also robust to restricting to genes detected in both species and to removing mouse
  sex-biased genes.

*What the robust pathways contain.* Each robust pathway was classed by its contributing members
(top 10% by T_spatial) zonated in either species:

| Category | Pathways |
|---|---:|
| Summary of mouse zonation (≥ 50% of members mouse-led) | 26 |
| Shared (Nitrogen metabolism, Aspirin ADME, Glutathione metabolism) | 3 |
| Mixed (Drug ADME, Mitochondrial Fatty Acid Beta-Oxidation) | 2 |
| Unresolved (fewer than 3 zonated members) | 5 |
| Human-led | 0 |

The mouse summaries include fatty-acid oxidation, peroxisome, bile acid, branched-chain amino-acid
and vitamin metabolism.

*Why.* A species × position statistic cannot recover a conserved core. A conserved gene enters
T_spatial only through an amplitude difference, and the interaction is about −0.8 to −0.9 times the
mouse gradient in our data and in both external human references (R4). The robust list is therefore
best read as "programs whose mouse zonation the human sections do not share".

*Probe sensitivity.* Eight robust pathways lose significance with equal-probe genes only:
- Glutathione metabolism;
- Steroid hormone biosynthesis;
- Chemical carcinogenesis;
- Retinol metabolism;
- Metabolism of xenobiotics by CYP450;
- Aspirin ADME;
- Formation of Cornified Envelope;
- Biological Oxidations.

Of the four pathways that notebook 39 placed "beyond global flattening", Glutathione metabolism and
Formation of Cornified Envelope fail this check and are dropped. Estrogen Response Early and Late
remain, but their members are mostly mouse-led.

*Peaks.* The 10 human-high pathways peak early (median 0.16) and the 26 mouse-high pathways late
(0.63). Most differences are graded and one-signed: 26 graded, 4 localised and 6 reversing, with a
median averaging loss of 0. The early/late split therefore mainly reflects mouse late-PT programs that
the human sections do not zonate.

*Earlier list.* The earlier 82 → 66 → 42 list is kept as a sensitivity analysis. It contains 27 of
the 36 robust and 17 of the 26 core pathways.

### R4 · Robustness and agreement with independent data

**Robustness of the pathway calls** (Figure 4a; Table S2).

*Six checks keep the effect.* Every robust pathway keeps at least 70% of its effect under each of
these:
- an anatomy-anchored coordinate, in which human sections are warped so that their segment
  transitions fall on the mouse transitions;
- cortical-like mouse S3 only;
- removal of genes flagged as ambient, neighbouring-segment or shared-exon;
- 5-bin covariate matching;
- equal probe counts.

*Significance.* The first two checks define "robust", so they keep significance (joint q ≤ 0.10).
The flagged-gene and 5-bin checks keep it too. With equal-probe genes only, 8 of the 36 lose
significance although they retain their effect (R3). The two gene-subset checks that define the core
set were passed (effect ratio ≥ 0.7 and q ≤ 0.10) by:
- 28 of 36 pathways for genes detected in both species;
- 32 of 36 after removing mouse sex-biased genes.

*Coordinate.* Pathway calls depend on the coordinate (R1): 22 of the 36 are also robust under DPT13
and 32 under SCF13-ED. The direction at the divergence peak never disagrees between coordinates. The
earlier notebook-12 analysis on notebook 03's DPT coordinate is historical and is not used as support.

**Specimen-level split-plot test** (Figure 4e). This test treats specimen × position pseudobulks as
replicates (Methods 6.18).
- *Genes.* It was calibrated under relabeling: 0.6% and 0.4% of genes at p ≤ 0.01, λ ≤ 0.87. Its gene
  ranking agreed only weakly with T_spatial (Spearman 0.21).
- *Pathways, pre-specified joint test.* It called nothing; decoys showed that its pseudobulk
  correlation correction overcorrects (decoy a = 1.10 against a median correction of 2.07).
- *Pathways, matched test alone.* It called 41 pathways, with 0 and 0 under relabeling. 15 of the
  36 robust pathways were among them, and the robust pathways' split-plot z exceeded that of uncalled
  pathways (median 2.6 against −0.5).

**External agreement is at the measurement level** (Figure 4c, d; Tables S3, S5). Donors are the
replicates. The human references were cortex snRNA (Census, 6 donors, primary; Lake/KPMP, 7 donors),
each compared with 12 male or 12 female mice.

*Genes.* Our top 10% of genes by T_spatial correlated with the external species × position contrast
above expression-matched chance, but about as well as the most highly expressed genes:

| Reference | Our top-10% genes | Matched random genes, median (97.5th percentile) | Highest-expressed genes |
|---|---:|---|---:|
| Census | 0.67 | 0.57 (0.62) | 0.64 |
| Lake/KPMP | 0.73 | 0.60 (0.66) | 0.72 |

*Pathways.*
- Robust pathways replicated better than uncalled ones. They replicated only marginally better than
  pathways called under relabeling, which are specimen artefacts by construction (Mann–Whitney
  p = 0.016 and 0.033 with male-mouse references; 0.19 and 0.21 with female).
- Only 17 of 36 replicated individually against Census, and 6 against Lake (BH ≤ 0.10 within the
  list), so the pre-specified criterion of at least 50% was not met.
- All three human datasets share a global component: the interaction is −0.81 to −0.90 times the
  mouse gradient (R² 0.60–0.82; Table S5). Once it is removed, 1–3 robust pathways replicate in any
  reference.
- Pathways called only by conventional screens (184) replicated their average difference 6 times at
  nominal p ≤ 0.05 and never at BH ≤ 0.10.
- These pathway comparisons treat overlapping, selected pathways as independent, so their p-values
  are nominal. We therefore claim agreement of gene measurements, not pathway-by-pathway replication.

*Genome-wide, against the cortex-only human atlas* (segment labels; Figure S7).
- Our human S3 − S1 contrasts correlated with Lake at Spearman 0.43 (6,787 genes), with 93.1%
  direction agreement for 335 strongly zonated genes.
- For the 10% most position-dependent genes, our interaction correlated with Lake minus male mouse
  snRNA at 0.70, with 96.7% direction agreement for 276 genes strong in both. With female mice as the
  reference the figures were 0.67 and 93.7%.

*Additional human donors with spatial data* (notebook 44). These were not informative. The authors'
labels in the Abedini et al. (2024) Visium samples failed the pre-specified S3 marker check: SLC5A1,
SLC7A13 and SLC22A7 did not rise in their PT_S3 spots. Descriptively, the S1-high human arms of UGT1A9,
UGT3A1, ACOX2 and SLC6A19, and the rise of RBP4, replicate in those samples.

### R5 · Gene programs placed differently along the PT

**How each claim is stated.** Every gene claim below carries its notebook 43 symmetric class, read
over S2 − S1 and S3c − early on reviewed labels, and says no more than that class allows:
- *conserved:* zonated in the same direction in both species;
- *reversal:* zonated in opposite directions;
- *mouse-only or human-only:* zonated in one species and flat in the other;
- *indeterminate:* at least one species can neither be called zonated nor confidently flat with two
  specimens.

"Flat in human" is used only where the human call is flat. A separate statement, marked *T_spatial*,
gives the gene's rank on the species-by-position interaction statistic. T_spatial is a ranking among
11,071 genes, not a gene-level test, and a high rank can arise from an amplitude difference alone.

**How each claim was checked.** External checks use reviewed segment labels against three references:
- the cortex-only human atlas (Lake et al. 2023; 7 donors, 4 female and 3 male);
- mouse snRNA by sex (12 male and 12 female donors);
- microdissected mouse segments from both sexes (Chen et al. 2021, 2023).

Pre-specified rules gave 30 headline, 5 supporting and 4 non-replicating genes of 39 claimed (Figure
S5). Rat segment proteomics (Limbutara et al. 2020) adds protein for the rodent side, quoted as
S3/S2 ratios because rat S1 samples carry less protein. Human protein evidence is presence only.

Standing caveat: one human donor, two male mice, human cortex only, different probe panels.

**Creatine synthesis: a large difference in amplitude** (Gatm and Gamt indeterminate; Figure 5).

- *Mouse arm.* Gatm (AGAT, the first step) falls steeply from S1 to S3 in mouse. The fall is 6.8 and
  6.3 log2 units in male and female snRNA and 12.1 and 11.0 in male and female microdissection. Even
  the cortical medullary-ray segment of mouse has 23.5 TPM, against 4,270 in S1.
- *Rodent protein.* In rat, AGAT protein falls at S3 (S3/S2 −3.7). Transamidinase activity is
  confined to S1 and S2 (Takeda et al. 1992).
- *Human arm.* Human GATM changes little and not consistently: −0.78 in our sections and +0.51 in
  Lake. It is neither confidently zonated nor confidently flat, so the class is indeterminate.
- *Defensible claim.* Mouse confines the first step of creatine synthesis to early PT, whereas human
  GATM extends along the PT. This is a large amplitude difference, not a human reversal.
- *Ranking.* Gatm ranks in the top 0.1% of genes by T_spatial.
- *Second step.* The second step, GAMT, rises toward late PT in rodents: mouse microdissection in both
  sexes, and rat protein (S3/S2 +1.6). Rat segment proteomics thus supports, at protein level, the
  separation of the two steps along the rodent tubule. Our Visium data do not resolve Gamt: it is
  indeterminate in both species (mouse flat in S3 − early), with unequal probes (3/4), and human GAMT
  is flat in the atlas. Figure 5 shows Gamt only as the creatine partner of Gatm.
- *Literature.* GATM is a human PT enzyme (Reichold et al. 2018), and the human kidney releases
  guanidinoacetate in vivo (Edison et al. 2007). The human axial profile and the species contrast
  were not found in the literature.

**Mitochondrial β-oxidation** (Acadm conserved; Acaa2 reversal; Figure 5).

- *Acadm (MCAD)* rises toward late PT in both species. The rise is much larger in mouse: snRNA +2.1 and
  +2.3, against +0.46 in Lake. The class is conserved, so the earlier claim "flat in human" is
  withdrawn.
- *Acaa2* (mitochondrial 3-ketoacyl-CoA thiolase) is confined to early PT in mouse (snRNA −4.0 and
  −4.2). It rises modestly in human: Lake +0.58, +0.95 in male donors. It is a reversal in our S2 − S1
  data.
- *Rodent protein.* Rat protein matches both mouse directions: Acadm S3/S2 +0.7 and Acaa2 −2.7.
- *State checks.* Neither gene is sex-biased in mouse PT, and both proteins are high in human tubules.
  Their mouse zonation is unchanged between ZT4 and ZT16, so it is not a time-of-day artefact.
  Whether fasting flattens it in human cannot be tested with public data.
- *Ranking.* Both rank in the top 1% by T_spatial.
- *Other fatty-acid oxidation genes.* Acsm3, Crot and Nudt19 rise in mouse but are indeterminate in
  our classes. Acsm3's late rise is moreover a male-mouse program. Of the peroxisomal genes, only
  Hsd17b4 is confidently mouse-only.
- *Literature.* No segment-resolved data were found for these genes. In rat, 3-hydroxyacyl-CoA
  dehydrogenase activity is roughly uniform along cortical segments (Le Hir et al. 1982; Bastin et al.
  1990), so activity assays would not see gene-specific zonation. An intravital study found that late
  mouse PT mobilises lipid droplets with lipase-rich lysosomes (Kaminska et al. 2026). Early- and
  late-PT mitochondrial fatty-acid oxidation capacities differ and change with 24-h fasting (Feola et
  al. 2026). PT PPARα is fasting-activated and drives renal fatty-acid oxidation (Aomura et al. 2026).

**Glutathione synthesis and drug handling** (Gclc and Gclm indeterminate; Gss conserved; Cyp2e1 not
classifiable; Figure 5 shows Gclc and Cyp2e1). This story was promoted to lead against the pre-specified rule (disclosed); under the
symmetric classes it is a supporting observation.

- *Glutathione synthesis.*
  - Gclc, Gclm and Gss rise toward late PT in mouse, in both sexes and in microdissection.
  - Rat protein agrees: S3/S2 +2.8, +1.9 and +2.4. Late placement is therefore supported on the
    rodent side.
  - In human the rise is small or absent (GCLC +0.61 in Lake). Both GCL subunits are indeterminate
    (Gclc also probe-imbalanced, 3 mouse against 6 human probes), and Gss is conserved.
  - Glutathione metabolism is a robust pathway, but it is probe-sensitive and shared rather than
    mouse-led (R3).
  - Axial glutathione biology differs between species. In rabbit PT, synthesis is highest in S1 and
    cellular GSH highest in S3 (Parks et al. 1998). In rat, glutathione is highest in convoluted and
    early straight PT (Brehe et al. 1976).
- *Cyp2e1.*
  - Mouse Cyp2e1 is confined to S1–S2 in our data and in microdissection. Rat protein peaks in S2
    (S3/S2 −5.5).
  - Mouse renal CYP2E1 is a known proximal-tubular, male and testosterone-regulated enzyme (Hu et al.
    1990; Speerschneider et al. 1995).
  - Human CYP2E1 mRNA is essentially absent, so it cannot be classified. The literature on human
    protein conflicts (Speerschneider et al. 1995; Akakpo et al. 2023), so we state low mRNA, not
    absent protein.
- *UGT1A9* is expressed only in human and highest in early PT (Lake −0.67). Its S1-high profile also
  appears in the additional Visium donors (R4). UGT1A9 is among the most abundant renal UGTs in human
  (Margaillan et al. 2015); its axial profile is not reported. It is not shown in Figure 5. UGT1A isoforms share exons
  2–5, so the probes must target the unique exon 1 for this to be a UGT1A9 call [VERIFY probe
  positions].
- *Reading.* The rodent placements of glutathione synthesis (late) and Cyp2e1 (S1–S2) are supported
  at transcript and protein level. The human side shows weaker or no gradients, not a different
  placement.

**Genes whose difference changes sign** (Dcxr, Ugt3a1 reversal; Figure 5).

- *Dcxr* (dicarbonyl/L-xylulose reductase) falls toward late PT in mouse and rises in human. It is a
  reversal in S3c − early and in S3 − early.
  - *Human rise.* It holds in the two human sources with valid labels: our sections +0.97 and Lake
    +1.53 (+1.59 in men, +1.48 in women).
  - *Mouse fall.* It holds in mouse snRNA (−0.46, −0.54) and microdissection (−1.94, −2.03). Rat
    protein agrees (S3/S2 −1.1).
  - *Other checks.* Probes are balanced. The Visium donors cannot test a rising gene, because their S3
    labels failed the marker check.
  - *Literature.* DCXR protein is known in mouse PT (Nakagawa et al. 2002), and DCXR loss tracks human
    CKD (Perco et al. 2019). Its axial distribution is not reported in either species.
- *Ugt3a1* rises toward late PT in mouse (snRNA +1.0 and +1.5) and falls in human (ours −0.62, Lake
  −0.98). It is a reversal in S3c − early, and its human S1-high arm replicates in the additional Visium
  donors. Rat protein does not quantify it.
- *Ranking.* Both rank in the top 1% by T_spatial.
- *Genes previously called reversals.* Earlier interaction-sign rules flagged more genes. The symmetric
  classes reassign them:

  | Gene | Class |
  |---|---|
  | Pah, Cyp24a1 | Conserved |
  | Igfbp4 | Mouse-only |
  | Acox2, Glyat, Nt5e, Cyp7b1 | Indeterminate |

**Expressed and zonated only in human** (Rbp4, Aox1). These genes do not pass the expression
requirement in mouse, so they are outside the symmetric classes. They are listed with the genes
expressed in one species only (Table S7).

- *RBP4* rises steeply toward late human PT (ours +2.46, Lake +3.46). The rise also appears in 7 of 8
  informative Visium samples.
  - Mouse Rbp4 is absent: 0.0 TPM in every microdissected segment of both sexes.
  - Rat RBP protein is found only in S1, where it is reabsorbed ligand. Rat in situ hybridisation
    placed RBP mRNA in outer-stripe S3 and also in perinephric fat (Makover et al. 1989).
  - The claim is therefore an mRNA claim, and ambient signal from perirenal fat is untested.
    Transcript diffusion between neighbouring Visium HD bins is a further alternative. Protein
    signals, such as the Human Protein Atlas staining, may include filtered RBP4 taken up by the
    tubule, as the rat protein pattern suggests.
- *AOX1* rises toward late human PT (Lake +1.13). It is supporting only (2 human against 3 mouse
  probes). Humans carry one AOX gene where mice express four (Terao et al. 2016).

**Organic anion uptake** (Slc22a6 indeterminate; Slc13a3 conserved).

- *OAT1 (Slc22a6)* peaks in S2 in mouse (970 TPM in the cortical medullary-ray segment) and in rat
  protein (S3/S2 −6.5). Human OAT1 mRNA shows no clear gradient. Human OAT1 protein is reported
  stronger in S1/S2 than in medullary-ray S3 (Breljak et al. 2016), so the protein map and our
  transcript call disagree.
- *NaDC3 (Slc13a3)* falls toward late PT in both species (conserved in S3c − early). In S2 − S1 it
  rises in mouse only (mouse-only in that contrast), so its union class is conserved but its early-PT
  behaviour differs. The v0.2 call "mouse-only" is withdrawn.
- *Whole-PT averages.* The module illustrates how a whole-PT comparison can hide regional differences.
  Slc22a6, Slc13a3 and Cyp24a1 have non-significant whole-PT DESeq2 differences but significant segment
  differences of opposite sign (Figure S6). Those S3 differences compare human cortex S3 with mouse S3
  that includes outer stripe.

**Supporting observations.**
- *Sterol synthesis.*
  - Mouse sterol genes rise toward late PT. Their zonation varies with time of day in mouse, and
    fasting changes them at whole-kidney scale by amounts comparable to the species gap.
  - Classes: Hmgcr is conserved, Ebp mouse-only, and Hmgcs1, Cyp51 and Lss indeterminate.
  - Reported as supporting only; the earlier rat argument is withdrawn.
- *Early-PT transporters in late PT.*
  - Slc9a3 is zonated only in human (rising S1→S2), Slc4a4 is conserved, and Slc6a19 is
    indeterminate.
  - The human S1-high arm of SLC6A19 replicates in the Visium donors.
  - Mouse B0AT1's S1–S2 restriction is known (Navarro Garrido et al. 2022).

**What did not hold.**
- *Serine synthesis is not human-specific.* Psat1 is conserved, though probe-imbalanced. Mouse Phgdh
  and Psat1 each have one probe against three in human, and mouse microdissection shows both
  S1-restricted, as in human.
- *Not replicated:* Gamt (in our data), Hadh, Ephx1 and Me1.

---

## 5. Discussion

**What is new.** To our knowledge this is the first comparison of human and mouse PT gene programs as
a function of position along the tubule, measured on intact tubules in tissue in both species. It
yields three results.

1. **When pathway screens are specific.** With two specimens per species, pathway screens that compare
   averages cannot be separated from specimen variation, in our cohort and in hundreds of external
   2 + 2 designs. Position-dependent screens are specific when the specimens within each group share
   their positional biology. They are not when specimens differ in age or injury severity. The
   split-plot argument motivates this but does not guarantee it; in our cohort the evidence is the
   small, incoherent within-species positional variation.
2. **Zonation is largely species-specific beyond a conserved core.**
   - Against same-species ceilings, human–mouse conservation of gene gradients is 0.17 for S1→S2 and
     about 0.28 for S3 versus early PT. The external atlases alone give 0.27–0.28.
   - With human labels transferred from a reference atlas on held-out genes, the index stays at or
     below 0.31 for both contrasts (0.33 for the region-matched S3c − early).
   - The shared component is concentrated in strong transport and segment-identity genes:
     transporters are conserved more often (77%) than zonated genes overall (51%).
   - Confidently species-only genes are few and not skewed toward mouse (41 against 76). Human
     gradients are smaller (S1→S2 spread 0.47–0.58 of mouse; 0.84 for the region-matched S3
     contrast).
   - The earlier reading that mouse zonates many more programs was a threshold artefact and is
     withdrawn.
3. **A short list of positional differences, at stated levels of evidence.**
   - *Reversals:* Dcxr, Ugt3a1 and Acaa2.
   - *Expressed and zonated only in human:* RBP4, and AOX1 as a supporting observation.
   - *Large amplitude difference:* creatine synthesis, which mouse confines to early PT.
   - *Rodent placements with weaker or no human gradients:* glutathione synthesis (late) and Cyp2e1
     (S1–S2).

   Rodent zonation of several programs was known from microdissection or enzyme activity. Rat protein
   now supports the rodent side of creatine, glutathione, β-oxidation and Dcxr. The human arm and the
   species contrast are new.

**What the spatial data add.** The central claim can be reproduced from public single-nucleus atlases
alone (external conservation index 0.27–0.28). Our data add three things:
- the first intact-tissue measurement of both species on one platform;
- region matching by physical depth, which overturned an earlier reading of the OAT module and lets
  S3 be compared within cortex;
- replication of measurements across structures and positions within each specimen, which the
  positional screen and the specimen-level checks use.

**Relation to prior work.**
- *Human atlases* label PT-S1/S2/S3 as discrete classes of dissociated nuclei (Lake et al. 2023).
- *Mouse microdissection* gives three anatomical points (Chen et al. 2021, 2023).
- *Rat segment proteomics* gives protein along the rodent tubule (Limbutara et al. 2020).
- *Cross-species atlases* align cell states (Klötzer et al. 2025). They find little overlap of
  single-cell expression changes but agreement at pathway level (Zhou et al. 2023).
- *Cross-species integration* is known to trade species mixing against biological conservation (Song
  et al. 2023). This is one reason we do not interpret the integrated embedding biologically and check
  the labels by reference transfer.
- *The mouse sex atlas* compares whole-PT sex differences with human, not positions (Chen S et al.
  2025).

The 2025 cross-species atlas is not open access, and its PT annotation depth could not be checked;
this is the main residual novelty risk.

**Why "flattening" is mostly a mouse-centric view.** Measured on mouse-zonated genes, human amplitude
looks like 0.15 of mouse; measured on human-zonated genes, it looks like 2–3. The reference-free
spread shows a moderate global reduction. It is about half for S1→S2 and smaller for the
region-matched S3 contrast. Within healthy reference donors, injury lowers amplitude further, and our
donor sits at the median of those donors.

**Implications for mouse models of PT physiology and toxicity.**
- *Transport.* For strong apical transport and segment-identity genes, mouse predicts where human PT
  zonates them more often than for other genes. Conserved zonation is not conserved usage, however:
  SLC34A3 carries a large share of phosphate transport in human but not in mouse (Fernandes et al.
  2026). OAT1 and OAT2 transcript profiles also disagree with human protein maps.
- *Metabolism.* Several rodent metabolic placements have weak or no counterpart in the human sections,
  including early-PT creatine synthesis, late-PT glutathione synthesis, β-oxidation steps and sterol
  synthesis. An effect confined to one mouse segment may be spread along the human PT, or weaker.
- *Physiological state.* Mouse PPARα and sterol programs respond strongly to fasting (Rakhshandehroo et
  al. 2009 for species differences in PPARα targets; Aomura et al. 2026). The renal tubular clock
  controls β-oxidation (Bignon et al. 2023). Laboratory mice are fed ad libitum and young, so part of
  any mouse–human difference in these programs may be physiological state.
- *Injury models.* Experimental injury models already differ in where injury falls along the nephron
  (Heyman et al. 2010). Our map adds that many programs are zonated differently.
- *Computational models.* Human nephron models borrow rodent axial parameters (Layton et al. 2019).
  Strong transporters support that choice more than metabolic genes do.
- *Caution.* All of these are mRNA-level observations.

**Limitations.**

*Design and specimens.*
- *Replication.* Two male control mice and two cortex sections from one male human donor. Specimens,
  not structures, are the replicates. Agreement between the two human sections is not donor
  replication, and no statement here is a population estimate for humans.
- *Sex.* All specimens are male, so species and sex are confounded. Female mouse references confirm
  most late mouse programs but not Acsm3, and female confirmation gives similar class counts.
- *Region.* Human sections are cortex only, so mouse outer-stripe S3 has no counterpart. Region checks
  (cortical-like mouse S3; the cortical medullary-ray mouse segment) changed the conclusion for the OAT
  module. With all mouse S3 included, human-only genes outnumber mouse-only ones (110 against 1).

*Labels and classification.*
- *Human segment labels.* The labels are cross-species cluster memberships. They have weak reference
  support (transfer agreement 0.71 against 0.87–0.92 in mouse), and human S2 is the least supported
  label. The conservation index is robust to transferred labels, but the S2 − S1 contrast is the
  noisiest.
- *Classification.* Most genes (6,307 of 7,407) cannot be classified with two specimens per species.
  The human confirmation atlas has fewer donors than the mouse one.

*Measurement.*
- *Probe panels and processing.* The two species were measured with different probe panels and Space
  Ranger versions. Genes with unequal probes are reported separately, and one apparent species
  difference (serine synthesis) was a panel artefact.
- *Coordinate.*
  - Within-segment order is weakly supported.
  - Human S1 ordering depends on read depth.
  - Human–mouse registration moves under perturbation, and the curve is numerically sensitive.
  - Pathway calls depend on the coordinate: 22 of 36 are robust on DPT13. Their directions do not.
- *External references.* External human positions are single-nucleus cluster annotations, not anatomy.
  Pathway-level replication beyond the shared global component is not shown.

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

The source, age and procurement of our donor are not recorded here [VERIFY].

*Decisions taken after seeing results* (disclosed):
1. The pre-registered coordinate rule selected DPT13. The downstream specificity gate was applied after
   DPT13 failed it on the matched-only test, although DPT13 is specific on the joint test. The authors
   kept scFates.
2. The 5% condition was added to the specificity rule after the first run.
3. NCDR replaced the pre-specified decoy FDP as the primary metric.
4. The joint test replaced the matched test after the matched test's relabeling behaviour was seen.
5. The drug-handling story was promoted to lead against the pre-specified rule. It is a supporting
   observation under the symmetric classes.
6. The notebook 39 diagnostics that turned "flattening" into "different genes" were post hoc.
7. In notebook 45, the age analysis of mouse null calls was post hoc.

**What would move this forward.**
- Several human donors with spatial data at single-cell resolution and anatomical segment
  annotation.
- Female animals and female human donors.
- Fasted-mouse or fed-human references.
- In situ validation of position-defining genes on human sections where S1 and S3 can be
  distinguished: GATM, ACAA2, DCXR, UGT3A1 and RBP4.

---

## 6. Methods

This section merges `docs/paper/methods_draft.md` (workstream 4), `docs/paper/methods_pathway_draft.md`
(workstream 2) and the corrections in `docs/paper/methods_corrections_ws1.md` (workstream 1). It
adds the zonation amplitude (notebook 39), zonation class (notebook 40) and biology verification
(notebook 38) methods. Parameters were read from the code.

Version 0.3 adds the revision analyses:
- segment-label validation (notebook 42; 6.25);
- same-species ceilings, the conservation index and symmetric classes (notebook 43; 6.26);
- physiological and tissue state, additional human donors and rodent protein (notebook 44; 6.27);
- many-draw specimen nulls and the AKI positive control (notebook 45; in 6.14);
- cross-fitted label sensitivity (notebook 46; 6.28);
- coordinate robustness of pathway calls (6.29);
- figure assembly (notebook 41; 6.30).

### 6.1 Study design and specimens

**Cohort.**
- Mouse: two control mouse kidney sections, `Ctrl1A2` and `Ctrl1A4`.
- Human: two sections of renal cortex from one donor, `HUK1_COR1` and `HUK1_MED1`. The second is
  named "MED" at source, but its segmentation contains cortex only, and it is analysed as cortex.

All four specimens express Y-linked transcripts and are treated as male [VERIFY against records].
The mouse sections are whole-kidney sections, so they include cortex and outer medulla. The two
ischaemia–reperfusion kidneys of the same mouse experiment (`IR2A2`, `IR2A4`) are used only in the
tissue-state analyses (notebooks 39 and 44) and as the positive control of notebook 45 (6.14).
Specimens are the units of replication.

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
- The human `*_v2.geojson` re-segmentations have no recorded provenance [VERIFY]. Their µm fields
  assume 0.5 µm per pixel, whereas the segmentation model works at 0.44068 µm per pixel (6.25)
  [VERIFY against scanner metadata].
- Mouse analyses use only the quality-controlled `*_kept_tubules_labeled_fine.geojson`. The rules
  and code of that upstream QC are not in the repository [VERIFY].
- Only polygon geometry is used downstream.

**Centroid check.** Every join between expression and segmentation recomputes the indexed polygon's
centroid and stops at any deviation above 1 px.

### 6.4 Tubule-by-gene matrices

Each 2-µm bin was assigned to the polygon covering its centre, in full-resolution image coordinates
(scale factor 1, shapely STRtree, at most one polygon per bin). Raw counts were summed per polygon.
The segmentation and the Visium HD image must share one coordinate frame. How they were registered
is not recorded [VERIFY]. Notebook 42 checked the result (6.25):
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

**What the labels are.**
- The S1/S2/S3 labels are memberships of joint cross-species Leiden clusters on the integrated
  embedding, named by marker review. They are not anatomical calls made separately in each species.
- How joint clustering merges or splits cell types across species depends on the integration
  strategy and its strength (Song et al. 2023).
- Notebook 42 tests the labels against distance to glomeruli and against expression-only reference
  transfer (6.25). Notebook 46 repeats the conservation index and class counts with transferred
  labels (6.28).
- No current notebook records why θ = 6 was chosen [VERIFY rationale].

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

**Outcome (disclosed).**
- The rule selected DPT13 (P5 0.80; P1 0.74 mouse, 0.29 human; gap 0.14).
- The scFates equal-depth refit failed (P5 0.73; gap 0.36).
- scFates itself would also have failed (P5 0.72; gap 0.29).
- The rule omitted the pre-registered downstream specificity gate: criterion 5 of rule C6, which
  requires the T_spatial screen to keep NCDR < 0.25 and fewer than 5% of pathways called under
  relabeling.
- After DPT13 was seen to fail that gate on notebook 36's matched-only test (NCDR 0.40), we applied
  the gate to every candidate. No candidate then qualified, and scFates was kept as primary.
- Notebook 37's full pipeline was then run on DPT13 and on the scFates equal-depth refit (SCF13-ED;
  6.29). Under notebook 37's primary joint test, DPT13 is specific: 45 species calls against 3 and 0
  relabeled (NCDR 0.03). The matched-only test in the same run gives 99 calls against 41 and 36
  (NCDR 0.39). DPT13 gives 26 robust and 15 core pathways.
- Whether DPT13 passes the gate therefore depends on the test. The authors kept scFates, the
  coordinate on which the analysis was designed. This is a decision, not an outcome of the
  pre-registered rule.

**Outcomes.**
- *DPT* on notebook 13's embedding (notebook 36): 96 calls; 56 of 66 robust and 38 of 42 core
  pathways from the earlier list retained; NCDR 0.40 (matched test).
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
  both, peak positions correlate with Spearman ρ = 0.71 (6.29).
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
- Count-based trajectory tests such as tradeSeq and condiments (Van den Berge et al. 2020; Roux de
  Bézieux et al. 2024) test conditions along a trajectory with cells as units. With two specimens
  per species, that would treat structures as replicates. We therefore use the gene statistics only
  as ranks, and place inference at the pathway level under specimen relabeling (6.14) and at the
  specimen level (6.18).
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
  second condition was added after the first run (disclosed).
- Decoy-estimated FDP (Elias & Gygi 2007) is FDP(t) = (mean relabeled count ≥ t) / (species count
  ≥ t), made monotone as a q-value. We report calls at FDP ≤ 10%, with 300-draw pathway-resampling
  intervals, and the minimum attainable FDP.
- Decoy FDP was the pre-specified primary R2 metric. It was replaced by NCDR because pathway-score
  models reach FDP 0 (disclosed).
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
  T_spatial. Detection passes if injury programs are called; specificity is judged by the
  relabeled ÷ AKI call ratio.
- The atlases have no coordinate, so the external nulls compare average-level screens with
  segment-step position-dependent screens. They do not test the smooth coordinate itself.

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
  pathways). This is a historical sensitivity only and is not cited as robustness support.
  Robustness to the coordinate is assessed instead by rerunning notebook 37 on DPT13 and SCF13-ED
  (6.29; Table S2 coordinate columns);
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

**Reference-axis (projection) ratio.**
- Both of our species gradients were regressed on an independent reference gradient. The ratio is
  cov(human, ref) / cov(mouse, ref).
- This ratio measures how much of each species' gradient lies along the reference axis. When the
  cross-species correlation is low (R3), it is a projection ratio and understates the human
  amplitude off that axis. The noise-corrected SD ratio (below) is the amplitude measure used in R3.
  (v0.2 called this the "dilution-free ratio".)
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

### 6.20 Zonation classes (notebook 40; superseded as primary by 6.26)

Notebook 40's classes are reported for the record (Figure S8, Table S7). The primary classes of v0.3
are notebook 43's symmetric classes (6.26), which reuse the calls below with amplitude-matched human
thresholds and external flat confirmation.

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
annotation. This check was added after the classes were computed (disclosed).

### 6.21 External segment-resolved datasets and verification of lead genes

**Datasets.**

| Dataset | Content |
|---|---|
| GSE150338 (Chen et al. 2021) | Microdissected male mouse PTS1–3; TPM replicates |
| GSE212213 (Chen et al. 2023) | Microdissected mouse PTS1–3, both sexes; mean TPM |
| GSE56743 (Lee et al. 2015) | Microdissected rat S1–S3; RPKM; excluded by notebook 43's reliability rule (6.26) |
| CELLxGENE Census 2025-11-08 (CZI Cell Science Program et al. 2025), mouse dataset `25818bf7…` (Chen S et al. 2025) | 12 male and 12 female donors with ≥ 50 nuclei per segment |
| Census human renal-cortex dataset `09b518f9…` (Acera-Mateos et al. 2026) | 6 donors; convoluted PT vs S3 |
| Lake et al. (2023) snRNA v1.5 (CELLxGENE `a12ccb9b…`) | Cortex nuclei (`region C`, cortex tissue), healthy donors, author labels PT-S1/S2/S3. 7 donors with ≥ 50 nuclei per segment (4 female, 3 male); 12 at ≥ 20 nuclei as sensitivity. AKI and CKD donors were used for tissue state. |

Datasets used only in the revision analyses are described with those analyses: the many-draw nulls
in 6.14 and the state, Visium and rat protein data in 6.27.

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
- Two relabelings screen only for gross nonspecificity. Notebook 45's many-draw nulls place them in
  an external distribution (6.14), but our cohort's own NCDR still rests on two draws.
- NCDR is not a false-discovery proportion (6.14).
- Noise-corrected correlations and the conservation index have gene-bootstrap intervals (1,000
  draws). They do not resample donors of our cohort, which has two per species.
- The symmetric classes, cross-fitted label checks and coordinate reruns are sensitivity analyses
  of fixed rules. None selects a model or a threshold.

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

### 6.25 Segment-label validation (notebook 42)

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
  - Human µm distances below are corrected by this factor (−13.5%) [VERIFY against scanner metadata].
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

### 6.26 Same-species ceilings, conservation index and symmetric classes (notebook 43)

Notebook 43 runs in `results/pt_revision_classes/<labels>/` (logic `43.revision_classes.1`). Its
protocol was pre-registered in the notebook header. Labels are a parameter: the reviewed labels are
primary, and transferred labels are sensitivities (6.28).

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
zonated calls (6.20), plus two changes:
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
added after the notebook ran [VERIFY: coordinator to confirm the count].

### 6.27 Physiological state, donor tissue state, additional human donors and rodent protein (notebook 44)

Notebook 44 runs in `results/pt_revision_state/` (logic `44.revision_state.1`). Its protocol was
recorded before it ran. Every input is public except the tissue-state placement, which uses our
structures.

**Physiological state.**
- *Programs* (fixed in advance): sterol/SREBP2 (12 genes), PPARα/fatty-acid oxidation (13),
  glutathione synthesis (4), tyrosine catabolism (6), the lead genes (Gatm, Gamt, Acaa2, Dcxr, Ugt3a1,
  Acox2), and structural markers (Slc5a2, Slc5a12, Slc7a13, Slc22a7, Slc34a1, Lrp2).
- *Whole-kidney state effects:*
  - 24-h fasting against fed, GSE267280 (one pooled column per condition, so a point estimate)
    [CITE: publication for GSE267280, if any];
  - ZT18 against ZT6, GSE277302 [CITE: publication for GSE277302, if any].

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

### 6.28 Cross-fitted label sensitivity (notebook 46)

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

### 6.29 Coordinate robustness of pathway calls

Notebook 37 was run unchanged through nbclient on two further coordinates, DPT13 and SCF13-ED.
Outputs are in `results/pt_pathway_final/<label>/`, and the summary is `d4_coordinate_runs.csv`.

For the union of the robust pathways under scFates and DPT13 (40 pathways),
`coordinate_robustness.csv` records:
- robust and core status under each coordinate;
- the joint z and peak position under each;
- sign agreement of the direction at the peak, wherever the pathway is a candidate;
- the notebook 43 category.

These columns form the coordinate part of Table S2.

### 6.30 Figure assembly (notebook 41)

Notebook 41 draws every main and supplementary figure from saved results. It fits no model and runs no
test.
- Descriptive aggregates that no notebook saves (segment pseudobulks, specimen bin means, reference
  gradients) are recomputed with the source notebooks' functions. Each must reproduce saved
  statistics before it is drawn.
- Key counts are asserted against their tables.
- A layout check reports overlapping text.

Outputs go to `results/paper_figures/`, with a panel-source manifest.

### 6.31 Data and code availability

[VERIFY: deposition of the Visium HD data and segmentation polygons (controlled access for human
tissue), code archive DOI and licence.] The public datasets are listed in 6.21.

---

---

## 7. Figures and tables

These legends are those of `SCRATCH/figures/figure_legends.md` (notebook 41, round 3, commit
e23a987). They were changed only where they conflict with the saved tables or with the v0.3 text.
Each change is marked "(v0.3)" and listed in the change log. Numbers in square brackets are read from
saved tables and must be re-read if notebook 43 is rerun with other labels.

"Section 3" in the sources means notebook 41's descriptive aggregates: pseudobulk means and bin means
only, with no model and no test. Each reproduces its source notebook's saved statistics before it is
drawn.

**Standing caveat for every figure.** Two control mice (M1, M2) and two cortex sections (H1, H2) from
one male human donor. H2 is labelled medulla but is cortex. All findings are descriptive.

Pathway names carry their library: (H) MSigDB Hallmark 2020, (K) KEGG 2019 Mouse, (R) Reactome 2022.
Enrichr's Hallmark set "Pperoxisome" is shown as "Peroxisome (H)". (v0.3: the round 3 legend also
listed (W) WikiPathways, but no WikiPathways set is among the 1,513 tested pathways, Table S2.)

### Main figures

**Figure 1 | Study design and a continuous cross-species PT coordinate.**

The coordinate is the declared primary coordinate (`primary_pt_coordinate.csv`, SCF13). It is checked
to equal, structure by structure, notebook 13's saved scFates coordinate.

**a,** Design.
- Visium HD 2-µm bins were collected from two control mouse sections and two human cortex sections
  from one donor.
- Panoptic segmentation of tubule cross-sections on H&E gave tubule-by-gene profiles for 26,839
  structures.
- Reviewed clustering over two integration passes kept 24,335 nephron structures. The PT subset has
  12,866 structures (S1 5,777, S2 3,806, S3 3,283).
- A nonbranching scFates curve on the integrated embedding gives the PT coordinate (0 to 1).

*(Notebook 13: `pass2_nephron_filter_counts.csv`, `cross_species_pt_scfates.h5ad`.)*

**b,** Left: the pass-2 Harmony embedding (dimensions 1 and 2) of mouse and human PT structures,
coloured by reviewed segment. The black line is a coordinate trace: the median embedding position in
30 equal-count bins of the scFates coordinate. The dot marks the start. scFates' own curve nodes were
not saved, so the trace stands in for the curve. Right: each specimen's segment medians (circles) and
reviewed segment transitions (ticks) on the coordinate.

*(Notebook 13, pass-2 embedding cache and `pt_segment_order_by_specimen.csv`; transitions from
notebook 37.)*

**c,** Agreement of the coordinate with external zonation.
- *P1:* the per-gene, fold-protected within-specimen correlation with the coordinate, correlated
  across genes with external S3 − early contrasts. Filled bars are the coordinate; open bars use
  segment labels in place of the coordinate.
- *P1w:* within-segment agreement with adjacent external contrasts. Coordinate bars are coloured by
  species, read-depth ordering is shown in grey and dots are specimens. Read-depth ordering is not
  uniformly null: it reaches ρ = −0.16 in human early PT.

*(Notebook 35: `p1_external_concordance.csv`, `p1w_within_segment_concordance.csv`; means from
notebook 36 `arm_metrics.csv`.)*

**d,** Identifiability checks.
- Top: absorption ratio for an injected human-only gradient (1 = fully kept).
- Middle: sign agreement with external data for externally species-specific genes. Filled points are
  the coordinate and open points segment labels; n is given per set.
- Bottom: count split, with the coordinate built on read half A and tested on half B (saved, A and B
  arms). It reports robust pathways retained from the earlier list of 66 and core pathways from the
  earlier list of 42.

*(Notebook 35: `p2_absorption.csv`, `positive_control_recovery.csv`, `g6_count_split.csv`.)*

---

**Figure 2 | With two specimens per species, average-level pathway screens are not specific;
position-dependent screens are, when specimens share their positional biology.**

**a,** The three balanced 2 + 2 contrasts of four specimens (species split and two relabelings), with
the split-plot error terms.

**b,** Pathways called (BH q ≤ 0.05) per strategy, on a symlog axis.
- *Markers:* the species split is filled. Black marks a strategy that passes the specificity rule
  (NCDR < 0.25 and fewer than 5% of pathways called under relabeling); grey marks one that fails. The
  two relabelings are open squares and diamonds.
- *Columns:* NCDR, null rate and the number of calls at decoy-estimated FDP ≤ 10%.

*(Notebook 37, Table S1 `table_S1_strategy_specificity.csv`.)* (v0.3: "Table S2" corrected to Table S1,
which is the strategy table; Table S2 is the pathway table.)

**c,** Decoy-estimated FDP as the threshold is lowered, for whole-PT and S1/S2/S3 GSEA, T_level,
S1/S2/S3 steps, and T_spatial (matched and joint tests). The horizontal line marks FDP = 10%.

*(Target–decoy summary of the saved per-pathway scores from notebook 37
`matched_and_joint_tests_all_partitions.csv` and notebook 31
`conventional_signed_gsea_all_partitions.csv`. Each curve's minimum FDP is checked against Table S1.)*
(v0.3: Table S2 → S1, as in b.)

**d,** Within-species controls: offset (T_level, light) and smooth positional (T_spatial, black)
pathway calls for the species split, mouse versus mouse, and human section versus section.

*(Notebook 37, `within_species_controls.csv`.)*

**e,** Mean decoy z² against pathway size (sextiles), with the fitted inflation lines
a·(1 + (m − 1)ρ) for T_level (grey) and T_spatial (black). Whole-PT DESeq2 gene counts are shown for
reference.

*(Notebook 37 `decoy_inflation_by_statistic.csv` and matched tests; DESeq2 from notebook 31.)*

**f–h,** The many-draw specimen null, at segment resolution. Public atlases have no coordinate, so
these panels compare whole-PT and segment GSEA with the segment-step joint screen, not the smooth
coordinate. Grey marks average-level screens (whole-PT GSEA, segment GSEA) and black the segment-step
joint screen.

**f,** Pathways called (BH q ≤ 0.05) per same-species 2 + 2 split. Bars span the 5th–95th percentile,
boxes the 25th–75th, and ticks mark the median. Pools are Census mouse male (300 splits) and female
(300), Lake/KPMP human cortex (105) and Census human cortex (45). Our cohort's two relabelings are
shown as stars. [Medians: whole-PT GSEA 22, 40, 34, 13; segment steps 11, 5, 0, 0.]

**g,** NCDR per cross-species draw (100 draws of 2 human + 2 mouse donors). Bars are medians [1.34,
0.69, 0.05]. The dashed line is notebook 31's NCDR threshold (0.25), and the share of draws meeting
notebook 31's full rule is printed above each strategy [19%, 21%, 80%]. Stars mark our cohort. Draws
with no species calls have no NCDR and are not shown.

**h,** AKI positive control (two AKI against two control mice, mouse-only design). Bars give the calls
for condition and open markers the two relabelings; the number is relabel ÷ condition. Positional
screens detect injury biology (3 injury Hallmarks; 13 of 14 injury markers rise most in S3) but fail
the 5% relabeling condition (0.35 for steps, 0.17 for T_spatial), because the two AKI mice differ in
positional injury severity.

**The specificity of position-dependent screens is empirical and conditional.** It holds when the
specimens of a group share their positional biology (adult human cortex; age-matched mice). It fails
when they do not: pubertal age in mice, and injury severity in AKI. NCDR measures vulnerability to
specimen noise, not the false-discovery proportion.

*(Notebook 45, `results/pt_revision_method/tables/`: `same_species_null_calls.csv`,
`cross_species_draws.csv`, `table_aki_strategy_calls.csv` and `anchor_our_cohort.csv`; every
distribution is checked against `table_null_call_distribution.csv` and
`table_cross_species_draws_summary.csv`. See `docs/paper/r2_revision_results.md`.)*

---

**Figure 3 | Beyond a conserved core, PT zonation is largely species-specific against same-species
reliability ceilings.**

**a,** Per-gene S1 → S2 zonation: the mean log2 pseudobulk ratio of the two mice against that of the
two human sections, for probe-balanced genes classified in S2 − S1 [7,391]. Points are coloured by the
symmetric class, in which human thresholds are matched to human amplitude and species-only calls need
a flat external atlas. Class counts are in the legend [conserved 109, reversal 37, mouse-only 31,
human-only 58].

*(Classes and counts: notebook 43, `gene_classes_symmetric.csv` and
`class_counts_by_setting_and_contrast.csv`. Means: section 3 under notebook 43's labels, reproducing
notebook 43's human–mouse correlations.)*

**b,** Left: noise-corrected correlations (r*) of gene gradients for every dataset pair with
replicates, rat excluded. Pairs are grouped as within human, within mouse and human–mouse, for
S2 − S1 and S3 − early. Bars are medians.

Right: for each contrast, the conservation index with its 95% gene-bootstrap interval (filled),
together with the cross-species r* (open), the human ceiling (r* of our human sections against Lake
cortex; orange triangle) and the mouse ceiling (r* of our mice against male mouse snRNA; blue
triangle). Two versions are shown: our data, and external data (Lake against male mouse snRNA).
- The dashed lines are notebook 43's pre-registered decision thresholds for the upper 95% bound of
  the index. Below 0.5 the wording is "largely different genes", between 0.5 and 0.8 "partly
  conserved", and at 0.8 or above "largely conserved". The three zones are labelled under the axis.
- Every upper bound lies below 0.5. [Index: S2 − S1 0.17 (0.12–0.22) ours and 0.27 external;
  S3 − early 0.28 ours and 0.28 external; ceilings about 0.65 human and 0.87 mouse.]

*(Notebook 43, `pairwise_gradient_correlations.csv` and `conservation_index.csv`.)*

**c,** Mouse-only and human-only gene counts (probe-balanced) for each setting and contrast. The cell
colour is log2(mouse-only ÷ human-only): blue means a mouse excess, orange a human excess. The bold
row is the revision's primary symmetric setting; the v1 rows are the superseded notebook 40 rules.
Under the symmetric rules the mouse-only excess disappears [S2 − S1: 31 mouse-only against 58
human-only].

*(Notebook 43, `class_counts_by_setting_and_contrast.csv`.)*

**d,** The 36 robust position-dependent pathways.
- *Bars:* each pathway's contributing members (top 10% by T_spatial) split into mouse-led, shared,
  human-led, mixed, not zonated and not classifiable, as a share of contributing members.
- *Grouping:* pathways are grouped by their label (★ core). The two "Fatty Acid Metabolism" sets are
  the Hallmark (H) and Reactome (R) sets. A pathway is a summary of mouse zonation, shared or
  human-led when ≥ 50% of its zonated members fall in that category. It is mixed otherwise, and
  unresolved when it has fewer than 3 zonated members.
- *Right:* label totals [26 summary of mouse zonation, 3 shared, 0 human-led, 2 mixed, 5 unresolved].

*(Notebook 43, `robust_pathways_reframed.csv`; pathway list from notebook 37.)*

---

**Figure 4 | Robustness and agreement with independent data.**

**a,** Effect retained (AUC effect ÷ original) for the 36 robust pathways under eight checks, with the
share ≥ 0.7 at right. Core pathways are dark and robust-only pathways light. The green bar marks
power-matched random removal of mouse S3 (median 5th percentile to median).

Filled points are also robust when notebook 37's analysis is rerun unchanged on the DPT13 coordinate
[22 of 36]; open points are robust on scFates only. Peak direction agrees between the two coordinates
for every pathway robust on either [40 of 40].

*(D4 rerun, `pt_revision_labels/coordinate_robustness.csv` and `coordinate_robustness_funnel.csv`.
Notebook 37, `table_S2_pathways.csv`; DPT row from notebook 12's saved DPT pathway effects.)*

(v0.3) The row labelled "DPT coordinate (notebooks 03, 12)" is the earlier analysis on notebook 03's
DPT coordinate, not DPT13. It is historical and is not used as robustness support (Methods 6.16).
Robustness to DPT13 is shown by the filled points and in Table S2's coordinate columns.

**b,** Reviewed segment transitions on the coordinate, per specimen. *(Notebook 37.)*

**c,** Gene-level agreement. Our centred species × segment interaction is plotted against the Census
donor-level Welch t of (human − male mouse)(S3 − early). Black points are the top 10% of genes by
T_spatial. The ρ values are notebook 37's saved values, computed on the centred contrast, with the
expression-matched random median for comparison.

*(Notebook 37, `gene_level_statistics.csv`, `external_replication_summary.csv`.)*

**d,** Pathway-level replication z with donors as replicates. The groups are robust and uncalled
pathways against Census, robust pathways against Lake/KPMP, and relabel-called pathways against
Census. Lake z for uncalled pathways is not saved. *(Notebook 37, Tables S2/S3.)*

**e,** Specimen-level split-plot gene test, Q–Q plot for the species split and both relabelings
(relabeling λ ≤ 0.87). *(Notebook 37, split-plot cache; species p-values checked against
`gene_level_statistics.csv`.)*

---

**Figure 5 | Lead genes against direct segment measurements.** The genes are Gatm, Gamt, Acadm, Acaa2,
Dcxr, Ugt3a1, Gclc and Cyp2e1, plus Slc7a13, a conserved apical transporter shown as a reference (no
claim). Gamt is shown only as the creatine partner of Gatm.

- *Left of each card:* fitted human (orange) and mouse (blue) curves along the PT coordinate, with each
  specimen's mean in 12 coordinate bins (circles M1 and H1, squares M2 and H2). Shading marks each
  reviewed segment's 10th–90th percentile. No curve is fitted for the figure.
- *Right:* each external dataset's S1, S2 and S3 values centred on its own S1 (log2):
  - Lake human cortex (7 donors);
  - mouse microdissection, male and female;
  - mouse snRNA, male and female.

  Rat is omitted because it fails notebook 43's reliability rule.
- *Card titles:* gene, story, the pre-specified verdict ("no claim" for Slc7a13), the symmetric
  zonation class from notebook 43 (union of S2 − S1 and S3c − early), and probe counts (mouse/human).
  [Classes: Gatm, Gamt and Gclc indeterminate; Acadm and Slc7a13 conserved; Acaa2, Dcxr and Ugt3a1
  reversal; Cyp2e1 not classifiable, because it is not expressed in human.]

*(Sources:*
- *curves: notebook 37's saved species fit, checked against `gene_level_statistics.csv` and identical
  to notebook 38's story-gene curves;*
- *specimen bins, microdissection TPM and mouse snRNA segment means: notebook 41 section 3,
  reproducing notebook 37's grid, notebook 38's detection, its saved TPM and its snRNA S3 − S1
  contrasts;*
- *Lake: notebook 38 `lake_cortex_segment_means.csv`;*
- *verdicts: `story_gene_evidence.csv`;*
- *probes: `probe_counts.csv`.)*

### Main table

**Table 1 | Conservation index and symmetric zonation classes.** (v0.3, new)
- *Part A:* the conservation index for each contrast and version (our data, external only,
  alternative ceilings), with 95% intervals, the cross-species r* and both ceilings.
- *Part B:* class counts (conserved, reversal, mouse-only, human-only, indeterminate, neither) for
  every setting and contrast. The settings are v1, A6, A7, A3b, A8, symmetric, and symmetric with A6
  or A7; the contrasts are S2 − S1, S3c − early, S3 − early and the union. Probe-balanced genes only.

*(Notebook 43: `conservation_index.csv`, `class_counts_by_setting_and_contrast.csv`.)*

### Supplementary figures

**Figure S1 | Coordinate limitations.**

**a,** Thinning to equal depth.
- Left: agreement of the equal-depth refit with the original order, per specimen and segment.
- Right: correlation of position with log library size before (dot) and after (arrow tip) thinning.
  Human S1 is depth-dependent.

**b,** Largest human − mouse segment shift across refits, for scFates and DPT.

**c,** Left: within-segment rank agreement between coordinates built from disjoint gene folds (bar:
median). Right: the largest position change after a 10⁻¹⁶ input change (numerical floor).

**d,** Gain in held-out gene prediction in the other specimen over the segment oracle (circles) and the
coverage oracle (diamonds), per coordinate arm and species.

**e,** Left: amplitudes of 38 conserved held-out genes in human against mouse. Filled points use the
coordinate and open points segment labels. Shares with smaller human amplitude are 76% (coordinate),
71% (labels) and 79% (external). Right: the human-high and mouse-high pathway counts at peak, **in the
earlier list of 66 robust pathways**, after adding positional noise to mouse. The dotted line marks
the estimated excess noise.

(v0.3) On the final 36 robust pathways the same procedure accounts for about 15% of the peak gap and
7–8% of the mouse-high count (R1).

*(Notebook 35: `p5_exposure_invariance`, `species_registration_across_refits`,
`p3a_gene_fold_within_segment_agreement`, `numerical_floor`, `a3_conserved_heldout_amplitudes`,
`a3b_noise_matching_summary`. Notebook 36: `p4_heldout_prediction`, `early_late_split_by_coordinate`,
`attenuation_by_arm`.)*

**Figure S2 | DPT specificity diagnostic (descriptive).**

**a,** Categories of called pathways for DPT and scFates. For the relabelings, the calls of both
relabelings are **summed**, not unioned: 77 = 41 + 36 for DPT, against a union of 57.

**b,** Within-species positional pathway calls (DPT: mouse 36, human 5; scFates: 16, 0).

**c,** Pooled within-segment IQR of each coordinate.

*(DPT diagnostic, `pt_reconstruction_v2/dpt_diagnostic/`.)*

(v0.3) These panels use the matched-only test. On notebook 37's joint test, DPT13 is specific: 45
species calls against 3 and 0 relabeled (NCDR 0.03; `d4_coordinate_runs.csv`, Methods 6.10).

**Figure S3 | Position-dependent pathways.**

**a,** Joint-test funnel (of 1,513 tested).

**b,** Median member Δ(s) (human − mouse, log-normalised) for the 36 robust pathways in 17 programs.
Ticks mark peaks and ★ marks core pathways. Annotations show the program, Census external replication
and the split-plot matched test.

**c,** Peak position against the median member Z at the peak.

**d,** The six largest programs.
- Top: member Δ(s) (grey) and the program mean, coloured by direction at the peak.
- Bottom: a **direction-matched driver**, which is the first shared driver in T_spatial order whose
  fitted Δ at the program's median peak has the program's sign. Its fitted curves are shown with each
  specimen's 12-bin means. When the mouse and human symbols differ, both are given (Cyp3a11 /
  CYP3A5).

**e,** Averaging loss by shape class.

*(Notebook 37: funnel, Table S2, Table S4, fitted-curve cache. Notebook 31: gene sets, checked against
the member counts and shape classes. Section 3: specimen bins.)*

**Figure S4 | Zonation amplitude.**

**a,** Human ÷ mouse amplitude against each reference: OLS with 95% interval (●), Deming (◇) and
orthogonal (□). The dashed line is the primary ratio r. External-only ratios are in green. (v0.3) With
a low cross-species correlation these are projection ratios onto each reference axis (Methods 6.19).

**b, c,** Both species' S2 − S1 gradients against male microdissection (b) and Lake healthy cortex (c),
with the saved OLS slopes. *(Gradients from section 3, which reproduces notebook 39's slopes and gene
counts.)*

**d,** Measurement checks: expression tertiles, log1p means and equal probe counts.

**e,** Tissue state: human ÷ mouse, AKI ÷ control mice, and low- and high-injury halves.

**f,** Lake donors' injury score against S2 − S1 amplitude. Spearman ρ is computed from the saved donor
table.

**g,** Lead genes against global flattening (dashed line: human = r × mouse).

**h, i,** Post hoc: noise-corrected SD ratios and disattenuated correlations.

*(Notebook 39 tables.)*

**Figure S5 | Lead-gene verification matrix.** One row per claimed gene and one column per
pre-specified test.
- *Cells:* dark ✓ passes; open × fails; grey not testable; blank not claimed.
- *Right:* the verdict, probe counts and annotations.
- *Counts:* 39 claimed genes, of which 30 are headline, 5 supporting and 4 not replicated.

(v0.3) The verdicts are notebook 38's and predate the symmetric classes. Several headline genes are
indeterminate under notebook 43 (for example Gatm, Gclm, Acsm3 and Acox2), so the text states each
gene's class (R5).

*(Notebook 38, `story_gene_evidence.csv`, `run_manifest.json`.)*

**Figure S6 | Opposite regional differences that a whole-PT comparison averages away (Slc22a6,
Slc13a3, Cyp24a1).**

**a,** Fitted human and mouse curves with each specimen's 12-bin means.

**b,** DESeq2 log2 fold changes (human − mouse) for whole PT and each segment (* padj ≤ 0.05).

In all three genes the whole-PT difference is not significant, but segment differences of opposite
sign are. **Caveat:** the S3 bars compare human cortex S3 with mouse S3 that includes outer stripe,
which has no human counterpart, so the S3 differences are partly regional. Slc22a7 is no longer shown,
because its whole-PT difference is itself significant.

*(Notebook 37 fitted-curve cache; DESeq2 from notebook 31, checked against notebook 38; bins from
section 3.)*

**Figure S7 | Segment-label agreement with the cortex atlas (notebook 38, genome-wide Lake check).**

**a,** Our human S3 − S1 against Lake cortex S3 − S1 (ρ, genes, direction agreement for strongly
zonated genes).

**b, c,** Our species × segment interaction against Lake − male mouse snRNA (b) and Lake − female mouse
snRNA (c), for all shared genes and for the top 10% by T_spatial (black).

Notebook 34's earlier Census-based version is superseded. No mouse within-species panel is shown.

*(Notebook 38, `gene_level_lake_contrasts.csv`; ρ and gene counts are rebuilt from it and match
`genome_wide_lake_check.csv`.)*

**Figure S8 | Superseded: notebook 40's zonation classes (Figure 3 of draft v0.1), kept for the
record.**

**a,** Per-gene S2 − S1 in each species, coloured by notebook 40's v1 class (65/20/219/44
conserved/reversal/mouse-only/human-only). The means come from section 3.

**b,** v1 class sizes under the primary rule and seven alternatives.

**c,** v1 representative genes, as specimen log2 CPM at S1, S2 and S3c (section 3).

**d,** Programs enriched in each v1 class, named by the class's best pathway in the program.

**e,** Robust pathways by v1 member class.

The mouse-only excess is a threshold effect of lower human amplitude (Figure 3c). *(Notebook 40
tables.)*

**Figure S9 | Cross-fitted label sensitivity of the conservation index (notebook 46).** Genes were
split into two halves. Segment labels were transferred using one half and the index was evaluated on
the other half, so labels and evaluation share no genes. The label settings are (i) human labels
transferred from Lake with mouse reviewed labels, (ii) human from Lake with mouse from
microdissection, and (iii) reviewed labels.

**a,** The conservation index (point) with its 95% interval and the upper bound marked, per setting
and gene half (half 0 → 1 means labels from half 0, evaluated on half 1), for S2 − S1, S3c − early
and S3 − early. The dashed line is the 0.5 threshold, and grey lines mark the reviewed-label index on
all genes, which notebook 46 reproduces exactly from notebook 43. [Largest index for S2 − S1 and
S3 − early 0.31, largest upper bound 0.38: the claim stands under every setting.] (v0.3) For the
region-matched S3c − early contrast, which the panel also shows, the largest index is 0.33 and the
largest upper bound 0.43, in setting (ii); both remain below 0.5.

**b,** Human label sets: S1/S2/S3 fractions and agreement with the reviewed labels [75–76% for
Lake-transferred labels].

*(Notebook 46, `pt_revision_classes/crossfit/tables/`: `crossfit_conservation_index.csv`,
`crossfit_label_sets.csv`, `reproduction_reviewed_all_genes_index.csv`, `decisions.csv`.)*

**Figure S10 | Anatomical and reference-transfer checks of the PT segment labels (notebook 42).**

**a,** Distance of each PT structure to the nearest glomeruli by reviewed label and species: median
and interquartile range. Filled points use morphological glomeruli (primary) and open points
expression-labelled glomeruli. [Spearman ρ of label order with distance: mouse 0.65, human 0.23.]
(v0.3) The mouse morphological glomerulus polygons fail notebook 42's pre-registered size rule
(1.3–1.4× the PT median area, against ≥ 3×), so the mouse distances are descriptive.

**b, c,** Reference transfer. Labels transferred from Lake/KPMP human cortex (b) and from male mouse
microdissection (c) are compared with our reviewed labels. Cells give counts and shading the row
share. [Agreement 71%, κ = 0.54 in human; 87%, κ = 0.80 in mouse.] Row totals equal the reviewed
label counts.

**Scale caveat:** the human distances assume 0.5 µm per pixel, against about 0.44 for the mouse
sections. If the human scale matches the mouse scale, the human distances are overstated by 13.5%.
The human pixel scale is unconfirmed [VERIFY against scanner metadata].

*(Notebook 42, `pt_revision_labels/`: `glomerulus_distance_by_label.csv`,
`transfer_confusion_human_lake_primary.csv`, `transfer_confusion_mouse_microdissection_primary.csv`,
`transfer_agreement.csv`, `summary.json`.)*

**Figure S11 | Physiological and tissue state (notebook 44).** (v0.3, proposed. Notebook 44 saves this
figure, but notebook 41 does not yet assemble it.)

**a,** Whole-kidney state effects. For each program, the median largest fasting or time-of-day effect
against the median |species × position interaction|. The dashed line marks the 0.5 decision ratio.

**b,** Mouse snRNA S3 − S1 at ZT4 against ZT16 for program genes.

**c,** Our donor on the Lake injury axis. Lake donors are coloured by disease; our two sections are
shown after mouse-calibrated platform correction.

**d,** Visium donors (Abedini et al. 2024), whose label check failed. Lead-gene S3 − S1 per sample
(post hoc centred) beside Lake and mouse snRNA values.

*(Notebook 44, `results/pt_revision_state/figures/fig_revision_state.pdf` and `tables/`.)*

### Supplementary tables

| Table | Content | Source |
|---|---|---|
| S1 | Strategy specificity: species and relabeled calls, NCDR, null rate, decoy calls with interval, minimum FDP | `pt_pathway_final/scfates/tables/table_S1_strategy_specificity.csv` |
| S2 | All 1,513 pathways: joint test, specificity margin, every robustness ratio, stability, program, shape and external replication. (v0.3) Coordinate columns for the 40 pathways robust under scFates or DPT13: robust and core status under SCF13, DPT13 ("DPT13-robust") and SCF13-ED; joint z and peak position under each; whether the direction at the peak agrees across coordinates; notebook 43 category. These columns still have to be merged into the table file. | `table_S2_pathways.csv`; `pt_revision_labels/coordinate_robustness.csv`; `d4_coordinate_runs.csv` |
| S3 | External replication per pathway | `table_S3_external_replication.csv` |
| S4 | The 17 programs with drivers and member pathways | `table_S4_programs.csv` |
| S5 | External replication summary and global component | `table_S5_global_component.csv`, `external_replication_summary.csv` |
| S6 | Coordinate arm comparison and decision rules (SCF13, SCF13-ED, DPT13, SCF-PT, CAC), including the post-hoc specificity gate and the D4 reruns | `pt_reconstruction_v2/nb36/arm_metrics.csv`, `c6_decision.csv`, `selection_rule.csv`, `downstream_summary.csv`; `pt_revision_labels/d4_coordinate_runs.csv` |
| S7 | Probe counts per gene in both panels; symmetric zonation classes for all genes; genes expressed in one species only; superseded notebook 40 classes | `pt_literature_deep/probe_counts.csv`; `pt_revision_classes/reviewed/tables/gene_classes_symmetric.csv`; `pt_zonation_classes/tables/gene_zonation_classes.csv`, `one_species_expressed_genes.csv` |
| S8 | Segmentation performance (PQ, SQ, DQ, merge and split rates) | `segmentation/model_summary_v2.txt` [VERIFY: not produced by a notebook] |
| S9 | Many-draw nulls: same-species null call distributions, cross-species draw summary, AKI strategy calls, coherence by pool | `pt_revision_method/tables/table_null_call_distribution.csv`, `table_cross_species_draws_summary.csv`, `table_aki_strategy_calls.csv`, `table_mechanism_coherence.csv` |
| S10 | Label validation: transfer agreement and confusion, glomerulus distance and attachment, spillover, registration QC; cross-fitted conservation index and class counts | `pt_revision_labels/*.csv`; `pt_revision_classes/crossfit/tables/*.csv` |
| S11 | State checks: whole-kidney state effects, time-of-day zonation, donor injury placement, Visium lead genes, rat segment proteome | `pt_revision_state/tables/*.csv` |

---

## 8. References

All references were checked in Europe PMC, Crossref or on the publisher page by workstream 4. Those
marked † were also verified by workstream 3 (`SCRATCH/ws3/literature_phase2.md`). The 14 references
added in v0.3 were checked by PubMed ID in Europe PMC.

Public datasets without a citable publication are identified by accession in Methods (GSE267280,
GSE277302) [CITE: publications, if any].

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

---

## 9. Remaining open items

The full list, with what each item needs and who can resolve it, is `open_items_v03.md` (proposed
replacement for `docs/paper/open_items.md`). The point-by-point referee response is
`referee_response_v01.md`. Every [PENDING D4] and [PENDING D8] marker of v0.2 is resolved in this
version. The markers still in the text:

**[VERIFY]: needs lab records or the user**
1. Specimens: mouse strain, age, supplier, housing, euthanasia; sex from records (6.1).
2. Human donor: source, age, kidney function, procurement, serial or separate blocks (6.1;
   Limitations).
3. Tissue processing, imaging, Visium HD workflow and sequencing (6.1, 6.2).
4. Ethics approvals (6.1).
5. Segmentation provenance: human `v2` segmentations, upstream mouse QC, H&E–Visium registration,
   whether to report PQ, inference software versions (6.3, 6.4, 6.24; Table S8).
6. **Pixel scale.** The human GeoJSON µm fields assume 0.5 µm per pixel; the segmentation runs at
   0.44068. Confirm against scanner metadata (R1, 6.3, 6.25, Figure S10).
7. HCOP download date (6.5); pathway library source and date (6.13); θ = 6 rationale (6.7).
8. AI-agent wording for the literature check (6.22).
9. Data and code deposition (6.31).

**[VERIFY]: computed outside a notebook**
10. The transporter shares (33 of 43 zonated transporters conserved, against 196 of 388 zonated
    genes) were counted by workstream 4 from `gene_classes_symmetric.csv` (R3, 6.26).
11. UGT1A9 probe positions relative to the unique exon 1 (R5).

**[CITE]**
12. OpenMidnight (Kaplan et al. 2025): model card only.
13. Publications for GSE267280 and GSE277302, if any (6.27).

**Figures (after round 3)**
14. Figure 5 draws the v0.3 gene set. Recommended: Gclm in place of Gclc, and possibly Rbp4
    (open items, D).
15. Figure 4a's historical DPT row should be replaced or supplemented by a DPT13 row.
16. Figure S11 (notebook 44 state figure) is proposed; notebook 41 does not yet assemble it.
    Figure 2f needs our relabelings on their strategy rows (`v02_changes.md` item 15).
17. Table S2's coordinate columns must be merged from `coordinate_robustness.csv`.
