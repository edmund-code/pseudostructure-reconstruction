# Citation check, references [1]–[21] (Introduction of manuscript_v1_intro_methods.md)

Checked 2026-11-03. Downloads are in `fulltext/dl_A/` (plain-text conversions in `dl_A/txt/`). No personal identifier was sent to any service. Section labels in quotes refer to the paper's own headings.

## Summary table

| Ref | Paper | Access route used | Full text read | Verdict | Proposed fix |
|---|---|---|---|---|---|
| 1 | Ståhl 2016 *Science* | Published PDF (6 pp.) on a university course-reading page; science.org 403, no PMC | Y (main text + legends; online SM not read) | SUPPORTED | none |
| 2 | Marx 2021 *Nat Methods* (Technology Feature) | nature.com paywall; no PMC, no preprint | **N — INACCESSIBLE, standfirst only** | NOT VERIFIABLE | Replace with the open editorial *Method of the Year 2020: spatially resolved transcriptomics*, Nat. Methods 18, 1 (2021), doi:10.1038/s41592-020-01042-x (read in full), or drop [2] |
| 3 | Moses & Pachter 2022 *Nat Methods* | bioRxiv preprint v2 (25 May 2021) PDF via CaltechAUTHORS; published version paywalled | Y (preprint only) | SUPPORTED (preprint) | none |
| 4 | Oliveira 2025 *Nat Genet* | Europe PMC XML PMC12165841 | Y | PARTLY | "whole sections" is wrong: capture area is 6.5 × 6.5 mm; see rewrite S1 |
| 5 | Halpern 2017 *Nature* | Europe PMC XML PMC5321580 (author manuscript) | Y | PARTLY (lobule yes; nephron not covered) | add a nephron ref; "along their length" → "along an axis"; see S2 |
| 6 | Moor 2018 *Cell* | cell.com 403, bioRxiv 429 (Cloudflare), institutional copy 403 | **N — INACCESSIBLE, abstract only** (Cell and preprint abstracts) | abstract supports villus; provisional | re-read before submission; see S2 |
| 7 | Curthoys & Moe 2014 *CJASN* | PMC HTML PMC4152816 | Y | PARTLY (transport, metabolism yes; injury no) | add injury source; see S3 |
| 8 | Lee 2015 *JASN* | PMC HTML PMC4625681 | Y (tables not rendered) | NOT SUPPORTED for this sentence (rat; no S1/S2/S3 differences discussed) | move [8] to the "repeated units" sentence (S2) |
| 9 | Chen 2021 *JASN* | PMC HTML PMC8017530 | Y | Intro: PARTLY (transport/metabolism yes, injury no). Methods: SUPPORTED | see S3 |
| 10 | Ransick 2019 *Dev Cell* | PMC HTML PMC6948019 (author manuscript) | Y | SUPPORTED | can also carry the S3-injury clause in S3 |
| 11 | Muto 2021 *Nat Commun* | Europe PMC XML PMC8044133 | Y | SUPPORTED | none |
| 12 | Lake 2023 *Nature* | Europe PMC XML PMC10356613 | Y | SUPPORTED | none |
| 13 | Kirita 2020 *PNAS* | Europe PMC XML PMC7355049 (SI Appendix not read) | Y | SUPPORTED | none |
| 14 | Marshall 2022 *iScience* | Europe PMC XML PMC8971939 | Y | NOT SUPPORTED (no PT segment states, no AKI model) | replace with Melo Ferreira 2021 *JCI Insight* 6, e147703, doi:10.1172/jci.insight.147703 (full text read, Europe PMC PMC8262485) |
| 15 | Dixon 2022 *JASN* | PMC HTML PMC8819997 | Y | SUPPORTED | none |
| 16 | Zhai 2003 *JASN* | Publisher full-text HTML (Ovid/LWW JASN, free) | Y (tables are images) | SUPPORTED | see S5 wording |
| 17 | Zhai 2006 *JASN* | Publisher full-text HTML (Ovid/LWW JASN, free) | Y (tables are images) | SUPPORTED | see S5 wording |
| 18 | Wang 2018 *Science* (STARmap) | PMC HTML PMC6339868 (author manuscript; SM not read) | Y | SUPPORTED | see S5 wording |
| 19 | Vickovic 2022 *Commun Biol* | Europe PMC XML PMC8837632 | Y | SUPPORTED | none |
| 20 | Schott 2024 *Cell* (Open-ST) | Published version (CC BY-NC) PDF, MDC repository edoc 24484 | Y | PARTLY (3D whole-transcriptome yes; it presents itself as low-cost, not laborious) | see S5 wording |
| 21 | Kiemen 2022 *Nat Methods* (CODA) | PMC HTML PMC10500590 (author manuscript) | Y | PARTLY (serial H&E reconstruction, no expression profiling) | regroup with serial-section tracing; see S5 |

### Proposed rewrites of the citing sentences

- **S1 (refs 1–4).** "Spatial transcriptomics (ST) technologies measure gene expression in intact tissue sections while keeping each measurement's position[1–3], and Visium HD captures the protein-coding transcriptome on a continuous lawn of 2-µm barcoded squares across a 6.5 × 6.5-mm capture area[4]." The product name is "Visium HD" (two words) in ref. 4.
- **S2 (refs 5, 6).** "Many organs are built from repeated units whose function changes continuously along an axis, such as the liver lobule[5], the intestinal villus[6] and the kidney nephron[8,10]." Lee 2015 and Ransick 2019 both describe the nephron as a repeated unit (quotes below). Moses 2022 [3] could also be cited here ("liver lobules, intestinal villi").
- **S3 (refs 7–9).** "Along the nephron, the proximal tubule runs through three subsegments, S1, S2 and S3, which differ in solute transport and metabolism[7,9] and in susceptibility to ischaemic injury[10,NEW]." NEW = Venkatachalam, M. A., Bernard, D. B., Donohoe, J. F. & Levinsky, N. G. Ischemic damage and repair in the rat proximal tubule: differences among the S1, S2, and S3 segments. *Kidney Int.* 14, 31–49 (1978), doi:10.1038/ki.1978.87, PMID 682423. **I checked only the abstract**, which says "S1 and S2 cells alike are only reversibly injured ... whereas S3 cells selectively undergo progressive cell injury and death"; read the full text before citing it. If no primary source is wanted, Ransick 2019 [10] alone carries the clause (full text, quoted below).
- **S4 (refs 10–15).** Keep the sentence and swap [14] Marshall 2022 for Melo Ferreira 2021. Alternatively keep Marshall and change "in acute injury" to "in acute injury and other kidney diseases", but Marshall still maps no PT segment states.
- **S5 (refs 16–21).** "Recovering that position physically requires tracing individual nephrons through hundreds to thousands of serial sections[16,17] or reconstructing tissue in three dimensions[18–21]; these approaches need many sections or assays per specimen and so far capture morphology[21], targeted gene panels[18] or transcriptomes over sub-millimetre depths[19,20], not whole organs."

## Factual details in our draft that these papers contradict or qualify

1. **Visium HD bin size and what Oliveira 2025 did.** The array is "2 × 2-µm squares ... directly adjacent to each other", 6.5 × 6.5 mm per capture area. Space Ranger 3.0 outputs raw 2-µm data plus 8- and 16-µm bins, and the authors "used the 8-µm binned data" for most analyses. They used 2-µm bins only with StarDist nuclear segmentation, to locate T cells. The assay is probe-based ("ligated probe pairs targeting the whole protein-coding transcriptome"), applied to FFPE human colorectal cancer. "Measure whole sections" overstates it: one capture area is 6.5 × 6.5 mm, and the authors say "across a large tissue section". Our Methods step of summing 2-µm bins inside H&E polygons is essentially their nuclear-segmentation custom binning, so ref. 4 can also be cited in Methods.
2. **Lee 2015 is rat.** The rats were "Male Sprague–Dawley". The paper treats PT segments as a group ("the proximal tubule segments stand out as having a unique profile") and does not discuss S1/S2/S3 differences in transport, metabolism or injury.
3. **Curthoys & Moe 2014** describe S1/S2 vs S3 differences in transport (SGLT2 vs SGLT1, earlier Pi and HCO3⁻ uptake) and in metabolism (acidosis-induced glutamine catabolism in the proximal convoluted tubule; SNAT3 moves from S3 to S1/S2 in chronic acidosis). They **do not discuss injury susceptibility** at all; "injur", "ischem" and "nephrotox" do not appear in the body.
4. **Zhai 2003 and 2006 did trace nephrons through serial sections.** In 2003 they cut about 800 serial 2.5-µm Epon sections per kidney (surface to inner stripe), traced 160 PTs (11 by hand on paper copies) and discarded paths leaving the matrix. In 2006 they cut 2,500 consecutive 2.5-µm sections per kidney, digitized every second one, placed lumen markers by hand, and traced 200 nephrons and CDs; 17 of 41 long-looped nephrons were incomplete because they left the block. **Neither paper states time or money cost**, so "laborious" is an inference from the section counts. The work is morphology only, in C57BL/6J males.
5. **Zhai 2003 qualifies our "three subsegments" sentence for mouse.** It found "no parameters enabled definition of distinct segments in this strain of mice" (ultrastructure). It adds that this does "not exclude a possible segmentation based on variations in enzyme/transporter-protein expression". In mouse, S1/S2/S3 are therefore defined by transporter and marker expression and location, not by ultrastructure. This supports our continuous framing and is worth a clause.
6. **Kirita 2020 and Dixon 2022 did map PT segment states in injury; Marshall 2022 did not.**
   - Kirita (snRNA-seq; male mice; 18-min bilateral clamp) has "injured S1/2", "injured S3", "severe injured PT" and FR-PTC states.
   - Dixon (Visium; **female** mice; a 34-min clamp matched 22 min in males) maps PTS1–3, InjPT and FR-PTC signatures over time.
   - Marshall (Slide-seqV2) studies a diabetic model (BTBR ob/ob), a UMOD proteinopathy model and nine human nephrectomies, one with "ischemia due to tumor compression". Its only PT label is "PCT".
   - Note that Kirita's clamp time (18 min) differs from the "22 minutes in male mice" that Dixon cites as their own male calibration.
7. **Halpern 2017 and Moor 2018 say nothing about the kidney.** Also, the lobule axis is radial (porto-central), not a "length".
8. **Open-ST does not call itself laborious.** It states "<€130 per 12 mm² capture area" and that "a single researcher can prepare 10–15 libraries in 3 days". The PDF renders € as "V" and µm as "mm". Its 3D block is 19 profiled sections over 350 µm of one lymph node.
9. **CODA profiles no expression.** It uses H&E only and lists spatial 'omics on intervening sections as future work. It is a serial-section reconstruction like Zhai, not "profiling tissue in three dimensions" in the molecular sense.
10. **Ransick 2019 points to a sex issue that matters for us.** S3 shows "the highest disparity in gene activity between the sexes", and S2/S3 sex-specific domains form "a tight, non-overlapping boundary at the S2-S3 junction". All our specimens are male, so this belongs in Methods or limitations whenever S3 genes are compared.

---

## [1] Ståhl, P. L. et al. *Science* 353, 78–82 (2016)

**Access.** Published Science PDF (6 pages; main text, figure legends and references), posted as a reading on a university course page (aliyxcia.github.io). science.org returned 403; there is no PMC copy. The online Supplementary Materials and Methods were not read.

**Citing sentence (Introduction).** "Spatial transcriptomics (ST) technologies enable gene expression characterization within intact tissue sections[1–3], and recent platforms such as VisiumHD measure whole sections at 2-µm resolution[4]."

**Verdict: SUPPORTED.**
- Main text, p.78: "introduce positional molecular barcodes in the complementary DNA (cDNA) synthesis reaction within the context of an intact tissue section before RNA-seq."
- Main text, p.78: "~200 million oligonucleotides in each of 1007 features, with a diameter of 100 µm and a center-to-center distance of 200 µm".
- Main text, p.82: "in contrast to when different regions of a tissue are dissected for analysis, the information for the whole section is maintained".

**Writing ideas.**
- *Introduction / Methods (how to state resolution).* It gives feature diameter, pitch and area together: "a diameter of 100 µm and a center-to-center distance of 200 µm, over an area of 6.2 mm by 6.6 mm". Give Visium HD the same treatment (2 µm squares, contiguous, 6.5 × 6.5 mm).
- *Introduction (the 2D framing).* It names the limitation we build on: "high-quality RNA-sequencing data with maintained two-dimensional positional information" (Abstract). ST is explicitly 2D per section.
- *Validation.* They benchmark sensitivity against an orthogonal gold standard: "The sensitivity of spatial transcriptomics was 6.9 ± 1.5% of single-molecule fluorescent in situ hybridization". This is a model for stating a calibrated, quantitative check.
- *Discussion.* "An individual spatial transcriptomics experiment thus serves as a permanent resource". This supports our pitch that reconstruction reuses existing ST data without new experiments.

## [2] Marx, V. *Nat. Methods* 18, 9–14 (2021)

**Access. INACCESSIBLE — standfirst only.** nature.com shows a paywall ("This is a preview of subscription content"), and the PDF link redirects to an institutional login. There is no PMC copy or preprint, and Europe PMC has no abstract. Only the standfirst is visible: "Nature Methods has crowned spatially resolved transcriptomics Method of the Year 2020." A Publisher Correction exists (doi:10.1038/s41592-021-01065-y).

**Citing sentence.** Same as [1].

**Verdict: NOT VERIFIABLE.** The standfirst supports only that spatial transcriptomics was named Method of the Year.

**Replacement, read in full.** *Method of the Year 2020: spatially resolved transcriptomics.* *Nat. Methods* 18, 1 (2021), doi:10.1038/s41592-020-01042-x. This is the open, unsigned editorial for the same special issue; nature.com HTML was read in full.
- Editorial, para 2: "these methods are able to elucidate single-cell heterogeneity and define cell types while also retaining spatial information."
- Editorial, para 3: "These approaches are often complementary and differ in their target coverage, spatial resolution and throughput".

Alternatively drop [2]: [1] and [3] carry the claim.

**Writing ideas (from the editorial).**
- *Introduction.* For framing the gap at organ scale: specialized cell types "remain poorly explored on the scale of whole tissues and organisms".
- *Introduction / Discussion.* Its two-way split of the field ("read out by microscopy ... Alternatively, RNA is captured in a way that retains spatial information") is a compact way to say our method is platform-agnostic.
- *Discussion (future work).* It anticipates "software tools for automated processing of the data and for dealing with huge datasets proving important". This is a citable motivation for a computational method.

## [3] Moses, L. & Pachter, L. *Nat. Methods* 19, 534–546 (2022)

**Access.** bioRxiv preprint v2 (posted 25 May 2021; 31-page PDF; CC BY-NC), downloaded from the CaltechAUTHORS record. bioRxiv itself returned 429 (Cloudflare), and the published version is paywalled. Differences I can see: the published abstract says "To comprehensively document ... we present a curated review", while the preprint says "To assess the ability and potential ... we present a curated database". The opening sentence is identical. I could not compare the published body.

**Citing sentence.** Same as [1].

**Verdict: SUPPORTED (preprint).**
- Section 4, array-based techniques: "Array based techniques have been applied to large areas of tissue ... Nevertheless, they do not have single cell resolution".
- Section 5, Future prospective: "most current era techniques require tissue sections."
- Abstract: "The function of many biological systems, such as embryos, liver lobules, intestinal villi, and tumors depends on the spatial organization of their cells."

**Writing ideas.**
- *Introduction (S2).* The abstract sentence above names liver lobules and intestinal villi together, so [3] can co-support the repeated-units sentence.
- *Introduction (3D gap).* "more efficient computational tools for aligning multiple sections ... should be developed to extend current era techniques to 3D" (Section 5). This is a direct call for computational routes to 3D; we offer one that needs no alignment.
- *Methods / Introduction (resolution wording).* They state resolution as spot diameter and pitch per platform: "the spot diameter of ST is 100 µm with spots 200 µm apart"; "Bead diameter is 10 µm in Slide-seq, and 2 µm in HDST".

## [4] Oliveira, M. F. et al. *Nat. Genet.* 57, 1512–1523 (2025)

**Access.** Europe PMC full-text XML, PMC12165841. I read the Main text, Results, Discussion, Methods and figure legends.

**Citing sentence.** "... recent platforms such as VisiumHD measure whole sections at 2-µm resolution[4]."

**Verdict: PARTLY SUPPORTED.** The 2-µm feature size is right. "Whole sections" is not: the capture area is 6.5 × 6.5 mm, and the authors analysed mostly 8-µm bins.
- Main: "~11,000,000 continuous 2-µm features in a 6.5 × 6.5-mm Visium HD capture area, compared to ~5,000 55-µm features with gaps".
- Results, Visium HD specifications: "the Space Ranger (version 3.0) pipeline outputs the raw 2-µm data and data binned at 8- and 16-µm resolution. Unless otherwise described, we used the 8-µm binned data".
- Discussion: "with micrometer-level precision across a large tissue section".

**Fix.** Use rewrite S1. If "whole section" is meant for our own data, say so only if our mouse and human sections fit inside one capture area, which needs checking against our slides.

**Writing ideas.**
- *Methods, Aggregating bins into units.* Our step has a direct precedent: "we followed the analysis guide 'Nuclei Segmentation and Custom Binning of Visium HD Gene Expression Data'"; "a gene-by-nucleus UMI count matrix" (Methods: Nuclear segmentation). Cite [4] there as well.
- *Introduction / Discussion (why segmentation-based aggregation).* They state the gap we fill: "truly accurate cell type assignment requires advanced cell segmentation algorithms that are not yet available" (Discussion).
- *Methods (stating bin choice).* On why 8 µm: "bins are small enough to provide single-cell scale for most cell types" (Methods: Visium HD). This is a model for justifying our choice of tubule cross-sections over fixed bins.
- *Discussion, limitations.* Plain phrasing to borrow: "A limitation of this study is the small sample size, preventing us from drawing broad conclusions." It fits our 2-vs-2 design.

## [5] Halpern, K. B. et al. *Nature* 542, 352–356 (2017)

**Access.** Europe PMC full-text XML, PMC5321580 (UKPMC Funders author manuscript): main text, Methods and Extended Data legends. The published surname is "Bahar Halpern"; an erratum exists, *Nature* 543, 742 (2017).

**Citing sentence (Introduction).** "Many organs are built from repeated three-dimensional units whose function changes continuously along their length, such as the liver lobule, the intestinal villus and the kidney nephron[5,6]."

**Verdict: PARTLY SUPPORTED.** The liver-lobule part is supported. The paper does not mention the kidney, and the lobule axis is radial, not a "length".
- Main text, para 2: "The liver is a heterogeneous tissue that consists of hepatocytes operating in repeating anatomical units termed lobules."
- Main text: "most liver genes exhibit spatially-graded rather than binary expression patterns."
- Main text: "around 50% of the expressed liver genes are non-randomly spatially zonated".

**Fix.** Use rewrite S2.

**Writing ideas.**
- *Discussion (continuous beats discrete).* Their exact rhetorical move is ours: "Our study challenges the traditional binary classification of liver into periportal and pericentral hepatocytes and reveals multiple roles for the intermediate lobule coordinates."
- *Results (non-monotonic genes).* "a new class of non-monotonic zonation profiles that peak at intermediate lobule layers". This is a model for presenting genes that peak inside a segment or at a boundary.
- *Validation.* They report the reconstruction's sensitivity to its inputs and validate it independently: "Our reconstruction accuracy approached saturation at 6 genes"; "using smFISH on 20 genes". Mirror this with held-out genes and oracles.
- *Discussion (generality).* "can be readily applied to extract spatial blueprints of other structured mammalian organs". It supports our generality claim (crypts, lobules).

## [6] Moor, A. E. et al. *Cell* 175, 1156–1167.e15 (2018)

**Access. INACCESSIBLE — abstract only.** Routes tried:
- cell.com HTML and the PDF that Semantic Scholar lists as "BRONZE" open (cell.com/article/S0092867418311644/pdf): 403 by both curl and WebFetch.
- ScienceDirect: 403.
- bioRxiv preprint 10.1101/261529: 429 (Cloudflare rate limit) by curl, WebFetch and the JATS link from the bioRxiv API.
- Weizmann institutional (cris.iucc.ac.il) file: 403.
- Europe PMC: no PMC copy and no preprint full text.

I read only the Europe PMC abstracts of the Cell version and of the preprint.

**Citing sentence.** Same as [5].

**Verdict (from the abstract only; provisional): supports the villus part.**
- Cell abstract: "The intestinal epithelium is a highly structured tissue composed of repeating crypt-villus units."
- Cell abstract: "uncovering broad zonation of enterocyte function along the villus."
- Preprint abstract: "We found that most enterocyte genes were zonated."

Re-read the full text before submission; it is open on cell.com from a normal browser.

**Writing ideas (abstract only).**
- *Discussion (generality).* "Our approach can be used to uncover zonation patterns in other organs when prior knowledge of landmark genes is lacking." This is the closest phrasing to our label-free claim. Note that Moor still needed landmark genes from microdissection, which we do not.
- *Results framing.* Their sequence of functional zones ("sequential expression of carbohydrates, peptides, and fat absorption machineries in distinct villus compartments") is a template for describing transporter handoffs along the PT.

## [7] Curthoys, N. P. & Moe, O. W. *Clin. J. Am. Soc. Nephrol.* 9, 1627–1638 (2014)

**Access.** PMC HTML, PMC4152816. I read the full review body and figure legends.

**Citing sentence (Introduction).** "Along the nephron, the proximal tubule alone runs through three subsegments, S1, S2 and S3, which differ in solute transport, metabolism and susceptibility to injury[7–9]."

**Verdict: PARTLY SUPPORTED.** Transport and metabolism are supported. Injury is not mentioned anywhere in the review.
- Glucose Transport: "The SGLT2 transporter is localized to the brush border membrane of the S1 and S2 segments of the proximal tubule where it extracts the bulk of the filtered glucose."
- Glucose Transport: SGLT1 "is preferentially expressed in the apical membrane of the S3 segment."
- Response to Acidosis: "SNAT3 is normally localized solely to the basolateral membrane of the S3 segment"; in chronic acidosis it is found "primarily in the basolateral membranes of the S1 and S2 segments".

**Fix.** Use rewrite S3.

**Writing ideas.**
- *Results, transporter handoffs (species and AKI).* Physiological rationale for a handoff: "the sequential positioning of the two transporters that favor capacity and affinity in the early and late tubule, respectively".
- *AKI / Discussion.* SNAT3's shift from S3 to S1/S2 in acidosis is a precedent for a condition-dependent axial move of a transporter domain, which is exactly the kind of switch-point shift we report.
- *Introduction (axial heterogeneity).* Figure 1 legend: "Inorganic phosphate (Pi) reabsorption is more avid in the earlier parts of the proximal tubule"; Fig. 1A plots solutes "along the length of the proximal tubule". Classical physiology already treats the PT as a continuum.
- *Introduction (scale).* The Summary has a crisp quantitative opener: "reabsorbing approximately 60%–70% of the water and NaCl ... and nearly all of the nutrients".

## [8] Lee, J. W., Chou, C.-L. & Knepper, M. A. *J. Am. Soc. Nephrol.* 26, 2669–2677 (2015)

**Access.** PMC HTML, PMC4625681. I read the text and figure legends; Tables 1–2 did not render as text.

**Citing sentence.** Same as [7].

**Verdict: NOT SUPPORTED for this sentence.** The study is rat. It profiles S1, S2 and S3 as separate microdissected segments but does not describe how they differ in transport, metabolism or injury; its Results treat "the proximal tubule segments" as one group.
- Concise Methods: "Male Sprague–Dawley rats weighing 200–250 g".
- Critical Enzymes in Metabolic Pathways: "the proximal tubule segments stand out as having a unique profile."
- Figure 2 legend: "Replicates from 14 different renal tubule segments were found to form six distinct clusters according to their anatomic and functional organization."

**Fix.** Move [8] to the repeated-units sentence (S2). It supports that clause directly. Introduction: "The mammalian kidney is made up of thousands of individual nephron units that consist of a glomerulus ... followed by a long epithelial tubule".

**Writing ideas.**
- *Introduction (discrete status quo).* Microdissection is the classic discrete comparator: "comprehensively profile gene expression in each of 14 renal tubule segments" (Abstract). Use it to frame segment-mean analyses.
- *Results, validation against microdissection.* "Manual microdissection allowed a high degree of sample purity" (Abstract). This is why microdissection data are a fair external reference for our ordering.
- *Introduction (scale).* The clustering of the three PT segments into one of six clusters is a citable sign that the S1–S3 differences are graded rather than categorical.

## [9] Chen, L., Chou, C.-L. & Knepper, M. A. *J. Am. Soc. Nephrol.* 32, 897–912 (2021)

**Access.** PMC HTML, PMC8017530. I read the full text and figure legends.

**Citing sentences.**
- (Introduction) Same as [7].
- (Methods) "Mouse pairs were confirmed when microdissected mouse PT segments[9] showed the early member higher in S1 than S3 and the late member the reverse, and human pairs likewise against single-cell PT segment profiles of six donors in CZ CELLxGENE Census[ref]."

**Verdict.** Introduction: PARTLY SUPPORTED (transport and metabolism, not injury). Methods: SUPPORTED.
- Results: "initial segment of the PT [PTS1], proximal straight tubule in cortical-medullary rays [PTS2], and last segment of the proximal straight tubule in the outer stripe of outer medulla [PTS3]".
- Results: "Some transcripts are exclusively localized to the PTS1 ... linked to glucose reabsorption (Slc5a2 [SGLT2]), lactate reabsorption (Slc5a12 ...), vitamin D metabolism".
- Results / Resource Web Page: segment TPMs are available at "https://esbl.nhlbi.nih.gov/MRECA/Nephron/ or Supplemental Table 1". This is what the Methods check needs.

Note for Methods: Chen's segments are defined by location ("PTS2" means straight tubule in medullary rays). The mice were "Male, 6- to 8-week-old, C57BL/6", and "Four to eight tubules were collected for each sample". State these where the confirmation step is described.

**Writing ideas.**
- *Results (subtle PT axis).* "The three PT segments (PTS1, PTS2, and PTS3) exhibited comparable expression patterns, which differed markedly from those of other tubule segments". Use it to motivate why within-PT gradients need a dedicated method.
- *Validation (independence of the reference).* Discussion: "Specific segments were identified on the basis of morphology, topology, and location in the kidney, and (unlike scRNA-seq) small-sample RNA-seq does not rely on particular marker genes". This makes it a reference that is not circular with marker-based labels.
- *Species comparison.* Discussion: "the methodology described here can be extended ... between strains or species". They flag species comparison as future work.

## [10] Ransick, A. et al. *Dev. Cell* 51, 399–413.e7 (2019)

**Access.** PMC HTML (NIH author manuscript), PMC6948019. I read the Summary, Introduction, Results, Discussion and STAR Methods.

**Citing sentence (Introduction).** "Single-cell and spatial studies have mapped these segment states in mouse and human kidneys and in acute injury[10–15], but how gene expression changes along the length of the tubule, within and between its segments, remains largely uncharacterized."

**Verdict: SUPPORTED** (mouse, healthy).
- Results, Sex-related diversity: "genes known to highlight the segmental organization of the proximal tubule ...: Spp2/Slc5a2 (S1), Slc34a1 (S1–2) and Slc22a7 (S3)".
- Results: "The highest disparity in gene activity between the sexes was in the S3 region."
- Results: "mice and humans both show a male-biased susceptibility to ischemia invoked kidney injury where the S3 region is particularly vulnerable to injury". This can carry the injury clause of S3.

**Writing ideas.**
- *Introduction (repeated units, S2).* "A mouse kidney contains around 14,000 nephrons ..., the human kidney 900,000–1,000,000"; "Each nephron comprises multiple segments with distinct physiological activities".
- *Results / Discussion (boundaries).* "S2 and S3 sex-specific expression domains form a tight, non-overlapping boundary at the S2-S3 junction". This is a known sharp boundary our axis should reproduce; cite it when reporting a step-like gene at S2/S3.
- *Methods / Discussion, limitations (sex).* All our specimens are male; S3 is the most sex-divergent region. One clause belongs in Methods.
- *Introduction (why not dissociation).* "male S3 cells were markedly underrepresented in the scRNA-seq profile (~10% of expected". This is a dissociation bias that in-situ ST avoids.

## [11] Muto, Y. et al. *Nat. Commun.* 12, 2190 (2021)

**Access.** Europe PMC full-text XML, PMC8044133.

**Citing sentence.** Same as [10].

**Verdict: SUPPORTED** (human cortex, segment states plus an injury-like PT state).
- Results: "snRNA-seq and snATAC-seq were performed on five healthy adult kidney samples".
- Results: "snATAC-seq was able to detect two subpopulations within the proximal tubule cluster, which likely represent the proximal convoluted tubule ... and the proximal straight tubule".
- Results: PT_VCAM1 "may represent an injured or regenerative subpopulation of cells."

**Writing ideas.**
- *Introduction / species comparison (discrete calls are weak in human).* "The delineation between S1/S2 and S3 was less clear in the snRNA-seq dataset". This matches our human S1/S2 diagnosis.
- *Methods / limitations (human cohort).* "our study focused on kidney cortex and did not include samples from the medulla" (Discussion). This mirrors our human slices, which are both cortex.
- *Discussion, limitations.* "The small sample size of this study does not adequately capture the expected heterogeneity of the general population." This is a ready template for the one-donor caveat.

## [12] Lake, B. B. et al. *Nature* 619, 585–594 (2023)

**Access.** Europe PMC full-text XML, PMC10356613: main text, Methods and Extended Data legends.

**Citing sentence.** Same as [10].

**Verdict: SUPPORTED** (human; single-cell, Slide-seq2 and Visium; healthy, AKI and CKD).
- Constructing a kidney cellular atlas: "spatial transcriptomic Slide-seq2 ... (6 individuals, 67 pucks) and Visium assays (22 individuals, 23 samples)".
- Reference and altered states: "Putative adaptive (successful or maladaptive tubular repair) states were primarily found within the proximal tubule (PT) and TAL subclasses in mouse and human kidneys."
- Methods, Integrating snCv3 and scCv3: "Cell types that could not be accurately resolved (PT-S1/PT-S2) were kept merged."

**Writing ideas.**
- *Introduction / species comparison.* The PT-S1/PT-S2 merge sentence is strong, citable evidence that discrete human PT segment labels are unstable.
- *Methods (stating resolution per platform).* "Slide-seq and Visium transcriptomic data at different resolutions (10 µm and 55 µm beads, respectively)".
- *AKI.* Their altered-state taxonomy ("cycling, transitioning, adaptive (successful or maladaptive repair) and degenerative") gives names to use if our AKI section touches injured-PT states.

## [13] Kirita, Y. et al. *Proc. Natl Acad. Sci. USA* 117, 15874–15883 (2020)

**Access.** Europe PMC full-text XML, PMC7355049. I read the main text, Materials and Methods and figure legends; the SI Appendix was not read.

**Citing sentence.** Same as [10].

**Verdict: SUPPORTED** (mouse snRNA-seq, segment states in acute injury).
- Proximal Tubule Responses to Acute Injury: "yielded three healthy subclusters (the S1, S2, and S3 segments of the proximal tubule), one repairing subcluster, and three injured subclusters".
- Same section: "we annotated these as “injured S1/2,” “injured S3,” and “severe injured PT,” respectively."
- Materials and Methods: "8- to 10-wk-old male mice"; the renal pedicle was clamped "for 18 min".

**Writing ideas.**
- *AKI (our injury model).* Their time course (4 h, 12 h, 2 d, 14 d, 42 d; n = 3 per group) and clamp time are the reference for stating ours. Report the clamp time and the post-injury day the same way.
- *Introduction.* "a subpopulation of injured proximal tubule (PT) epithelia (the segment most susceptible to injury)". This is citable for PT vulnerability, but not S3-specific.
- *Introduction / Results (segment-specific transporters).* "Each nephron segment performs unique reabsorptive and secretory functions ... reflected by segment-specific expression of all detected solute-linked carriers". This is the discrete view we contrast with.

## [14] Marshall, J. L. et al. *iScience* 25, 104097 (2022)

**Access.** Europe PMC full-text XML, PMC8971939.

**Citing sentence.** Same as [10].

**Verdict: NOT SUPPORTED** for "mapped these segment states ... in acute injury".
- The mouse models are a diabetic kidney disease model and a toxic UMOD proteinopathy model. The human cohort is nine nephrectomies, one with ischaemic injury from tumour compression.
- The only PT label used is "PCT"; no S1/S2/S3 or injured-PT states are mapped.
- Summary: "in a mouse model of diabetic kidney disease, we detected changes in the cellular organization of the spatially restricted kidney filter".
- Results: "proximal convoluted tubules (PCTs) were enriched in cortex".
- Results: "one sample displayed clear evidence of injury (ischemia due to tumor compression".

**Replacement, full text read.** Melo Ferreira, R. et al. Integration of spatial and single-cell transcriptomics localizes epithelial cell-immune cross-talk in kidney injury. *JCI Insight* 6, e147703 (2021). doi:10.1172/jci.insight.147703, PMID 34003797, PMCID PMC8262485 (Europe PMC XML read).
- Visium on a human reference nephrectomy and on mouse IRI and CLP models at 6 h, with PT S1, S2 and S3 mapped.
- Results: "Both a cortical PT S3 cluster (PT-S3-C) and a medullary PT S3 cluster (PT-S3-OS) were identified".
- Abstract: "Atf3 was identified as a chemotactic factor in S3 proximal tubules."
- Caveat: Fig. 4 legend says "n = 1 murine sample per model".

**Writing ideas (Marshall, if kept elsewhere).**
- *Methods / Introduction (resolution).* "a feature (or bead) size of 10 μm most commonly capturing one to two cells/feature and no more than three cells/feature". This is a model for saying what a feature contains.
- *Discussion, limitations.* "Slide-seqV2 data have limited power for the detection of genes with low expression levels at near-single-cell resolution". The same applies to 2-µm bins before aggregation.
- *Results (aggregating to structures).* They filtered PCT, DCT and CD beads by k-nearest neighbours "selecting parameters that maintained beads contained within tubule structures". This is a precedent for structure-level aggregation.

## [15] Dixon, E. E. et al. *J. Am. Soc. Nephrol.* 33, 279–289 (2022)

**Access.** PMC HTML, PMC8819997.

**Citing sentence.** Same as [10].

**Verdict: SUPPORTED** (Visium, female mouse AKI time course, PT segment and injured states).
- Results: "a bilateral renal pedicle clamp time of 34 minutes in female mice induced an equivalent degree of injury as 22 minutes in male mice."
- Figure 4 legend: the cell types include "injured proximal tubule cells (InjPT), failed repair proximal tubule cells (FR-PTC) ... proximal tubule segments 1–3 (PTS1–3)".
- Results (Giotto section): Cryab "in the S3 segment ... was upregulated with injury and downregulated with repair".

**Writing ideas.**
- *AKI (regional responses).* "region-specific injury and repair transcriptional responses" (Abstract, Conclusions). This is the closest prior claim to our within-segment AKI shifts; we resolve it along the axis instead of by region.
- *Discussion, limitations (sections as replicates).* "we anticipate challenges in normalizing biologic replicate variability attributed to tissue slice variation". This applies directly to our section-level design.
- *Introduction / Methods (resolution).* "resolution is limited to spots encompassing 10–50 cells due to the distance between oligo-embedded spots". This is the contrast with Visium HD.
- *AKI (sex).* "The vast majority of murine AKI studies to date have used male mice because they are more susceptible to ischemia-reperfusion injury". Our male AKI mice are the standard choice, but say so.

## [16] Zhai, X. Y. et al. *J. Am. Soc. Nephrol.* 14, 611–619 (2003)

**Access.** Publisher full-text HTML (JASN on Ovid/LWW, free access): Abstract, Materials and Methods, Results, Discussion and figure legends. Tables are images and were not read.

**Citing sentence (Introduction).** "Recovering that position physically requires tracing individual nephrons through serial sections[16,17] or profiling tissue in three dimensions[18–21], approaches that remain laborious and are rarely paired with transcriptome-wide measurements across whole organs."

**Verdict: SUPPORTED** (serial-section tracing). "Laborious" is implied, not stated: no time or cost is given.
- Materials and Methods: "Serial sections, 2.5-μm-thick, were generated from the surface of the kidneys to the inner stripe ... about 800 sections from each animal".
- Manual Tracing: tubules "were traced through about 800 hard copies made from the digitized images."
- Computer-Assisted Digital Tracing: "If a tubular path exited the boundary of the matrix, it was discarded."

**Writing ideas.**
- *Introduction / Discussion (why continuous).* The mouse PT has no ultrastructural segments: "Ultrastructural analyses reveal no obvious morphologic segmentation of the proximal tubule". The caveat "do not exclude a possible segmentation based on variations in enzyme/transporter-protein expression" is exactly what our axis tests.
- *Introduction (cost, quantified).* State the burden in numbers: about 800 sections per kidney for 160 PTs.
- *Methods (smoothing precedent).* They smoothed traced paths with a moving average to remove "a zig-zag appearance never observed in longitudinal sectioned proximal tubules". This is a precedent for smoothing ordered positions.
- *Validation.* Computer tracing was validated against manual tracing: "no significant difference between the two methods". This is a model for an oracle comparison.

## [17] Zhai, X. Y. et al. *J. Am. Soc. Nephrol.* 17, 77–88 (2006)

**Access.** Publisher full-text HTML (JASN on Ovid/LWW): Abstract, Materials and Methods, Results, Discussion and figure legends. Table 1 is an image and was not read.

**Citing sentence.** Same as [16].

**Verdict: SUPPORTED.** No time or cost is stated.
- Preparation of Renal Tissue: "From each kidney, a total of 2500, 2.5-μm-thick consecutive sections were obtained from the surface to the papillary tip".
- Digital Tracing: "A marker was placed manually in the tubule lumen along the course of the tubule."
- Results: "In total, 30 CD, 127 SLN, and 41 LLN were traced. Of these, 24 LLN were complete, whereas 17 LLN were incomplete".

**Writing ideas.**
- *Introduction (cost, quantified).* "2,500 sections per kidney, every second one digitized, markers placed by hand, 200 nephrons" is a concrete contrast to one section per specimen.
- *Results (segment lengths as a prior).* PT, thin limb, TAL and DCT lengths per nephron type give an anatomical yardstick if we report relative positions along the nephron. Relative segment lengths "are constant for all SLN".
- *Discussion (gradual transitions).* "The epithelia along DTL of the LLN transformed gradually from one type to another". This is morphological support that segments grade into each other.
- *Discussion, limitations wording.* "it should be emphasized that the findings of this study are based on C57/BL/6J mice only." This is a template for stating strain scope.

## [18] Wang, X. et al. *Science* 361, eaat5691 (2018)

**Access.** PMC HTML (author manuscript), PMC6339868: Abstract, Results, Discussion, Methods summary and figure legends. The Supplementary Materials were not read.

**Citing sentence.** Same as [16].

**Verdict: SUPPORTED** (3D in situ profiling with a targeted panel in thick tissue).
- Adapting STARmap to thick tissue blocks: "Using a curated gene set including 23 cell-type markers and five ARGs read out over seven cycles".
- Same section: "more than 30,000 cells across volumes spanning all six layers and the corpus callosum."
- Same section: "STARmap successfully detected yellow fluorescent protein (YFP) mRNA across 150 µm of tissue thickness". This is a validation sample, not the main volume.

**Writing ideas.**
- *Introduction (3D needs).* "a short-range self-clustering organization ... that could only be accurately measured in three but not in two dimensions". Use it for the value of 3D, then note our structural axis is recovered from 2D.
- *Introduction (trade-off).* Thin sections allow 160–1,020 genes; thick volumes used 28. This is the gene-number versus volume trade-off that motivates a computational route.
- *Introduction.* "methods of 3D spatial transcriptomic analysis in thick tissues have long been sought". This phrasing establishes a long-standing need.

## [19] Vickovic, S. et al. *Commun. Biol.* 5, 129 (2022)

**Access.** Europe PMC full-text XML, PMC8837632.

**Citing sentence.** Same as [16].

**Verdict: SUPPORTED** (3D, transcriptome-wide, one array per consecutive section, biopsies only).
- 3D spatial profiling of RA synovia: "we profiled 27 tissue sections and a total of 17,117 tissue spots by ST".
- Methods, Spatial Transcriptomics reactions: "Tissues were cryosectioned at 7 µm thickness. Each section was ... placed onto an individual array".
- Methods, Image registration: "As the spatially resolved data are of the restricted resolution, the data were interpolated using the akima package".

**Writing ideas.**
- *Introduction (cost per z-level).* Every z-level costs one array, and six biopsies took 27 arrays. This is the per-specimen burden our single-section approach avoids.
- *Introduction / Methods (feature density).* "Tissue sections are placed on a glass slide covered with 1000–2000 features". This is a model for stating the platform's limit next to Visium HD's.
- *Methods (study design wording).* They describe a prospective power analysis: "sampling at least 3 patients in 4 replicate sections ... would provide sufficient power". This is a useful model if we justify sections versus specimens.

## [20] Schott, M. et al. *Cell* 187, 3953–3972.e26 (2024)

**Access.** Published version (open access, CC BY-NC), PDF from the MDC institutional repository (edoc.mdc-berlin.de, eprint 24484). I read the full article including STAR Methods. The PDF text renders µm as "mm" and € as "V"; the quotes below are corrected.

**Citing sentence.** Same as [16].

**Verdict: PARTLY SUPPORTED.** It is 3D, whole-transcriptome, built from serial sections, which fits "profiling tissue in three dimensions" and "rarely ... across whole organs" (one 350-µm lymph-node block). But it presents itself as cheap and scalable, which conflicts with "laborious".
- Open-ST generates 3D virtual tissue blocks: "we obtained 10 µm sections spanning 350 µm of tissue depth. Our experimental design consisted of profiling sections with Open-ST (19), H&E staining only (11)".
- Introduction: "a total cost of <€130 per 12 mm² capture area"; "a single researcher can prepare 10–15 libraries in 3 days".
- Limitations: "The discrepancy in z- versus x-y resolution (10 µm vs. 0.6 µm) presents a complementary limitation".

**Fix.** Use rewrite S5, which drops "laborious" in favour of "many sections or assays per specimen ... sub-millimetre depths".

**Writing ideas.**
- *Introduction (gap statement, 2024).* "no end-to-end platform currently exists to generate and computationally analyze ST in 3D" (Introduction).
- *Introduction / Discussion (what 3D buys).* Summary: "the 3D reconstruction ... revealed spatially contiguous structures not visible in 2D". We claim a structural axis without 3D data, so we should say what we cannot recover, such as neighbourhoods along z.
- *Discussion, limitations (cost framing).* Their honest limitation is a model: "requires an initial financial investment ... and is therefore not well-suited for low-throughput experiments".

## [21] Kiemen, A. L. et al. *Nat. Methods* 19, 1490–1499 (2022)

**Access.** PMC HTML (author manuscript), PMC10500590: Abstract, Introduction, Results, Discussion, Online Methods and Extended Data legends.

**Citing sentence.** Same as [16].

**Verdict: PARTLY SUPPORTED.** CODA reconstructs 3D tissue from serial H&E with deep-learning labels and measures no gene expression. It supports "laborious" (thousands of sections) and "rarely paired with transcriptome-wide measurements". It belongs with serial-section reconstruction rather than molecular "profiling in three dimensions".
- Introduction: "We analyzed 4,114 H&E sections to reconstruct thirteen samples of up to 3.5cm3".
- Discussion: "The sectioning of tissue causes unpredictable warping, requiring sophisticated registration techniques."
- Discussion: "future addition of IHC labeling, spatial ‘omics’, and gene expression imaging to the intervening sections will increase the number of labels".

**Fix.** Use rewrite S5, which lists CODA as the morphology-only example.

**Writing ideas.**
- *Introduction (cost of serial sectioning, quantified).* "time-consuming manual annotations and costly immunohistochemical (IHC) labeling and mass spectrometry have been required" (Introduction).
- *Methods / validation (sub-sampling test).* They quantify what is lost when sections are skipped: "<5% error in 3D cell count and tissue composition extrapolation when skipping up to two serial sections". This is a model for our held-out and subsampling checks.
- *Methods (H&E deep-learning segmentation).* CODA is a precedent for multi-class H&E semantic segmentation feeding a structural analysis, and it states timing gains: "on average 3-fold faster than the other techniques".

---

## Notes on earlier entries in intro_research.md (section 3)

- Entries reused only where labelled FULL TEXT, and every quote above was re-read from the downloaded text: Oliveira, Lake, Muto, Kirita, Marshall, Chen, Lee, Ransick, Curthoys, Vickovic, Wang, Kiemen, Schott, Melo Ferreira.
- Entries labelled ABSTRACT in the notes (Ståhl, Moses, Zhai 2003, Zhai 2006, Dixon) are now read in full here, except Moses (preprint only).
- Marx stays at standfirst only; Moor stays at abstract only.
- Correction to the earlier Marshall note: it is not a mouse AKI study, and its PT annotation is PCT only.
- Correction to the earlier Curthoys note ("adequate as the PT-segment review"): it carries transport and metabolism, not injury.
