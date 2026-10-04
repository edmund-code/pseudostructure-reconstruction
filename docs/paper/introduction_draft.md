# Introduction (draft, workstream 4)

Status: draft written 2026-10-04. Every factual statement below cites a source whose title and
abstract were checked in Europe PMC during this session. Markers: **[CITE]** = a claim not yet
sourced; **[WS1]**, **[WS2]**, **[WS3]** = depends on that workstream's final result. Numbers from
our own data are deliberately left out until workstreams 1–3 freeze them.

---

**The proximal tubule matters for translation.** The proximal tubule (PT) does much of the
kidney's work:

- It reabsorbs about 60–70% of filtered water and NaCl, more of the filtered bicarbonate, and nearly
  all filtered nutrients. It also hosts active solute secretion, hormone production and most of the
  kidney's metabolic functions [1].
- During fasting and stress, the kidney supplies up to 40% of systemic glucose production through
  gluconeogenesis. This capacity is impaired in PT cells during acute kidney injury (AKI) [2].
  Defective fatty-acid oxidation in tubular epithelium contributes to fibrosis [3].
- The organic anion transporters of the SLC22 family handle common drugs, toxins, nutrients and
  putative uremic toxins [4].

The PT is regarded as the primary target of injury in AKI and chronic kidney disease [5], and as
the main site of drug-induced kidney injury [6]. Much of what is known about these processes comes
from rodents. Yet no intervention has been shown to prevent AKI or speed recovery in patients,
which has raised the question of whether animal models predict human responses [7]. Rodent AKI
models do not fully recapitulate the human disease [8], and cross-species differences in renal drug
transport and metabolism contribute to poor translation of preclinical pharmacokinetic and
toxicity data [9,10]. Knowing where human and mouse PT differ is therefore practically important
for drug development and for interpreting mouse injury models.

**The PT is not uniform, and its zonation differs between species.** The PT is divided along its
length into S1, S2 and S3 segments, and functions are distributed along this axis:

- In mouse, the glucose transporter SGLT2 is found in the early PT [11].
- In isolated human PT segments, S2 and S3 make more glucose from lactate than S1 [12].
- After ischaemia in mice, tubular damage is most severe in S3 of the outer stripe of the outer
  medulla [13].

Segment-resolved transcriptomes of microdissected tubules exist for rat [14] and mouse [15].
Single-cell atlases now describe PT diversity in mouse [16] and human [17] kidney.

Several species differences map onto this axis. In human kidney, the organic anion transporters
OAT1–3 stain more strongly in S1/S2 than in S3 within the medullary rays, and S3 in the outer
stripe is unstained [18]. In mouse, sexually dimorphic gene activity maps predominantly to PT
segments, whereas only a limited set of genes shows conserved sex-linked regulation in human
kidney [19]. Sex-biased spatial patterns in mouse extend to the cortex and the outer stripe [20].

These observations come from different studies, assays and segment definitions. To our knowledge,
no study has compared human and mouse PT gene programs position by position along the S1→S3 axis
[WS3: confirm against the final literature check].

**Comparing matched cell types or segments cannot show where along the PT the species differ.**
Cross-species single-cell atlases align conserved cell states between human and rodent kidneys
[21]. This answers which cell types correspond, not where along a segment a program changes.
Dissociation also discards the spatial relationships between cells [22]. Assigning cells or
structures to S1, S2 and S3 turns the PT into three bins, so any difference within a bin, or a
difference in where an expression boundary lies, is averaged into a constant offset per segment
[CITE: evidence that PT gene expression varies continuously within and across S1–S3, rather than in
discrete steps].

Continuous reconstructions in other organs show what such averaging can hide:

- In mouse liver, about half of all genes are zonated along the lobule axis, and many of the
  profiles are non-monotonic, peaking in the middle layers rather than at either end [23].
- Enterocyte functions are broadly zonated along the intestinal villus [24].

**Spatial transcriptomics and pseudospace reconstruction.** Spatial transcriptomics keeps tissue
context. In kidney it has been used to map cell types and injury niches in human and mouse,
combined with single-nucleus data [25] or using Slide-seqV2 [26], and in the human kidney atlas
[17]. Visium HD measures the
whole transcriptome at single-cell-scale resolution on tissue sections [1a].

A tissue section, however, cuts each tubule at an arbitrary point, so a section gives no direct
reading of where a structure lies along the nephron. That position has to be reconstructed. There
are several precedents:

- Landmark genes were used to infer lobule coordinates in liver [23] and villus coordinates in
  intestine [24].
- Spatial arrangements have been recovered from expression alone [22].
- Recent methods unroll tissue architecture from spatial data to expose zonated programs [27].
- Trajectory inference orders cells along continuous processes, with dozens of methods
  benchmarked [28]. These include diffusion pseudotime [29] and principal-graph methods such as
  scFates [30].

Applied to segmented PT cross-sections, the same idea yields a continuous "pseudospace" coordinate
along which human and mouse expression can be compared.

**Few specimens make species comparisons hard to separate from specimen differences.**
Cells or structures from one specimen are not independent replicates. Methods that ignore
variation between biological replicates are biased, and can report hundreds of differentially
expressed genes where no biological difference exists [31,32]. Aggregating counts to one
pseudobulk per specimen is conservative but underpowered relative to mixed models [32]. Both
approaches have been benchmarked for multi-sample designs [33]. Lamian extended sample-aware testing to pseudotime. When its authors
randomly partitioned samples into two groups with no expected difference, methods that ignore
sample-level variation reported thousands of false-positive differential genes [34].

Pathway analysis has a related problem. Competitive gene-set tests assume that genes are
independent. Correlation among genes in a set inflates their false-positive rate [35,36].

Human tissue is scarce, and spatial studies often have only a few sections per species. In that
setting, an average species difference estimated from two specimens per species is confounded with
specimen-level variation. A useful analysis must therefore show that its calls are specific to the
species contrast, not merely that they are numerous.

**This study.** We profiled two control mouse kidney sections and two cortex sections from one
human donor with Visium HD. We segmented individual tubule cross-sections from the paired H&E
images with a panoptic segmentation model, and aggregated 2-µm bins into tubule-level profiles in
a shared one-to-one ortholog space. After integrating the specimens and reviewing cluster labels,
we reconstructed a continuous PT coordinate across species [WS1: final method and its
justification].

Along this coordinate we asked where, not whether, human and mouse PT programs differ. The test
statistic is a gene-level species-by-position term from nested regression models; pathways are
tested with a covariate-matched competitive rank test. To check specificity, we reran every
strategy on balanced relabelings of the four specimens, which contain no species difference.
Average-level pathway screens reported as many pathways for these relabelings as for species,
whereas the smooth position-dependent screen was specific [WS2: final numbers]. We then tested the
position-dependent differences against coordinate misregistration, the cortex-only human
sampling, genes detected in one species only, mouse sex bias and the choice of trajectory method.
Finally, we compared them with independent segment-resolved data from microdissected tubules and
single-nucleus atlases [WS2/WS3: final replication numbers and the biology to highlight].

---

## References (all checked in Europe PMC on 2026-10-04)

1. Curthoys NP, Moe OW. Proximal tubule function and response to acidosis. *Clin J Am Soc Nephrol* 9:1627–1638 (2014). https://pubmed.ncbi.nlm.nih.gov/23908456/

   1a. Oliveira MF et al. High-definition spatial transcriptomic profiling of immune cell populations in colorectal cancer. *Nat Genet* 57:1512–1523 (2025). https://pubmed.ncbi.nlm.nih.gov/40473992/ (introduces Visium HD; renumber when merging.)
2. Legouis D et al. Altered proximal tubular cell glucose metabolism during acute kidney injury is associated with mortality. *Nat Metab* 2:732–743 (2020). https://pubmed.ncbi.nlm.nih.gov/32694833/
3. Kang HM et al. Defective fatty acid oxidation in renal tubular epithelial cells has a key role in kidney fibrosis development. *Nat Med* 21:37–46 (2015). https://pubmed.ncbi.nlm.nih.gov/25419705/
4. Nigam SK et al. The organic anion transporter (OAT) family: a systems biology perspective. *Physiol Rev* 95:83–123 (2015). https://pubmed.ncbi.nlm.nih.gov/25540139/
5. Chevalier RL. The proximal tubule is the primary target of injury and progression of kidney disease: role of the glomerulotubular junction. *Am J Physiol Renal Physiol* 311:F145–F161 (2016). https://pubmed.ncbi.nlm.nih.gov/27194714/
6. Zuidervaart SA et al. Rat and human primary proximal tubule epithelial cells allow for cross-species comparisons of drug-induced nephrotoxicity in vitro. *Arch Toxicol* 100:2015–2027 (2026). https://pubmed.ncbi.nlm.nih.gov/41699311/ (Abstract states that PT is the primary site of drug-induced injury and that rat-to-human translation is hampered by species differences. A review would be a stronger citation for the first point [CITE].)
7. de Caestecker M et al. Bridging Translation by Improving Preclinical Study Design in AKI. *J Am Soc Nephrol* 26:2905–2916 (2015). https://pubmed.ncbi.nlm.nih.gov/26538634/
8. Fu Y et al. Rodent models of AKI and AKI-CKD transition: an update in 2024. *Am J Physiol Renal Physiol* 326:F563–F583 (2024). https://pubmed.ncbi.nlm.nih.gov/38299215/
9. Thakur A et al. Sex and the Kidney Drug-Metabolizing Enzymes and Transporters: Are Preclinical Drug Disposition Data Translatable to Humans? *Clin Pharmacol Ther* 116:235–246 (2024). https://pubmed.ncbi.nlm.nih.gov/38711199/
10. Basit A et al. Kidney Cortical Transporter Expression across Species Using Quantitative Proteomics. *Drug Metab Dispos* 47:802–808 (2019). https://pubmed.ncbi.nlm.nih.gov/31123036/
11. Vallon V et al. SGLT2 mediates glucose reabsorption in the early proximal tubule. *J Am Soc Nephrol* 22:104–112 (2011). https://pubmed.ncbi.nlm.nih.gov/20616166/
12. Conjard A et al. Gluconeogenesis from glutamine and lactate in the isolated human renal proximal tubule: longitudinal heterogeneity and lack of response to adrenaline. *Biochem J* 360:371–377 (2001). https://pubmed.ncbi.nlm.nih.gov/11716765/
13. Kim J et al. Intra-renal slow cell-cycle cells contribute to the restoration of kidney tubules injured by ischemia/reperfusion. *Anat Cell Biol* 44:186–193 (2011). https://pubmed.ncbi.nlm.nih.gov/22025970/ (A classic reference on S3 vulnerability to ischaemia would be preferable [CITE].)
14. Lee JW et al. Deep Sequencing in Microdissected Renal Tubules Identifies Nephron Segment-Specific Transcriptomes. *J Am Soc Nephrol* 26:2669–2677 (2015). https://pubmed.ncbi.nlm.nih.gov/25817355/
15. Chen L et al. A Comprehensive Map of mRNAs and Their Isoforms across All 14 Renal Tubule Segments of Mouse. *J Am Soc Nephrol* 32:897–912 (2021). https://pubmed.ncbi.nlm.nih.gov/33769951/
16. Ransick A et al. Single-Cell Profiling Reveals Sex, Lineage, and Regional Diversity in the Mouse Kidney. *Dev Cell* 51:399–413 (2019). https://pubmed.ncbi.nlm.nih.gov/31689386/
17. Lake BB et al. An atlas of healthy and injured cell states and niches in the human kidney. *Nature* 619:585–594 (2023). https://pubmed.ncbi.nlm.nih.gov/37468583/
18. Breljak D et al. Distribution of organic anion transporters NaDC3 and OAT1-3 along the human nephron. *Am J Physiol Renal Physiol* 311:F227–F238 (2016). https://pubmed.ncbi.nlm.nih.gov/27053689/
19. Xiong L et al. Direct androgen receptor control of sexually dimorphic gene expression in the mammalian kidney. *Dev Cell* 58:2338–2358 (2023). https://pubmed.ncbi.nlm.nih.gov/37673062/
20. Chen S et al. Multi-omic and spatial analysis of mouse kidneys highlights sex-specific differences in gene regulation across the lifespan. *Nat Genet* 57:1213–1227 (2025). https://pubmed.ncbi.nlm.nih.gov/40259083/
21. Klötzer KA et al. Analysis of individual patient pathway coordination in a cross-species single-cell kidney atlas. *Nat Genet* 57:1922–1934 (2025). https://pubmed.ncbi.nlm.nih.gov/40775269/
22. Nitzan M et al. Gene expression cartography. *Nature* 576:132–137 (2019). https://pubmed.ncbi.nlm.nih.gov/31748748/
23. Halpern KB et al. Single-cell spatial reconstruction reveals global division of labour in the mammalian liver. *Nature* 542:352–356 (2017). https://pubmed.ncbi.nlm.nih.gov/28166538/
24. Moor AE et al. Spatial Reconstruction of Single Enterocytes Uncovers Broad Zonation along the Intestinal Villus Axis. *Cell* 175:1156–1167 (2018). https://pubmed.ncbi.nlm.nih.gov/30270040/
25. Melo Ferreira R et al. Integration of spatial and single-cell transcriptomics localizes epithelial cell-immune cross-talk in kidney injury. *JCI Insight* 6:e147703 (2021). https://pubmed.ncbi.nlm.nih.gov/34003797/
26. Marshall JL et al. High-resolution Slide-seqV2 spatial transcriptomics enables discovery of disease-specific cell neighborhoods and pathways. *iScience* 25:104097 (2022). https://pubmed.ncbi.nlm.nih.gov/35372810/
27. Guo J et al. SMURF: soft-segmentation for single-cell reconstruction and topological analysis of spatial transcriptomic data. *Nat Commun* 17:7990 (2026). https://pubmed.ncbi.nlm.nih.gov/42350392/
28. Saelens W et al. A comparison of single-cell trajectory inference methods. *Nat Biotechnol* 37:547–554 (2019). https://pubmed.ncbi.nlm.nih.gov/30936559/
29. Haghverdi L et al. Diffusion pseudotime robustly reconstructs lineage branching. *Nat Methods* 13:845–848 (2016). https://pubmed.ncbi.nlm.nih.gov/27571553/
30. Faure L et al. scFates: a scalable python package for advanced pseudotime and bifurcation analysis from single-cell data. *Bioinformatics* 39:btac746 (2023). https://pubmed.ncbi.nlm.nih.gov/36394263/
31. Squair JW et al. Confronting false discoveries in single-cell differential expression. *Nat Commun* 12:5692 (2021). https://pubmed.ncbi.nlm.nih.gov/34584091/
32. Zimmerman KD et al. A practical solution to pseudoreplication bias in single-cell studies. *Nat Commun* 12:738 (2021). https://pubmed.ncbi.nlm.nih.gov/33531494/
33. Crowell HL et al. muscat detects subpopulation-specific state transitions from multi-sample multi-condition single-cell transcriptomics data. *Nat Commun* 11:6077 (2020). https://pubmed.ncbi.nlm.nih.gov/33257685/
34. Hou W et al. A statistical framework for differential pseudotime analysis with multiple single-cell RNA-seq samples. *Nat Commun* 14:7286 (2023). https://pubmed.ncbi.nlm.nih.gov/37949861/
35. Goeman JJ, Bühlmann P. Analyzing gene expression data in terms of gene sets: methodological issues. *Bioinformatics* 23:980–987 (2007). https://pubmed.ncbi.nlm.nih.gov/17303618/
36. Wu D, Smyth GK. Camera: a competitive gene set test accounting for inter-gene correlation. *Nucleic Acids Res* 40:e133 (2012). https://pubmed.ncbi.nlm.nih.gov/22638577/

## Notes for the coordinator on claims and sources

- Statements about Squair, Zimmerman and Crowell are paraphrased from their abstracts. The claim
  that pseudobulk is "conservative and underpowered relative to mixed models" is Zimmerman's. The
  Lamian random-partition null was checked in the open-access full text (PMC10638410).
- I did not use CYP2E1 as an example of a species difference. The literature conflicts: one 2023
  report localises CYP2E1 to human PT (PMID 38042273), whereas an older study describes
  male-specific renal CYP2E1 activity in mouse (PMID 7839370). Workstream 3 should decide.
- Remaining [CITE] items: (a) continuous rather than stepwise PT zonation; (b) a classic review for
  S3 ischaemic vulnerability; (c) a review stating that PT is the main site of drug-induced injury.
