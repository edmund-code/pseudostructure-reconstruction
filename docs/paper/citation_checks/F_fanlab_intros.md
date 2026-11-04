# How the Fan lab (JEFworks) writes Introductions

Compiled 2026-11-03. These are notes on writing style. They are not draft prose for our paper.

## How this was checked

- Each paper was found with a Europe PMC REST search on its title, then fetched by PMCID. Every download went into its own directory under `scratchpad/fanlab/<paper>/`.
- Text was extracted with `python3 -I` scripts kept in `scratchpad/fanlab_scripts/`:
  - `jats_intro.py` reads JATS XML.
  - `pmc_html_intro.py` reads PMC HTML.
- Word counts strip citation markers, both numbered and author–year. They exclude figure captions and include numbers and URLs. The HTML-derived count for VeloViz may be off by about 2 words, because two author surnames in split citation links survived the stripping.
- Every quotation below was checked by script as an exact substring of the extracted text. Superscript citations are kept as [n]. Author–year citations are shown as printed.
- No personal identifiers were sent anywhere except as part of title or author search queries.

| # | Paper | Source read | Full text? |
|---|---|---|---|
| 1 | PAGODA, Fan 2016 *Nat Methods* 13:241–244 | Europe PMC XML of the NIH author manuscript (NIHMS746594, PMC4772672) | FULL TEXT, author-manuscript version |
| 2 | HoneyBADGER, Fan 2018 *Genome Res* 28:1217 | Europe PMC XML (PMC6071640), published version | FULL TEXT |
| 3 | MERINGUE, Miller 2021 *Genome Res* 31:1843 | Europe PMC XML (PMC8494224), published version | FULL TEXT |
| 4 | STdeconvolve, Miller 2022 *Nat Commun* 13:2339 | Europe PMC XML (PMC9055051), published version | FULL TEXT |
| 5 | VeloViz, Atta 2022 *Bioinformatics* 38:391 | PMC HTML of the published version (PMC8723140). See note | FULL TEXT |
| 6 | SEraster, Aihara 2024 *Bioinformatics* btae412 | Europe PMC XML (PMC11226864), published version | FULL TEXT |
| 7 | STcompare, Clifton 2026 *Bioinformatics* btag644 | Europe PMC XML (PMC13577157), published version | FULL TEXT |
| 8a | STARIT, Velazquez … Fan 2026 *PLoS Comput Biol* (extra) | Europe PMC XML (PMC13657327), published version | FULL TEXT |
| 8b | Gene count normalization in im-SRT, Atta … Fan 2024 *Genome Biol* 25:153 (extra) | Europe PMC XML (PMC11167774), published version | FULL TEXT |

Notes on sources:

- **VeloViz.** Europe PMC `fullTextXML` returned HTTP 500, and so did the preprint (PPR275603). NCBI efetch returned front matter only, with the comment "The publisher of this article does not allow downloading of the full text in XML form". The NCBI BioC OA API had no record. The PMC HTML article page carries the full published text.
- **STcompare.** I read the published Bioinformatics version. The bioRxiv preprint (PPR1125133) was not readable through Europe PMC (HTTP 500), so I did not compare the two versions.
- **PAGODA.** The pages are 241–244, and the main text of the author manuscript is about 2,100 words (Introduction 477, Results 1,502, Discussion 122), which is the size of a Brief Communication. The XML labels it "Article".
- **Choice of extras.** I picked the two 2023–2026 Fan last-author papers that frame a new representation or a methodological problem. I did not read:
  - the 2026 eLife Xenium off-target probe paper;
  - the 2025 scatterbar application note;
  - the 2026 Genome Biology cold-ischemia kidney paper, which is applied work.

---

## 1. PAGODA. Fan et al. 2016, *Nature Methods*

- **Source:** NIH author manuscript, Europe PMC XML. FULL TEXT.
- **Size:** 4 paragraphs, 477 words (57 / 164 / 163 / 93).

**Paragraph roles**

1. Context. scRNA-seq is unbiased, but analysis is hard because only some cell-to-cell differences matter.
2. Limits of PCA, GP-LVM, tSNE and clustering. Technical noise restricts them to highly expressed genes. Cross-cutting signals such as the cell cycle mask the processes of interest. Partitioning yields a single structure.
3. Tool, its basis, a rationale borrowed from GSEA, a worked example (Neurod1), and a one-sentence result on published data.
4. The biological system: what may drive heterogeneity in mouse NPCs. New data were generated, and a one-sentence result follows.

**How prior methods are criticised**

- Criticism is mostly by category: PCA and "other dimensionality reduction approaches", "commonly-used clustering approaches", "Partitioning methods".
- Named tools (GP-LVM, tSNE, BackSPIN, k-means) appear only as members of a category.
- The tone is mild and mechanistic: methods fail because of what the data contain, not because they are badly built.
- Limitation verbs: "pose difficulties for", "is often restricted to", "can mask", "posing a challenge for", "may, for example, choose to classify cells first based on the cell cycle phase", "difficult to discern by the existing heterogeneity analysis approaches".

**How the gap is framed**

- Mainly conceptual: cells admit several valid, cross-cutting classifications, while clustering seeks one partition.
- Secondarily a data limit: technical noise.

**How the tool is introduced**

> "Here, we describe an alternative approach for analyzing transcriptional heterogeneity called PAGODA that aims to detect all statistically-significant ways in which measured cells can be classified."

What follows it:

- one sentence on the mechanism (coordinated variability of annotated and de novo gene sets);
- an analogy to gene set testing in differential expression;
- a concrete example;
- a results sentence on published datasets;
- a separate paragraph on new NPC data and its result.

There is no URL and no generality claim beyond "previously published datasets".

**Example, figure callout, analogy**

- Running example: yes. Cell cycle or metabolic state masking tissue signals, and Neurod1 versus a coordinated neuronal-differentiation set.
- Figure callout: no.
- Analogy: yes, the GSEA rationale ("A similar rationale can be applied…").

**Discussion, first paragraph (the Discussion has one paragraph):** it opens with an analogy, "Just like organisms as a whole, individual cells can be classified according to a variety of meaningful criteria." It then states that an optimal single partition "is unlikely to be fully informative".

**Characteristic quotes**

- "Nevertheless, analysis of scRNA-seq data remains challenging, as measurements expose numerous differences between cells, only some of which may be relevant for system-level functions."
- "Such cross-cutting transcriptional features represent alternative ways to classify cells, posing a challenge for the commonly-used clustering approaches that aim to reconstruct a single subpopulation structure[5,8,9,11]."
- "A similar rationale can be applied in the context of heterogeneity analysis."
- "Examining previously published datasets, we illustrate that PAGODA recovers known subpopulations and reveals additional subsets of cells in addition to providing important insights about the relationships amongst the detected subsets."

## 2. HoneyBADGER. Fan et al. 2018, *Genome Research*

- **Source:** Europe PMC XML, published version. FULL TEXT.
- **Size:** 3 paragraphs with no heading (they come before "Results"), 367 words (171 / 93 / 103).

**Paragraph roles**

1. Biological and clinical motivation. Intratumor heterogeneity matters. How genetic alterations act on transcription is "not well understood". Linking subclones to transcriptional states "would be helpful", and the paragraph lists clinical payoffs.
2. Technology limit. Bulk data are insufficient, single cells are needed, and joint DNA and RNA profiling of one cell "remains challenging".
3. Existing computational routes:
   - bulk CNV callers that use coverage and allelic imbalance;
   - scRNA-seq expression-deviation CNV calls, which detect only "chromosome-scale CNVs".

   The paragraph then gives the allele idea and the tool sentence.

**How prior methods are criticised**

- Prior methods are cited by author (Patel 2014, Tirosh 2016) rather than by tool name.
- The limitation is implied by scope, not stated: they illustrate "chromosome-scale CNVs", while the tool targets "focal CNV and … LOH". There is no explicit verb of failure.
- Tone: very mild.

**How the gap is framed**

- A data and technology limit: no simultaneous DNA and RNA from a single cell.
- A resolution limit of current inference.

**How the tool is introduced**

> "Here, we propose a computational approach called HoneyBADGER to quantitatively infer the presence of subclone-specific focal CNV and loss-of-heterozygosity (LOH) events in individual cells using allele and expression information from scRNA-seq data."

It is the last sentence of the Introduction. There is no results summary, no data list and no URL. The acronym is not expanded. The first Results subsection then demonstrates the remaining technical obstacle: "Prevalence of monoallelic detection in scRNA-seq data presents challenges".

**Example, figure callout, analogy:** none. Bulk CNV methods serve as the precedent ("Likewise, analysis of allelic imbalance may also be informative…").

**Discussion, first paragraph (the Discussion has one paragraph):**

- It opens "Altogether, our results demonstrate…".
- It recaps the gains in sensitivity and precision over an expression-based approach.
- It makes a protocol-generality claim: "suitable for both full-transcript-length and 3′-tagging scRNA-seq protocols".
- It ends on biological interplay.

**Characteristic quotes**

- "In that regard, the ability to examine the transcriptional states of genetically distinct intratumoral subclones would be helpful for evaluating the likely functional impact of associated subclonal mutations."
- "However, integration with genotypic information at the single-cell level is necessary to establish correspondence between transcriptionally distinct subpopulations and genetic subclones."
- "At the same time, simultaneous unbiased assessment of DNA and RNA from an individual cell remains challenging (Dey et al. 2015; Macaulay et al. 2015; Wang et al. 2017)."
- "Likewise, analysis of allelic imbalance may also be informative about CNVs in the context of scRNA-seq."

## 3. MERINGUE. Miller et al. 2021, *Genome Research*

- **Source:** Europe PMC XML, published version. FULL TEXT.
- **Size:** 3 paragraphs with no heading, 638 words (112 / 145 / 381).

**Paragraph roles**

1. Biological context. Spatial context matters (brain, development, cancer). Dissociation loses it, and spatial organisation "remains an open question".
2. Technologies, imaging-based and capture-based with 10–100 µm pixels, and what analysis must therefore do: be scalable, and "accommodate the nonuniform cell density" from both biological and technical causes.
3. Named prior methods, each given one neutral "Briefly" sentence on its mechanism: SpatialDE, Trendsceek, SPARK. The paragraph then:
   - names their shared assumption (direct use of cell distances);
   - gives two failure modes, with Supplemental Fig. S1A and S1B;
   - adds SpaOTsc, which lacks a significance test;
   - says 3D or multi-section analysis "have yet to be shown";
   - gives the tool sentence and a one-sentence result.

**How prior methods are criticised**

- Named tools are summarised fairly, one by one, before one shared, structural limitation is drawn.
- Tone: conditional and polite.
- Limitation verbs: "could present a challenge for analyses within tissues where…", "may give rise to apparent spatial aggregation", "could also potentially mask", "does not provide a statistical framework to distinguish…", "have yet to be shown".

**How the gap is framed**

- Tools fail because of an assumption: distance-based statistics confound density with expression.
- The new technology creates the demand.
- 3D is an unmet capability.

**How the tool is introduced**

> "Here, we developed MERINGUE, a density-agnostic method for identifying spatial gene expression heterogeneity using spatial autocorrelation and cross-correlation analyses."

It is followed by one sentence of results and generality: "Using a variety of spatially resolved transcriptomic data sets, we show that MERINGUE is able to identify biologically relevant spatial gene expression patterns in both 2D and 3D in a manner that is independent of cell density." The URL is given only in "Software availability".

**Example, figure callout, analogy**

- Running example: a hypothetical failure (dense regions creating false aggregation).
- Figure callout: Supplemental Fig. S1A and S1B.
- Analogy: none.
- It also gives a formal definition of the target quantity ("spatial gene expression heterogeneity, defined as an uneven, aggregated, or patterned spatial distribution of gene expression magnitudes").

**Discussion, first paragraph**

- It opens "Spatially resolved transcriptomic measurements demand computational approaches…" and "Here, we presented MERINGUE as a density-agnostic approach…".
- It then recaps validation across platforms in 2D and 3D, the MERFISH biology, scalability, and robustness to density.

**Characteristic quotes**

- "Thus, how these subpopulations of cells are organized in space and how they may interact with each other remains an open question in many systems."
- "Thus, each method directly incorporates information regarding cell distances, which could present a challenge for analyses within tissues where cells are distributed with nonuniform densities."
- "It is, however, important to identify variations in gene expression magnitudes across cells that do not arise from variations in cellular density."
- "Furthermore, cells in tissues inherently exist in a three-dimensional context, yet computational approaches capable of taking into consideration z-axis information, often at differing length scales such as multiple noncontiguous tissue sections, have yet to be shown."

## 4. STdeconvolve. Miller et al. 2022, *Nature Communications*

- **Source:** Europe PMC XML, published version. FULL TEXT.
- **Size:** "Introduction" heading, 3 paragraphs, 642 words (64 / 280 / 298). Fig. 1, an overview schematic, sits inside the section.

**Paragraph roles**

1. Context and data limit in four sentences. Pixels are mixtures of cells, which "hinders" cell-type-specific analysis.
2. Named reference-based tools (SPOTlight, RCTD, SpatialDWLS), each given one sentence on its mechanism. Their shared dependence on a reference follows, then the reasons a reference may be absent or mismatched, with a concession that atlases help. The paragraph concludes that a reference-free route is needed.
3. Tool sentence with URL and Fig. 1, followed by:
   - the LDA analogy from natural language processing;
   - four "unique aspects of ST data" that the method leverages;
   - feature selection;
   - the a-priori choice of K.

**How prior methods are criticised**

- Named tools are summarised neutrally. The criticism targets the input they "rely on", not their performance.
- Concession: atlases "may help alleviate such limitations particularly for healthy tissues".
- Limitation phrases: "rely on the availability of a suitable single-cell reference, which may present limitations if…", "may still result in systematically different gene expression quantifications", "may result in missing or inconsistent cell types", "affecting reference-based deconvolution accuracy".

**How the gap is framed**

- A data limit: multi-cellular pixels.
- A dependency of current tools: they need an external reference.

**How the tool is introduced**

> "Here, we developed STdeconvolve (available at https://github.com/JEFworks-Lab/STdeconvolve and as Supplementary Software) as a reference-free, unsupervised approach for deconvolving multi-cellular pixel-resolution ST data (Fig. 1)."

It is followed by the mechanism and an explicit statement of the modelling assumptions, framed as properties of the data that the method exploits. It also states a user-set parameter (K) and how K is guided. There is **no results summary** in the Introduction; the paragraph ends on K.

**Example, figure callout, analogy**

- Running example: none.
- Figure callout: Fig. 1.
- Analogy: yes, a sustained one. Words, documents and topics map to genes, pixels and cell types.

**Discussion, first paragraph**

- It restates the need and names the competitors' dependency ("As we have shown, this reliance…").
- It then says "Here, we have presented STdeconvolve…" and recaps every Results section in order: accuracy, competitiveness when a reference is ideal, clustering versus deconvolution, perturbation, breast cancer.
- Limitations sit in a later paragraph: "its use of LDA modeling relies on several underlying assumptions, which may present limitations when these assumptions are not satisfied".

**Characteristic quotes**

- "This lack of single-cell resolution hinders the characterization of cell-type-specific spatial organization and gene expression variation."
- "As such, these approaches rely on the availability of a suitable single-cell reference, which may present limitations if such a reference does not exist due to budgetary, technical[6], or biological limitations[7]."
- "As such, a reference-free deconvolution approach provides an alternative strategy for deconvolving cell types when an appropriate reference is not available."
- "In the context of ST data, given a count matrix of gene expression in multi-cellular ST pixels, STdeconvolve applies LDA to infer the putative transcriptional profile for each cell type and the proportional representation of each cell type in each multi-cellular ST pixel (“Methods”)."
- "Briefly, these unique aspects of ST data include (i) the limited number of cells and cell types represented in each ST pixel, (ii) the limited impact of batch effects on the measured gene expression across pixels, (iii) the large number of pixels compared to cell types, and (iv) the likely heterogeneity of cell-type proportional distribution across pixels in tissues."

## 5. VeloViz. Atta et al. 2022, *Bioinformatics*

- **Source:** PMC HTML of the published version. FULL TEXT.
- **Size:** "1 Introduction", 3 paragraphs, about 339 words (78 / 168 / 93).

**Paragraph roles**

1. Context. scRNA-seq is a static snapshot. Trajectories are inferred from continua, and RNA velocity predicts future states.
2. Current practice projects velocity onto PCA, t-SNE, UMAP or diffusion-map embeddings. These embeddings disagree with each other and break when intermediate states are missing. A concession follows: non-visual driver-gene analyses exist, but visualisation "remains an important approach".
3. Tool sentence and a two-sentence results summary.

**How prior methods are criticised**

- By category of embedding. Tools are named only as the embeddings themselves.
- Limitation verbs: "can yield different representations of the underlying cellular trajectory", "may be unable to capture global relationships … thereby hindering downstream interpretation".

**How the gap is framed**

- Tools fail under a data limit: missing or rare intermediate states.
- Tools are also inconsistent with each other.

**How the tool is introduced**

> "Here, we developed VeloViz to visualize cellular trajectories by incorporating information about each cell’s predicted future transcriptional state inferred from RNA velocity analysis."

It is followed by "Using both real and simulated data … we demonstrate that VeloViz embeddings are better able to consistently capture underlying cellular trajectories across diverse trajectory topologies compared with other evaluated methods". A second sentence covers missing intermediates. There is no URL.

**Example, figure callout, analogy:** none.

**Discussion, first paragraph**

- It opens "To facilitate better visual representation…, we developed VeloViz…".
- It recaps the results and warns about the misreading that t-SNE and UMAP invite ("erroneous interpretation that these cells are part of biologically distinct subpopulations").
- The next paragraph opens "However, several limitations of VeloViz should be considered…".

**Characteristic quotes**

- "Current technologies for high-throughput single-cell transcriptomics profiling provide a static snapshot of the transcriptional states of individual cells."
- "However, these approaches can yield different representations of the underlying cellular trajectory."
- "Although alternative non-visual methods such as identifying dynamic driver-genes have been developed to help interpret information from RNA velocity analysis (Bergen et al., 2020), visual representation of cellular trajectories remains an important approach to understanding the overall relationships between cell states."
- "Likewise, given simulated and sparsely sampled cellular trajectories with missing intermediate cell states, we find that VeloViz embeddings are able to more robustly retain the overall cell state relationships in the underlying trajectories compared with other evaluated methods."

## 6. SEraster. Aihara et al. 2024, *Bioinformatics*

- **Source:** Europe PMC XML, published version. FULL TEXT.
- **Size:** "1 Introduction", 3 paragraphs, 349 words (128 / 96 / 125).

**Paragraph roles**

1. Context, with downstream tool categories listed: SVG tools and cell-type interaction tools. Their compute grows with N, so preprocessing is "needed".
2. Precedent in scRNA-seq: subsampling, metacells, and self-organising maps for spatial data. The precedents are presented positively as solutions to a similar problem.
3. Tool sentence, implementation (sparse matrices, parallelism), benchmarks, availability (two GitHub URLs) and ecosystem fit (SpatialExperiment, Seurat).

**How prior methods are criticised**

- By category, and only on scale: "many of these computational tools have runtime and memory requirements that increase with the number of single cells or spatial points analyzed".
- The tone is neutral. SEraster is positioned as complementary to existing tools, not as a competitor.

**How the gap is framed:** data scale outgrowing tools, which is a technology-driven limit.

**How the tool is introduced**

> "Here, to enhance the scalability of spatial omics data analysis, we developed a preprocessing framework called SEraster to aggregate spatially proximal cells into pixels using rasterization prior to downstream analysis."

It is followed by:

- implementation details;
- "We benchmarked the performance of SEraster on downstream spatial omics data analyses including identifying SVGs and cell-type co-enrichment…";
- the URLs and package infrastructures.

**Example, figure callout, analogy**

- Running example: none.
- Figure callout: none.
- Analogy: a historical one. "A similar scalability challenge previously emerged for … scRNA-seq".

**Discussion, first paragraph**

- It opens with a context sentence, then "we developed SEraster to use rasterization…".
- It recaps runtime, accuracy, comparison with down-sampling, and multi-resolution co-enrichment.
- A later paragraph opens "Several limitations of SEraster should be considered…".

**Characteristic quotes**

- "However, many of these computational tools have runtime and memory requirements that increase with the number of single cells or spatial points analyzed, presenting challenges as technologies continue to improve and researcher apply them to generate large-scale spatial omics data with millions of spatial points." (The typo "researcher" is in the original.)
- "Therefore, preprocessing to streamline such large-scale spatial omics data analyses is needed."
- "A similar scalability challenge previously emerged for the analysis of single-cell RNA sequencing (scRNA-seq) datasets."
- "SEraster further implements sparse matrix representations and parallel processing for enhanced efficiency."

## 7. STcompare. Clifton et al. 2026, *Bioinformatics* (bioRxiv 2025.11.21)

- **Source:** Europe PMC XML, published version. FULL TEXT.
- **Size:** "1 Introduction", 3 paragraphs, 429 words (92 / 156 / 181). Fig. 1, an 8-panel overview, sits inside the section.

**Paragraph roles**

1. Context and motivation: finding genes whose spatial pattern changes in disease against control.
2. Two categories of existing analysis, differential gene expression (DGE) and spatially variable genes (SVG), and the logical reasons they cannot detect a change in pattern. Each reason is spelled out as a false-negative and false-positive case.
3. A toy simulation with test statistics: genes A, B and C, radial patterns, Wilcoxon and Moran's I values (Fig. 1A). It proves the logical point from paragraph 2.

**How prior methods are criticised**

- By category (DGE, SVG), not by tool name.
- The criticism is logical: these methods answer a different question.
- Tone: firmer than elsewhere in the lab.
- Limitation phrases: "does not consider the spatial information", "fail to characterize changes in spatial expression patterns", "significance cannot be used to distinguish…", "do not reflect changes or lack of changes in spatial expression patterns".

**How the gap is framed:** conceptual. Existing tests target a different estimand.

**How the tool is introduced**

- It is **not introduced in the Introduction.** The Introduction ends on the toy-simulation conclusion.
- The tool is named in the first sentence of section 2, Methods:

  > "To enable comparative analysis for ST datasets that considers spatial information in structurally matched tissues, we developed STcompare, available as an R package at https://github.com/JEFworks-Lab/STcompare."

- That sentence is followed by the two tests and the preprocessing prerequisites: alignment, including the lab's STalign, and rasterization with SEraster.

**Example, figure callout, analogy**

- Running example: yes, genes A, B and C, which recur in Fig. 1C–D.
- Figure callout: Fig. 1A.
- Analogy: none.

**Discussion, first paragraph**

- It opens "Comparative spatial transcriptomics data analysis is needed…" and "Here, we presented STcompare…".
- It describes both tests and their rationale.
- A later paragraph says "we acknowledge that there are limitations". It names the Gaussian-kernel assumption and how a user would detect that the assumption fails.

**Characteristic quotes**

- "However, this type of DGE analysis does not consider the spatial information."
- "As such, a gene can be identified as not differentially expressed across two datasets by having the same mean gene expression despite having distinct spatial expression patterns."
- "As such, existing DGE and SVG analyses fail to characterize changes in spatial expression patterns."
- "To demonstrate these cases, we simulated ST datasets where one pair of genes, A and B, had different radial expression patterns but the same mean expression, while the other pair, A and C, had same radial expression patterns but different mean expression (Fig. 1A)."
- "These results demonstrated that bulk DGE and SVG analysis results do not reflect changes or lack of changes in spatial expression patterns."

## 8a. STARIT. Velazquez … Fan 2026, *PLoS Computational Biology* (extra)

- **Source:** Europe PMC XML, published version. FULL TEXT.
- **Size:** "1 Introduction", 4 paragraphs, 617 words (142 / 192 / 245 / 38).

**Paragraph roles**

1. A technology explainer: how imaging-based spatial transcriptomics (imSRT) decodes transcripts and how cells are segmented.
2. Conventional practice: per-cell counts (Fig. 1A). These approaches "inherently neglect" subcellular localisation. The paragraph says why localisation matters, gives a disease example, notes that frameworks "remain limited", and states the hypothesis.
3. Tool sentence (Fig. 1B), then:
   - a step-by-step mechanism;
   - positioning as an "extension" of prior rasterisation to the subcellular scale;
   - the downstream pipeline.
4. One-sentence results summary.

**How prior methods are criticised**

- By practice ("conventional analysis approaches"), with no tool named.
- Limitation phrases: "inherently neglect", "remain limited", "conventional count-based approaches fail to resolve".

**How the gap is framed:** a representation problem, which the paper calls a "knowledge gap". The standard data structure discards information.

**How the tool is introduced**

> "As an alternative way to represent imSRT data to incorporate subcellular molecular spatial information and address this knowledge gap, we introduce STARIT (Spatial Transcriptomics As Rasterized Image Tensors) (Fig 1B)."

It is followed by the mechanism and lineage, then the final paragraph: "Using simulated and real imaging-based spatial transcriptomics datasets, we show that STARIT can distinguish cellular variation … that conventional count-based approaches fail to resolve while also recovering expected cell-type-specific transcriptional heterogeneity in real datasets."

**Example, figure callout, analogy**

- Running example: a disease example (synaptic RNA mislocalisation in Alzheimer's disease).
- Figure callouts: Fig. 1A and Fig. 1B.
- Analogy: count analysis is "analogous to conventional single-cell RNA-sequencing data analyses".

**Discussion, first paragraph**

- It restates the gap, then "We developed STARIT…".
- It recaps simulations (robustness, trajectories) and real data (osmFISH, bacterial MERFISH).
- A later paragraph opens "Despite these advances, several limitations remain". It states the core assumption explicitly: "STARIT assumes that the observed two-dimensional transcript pattern snapshot contains enough information…".

**Characteristic quotes**

- "Conventionally, analysis of such imSRT data has relied on assigning and tallying the number of transcripts detected for each gene per segmented cell (Fig 1A)."
- "However, these conventional analysis approaches inherently neglect the subcellular heterogeneity of these transcripts within cells such as potential nuclear, perinuclear, cytoplasmic, polarized, and other compartment-associated localizations."
- "While imSRT offers transcriptome-wide profiling to map subcellular RNA heterogeneity, computational frameworks capable of systematically incorporating such information remain limited."
- "Such use of rasterization to create image-based representations has previously been explored in spatial omics applications at the cellular level and hence STARIT represents an extension of these similar principles to the subcellular level [7,8]."

## 8b. Gene count normalization in im-SRT. Atta … Fan 2024, *Genome Biology* (extra)

This is an evaluation paper, not a tool paper.

- **Source:** Europe PMC XML, published version. FULL TEXT.
- **Size:** "Background", 4 paragraphs, 653 words (91 / 146 / 174 / 242).

**Paragraph roles**

1. Technology context: im-SRT is aggregated to cells by segmentation.
2. A catalogue of normalization methods (library size, scTransform, DESeq2, TMM, cell volume) and the data each was designed for.
3. Gene panels are pre-selected and often skewed toward cell types of interest.
4. Study aim (Fig. 1), then:
   - an "ideal" invariance criterion, with a hypothetical gene example;
   - "However, we show…" results;
   - generalisation across technologies and tissues;
   - mitigation with larger panels;
   - a caution and a decision tree.

**How prior methods are criticised**

- Named methods are introduced neutrally, with only the note that they "were designed for" other data.
- The failure is shown, not asserted: "may introduce errors in downstream analyses".
- Tone: cautionary ("we caution").

**How the gap is framed:** a broken assumption. The experimental design (skewed panels) breaks the assumptions behind count-based scaling.

**How the work is introduced**

> "Here, we investigate how different gene count normalization procedures may impact downstream transcriptomic analyses for im-SRT data with different gene panels (Fig. 1)."

It is followed by an ideal-world criterion, results, a generality claim ("We further demonstrate that these concerns generalize across diverse im-SRT technologies and tissues") and a practical deliverable (the decision tree).

**Example, figure callout, analogy**

- Running example: a hypothetical gene that is differentially expressed in one panel should stay so in another.
- Figure callout: Fig. 1.
- Analogy: none.

**Discussion, first paragraph**

- It opens "Here, we investigate…" and recaps the findings.
- It ends on generality: "consistent across multiple simulated skewed gene panels across multiple im-SRT data sets spanning several im-SRT technologies and tissues, demonstrating that the results are independent of upstream, technology-specific technical factors".
- Later paragraphs list each method's assumptions one by one, for example "library size normalization assumes…" and "DESeq2 assumes…".

**Characteristic quotes**

- "Ideally, biological conclusions drawn from im-SRT studies in the same biological system would be robust to experimental design choices such as gene panel composition."
- "As such, these gene panels may be skewed to enable more thorough characterization of the cell types of interest at the expense of other cell types present in the tissue."
- "We further demonstrate that these concerns generalize across diverse im-SRT technologies and tissues."
- "Overall, we caution that the choice of normalization method and gene panel may impact the biological interpretation of the im-SRT data, and present readers with a decision tree to help guide the choice of gene expression normalization method for their im-SRT data."

---

## Synthesis

**What is consistent**

- **Structure.** 3–4 paragraphs, 340–650 words (median 477), in one funnel:
  - *X technologies enable Y*;
  - analysis "demands" or "is needed";
  - prior work summarised fairly, one sentence per tool or category;
  - "As such" or "Thus", then one shared limitation;
  - *Here, we developed X, a [property] method for [task] using [mechanism]*;
  - zero to two sentences of results.

  The tool appears only in the last or second-to-last paragraph; STcompare defers it to Methods.
- **Tone.**
  - Hedged: "may" appears 5–8 times in PAGODA, STdeconvolve and normalization, but never in SEraster or STcompare.
  - Criticism targets an input or assumption that a whole category shares, not a named tool's accuracy: "rely on the availability of a suitable single-cell reference", "directly incorporates … cell distances", "does not consider the spatial information", "inherently neglect".
  - Concessions are common ("may help alleviate").
  - No Introduction uses "novel", "state-of-the-art", or "first" in the priority sense.
  - The method is "an alternative approach" or "an alternative strategy".
- **Generality** is claimed by listing data: "in both 2D and 3D", "real and simulated", "diverse im-SRT technologies and tissues".
- **Assumptions** enter the Introduction in one of three forms:
  - data properties the method leverages (STdeconvolve's (i)–(iv));
  - a precise definition of the target (MERINGUE);
  - an ideal invariance (normalization).

  Limitations wait for a Discussion paragraph (*Several limitations of X should be considered*; "Despite these advances, several limitations remain").
- **Ending.** *we show/demonstrate that X…*. The first Discussion paragraph restates the need, says *Here, we presented X*, and recaps the Results in order.

**Differences by venue** (2–3 papers each; tendencies only).

- **Nature Methods and Nature Communications** (PAGODA, STdeconvolve):
  - a headed Introduction and numbered citations;
  - a cross-domain analogy that explains the model (GSEA rationale; LDA's words, documents, topics);
  - "(Fig. 1)" on the tool sentence;
  - a tool paragraph of 160–300 words;
  - in PAGODA, a biology paragraph that sets up discovery in new data.
- **Genome Research** (HoneyBADGER, MERINGUE):
  - no heading, and author–year citations;
  - HoneyBADGER stops at the tool sentence and shows the obstacle in Results;
  - MERINGUE gives the most detailed tool-by-tool critique.
- **Bioinformatics** (VeloViz, SEraster, STcompare):
  - the shortest Introductions (340–430 words);
  - practical gaps: scale, visualisation, comparison;
  - URLs and package compatibility can sit in the Introduction (SEraster);
  - STcompare proves its gap with a toy simulation that reports statistics.
