# Citation check B: references [22]–[40] (closest prior art)

Draft checked: `docs/paper/manuscript_v1_intro_methods.md` (Introduction and Methods, v1 of 2026-11-03). Citing sentences are taken from `fulltext/claims.md`. Every reference was read in full text, from Introduction through Methods and figure legends; where only a preprint or the authors' manuscript source could be read, the table says so. Downloads are in `fulltext/dl_B/`, and the parsing scripts are in `scripts_B/`, run with `python -I`. Every verbatim quote below was checked by script against the downloaded text (`scripts_B/verify_quotes.py`). Every source quote was found; the only misses were our own manuscript text and proposed wording. The earlier notes in `intro_research.md` were re-verified against these fresh copies rather than reused.

## Summary table

| Ref | Access route actually used | Full text read | Verdict | Proposed fix |
|---|---|---|---|---|
| [22] Traag 2019 (Leiden) | Europe PMC OA fullTextXML, PMC6435756 | Y | Methods: SUPPORTED. Intro: PARTLY (a general network-community algorithm; no cells, spots or kidney) | Fix A; Methods unchanged |
| [23] Hao 2021 (Seurat v4) | Europe PMC fullTextXML, PMC8238499 (author manuscript with STAR Methods) | Y | PARTLY (clusters dissociated PBMCs on a WNN graph by SLM, then annotates by markers; not spatial, not kidney) | Fix A. Optional: cite Scanpy (Wolf 2018), which the pipeline actually uses |
| [24] Zhao 2021 (BayesSpace) | NIH author manuscript NIHMS1768714 via NCBI efetch, PMC8763026 (Europe PMC: HTTP 500) | Y | PARTLY (each spot gets one of q clusters under a Potts prior on the spot lattice; clusters, not segment labels; no kidney) | Fix A |
| [25] Hu 2021 (SpaGCN) | Research Square preprint rs-119776/v1 PDF, via the Semantic Scholar openAccessPdf link. Nat Methods version paywalled; bioRxiv rate-limited (Cloudflare 1015) | Y (preprint only) | PARTLY (domains from a graph built on x/y plus histology; no kidney) | Fix A |
| [26] Long 2023 (GraphST) | Europe PMC OA fullTextXML, PMC9977836 | Y | PARTLY (spot kNN graph from Euclidean x/y; one domain per spot; no kidney) | Fix A |
| [27] Chitra 2025 (GASTON) | NIH author manuscript NIHMS2093871 via NCBI efetch, PMC12370209 (Europe PMC: HTTP 500) | Y | SUPPORTED, Intro and Methods (isodepth is a smooth function of spot x/y; GLM-PC input; 30 random initializations) | None. Optional: "as a smooth function of tissue x/y position" |
| [28] Ma 2022 (Belayer) | NIH author manuscript NIHMS1837222 via NCBI efetch, PMC9814896 (Europe PMC: HTTP 500) | Y | SUPPORTED (relative depth Φ(x, y) by conformal maps; curved layers need approximate boundaries as input) | None. Optional: "in layered tissues" |
| [29] Miao 2026 (CoPro) | Europe PMC fullTextXML, PMC13131544 (bioRxiv preprint v1, 2026-04-21; not peer reviewed) | Y (preprint) | PARTLY. Segment order is an input only in its supervised kidney run (ordinal PTS1 = 1 … PC/IC = 7; Kendall τ 0.74–0.78). CoPro needs tissue x/y and cell-type labels in every mode | Fix B |
| [30] Clifton 2023 (STalign) | Europe PMC OA fullTextXML, PMC10709594 | Y | SUPPORTED | None (see flag 2 on our first sentence) |
| [31] Clifton 2026 (STcompare) | Europe PMC OA fullTextXML, PMC13577157 (published version) | Y | SUPPORTED ("dependent on the datasets under comparison having inherent structural similarity") | None. Optional AKI clause (flag 5) |
| [32] Trapnell 2014 (Monocle) | Europe PMC XML of NIH author manuscript NIHMS570726 (main text and legends), plus its Supplementary Figures PDF (PMC Open Data S3). Online Methods paywalled | Y (main text); N (Online Methods) | SUPPORTED except "alone" (path direction is left to the user) | Fix C |
| [33] Haghverdi 2016 (DPT) | bioRxiv preprint v2 (2016-03-02) JATS, including Online Methods. Published text not read: no PMCID, MDC accepted version HTTP 401, Semantic Scholar HTTP 429. Published abstract read (Europe PMC) | Y (preprint only) | Intro: PARTLY (needs a root cell). Methods: SUPPORTED; "15 diffusion components" is a Scanpy setting, while the paper's DPT is full-rank | Fix C; Methods fix M1 |
| [34] Street 2018 (Slingshot) | Europe PMC OA fullTextXML, PMC6007078, including Table 1 | Y | Intro: PARTLY (needs clusters and a start cluster). Methods: SUPPORTED | Fix C; Methods fix M2 |
| [35] Wolf 2019 (PAGA) | Europe PMC OA fullTextXML, PMC6425583 (Additional file 1 not read) | Y | PARTLY (orders on top of a partition, by distance from a root cell) | Fix C |
| [36] Saelens 2019 (benchmark) | Authors' manuscript source `manuscript/paper.Rmd`, GitHub dynverse/dynbenchmark. bioRxiv v1 blocked (Cloudflare 1015 ×4); Ghent copy "UGent only". Published abstract read | Y (authors' manuscript source, not the typeset paper) | PARTLY ("A large portion of the methods requires some form of prior information (e.g. a start cell)") | Fix C |
| [37] Niec 2022 (SpaceFold) | NIH author manuscript NIHMS1809914 via NCBI efetch, PMC9271639 (Europe PMC: HTTP 500) | Y | SUPPORTED (BayesPrism cell-type fractions projected to 1D by PHATE; no x/y; qualitative validation only) | Optional: "Visium spots along intestinal crypts" |
| [38] Karin 2023 (scPrisma) | Europe PMC OA fullTextXML, PMC10635821 | Y | PARTLY (no testis data; it supplies the cyclic template that [39] used; linear templates also exist, so the cycle is the user's choice) | Fix D |
| [39] Chakravorty 2026 | NIH author manuscript via NCBI efetch, PMC13326559 (Europe PMC: HTTP 500) | Y | SUPPORTED (hand-segmented tubule aggregates ordered by scPrisma's de novo cyclic pipeline) | Fix D (split the [38]/[39] roles) |
| [40] Asami 2026 | NIH author manuscript NIHMS2193327 via NCBI efetch, PMC13379681 (Europe PMC: HTTP 500). First full-text read; the earlier notes had the abstract only | Y | SUPPORTED. Caution: its Discussion names **kidney nephrons** as a target for "spatial pseudotiming" | Fix F (one Discussion sentence) |
| "These approaches rely on a change in cell-type composition, a cyclic topology or three-dimensional data" | — | — | Holds for [37]–[40] taken separately, but the referent can be read as including [32]–[36], which rely on none of the three | Fix E |

## Proposed wording (exact)

**Fix A: clustering sentence ([22]–[26]).** None of the five papers mentions kidney, renal or nephron (text search of all five full texts). All five output clusters or domains; segment names come from a separate annotation step. The three spatial-domain methods push neighbouring spots in the tissue plane into the same group. That prior is exactly what fails when neighbouring cross-sections belong to different nephrons, and the current sentence does not use it.
> "As such, ST analyses of the kidney typically discretize[kidney refs among 10–15 that cluster spots or cells]. Clustering assigns each cell, spot or cross-section to one of a few discrete groups that are then annotated as segments[22,23], spatial-domain methods additionally encourage neighbouring spots in the tissue plane to share a group[24–26], and comparisons between conditions or species then contrast group means."

Shorter alternative: "Clustering and spatial-domain methods assign each spot or cross-section to a discrete cluster or domain[22–26], which is then annotated as a segment, and comparisons between conditions or species then contrast segment means."

**Fix B: CoPro ([29]).**
> "Methods that model continuous expression gradients, such as GASTON and Belayer, build their coordinate from tissue x/y positions[27,28], and CoPro, which also uses tissue coordinates and cell-type labels, was given the known order of nephron segments as an input when applied to the kidney[29]."

Shorter alternative: "Methods that model continuous expression gradients build their coordinate from tissue x/y positions[27–29]; applied to the kidney, CoPro was also given the order of nephron segments as an input[29]."

**Fix C: trajectory inference ([32]–[36]).** "From expression alone" overstates. DPT and PAGA need a root cell, Slingshot needs clusters and a start cluster, and Monocle leaves the direction to the user. Saelens states that "a large portion" of methods need prior information. Our own curve is also oriented by a marker-axis root, so "alone" does not set up a contrast with our method.
> "Trajectory inference orders cells along continuous processes from their expression profiles, usually given a start cell or cluster that sets the direction[32–36], and related ideas have ordered repeated structures: …"

**Fix D: testis clause ([38], [39]).**
> "… testis tubules around the known cycle of the seminiferous epithelium[39], using a cyclic topological template[38], …"

**Fix E: the "These approaches" sentence ([37]–[40]).**
> "These repeated-structure approaches rely on a change in cell-type composition, a known cyclic topology or three-dimensional capture of whole structures."

("Three-dimensional capture of whole structures" is more exact than "three-dimensional data". Asami's ordering itself uses only per-follicle counts; the 3D imaging is what ensures each follicle is complete. Partial follicles were excluded.)

**Fix F: Asami in the Discussion ([40]).** Asami et al. write: "This approach should be broadly applicable to systems where similar structures develop asynchronously, including kidney nephrons, islets of Langerhans, and lung alveoli." Proposed Discussion sentence:
> "Spatial pseudotiming orders whole structures that develop asynchronously by developmental age[40]; pseudostructure reconstruction instead orders cross-sections of mature structures by position along them, from 2D sections in which no structure is captured whole."

**Methods fix M1: DPT ([33]).**
> "…were ordered by diffusion pseudotime[33] as implemented in Scanpy 1.11.5[ref Wolf 2018] (max(30, √n) neighbours; 15 diffusion components; rooted at the same root unit), …"

**Methods fix M2: Slingshot ([34]).** Per `results/paper_figures/fig3/baselines/methods_notes.md`, Slingshot was given the root's cluster (`start.clus`), not the root unit. It produced 5–6 lineages, averaged with `slingAvgPseudotime`.
> "Slingshot 2.14.0[34] on the pass-2 Leiden clusters, starting from the root unit's cluster, with each unit's per-lineage pseudotimes averaged by its lineage weights (R 4.4.3)"

## Other flags from the full texts

1. **Novelty boundaries.**
   - CoPro already compares samples on a shared expression-derived axis without spatial registration ("projecting all cells from different samples onto a biologically defined axis derived from shared gene expression programs"). Asami already compares genotypes on one jointly computed pseudotime. Our shared-axis comparison should not be presented as new in itself. What is new is building the axis from expression alone, with no x/y, labels, reference program or whole-structure capture.
   - GASTON already "pools" information across spatially separated locations of equal isodepth, but only through a smooth field over one slice.
2. **Our first sentence is almost verbatim STalign's.** STalign's abstract opens "Spatial transcriptomics (ST) technologies enable high throughput gene expression characterization within thin tissue sections." Our abstract and Introduction open with nearly the same words. As lab house style this is reasonable, but the authors should be aware of it.
3. **"The ends of the structure are not marked"** (Introduction) sits uneasily with our marker-axis root (flagged by the [32]–[36] check). Consider "no landmark fixes the ends", or say that orientation comes from an early→late marker axis.
4. **Uniform-sampling caveat.** Chakravorty can justify evenly spaced angles because tubule stages are spatially random. A kidney section samples nephron positions roughly in proportion to the tubule length lying in the plane (convoluted S1/S2 versus straight S3). Say that t is an order-preserving position, not physical length. Asami's limitation, "pseudotime, which may not have a linear relationship with real time", is the model sentence.
5. **STcompare is the direct kidney AKI precedent** (Visium, 24 h IRI against sham, 200-µm pixels, compartments cortex / outer medulla / inner medulla). It reports PT transporters (Acsm2, Slc34a1) "remained localized to the cortex, but expression was greatly reduced". The AKI section could say that matched-location comparison resolves the kidney at compartment scale, while we compare positions along the PT.
6. **Validation bar set by the repeated-structure precedents.** SpaceFold was validated only qualitatively (expected cell-type order, zone markers, SOX9). Chakravorty used composition, published stage frequencies and DAPI morphology. Asami correlated pseudotime with follicle length (Spearman ρ = 0.94 and 0.90). None compares its axis with a measured position along the structure, so our segment-order ρ and the held-out-gene and downsampling stability tests are a fair point of comparison.

---

# Per-reference sections

## [22] Traag, Waltman & van Eck 2019, From Louvain to Leiden. *Sci. Rep.* 9, 5233

**1. Access.** Europe PMC fullTextXML, PMC6435756 (open access). Read: Introduction, Louvain/Leiden algorithm sections, Guarantees, Experimental Analysis (benchmark and empirical networks), Discussion and figure legends. The Supplementary Information (proofs) was not read.

**2. Our citing sentences.**
- Intro: "Clustering and spatial-domain methods assign each spot or cross-section to a segment label[22–26], and comparisons between conditions or species then contrast segment means."
- Methods: "Units were clustered by Leiden[22] (resolution 2, 30 nearest neighbours, 2,000 highly variable genes), and each cluster was labelled from curated marker panels ..."

**3. Verdict.**
- Methods: **SUPPORTED.** Leiden is the algorithm we ran, and γ is its resolution parameter.
- Intro: **PARTLY SUPPORTED.** This is a general community-detection paper. Its empirical networks are Amazon, DBLP, IMDB, Live Journal, Web of Science and Web UK. It never mentions cells, spots or gene expression. As a citation for the clustering algorithm it is standard and acceptable. It does not support "assign each spot or cross-section to a segment label". Fix: Fix A, where [22] sits with "Clustering assigns each cell, spot or cross-section to one of a few discrete groups".
- Our wording does not claim the connectivity guarantee. If a later draft does, the paper's own hedge should go with it (third quote).

Quotes:
- Abstract: "We prove that the Leiden algorithm yields communities that are guaranteed to be connected."
- Introduction: "Higher resolutions lead to more communities, while lower resolutions lead to fewer communities."
- Experimental Analysis > Badly connected communities: "The Leiden algorithm guarantees all communities to be connected, but it may yield badly connected communities."

**4. Writing ideas.**
- "Higher resolutions lead to more communities, while lower resolutions lead to fewer communities." Use it in **Methods** to explain why resolution 2 over-splits before the marker-panel labelling merges clusters.
- "nodes in a community in biological or neurological networks are often assumed to share similar functions or behaviour". Use it in the **Introduction or Discussion** to name the assumption behind contrasting segment means: members of a cluster are treated as equivalent.
- The paper pairs each guarantee with its limit: "guarantees all communities to be connected, but it may yield badly connected communities". Use the same pattern in **Methods** for our identification statement: what the curve guarantees (an order along one axis) and what it does not (absolute direction, endpoints).
- "Such a modular structure is usually not known beforehand." Use a parallel in the **Introduction**: position along the nephron is not known beforehand either.

---

## [23] Hao, Y. et al. 2021, Integrated analysis of multimodal single-cell data. *Cell* 184, 3573–3587.e29

**1. Access.** Europe PMC fullTextXML, PMC8238499 (author manuscript including the STAR Methods). Read: Introduction, all Results sections, Discussion and Limitations, the STAR Methods sections on WNN, clustering and annotation, differential abundance and DE across conditions, and the gradient analysis, plus the figure legends.

**2. Our citing sentence.** Intro: "Clustering and spatial-domain methods assign each spot or cross-section to a segment label[22–26], and comparisons between conditions or species then contrast segment means."

**3. Verdict: PARTLY SUPPORTED.**
- The paper clusters **single dissociated cells** (CITE-seq PBMCs and bone marrow) on a weighted-nearest-neighbour graph. It annotates the clusters by markers (57 clusters in three tiers), then compares conditions within clusters (vaccination time points).
- That supports the generic claim that clustering assigns each unit to a discrete annotated group and that comparisons then contrast groups.
- It is not spatial, not kidney and not segment-based. Its clustering is SLM (Waltman & van Eck 2013), not Leiden.
- It is the Seurat v4 (WNN) paper, so it is an acceptable Seurat citation. Its subject, however, is multimodal integration.
- Fix: see Fix A (top of this report). Optionally replace [23] with Scanpy (Wolf, Angerer & Theis 2018), which matches our pipeline and is already a Methods placeholder (I did not read it).

Quotes:
- STAR Methods > Preprocessing (HVTN PBMC) > Clustering and annotation: "This distance represents the edge weight in a shared nearest neighbor graph (SNN), which we used as input to the graph-based smart local moving (SLM) algorithm"
- Results > Mapping query datasets: "These datasets are typically processed using a workflow that consists of unsupervised clustering, which assumes minimal prior knowledge and is ideally suited for cell type discovery."
- Figure 6A legend: "For each of our level 2 annotated cell clusters, we calculated the number of differentially expressed genes between unvaccinated (day 0) and day 3 samples"

**4. Writing ideas.**
- "NK cells were subdivided into five clusters, representing variation across a continuous landscape." This is the authors of a clustering tool admitting that clusters cut a continuum. Use it in the **Introduction** (the "Segment labels cannot represent a gene that changes within a segment" sentence) or the **Results overview**, where clustering obscures within-segment gradients.
- "the traditional division of NK cells into CD56-bright and CD56-dim categories represents a broader continuum defined by CD16 expression". It parallels S1/S2/S3 as discrete names for a continuum. Use it in the **Discussion**.
- "typically processed using a workflow that consists of unsupervised clustering, which assumes minimal prior knowledge and is ideally suited for cell type discovery". It credits clustering before stating its limit. Use it in the **Introduction**.
- "We merged clusters that did not exhibit clear evidence of separation". It is a model for stating plainly how our resolution-2 Leiden clusters were collapsed into segment labels, in **Methods**.

---

## [24] Zhao, E. et al. 2021, Spatial transcriptomics at subspot resolution with BayesSpace. *Nat. Biotechnol.* 39, 1375–1384

**1. Access.** NIH author manuscript NIHMS1768714, fetched as JATS through NCBI efetch (db=pmc, PMC8763026). Europe PMC fullTextXML returned HTTP 500. Read: Introduction (unheaded), all Results, Discussion, Methods (spatial clustering model, enhanced resolution, simulations) and the figure legends.

**2. Our citing sentence.** As for [23].

**3. Verdict: PARTLY SUPPORTED.**
- The paper assigns each Visium/ST spot to one latent cluster z_i ∈ {1,…,q}, with a Potts prior that pushes lattice neighbours into the same cluster.
- It uses spot positions through the square or hexagonal lattice neighbourhood.
- It is evaluated against manual DLPFC layer annotations and on melanoma, IDC, ovarian and SCC samples. There is no kidney.
- The outputs are clusters. Layer or segment names come from the annotation used as ground truth.
- Fix: see Fix A (top of this report).

Quotes:
- Introduction: "BayesSpace enables spatial clustering by modeling a low-dimensional representation of the gene expression matrix and encouraging neighboring spots to belong to the same cluster via a spatial prior"
- Methods > Spatial clustering model: "ST and Visium spots are arranged on square and hexagonal lattices, which provide a natural way to define a neighborhood structure"
- Methods > Spatial clustering model: "zi ∈ {1,…,q} denotes the latent cluster that i belongs to"

**4. Writing ideas.**
- "encouraging neighboring spots to belong to the same cluster via a spatial prior". This states the prior exactly. Use it in the **Introduction** or the **comparison with other methods**: the prior is wrong where neighbouring cross-sections belong to different nephrons.
- "The number of clusters q is determined by prior biological knowledge when available or otherwise by the elbow of the pseudo-log-likelihood plot". Discrete methods must be told how many groups to find, and our curve needs no k. Use it in the **comparison with other methods**.
- "we use the adjusted Rand index (ARI) to quantify similarity between cluster labels and manual annotations, which are considered the ground truth". Use it in **validation**: name what we score against (segment order) and say it is a reference, not anatomical truth.
- "Clustering is an important step in the analysis of such data that allows downstream analyses, such as cell type or tissue annotation and differential expression". It credits clustering before the limit, for the **Introduction**.

---

## [25] Hu, J. et al. 2021, SpaGCN. *Nat. Methods* 18, 1342–1351

**1. Access.**
- **Preprint only:** Research Square rs-119776/v1 (posted 10 Dec 2020, CC BY), 41-page PDF, from the Semantic Scholar `openAccessPdf` link. Read: Abstract, Introduction, all Results sections, Discussion and the full Methods (graph construction, distance, GCN, iterative clustering, SVG detection, Moran's I).
- The published version at nature.com is paywalled; one fetch returned only the abstract and "Access options". There is no PMC copy.
- The bioRxiv preprint (10.1101/2020.11.30.405118) returned Cloudflare 1015 rate-limit errors.
- The published abstract matches the preprint except that it reports "seven SRT datasets" (the preprint has five). The two added datasets were not checked.

**2. Our citing sentence.** As for [23].

**3. Verdict: PARTLY SUPPORTED** (preprint).
- SpaGCN builds a spot graph from (x, y) coordinates plus a histology-derived third coordinate, aggregates expression over spatial neighbours, and clusters the spots. Each cluster is a spatial domain.
- Datasets: mouse olfactory bulb, mouse posterior brain, DLPFC, pancreatic cancer and MERFISH hypothalamus. There is no kidney.
- It outputs domains, not segment labels.
- Fix: see Fix A (top of this report).

Quotes:
- Results > Overview of SpaGCN: "Each cluster is considered as a spatial domain from which SpaGCN then detects SVGs that are enriched in a domain by differential expression analysis"
- Methods > Data preprocessing: "with N samples and D genes, along with the (x, y) 2-dimensional spatial coordinates of each sample"
- Methods > Calculation of distance between two vertices: "SpaGCN considers two spots to be close if and only if 1) the two spots are physically close, and 2) they have similar pixel features"

**4. Writing ideas.**
- "an important step is identifying spatial domains defined as regions that are spatially coherent in both gene expression and histology". This is a clean definition of the category we contrast with, for the **Introduction**.
- "Traditional clustering methods such as K-means and Louvain's method can only take gene expression data as input, and the resulting clusters may not be contiguous". Use it in the **comparison with other methods**: in the kidney the cross-sections of one segment really are not contiguous in the plane, so spatial contiguity is the wrong goal. Our ordering uses expression only, on purpose.
- "provide a comprehensive picture on the spatial gradients in gene expression in tissue" (Introduction). SpaGCN calls domain-level enrichment "gradients". Use it in the **Discussion** to separate domain contrasts from position along a structure, so a reviewer cannot say domain methods already cover gradients.
- The Results open with "Overview of SpaGCN and evaluation", then go dataset by dataset, each with a known layered ground truth. That is a model for our **Results overview** followed by per-dataset sections.

---

## [26] Long, Y. et al. 2023, Spatially informed clustering, integration, and deconvolution of spatial transcriptomics with GraphST. *Nat. Commun.* 14, 1155

**1. Access.** Europe PMC fullTextXML, PMC9977836 (open access). Read: Introduction, all Results, Discussion, the full Methods (graph construction, contrastive learning, spatial domain assignment and refinement, integration, deconvolution, baselines and evaluation) and the figure legends.

**2. Our citing sentence.** As for [23].

**3. Verdict: PARTLY SUPPORTED.**
- GraphST builds a spot kNN graph from Euclidean distances in tissue coordinates and learns spatially smoothed embeddings. It then clusters the spots with mclust, one domain per spot, with an optional step that reassigns each spot to the majority label within radius r.
- Datasets: DLPFC, mouse olfactory bulb, hippocampus, embryo, breast cancer and lymph node. There is no kidney.
- It outputs domains. The number of clusters is set to the number of manual annotations.
- Fix: see Fix A (top of this report).

Quotes:
- Introduction: "Within the spatial transcriptomics analysis workflow, assigning capture spots to spatial domains with unsupervised clustering is an essential task."
- Methods > Graph construction: "for a given spot, its neighbors are defined by its proximity to other spots based on the Euclidean distance computed from the spatial location information."
- Methods > Spatial domain assignment via clustering and refinement: "we use the reconstructed spatial gene expression Hs from the decoder (Fig. 1A) with the nonspatial assignment algorithm, mclust[55], to cluster the spots into different spatial domains."

**4. Writing ideas.**
- "cells are organized into groups of similar cells physically clustered together" (Introduction, first sentence). This is the premise our second paragraph overturns for repeated structures. Use it in the **Introduction**.
- "For tissue slices with manual annotation, we set the number of clusters to be the same as the ground truth." Use it in the **comparison with other methods**: discrete baselines are usually handed the true number of groups, and our method is not.
- "batch correction methods developed for scRNA-seq, such as Harmony and scVI, are not suitable as they only consider gene expression". This anticipates a reviewer objection to our Harmony step. In **Methods or the Discussion**, say that ignoring tissue position is deliberate here, because tissue position does not encode position along the nephron.
- "Multiple tissue sections can also be obtained via serial sectioning of the organ of interest, yielding serial tissue slices that capture 3D information". Use it in the **Introduction**, in the 3D/serial-section sentence: the alternative we avoid.

---

## [27] Chitra, U. et al. Mapping the topography of spatial gene expression with interpretable deep learning. *Nat. Methods* 22, 298–309 (2025)

Verified from the XML: vol. 22, issue 2, pp. 298–309; epub 2025-01-23; 7 authors (Chitra U … Raphael BJ); doi 10.1038/s41592-024-02503-3; PMID 39849132.

**1. Access.** NIH author manuscript NIHMS2093871, JATS from NCBI efetch (db=pmc, PMC12370209). The Europe PMC fullTextXML endpoint returned HTTP 500. I did not read the copy-edited publisher version, so the wording may differ slightly from print. Read: Introduction, all Results sections, Discussion, full Methods (model, conservative vector fields, MLE steps 1–2, training and implementation, model selection, isodepth scaling, comparison methods including SpaceFold), and the legends of Figures 1–6 and Extended Data Figures 1–8.

**2. Our citing sentences.**
- Intro: "Methods that model continuous expression gradients, such as GASTON and Belayer, build their coordinate from tissue x/y positions[27,28], and CoPro takes the order of nephron segments as an input[29]."
- Methods: "The same corrected embedding and root were ordered by diffusion pseudotime[33] (...), Slingshot 2.14.0[34] on the pass-2 clusters (R 4.4.3), scFates 1.2.5 (principal tree)[ref Faure 2023] and GASTON[27] (GLM-PCA input, 30 restarts per section), and each ordering was scored by the same ρ."

**3. Verdict: SUPPORTED (both sentences).** GASTON learns its isodepth d(x,y) as a continuously differentiable function of the spot coordinates of one slice. The input is the count matrix plus the N×2 location matrix, and the network maps coordinates to expression. The Methods sentence matches their implementation: they solve the problem on the top GLM-PCs, with 30 random initializations, keeping the best likelihood.
- Results > GASTON learns gene expression topography: "GASTON trains a neural network to learn a composite function (f∘d)(x,y) from spatial coordinates to gene expression features, where the isodepth d(x,y) corresponds to an interpretable hidden layer"
- Methods > Spatial Topography Problem: "Given SRT data (A, S) and a number P of spatial domains, find a continuously differentiable function d:R2→R and a piecewise linear function h(w)"
- Methods > Training and implementation: "Because of the non-convexity of (10), we use 30 random initializations and select the solution with the largest likelihood."

Optional sharpening. GASTON's coordinate is a **smooth** field over the tissue plane, so nearby locations must have nearby isodepth. Our paragraph 2 argues that neighbouring cross-sections belong to different nephrons. Saying "build their coordinate as a smooth function of tissue x/y position[27,28]" makes the objection specific rather than generic. It is also fair to GASTON, which explicitly allows the same isodepth at separated locations ("the spatial domains Rp need not be contiguous"). So the problem is smoothness, not contiguity.

**4. Writing ideas.**
- **Introduction (concept before tool).** GASTON introduces the coordinate one paragraph before the tool: "We introduce gene expression topography, a fundamentally different approach to modeling spatial variation in gene expression." It then defines the isodepth as "a 1-dimensional coordinate over the tissue slice which describes both the arrangement of spatial domains and the relative position of each spatial location". We could define the pseudostructure position t the same way (one coordinate that carries both segment order and position within a segment) before "To address this, we develop".
- **Introduction (segment-label limitation).** "most existing SRT algorithms assume that mean gene expression is constant within a spatial domain" (Methods). This directly supports our sentence "Segment labels cannot represent a gene that changes within a segment"; [27] could be cited there too.
- **Comparison with other methods / Discussion.** "existing approaches are unable to model continuous variation in physical space" (Discussion, about trajectory inference). Use the mirror image: our t is a position along a structure, neither physical (x/y) space nor process time. Pair it with GASTON's stated assumption, "we assume that all (spatially varying) genes share the same vector field of spatial gradients", to explain why a per-section isodepth in kidney follows the cortico-medullary axis rather than nephron position (our GASTON baseline).
- **Discussion (limitations as stated assumptions).** "the current derivation of isodepth by GASTON relies on two simplifying assumptions that may require adjustment for specific applications". This is a model for stating our identification limits (one unbranched curve; direction set by an external orientation; endpoints not marked).

---

## [28] Ma, C., Chitra, U., Zhang, S. & Raphael, B. J. Belayer: modeling discrete and continuous spatial variation in gene expression from spatially resolved transcriptomics. *Cell Syst.* 13, 786–797.e13 (2022)

Verified from the XML: vol. 13, issue 10, pp. 786–797.e13 (2022-10-19); 4 authors as cited; doi 10.1016/j.cels.2022.09.002.

**1. Access.** NIH author manuscript NIHMS1837222 via NCBI efetch (PMC9814896). Europe PMC returned HTTP 500. Read: Introduction, Results (algorithm, simulations, DLPFC, skin wound, somatosensory cortex), Discussion, STAR Methods (layered-tissue model, axis-aligned problem and its DP, conformal maps, the approximate-boundary and line-boundary cases, parameter selection, GLM-PCA, GP relationship).

**2. Our citing sentence.** Intro: "Methods that model continuous expression gradients, such as GASTON and Belayer, build their coordinate from tissue x/y positions[27,28], and CoPro takes the order of nephron segments as an input[29]."

**3. Verdict: SUPPORTED.** Belayer models expression as a piecewise linear function of a relative depth Φ(x, y), computed from spot coordinates through conformal maps or the heat equation. For curved layers it also needs approximate layer boundaries (from histology or prior anatomy) as input.
- Results > Belayer Algorithm: "Given an N × G transcript count matrix A ... and spatial coordinate matrix S = [si] where si = (xi, yi) are the coordinates of spot i"
- Abstract: "Belayer models gene expression as a piecewise linear function of the relative depth of a tissue layer ... We use conformal maps to model relative depth"
- Results > Belayer Algorithm: "given approximate layer boundaries Γ^ℓ of the tissue slice T, such as from histological images or prior anatomical knowledge, we estimate the real part Re Φℓ"

Optional: "...such as GASTON and Belayer, build their coordinate from tissue x/y positions within one section[27,28]". Or, for Belayer alone, "in layered tissues". Belayer defines "layered tissues, a global model of tissue organization for tissues that consist of consecutive layers of cell types", which a field of convoluted nephrons is not.

**4. Writing ideas.**
- **Introduction (discrete vs continuous framing).** "Most current computational models for analyzing spatial variation in gene expression in SRT data make one of two distinct modeling assumptions." The next sentences name each assumption and its cost. This is a compact template for our "As such, ST analyses of the kidney typically discretize" paragraph.
- **Introduction (within-segment limitation).** "also assume that gene expression is constant within each cluster, and thus do not account for continuous spatial variation". Second support (with [27]) for "Segment labels cannot represent a gene that changes within a segment".
- **Methods / Results (our clusters + curve model).** "Piecewise functions allow for discontinuities in expression ... while continuous functions model gradients of gene expression within a tissue layer". This is the precedent for our "clusters + curve" negative-binomial model and for separating boundary steps from within-segment gradients.
- **Discussion (geometry limits).** "extend the definition of a layered tissue to account for more complex tissue geometries such as ... 3-D tissues". The authors list complex 3D geometry as unsolved. Kidney nephrons are such a case, which positions our method as covering where layered models stop.

---

## [29] Miao, Z. et al. Dissecting the coordinated progression of cell states in spatial transcriptomics with CoPro. Preprint at *bioRxiv* https://doi.org/10.64898/2026.04.17.719309 (2026)

Verified from the XML: bioRxiv, elocation 2026.04.17.719309, epub 2026-04-21; 11 authors (Miao Z … Zhang NR); PMID 42079111.

**1. Access.** Europe PMC fullTextXML (PMC13131544), the bioRxiv preprint v1. It is not peer reviewed. Read: Introduction, all Results (overview, simulations, colon ×3, healthy liver, aging liver, kidney), Discussion, and the full Methods (input data, spatial kernel, skrCCA, multi-sample framework, bandwidth, significance testing, projection, comparison methods GASTON/SpaceFlow/ONTraC, and the kidney seqFISH processing, supervised analysis and transfer).

**2. Our citing sentence.** Intro: "...and CoPro takes the order of nephron segments as an input[29]."

**3. Verdict: PARTLY SUPPORTED.**

What CoPro does (task point c):
- **Tissue coordinates: required in every mode.** The input is expression plus spatial coordinates per annotated cell type. A Gaussian kernel on Euclidean cell–cell distances is the core of its objective. Its significance test shuffles spatial patches while holding coordinates fixed.
- **Cell-type labels: required in every mode.** Each cell or pixel must first be assigned to a broad cell type.
- **Nephron segment order: an input in the kidney analysis only.** That analysis ran CoPro in supervised mode on seqFISH data: 1,298 genes, three control mice, 63,740 tubular and 18,185 vascular cells. Eleven tubular cell types were given ordinal labels (PTS1 = 1 ... PC = IC = 7). The labels were regressed onto 40 PCs by OLS to fix the tubular axis, and CoPro then found the vascular program that co-varies with it through the spatial kernel. Recovery of the known order was scored by Kendall τ = 0.74–0.78. The target was tubule–vasculature co-progression, not gradients within the PT.
- Colon, brain and liver were run unsupervised, with no order input. The liver run used a single-cell-type mode on 8-µm imputed super-pixels and still used pixel locations.

So the sentence is accurate for CoPro's kidney use. As written, though, it suggests a general requirement and omits the stronger, mode-independent differences (x/y and cell-type labels).

Quotes:
- Methods > Input data: "Spatial transcriptomic data were represented as a normalized gene-expression matrix and spatial coordinates for each annotated cell type."
- Methods > Supervised CoPro analysis: "To define the reference tubular axis, we assigned ordinal labels reflecting the known anatomical ordering of nephron segments: PTS1 = 1, PTS2 = 2, PTS3 = 3"
- Results > Overview of CoPro framework: "In the unsupervised setting, it identifies axes of spatial coordination de novo, using only the expression matrix and the spatial position of cells."

Proposed replacement (keeps the sentence structure):
> "Methods that model continuous expression gradients, such as GASTON and Belayer, build their coordinate from tissue x/y positions[27,28], and CoPro, which also uses tissue coordinates and cell-type labels, was given the known order of nephron segments as an input when applied to the kidney[29]."

Shorter alternative:
> "Methods that model continuous expression gradients build their coordinate from tissue x/y positions[27–29]; applied to the kidney, CoPro was also given the order of nephron segments as an input[29]."

Novelty note for the Discussion. CoPro also compares samples on a shared axis without spatial registration ("projecting all cells from different samples onto a biologically defined axis derived from shared gene expression programs"). Our cross-group comparison on a shared t should not be presented as the first registration-free cross-sample axis. What differs is that ours is built from expression alone, without x/y, labels or a reference program.

**4. Writing ideas.**
- **Introduction (discretization).** "it may impose artificial categorization on what might be continuous spatial transitions" (Introduction, on spatial-domain methods). This is a measured way to state our segment-label critique.
- **Introduction (comparative methods; sets up STalign/STcompare).** "many rely on spatial registration, assuming that tissue regions with similar histological appearance correspond to biologically comparable states". This states the assumption that fails for nephron cross-sections.
- **Species / AKI comparison and Discussion (shared axis, its cost).** "the transfer of gene weights from a reference condition to additional samples defines a common biological axis across samples without relying on spatial registration". The honest limit comes later: "the approach may obscure tissue-specific patterns that exist in individual samples but deviate from the reference program". Our joint axis across species and injury carries the same risk, and a Discussion sentence modelled on this would pre-empt the reviewer.
- **Methods (null calibration).** "This patch-based design preserves within-cell-type spatial autocorrelation while breaking cross-cell-type spatial coupling, thereby generating an appropriate null distribution". This is a template for saying what our balanced specimen relabeling preserves (within-section structure) and what it breaks (the group label).

---

## [30] Clifton, K. et al. STalign: alignment of spatial transcriptomics data using diffeomorphic metric mapping. *Nat. Commun.* 14, 8123 (2023)

Verified from the XML: vol. 14, article 8123, epub 2023-12-08; 9 authors (Clifton K … Fan J); PMC10709594.

**1. Access.** Europe PMC OA fullTextXML (PMC10709594). Read: Abstract, Introduction, all Results (overview; within- and across-technology alignment; 3D CCF; diverse tissues), Discussion, Methods (datasets, application, evaluation metrics, statistics, implementation) and figure legends.

**2. Our citing sentence.** Intro: "Comparative methods align sections to compare expression at matched locations[30,31], which requires sections that share structure in the tissue plane."

**3. Verdict: SUPPORTED.** STalign aligns sections through cell-density images (LDDMM), so that expression can be compared at matched locations. It relies on visible shared structure, and it can handle sections that only partly match.
- Introduction: "Comparative analysis of ST datasets at matched spatial locations across tissues, individuals, and samples provides the opportunity to interrogate spatial gene expression and cell-type compositional variation"
- Results > diverse tissues: "STalign relies on variations in cell densities that generally form visible structures that can be used for alignment."
- Discussion: "alignment accuracy at the resolution of single cells is limited by the fact that there is generally no one-to-one correspondence between cells across samples, particularly for complex tissues." (The same holds for tubule cross-sections: no cross-section in one kidney corresponds to one in another.)

**4. Writing ideas.**
- **Abstract / Introduction first sentence (overlap to be aware of).** STalign's abstract opens "Spatial transcriptomics (ST) technologies enable high throughput gene expression characterization within thin tissue sections." Our abstract and Introduction open with almost the same sentence. That is reasonable lab house style, but the authors should know the overlap is near-verbatim. The next sentence, "However, comparing spatial observations across sections, samples, and technologies remains challenging", has the same "However, ..." shape as our abstract.
- **Introduction (tool entry).** "To address these challenges, we present an approach called STalign that builds on recent developments in Large Deformation Diffeomorphic Metric Mapping". Our "To address this, we develop pseudostructure reconstruction, which ..." already follows this. Their version names the mathematical basis in the same sentence, which suits our "elastic principal curve".
- **Results overview / validation.** "We show that STalign achieves high gene expression and cell-type correspondence across matched spatial locations that is significantly improved over landmark-based affine alignments." Claim, metric and comparator in one sentence, as we need for "orders segments as well as the best tree-based trajectory method".
- **Discussion (positioning and limits).** "STalign complements gene-expression-based approaches for sample alignment by focusing on the real space rather than a higher-order transcriptomic manifold." Our method is the reverse case, so this sentence lets us frame it as complementary to the lab's own alignment tool. Their limits paragraph ("As LDDMM relies on optimization using gradient descent, the resulting alignment solution may converge on local minima") is a model for stating that the elastic curve fit depends on initialization.

---

## [31] Clifton, K. et al. STcompare: comparative spatial transcriptomics data analysis of structurally matched tissues to characterize differentially spatially patterned genes. *Bioinformatics* 42, btag644 (2026)

Verified from the XML: vol. 42, issue 9, btag644, published 2026-08-26; 7 authors (Clifton K, Jiang V, dos Santos Peixoto R, Singh S, Matsuura R, Rabb H, Fan J; the eighth contributor in the XML is the Associate Editor). "Clifton, K. et al." is correct.

**1. Access.** Europe PMC OA fullTextXML (PMC13577157), the published version. Read: Abstract, Introduction, Methods, all Results (simulations, MERFISH brain replicates, kidney AKI), Discussion and the Figure 1–2 legends. The Supplementary Methods were not read.

**2. Our citing sentence.** As for [30].

**3. Verdict: SUPPORTED.**
- Methods: "To use STcompare, ST datasets should be preprocessed to ensure matched spatial locations for comparison such as through alignment of shared tissue structures"
- Discussion: "Both the spatial correlation test and spatial fold-change test of STcompare are dependent on the datasets under comparison having inherent structural similarity."
- Results 3.3 (kidney): "The two tissue sections exhibit strong structural correspondence." The comparison was one 24-h ischemia–reperfusion Visium section against one sham section, aligned by STalign on cortex / outer medulla / inner medulla labels and rasterized to 200-µm pixels.

Note for our AKI section. STcompare is the closest kidney AKI precedent, and it works at the scale of cortex and medulla. Its PT findings are compartment-level: of Acsm2 and Slc34a1 it reports that "these cell-type specific genes remained localized to the cortex, but expression was greatly reduced". A clause saying that matched-location comparisons resolve the kidney at the scale of its compartments, while we compare positions along the PT, would sharpen what is new without criticising the lab's own tool.

**4. Writing ideas.**
- **Introduction (gap stated as a capability).** "Existing methods generally fail to distinguish changes in spatial patterning by focusing only on changes in gene expression magnitude" (Abstract), then "To address these limitations, we develop STcompare". This is the lab's current gap → tool sentence pair.
- **Methods / Results (calibration).** STcompare first shows the naive test failing ("We observed a high FPR of 0.50 ... well beyond the expected FPR of 0.05"), then introduces its fix ("STcompare generates empirical null distributions that model spatial autocorrelation"). We could show the same two steps for group comparisons: an uncalibrated per-tubule test, then the balanced specimen relabeling.
- **Validation (replicates as the negative control).** "Because these are biological replicates, we hypothesized that genes should be significantly positively correlated with high similarity". This is a template for stating what the control-versus-control comparison should show before showing it.
- **Discussion (the caveat that maps onto ours).** "when structural alignment is poor, the interpretation of STcompare results may become challenging, as any statistically significant differences could reflect underlying alignment error". For us, a shift in t between groups could reflect axis error. The held-out-gene axes and the relabeling null are our answer, and a parallel sentence would make that explicit.

---

## [32] Trapnell, C. et al. The dynamics and regulators of cell fate decisions are revealed by pseudotemporal ordering of single cells. *Nat. Biotechnol.* 32, 381–386 (2014).

**1. Access.**
- Read: Europe PMC `fullTextXML` for PMC4122333, the NIH author manuscript NIHMS570726. Its title is the manuscript title, "Pseudo-temporal ordering of individual cells reveals dynamics and regulators of cell fate decisions". It includes the main text and figure legends.
- Also read: the Supplementary Figures PDF (NIHMS570726-supplement-1.pdf, 17 pages of figures and legends), downloaded from the PMC Open Data S3 bucket (`pmc-oa-opendata/PMC4122333.1/`). The Europe PMC supplementary-files endpoint returned only one image, and the pmc.ncbi.nlm.nih.gov link served a proof-of-work page.
- Also read: the published abstract, from Europe PMC.
- **NOT read:** the Online Methods and the "Supplemental Methods". They are not in the author manuscript, and nature.com shows "This is a preview of subscription content". How the ordering genes and the path direction were chosen is therefore unverified from the primary source.

**2. Citing sentence (Introduction).** "Trajectory inference orders cells along continuous processes from expression alone[32–36], and related ideas have ordered repeated structures: …"

**3. Verdict: SUPPORTED** for "orders cells along continuous processes from expression". "Alone" is covered by Fix C.
- Main text: "we developed a novel unsupervised algorithm, Monocle, which re-ordered the cells to maximize the transcriptional similarity between successive pairs (Fig 2a)."
- Main text: "Importantly, Monocle is unsupervised and needs no prior knowledge of specific genes that distinguish cell fates, and is thus suitable for studying a wide array of dynamic biological processes."
- Caveat, from a secondary source (Slingshot, Background, FULL TEXT): "The direction of this path and the number of branching events are left to the user". Slingshot's Table 1 lists Monocle's supervision as "Path direction, number of lineages". The published abstract also says the data were "collected at multiple time points".

**4. Writing ideas.**
- "In essence, a single-cell RNA-Seq experiment might constitute a time-series, with each cell representing a distinct time point along a continuum." Serves the **Introduction**: the direct analogue of "each cross-section is a sample of one shared structure".
- "Drawing conclusions from a group of individuals based on the properties of their average is a hazardous practice because the average can mask important trends among the individuals". Serves the **Introduction**, at the point about contrasting segment means. A further phrase serves the **Results** (continuous vs clusters): "simply comparing gene expression levels between groups of cells collected on different days masked changes in key transcriptional regulators".
- "Monocle placed subsets as small as 50 cells in pseudo-temporal order highly similar (spearman >= 0.8) to their relative order within the full data set." Serves **validation**: a template for reporting our tubule-downsampling stability as a correlation with the full-data order.
- "genes that play active roles at the early and late stages … showed pseudo-temporal kinetics that were highly consistent with expectations … confirming the accuracy of the ordering". Serves the **Results overview / validation**: known markers confirm the order, as our segment-order recapitulation does.

---

## [33] Haghverdi, L., Büttner, M., Wolf, F. A., Buettner, F. & Theis, F. J. Diffusion pseudotime robustly reconstructs lineage branching. *Nat. Methods* 13, 845–848 (2016).

**1. Access.**
- Read: the **bioRxiv preprint v2** (doi 10.1101/041384, posted 2016-03-02), as JATS XML from `biorxiv.org/content/early/2016/03/02/041384.source.xml`. It is complete, with main text, Box 1, figure legends and Online Methods.
- **The published Brief Communication was NOT read.** It has no PMCID. OpenAlex lists an MDC repository accepted-version PDF, which returned HTTP 401. The Semantic Scholar `openAccessPdf` query returned HTTP 429 three times. The only published text read is the abstract, from Europe PMC, and it is worded differently from the preprint's, so other published wording may differ too.

**2. Citing sentences.**
- Introduction: "Trajectory inference orders cells along continuous processes from expression alone[32–36], …"
- Methods: "The same corrected embedding and root were ordered by diffusion pseudotime[33] (Scanpy 1.11.5; max(30, √n) neighbours, 15 diffusion components), Slingshot 2.14.0[34] on the pass-2 clusters (R 4.4.3), scFates 1.2.5 (principal tree)[ref Faure 2023] and GASTON[27] …"

**3. Verdicts.**
- **Introduction: PARTLY SUPPORTED.** The order is computed from expression, but only relative to a user-given root cell.
  - Published abstract: "The temporal order of differentiating cells is intrinsically encoded in their single-cell expression profiles."
  - Main text (preprint): "Fixing a known root cell r as start of the dynamical process of interest, we define the pseudotime of cell x as dpt(x,r)."
  - Online Methods, "Overview of DPT algorithm" (preprint): the inputs are "The n by G data matrix" and "One (or several) root cell(s)."
  - Fix: see Fix C (top of this report).
- **Methods: SUPPORTED.** Our use of a single root unit matches DPT's definition of pseudotime as the distance from a root cell, and the preprint allows several roots, whose pseudotimes are averaged. One precision issue remains: the paper's DPT uses the full transition matrix, while Scanpy's `dpt` works on the top `n_dcs` diffusion components. Our "15 diffusion components" is a Scanpy setting, not part of the method in [33].
  - Main text (preprint): "Consequently, DPT is based on the full rank T rather than a low rank approximation."
  - Proposed wording: "…were ordered by diffusion pseudotime[33] as implemented in Scanpy 1.11.5[ref Wolf 2018] (max(30, √n) neighbours; 15 diffusion components; rooted at the same root unit), …"

**4. Writing ideas.**
- "we counted how often a pseudotime puts a cell from a later temporal sorting before an earlier one (measured by Kendall rank correlation τ)." Serves **validation / comparison with other methods**: a one-sentence model for defining our segment-order ρ in words.
- Box 1: "Pseudotimes are proxies for universal time". Serves the **Discussion / Methods**: a template for saying that t is a proxy for position along the structure, not a measured distance.
- "Even in the absence of branching, DPT is significantly more robust with respect to noise in low-density regions and cell outliers than existing methods". Serves the **comparison with other methods**: how to phrase a robustness claim for a single-path case.
- Fig. 2f legend: "Comparisons of Kendall rank correlation on bootstrap samples (n=100 bootstrap runs, downsampling to maximal 1800 cells…)". Serves **validation**: report method-versus-method stability on repeated subsamples.

---

## [34] Street, K. et al. Slingshot: cell lineage and pseudotime inference for single-cell transcriptomics. *BMC Genomics* 19, 477 (2018).

**1. Access.** Europe PMC open-access `fullTextXML` for PMC6007078. I read the Background, Results and discussion, Conclusions, Methods, the figure legends, and Table 1 (parsed from its XML rows). Additional file 1 was not read.

**2. Citing sentences.**
- Introduction: "Trajectory inference orders cells along continuous processes from expression alone[32–36], …"
- Methods: "… Slingshot 2.14.0[34] on the pass-2 clusters (R 4.4.3) …"

**3. Verdicts.**
- **Introduction: PARTLY SUPPORTED.** Slingshot needs a clustering and a start cluster. It also states that every method needs a user-identified starting state.
  - Methods: "Like other methods (TSCAN, Waterfall, Monocle 2), Slingshot requires the user to identify the initial cluster or root node."
  - Methods: "all existing lineage inference methods explicitly or implicitly make the assumption that a starting state can be identified by the user".
  - Fix: see Fix C (top of this report).
- **Methods: SUPPORTED.** Running Slingshot on clusters is exactly its design: "Slingshot assumes that the n cells have been partitioned into K clusters, potentially corresponding to distinct cellular states." (Methods). Two implementation details are worth stating, both taken from `results/paper_figures/fig3/baselines/methods_notes.md`:
  - The run used `start.clus` = the root unit's pass-2 cluster.
  - The run used `slingAvgPseudotime`, which averages the per-lineage pseudotimes with the cells' lineage weights. This averaging is a package function and is not described in the 2018 paper; the paper describes per-lineage pseudotimes and cell weights.
  - The pipeline produced 5–6 lineages on our data.
  - Proposed wording: "Slingshot 2.14.0[34] on the pass-2 Leiden clusters, starting from the root unit's cluster, with each unit's per-lineage pseudotimes averaged by its lineage weights (R 4.4.3)". This also corrects "the same … root": Slingshot was given the root's cluster, not the root unit.

**4. Writing ideas.**
- Results: "principal curves give very robust pseudotimes when there is a single lineage." Results: "The principal curve approach was the most stable method, but on more complex datasets, it has the obvious limitation of only characterizing a single lineage." These serve the **Results / Discussion**: they support fitting one unbranched curve, and they phrase the matching limitation honestly.
- Methods: "Project all data points onto the curve and calculate the arc length from the beginning of the curve to each point's projection." Serves **Methods**: standard wording for our position t, already close to our draft.
- Fig. 2: "how consistently it orders cells over 50 subsamples (bootstrap samples with duplicates removed)". Serves **validation**: a subsampling-stability design and how to state it.
- Background: "because different lineage inference methods make different upstream choices, they can be difficult to compare." Serves the **comparison with other methods**: it justifies running every method on the same corrected embedding and root.

---

## [35] Wolf, F. A. et al. PAGA: graph abstraction reconciles clustering with trajectory inference through a topology preserving map of single cells. *Genome Biol.* 20, 59 (2019).

**1. Access.** Europe PMC open-access `fullTextXML` for PMC6425583. I read the Background, all Results subsections, Conclusions, Methods and the figure legends. Additional file 1 (the Notes) was not read.

**2. Citing sentence (Introduction).** "Trajectory inference orders cells along continuous processes from expression alone[32–36], …"

**3. Verdict: PARTLY SUPPORTED.** PAGA does order cells along a continuous variable, but on top of a partition, and by distance from a root cell.
- Background: "inferring pseudotemporal orderings or trajectories of cells [2–4] assumes that data lie on a connected manifold and labels cells with a continuous variable—the distance along the manifold."
- Fig. 1 legend: "we order cells within each partition according to their distance from a root cell."
- Methods: "PAGA works with any underlying clustering algorithm or experimentally generated groupings of observations." Pseudotime is "an extended version of diffusion pseudotime (DPT) … that accounts for disconnected graphs."
- Fix: see Fix C (top of this report).

**4. Writing ideas.**
- Background: "Clustering assumes that data is composed of biologically distinct groups such as discrete cell types or states and labels these with a discrete variable—the cluster index." Serves the **Introduction**: the cleanest one-line framing of "discretize" against "continuous".
- Results: "graph clustering alone is, as any clustering method, an ill-posed problem in the sense that many highly degenerate quasi-optimal clusterings exist". Serves the **Introduction / Discussion**: why segment labels from clustering are not a ground truth for boundaries.
- Results: "Sampling widely varying parameters, which leads to widely varying clusterings, we find that the inferred abstraction … is much more robust than the underlying graph clustering algorithm". Serves **validation**: phrasing for robustness to upstream parameters, such as our half-gene axes.
- Results: "Even though the connections in PAGA graphs often correspond to actual biological trajectories, this is not always the case." Serves the **Discussion**: the model for saying that an expression-derived order is not anatomical truth.

---

## [36] Saelens, W., Cannoodt, R., Todorov, H. & Saeys, Y. A comparison of single-cell trajectory inference methods. *Nat. Biotechnol.* 37, 547–554 (2019).

**1. Access.**
- Read: the **authors' manuscript source**, `manuscript/paper.Rmd` in the GitHub repository dynverse/dynbenchmark (raw.githubusercontent.com). Its header says "Published in Nat Biotechnol 37, 547–554 (2019)". Its text commits run 2018-11-13 to 2019-01-09, and a render update followed in 2021. I read the Abstract, Introduction, Results (accuracy, scalability, stability, usability), Discussion and Online Methods (method input and prior information, datasets, metrics).
- Caveats on this source:
  - It is not the typeset article.
  - Its abstract is worded differently from the published one, so other wording may also differ.
  - Its numbers are R inline code (`` `r n_methods_evaluated` ``), so it states no counts. The counts come from the published abstract (Europe PMC): ">70" tools, "45 of these methods on 110 real and 229 synthetic datasets".
- The bioRxiv preprint v1 (10.1101/276907, 2018-03-05, titled "A comparison of single-cell trajectory inference methods: towards more accurate and robust tools"; its counts may differ from the published ones) could not be fetched. Both the JATS and the PDF returned Cloudflare error 1015 on four tries. The Ghent University repository copy (hdl 1854/LU-8619545) is "UGent only".

**2. Citing sentence (Introduction).** "Trajectory inference orders cells along continuous processes from expression alone[32–36], …"

**3. Verdict: PARTLY SUPPORTED.** The benchmark supports "orders cells along a trajectory by expression", but it says explicitly that many methods need prior information.
- Introduction: "trajectory inference (TI) methods, also called pseudotime analysis, which order cells along a trajectory by their overall similarity in -omics datasets".
- Online Methods, Method input: "A large portion of the methods requires some form of prior information (e.g. a start cell) in order to be executable."
- Online Methods, Method input: "Prior information can be supplied as a starting cell from which the trajectory will originate, a set of important marker genes, or even a grouping of cells into cell states."
- Fix: see Fix C (top of this report).

**4. Writing ideas.**
- Results, Accuracy: "Real datasets were classified as "gold standard" if the reference trajectory was not extracted from the expression data itself, such as via cellular sorting or cell mixing". Serves **validation**. It gives a vocabulary for saying honestly that our segment labels are an expression-derived ("silver") reference, and that DPT and scFates are comparators, not truth.
- Discussion: "it is critical that a trajectory, and the downstream results and/or hypotheses originating from it, are confirmed by multiple TI methods." Serves the **comparison with other methods**: the rationale for reporting results that hold on both our curve and the DPT and scFates baselines.
- Results, Accuracy: "methods tend to be either too optimistic or too pessimistic regarding the complexity of the topology in the data." Serves the **Results (comparison with other methods)**: context for the tree methods that branched on our unbranched structure (Slingshot gave 5–6 lineages; scFates gave 4–7 tips, per methods_notes).
- Results, Stability: "we executed each method on 10 different subsamples of the datasets (95% of the cells, 95% of the features), and calculated the average similarity between each pair of models". Serves **validation**: a compact template for our downsampling and half-gene stability tests. A further phrase serves the **Discussion**, on choosing the root from a marker axis: "incorrect or noisy prior information can bias the trajectory towards current knowledge".

---

## [37] Niec, R. E. et al. Lymphatics act as a signaling hub to regulate intestinal stem cell activity. *Cell Stem Cell* 29, 1067–1082.e18 (2022)

Verified from the XML: vol. 29, issue 7, pp. 1067–1082.e18 (ppub 2022-07-07); 15 authors (Niec RE … Fuchs E).

**1. Access.** NIH author manuscript NIHMS1809914 via NCBI efetch (PMC9271639). Europe PMC returned HTTP 500. Read: Introduction, Results (all sections, especially "The in vivo spatial maps of cell types" and "Cell type and expression cartography"), Discussion (including "The Power and Universality of Our Computational Strategy" and Limitations), the STAR Methods sections on Visium, BayesPrism and SpaceFold cartography, and the Figure 5 legend.

**2. Our citing sentence.** Intro: "...and related ideas have ordered repeated structures: intestinal crypt positions from deconvolved cell-type composition[37], ..." plus "These approaches rely on a change in cell-type composition, ...".

**3. Verdict: SUPPORTED.** SpaceFold describes each Visium spot by its BayesPrism cell-type fractions (about 30 types, from a matched scRNA-seq reference). It z-scores the fractions and projects them onto one dimension with PHATE (knn = 10, ndim = 1, α = 40), pooling hundreds of crypts (about 4 spots per crypt). No x/y is used, and the stated assumption is that composition changes gradually along the crypt.
- STAR Methods > Projecting spots onto a 1D spatial axis: "we assume that two spots with similar cell type fractions are in a similar spatial position along the crypt, and moreover these cell type fractions change gradually along the crypt axis"
- STAR Methods: "We used the cell-type fraction inferred by BayesPrism as input, and used PHATE (Moon et al., 2019) to project this vector onto a 1D axis."
- Results: "We hypothesized that, as spatial spots were randomly sampled along the stereotypical axis, the cell type composition of each spatial spot ... would inform the relative physical coordinates of cells"

Precision (optional). SpaceFold orders spots sampled along crypts, not crypts as wholes; "intestinal crypt positions" is a fair shorthand. If a reviewer could read "ordered repeated structures" literally, use "Visium spots along intestinal crypts from deconvolved cell-type composition[37]".

Validation was qualitative: expected cell-type order, villus zone markers (Moor 2018) and SOX9 immunofluorescence. There is no quantitative comparison with a measured position along the crypt, so our validation design can claim to go further.

**4. Writing ideas.**
- **Introduction / Results overview (the repeated-unit premise).** "SpaceFold took advantage of the hundreds of highly stereotyped structures of the intestinal epithelium, each comprised of repeating crypt-villus units". The Discussion adds "pooling information across many crypts to elucidate robust shared biological features". These are the closest published wording to our "treats every cross-section of a repeated structure as a sample of one shared structure".
- **Discussion (identical-unit assumption).** "Intestinal tissue is known to show heterogeneous gene expression across disparate tissue segments ... which may violate the assumption of identical repetitive units." SpaceFold handled this by running regions separately. We should state the parallel risk (superficial vs juxtamedullary nephrons; injured vs control) and how a shared axis is checked.
- **Discussion (generality claim).** "it can be easily generalized to other tissue types of repetitive structures" (STAR Methods), with "e.g. liver lobules" in the Discussion. This is a precedent for our generality claim (crypts, lobules), with the difference that ours needs no deconvolution reference.
- **Comparison with other methods.** Their result statement, "1D pseudo-space inferred by SpaceFold recapitulated the expected cell type distributions along SI crypt-villus and LI crypt axes", shows that the precedent was validated only qualitatively. It sets the bar our ρ-against-segment-order validation clears.

---

## [38] Karin, J., Bornfeld, Y. & Nitzan, M. scPrisma infers, filters and enhances topological signals in single-cell data using spectral template matching. *Nat. Biotechnol.* 41, 1645–1654 (2023)

Verified from the XML: vol. 41, issue 11, pp. 1645–1654; epub 2023-02-27; 3 authors as cited.

**1. Access.** Europe PMC OA fullTextXML (PMC10635821). Read: the introduction (unheaded "Main"), all Results (template matching, HeLa, liver, Chlamydomonas, SCN, generalized templates), Discussion, and the Methods (cyclic and linear spectra, preprocessing, general-case algorithm, reconstruction problems 1–2, gene inference, filtering and enhancement).

**2. Our citing sentence.** Intro: "...testis tubules around the known cycle of the seminiferous epithelium[38,39], ..." plus "These approaches rely on ... a cyclic topology ...".

**3. Verdict: PARTLY SUPPORTED.** scPrisma analyses no testis or seminiferous data. Its datasets are HeLa, liver hepatocytes, Chlamydomonas, SCN, pancreas, hESC and a Slide-seqV2 hippocampus. It is the method Chakravorty et al. [39] used ("de-novo cyclic pipeline in scPrisma"). The topology is a template the user chooses, and linear templates are derived alongside cyclic ones. In the main text, the reconstruction (ordering) step is shown for cyclic signals, while linear templates are used for enhancement and filtering.
- Methods > scPrisma general-case algorithm: "Choose the desired topology (for example, periodic/linear). Calculate the theoretical covariance eigenvectors and eigenvalues."
- Results > Spectral template matching: "Circulant matrices have closed-form formula for their eigenvectors and eigenvalues ... which we used to estimate the ordering of cells along the cyclic topology."
- Introduction: "it enables integration of diverse types of prior knowledge (such as low-resolution temporal ordering) but does not rely on it and can be used for de novo analyses."

Proposed replacement for the clause:
> "testis tubules around the known cycle of the seminiferous epithelium[39], using a cyclic topological template[38]"

or keep "[38,39]" but say "with a cyclic topological prior".

**4. Writing ideas.**
- **Introduction / Discussion (global shape, not markers).** "Relying on such global topological features can provide a more robust and generalizable approach than relying on specific features such as marker genes". It argues for ordering by global shape (one curve) rather than marker panels, which is our case for not using segment labels.
- **Methods (topology as a declared prior).** "Choose the desired topology (for example, periodic/linear)" is the first line of their algorithm. Our open, unbranched curve is the same kind of choice and should be stated as a prior, not hidden.
- **Validation.** "evaluating results obtained by scPrisma can be done by comparison to partial prior knowledge, such as marker genes known to be related (or unrelated) to recovered biological signals". This frames our evaluation against segment labels and known S1/S3 markers as comparison with partial prior knowledge, not ground truth.
- **Discussion (overfitting to the template).** "scPrisma does not converge to a topology that is not a reflection of a strong-enough signal in the data". This is a model for the claim we need to make and test: the curve does not impose an order where expression carries none (our downsampling and half-gene stability).

---

## [39] Chakravorty, A. et al. The temporal architecture of the seminiferous epithelial cycle revealed by spatial transcriptomics. *Cell* 189, 4122–4139.e24 (2026)

Verified from the XML: vol. 189, issue 13, pp. 4122–4139.e24; epub 2026-05-15; 8 authors (Chakravorty A … Cai L).

**1. Access.** NIH author manuscript via NCBI efetch (PMC13326559). Europe PMC returned HTTP 500. Read: Introduction, Results (all sections, especially "Tubule-level analysis reveals the temporal architecture"), Discussion and Limitations, the STAR Methods (tubule segmentation, preprocessing, "Determining order of tubules and Tubule Angle", mapping to stages, population dynamics, modified Moran's I, synchronization), and the legends of Figure 2 and Figure S3.

**2. Our citing sentence.** As for [38].

**3. Verdict: SUPPORTED.**

What is ordered, by what, with what prior and how it is validated (task point b):
- **What is ordered.** Whole seminiferous-tubule cross-sections, hand-segmented in FIJI. The untreated set is 337 cross-sections from two mice, plus 317 from a third; perturbed sets were ordered separately.
- **By what.** seqFISH+ transcripts (2,653 genes) summed per tubule, normalized and log-scaled, then ordered by scPrisma's de novo cyclic reconstruction. The order was converted to a "Tubule Angle" spread evenly over 0–360°.
- **Prior.** A closed cycle (stages I–XII). The even spacing of angles is justified by showing that stages are spatially random in the section (cyclic Moran's I). Tissue x/y enters only in that check, not in the ordering.
- **Validation.** Sequential germ-cell composition along the angle; mapping to stages from published stage frequencies (3,000 tubules from 12 mice); independent DAPI nuclear morphology; and a simulation showing that 337 samples recover the dynamics. Repeat runs checked that the order is stable. There is no comparison with a measured time.

Quotes:
- STAR Methods > Determining order of tubules: "The order of tubules in both Untreated and Busulfan conditions was determined using the de-novo cyclic pipeline in scPrisma. Prior to running the pipeline, counts were normalized and log-scaled."
- Results: "After normalization and scaling (Methods), tubules organized into a circular transcriptional topology in principal component (PC) space (Figure 2C)."
- Results: "Using published stage durations, we mapped Tubule Angles to classical stages I–XII (Figures S3J,S3N). DAPI nuclear morphology independently validated these assignments"

**4. Writing ideas.**
- **Introduction / Results overview (each cross-section is a sample).** "Since each tubule cross-section represents a time point in the seminiferous epithelial cycle, we resolved tubule ordering to reconstruct the cycle's temporal architecture". This is the closest published parallel to our premise, with position replacing time.
- **Results (sampling argument and its kidney caveat).** "each tubule cross-section represents a random time point along the seminiferous epithelial cycle, enabling unbiased sampling". We cannot claim uniform sampling: a section samples nephron positions roughly in proportion to how much tubule length lies in the plane (convoluted S1/S2 in cortex versus straight S3). Our t is therefore a rank-like position, not physical length. Saying this in Results or the Discussion pre-empts the obvious question.
- **Methods (robustness).** "To ensure robustness of the cyclic ordering, we repeated the pipeline multiple times for each dataset and compared the resulting cyclic orders for consistency." A plain one-line robustness statement, a model for introducing our downsampling and half-gene stability.
- **Discussion (2D limit).** "our 2D analysis may underrepresent cell-cell interactions along the tubule's longitudinal axis; three-dimensional reconstruction would provide more spatial context." This supports our framing that a 2D snapshot of repeated tubules can stand in for 3D tracing, with an honest limit attached.

---

## [40] Asami, S., Yin, C., Fan, J., Garza, L. A. & Kalhor, R. Four-dimensional molecular mapping from a spatial snapshot reveals the dynamics of hair follicle organogenesis. *Cell* 189, 5653–5672.e14 (2026)

Verified from the XML: vol. 189, issue 18, pp. 5653–5672.e14; epub 2026-07-01, ppub 2026-09-03; 5 authors as cited (Kalhor R senior; Fan J third); doi 10.1016/j.cell.2026.06.014.

**1. Access.** NIH author manuscript NIHMS2193327 via NCBI efetch (PMC13379681). Europe PMC returned HTTP 500. The earlier notes held only the abstract; this is the first full-text read. Read: Introduction; Results ("3DEEP-HybISS captures 3D hair follicles", "Spatial snapshot animates hair follicle development", "Virtual 3D dissection", the nude-mouse comparison); the full Discussion including "Spatial pseudotiming provides a new paradigm" and Limitations; and the STAR Methods on follicle isolation and linearization, pseudotime inference, secondary/tertiary classification, squared-error analysis, and KO-vs-WT maturation, stability, composition and differential expression.

**2. Our citing sentence.** Intro: "...and hair follicles by developmental age from three-dimensional profiles[40]." plus "These approaches rely on ... three-dimensional data."

**3. Verdict: SUPPORTED.**

What is ordered, by what, with what prior and how it is validated (task point b):
- **What is ordered.** Whole hair follicles imaged in 3D in sections up to 400 µm thick (3DEEP-HybISS). There were 371 complete follicles in the 85-gene wild-type sample, 261 wild-type follicles in the 228-gene samples, and 262 Foxn1 nude follicles. Follicles were traced by hand every 10th z-plane and interpolated. Partial follicles at field edges were excluded "due to incomplete capture or uncertain 3D organization".
- **By what.** Transcripts were summed per follicle into a follicle-by-gene count table and normalized (scTransform), then reduced by PCA. Slingshot gave a 0–1 pseudotime for the 85-gene data. For the 228-gene data, PC1 of a joint WT + KO PCA was used as the pseudotime proxy.
- **Prior.** A single, unbranched developmental trajectory: follicles are induced asynchronously and grow continuously, so one snapshot spans the trajectory.
- **Validation.** Spearman correlation with follicle **length**, a physical proxy for stage (ρ = 0.94 and 0.90; 0.78 in KO). Separating secondary from tertiary follicles halved the error of the length–pseudotime fit. The two independent panels agreed, and follicle structures matched immunohistochemistry and smFISH.
- **Comparison across genotypes.** One axis was computed jointly for WT and KO. KO follicles were 8% behind WT in pseudotime at equal length, and the residual variance was 30% higher. Mixed models treated each animal as a group to avoid pseudoreplication.

Quotes:
- Results > Spatial snapshot animates hair follicle development: "We thus used Slingshot to order follicles along a developmental trajectory (Figure 2C), creating a whole-organ pseudotime ranging from 0 to 1."
- Results: "The strong correlation between follicle pseudotime and length (Spearman's ρ = 0.94, p < 6e-178) indicates that pseudotime accurately orders follicles by developmental stage"
- Results: "ordering the 3D transcriptomes of follicles by developmental stage—i.e., spatial pseudotiming—transformed our spatial snapshot into a 4D molecular map of organogenesis"

**Caution for our Introduction and Discussion.** Asami's Discussion names the kidney as a target: "This approach should be broadly applicable to systems where similar structures develop asynchronously, including kidney nephrons, islets of Langerhans, and lung alveoli." It means developmental time across whole structures, not position along one mature structure. A reviewer who knows the paper (Jean Fan is a co-author) may still raise it. One Discussion sentence would settle it: "Spatial pseudotiming orders whole structures that develop asynchronously by developmental age[40]; pseudostructure reconstruction instead orders cross-sections of mature structures by position along them, from 2D sections in which no structure is captured whole." Asami's own Introduction supports the 2D point: thin 2D sections "perform best for organs large enough for serial sectioning (e.g., the heart) or those repetitive along at least one axis (e.g., a nerve cord)".

**4. Writing ideas.**
- **Discussion (naming the idea and an analogy).** "A key conceptual advance of this work is the recognition that asynchronously developing cell groups within a single sample can be computationally ordered to capture developmental time". The cryo-EM analogy follows: "spatial pseudotiming is to serial sampling what cryogenic electron microscopy (Cryo-EM) is to X-ray crystallography". We need our own one-line concept statement. A parallel analogy (reconstruction from many randomly oriented copies) would fit our "every cross-section is a sample of one shared structure" and is strongest if it cites the precedent rather than echoing it.
- **Validation.** "The strong correlation between follicle pseudotime and length ... indicates that pseudotime accurately orders follicles by developmental stage". This is the form of a validation sentence: an independent physical or anatomical proxy, a rank correlation, and the conclusion. Ours is segment order (ρ), with DPT/scFates as comparators.
- **AKI / species comparison.** "Whole-organ pseudotime, calculated jointly for WT and KO, successfully ordered KO follicles ... revealing that KO follicles follow the same developmental trajectory as WT." It is followed by the lag and instability results, and in Methods: "Hair follicles from the same sample ... were treated as a group to prevent counting them as completely independent". This is a ready model for our shared axis across control and AKI and for "Sections are the replication units".
- **Discussion (limits).** "the dynamics are based on pseudotime, which may not have a linear relationship with real time". The direct analogue for us: t need not be linear in physical length along the nephron.

---

# Cross-reference

## Ambiguous referent: "These approaches rely on a change in cell-type composition, a cyclic topology or three-dimensional data."

Each part holds for its precedent: composition for SpaceFold (assumption quoted under [37]); a cycle for Chakravorty/scPrisma; complete 3D structures for Asami (partial follicles excluded). The sentence directly follows "Trajectory inference orders cells along continuous processes from expression alone[32–36], and related ideas have ordered repeated structures...". "These approaches" can therefore be read as including trajectory inference, which relies on none of the three. Proposed:
> "These repeated-structure approaches rely on a change in cell-type composition, a known cyclic topology or three-dimensional capture of whole structures."

"Three-dimensional capture of whole structures" is more exact than "three-dimensional data". Asami's ordering itself uses only per-follicle counts; the 3D imaging is what ensures every follicle is complete.
