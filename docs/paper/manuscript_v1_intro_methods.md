# Pseudostructure reconstruction of kidney nephrons enables comparative 3D analysis of molecular gradients

Draft v1 of the Introduction and Methods (2026-11-03; Methods restructured the same day after the lab's pattern:
a short overview whose run-in paragraphs are the method's components, then gene models, handoffs and group
comparisons, one subsection per dataset, evaluation, other methods, statistics). Written against the pipeline as it now runs: H&E segmentation
with a frozen OpenMidnight encoder (no SAM2), the reconstruction package's elastic principal curve (diffusion
pseudotime is a baseline only), and negative-binomial gene models on per-gene held-out axes. Bracketed items
`[...]` are facts we do not hold in the repo and need from collaborators. Numbers come from the frozen protocols in
`docs/paper/fig3-*`, `fig4-*`, `fig5-*` and `fig6-*`.

## Abstract (agreed 2026-11-02, 145 words)

Spatial transcriptomics (ST) technologies enable gene expression characterization within tissue sections.
However, sections of repeating three-dimensional structures such as kidney nephrons capture scattered
cross-sections whose position along the structure is lost. To address this, we develop pseudostructure
reconstruction to segment such cross-sections from histology and order them along a shared axis using paired ST
data, without segment labels. We apply pseudostructure reconstruction to VisiumHD data of mouse and human kidneys
to recapitulate expected nephron segment order and resolve continuous gradients within proximal tubule segments
that clustering obscures. Comparing species, we identify pathways and paired solute transporters whose proximal
tubule positions differ between mouse and human. Finally, we compare control and acutely injured mouse kidneys,
identifying injury-associated shifts in transporter switch points and pathway responses within segments.
Overall, we anticipate pseudostructure reconstruction will enable fine-resolution analysis of gene expression
gradients in structures and tissues of interest.

## Introduction

Spatial transcriptomics (ST) technologies enable gene expression characterization within intact tissue
sections[1–3], and recent platforms such as VisiumHD measure whole sections at 2-µm resolution[4]. Many organs are
built from repeated three-dimensional units whose function changes continuously along their length, such as the
liver lobule, the intestinal villus and the kidney nephron[5,6]. Along the nephron, the proximal tubule alone runs
through three subsegments, S1, S2 and S3, which differ in solute transport, metabolism and susceptibility to
injury[7–9]. Single-cell and spatial studies have mapped these segment states in mouse and human kidneys and in
acute injury[10–15], but how gene expression changes along the length of the tubule, within and between its
segments, remains largely uncharacterized.

However, a two-dimensional section cuts each nephron many times. A single kidney section therefore holds thousands
of tubule cross-sections, each from an unknown position along its nephron, and neighbouring cross-sections often
belong to different nephrons, so their location in the tissue plane does not reveal their position along the
structure. Recovering that position physically requires tracing individual nephrons through serial sections[16,17]
or profiling tissue in three dimensions[18–21], approaches that remain laborious and are rarely paired with
transcriptome-wide measurements across whole organs.

As such, ST analyses of the kidney typically discretize. Clustering and spatial-domain methods assign each spot or
cross-section to a segment label[22–26], and comparisons between conditions or species then contrast segment
means. Segment labels cannot represent a gene that changes within a segment, peaks at a segment boundary, or hands
its function to a partner transporter part-way along the tubule. Methods that model continuous expression
gradients, such as GASTON and Belayer, build their coordinate from tissue x/y positions[27,28], and CoPro takes the
order of nephron segments as an input[29]. Comparative methods align sections to compare expression at matched
locations[30,31], which requires sections that share structure in the tissue plane. Trajectory inference orders
cells along continuous processes from expression alone[32–36], and related ideas have ordered repeated structures:
intestinal crypt positions from deconvolved cell-type composition[37], testis tubules around the known cycle of the
seminiferous epithelium[38,39], and hair follicles by developmental age from three-dimensional profiles[40]. These
approaches rely on a change in cell-type composition, a cyclic topology or three-dimensional data. In the proximal
tubule, the epithelial cell type is shared along the axis, the ends of the structure are not marked, and each
section samples cross-sections at unknown positions.

To address this, we develop pseudostructure reconstruction, which treats every cross-section of a repeated
structure as a sample of one shared structure and recovers its position from its expression, much as pseudotime
recovers progress through a process. Pseudostructure reconstruction segments tubule cross-sections from H&E
histology with a frozen pathology foundation model[41–43] and a trainable decoder, sums the paired VisiumHD bins
within each cross-section, harmonizes sections[44], and fits a single elastic principal curve[45,46] through all
cross-sections, placing each at a position t between 0 and 1 without segment labels or tissue coordinates. Genes
are then modelled along t with negative-binomial generalized additive models[47], each gene on an axis built without
it, and comparisons between groups are calibrated against balanced relabelings of the specimens.

We apply pseudostructure reconstruction to VisiumHD data of four mouse and two human kidney sections. Without
segment labels, it recapitulates the order of nephron segments in every section, remains stable when tubules are
downsampled or the axis is built from half of the genes, and orders segments as well as the best tree-based
trajectory method while fitting one unbranched curve. Within
the proximal tubule, it resolves genes that change inside a segment or peak at a segment boundary, patterns that
segment clustering obscures. Comparing species on a shared axis, it identifies genes whose proximal tubule gradient
reverses between mouse and human, and paired solute transporters whose switch point sits at a different position
along the human proximal tubule, later for glucose and amino-acid transporters. Finally, comparing control and acutely injured mouse kidneys, it identifies transporter switch
points that move upstream as injured S3 loses its late transporters, and gradients within segments that injury
creates or erases. Pseudostructure reconstruction is available as an open-source Python package at [URL].

**References cited in the Introduction** (verified against Europe PMC, Crossref or the publisher; numbering
continues into Results and Methods)

1. Ståhl, P. L. et al. Visualization and analysis of gene expression in tissue sections by spatial transcriptomics. *Science* 353, 78–82 (2016).
2. Marx, V. Method of the Year: spatially resolved transcriptomics. *Nat. Methods* 18, 9–14 (2021).
3. Moses, L. & Pachter, L. Museum of spatial transcriptomics. *Nat. Methods* 19, 534–546 (2022).
4. Oliveira, M. F. et al. High-definition spatial transcriptomic profiling of immune cell populations in colorectal cancer. *Nat. Genet.* 57, 1512–1523 (2025).
5. Halpern, K. B. et al. Single-cell spatial reconstruction reveals global division of labour in the mammalian liver. *Nature* 542, 352–356 (2017).
6. Moor, A. E. et al. Spatial reconstruction of single enterocytes uncovers broad zonation along the intestinal villus axis. *Cell* 175, 1156–1167.e15 (2018).
7. Curthoys, N. P. & Moe, O. W. Proximal tubule function and response to acidosis. *Clin. J. Am. Soc. Nephrol.* 9, 1627–1638 (2014).
8. Lee, J. W., Chou, C.-L. & Knepper, M. A. Deep sequencing in microdissected renal tubules identifies nephron segment-specific transcriptomes. *J. Am. Soc. Nephrol.* 26, 2669–2677 (2015).
9. Chen, L., Chou, C.-L. & Knepper, M. A. A comprehensive map of mRNAs and their isoforms across all 14 renal tubule segments of mouse. *J. Am. Soc. Nephrol.* 32, 897–912 (2021).
10. Ransick, A. et al. Single-cell profiling reveals sex, lineage, and regional diversity in the mouse kidney. *Dev. Cell* 51, 399–413.e7 (2019).
11. Muto, Y. et al. Single cell transcriptional and chromatin accessibility profiling redefine cellular heterogeneity in the adult human kidney. *Nat. Commun.* 12, 2190 (2021).
12. Lake, B. B. et al. An atlas of healthy and injured cell states and niches in the human kidney. *Nature* 619, 585–594 (2023).
13. Kirita, Y., Wu, H., Uchimura, K., Wilson, P. C. & Humphreys, B. D. Cell profiling of mouse acute kidney injury reveals conserved cellular responses to injury. *Proc. Natl Acad. Sci. USA* 117, 15874–15883 (2020).
14. Marshall, J. L. et al. High-resolution Slide-seqV2 spatial transcriptomics enables discovery of disease-specific cell neighborhoods and pathways. *iScience* 25, 104097 (2022).
15. Dixon, E. E., Wu, H., Muto, Y., Wilson, P. C. & Humphreys, B. D. Spatially resolved transcriptomic analysis of acute kidney injury in a female murine model. *J. Am. Soc. Nephrol.* 33, 279–289 (2022).
16. Zhai, X. Y. et al. Digital three-dimensional reconstruction and ultrastructure of the mouse proximal tubule. *J. Am. Soc. Nephrol.* 14, 611–619 (2003).
17. Zhai, X. Y. et al. Three-dimensional reconstruction of the mouse nephron. *J. Am. Soc. Nephrol.* 17, 77–88 (2006).
18. Wang, X. et al. Three-dimensional intact-tissue sequencing of single-cell transcriptional states. *Science* 361, eaat5691 (2018).
19. Vickovic, S. et al. Three-dimensional spatial transcriptomics uncovers cell type localizations in the human rheumatoid arthritis synovium. *Commun. Biol.* 5, 129 (2022).
20. Schott, M. et al. Open-ST: high-resolution spatial transcriptomics in 3D. *Cell* 187, 3953–3972.e26 (2024).
21. Kiemen, A. L. et al. CODA: quantitative 3D reconstruction of large tissues at cellular resolution. *Nat. Methods* 19, 1490–1499 (2022).
22. Traag, V. A., Waltman, L. & van Eck, N. J. From Louvain to Leiden: guaranteeing well-connected communities. *Sci. Rep.* 9, 5233 (2019).
23. Hao, Y. et al. Integrated analysis of multimodal single-cell data. *Cell* 184, 3573–3587.e29 (2021).
24. Zhao, E. et al. Spatial transcriptomics at subspot resolution with BayesSpace. *Nat. Biotechnol.* 39, 1375–1384 (2021).
25. Hu, J. et al. SpaGCN: integrating gene expression, spatial location and histology to identify spatial domains and spatially variable genes by graph convolutional network. *Nat. Methods* 18, 1342–1351 (2021).
26. Long, Y. et al. Spatially informed clustering, integration, and deconvolution of spatial transcriptomics with GraphST. *Nat. Commun.* 14, 1155 (2023).
27. Chitra, U. et al. Mapping the topography of spatial gene expression with interpretable deep learning. *Nat. Methods* 22, 298–309 (2025).
28. Ma, C., Chitra, U., Zhang, S. & Raphael, B. J. Belayer: modeling discrete and continuous spatial variation in gene expression from spatially resolved transcriptomics. *Cell Syst.* 13, 786–797.e13 (2022).
29. Miao, Z. et al. Dissecting the coordinated progression of cell states in spatial transcriptomics with CoPro. Preprint at *bioRxiv* https://doi.org/10.64898/2026.04.17.719309 (2026).
30. Clifton, K. et al. STalign: alignment of spatial transcriptomics data using diffeomorphic metric mapping. *Nat. Commun.* 14, 8123 (2023).
31. Clifton, K. et al. STcompare: comparative spatial transcriptomics data analysis of structurally matched tissues to characterize differentially spatially patterned genes. *Bioinformatics* 42, btag644 (2026).
32. Trapnell, C. et al. The dynamics and regulators of cell fate decisions are revealed by pseudotemporal ordering of single cells. *Nat. Biotechnol.* 32, 381–386 (2014).
33. Haghverdi, L., Büttner, M., Wolf, F. A., Buettner, F. & Theis, F. J. Diffusion pseudotime robustly reconstructs lineage branching. *Nat. Methods* 13, 845–848 (2016).
34. Street, K. et al. Slingshot: cell lineage and pseudotime inference for single-cell transcriptomics. *BMC Genomics* 19, 477 (2018).
35. Wolf, F. A. et al. PAGA: graph abstraction reconciles clustering with trajectory inference through a topology preserving map of single cells. *Genome Biol.* 20, 59 (2019).
36. Saelens, W., Cannoodt, R., Todorov, H. & Saeys, Y. A comparison of single-cell trajectory inference methods. *Nat. Biotechnol.* 37, 547–554 (2019).
37. Niec, R. E. et al. Lymphatics act as a signaling hub to regulate intestinal stem cell activity. *Cell Stem Cell* 29, 1067–1082.e18 (2022).
38. Karin, J., Bornfeld, Y. & Nitzan, M. scPrisma infers, filters and enhances topological signals in single-cell data using spectral template matching. *Nat. Biotechnol.* 41, 1645–1654 (2023).
39. Chakravorty, A. et al. The temporal architecture of the seminiferous epithelial cycle revealed by spatial transcriptomics. *Cell* 189, 4122–4139.e24 (2026).
40. Asami, S., Yin, C., Fan, J., Garza, L. A. & Kalhor, R. Four-dimensional molecular mapping from a spatial snapshot reveals the dynamics of hair follicle organogenesis. *Cell* 189, 5653–5672.e14 (2026).
41. Oquab, M. et al. DINOv2: learning robust visual features without supervision. *Trans. Mach. Learn. Res.* (2024).
42. Karasikov, M. et al. Training state-of-the-art pathology foundation models with orders of magnitude less data. In *Medical Image Computing and Computer Assisted Intervention – MICCAI 2025*, LNCS 15967, 573–583 (Springer, 2025).
43. Kaplan, D. et al. How to train a state-of-the-art pathology foundation model with $1.6k. Technical Report SOPHONT-TR-2025-001 (Sophont, 2025). https://doi.org/10.5281/zenodo.20711012
44. Korsunsky, I. et al. Fast, sensitive and accurate integration of single-cell data with Harmony. *Nat. Methods* 16, 1289–1296 (2019).
45. Gorban, A. N. & Zinovyev, A. Y. Elastic principal graphs and manifolds and their practical applications. *Computing* 75, 359–379 (2005).
46. Albergante, L. et al. Robust and scalable learning of complex intrinsic dataset geometry via ElPiGraph. *Entropy* 22, 296 (2020).
47. Van den Berge, K. et al. Trajectory-based differential expression analysis for single-cell sequencing data. *Nat. Commun.* 11, 1201 (2020).

Methods references continue the numbering: Ronneberger 2015 (U-Net), Kirillov 2019 (panoptic quality), Meyer &
Beucher 1990 (watershed), Yates 2021 (HCOP), Wolf 2018 (Scanpy), Hastie & Stuetzle 1989 (principal curves),
Chen 2021 (microdissection), CZ CELLxGENE Census, MSigDB/Reactome/KEGG and Enrichr. Full entries:
`scratchpad` research notes, to be moved into a reference manager.

## Methods

### Pseudostructure reconstruction overview

Pseudostructure reconstruction takes hematoxylin and eosin (H&E) whole-slide images with registered VisiumHD 2-µm
bins from one or more sections and returns, for every segmented cross-section of a repeated structure (a unit),
its summed expression profile and a position t ∈ [0, 1] along the structure, shared by all sections. Positions
come from expression alone: segment labels and tissue coordinates never enter the placement of a unit. Every step
after segmentation is implemented in the pseudostructure reconstruction package (`pseudospace_reconstruction`
0.1.0, commit 59da04b).

**Tubule segmentation: annotations.** Tubules were segmented from H&E whole-slide images with a deep-learning model
trained on partial manual annotations. Structures were outlined as polygons in QuPath [version] on seven mouse
kidney sections: six from a separate nephrectomy study (three 24 h and three 4 weeks after nephrectomy; not
VisiumHD data) and one control section (Ctrl1A2) that is also analysed here. These contained 2,422 annotated
objects: 1,438 proximal tubules, 457 distal tubules, 336 collecting ducts and 191 glomeruli. For human data, one
cortical section (HUK1_COR1, also analysed here) contained 1,547 objects: 1,443 tubules and 104 glomeruli.
Structures that were not annotated (vessels, interstitium, partially visible tubules) and tissue outside annotated
regions did not contribute to the loss; 52 annotated background regions served as negative examples. Mouse images
were analysed at 0.44 µm per pixel and human images at 0.5 µm per pixel, without resampling or stain
normalization. Separate mouse and human models were trained on patches of 512 × 512 and 640 × 640 pixels. Five
nephrectomy sections were used for training and one was withheld for validation; the control section was divided
into alternating 2,048-pixel blocks for training and validation, and the human section likewise into
non-overlapping training and validation blocks.

**Tubule segmentation: model and training.** Each model combines a frozen pathology foundation model, OpenMidnight
(an open replication of Midnight, a DINOv2 ViT-g/14 with registers and 1.13 billion parameters)[41–43], with a
trainable full-resolution convolutional stem that supplies detail below the 14-pixel patch size, and a U-Net-style
decoder[ref] (7.31 million trainable parameters in total). The decoder has three outputs per pixel: a class map
(background, proximal tubule, distal tubule, collecting duct and glomerulus in mouse; background, tubule and
glomerulus in human); a structure map with four categories (background, object interior, edge facing background,
edge facing another object; a pixel is an edge when a different label occurs within a 3-pixel radius); and a
glomerular centre heatmap. Because annotations were partial, the loss was computed only over annotated objects, a
5-pixel margin around them, annotated background and off-tissue pixels, and combined cross-entropy and soft Dice
terms for the class and structure maps with mean-squared error for the centre map. Models were trained with AdamW
(learning rate 3 × 10⁻⁴, weight decay 0.03, batch size 8 for mouse and 6 for human), with a two-epoch linear
warm-up followed by cosine decay over 60 epochs, geometric, colour and stain augmentation, and oversampling of rare
classes. The checkpoint with the highest validation panoptic quality (PQ)[ref Kirillov 2019] was kept. Each model
was trained once on one NVIDIA RTX 6000 Ada GPU, in about 3.5 h (mouse) and 2.5 h (human).

**Tubule segmentation: instances and accuracy.** Objects were separated by marker-controlled watershed[ref]: seeds
of high interior probability were flooded up to the edge probabilities. Each object took the majority class of
its interior, glomeruli with more than one centre peak were split, and objects below a class-specific minimum area
were discarded. Whole slides were processed in tiles with 25% overlap, averaging predictions over the eight
rotations and reflections of each tile. Accuracy was assessed on annotated objects of the withheld validation
regions, matching predictions one-to-one to annotations at intersection over union (IoU) above 0.5. In mouse, PQ
was 0.880 (segmentation quality 0.914, detection quality 0.963; 98 objects in five regions of the withheld
section); 1.0% of annotated objects were split and 6.1% merged. The human model reached a validation PQ of 0.824.
[See note A: reconcile with the Fig. 2 evaluation.]

**Aggregating bins into units.** Every in-tissue 2-µm bin was assigned to the polygon covering its centre, and its
counts were summed within each polygon, giving one raw-count profile per unit. Morphometrics (area, perimeter,
equivalent diameter, solidity, circularity, axis lengths) were recorded for every polygon.

**Quality control.** Units were kept when they had at least 500 counts and at least 90% of their polygon area
covered by bins. Within each section and segmentation class (at least 50 units per group), units more than three
median absolute deviations from the median of log total counts, log counts per bin or the residual of detected
genes on total counts were removed. Merged neighbouring tubules would mix expression from different segments, so
every unit was also screened for merged objects: conflicting identity (marker panels of two non-adjacent nephron
families score high; flagged, not removed); a split test (the unit's bins are split along their best spatial cut,
the halves are scored for identity, and the separation is tested against 99 label shuffles at α = 0.05; halves of
different, non-adjacent identity were removed, while adjacent ones, likely genuine segment transitions, were kept);
and morphology (solidity below 0.75, or a narrow neck between two lobes). Glomeruli were exempt from the merge
tests. Counts were normalized to 10⁴ per unit and log-transformed for embedding; all gene models use raw counts.

**Sections from two species.** Human genes were mapped to one-to-one mouse orthologs (HCOP[ref]). A pair was kept
when its mouse symbol is measured in every mouse section and its human symbol in every human section; human genes
were renamed to the mouse symbols and every non-shared gene was dropped before quality control, so both species
enter the embedding with the same features. Because the two probe panels differ, expression levels are never
compared across species, only patterns along the axis.

**Embedding, section correction and labels.** The kept units of all sections were embedded jointly in two passes
(Scanpy 1.11.5[ref]; R harmony 2.0.5[44]). In each pass, highly variable genes (mean 0.0125–3, normalized
dispersion ≥ 0.5, intersected across sections) gave 50 principal components, corrected for section with Harmony
(θ = 6, λ = 1, 30 iterations). Units were clustered by Leiden[22] (resolution 2, 30 nearest neighbours, 2,000
highly variable genes), and each cluster was labelled from curated marker panels (PT-S1, PT-S2, PT-S3, thin limb,
thick ascending limb, macula densa, DCT, CNT and collecting-duct segments; best panel score ≥ 0.5 and margin ≥ 0.1
over the runner-up). Non-tubule clusters were dropped and the nephron units re-embedded (pass 2) with genes
detected in at least 5% of units and with at least 20 counts. Labels were reviewed cluster by cluster and none was
overridden. Labels are not inputs to the curve.

**Elastic principal curve.** An elastic principal curve[45,46] was fitted through all pass-2 units in the 50
corrected components. The curve is a chain of k nodes minimizing

U = (1/N) Σᵢ ‖xᵢ − n(xᵢ)‖² + λ Σ_edges ‖nⱼ − nⱼ₊₁‖² + μ Σ_inner nodes ‖nⱼ − (nⱼ₋₁ + nⱼ₊₁)/2‖²,

where n(xᵢ) is the node nearest to unit i; the terms measure fit, stretching and bending. The curve was initialized
as three nodes on the first principal component and grown by bisecting the edge whose bisection gave the lowest
energy, relaxing node positions after each addition, until k = 50 (λ = 0.01, μ = 0.1). Each unit was projected
onto its nearest edge, and its position is the path length from the start of the curve to its projection; units
beyond an end node were placed on the straight continuation of the end edge, up to twice that edge's length. The
start is the curve end nearest the root unit, the medoid of the proximal-tubule cluster lowest on a marker axis
(mean score of early genes Cyp2e1, Gatm, Slc13a3, Slc22a6, Slc5a12, Slc5a2 minus that of late genes Aqp2, Aqp4,
Atp6v0d2, Rhcg, Slc14a2, Wnt7b), which decides only which end is 0; positions were then scaled to [0, 1]. The fit is
deterministic (seed 0). The curve code is a port of elpigraph-python 0.3.4 reduced to the principal-curve case.

### Gene models along the structure

**A position for each gene that the gene did not help build.** A gene that defines the axis could look graded
along it by construction. The genes were therefore split at random into halves five times, and stages 03–04 were
rerun on each half alone, with the other half withdrawn from the highly-variable-gene pool. Each gene's position for
a unit is the mean of the five axes built without that gene; a unit needs a position on at least three of the five.

**PT coordinate x.** The curve runs through every nephron unit, so the proximal tubule (PT) is a stretch of it.
Within the PT (units labelled PT-S1, PT-S2 or PT-S3), x is a unit's rank among the PT units of a reference group of
sections, divided by their number, so x runs in equal numbers of reference tubules. Units of the other group are
mapped through the reference units' empirical distribution (linear interpolation) and clipped to [0, 1].

**Boundaries and cores.** b₁₂ is the threshold on x that misclassifies the fewest reference units between the
PT-S1 and PT-S2 labels, and b₂₃ likewise for PT-S2 and PT-S3, computed on each axis. A segment's core is its x range
more than 0.06 from any boundary.

**Models.** For each gene and section, counts were modelled by negative-binomial generalized linear models with a
log(library size) offset, the library being the unit's total over non-mitochondrial genes: a curve (cubic B-spline
in x with 6 degrees of freedom, knots at reference x quantiles, the same basis in every section), clusters (the
PT-S1/S2/S3 label factor), clusters + curve, and bulk (an intercept)[47]. The dispersion was estimated per gene and
section by maximum likelihood under clusters + curve and held fixed for the nested fits, which were fitted by
iteratively reweighted least squares. Curves are reported as log2 counts per 10⁴ on 201 grid points. Mitochondrial
genes were excluded, and genes were dropped when their rate in any neighbouring segment group (ascending limb, DCT,
collecting duct) exceeded twice their PT rate (neighbour filter).

**Curve against clustering.** The likelihood ratio of clusters + curve against clusters measures what position adds
beyond the segment labels. Its null was a decoy that shuffles x within each label in each section and refits. For
each section, the threshold was the smallest at which (decoy genes above it) / (real genes above it) ≤ 0.05, and a
gene passed when it cleared the threshold in both sections of a group. Passing genes were classed as a peak between
segments (global maximum within 0.06 of a boundary, at least 0.5 log2 above both flanking cores, at the same
boundary in both sections) or a change within a segment (mean of the last third of a core minus its first third
≥ 1 log2, with the same sign and segment in both sections).

### Transporter handoffs

Within each solute-carrier and ABC transporter family (HGNC groups 752 and 417)[ref HGNC], a gene was early in a
species when its end contrast (mean of the last 30% of its curve minus the first 30%) was ≤ −0.5 log2 in both
sections of that species, and late when ≥ +0.5 in both; every early × late pair of a family is a candidate handoff.
Mouse pairs were confirmed when microdissected mouse PT segments[9] showed the early member higher in S1 than S3 and
the late member the reverse, and human pairs likewise against single-cell PT segment profiles of six donors in CZ
CELLxGENE Census[ref]. The confirmed pairs (50 mouse, 16 human; 57 distinct) were fixed before any species or injury
comparison. For each pair and section, both curves were min–max scaled over the grid points inside that section's
2.5–97.5% range of x, and the switch x* is the first point where the early profile minus the late one turns from
positive to non-positive (linear interpolation). Scaling makes x* invariant to a per-gene offset and scale of the
log curve, so neither probe efficiency nor overall gradient amplitude moves it. Its uncertainty came from 200
bootstraps that resample units within each section and refit both curves at fixed dispersion.

### Comparing groups of sections

**Design and calibration.** Each comparison has two groups of two sections. For a per-section quantity Q, the group
contrast is mean(Q, group 2) − mean(Q, group 1), and its four pairs are the cross differences of single sections.
Two balanced relabelings, which put one section of each group into each arm, give the same statistic with the group
effect cancelled; they measure section-to-section variation. A contrast is specific when both relabelings are at
most 0.25 of it in absolute value. Every call requires the four pairs to share its sign, and every list is
calibrated by applying the full rule to the two relabelings (or to the decoy, where stated): the estimated
false-call share is (null calls) / (real calls), and a list above 0.10 is reported as not calibrated.

**Changes within a segment.** For a segment core, W is the curve's mean over the last third minus its mean over the
first third. A gene's response changes inside a segment when |W| of the group contrast is ≥ 1 log2, all four pairs
have that sign with |W| ≥ 0.5, the contrast is specific, the same change is seen in raw counts of units carrying that
segment's label (summed counts, pseudo-count 0.5; ≥ 0.5 log2, all four pairs agreeing), and each third of the core
holds at least 25 such units in every section. This list was calibrated with the decoy.

**Gradients that differ in direction or steepness.** For each section, E is the curve's mean over the S3 core minus
its mean over the S1 core (over grid points within the section's sampled range), and R the same from raw counts of
PT-S3-labelled units in the S3 core and PT-S1-labelled units in the S1 core. A gradient reverses when both sections
of one group share the sign of E, both of the other group have the opposite sign, |E| ≥ 0.5 in all four, and R
agrees in sign with E in each. A gradient is steeper in group 2 when E has one sign in all four sections, |E| ≥ 1 in
both group-2 sections, |E(group 2)| − |E(group 1)| ≥ 0.5 for all four pairs, and R shows the same.

**Switch shifts.** A pair's switch shifted when it has a switch in all four sections, all four pair differences in
x* share a sign, |Δx*| ≥ 0.05, Δx* is specific, and its 95% bootstrap interval excludes 0.

**Pathways.** Gene sets came from MSigDB Hallmark 2020, Reactome 2022 and KEGG 2019 Mouse (Enrichr libraries)[refs],
restricted to screened genes with at least 10 members. A set's statistic is the mean of its members' per-gene
statistic. Each set was compared with 1,000 random sets of the same size drawn from the screened genes (two-sided
empirical p) and Benjamini–Hochberg corrected within each contrast; a call needs q ≤ 0.05, four agreeing pairs and
specificity.

**Gradient retention.** Genes with a change of at least 0.75 log2 in one section were selected, and each other
section's median signed change over the median selected change was reported. Selecting in one section and reading
another avoids regression to the mean.

### Analysis of control and injured mouse kidneys

**Tissue and data.** Kidneys from two control mice and two mice after ischemia–reperfusion injury (IRI) were
profiled, one section per mouse (Ctrl1A2, Ctrl1A4, IR2A2, IR2A4). [Strain, age, IRI protocol (clamp duration, uni-
or bilateral, time to harvest), fixation and embedding, animal protocol number.] All four mice are male (Ddx3y,
Eif2s3y, Kdm5d and Uty are expressed). Sections were processed with the VisiumHD [mouse probe set, version; 19,059
genes], H&E-stained, imaged at [scanner, objective] and processed with Space Ranger [version] to 2-µm bins.

**Reconstruction.** The four sections were reconstructed together. Of 64,272 segmented polygons, 43,984 contained
bins, 36,445 passed quality control and 32,939 nephron units were placed on the curve; PT units numbered 2,780,
2,745, 3,937 and 5,551 (Ctrl1A2, Ctrl1A4, IR2A2, IR2A4).

**Gene models.** The control mice are the reference group for x; fewer than 0.1% of IR units were clipped. The
comparison of curve and clustering used the control sections only, screening genes detected in at least 10% of PT
units in each control mouse. The control–AKI comparison screened genes detected in at least 10% of PT units in each
mouse of at least one condition, so that genes switched on or off by injury are included, with the neighbour filter
applied within each condition. Pathway statistics were W per segment core.

**Placement of injured tubules.** To test whether IR units are placed by injury genes rather than by position,
injury-responsive genes were defined from whole-PT pseudobulk (all four IR × control pairs of the same sign and at
least 0.5 log2; detected in at least 5% of PT units in some section), and the 3,812 such genes were withdrawn from
the highly-variable-gene pool in a rerun of stages 03–04. Each section's PT positions on this injury-blind axis were
compared with the reference by Spearman ρ, over all PT units and within each label, against the ρ between the two
halves of the gene-half splits.

### Analysis of human and mouse kidneys

**Tissue and data.** Two sections of healthy kidney cortex from one male donor (HUK1_COR1 and HUK1_MED1; both are
cortex, the "MED" name does not denote medulla; DDX3Y, KDM5D, UTY and EIF1AY are expressed) were profiled with the
VisiumHD [human probe set, version; 18,132 genes]. [Source of the tissue, donor age, IRB protocol and consent.] They
were analysed together with the two control mouse sections, whose bins were reprocessed in the ortholog space.

**Reconstruction.** The ortholog space held 15,567 genes. Of 51,155 segmented polygons, 31,533 contained bins,
27,395 passed quality control and 25,048 nephron units (13,539 mouse, 11,509 human) were placed on the curve; PT
units numbered 2,644, 2,650, 3,194 and 3,459 (Ctrl1A2, Ctrl1A4, HUK1_COR1, HUK1_MED1).

**Gene models.** The mice are the reference group for x; up to 3.2% of human units, all at the start of the PT, were
clipped. Boundaries come from the mouse units only, because most human units in the S2 range of x are labelled
PT-S1. Genes were screened when detected in at least 10% of PT units in all four sections (a comparison of
patterns needs the gene in both species), with the neighbour filter applied within each species. The steeper-in-mouse
category was counted but not called, because the human sections are flatter overall and the two cannot be told
apart. For pathways, the per-gene statistic was the distance from the genome-wide line of E(human) on E(mouse),
fitted by least squares for each human × mouse pair and averaged over the four pairs, which removes the overall
flattening of human gradients. Switch positions were also re-expressed as quantiles among the switch points of all
graded genes in the same section; a shift that persists on that scale is not explained by a smooth warp between the
two species' axes.

### Evaluation of the axis

**Segment order.** Each unit's label was mapped to seven groups in nephron order (PT-S1 < PT-S2 < PT-S3 <
ascending limb < DCT < CNT < collecting duct), and the Spearman ρ between t and the group rank was computed per
section, against a null that shuffled t within the section 1,000 times. The mouse sections were scored in the
control/AKI reconstruction and the human sections in the human/mouse reconstruction.

**Robustness.** Stages 03–04 were rerun on perturbed inputs and each rerun scored by Spearman ρ with the reference t
over shared units: tubule downsampling to 10–100% of units per section (three seeds per fraction); binomial thinning
of every count to 10–100% of reads, followed by renormalization; and the five gene-half splits, whose two halves were
also scored against each other. Reruns used four BLAS threads, because the Leiden partition depends on the thread
count.

**Cortical depth.** On the mouse sections, each PT unit's distance to the kidney capsule was measured on the H&E
image after masking the hilum and cracks, and correlated with t. Depth is independent of expression and served as a
consistency check only; it never selected a parameter. The human sections have no capsule (every edge is a cut face)
and were not scored.

### Analysis with other methods

The same corrected embedding and root were ordered by diffusion pseudotime[33] (Scanpy 1.11.5; max(30, √n)
neighbours, 15 diffusion components), Slingshot 2.14.0[34] on the pass-2 clusters (R 4.4.3), scFates 1.2.5 (principal
tree)[ref Faure 2023] and GASTON[27] (GLM-PCA input, 30 restarts per section), and each ordering was scored by the
same ρ. [Exact settings: `results/paper_figures/fig3/baselines/methods_notes.md`.]

### Statistics and reproducibility

Sections are the replication units: two per group, with the two human sections from one donor, so human results
describe that donor. No statistical method was used to predetermine sample size, and no data were excluded beyond
the quality-control rules above. Every screen's rules, thresholds and nulls were written down before its first fit;
checks added afterwards are reported as post hoc. Random seeds are fixed in the code. [Runtime and hardware for
stages 01–04 and the gene models.]

### Data availability

[VisiumHD data and H&E images: repository and accession. Segmentation annotations and model weights: where.]

### Code availability

[Pseudostructure reconstruction package and analysis code: GitHub URLs, Zenodo DOI and jef.works documentation.]

---

## Notes for the authors

- **A. Segmentation numbers disagree across our documents.** The paragraph above keeps the numbers from the
  current manuscript draft (mouse PQ 0.880 on 98 objects). The handoff deck gives mouse PQ 0.727 on 1,298 objects,
  and the Fig. 2 ablation gives 0.762 pooled over 454 objects with test-time augmentation. One evaluation should be
  chosen and used in Methods, Fig. 2 and the text.
- **B. Unit counts.** The figures are computed on our local run (64,272 mouse polygons → 32,939 units on the curve;
  human vs mouse 25,048). The posters used the v4fix segmentations (64,080 → 32,785). Methods uses the local run.
- **C. Outdated in the current draft:** "Segment Anything Model 2 (SAM2) encoder" (the encoder is OpenMidnight with
  a convolutional stem); "applying diffusion pseudotime" (the axis is the elastic principal curve; diffusion
  pseudotime is a baseline); "pseudospace" in headings (now "pseudostructure"); "progressive, segment-dependent
  decoupling of transport functions and cell-cycle activation" (the AKI results are switch shifts and
  within-segment gradients; the transport-loss vs cell-cycle program claim failed its relabeling control).
- **D2. Abstract, species sentence.** The calibrated Fig. 5 pathway screen calls no set, so "we identify pathways
  and paired solute transporters whose proximal tubule positions differ" should become "we identify genes and paired
  solute transporters whose proximal tubule positions differ" (the Introduction already says this).
- **E. Prior art to cite and distinguish** (full-text reads, 2026-11-03): SpaceFold (Niec 2022; axis from
  deconvolved cell-type composition), Chakravorty 2026 (hand-segmented testis tubules ordered around a known cycle),
  Asami 2026 (Jean Fan co-author; hair follicles ordered by developmental age from 3D profiles), CoPro (needs x/y,
  takes segment order as input), GASTON/Belayer (x/y). "First to order cross-sections" must not be claimed.
- **D. Title grammar:** "Pseudostructure reconstruction of kidney nephrons enables …" (singular verb).
