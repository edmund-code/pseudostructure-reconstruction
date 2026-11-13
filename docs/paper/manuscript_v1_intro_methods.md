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
data, without segment labels. We apply pseudostructure reconstruction to Visium HD data of mouse and human kidneys
to recapitulate expected nephron segment order and resolve continuous gradients within proximal tubule segments
that clustering obscures. Comparing species, we identify pathways and paired solute transporters whose proximal
tubule positions differ between mouse and human. Finally, we compare control and acutely injured mouse kidneys,
identifying injury-associated shifts in transporter switch points and pathway responses within segments.
Overall, we anticipate pseudostructure reconstruction will enable fine-resolution analysis of gene expression
gradients in structures and tissues of interest.

## Introduction

Spatial transcriptomics (ST) technologies enable gene expression characterization within tissue sections while
retaining each measurement's position[1–3], and Visium HD captures the protein-coding
transcriptome on a continuous grid of 2-µm squares across a 6.5 × 6.5-mm capture area[4]. Many organs
are built from repeated units whose function changes continuously along an axis, such as the liver lobule[5],
the intestinal villus[6] and the kidney nephron[7,8]. Along the nephron, the proximal tubule
runs through three subsegments, S1, S2 and S3, which differ in solute transport and metabolism[9,10]
and in susceptibility to ischaemic injury[8]; in the mouse, they are distinguished by transporter expression
rather than by ultrastructure[11]. Single-cell and spatial studies have mapped these segment states in mouse
and human kidneys and in acute injury[8,12–16], but gene
expression along the tubule is known mostly at the resolution of these three segments.

However, a two-dimensional section cuts each nephron many times. A single kidney section therefore holds thousands
of tubule cross-sections, each from an unknown position along its nephron, and neighbouring cross-sections often
belong to different nephrons, so their location in the tissue plane does not reveal their position along the
structure. Recovering that position physically requires tracing individual nephrons through hundreds to thousands
of serial sections[11,17,18] or profiling tissue in three dimensions[19–21];
these approaches need many sections or assays per specimen and have so far captured morphology, targeted gene
panels or transcriptomes over sub-millimetre depths rather than whole organs.

As such, ST analyses of the kidney typically summarize expression by segment or cell type[15,16].
Clustering assigns each cell, spot or cross-section to one of a few discrete groups that are then annotated as
segments[22,23], spatial-domain methods additionally encourage neighbouring spots in the tissue plane to
share a group[24–26], and comparisons between conditions or species then contrast group means.
Discrete groups cannot represent a gene that changes within a segment, peaks at a segment boundary, or hands its
function to a partner transporter part-way along the tubule, and a spatial prior misleads where neighbouring
cross-sections belong to different nephrons. Methods that model continuous expression gradients build their
coordinate from tissue x/y positions[27–29]; applied to the kidney, CoPro was also given the
order of nephron segments as an input[29]. Comparative methods align sections to compare expression at matched
locations[30,31], which requires sections that share structure in the tissue plane. Trajectory
inference orders cells along continuous processes from their expression profiles, usually given a start cell or
cluster that sets the direction[32–36], and related ideas have
ordered repeated structures: intestinal crypt positions from deconvolved cell-type composition[37], testis
tubules around the known cycle of the seminiferous epithelium[38] using a cyclic topological
template[39], and hair follicles by developmental age, each captured whole in three dimensions[40].
These repeated-structure approaches rely on a change in cell-type composition, a known cyclic topology or
three-dimensional capture of whole structures. In the proximal tubule, the epithelial cell type is shared along the
axis, no landmark marks where a cross-section lies, and each section samples cross-sections at unknown positions.

To address this, we develop pseudostructure reconstruction, which treats every cross-section of a repeated
structure as a sample of one shared structure and recovers its position from its expression, much as pseudotime
recovers progress through a process. Pseudostructure reconstruction segments tubule cross-sections from H&E
histology with a frozen pathology foundation model[41–43] and a trainable decoder, sums
the paired Visium HD bins within each cross-section, harmonizes sections[44], and fits a single elastic
principal curve[45–47] through all cross-sections, placing each at a position t
between 0 and 1 without segment labels or tissue coordinates. Genes are then modelled along t by negative-binomial
regression on a spline basis[48,49], each gene on an axis built without it, and comparisons
between groups are calibrated against balanced relabelings of the specimens.

We apply pseudostructure reconstruction to Visium HD data of four mouse and two human kidney sections. Without
segment labels, it recapitulates the order of nephron segments in every section, remains stable when tubules are
downsampled or the axis is built from half of the genes, and orders segments as well as the best tree-based
trajectory method while fitting one unbranched curve. Within
the proximal tubule, it resolves genes that change inside a segment or peak at a segment boundary, patterns that
segment clustering obscures. Comparing species on a shared axis, it identifies genes whose proximal tubule gradient
reverses between mouse and human, and paired solute transporters whose switch point sits at a different position
along the human proximal tubule, later for glucose and amino-acid transporters. Finally, comparing control and acutely injured mouse kidneys, it identifies transporter switch
points that move upstream as injured S3 loses its late transporters, and gradients within segments that injury
creates or erases. COAST is available as an open-source Python package at github.com/edmundtsou/COAST.

## Methods

### Pseudostructure reconstruction overview

Pseudostructure reconstruction takes hematoxylin and eosin (H&E) whole-slide images with registered Visium HD 2-µm
bins from one or more sections and returns, for every segmented cross-section of a repeated structure (a unit),
its summed expression profile and a position t ∈ [0, 1] along the structure, shared by all sections. Positions
come from expression alone: segment labels and tissue coordinates never enter the placement of a unit, and a
marker axis only orients the curve. t orders units along the structure; it is not physical length, because a
section samples each part of the structure in proportion to the length it lays in the plane. Every step
after segmentation is implemented in COAST (Continuous Ordering of Anatomical Structures from Transcriptomics;
github.com/edmundtsou/COAST, commit db5502c); the analyses of this paper call it through
`analysis/paper/coast_fit.py`.

**Tubule segmentation: annotations.** Tubules were segmented from H&E whole-slide images with a deep-learning model
trained on partial manual annotations. Structures were outlined as polygons in QuPath [version] on seven mouse
kidney sections: six from a separate nephrectomy study (three 24 h and three 4 weeks after nephrectomy; not
Visium HD data) and one control section (Ctrl1A2) that is also analysed here. These contained 2,422 annotated
objects: 1,438 proximal tubules, 457 distal tubules, 336 collecting ducts and 191 glomeruli. For human data, one
cortical section (HUK1_COR1, also analysed here) contained 1,547 objects: 1,443 tubules and 104 glomeruli.
Structures that were not annotated (vessels, interstitium, partially visible tubules) and tissue outside annotated
regions did not contribute to the loss; 52 annotated background regions served as negative examples. Mouse images
were analysed at 0.44 µm per pixel and human images at 0.5 µm per pixel, without resampling or stain
normalization. Separate mouse and human models were trained on patches of 512 × 512 and 640 × 640 pixels. Five
nephrectomy sections were used for training and one was withheld for validation; the control section was divided
into alternating 2,048-pixel blocks for training and validation, and the human section likewise into
non-overlapping training and validation blocks.

**Tubule segmentation: model and training.** Each model combines a frozen pathology foundation model, OpenMidnight[43], an open replication of
Midnight[42] trained on about 12,000 TCGA H&E slides: a DINOv2[41] ViT-g/14 with register tokens[50] (1.13
billion parameters in the released checkpoint). It is joined to a
trainable full-resolution convolutional stem that supplies detail below the 14-pixel patch size, and a U-Net-style
decoder[51] (7.31 million trainable parameters in total). The decoder has three outputs per pixel: a class map
(background, proximal tubule, distal tubule, collecting duct and glomerulus in mouse; background, tubule and
glomerulus in human); a structure map with four categories (background, object interior, edge facing background,
edge facing another object; a pixel is an edge when a different label occurs within a 3-pixel radius); and a
glomerular centre heatmap. Because annotations were partial, the loss was computed only over annotated objects, a
5-pixel margin around them, annotated background and off-tissue pixels, and combined cross-entropy and soft Dice
terms for the class and structure maps with mean-squared error for the centre map. Models were trained with AdamW
(learning rate 3 × 10⁻⁴, weight decay 0.03, batch size 8 for mouse and 6 for human), with a two-epoch linear
warm-up followed by cosine decay over 60 epochs, geometric, colour and stain augmentation, and oversampling of rare
classes. The checkpoint with the highest validation class-agnostic panoptic quality (PQ)[52] was kept. Each model
was trained once on one NVIDIA RTX 6000 Ada GPU, in about 3.5 h (mouse) and 2.5 h (human).

**Tubule segmentation: instances and accuracy.** Objects were separated by marker-controlled watershed[53]: seeds
of high interior probability were flooded up to the edge probabilities. Each object took the majority class of
its interior, glomeruli with more than one centre peak were split, and objects below a class-specific minimum area
were discarded. Whole slides were processed in tiles with 25% overlap, averaging predictions over the eight
rotations and reflections of each tile. Accuracy was assessed on annotated objects of the withheld validation
regions, matching predictions one-to-one to annotations at intersection over union (IoU) above 0.5. In mouse, PQ
was 0.880 (segmentation quality 0.914, recognition quality 0.963; 98 objects in five regions of the withheld
section); 1.0% of annotated objects were split and 6.1% merged. The human model reached a validation PQ of 0.824.
[See note A: reconcile with the Fig. 2 evaluation.]

**Aggregating bins into units.** Every in-tissue 2-µm bin was assigned to the polygon covering its centre, and its
counts were summed within each polygon, giving one raw-count profile per unit, as Visium HD 2-µm bins have been
summed within segmented nuclei[4]. Morphometrics (area, perimeter,
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
tests (`coast.pp.filter_units`, `coast.pp.doublets`). Counts were normalized to 10⁴ per unit and log-transformed for
embedding; all gene models use raw counts.

**Sections from two species.** Human genes were mapped to one-to-one mouse orthologs from HCOP[54]: pairs
supported by at least three of its orthology resources, kept when each gene is the other's best-supported
partner. A pair was used when its mouse symbol is measured in every mouse section and its human symbol in every
human section; human genes
were renamed to the mouse symbols and every non-shared gene was dropped before quality control, so both species
enter the embedding with the same features. Because the two probe panels differ, expression levels are never
compared across species, only patterns along the axis.

**Embedding, section correction and labels.** The kept units of all sections were embedded jointly in two passes
(Scanpy 1.11.5[55]; harmonypy 2.1.0, the Python port of Harmony[44]). In each pass, genes detected in at least 5% of
units and with at least 20 counts (plus the marker-panel genes) were kept, and highly variable genes (Scanpy's
Seurat flavour: mean 0.0125–3, normalized dispersion ≥ 0.5; ranked within each section, and kept only when variable
in every condition or species; mitochondrial and ribosomal genes excluded) gave 50 principal components, corrected
for section with Harmony (θ = 6, λ fixed at 1, at most 30 rounds). Units were clustered by Leiden[22] (resolution 2)
on the 25 exact nearest neighbours in the corrected space, and each cluster was labelled from curated marker panels
(PT-S1, PT-S2, PT-S3, thin limb, thick ascending limb, macula densa, DCT, CNT and collecting-duct segments; best mean
panel score ≥ 0.5, control genes drawn from 2,000 highly variable genes and the panel genes; the margin over the
runner-up is reported, not thresholded). Non-tubule clusters were dropped and the nephron units re-embedded (pass 2)
by the same steps. Labels are not inputs to the curve. [Labels were assigned automatically; cluster-by-cluster review
pending.]

**Elastic principal curve.** An elastic principal curve[46,47] was fitted through all pass-2 units in the 50
corrected components. The curve is a chain of k nodes minimizing

U = (1/N) Σᵢ ‖xᵢ − n(xᵢ)‖² + λ Σ_edges ‖nⱼ − nⱼ₊₁‖² + μ Σ_inner nodes ‖nⱼ − (nⱼ₋₁ + nⱼ₊₁)/2‖²,

where n(xᵢ) is the node nearest to unit i; the terms measure fit, stretching and bending. COAST fits it with
elpigraph-python 0.3.4. The curve was initialized
as three nodes on the first principal component (the elpigraph-python default) and grown by bisecting the edge whose
bisection gave the lowest energy[47], relaxing node positions after each addition, until k = 50 (λ = 0.01, μ = 0.1). Each unit was projected
onto its nearest edge, and its position is the path length from the start of the curve to its projection (the projection index of
principal curves[45]); units
beyond an end node were placed on the straight continuation of the end edge, up to twice that edge's length
(our addition to the curve projection). The
direction was chosen so that position rises with a marker axis (mean score of late genes Aqp2, Aqp4, Atp6v0d2, Rhcg,
Slc14a2, Wnt7b minus that of early genes Cyp2e1, Gatm, Slc13a3, Slc22a6, Slc5a12, Slc5a2; Spearman correlation),
which decides only which end is 0; positions were then scaled to [0, 1]. The fit is deterministic (seed 0).

### Gene models along the structure

**A position for each gene that the gene did not help build.** A gene that defines the axis could look graded
along it by construction. The genes were therefore split at random into halves five times, and both passes and the curve were
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
PT-S1/S2/S3 label factor), clusters + curve, and bulk (an intercept). The spline follows tradeSeq[49] (six cubic basis functions, knots at quantiles, a
library-size offset) but carries no smoothing penalty. The dispersion was estimated per gene and section by maximum
likelihood under clusters + curve and held fixed for the nested fits, which were fitted by iteratively reweighted
least squares and compared by likelihood ratio[48]. Curves are reported as log2 counts per 10⁴ on 201 grid points. Mitochondrial
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

Within each transporter family (the subgroups of HGNC gene groups 752, solute carriers, and 417, ATP binding
cassette transporters[56]; downloaded 2026-10-12), a gene was early in a
species when its end contrast (mean of the last 30% of its curve minus the first 30%) was ≤ −0.5 log2 in both
sections of that species, and late when ≥ +0.5 in both; every early × late pair of a family is a candidate handoff.
Mouse pairs were confirmed when microdissected mouse PT segments[10] showed the early member higher in S1 than S3 and
the late member the reverse, and human pairs likewise against single-cell profiles of proximal convoluted tubule (S1 and S2) and S3 cells
from the six donors with at least 50 cells of each in CZ CELLxGENE Census (release 2025-11-08)[57]. The confirmed pairs (50 mouse, 16 human; 57 distinct) were fixed before any species or injury
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

**Pathways.** Gene sets came from MSigDB Hallmark 2020[58], Reactome 2022[59] and KEGG 2019
Mouse[60,61], as distributed by Enrichr[62,63]; human gene symbols were mapped
to mouse through the one-to-one ortholog map, and sets were restricted to screened genes with at least 10 members. A set's statistic is the mean of its members' per-gene
statistic. Each set was compared with 1,000 random sets of the same size drawn from the screened genes (two-sided
empirical p) and Benjamini–Hochberg corrected[64] within each contrast; a call needs q ≤ 0.05, four agreeing pairs and
specificity.

**Gradient retention.** Genes with a change of at least 0.75 log2 in one section were selected, and each other
section's median signed change over the median selected change was reported. Selecting in one section and reading
another avoids regression to the mean.

### Analysis of control and injured mouse kidneys

**Tissue and data.** Kidneys from two control mice and two mice after ischemia–reperfusion injury (IRI) were
profiled, one section per mouse (Ctrl1A2, Ctrl1A4, IR2A2, IR2A4). [Strain, age, IRI protocol (clamp duration, uni-
or bilateral, time to harvest), fixation and embedding, animal protocol number.] All four mice are male (Ddx3y,
Eif2s3y, Kdm5d and Uty are expressed). Sections were processed with the Visium HD [mouse probe set, version; 19,059
genes], H&E-stained, imaged at [scanner, objective] and processed with Space Ranger [version] to 2-µm bins.

**Reconstruction.** The four sections were reconstructed together. Of 64,080 segmented polygons, 43,761 contained
bins, 36,610 passed quality control and 31,928 nephron units were placed on the curve; PT units numbered 2,864,
2,757, 3,385 and 5,104 (Ctrl1A2, Ctrl1A4, IR2A2, IR2A4).

**Gene models.** The control mice are the reference group for x; at most 0.24% of IR units were clipped. The
comparison of curve and clustering used the control sections only, screening genes detected in at least 10% of PT
units in each control mouse. The control–AKI comparison screened genes detected in at least 10% of PT units in each
mouse of at least one condition, so that genes switched on or off by injury are included, with the neighbour filter
applied within each condition. Pathway statistics were W per segment core.

**Placement of injured tubules.** To test whether IR units are placed by injury genes rather than by position,
injury-responsive genes were defined from whole-PT pseudobulk (all four IR × control pairs of the same sign and at
least 0.5 log2; detected in at least 5% of PT units in some section), and the 3,635 such genes were withdrawn from
the highly-variable-gene pool in a rerun of both passes and the curve. Each section's PT positions on this injury-blind axis were
compared with the reference by Spearman ρ, over all PT units and within each label, against the ρ between the two
halves of the gene-half splits.

### Analysis of human and mouse kidneys

**Tissue and data.** Two sections of healthy kidney cortex from one male donor (HUK1_COR1 and HUK1_MED1; both are
cortex, the "MED" name does not denote medulla; DDX3Y, KDM5D, UTY and EIF1AY are expressed) were profiled with the
Visium HD [human probe set, version; 18,132 genes]. [Source of the tissue, donor age, IRB protocol and consent.] They
were analysed together with the two control mouse sections, whose bins were reprocessed in the ortholog space.

**Reconstruction.** The ortholog space held 15,567 genes. Of 51,176 segmented polygons, 31,524 contained bins,
27,614 passed quality control and 25,352 nephron units (13,920 mouse, 11,432 human) were placed on the curve; PT
units numbered 2,712, 2,740, 3,140 and 3,392 (Ctrl1A2, Ctrl1A4, HUK1_COR1, HUK1_MED1).

**Gene models.** The mice are the reference group for x; at most 0.06% of human units were clipped. Boundaries come from the mouse units only, because most human units in the S2 range of x are labelled
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

The same corrected embedding and root were ordered by diffusion pseudotime[33] as implemented in Scanpy 1.11.5[55]
(max(30, √n) neighbours, 15 diffusion components, rooted at the same root unit), Slingshot 2.14.0[34] on the pass-2
Leiden clusters (R 4.4.3; starting from the root unit's cluster, with each unit's per-lineage pseudotimes averaged by
its lineage weights), scFates 1.2.5[65] (SimplePPT principal tree[66]) and GASTON[27] (GLM-PCA input, 30 restarts per section), and each ordering was scored by the
same ρ. [Exact settings: `results/paper_figures/fig3/baselines/methods_notes.md`.]

### Statistics and reproducibility

Sections are the replication units: two per group, with the two human sections from one donor, so human results
describe that donor. All specimens are male; in mouse, S3 is the most sex-divergent part of the PT[8], so
mouse S3 genes describe male mice. No statistical method was used to predetermine sample size, and no data were excluded beyond
the quality-control rules above. Every screen's rules, thresholds and nulls were written down before its first fit;
checks added afterwards are reported as post hoc. Random seeds are fixed in the code. [Runtime and hardware for
stages 01–04 and the gene models.]

### Data availability

[Visium HD data and H&E images: repository and accession. Segmentation annotations and model weights: where.]

### Code availability

[COAST: github.com/edmundtsou/COAST; analysis code: GitHub URL, Zenodo DOI and jef.works documentation.]

## References

Numbered by first citation (Introduction, then Methods). Every reference was checked against its full text on
2026-11-03 unless marked; the claim-by-claim record is `docs/paper/citation_ledger.md`.

1. Ståhl, P. L. et al. Visualization and analysis of gene expression in tissue sections by spatial transcriptomics. *Science* 353, 78–82 (2016).
2. Method of the Year 2020: spatially resolved transcriptomics. *Nat. Methods* 18, 1 (2021). doi:10.1038/s41592-020-01042-x
3. Moses, L. & Pachter, L. Museum of spatial transcriptomics. *Nat. Methods* 19, 534–546 (2022). **[bioRxiv preprint read]**
4. Oliveira, M. F. et al. High-definition spatial transcriptomic profiling of immune cell populations in colorectal cancer. *Nat. Genet.* 57, 1512–1523 (2025).
5. Halpern, K. B. et al. Single-cell spatial reconstruction reveals global division of labour in the mammalian liver. *Nature* 542, 352–356 (2017).
6. Moor, A. E. et al. Spatial reconstruction of single enterocytes uncovers broad zonation along the intestinal villus axis. *Cell* 175, 1156–1167.e15 (2018).
7. Lee, J. W., Chou, C.-L. & Knepper, M. A. Deep sequencing in microdissected renal tubules identifies nephron segment-specific transcriptomes. *J. Am. Soc. Nephrol.* 26, 2669–2677 (2015).
8. Ransick, A. et al. Single-cell profiling reveals sex, lineage, and regional diversity in the mouse kidney. *Dev. Cell* 51, 399–413.e7 (2019).
9. Curthoys, N. P. & Moe, O. W. Proximal tubule function and response to acidosis. *Clin. J. Am. Soc. Nephrol.* 9, 1627–1638 (2014).
10. Chen, L., Chou, C.-L. & Knepper, M. A. A comprehensive map of mRNAs and their isoforms across all 14 renal tubule segments of mouse. *J. Am. Soc. Nephrol.* 32, 897–912 (2021).
11. Zhai, X. Y. et al. Digital three-dimensional reconstruction and ultrastructure of the mouse proximal tubule. *J. Am. Soc. Nephrol.* 14, 611–619 (2003).
12. Muto, Y. et al. Single cell transcriptional and chromatin accessibility profiling redefine cellular heterogeneity in the adult human kidney. *Nat. Commun.* 12, 2190 (2021).
13. Lake, B. B. et al. An atlas of healthy and injured cell states and niches in the human kidney. *Nature* 619, 585–594 (2023).
14. Kirita, Y., Wu, H., Uchimura, K., Wilson, P. C. & Humphreys, B. D. Cell profiling of mouse acute kidney injury reveals conserved cellular responses to injury. *Proc. Natl Acad. Sci. USA* 117, 15874–15883 (2020).
15. Melo Ferreira, R. et al. Integration of spatial and single-cell transcriptomics localizes epithelial cell-immune cross-talk in kidney injury. *JCI Insight* 6, e147703 (2021).
16. Dixon, E. E., Wu, H., Muto, Y., Wilson, P. C. & Humphreys, B. D. Spatially resolved transcriptomic analysis of acute kidney injury in a female murine model. *J. Am. Soc. Nephrol.* 33, 279–289 (2022).
17. Zhai, X. Y. et al. Three-dimensional reconstruction of the mouse nephron. *J. Am. Soc. Nephrol.* 17, 77–88 (2006).
18. Kiemen, A. L. et al. CODA: quantitative 3D reconstruction of large tissues at cellular resolution. *Nat. Methods* 19, 1490–1499 (2022).
19. Wang, X. et al. Three-dimensional intact-tissue sequencing of single-cell transcriptional states. *Science* 361, eaat5691 (2018).
20. Vickovic, S. et al. Three-dimensional spatial transcriptomics uncovers cell type localizations in the human rheumatoid arthritis synovium. *Commun. Biol.* 5, 129 (2022).
21. Schott, M. et al. Open-ST: high-resolution spatial transcriptomics in 3D. *Cell* 187, 3953–3972.e26 (2024).
22. Traag, V. A., Waltman, L. & van Eck, N. J. From Louvain to Leiden: guaranteeing well-connected communities. *Sci. Rep.* 9, 5233 (2019).
23. Hao, Y. et al. Integrated analysis of multimodal single-cell data. *Cell* 184, 3573–3587.e29 (2021).
24. Zhao, E. et al. Spatial transcriptomics at subspot resolution with BayesSpace. *Nat. Biotechnol.* 39, 1375–1384 (2021).
25. Hu, J. et al. SpaGCN: integrating gene expression, spatial location and histology to identify spatial domains and spatially variable genes by graph convolutional network. *Nat. Methods* 18, 1342–1351 (2021). **[Research Square preprint read]**
26. Long, Y. et al. Spatially informed clustering, integration, and deconvolution of spatial transcriptomics with GraphST. *Nat. Commun.* 14, 1155 (2023).
27. Chitra, U. et al. Mapping the topography of spatial gene expression with interpretable deep learning. *Nat. Methods* 22, 298–309 (2025).
28. Ma, C., Chitra, U., Zhang, S. & Raphael, B. J. Belayer: modeling discrete and continuous spatial variation in gene expression from spatially resolved transcriptomics. *Cell Syst.* 13, 786–797.e13 (2022).
29. Miao, Z. et al. Dissecting the coordinated progression of cell states in spatial transcriptomics with CoPro. Preprint at *bioRxiv* https://doi.org/10.64898/2026.04.17.719309 (2026). **[preprint (not peer reviewed)]**
30. Clifton, K. et al. STalign: alignment of spatial transcriptomics data using diffeomorphic metric mapping. *Nat. Commun.* 14, 8123 (2023).
31. Clifton, K. et al. STcompare: comparative spatial transcriptomics data analysis of structurally matched tissues to characterize differentially spatially patterned genes. *Bioinformatics* 42, btag644 (2026).
32. Trapnell, C. et al. The dynamics and regulators of cell fate decisions are revealed by pseudotemporal ordering of single cells. *Nat. Biotechnol.* 32, 381–386 (2014). **[main text read; Online Methods paywalled]**
33. Haghverdi, L., Büttner, M., Wolf, F. A., Buettner, F. & Theis, F. J. Diffusion pseudotime robustly reconstructs lineage branching. *Nat. Methods* 13, 845–848 (2016). **[accepted author manuscript read; copyedited version paywalled]**
34. Street, K. et al. Slingshot: cell lineage and pseudotime inference for single-cell transcriptomics. *BMC Genomics* 19, 477 (2018).
35. Wolf, F. A. et al. PAGA: graph abstraction reconciles clustering with trajectory inference through a topology preserving map of single cells. *Genome Biol.* 20, 59 (2019).
36. Saelens, W., Cannoodt, R., Todorov, H. & Saeys, Y. A comparison of single-cell trajectory inference methods. *Nat. Biotechnol.* 37, 547–554 (2019). **[authors' manuscript source read]**
37. Niec, R. E. et al. Lymphatics act as a signaling hub to regulate intestinal stem cell activity. *Cell Stem Cell* 29, 1067–1082.e18 (2022).
38. Chakravorty, A. et al. The temporal architecture of the seminiferous epithelial cycle revealed by spatial transcriptomics. *Cell* 189, 4122–4139.e24 (2026).
39. Karin, J., Bornfeld, Y. & Nitzan, M. scPrisma infers, filters and enhances topological signals in single-cell data using spectral template matching. *Nat. Biotechnol.* 41, 1645–1654 (2023).
40. Asami, S., Yin, C., Fan, J., Garza, L. A. & Kalhor, R. Four-dimensional molecular mapping from a spatial snapshot reveals the dynamics of hair follicle organogenesis. *Cell* 189, 5653–5672.e14 (2026).
41. Oquab, M. et al. DINOv2: learning robust visual features without supervision. *Trans. Mach. Learn. Res.* (2024). https://openreview.net/forum?id=a68SUt6zFt
42. Karasikov, M. et al. Training state-of-the-art pathology foundation models with orders of magnitude less data. In *Medical Image Computing and Computer Assisted Intervention – MICCAI 2025*, LNCS 15967, 573–583 (Springer, 2025).
43. Kaplan, D., Grandhi, R. S., Lane, C., Warner, B., Abraham, T. M. & Scotti, P. S. How to train a state-of-the-art pathology foundation model with $1.6k. Technical Report SOPHONT-TR-2025-001 (Sophont, 2025). https://doi.org/10.5281/zenodo.20711012
44. Korsunsky, I. et al. Fast, sensitive and accurate integration of single-cell data with Harmony. *Nat. Methods* 16, 1289–1296 (2019).
45. Hastie, T. & Stuetzle, W. Principal curves. *J. Am. Stat. Assoc.* 84, 502–516 (1989).
46. Gorban, A. N. & Zinovyev, A. Y. Elastic principal graphs and manifolds and their practical applications. *Computing* 75, 359–379 (2005).
47. Albergante, L. et al. Robust and scalable learning of complex intrinsic dataset geometry via ElPiGraph. *Entropy* 22, 296 (2020).
48. McCarthy, D. J., Chen, Y. & Smyth, G. K. Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation. *Nucleic Acids Res.* 40, 4288–4297 (2012).
49. Van den Berge, K. et al. Trajectory-based differential expression analysis for single-cell sequencing data. *Nat. Commun.* 11, 1201 (2020).
50. Darcet, T., Oquab, M., Mairal, J. & Bojanowski, P. Vision transformers need registers. In *International Conference on Learning Representations* (2024). arXiv:2309.16588
51. Ronneberger, O., Fischer, P. & Brox, T. U-Net: convolutional networks for biomedical image segmentation. In *Medical Image Computing and Computer-Assisted Intervention – MICCAI 2015*, LNCS 9351, 234–241 (Springer, 2015).
52. Kirillov, A., He, K., Girshick, R., Rother, C. & Dollár, P. Panoptic segmentation. In *Proc. IEEE/CVF Conference on Computer Vision and Pattern Recognition* 9396–9405 (2019).
53. Beucher, S. & Meyer, F. The morphological approach to segmentation: the watershed transformation. In *Mathematical Morphology in Image Processing* (ed. Dougherty, E. R.) 433–481 (Marcel Dekker, 1993).
54. Yates, B., Gray, K. A., Jones, T. E. M. & Bruford, E. A. Updates to HCOP: the HGNC comparison of orthology predictions tool. *Brief. Bioinform.* 22, bbab155 (2021).
55. Wolf, F. A., Angerer, P. & Theis, F. J. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol.* 19, 15 (2018).
56. Seal, R. L. et al. Genenames.org: the HGNC resources in 2023. *Nucleic Acids Res.* 51, D1003–D1009 (2023).
57. CZI Cell Science Program et al. CZ CELLxGENE Discover: a single-cell data platform for scalable exploration, analysis and modeling of aggregated data. *Nucleic Acids Res.* 53, D886–D900 (2025).
58. Liberzon, A. et al. The Molecular Signatures Database (MSigDB) hallmark gene set collection. *Cell Syst.* 1, 417–425 (2015).
59. Gillespie, M. et al. The reactome pathway knowledgebase 2022. *Nucleic Acids Res.* 50, D687–D692 (2022).
60. Kanehisa, M. & Goto, S. KEGG: Kyoto Encyclopedia of Genes and Genomes. *Nucleic Acids Res.* 28, 27–30 (2000).
61. Kanehisa, M., Sato, Y., Furumichi, M., Morishima, K. & Tanabe, M. New approach for understanding genome variations in KEGG. *Nucleic Acids Res.* 47, D590–D595 (2019).
62. Kuleshov, M. V. et al. Enrichr: a comprehensive gene set enrichment analysis web server 2016 update. *Nucleic Acids Res.* 44, W90–W97 (2016).
63. Xie, Z. et al. Gene set knowledge discovery with Enrichr. *Curr. Protoc.* 1, e90 (2021).
64. Benjamini, Y. & Hochberg, Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J. R. Stat. Soc. B* 57, 289–300 (1995).
65. Faure, L., Soldatov, R., Kharchenko, P. V. & Adameyko, I. scFates: a scalable python package for advanced pseudotime and bifurcation analysis from single-cell data. *Bioinformatics* 39, btac746 (2023).
66. Mao, Q., Yang, L., Wang, L., Goodison, S. & Sun, Y. SimplePPT: a simple principal tree algorithm. In *Proc. 2015 SIAM International Conference on Data Mining* 792–800 (2015).

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
- **F. From the full-text citation check (2026-11-03; record in `docs/paper/citation_ledger.md`).**
  - Our first sentence is close to the opening of STalign's abstract (same lab: "Spatial transcriptomics (ST)
    technologies enable high throughput gene expression characterization within thin tissue sections").
  - Asami 2026 (Jean Fan co-author) names kidney nephrons as a target for its "spatial pseudotiming". The Discussion
    should separate the two: it orders whole structures that develop asynchronously by developmental age; we order
    cross-sections of mature structures by position along them, from 2D sections in which no structure is whole.
  - CoPro and Asami already compare samples on a shared expression-derived axis, so a shared axis is not new by
    itself; what is new is an axis built with no tissue coordinates, labels, reference program or whole-structure
    capture.
  - None of the repeated-structure precedents compares its axis with a measured position along the structure
    (Asami correlates with follicle length, ρ 0.94/0.90), so our segment-order ρ and stability tests are the fair
    comparison.
  - All 66 references have now been read in full (for DPT, the accepted author manuscript; for Moses 2022, SpaGCN and
    CoPro, preprints).
- **D. Title grammar:** "Pseudostructure reconstruction of kidney nephrons enables …" (singular verb).
