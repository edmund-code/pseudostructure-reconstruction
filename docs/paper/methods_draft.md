# Methods (draft, workstream 4)

Status: draft written 2026-10-04 from the code of notebooks 01, 03, 12, 13, 31, 33 and 34, the
`pseudospace/` and `segmentation/` packages, and saved result tables. Parameters were read from
code, not from memory. Counts come from the saved runs named in brackets.

Markers used in this draft:

- **[VERIFY]**: not recorded in the repository or its outputs. A person with the lab records must
  confirm or supply it.
- **[WS1]**, **[WS2]**, **[WS3]**: depends on a decision by workstream 1 (reconstruction method),
  2 (final pathway reporting set) or 3 (biology and literature).
- **[CITE]**: needs a reference that has not been verified.

References are numbered and listed at the end. Each one was checked against Europe PMC, Crossref
or the publisher page.

---

## Study design and specimens

The study compares healthy proximal tubule (PT) in two species. We used two control mouse
kidney sections (specimens `Ctrl1A2` and `Ctrl1A4`) and two sections of renal cortex from one
human donor (`HUK1_COR1` and `HUK1_MED1`). The second human section is named "MED" at source,
but its processed segmentation contains cortex only. It is analysed and described as cortex
throughout. All four specimens are male: Y-linked transcripts are detected in every specimen
[VERIFY against donor and animal records].
The mouse sections are whole-kidney sections that include cortex and outer medulla. The human
sections contain cortex only, so mouse S3 tubules in the outer stripe of the outer medulla have no
human counterpart (see *Sampled-region check*). Specimens, not tubules, are the units of
replication. With two specimens per species, every species comparison in this study is
descriptive.

- Mouse: strain, age, supplier, housing, diet and euthanasia method [VERIFY]. The task brief
  describes the mice as C57BL/6-type; no file in the repository records the strain [VERIFY]. The
  same experiment includes two ischaemia–reperfusion kidneys (`IR2A2`, `IR2A4`) that this paper
  does not analyse.
- Human: tissue source (for example tumour-free nephrectomy tissue or a deceased-donor kidney),
  donor age, kidney function, and whether the two sections are serial sections or come from
  separate blocks [VERIFY].
- Tissue processing: fixation and embedding (FFPE or fixed frozen), section thickness, H&E
  staining protocol and imaging system [VERIFY]. A CytAssist image is present for the human
  sections, which indicates a CytAssist-based workflow [VERIFY].
- Ethics: animal-protocol approval and human-tissue consent or IRB approval [VERIFY].

## Visium HD data generation and preprocessing

Sections were profiled with Visium HD spatial gene expression, a probe-based assay at 2-µm bin
resolution [1]. Library preparation, sequencing platform and read depth are not recorded here
[VERIFY].

Raw data were processed with Space Ranger. The two species used different versions and
references. These values were read from the output metadata:

| | Mouse (`Ctrl1A2`, `Ctrl1A4`) | Human (`HUK1_COR1`, `HUK1_MED1`) |
|---|---|---|
| Space Ranger | 3.1.1 | 4.0.1 |
| Probe set | Visium Mouse Transcriptome Probe Set v2.0 | Visium Human Transcriptome Probe Set v2.1.0 |
| Genome reference | mm10 | GRCh38 (GRCh38-2024-A) |
| Chemistry | Visium HD v1 | Visium HD H1 slide, HD probe-based v1 |
| Features in the filtered matrix | 19,059 | 18,132 |

We used the 2-µm square-bin output (`square_002um`): the filtered feature-barcode matrix and the
tissue-position table. The workflow confirms the bin size from the barcode prefix (`s_002um_`),
so a coarser 8-µm or 16-µm bin cannot be used by mistake. Only bins flagged as under tissue were
kept.

## Tubule segmentation

**Model.** Tubules and glomeruli were segmented on H&E whole-slide images with `kidney_panoptic`,
a panoptic instance segmentation model developed for this project (repository `segmentation/`,
checkpoint "v4").

- *Encoder.* A frozen pathology foundation model: OpenMidnight, a DINOv2 ViT-g/14 with registers
  that was fine-tuned on H&E images [2,3]. Tokens were taken from blocks 9, 19, 29 and 39 at an
  input size of 518 px.
- *Decoder.* A trainable convolutional stem on native pixels (strides 1–16) feeds a decoder that
  fuses the tokens at stride 16 and outputs predictions at stride 1. It has 7.31 million trainable
  parameters, uses GroupNorm and dropout 0.10.
- *Heads.* Three heads: a 5-class semantic head (background, proximal, distal, collecting-duct
  tubule, glomerulus); a 4-class boundary head (background, interior, boundary to background,
  boundary to another instance); and a centre heat map for glomeruli.
- *Working resolution.* 0.44068 µm per pixel.

**Instance decoding.** Instances were produced by marker-controlled watershed:

1. Foreground is 1 − p(background) ≥ 0.50.
2. Seeds are connected components of p(interior) ≥ 0.40 with area ≥ 24 px².
3. The watershed floods outward from the seeds over the boundary probability.
4. Each instance takes the majority semantic class of its pixels, with the boundary band excluded.
5. Instances below a per-class minimum area are removed (px²): proximal 1,950; distal 1,000;
   collecting duct 950; glomerulus 2,300.

Glomeruli containing two or more centre peaks were split further. Fragment merging and
semantic splitting were disabled.

**Whole-slide inference.**

- Tiles of 512 px with 128 px overlap were blended with a Gaussian window (σ = 0.25 × tile).
- Decode windows were 2,048 px with 256 px overlap. Instances crossing window borders were
  stitched at IoU ≥ 0.5.
- A saturation-based tissue mask removed instances that lay more than half off tissue.
- Reported inference used D4 test-time augmentation: the probability maps from the 8 rotations
  and reflections were averaged.
- Polygons were simplified with a 1-px tolerance and capped at 400 vertices.

**Training and model selection.**

- *Data.* Partially annotated mouse H&E slides: 6,475 patches of 512 px (stride 256) from 7 slide
  regions of nephrectomy and control kidneys. Part of `Ctrl1A2` was in the training data.
- *Ignore mask.* Annotated instances were trusted foreground. A 5-px ring around each instance
  and all off-tissue glass were trusted background. Unannotated tissue was ignored in the loss.
- *Loss.* Cross-entropy plus Dice on the semantic and boundary heads. Boundary pixels were
  up-weighted (×3) and instance-to-instance boundary pixels more (×6). Mean squared error was
  used for the centre head.
- *Optimisation.* AdamW, learning rate 3 × 10⁻⁴ with cosine decay to 10⁻⁶ and 2 warm-up epochs,
  weight decay 0.03, batch size 8, mixed precision, gradient clipping at 1.0. Augmentation
  covered geometry, stain (HED), blur, noise, JPEG artefacts, cutout, elastic deformation and
  copy-paste of rare classes.
- *Checkpoint selection.* The checkpoint was chosen by pooled, class-agnostic panoptic quality
  (PQ) [4] on the validation set, never by validation loss. Ties within 0.005 PQ went to the lower
  merge rate. The selected checkpoint is epoch 35 of 54.
- *Performance.* On five 4,096-px regions of a fully held-out nephrectomy slide (98 objects), PQ
  was 0.837 without test-time augmentation and 0.880 with D4 test-time augmentation [VERIFY
  whether to report these numbers, which come from `segmentation/model_summary_v2.txt`].
- *Known limitations.* Predicted objects are smaller than the annotations. The most common
  error is merging of touching tubules of the same type.

**Human segmentation.** The human sections use separate re-segmentations
(`HUK1_COR1_v2.geojson`, `HUK1_MED1_v2.geojson`). The model was trained only on mouse slides. How
the human segmentations were produced (model, checkpoint, test-time augmentation, any manual
editing), and whether they were validated, is not recorded [VERIFY].

**Quality control of mouse segmentations.** Mouse analyses use only the quality-controlled files
`*_kept_tubules_labeled_fine.geojson`. In these files, structures with low support or
mis-assignment were removed upstream. Each remaining feature carries `kept_tubule`,
`doublet_tier1_flag`, `coarse_class` and `segment_class` attributes. The QC procedure that
produced these files is not in the repository [VERIFY: describe the rules and supply the code].
Only polygon geometry is used downstream. The fine labels in the GeoJSON are not used to label PT
segments (see *Integration, clustering and reviewed labels*).

**Centroid verification.** Every expression matrix records each structure's polygon index and
centroid at build time. Each analysis that joins expression back to a segmentation recomputes
the centroid of the indexed polygon and stops if any deviation exceeds 1 px. This guards against
pairing a matrix with a different version of the segmentation file.

## Tubule-by-gene matrices

Each 2-µm bin was assigned to a segmented structure if its centre, in full-resolution image pixel
coordinates (scale factor 1), was covered by the structure's polygon. The search used a shapely
STRtree spatial index, and each bin was assigned to at most one structure. Raw counts of the
assigned bins were summed per structure, giving one structure-by-gene count matrix per specimen.
Each matrix keeps the polygon index, centroid and number of contributing bins (`n_spots`) for every
structure. The segmentation polygons must be in the same full-resolution coordinate frame as
Space Ranger's tissue positions [VERIFY how the H&E image used for segmentation was registered to
the Visium HD image].

Matrix sizes (polygons; polygons with at least one bin): `Ctrl1A2` 7,624 (7,624); `Ctrl1A4` 7,002
(7,002); `HUK1_COR1` 8,747 (6,416); `HUK1_MED1` 14,186 (7,223). The human segmentations contain many
polygons that receive no bins [VERIFY: for example, a segmented area larger than the capture area].
These polygons are removed by the structure filters below.

## Cross-species ortholog space

Human genes were converted to mouse-symbol space using the HGNC Comparison of Orthology
Predictions (HCOP) fifteen-column human–mouse table [5]. The table was downloaded on [VERIFY date]
with checksum [VERIFY].

- *Which pairs were kept.* Pairs supported by fewer than three source databases were removed.
  Each remaining pair was scored as 10 × (number of supporting databases) + 1 if the two symbols
  are identical. A pair was kept only if each gene was the other's highest-scoring partner
  (mutual best). This gave a bijective map of 17,449 human–mouse pairs.
- *Mapping.* Human counts were moved to mouse-symbol columns through a sparse transformation
  matrix. No reviewed overrides were used.
- *Measured versus structural zeros.* The map creates a mouse column for every accepted pair,
  whether or not the human probe panel measured that gene. A column that no human feature feeds is
  a structural zero, not an observed zero. Each gene therefore carries a flag recording whether
  every input matrix measured it. Of the accepted pairs, 15,567 were measured in all four specimens.
- *Where the flag is applied.* Gene-level analyses are restricted to measured genes. For the
  integration and clustering step, the flag masks the analysis gene set but no columns are
  removed from the object. Removing columns would change Scanpy's Seurat-flavour HVG binning, and
  with it the embedding and the reviewed clustering.

## Structure filters and normalisation for integration

The integration and clustering below are reproduced from raw counts by notebook 13, which applies
notebook 03's settings.

**Structure filters.**

1. *Mouse.* Mouse structures come from the QC'd segmentations and get no gene-count threshold.
2. *Human.* Human structures need at least 100 detected ortholog genes.
3. *Low-support tail, per species.* Among structures passing rules 1–2, the lowest 5% of `n_spots`
   was removed separately within each species. We used a per-species quantile rather than one
   absolute threshold for two reasons: the species differ about two-fold in bins per structure,
   and an absolute threshold unbalances the species counts.

This kept 26,839 of 37,559 structures: `Ctrl1A2` 7,231; `Ctrl1A4` 6,676; `HUK1_COR1` 6,098;
`HUK1_MED1` 6,834.

**Gene filter.** Genes were kept if detected in at least 5% of retained structures and with at least
20 total counts. This left 10,076 genes, of which 9,943 were measured in all inputs.

**Normalisation.** Counts were scaled to 10,000 per structure over the retained genes and
log-transformed with log1p. Mitochondrial (`mt-`) and ribosomal (`Rpl`, `Rps`, `Mrpl`, `Mrps`) genes
stayed in the object but were excluded from integration features. The 133 retained genes that lack
a human feature cannot be selected as integration features, because selection requires
variability in both species.

## Integration, clustering and reviewed labels (pass 1)

**Integration features.** Highly variable genes (HVGs) for integration were selected separately
within each species. The selection was Seurat-flavour and batch-aware by specimen, with mean
between 0.0125 and 3 and normalised dispersion ≥ 0.5. Only genes selected in both species were
kept: 732 genes.

**Integration.**

- Principal component analysis (PCA): 50 components on the log-normalised values of these genes
  (Scanpy [6]; random state 0).
- Harmony [7] through rpy2, with R harmony 2.0.5 (enforced at run time). Batch variable = specimen;
  θ = 6, λ = 1, at most 30 iterations, τ = 0, R seed 0.

Sample correction and species are confounded. Harmony can therefore remove species biology along
with specimen effects; this is why the integrated embedding is used only to group and order
structures. All gene-level comparisons use expression rebuilt from raw counts (*PT analysis set*).

**Clustering.**

- A 30-nearest-neighbour graph was built on the 50 Harmony dimensions.
- Leiden clustering [8]: igraph implementation, resolution 0.7, 2 iterations, random state 0.
  This gave 11 clusters.
- The partition is pinned by a fingerprint: number of structures, resolution, neighbours, seed,
  number of clusters, and a SHA-1 hash of the cluster memberships. A run whose fingerprint differs
  stops before any labels are applied.

**Labels.** Cluster labels were assigned by manual review of marker detection and differential
expression. Labels came from a shared vocabulary and were assigned only at the granularity the
markers support. The reference panels were:

- PT-S1: Slc5a2, Slc5a12, Gatm, Lrp2, Cubn, Slc34a1
- PT-S2: Slc22a6, Slc13a3, Cyp2e1
- PT-S3: Slc22a7, Slc7a13, Cyp7b1, Slc6a18, Acsm3
- Thin limbs, thick ascending limb, distal convoluted tubule, connecting tubule and collecting
  duct, podocytes, vessel, stroma, smooth muscle and immune panels

The detection fractions used in the review are saved with the run. For the PT clusters, Slc5a2 was
detected in 0.91 of S1 structures, Slc22a6 in 0.94 of S2, and Slc22a7 and Slc7a13 in 0.74 and 0.78
of S3.

The final map was PT-S1, PT-S2, PT-S3, two ascending-limb clusters (AL), distal convoluted tubule
(DCT), two collecting-system clusters (CNT_CD and CCD), glomerulus, smooth muscle, and one
unresolved cluster. The unresolved cluster had flat panel scores and low gene counts and is
labelled `Unassigned`. No cluster had thin-limb markers.

Notebook 13's labels were checked against notebook 03's reviewed memberships. The two partitions
have different membership hashes and slightly different cluster sizes [WS1/coordinator: report
their concordance; see open items].

## Nephron re-integration (pass 2) and PT subset

Glomerulus, smooth-muscle and `Unassigned` clusters (2,504 structures) were removed. That left
24,335 tubular structures: PT 12,866; AL 6,731; CNT_CD 3,363; DCT 1,375. The first-pass labels
were kept fixed. Integration was repeated on these structures with the same settings: species-
intersected HVGs (752 genes), 50 PCs, Harmony, and a 30-neighbour graph.

The PT subset was taken after this second pass. It contains 12,866 structures:

| Specimen | S1 | S2 | S3 |
|---|---|---|---|
| `Ctrl1A2` | 938 | 976 | 938 |
| `Ctrl1A4` | 873 | 1,017 | 918 |
| `HUK1_COR1` | 1,890 | 908 | 640 |
| `HUK1_MED1` | 2,076 | 905 | 787 |

These are the S1/S2/S3 labels used by every downstream analysis. They come from clustering, so they
mark coarse regions rather than measured anatomical boundaries.

## PT coordinate reconstruction

> **[WS1] Placeholder.** This section describes the current coordinate. Workstream 1 may replace
> it with a problem-specific reconstruction. If so, rewrite this section and update every
> downstream count that depends on the coordinate: common support, the PT analysis set,
> transition positions, and the pathway lists.

**Fit.** The PT coordinate is a nonbranching principal curve fitted with scFates 1.2.5 [9], which
uses ElPiGraph [10].

- *Input space.* The first five pass-2 Harmony dimensions of the PT structures.
- *Curve.* `tl.curve` with 30 nodes and seed 0.
- *Topology check.* The fitted graph had to be a single path with two tips and no forks;
  otherwise the workflow stops.

**Orientation.** Orientation used a marker axis built from early PT genes (Slc5a2, Slc5a12, Gatm)
and late PT genes (Slc22a7, Slc7a13, Cyp7b1). Each gene's log-normalised expression was z-scored
within species. The axis is the mean late z-score minus the mean early z-score.

The curve was rooted at the tip with the highest early-marker score (the negative of the axis).
Pseudotime was computed with one seeded projection of structures onto the curve (`n_map = 1`).
It was then min–max scaled to [0, 1]. The marker axis sets only the root. It does not shape the
curve.

**Checks applied before the coordinate was saved.** In every specimen, the median coordinate of
S3 structures had to exceed the S1 and S2 medians. The observed medians were:

| | S1 | S2 | S3 |
|---|---|---|---|
| Mouse | 0.16–0.17 | 0.37–0.40 | 0.76 |
| Human | 0.19 | 0.46–0.47 | 0.77 |

**DPT comparator.** As a comparator, Scanpy diffusion pseudotime (DPT) [11] was fitted to the same
five dimensions. The settings were a 30-neighbour graph, a diffusion map, and a root structure. The
root was chosen among the 1% of structures most strongly assigned to the scFates root node, as the
one with the highest early-marker score.

Per-specimen Spearman correlations between the two coordinates were 0.88–0.91. A nonbranching
curve over the full tubular nephron was also fitted, as a topology diagnostic only. External
checks of the coordinate are described under *External segment-resolved datasets*.

## PT analysis set, common support and gene universe

**Common support.** Gene models were fitted on the intersection of the four specimens' 1st–99th
percentile coordinate intervals. In the current run this is [0.018, 0.864] on the scFates scale.
The interval keeps 12,272 of the 12,866 PT structures:

| Specimen | Structures | S1 | S2 | S3 |
|---|---|---|---|---|
| `Ctrl1A2` | 2,552 | 823 | 958 | 771 |
| `Ctrl1A4` | 2,685 | 821 | 1,014 | 850 |
| `HUK1_COR1` | 3,355 | 1,845 | 904 | 606 |
| `HUK1_MED1` | 3,680 | 2,041 | 902 | 737 |

**Evaluation grid.** Fitted curves were evaluated on 61 equally spaced positions. The grid runs
from the largest specimen minimum to the smallest specimen maximum. Within ±8% of the interval
around every grid position, each specimen had to contribute at least 15 structures.

**Expression rebuilt from raw counts.** For the retained structures, expression was rebuilt from
the original count matrices, starting from all 17,449 accepted ortholog pairs. The saved
coordinate is used as is and nothing upstream is refitted.

- *Library size.* Each structure's library size is the sum of counts over the 15,567 orthologs
  measured in every specimen. Expression is log1p(counts / library size × 10⁴).
- *Eligible genes.* A gene is eligible if it was measured in every specimen and its mean
  detection fraction across the four specimens (equal weights) is at least 2%. A gene may be
  undetected in one species and still be eligible. A gene missing from an input cannot be.

This gave 11,071 eligible genes. The normalisation denominator is set before the detection filter,
so changing the detection threshold does not change the scaling of the genes that remain.

## Nested gene models and species-by-position statistics

**Models.** For each eligible gene, three nested Gaussian regression-spline models were fitted to
log-normalised expression y as a function of coordinate s:

- M₀: y = f(s) + specimen contrasts
- M_level: y = f(s) + β·human + specimen contrasts
- M_full: y = f(s) + β·human + human × g(s) + specimen contrasts

**Basis.** f and g share a fixed cubic B-spline basis: 6 columns plus an intercept, with 3 internal
knots at the 25th, 50th and 75th percentiles of the pooled coordinate and boundaries at 0 and 1.
The basis is unpenalised. Because all three models use the same basis, they nest exactly.

**Specimen terms and weights.** Specimen intercepts enter as contrasts that sum to zero within
each species, so the species term stays identifiable. Models were fitted by weighted least
squares. Weights give each species half of the total weight and split it equally between that
species' specimens: wᵢ = n / (2 · J_s · n_j), where J_s is the number of specimens in species s and
n_j the number of structures in specimen j.

**Test statistics.** From weighted residual sums of squares (SSE), each nested comparison gives a
partial F statistic:

T = [(SSE_reduced − SSE_full) / Δp] / [SSE_full / (n − p_full)]

| Statistic | Comparison | Δp | Question |
|---|---|---|---|
| T_level | M₀ vs M_level | 1 | Constant species offset |
| T_spatial | M_level vs M_full | 6 | Change of the species difference with position |
| T_total | M₀ vs M_full | 7 | Either |

These F statistics treat structures as observations. They are used only to rank genes, never as
donor-level p-values.

**Local species difference.** M_full gives the fitted difference δ(s) = μ̂_human(s) − μ̂_mouse(s).
Its standard error is a heteroskedasticity-robust HC3 sandwich estimate [12]. HC3 is used because
the balancing weights are not inverse-variance weights. The working signal-to-noise ratio is
Z(s) = δ(s)/SE(s). A signed level effect, Z_level, is the human coefficient of M_level divided by
its HC3 standard error.

T_spatial measures changes in amplitude as well as changes in where a gene is expressed. It is not
a pure relocation statistic.

## Pathway libraries and the covariate-matched rank-AUC test

**Libraries.** Three libraries were used: Reactome 2022 (1,816 terms) [13], MSigDB Hallmark 2020
(50 terms) [14] and KEGG 2019 Mouse (304 terms) [15]. They were used in the Enrichr library
format [16]; the download source and date are not recorded [VERIFY].

**Membership.** Members already in mouse-symbol form were matched directly (case-insensitively).
Other members were mapped through the accepted ortholog table. A pathway was tested if it had
10–300 eligible members. That gave 1,513 pathways: Reactome 1,181, Hallmark 50, KEGG 282. The
statistic was chosen before any result was seen.

**Rank-AUC.** For each pathway P and each gene statistic T,

AUC_P = Pr(T_member > T_nonmember) + ½ · Pr(tie), with effect = AUC_P − 0.5,

computed over all eligible genes. Because the statistics are unsigned, enrichment means more
species difference, not higher expression in one species.

**Matched null.**

- Each gene was described by three covariates: equal-specimen mean log-normalised expression,
  equal-specimen detection fraction, and positional coverage. Positional coverage is the fraction
  of 10 equal-width coordinate bins in which the gene is detected, averaged across specimens.
- Each covariate was cut into tertiles, with tied values kept together. This gives up to 27 joint
  strata.
- For each pathway, 9,999 random gene sets were drawn (seed 12). Each set had the pathway's exact
  number of genes in every stratum, drawn without replacement from all eligible genes in that
  stratum.
- The one-sided empirical p-value is (1 + #{AUC_null ≥ AUC_observed}) / 10,000.
- Benjamini–Hochberg (BH) correction [17] was applied separately for each statistic, across all
  1,513 pathways from all three libraries.

**Discovery rule.** A spatial candidate has positive effect and empirical q ≤ 0.05 for T_spatial.
The same null with matching removed was run as a diagnostic.

The matched null addresses differences in how well genes are measured. It does not preserve
correlation between genes and does not supply donor replication.

## Correlation-adjusted competitive test

To check whether a pathway's rank evidence survives correlation among its genes, we computed a
second, CAMERA-style p-value [18] for T_spatial. The pathway score itself is unchanged.

1. Take residuals from M_full. Within each specimen, compute the mean Pearson correlation over all
   distinct pairs of member genes.
2. Average the four specimen means equally on the Fisher z scale, transform back, and set negative
   values to zero.
3. With m members and n = G − m non-members, U = m·n·AUC. Its variance is the CAMERA rank-sum
   variance with that correlation, multiplied by a tie correction, as in limma's
   `rankSumTestWithCorrelation`. Z uses a continuity correction of 0.5, and p is from the upper
   normal tail.
4. BH correction spans all 1,513 pathways.

Pathways whose correlation could not be estimated were marked unavailable; they were never assigned
zero correlation. This test is an exploratory adaptation to unsigned F statistics, not a direct
application of `camera` or `cameraPR`.

## Local divergence curves and positional descriptions

At each grid position, all eligible genes were ranked by |Z(s)| and by Z(s). For each pathway:

- Local divergence: D_P(s) = AUC_P(|Z(s)|) − 0.5
- Relative direction: S_P(s) = 2 · [AUC_P(Z(s)) − 0.5]

S is relative to the other genes, so absolute direction was read from the median member Z and the
fitted gene curves. The peak position is argmax D_P(s). Positional summaries are descriptive and
were computed after discovery: peak position, half-height width, and early, mid and late shares
over thirds of the interval. No significance is claimed at individual grid positions.

## Sensitivity runs and the final pathway list

> **[WS2] Placeholder.** Workstream 2 is finalising the reporting set and may change the filters
> or their order. The funnel below describes the current notebook 31 and 33 run.

**Seven planned sensitivity runs.** Every gene was refitted and the full pathway family retested in
each run, with the coordinate and grid fixed:

1. Omit `Ctrl1A2`, with the weights rebalanced over the remaining specimens.
2. Omit `Ctrl1A4`, likewise.
3. Omit `HUK1_COR1`, likewise.
4. Omit `HUK1_MED1`, likewise.
5. Spline basis of 4 instead of 6.
6. Spline basis of 8 instead of 6.
7. Detection threshold 5% instead of 2%. Membership, strata and background are updated.

A pathway is retained in a run if its effect stays positive with matched q ≤ 0.05, with the same
9,999 draws. The leave-one-section-out runs measure section sensitivity, not resampling of human
donors.

**Final list.** Spatial candidates passed three filters, in this order:

1. *Not specimen-sensitive.* Not called for T_spatial (positive effect, normal-approximation
   q ≤ 0.05) in either balanced relabeling or either within-species comparison (see *Specimen
   relabeling control*).
2. *Stable.* Retained in at least 5 of the 7 planned runs.
3. *Correlation-supported.* Correlation-adjusted q ≤ 0.05.

Pathways passing all three are called *confident* pathways.

**Labels added after the list was frozen.**

- *Robust*: the pathway is still called after both the registration check and the sampled-region
  check.
- *Core*: robust, and also called using only genes detected in both species, and also called after
  removing mouse sex-biased genes.
- *Themes*: assigned by fixed keyword rules on pathway names. They are used for reading only.

## Conventional comparators

Every comparator used the same structures, coordinate, gene universe and pathway family as the
primary analysis.

**Pseudobulk.**

- Raw counts were summed by specimen × reviewed segment. Whole-PT pseudobulks are the sum over the
  three segments.
- Genes were kept if they had at least 10 counts in at least 2 of the 4 specimen pseudobulks within
  each of S1, S2 and S3. That left 10,933 genes and 1,508 testable pathways (10–300 members).
- Differential expression used PyDESeq2 0.5.4 [19,20] with design `~species`, Wald tests and no
  Cook's-distance refit, for whole-PT and for each segment.

**Pathway tests on the pseudobulk.** Signed Wald statistics were tested with GSEApy 1.3.1 [21,22]
preranked GSEA: weight 1, multilevel p-values, seed 12, 10–300 members. BH correction spanned all
pathways for whole-PT, and all pathway × segment tests for the segment analysis. A call was q ≤ 0.05
in either direction; for segments, the minimum q across the three segments was used.

**Other screens.**

- *Pathway-score GAMs.* For each structure, the score was the mean of member genes' standardised
  expression, with the specimen balancing weights. The nested models were then fitted to the scores
  and tested with structure-level F distributions, followed by BH correction.
- *Step models.* Gene-level models that replace the smooth interaction with steps, tested with the
  matched rank-AUC. One version uses species × S1/S2/S3 indicators beyond the offset
  (T_discrete_spatial). The other uses species × 6 equal-width coordinate bins beyond the offset.
- *Other set tests on T_spatial.* T_spatial was also tested against size-only random sets and with
  unweighted preranked GSEA (weight 0, multilevel).

## Specimen relabeling control and the negative-control discovery ratio

A screen that reports pathways when there is no species difference to find is measuring specimen
heterogeneity or pseudoreplication, not species biology. To measure this, every strategy was rerun
on the two balanced relabelings of the four specimens:

- {`HUK1_COR1` + `Ctrl1A2`} versus {`HUK1_MED1` + `Ctrl1A4`}
- {`HUK1_COR1` + `Ctrl1A4`} versus {`HUK1_MED1` + `Ctrl1A2`}

This follows the random-partition control used in Lamian [23]. Structures, coordinate, weights,
genes, pathway family and test were unchanged.

In the relabeled gene models, each group has half the weight and each specimen a quarter. Specimen
intercepts enter as within-group contrasts, which absorb the species offset. Species × spline
basis terms enter every model as nuisance, so the tested term is a difference between specimens of
the same species. Relabeled DESeq2 models used the design `~species + group`.

Two within-species contrasts were also run with the same models: mouse versus mouse, and human
section versus human section.

**Calls and summary.**

- Matched-test calls in the specificity analysis used a normal approximation from the same matched
  null: z = (AUC − mean null AUC) / SD of null AUC, then BH within each partition and statistic. A
  call needed a positive effect and q ≤ 0.05. We used the approximation because the empirical
  p-value has a floor of 10⁻⁴.
- The negative-control discovery ratio (NCDR) is the mean number of relabeled calls divided by the
  number of species calls.
- A strategy was called specific if NCDR < 0.25 *and* the relabeled groupings averaged fewer than
  5% of tested pathways. The second condition was added after the first run, because a
  pathway-score model passed the ratio while calling about 280 pathways under relabeling. The
  change affects no gene-statistic decision.

With only two relabelings, this control screens for gross nonspecificity. It is not a calibrated
false-discovery rate.

## Direction-reversing genes

A gene was flagged as reversing if its fitted species difference had opposite signs at two places
along the PT. Each place needed at least 3 contiguous grid points with |Z| ≥ 3 and |δ| ≥ 0.25
log-normalised units. This is a display rule, not a test.

A reversal was confirmed without the coordinate if segment-level DESeq2 called the gene
significantly higher in human in one segment and significantly higher in mouse in another (both
padj ≤ 0.05).

## Robustness checks

**Registration.** A species shift in the coordinate would by itself produce position-dependent
species differences. Four checks address this:

1. *Segment transitions.* Within each specimen, a logistic regression (C = 10⁴) of the later
   segment label on the coordinate was fitted for S1/S2 and for S2/S3. The transition is the
   position where both labels are equally likely.
2. *Global registration.* We used 500 genes with clear positional structure in both species: the
   largest minimum across species of the curve's SD, among genes detected in at least 10% of each
   species' structures. Human curves were transformed by s → a + b·s, with a from −0.12 to 0.12 (25
   steps) and b from 0.85 to 1.15 (16 steps). For each transform we computed the mean Pearson
   correlation with the mouse curves over at least 80% overlap, on a 121-point grid.
3. *Anchored re-analysis.* Each human section was warped piecewise-linearly so that its S1→S2
   and S2→S3 transitions landed on the mean mouse transitions. The ends were fixed at 0 and 1, and
   mouse positions were unchanged. All genes were then refitted and retested against the original
   strata (9,999 draws, seed 12).
4. *Per-gene shift.* Each reversing gene's human curve was given its best shift (−0.15 to 0.15,
   31 steps, at least 70% overlap), chosen to minimise the variance of the difference curve. The
   reversal rule was then re-applied.

**Sampled region.** Glomeruli occur only in cortex, so distance to glomeruli serves as a depth
proxy. A mouse structure's depth was its mean distance to the three nearest glomerulus-labelled
structures in the same specimen. Glomerulus labels were taken from notebook 03's pass-1
annotation [see open items]. Centroids of PT structures and glomeruli were verified against the
fine GeoJSON before use.

A mouse S3 structure was called cortical-like if its depth was no greater than the 95th percentile
of depth among the same specimen's S1 structures. Deeper mouse S3 structures were removed, and the
models and pathway screen were rerun.

**Detection in both species.** A gene that one species barely detects behaves differently from a
probe-efficiency difference. A probe-efficiency difference shifts log expression by a constant,
which cancels in T_spatial; a near-absent gene does not cancel. The matched test was therefore
repeated on genes detected in at least 5% of structures in each species, with pathways re-formed
on those genes (10–300 members) and the original strata.

**Mouse sex bias.** All specimens are male. We used the per-segment female-versus-male mouse PT
differential-expression tables of Xiong et al. [24] (repository `LingyunXiong/Kidney_SexDiff`, MIT
licence; files `ArchRpeaks_WNN_PT-S{1,2,3}_fWTvsmWT_DEGs_DESeq2Norm_SSeq.csv`; commit [VERIFY]). A
gene was sex-biased if |log2 fold change| ≥ 1 and padj ≤ 0.05 in any PT segment. This gave 511
male-biased and 700 female-biased genes in the eligible universe. All sex-biased genes were removed
and the matched test was repeated.

**Trajectory method.** An earlier run of the same pathway analysis used notebook 03's PT
coordinate instead of the scFates coordinate. That run is notebook 12 logic version 5, with results
in `results/pt_pathway_remodeling/`. The notebook 03 coordinate is Scanpy DPT recomputed on the PT
subset of the pass-2 Harmony embedding, rooted near the low end of the S1→S3 marker axis and
oriented and scaled by that axis (5th–95th percentile). Its common support, gene universe (11,106
genes) and pathway family (1,516) differ slightly from the current run. Agreement between the two
runs was compared at the level of pathway calls [WS1: rerun under the final coordinate from a
committed notebook; see open items].

## External segment-resolved datasets

Our PT profiles were compared with four public datasets that measured PT segments directly.

**Datasets.**

- *GSE150338* [25]: microdissected mouse PTS1, PTS2 and PTS3; male C57BL/6 mice aged 6–8 weeks;
  bulk RNA-seq, TPM per replicate.
- *GSE56743* [26]: microdissected rat S1, S2 and S3; male Sprague-Dawley rats; bulk RNA-seq, RPKM
  per replicate.
- *CELLxGENE Census* [27], version 2025-11-08, accessed with `cellxgene_census` 1.18.0:
  - Mouse dataset `25818bf7-e2a7-41ec-8ff2-bc369c0ff4f5` [28]: 10x multiome nuclei annotated as PT
    segment 1, 2 and 3.
  - Human dataset `09b518f9-da64-44cc-aec8-70a89d55611f` [29]: renal-cortex benchmarking dataset
    annotated as proximal convoluted tubule or PT segment 3. No human study has microdissected
    S1–S3, so human positional data come from these annotations.

**Processing.**

- *Microdissection data.* Values were transformed to log2(x + 1) and averaged over replicates
  within each segment. "Early" is the mean of S1 and S2. Rat symbols were matched to our mouse
  symbols case-insensitively, which is an approximation.
- *Census data.* Primary-data cells of the PT cell types were extracted for the genes in our
  ortholog universe. Raw counts were summed per donor × sex × cell type × developmental stage.
  Each sum was normalised by the cells' total counts over all genes and transformed to
  log2(CPM + 1). Embryonic donors and donor × segment groups with fewer than 50 cells were dropped.
  A donor was kept only if it had every segment needed. This left 12 female and 12 male mice, and
  4 female and 2 male human donors. Mouse early PT was pooled from S1 and S2 counts.
- *Our data.* For each specimen, we averaged log-normalised expression within S1, S2 and S3, and
  within early PT (S1 and S2 structures pooled). Species values are the mean of specimen means.

**Comparisons.**

- *Zonation within each species.* Our late − early contrast (S3 minus early) was compared with each
  external contrast by Spearman correlation across genes. Genes had to be detected in at least 10%
  of our structures in that species and reach a mean log2 value ≥ 1 in either external group.
  Direction agreement was computed for genes with |external contrast| ≥ 1.
- *Species × position.* Our interaction, (human late − early) − (mouse late − early), was compared
  with the external interaction from human snRNA minus male-mouse snRNA. The comparison was
  repeated with female mice and with microdissected male mice as the mouse reference. Genes had to
  be detected in at least 10% of structures in both species and reach external human level ≥ 1.
  Results are reported for all such genes and for the top 10% by T_spatial. Direction agreement used
  genes with |external| ≥ 1 and |ours| ≥ 0.25.
- *Pathway level.* For each confident pathway with at least 5 members detected in both species, we
  compared the median member interaction. A pathway counted as informative when both values had
  |median| ≥ 0.1.

## Literature check

> **[WS3] Placeholder.** Workstream 3 is extending this check. The current procedure is described
> in `docs/results/pt-pathway-literature-novelty.md`.

Each confident pathway (with up to eight driver genes), each direction-reversing gene, and each
pathway called only by conventional screens was searched in titles and abstracts through Europe
PMC. Queries contained only gene, pathway and generic terms. Each finding was assigned a status:

- *Known*: the species difference and its PT position are both reported.
- *Partly known*: one part is reported.
- *Not found*: nothing after at least two targeted queries.
- *Contradicts*: the literature reports the opposite.

"Not found" means not found at the abstract level. It is not proof of novelty.

## Statistics and multiple testing

**Replication.** Two mouse specimens and two sections from one human donor provide no estimate of
variation among human donors. All inference is therefore conditional on the observed specimens and
is reported as descriptive.

**How each statistic is used.**

- Structure-level statistics (partial F, HC3 standard errors, Z(s)) rank and describe genes. They
  are not donor-level tests.
- Pathway-level tests are competitive tests of a gene statistic against matched genes. They are not
  tests of a species effect in a population.

**Multiple testing.** Correction used Benjamini–Hochberg [17] at q ≤ 0.05 in every family:

| Family | Scope of BH correction |
|---|---|
| Matched rank-AUC | Each statistic separately, across all pathways in the three libraries |
| Correlation-adjusted test | All pathways |
| Signed GSEA | Whole-PT across pathways; segments across pathway × segment |
| Specificity (normal approximation) | Each partition × statistic |
| Sensitivity runs | The full pathway family, rerun per run |

With 9,999 draws, the smallest empirical p-value is 10⁻⁴. Monte Carlo standard errors are
reported.

**Descriptive rules, not tests.** Themes, peak positions, direction-reversal flags, robustness
labels and external concordance.

**Two relabelings.** These can detect only gross nonspecificity.

**Sex.** Because all specimens are male, species and sex cannot be separated in our data. Female
external references were used only as a check.

## Software and reproducibility

**Analysis environment.** Python 3.11.16, Scanpy 1.11.5 [6], AnnData 0.12.19, NumPy 1.26.4, SciPy
1.17.1, pandas 2.3.3, statsmodels 0.15.0, patsy 1.0.3, scikit-learn 1.9.0, python-igraph 1.0.0,
leidenalg 0.12.0, scFates 1.2.5 [9], PyDESeq2 0.5.4 [20], GSEApy 1.3.1 [22], shapely 2.1.2, rpy2 3.6.7,
R 4.5.3 with harmony 2.0.5 [7].

These are the installed versions. Notebook 12's manifest records NumPy, SciPy, pandas, AnnData,
patsy, statsmodels and GSEApy at these versions. The stage-cache keys record Scanpy 1.11.5 and
harmony 2.0.5. Notebook 13 enforces scFates 1.2.5 and harmony 2.0.5 at run time.

**Segmentation environment.** Python 3.10 and PyTorch 2.5.1 (CUDA 12.1), with scikit-image [30],
OpenCV and OpenSlide [VERIFY versions used for the mouse inference run].

**Reproducibility.**

- Expensive stages are cached under keys built from parameters, input checksums and code digests,
  and the logic version of each notebook is recorded in its manifest.
- Every analysis is a committed, output-free notebook that takes `--data-root` and
  `--results-root`.
- Census extraction was run with a script that is not yet committed [see open items].

## Data and code availability

[VERIFY/coordinator: deposition plan for the Visium HD data (for example GEO, or controlled access
for the human tissue); segmentation polygons; code archive DOI and licence.] The public datasets
are listed above with their accession numbers.

---

## References

1. Oliveira MF et al. High-definition spatial transcriptomic profiling of immune cell populations in colorectal cancer. *Nat Genet* 57:1512–1523 (2025). https://pubmed.ncbi.nlm.nih.gov/40473992/
2. Oquab M et al. DINOv2: Learning Robust Visual Features without Supervision. arXiv:2304.07193 (2023). https://arxiv.org/abs/2304.07193
3. Kaplan et al. (Sophont). OpenMidnight model card and blog post, "How to Train a State-of-the-Art Pathology Foundation Model with $1.6k" (2025). https://huggingface.co/SophontAI/OpenMidnight [CITE: no peer-reviewed or arXiv record found; confirm the preferred citation]
4. Kirillov A et al. Panoptic Segmentation. CVPR 2019; arXiv:1801.00868. https://arxiv.org/abs/1801.00868
5. Yates B et al. Updates to HCOP: the HGNC comparison of orthology predictions tool. *Brief Bioinform* 22:bbab155 (2021). https://pubmed.ncbi.nlm.nih.gov/33959747/
6. Wolf FA et al. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol* 19:15 (2018). https://pubmed.ncbi.nlm.nih.gov/29409532/
7. Korsunsky I et al. Fast, sensitive and accurate integration of single-cell data with Harmony. *Nat Methods* 16:1289–1296 (2019). https://pubmed.ncbi.nlm.nih.gov/31740819/
8. Traag VA et al. From Louvain to Leiden: guaranteeing well-connected communities. *Sci Rep* 9:5233 (2019). https://pubmed.ncbi.nlm.nih.gov/30914743/
9. Faure L et al. scFates: a scalable python package for advanced pseudotime and bifurcation analysis from single-cell data. *Bioinformatics* 39:btac746 (2023). https://pubmed.ncbi.nlm.nih.gov/36394263/
10. Albergante L et al. Robust and Scalable Learning of Complex Intrinsic Dataset Geometry via ElPiGraph. *Entropy* 22:296 (2020). https://pubmed.ncbi.nlm.nih.gov/33286070/
11. Haghverdi L et al. Diffusion pseudotime robustly reconstructs lineage branching. *Nat Methods* 13:845–848 (2016). https://pubmed.ncbi.nlm.nih.gov/27571553/
12. MacKinnon JG, White H. Some heteroskedasticity-consistent covariance matrix estimators with improved finite sample properties. *J Econometrics* 29:305–325 (1985). https://doi.org/10.1016/0304-4076(85)90158-7
13. Gillespie M et al. The reactome pathway knowledgebase 2022. *Nucleic Acids Res* 50:D687–D692 (2022). https://pubmed.ncbi.nlm.nih.gov/34788843/
14. Liberzon A et al. The Molecular Signatures Database (MSigDB) hallmark gene set collection. *Cell Syst* 1:417–425 (2015). https://pubmed.ncbi.nlm.nih.gov/26771021/
15. Kanehisa M et al. New approach for understanding genome variations in KEGG. *Nucleic Acids Res* 47:D590–D595 (2019). https://pubmed.ncbi.nlm.nih.gov/30321428/
16. Kuleshov MV et al. Enrichr: a comprehensive gene set enrichment analysis web server 2016 update. *Nucleic Acids Res* 44:W90–W97 (2016). https://pubmed.ncbi.nlm.nih.gov/27141961/
17. Benjamini Y, Hochberg Y. Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. *J R Stat Soc B* 57:289–300 (1995). https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
18. Wu D, Smyth GK. Camera: a competitive gene set test accounting for inter-gene correlation. *Nucleic Acids Res* 40:e133 (2012). https://pubmed.ncbi.nlm.nih.gov/22638577/
19. Love MI et al. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol* 15:550 (2014). https://pubmed.ncbi.nlm.nih.gov/25516281/
20. Muzellec B et al. PyDESeq2: a python package for bulk RNA-seq differential expression analysis. *Bioinformatics* 39:btad547 (2023). https://pubmed.ncbi.nlm.nih.gov/37669147/
21. Subramanian A et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *PNAS* 102:15545–15550 (2005). https://pubmed.ncbi.nlm.nih.gov/16199517/
22. Fang Z et al. GSEApy: a comprehensive package for performing gene set enrichment analysis in Python. *Bioinformatics* 39:btac757 (2023). https://pubmed.ncbi.nlm.nih.gov/36426870/
23. Hou W et al. A statistical framework for differential pseudotime analysis with multiple single-cell RNA-seq samples. *Nat Commun* 14:7286 (2023). https://pubmed.ncbi.nlm.nih.gov/37949861/
24. Xiong L et al. Direct androgen receptor control of sexually dimorphic gene expression in the mammalian kidney. *Dev Cell* 58:2338–2358 (2023). https://pubmed.ncbi.nlm.nih.gov/37673062/
25. Chen L et al. A Comprehensive Map of mRNAs and Their Isoforms across All 14 Renal Tubule Segments of Mouse. *J Am Soc Nephrol* 32:897–912 (2021). https://pubmed.ncbi.nlm.nih.gov/33769951/
26. Lee JW et al. Deep Sequencing in Microdissected Renal Tubules Identifies Nephron Segment-Specific Transcriptomes. *J Am Soc Nephrol* 26:2669–2677 (2015). https://pubmed.ncbi.nlm.nih.gov/25817355/
27. CZI Cell Science Program et al. CZ CELLxGENE Discover: a single-cell data platform for scalable exploration, analysis and modeling of aggregated data. *Nucleic Acids Res* 53:D886–D900 (2025). https://pubmed.ncbi.nlm.nih.gov/39607691/
28. Chen S et al. Multi-omic and spatial analysis of mouse kidneys highlights sex-specific differences in gene regulation across the lifespan. *Nat Genet* 57:1213–1227 (2025). https://pubmed.ncbi.nlm.nih.gov/40259083/
29. Acera-Mateos M et al. Systematic evaluation of single-cell multimodal data integration enhances cell type resolution and discovery of clinically relevant states in complex tissues. *Genome Biol* 27:64 (2026). https://pubmed.ncbi.nlm.nih.gov/41821037/ The Census collection cites the preprint https://doi.org/10.1101/2025.03.06.637075, which has the same first author and dataset description [VERIFY that the published version is the right citation].
30. van der Walt S et al. scikit-image: image processing in Python. *PeerJ* 2:e453 (2014). https://pubmed.ncbi.nlm.nih.gov/25024921/
