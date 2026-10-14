# Making better use of a continuous axis: what the literature offers (review, 2026-10-14)

**Question (user, 2026-10-14).** Are we making the best use of the continuous PT coordinate: for visualising genes,
combining genes into pathways, classifying patterns, and continuous pathway analysis? Have trajectory-inference,
time-course or spatial-zonation papers already built better ways?

**Method.** Four parallel literature searches:
- A: gene dynamics;
- B: continuous pathway analysis;
- C: spatial zonation and gradients;
- D: comparing trajectories across conditions and species.

Every PMID was checked against PubMed (title, year, journal). The detailed tables, with about 140 papers and the
software, replicate handling and fit to our design of each, are in `docs/research/continuous-analysis/`:
- `a-gene-dynamics.md`;
- `b-continuous-pathways.md`;
- `c-spatial-gradients.md`;
- `d-comparing-trajectories.md`.

This page is the synthesis.

## Bottom line

**Yes, the field has better tools than we use. They fall into five families, and the four searches converged on the
same few ideas independently.**

**No published method fits our design as it stands.** Our design combines:
- a reconstructed one-dimensional anatomical axis;
- aggregated tubule units with an exposure offset;
- 2 vs 2 specimens plus one donor;
- species-specific probe panels.

So each idea has to be adapted. Our specimen-level calibration (balanced relabeling) stays.

**Two things the literature confirms we already do right:**
- Specimens, not units, are the replicates. Lamian (PMID 37949861) ran cell-level trajectory methods on a null
  dataset and got thousands of false positives (5,822–8,783 genes).
- Our logistic edge fit (position and 10–90% width) is the published switchde model (PMID 28011787).

## Convergent recommendations, ranked

### 1. Classify every gene by competing models, including a "step at the label boundaries" model (A, C, D)

**Sources.**
- Droin et al. 2021 (PMID 33432202): liver, nested count models with a per-mouse intercept, five classes by BIC
  weights.
- switchde (PMID 28011787) and ImpulseDE2 (PMID 30102402).
- GASTON (PMID 39849132): piecewise-linear expression with continuous gradients and jumps.

**What it is.**
- Fit a small ladder of NB models with specimen terms for each gene:
  1. flat;
  2. smooth gradient;
  3. sigmoid switch (position, width);
  4. transient bump (peak, width);
  5. a **step function at the known S1/S2 and S2/S3 transitions**.
- Assign classes with BIC weights, so each gene gets P(class).

**Why it matters for us.** This is the per-gene test of "continuous vs discrete": does a step at the labels explain
the gene as well as a smooth or free-position model? Part A found **no published method that asks this gene by
gene**. It also gives the transient class that our edge fits lack.

**Figures.**
- A class atlas (exemplar curves per class).
- A scatter of within-segment slope against boundary jump (GASTON-style), coloured by class.
- Heatmaps split into class blocks and sorted by onset or peak within each block. This removes the argmax artefact
  that piles monotone genes at the two ends.

### 2. Positional information in bits, per gene set, against the 1.58-bit ceiling of three labels (B)

**Sources.**
- Dubuis et al. 2013 (PMID 24089448) and Tkačik et al. 2015 (PMID 25361898): positional information and
  positional error.
- Cross-specimen decoding as in Phoenix.
- TIPS (Zheng et al. 2021, PMID 34370020) is a crude single-cell precedent. It rebuilds the trajectory from each
  pathway's genes and correlates the result with the global pseudotime. Its score is sign-agnostic, so it keeps
  pathways whose members move in opposite directions. It treats cells as replicates, and its reference axis shares
  genes with the pathways. Bits decoded across specimens from held-out genes are the calibrated version of the same
  question.

**What it is.**
- Decode position from a pathway's genes.
- Train on one specimen and test on the other; then train on mouse and test on human.
- Report bits and the positional error σx(x) along the PT.
- Three discrete labels carry at most log2 3 = 1.58 bits.

**Why it matters for us.**
- It is a quantitative, per-pathway "continuous beats discrete" number.
- **No application to an adult tissue axis was found.**
- It needs genes held out from building the curve, or count splitting (Neufeld 2023, PMID 36511385), to avoid
  circularity.

**Figures.**
- σx(x) curves for the leading pathways, with the label bands behind them.
- Bits per pathway with a reference line at 1.58.
- A transfer matrix (specimen to specimen, mouse to human).

### 3. Place pathways by where their genes change, not by a median curve (B, C)

**Sources.**
- GeneSwitches (PMID 32058565): a pathway timeline.
- PSEA (PMID 26955841): enrichment of positions against background.
- TcGSA (PMID 26111374): sub-trends inside a set.
- Halpern et al. 2017 (PMID 28166538) and Ben-Moshe et al. 2019 (PMID 31535084): centre of mass and enzyme order in
  liver.
- Hensman et al. 2013 (PMID 23962281): a hierarchical GP for pathway curves with gene and replicate levels.

**What it is.**
- For each pathway, test whether its members' onset, peak or centre-of-mass positions sit earlier or later than
  background, against a matched random-set null with specimen agreement.
- Before summarising a pathway into one curve, check for members that rise and fall in opposite directions, so the
  median does not cancel them.
- Optionally, fit a hierarchical pathway curve with a band and a coherence measure.

**Why it matters for us.**
- It fixes the weakness that our median program curve can flatten heterogeneous pathways.
- It answers "where", which is the paper's point, rather than only "whether".

**Figures.**
- A ridgeline pathway timeline (density of member onsets per pathway), with the label bands behind it.
- A forest plot of human − mouse or IR − control shifts.
- A sub-trend heatmap for heterogeneous pathways.

### 4. Locate where along the PT injury acts, and compare species without the probe-panel offset (C, D)

**Sources.**
- Weiss et al. 2026 (PMID 41721135): human liver zonation; per-gene slope with a patient random intercept; gradients
  lost in disease.
- Yakubovsky et al. 2026 (PMID 41986723): Visium HD human liver; per-sample centre of mass; human-vs-mouse ortholog
  scatter of zonation bias.
- OrganAxis (Yayon et al. 2024, PMID 39567784): shape similarity plus a bin × group interaction.
- Lamian (PMID 37949861): mean shift vs trend difference.
- DEtime / BGP: divergence points.
- Rayon et al. 2020 (PMID 32943498): a global time warp first, then gene-level shifts.

**What it is.**
- **Level-free summaries per gene and specimen:** slope, centre of mass, S3/S1 bias. These carry no probe-panel
  offset, which suits human vs mouse.
- **For IR:** a difference curve Δ(x) with a specimen-based band, plus the position where IR and control first
  diverge.
- **Before any expression comparison:** check unit density along the axis per specimen, to separate "IR removes S3
  units" from "IR changes S3 expression".

**Why it matters for us.**
- Part D found **no kidney study modelling AKI along a continuous anatomical PT axis**; liver has done it.
- Yakubovsky 2026 is the closest design precedent: Visium HD, per-sample centre of mass, human vs mouse.

**Figures.**
- A human-vs-mouse slope or centre-of-mass scatter, with quadrant counts (conserved S1, conserved S3, reversed).
- An effect-location heatmap of Δ(x), sorted by where |Δ| peaks.
- A density of divergence points along the PT.
- Unit-density ridges per specimen.

### 5. Genome-wide change density along the axis (A, B, C)

**Sources.**
- Trendy (PMID 30326833): breakpoint density.
- GAM derivatives with simultaneous bands (Simpson 2018, doi 10.3389/fevo.2018.00149).
- Adjacency-constrained clustering of ordered bins.

**What it is.**
- Count at each position the genes whose derivative band excludes zero, or whose edges fall there.
- Ask whether change concentrates at the S1/S2 and S2/S3 transitions or is spread along the PT.

**Figures.**
- A change-density track, with the label transitions, above a pathway × window enrichment heatmap.

### 6. Learn smooth programs across specimens and species, with warping (B)

**Sources.**
- MEFISTO (Velten et al. 2022, PMID 35027765): smooth factors along a covariate, per-group "sharedness" between
  groups, and built-in axis alignment.
- Alternatives: TEMPTED (a tensor across specimens) and CoGAPS/projectR (learn programs in mouse, project human).

**Use.** An exploratory complement to library pathways, annotated afterwards. Restrict to genes that are measured in
both species and held out from building the curve.

## Cautions the literature raises

- **Cell- or unit-level tests along trajectories over-call.** Lamian's null experiment shows this, as does Squair et al.
  2021 (PMID 34584091). Most spatially variable gene tests give inflated P values; SPARK and SPARK-X are the
  exceptions (benchmark in `c-spatial-gradients.md`).
- **Continuous models do not automatically beat pairwise tests.** On short series, time-course tools lost to pairwise
  tests (Spies et al. 2019, PMID 29028903). This bears directly on any "continuous beats discrete" claim: show it with
  calibrated comparisons, not by assumption. Our notebook 74 found continuous calling only a few more sets than steps.
- **Testing along an axis estimated from the same counts is circular** (Neufeld et al. 2023, PMID 36511385). Hold out
  genes (our fold protection) or split counts.
- **Medians and averages of set members can cancel** opposite sub-trends (TcGSA). Check within-set heterogeneity
  before summarising.

## Where our project would be new (gaps no paper fills)

1. **A continuous S1→S3 PT transcriptome atlas in tissue.** Kidney work is discrete (clusters, segments), spans the
   cortico-medullary axis, or is imaging-based. Notebooks 81, 83 and 84 already showed what continuous adds over
   these.
2. **AKI along a continuous anatomical PT axis.** Kidney AKI work uses regions, clusters or injury-state pseudotime.
3. **A per-gene test of continuous vs discrete against the known segment boundaries**, with specimens as replicates.
4. **Pathway positional information in bits on an adult tissue axis**, compared with the labels' ceiling.
5. **Small-sample calibration of curve comparisons** (2 vs 2 plus one donor). Our balanced relabeling has no
   published equivalent.
6. **Cross-species comparison under species-specific probe panels**, using level-free shape summaries.

## Suggested order of adoption

Each step needs a frozen protocol before any fit, under the project's rules.

1. **Recommendation 1** (model-ladder classes, including the boundary-step model). It reuses our NB machinery and
   directly tests the paper's core claim gene by gene.
2. **Recommendation 3** (pathway timeline and position enrichment, with sub-trend checks). It replaces the median
   program curve and gives the "where" figure.
3. **Recommendation 2** (positional information in bits). It is the most novel, and is quantitative for methods
   reviewers.
4. **Recommendation 4** (level-free human-vs-mouse summaries; IR effect-location map with unit-density check).
5. **Recommendations 5 and 6** as supporting analyses.
