# Making better use of a continuous axis: what the literature offers (review, 2026-10-14)

**Question (user, 2026-10-14).** Are we making the best use of the continuous PT coordinate: for visualising genes,
combining genes into pathways, classifying patterns, and continuous pathway analysis? Have trajectory-inference,
time-course or spatial-zonation papers already built better ways?

**Method.** Four parallel literature searches:
- A: gene dynamics;
- B: continuous pathway analysis;
- C: spatial zonation and gradients;
- D: comparing trajectories across conditions and species.

The first version (commit 734116d) was written partly from abstracts and partly from searches inside full texts. On
the user's instruction ("read the full paper for everything you find, abstract often isn't enough"), **every cited
paper was then read in full** on 2026-10-14, and each claim was checked against the text. The checks gave the same
answer in many places but changed several recommendations; the main changes are listed below.

Where the material lives:
- One note per paper, with its read status, what it actually does and a verdict on each claim the review makes about
  it: `docs/research/continuous-analysis/full-text-notes/`.
- The detailed tables (each entry tagged with its read status):
  - `a-gene-dynamics.md`;
  - `b-continuous-pathways.md`;
  - `c-spatial-gradients.md`;
  - `d-comparing-trajectories.md`.

This page is the synthesis.

## Read status (137 works)

- **Full text read: 118.** These are journal versions or author manuscripts, with supplements where they were open.
  Three were read as the authors' own versions: Aach & Church 2001, Goeman & Bühlmann 2007 and Chechik & Koller 2009.
- **Preprint read instead of the journal version: 8.**
  - GeneSwitches, tcgsaseq, Neufeld (count splitting), Phoenix;
  - Moor 2018, Ben-Moshe 2022, Saelens 2019, Kanton 2019.
- **Only partly read: 6.**
  - Rayon 2020 and Monocle 2014: the methods supplements were not accessible.
  - Tirosh 2016: the control-gene-set score is described only in the supplement, which could not be downloaded.
  - ImpulseDE 2017: the supplement was not accessible.
  - maSigPro 2006 and Heinonen 2015: read only through a web page summary, not verbatim.
- **Abstract only, full text not accessible: 5.** Claims about these rest on the abstract or on the software
  documentation, and are marked so.
  - Zaslaver 2004.
  - Harnik 2024 (plus its Extended Data legends).
  - Heitz 2025.
  - cellAlign (README and code).
  - Hotspot (documentation and code).
- **Software cited without a paper:** TCseq and SiZer, described from their documentation.

## What the full-text read changed

1. **Recommendation 1 (model-ladder classes).**
   - Droin 2021 fitted Gaussian mixed models to per-mouse layer means, not count models. Their BIC weights come from
     80 aggregated points per gene.
   - We must fit the ladder on per-specimen binned profiles, or calibrate the calls by specimen relabeling.
   - GASTON has no test for a jump. Our boundary-step model needs its own test.
   - Two papers weaken the novelty claim:
     - SPARK-X already tests genes for spatial signal left after conditioning on labels (one sample, unit-level).
     - Joost 2016 compared continuous and cluster versions of one NB model and found a tie (cells as units).
2. **Recommendation 2 (positional information in bits).**
   - The decoders in Dubuis, Tkačik and Petkova need a position measured independently of expression, many embryos at
     each position, and levels on one absolute scale. We have none of these.
   - Our "bits" would measure agreement with the reconstructed coordinate, not with anatomy.
   - "Train on mouse, test on human" is invalid under species-specific probe panels.
   - log2 3 = 1.58 is only an upper bound. The right comparator is the information the labels actually carry.
   - "Cross-specimen decoding as in Phoenix" was a wrong attribution. Phoenix cross-validates over cells within one
     dataset.
3. **Recommendation 3 (pathway position).**
   - Random-gene-set nulls are anti-conservative when set members are correlated. Final calls need specimen-level
     calibration.
   - The signed eigengene (PLAGE, WGCNA) does not cancel opposite-direction members, unlike our median. Hensman's
     hierarchical curve does cancel them.
   - The ordered-enzyme claim came from description, not a test. Chechik 2008 found backward orders too.
4. **Recommendation 4 (species and injury comparisons).**
   - Yakubovsky 2026 used standard Visium, not Visium HD (Visium HD was only a 3-donor validation). Zwick 2024 is an
     equal precedent.
   - Rayon 2020 did not test gene-level shifts after its global warp. Bar-Joseph 2003, Kanton 2019 and
     Cardoso-Moreira 2019 did.
   - Weiss 2026's slopes and OrganAxis's similarity are cell-level or not level-free.
   - The IR difference curve needs a condition term, and its band must be calibrated by relabeling.
5. **Recommendation 5 (change density).**
   - Penalised smooths blur steps, which biases the result toward "continuous" (Simpson 2018).
   - Trendy's breakpoints miss genes with a single gradient.
6. **Recommendation 6 (MEFISTO).**
   - Its warping removes timing differences between groups, so it must stay off for position comparisons.
   - It has no tests and no count model.
7. **Cautions.**
   - Count splitting is valid only for Poisson counts. Our held-out gene folds are the supported fix.
   - Spies 2019 showed time-course tools losing to pairwise tests only on 4-point step-like series.
   - Most spatially-variable-gene tests are miscalibrated in both directions, not only inflated.
8. **Our own methods.**
   - Our edge fit shares switchde's position and steepness parameters, but is not the same model.
   - Lamian already permutes sample labels. The paper states no minimum, but its smallest evaluated design is 4 vs 4.

## Bottom line

**Yes, the field has better tools than we use. They fall into five families. The four searches converged on the same
few ideas independently, and the full-text read kept all five, with changes to how each must be run.**

**No published method fits our design as it stands.** Our design combines:
- a reconstructed one-dimensional anatomical axis;
- aggregated tubule units with an exposure offset;
- 2 vs 2 specimens plus one donor;
- species-specific probe panels.

Every method we found either treats cells or units as replicates, or needs many more animals than we have. So each
idea has to be adapted, and our specimen-level calibration (balanced relabeling) stays the arbiter.

**Three things the literature supports in what we already do:**
- **Specimens, not units, are the replicates.**
  - Lamian (PMID 37949861) ran cell-level trajectory methods on a null dataset and got 5,822–8,783 false genes.
    Per-sample pseudobulk limma gave 0.
  - Zimmerman et al. 2021 found mixed models already inflated at 5 animals per group, and many genes failing to fit
    at 2. With few animals they advise a fixed effect per animal, and any permutation must shuffle whole animals,
    which our relabeling does.
- **Held-out gene folds against circularity.** Neufeld et al. (PMID 36511385) show gene splitting gives valid
  p-values. Training on one specimen and projecting another onto the curve ("cell splitting") is still circular.
- **Our logistic edge fit uses switchde's position and steepness parameters** (PMID 28011787; 10–90% width =
  ln 81/|k|).
  - switchde fixes the low state at zero, uses a Gaussian likelihood on log expression and treats cells as
    replicates.
  - Ours has a free baseline and amplitude, is fitted to binned profiles, and adds a two-edge model chosen by BIC
    (`pseudospace/transition_sharpness.py`).

## Convergent recommendations, ranked

### 1. Classify every gene by competing models, including a "step at the label boundaries" model (A, C, D)

**Sources.**
- **Droin et al. 2021** (PMID 33432202), liver:
  - nested Gaussian mixed models on log2 layer means (8 layers × 10 mice) with a mouse random intercept;
  - five space × time classes chosen by BIC, with BIC weights as class probabilities;
  - spatial shapes are polynomials of degree ≤ 2.
- **Other shape models:**
  - sigmoid: switchde (PMID 28011787);
  - impulse: ImpulseDE2 (PMID 30102402), whose Fig. 5 is a heatmap split by class and sorted by peak;
  - piecewise linear: GASTON (PMID 39849132). It tests only the slope within a domain; jump and slope calls are
    percentile cut-offs.
- **Published baseline to beat:** Monocle 3 (PMID 30787437) calls genes activated, repressed or transient by where the
  maximum of the scaled curve falls.
- **Published comparator:** SPARK-X (PMID 34154649) tests each gene for spatial pattern after adjusting for discrete
  labels; about half its hits disappeared after adjustment. It works in one sample, with unit-level permutation.

**What it is.**
- Fit a small ladder of models for each gene:
  1. flat;
  2. smooth gradient;
  3. sigmoid switch (position, width);
  4. transient bump (peak, width);
  5. a **step function at the known S1/S2 and S2/S3 transitions**.
- Fit on per-specimen binned profiles, as Droin fitted per-mouse layer means, or calibrate the class calls with
  balanced specimen relabeling. BIC weights from unit-level fits would sit near 0 or 1, because they treat units as
  replicates.
- The step-versus-smooth comparison needs our own test and calibration; no paper supplies one.
- **Positive controls for the step class:** genes with a "tight, non-overlapping boundary at the S2-S3 junction" by
  RNAscope (Ransick 2019, PMID 31689386). These are male-specific genes, and all our specimens are male.

**Why it matters for us.**
- It is the per-gene test of "continuous vs discrete": does a step at the labels explain the gene as well as a smooth
  or free-position model?
- **What is new:** doing this with specimens as replicates, an offset count model, and the anatomical segment
  boundaries.
- **What is not new:** testing for signal beyond discrete labels (SPARK-X), and comparing continuous with cluster
  models in one NB framework. Joost 2016 (PMID 27641957) did the latter and found the two explained the data equally
  well, with cells as units.
- It also gives the transient class that our edge fits lack.

**Figures.**
- A class atlas (exemplar curves per class).
- A within-segment slope against boundary jump scatter, coloured by class. This is our own figure; GASTON reports
  class counts.
- Heatmaps split into class blocks and sorted by onset or peak within each block, as in ImpulseDE2 Fig. 5. This
  removes the argmax artefact that piles monotone genes at the two ends.

### 2. Positional information per gene set, against the information the labels carry (B)

**Sources.**
- **Dubuis et al. 2013** (PMID 24089448), fly embryo gap genes:
  - continuous levels carried about 4.1 bits;
  - the best on/off code carried 2.92 ± 0.03 bits.
  - So discretising lost information.
- **Tkačik et al. 2015** (PMID 25361898): the optimal decoder and positional error σx(x). It stops at 4 genes, and
  even that needed Gaussian approximations, sample-size extrapolation and a regularised covariance.
- **Petkova et al. 2019** (PMID 30712870): the same decoder in mutants. It requires:
  - position measured independently of expression;
  - many embryos at each position: the decoder used 38 wild-type embryos, and the noise covariances came from data sets
    of 24–102 embryos;
  - levels on one absolute scale. Rescaling profiles to 0–1 made predictions much worse.
  - Collapsing graded profiles into on/off domains lost precision.
- **Phoenix** (Halperin 2026, PMID 42373304; bioRxiv preprint read):
  - predicts pseudotime from a pathway's genes with a random forest;
  - compares the accuracy with 200 random gene sets of the same size, so it is a competitive test;
  - cross-validates over cells within one dataset.
- **TIPS** (Zheng et al. 2021, PMID 34370020) is a cruder single-cell precedent.
  - It rebuilds the trajectory from each pathway's genes and correlates it with the global pseudotime.
  - The score ignores direction, so pathways whose members move oppositely are kept.
  - It treats cells as replicates, and its reference axis shares genes with the pathways.
  - Its simulation favours the method by construction, random gene lists also reach high |r|, and results depend on
    the trajectory algorithm.

**What it is.**
- Within one species, decode the reconstructed coordinate from a pathway's **held-out genes**. Train on one specimen
  and test on the other; this is our adaptation, since none of the sources does it.
- Report decoding error, and bits where they can be estimated, relative to the consensus coordinate. This is not
  anatomical position: our coordinate is built from the transcriptome.
- **Comparators:**
  - the information the S1/S2/S3 labels actually carry about the same coordinate, estimated from the same units;
  - a discretised (on/off) version of the pathway, as in Dubuis.
  - log2 3 = 1.58 bits is only an upper bound for three labels.
- Joint bits for 20–200-gene sets would be strongly biased upward, so use decoding error across specimens rather than
  a plug-in joint information.
- **Cross-species transfer** needs level-free features, such as within-specimen standardised or ranked expression. A
  decoder trained on mouse levels cannot be applied to human levels measured on another probe panel.

**Why it matters for us.**
- It is a quantitative, per-pathway "continuous vs discrete" number. Dubuis and Petkova show the principle, that
  discretising loses positional information, in embryos.
- **No application to an adult tissue axis was found.**
- **Caveat:** in the embryo papers, noise is variability across embryos at a known position. Ours would be
  unit-to-unit noise around a reconstructed coordinate.

**Figures.**
- σx(x) curves for the leading pathways, with the label bands behind them.
- Information per pathway, with reference lines for the labels' measured information and the 1.58 upper bound.
- A specimen-to-specimen transfer matrix.

### 3. Place pathways by where their genes change, not by a median curve (B, C)

**Sources (set position).**
- **GeneSwitches** (PMID 32058565; preprint and README read; the journal version is paywalled): pathway time is the
  median of its members' switch times, shown as a ridge plot. That rule is stated in the README, not the preprint.
- **PSEA** (PMID 26955841): tests a set's member positions against background.
  - Its benchmark shows fixed windows plus Fisher's test miss sets centred near bin edges. This is a ready
    "continuous beats discrete" precedent.
- **Chechik & Koller 2009** (PMID 19193146): compares each GO category's onset times with its sibling categories
  (Wilcoxon), on a linear axis.
- **Ben-Moshe et al. 2019** (PMID 31535084): mean centre of mass of 146 miR-122 targets against 1,000 random sets of the
  same size (p = 0.018).
- **Halpern et al. 2017** (PMID 28166538):
  - sorts genes by centre of mass;
  - shows bile-acid enzymes peaking in reaction order, descriptively;
  - its formal test is a discordance score for connected enzyme pairs against random pairs, which ignores order.
- **Chechik et al. 2008** (PMID 18953355): found forward activation most often, but also backward orders depending on
  condition. A reaction-order test along the PT should be two-sided, shuffling profiles on a fixed reaction graph.
- **Timing by pathway:**
  - DEtime ranks GO terms by when half their genes have diverged.
  - Waterfall times pathways by averaging members' on/off states.
  - Pseudotimecascade runs GO enrichment on successive switch-point groups, within each direction class.

**Sources (summarising a set without cancellation).**
- **Signed eigengene** (PLAGE, PMID 16156896; WGCNA, PMIDs 19114008 and 18031580):
  - genes are standardised before the SVD, so weights follow co-variation, not variance;
  - WGCNA fixes the sign by correlation with the set mean; PLAGE's sign is arbitrary, so we must orient it the same way;
  - rising and falling members load with opposite signs and do not cancel. A second, independent sub-trend goes to
    the next component.
- **TcGSA** (PMID 26111374): its sub-trend step clusters the members' fitted curves, picks the number with the gap
  statistic and shows one heatmap row per trend.
  - Borrow only that step.
  - Its test is self-contained, and member heterogeneity alone makes a set significant. FUNNEL's simulation found it
    calling all 90 sets, including 72 with no signal.
- **Hensman et al. 2013** (PMID 23962281): a hierarchical GP with cluster, gene and replicate levels. Its top level is
  an inferred cluster, not a pathway. Its cluster curve follows the members' mean, so opposite members still cancel.

**What it is.**
- For each pathway, test whether its members' onset, peak or centre-of-mass positions sit earlier or later than
  background.
  - Use random sets matched on size and expression for the statistic.
  - Do not take its gene-sampling p-value as final: it is anti-conservative when members are correlated (PSEA, CAMERA).
  - Final calls need agreement between specimens or balanced relabeling.
- Before summarising a pathway into one curve, check for sub-trends. Then summarise with a signed eigengene, or with
  medians split by direction, instead of the plain median.

**Why it matters for us.**
- It fixes the weakness that our median program curve can flatten heterogeneous pathways.
- It answers "where", which is the paper's point, rather than only "whether".

**Figures.**
- A ridgeline pathway timeline (density of member onsets per pathway), with the label bands behind it.
- A forest plot of human − mouse or IR − control shifts.
- A sub-trend heatmap for heterogeneous pathways.
- A temporal-GO heatmap: GO terms × successive switch-point groups.

### 4. Locate where along the PT injury acts, and compare species without the probe-panel offset (C, D)

**Sources (species).**
- **Yakubovsky et al. 2026** (PMID 41986723):
  - standard Visium of human liver; Visium HD was only a 3-donor validation;
  - per-sample centre of mass;
  - a human-vs-mouse ortholog scatter of a within-species log end-ratio (Spearman 0.15).
- **Zwick et al. 2024** (PMID 38321203): a human-vs-mouse scatter of per-gene centre of mass along the gut, with 2 mice
  and 2 donors (Spearman 0.29).
- Neither study faces a species-specific probe panel or models specimens.
- **Aizarani et al. 2019** (PMID 31292543):
  - 68% of zonated hepatocyte genes correlated positively between human and mouse;
  - but the same pipeline matched the mouse's own landmark-based profiles for only 80% of genes.
  - So a species comparison needs a within-species baseline computed the same way.
- **Weiss et al. 2026** (PMID 41721135):
  - per-protein linear slopes on robust-scaled intensities, with a patient random intercept only, so slope tests
    treat cells as replicates;
  - a noisier condition shrinks the slope and looks like "gradient lost";
  - a linear slope misses genes that peak mid-axis.
- **OrganAxis** (PMID 39567784):
  - its cosine similarity is computed on uncentred median profiles, so it is not level-free;
  - its bin × age ANOVA has no donor term.

**Sources (warping before gene tests).**
- **Precedents for "one global warp, then gene-level tests":**
  - Bar-Joseph 2003 (PMID 12934016): a global linear shift-and-stretch, then a per-gene test.
  - Kanton 2019 (PMID 31619793; preprint read): dynamic time warping, then a spline F-test on the aligned axis, after
    trimming to the shared extent.
  - Cardoso-Moreira 2019 (PMID 31243369): a global stage map, then local deviations.
- **Rayon et al. 2020** (PMID 32943498): one global mouse-to-human time-scaling factor (about 2.5) that genes and
  clusters agree with. It does not analyse gene-level shifts after it. Its methods supplement was not read; the
  time-factor details come from a summary of the preprint, not verified verbatim.
- **Aach & Church 2001** (PMID 11395426): warps fitted from small gene groups are unstable. Per-gene shifts need a
  stability check.
- **For us**, a uniform stretch cannot cause shifts on a per-species normalised axis. The risks are a non-uniform warp
  (different S1/S2/S3 proportions) or a partial one (different extent covered).

**Sources (injury).**
- **Lamian** (PMID 37949861):
  - splits a mean shift from a trend difference;
  - is Gaussian on imputed log expression;
  - permutes sample labels, with 4 vs 4 as its smallest design.
- **Divergence points:** DEtime (PMID 27288495) and BGP (PMID 29843817) assume identical curves before the divergence.
  - They apply directly to IR vs control; for species they need centred curves first.
  - They are one-sided, so run them from each end of the PT.
  - DEtime fits better than BGP when the condition labels are known.
- **TrendCatcher** (PMID 35175937) has a per-interval group-label permutation test of where two curves part.
- **HGAM** (Pedersen 2019, PMID 31179172):
  - a difference curve needs a condition-by-position smooth (`s(x, by = condition)`) plus specimen smooths;
  - its intervals depend on smoothing parameters estimated from 2 specimens, so they must be calibrated;
  - NB families work in `gam()`/`bam()` but not `gamm()`, and Pedersen advises against likelihood-ratio p-values
    between model structures.
- **Unit density along the axis:**
  - Ben-Moshe 2022 (PMID 35659879; preprint read) bins cells into zone tertiles. But it plots the coordinate's
    distribution at each time point, and compares injured cells with controls resampled to the same distribution.
    It also shows injury reprogramming the landmark genes that build the coordinate.
  - Milo (PMID 34594043) fits an NB model to cell counts in each neighbourhood of the expression graph, with samples
    as replicates and a weighted spatial FDR. Applying it to windows along our axis is our adaptation.
  - condiments (PMID 38280860, supplement) has a liver case where a treatment shifted hepatocytes along the
    central–portal axis: a change in where cells sit, not in their expression.
- **Tritschler et al. 2019** (PMID 31249007): projecting a new condition onto an embedding built from controls can hide
  states that lie off it.

**What it is.**
- **Species:**
  - compute level-free summaries per gene and specimen: slope with our offset, centre of mass, S3/S1 bias;
  - test them with balanced relabeling;
  - report mouse-vs-mouse agreement, computed the same way, next to human vs mouse.
- **IR:**
  - a difference curve Δ(x) from a condition-by-position model with specimen terms;
  - its band calibrated by balanced relabeling;
  - the position where IR and control first diverge.
- **Before any IR expression comparison:**
  - count units per specimen in overlapping windows, to separate "IR removes S3 units" from "IR changes S3 expression";
  - compare at matched coordinate distributions;
  - check IR unit positions on a coordinate built without injury-responsive genes, and report each IR unit's distance
    from the curve.

**Why it matters for us.**
- **No kidney study models AKI along a continuous anatomical PT axis.**
  - Kirita 2020 uses clusters and injury-state trajectories, and reports failed-repair PT cells scattered throughout
    the PT. That is a prediction we can test along the axis.
  - Melo Ferreira 2021 and Dixon 2022 have one mouse, or one section, per model or time point.
  - Lake 2023 builds injury-state pseudotime.
- Liver has done this, but in a coarse form: Ben-Moshe 2022 uses three zones with cell-level tests, and Weiss 2026
  uses cell-level slopes.

**Figures.**
- A human-vs-mouse slope or centre-of-mass scatter, with quadrant counts and the mouse-vs-mouse baseline.
- An effect-location heatmap of Δ(x), sorted by where |Δ| peaks.
- A density of divergence points along the PT.
- Unit-density ridges per specimen.

### 5. Genome-wide change density along the axis (A, B, C)

**Sources.**
- **Simpson 2018** (doi 10.3389/fevo.2018.00149):
  - derivatives of GAM fits with simultaneous intervals;
  - it warns that a GAM with one smoothness parameter is ill-suited to step-like change, and that ignoring
    dependence makes periods of change anti-conservative;
  - it also propagates uncertainty in an estimated axis, by refitting on 100 draws of that axis. That is a recipe for
    our elastic curve.
- **Trendy** (PMID 30326833): its breakpoints are changes of slope. A gene with one gradient along the whole PT has no
  breakpoint.

**What it is.**
- Count, at each position, the genes whose derivative band excludes zero, with bands from specimen-level variation.
- Report separately the share of genes with a single gradient across the PT; they change everywhere but have no edge.
- Before concluding that change is spread along the PT, compare with an explicit step model (Recommendation 1) or
  use adaptive smoothing. Penalised smooths spread sharp boundaries over neighbouring positions.

**Figures.**
- A change-density track, with the label transitions, above a pathway × window enrichment heatmap.

### 6. Learn smooth programs across specimens (B), exploratory only

**Sources.**
- **MEFISTO** (Velten et al. 2022, PMID 35027765):
  - smooth factors along a covariate;
  - for each factor, a group-by-group correlation matrix and a "sharedness" score;
  - optional axis alignment.
- **Limits:**
  - no tests;
  - Gaussian likelihood, with no count model or offset;
  - group-wise centring hides level changes.
- **Its alignment removes timing differences between groups** (supplement Sec. 8.3).
  - Keep it off when comparing positions.
  - Use it only to estimate a global map between axes.
- **Alternatives:**
  - CoGAPS/projectR (PMID 31121116) fits least squares on the target's levels, so probe-panel offsets enter the
    projected weights. Centre each gene within each specimen first, and compare shapes within a target.
  - TEMPTED (PMID 39696594) gives each specimen only a scalar loading on a shared curve, so it cannot show a shift in
    where a program changes. It is a weak fit, kept only as a minor alternative.

**Use.** An exploratory complement to library pathways, annotated afterwards. Restrict it to genes measured in both
species and held out from building the curve.

## Cautions the literature raises

- **Cell- or unit-level tests over-call.**
  - Lamian's null experiment shows this along trajectories.
  - Squair et al. 2021 (PMID 34584091) show it for between-condition tests, including spot-level spatial tests and
    the GO enrichments built on them.
  - Mixed models are not the fix at our sample size (Zimmerman 2021).
- **Spatially-variable-gene tests are miscalibrated in both directions.**
  - In one benchmark (Li 2025, PMID 40968359), six methods, including nnSVG and SpatialDE, were too conservative, and
    four were inflated. SPARK and SPARK-X were calibrated.
  - The issue that matters for us is that none of them models specimens.
- **Continuous is not automatically better.**
  - Spies et al. 2019 (PMID 29028903): on 4-point step-like series, pairwise tests beat most time-course tools, mainly
    because those tools gave false positives. On 8 or more points, spline tools beat pairwise tests.
  - Joost 2016 found a tie between continuous and cluster versions of one model.
  - Weiss 2026's "397 more proteins with a continuous test" compares a 1-df test with a 19-df cell-level ANOVA,
    without calibration.
  - Our axis is dense, so the lesson is calibration: compare continuous and S1/S2/S3 tests under the same specimen
    relabeling. Notebook 74 found continuous calling only a few more sets than steps.
- **Smoothing biases toward "continuous".**
  - Penalised smooths blur steps (Simpson 2018).
  - Reconstruction methods that average each unit over nearby positions make any step look gradual. Hinze 2021's
    "smooth transition" in the collecting duct comes from such averaged profiles.
  - Any "gradual, not stepped" claim needs an explicit step model as the comparator.
- **Testing along an axis estimated from the same counts is circular** (Neufeld et al., PMID 36511385).
  - Hold out genes (our fold protection).
  - Count splitting is valid only for Poisson counts; with NB overdispersion and large means the two halves stay
    correlated.
  - Training on one specimen and placing another on the curve does not remove the circularity.
- **Medians and averages of set members can mislead.**
  - Opposite sub-trends cancel (TcGSA).
  - GSVA's default score cancels mixed-direction sets.
  - Spectra (PMID 37735262) shows a plain set average is biased when sets overlap.
  - Check within-set heterogeneity before summarising.

## Where our project would be new (gaps no paper fills)

1. **A continuous S1→S3 PT transcriptome atlas in tissue.**
   - Kidney work is discrete (clusters, segments), spans the cortico-medullary axis, or is imaging-based.
   - Hinze 2021 ordered collecting-duct cells, and found the PT carried the least regional signal.
   - Rood 2019 describes an expression-learned "pseudospace" in general, but not for repeated tubule units.
   - Notebooks 81, 83 and 84 already showed what continuous adds over these.
2. **AKI along a continuous anatomical PT axis.** Kidney AKI work uses regions, clusters, single-animal sections or
   injury-state pseudotime.
3. **A per-gene test of continuous vs discrete against the known segment boundaries, with specimens as
   replicates.** This is qualified: SPARK-X conditions on labels in one sample, and Joost 2016 compared continuous and
   cluster models with cells as units.
4. **Positional information per pathway on an adult tissue axis.** It must be stated relative to the reconstructed
   coordinate, and compared with the information the labels carry.
5. **Small-sample calibration of curve comparisons (2 vs 2 plus one donor).** Lamian permutes sample labels, but its
   smallest evaluated design is 4 vs 4, and nothing is calibrated at 2 vs 2. Balanced relabeling there is our own.
6. **Cross-species comparison under species-specific probe panels.** It uses level-free shape summaries, with a
   within-species baseline. The closest precedents (Yakubovsky 2026, Zwick 2024) do not face a panel offset.

## Papers to read first (user request, 2026-10-14)

These are ranked by how close each is to our design, or by how much we would borrow from it. All of them were read in
full, except where a preprint is noted. Each note in `full-text-notes/` has the details.

**Closest to what we do.**
1. **Joost et al. 2016, *Cell Systems*** (PMID 27641957).
   - Like ours, its spatial axis (hair follicle) is built from expression alone.
   - It asks "clusters or continuous axes?" with one NB model, and finds a tie in explained molecules. Cells were
     pooled from 19 mice.
   - This is the benchmark our "continuous beats discrete" claim has to clear.
   - Its Fig. 4E (density of each population along the axis) is a figure to copy.
2. **Droin et al. 2021, *Nature Metabolism*** (PMID 33432202). The template for Recommendation 1:
   - nested models and five classes with BIC weights;
   - one intercept per mouse;
   - a condition (time of day) that varies between animals, as IR does for us.
3. **Hinze et al. 2021, *JASN*** (PMID 33239393).
   - The closest kidney precedent: a 1-D axis built from expression (collecting duct), a "continuous, not subtypes"
     conclusion, and a knockout read as a flattened gradient.
   - Its statistics are weak (n = 2; P values over reconstructed positions).
   - It warns that ordering by similarity may fail under ischaemia.
4. **Zwick et al. 2024, *Nature Cell Biology*** (PMID 38321203).
   - Our scale (2 mice, 2 human donors) and our question, along the small intestine.
   - It compares human and mouse by per-gene centre of mass.
   - It concludes "five discrete domains" from a gap statistic that was never tested against a smooth-gradient null.
     Reviewers may cite it against us, so read it closely.
5. **Weiss et al. 2026, *Nature Metabolism*** (PMID 41721135).
   - Human liver zonation in health and in disease (spatial proteomics), with per-protein linear slopes, and gradients
     lost in disease.
   - It is the nearest analogue of "IR flattens the S1-high program".
   - Take its figures, not its cell-level slope tests.

**Ideas to take for specific recommendations.**
6. **Kanton et al. 2019, *Nature*** (PMID 31619793; bioRxiv preprint read). The template for human vs mouse along one
   axis:
   - trim to the shared extent;
   - align;
   - spline test of shared vs species-specific curves;
   - a same-species split to estimate false positives, close to our relabeling.
7. **Ben-Moshe et al. 2022, *Cell Stem Cell*** (PMID 35659879; bioRxiv preprint read). Injury along an axis:
   - plot the coordinate's distribution per condition;
   - compare injured units with controls matched on position;
   - beware injury rewriting the genes that build the coordinate.
8. **Hou et al. 2023 (Lamian), *Nature Communications*** (PMID 37949861).
   - The null experiment that supports specimens as replicates.
   - The mean-shift vs trend split that we can rebuild in our NB framework.
9. **Dubuis et al. 2013, *PNAS*** (PMID 24089448) **and Petkova et al. 2019, *Cell*** (PMID 30712870).
   - Positional information in bits.
   - Continuous levels against an on/off code.
   - Why their decoder needs things we lack.
10. **Yakubovsky et al. 2026, *Nature*** (PMID 41986723). The figure vocabulary for species comparison:
    - centre-of-mass-sorted heatmaps;
    - ortholog scatter with conserved and reversed anchor genes;
    - an honest, weak cross-species correlation.
11. **PSEA (Zhang et al. 2016, *J Biol Rhythms*, PMID 26955841) and Chechik & Koller 2009 (*J Comput Biol*, PMID
    19193146).** Short reads for Recommendation 3:
    - testing set positions against background;
    - a benchmark showing fixed windows miss sets near window edges;
    - "edge = position of half amplitude".
12. **GASTON (Chitra et al. 2025, *Nature Methods*, PMID 39849132).** Its "slope inside a segment vs jump at a
    boundary" split maps onto S1/S2/S3. We must add the jump test it lacks.

**For defending the analysis to reviewers (skim).**
- Spies 2019 (continuous vs pairwise benchmark).
- Neufeld et al. (circularity; gene splitting).
- Zimmerman 2021 and Squair 2021 (replicate unit).
- Simpson 2018 (derivatives; steps; axis uncertainty).
- SPARK-X (the label-adjusted comparator).
- Aizarani 2019 (human vs mouse liver; within-species baseline).
- Ransick 2019 and Kirita 2020 (kidney priors and positive controls). A marker counts as a control only if it is on
  our mouse panel and the time point is comparable; Havcr1, Lcn2, Mki67 and Top2a are reported absent from it.

## Suggested order of adoption

Each step needs a frozen protocol before any fit, under the project's rules.

1. **Recommendation 1** (model-ladder classes, including the boundary-step model). Fit on per-specimen binned
   profiles, and calibrate the calls by relabeling. It reuses our NB machinery and tests the paper's core claim gene
   by gene.
2. **Recommendation 3** (pathway timeline and position enrichment, with a sub-trend check and a signed-eigengene
   summary). It replaces the median program curve and gives the "where" figure.
3. **Recommendation 2** (within-species decoding of the coordinate from held-out pathway genes). It is the most
   novel, and gives a quantitative number for methods reviewers.
4. **Recommendation 4** (level-free species summaries with a mouse-vs-mouse baseline; an IR effect-location map after
   the unit-density check).
5. **Recommendations 5 and 6** as supporting analyses.

## Run in notebook 90 (2026-10-14)

Recommendations 1–5 and the IR part of recommendation 4 were run in notebook 90, one figure each. Protocol and results
are in `docs/results/pt-continuous-axis-experiments.md`. Recommendation 6 (MEFISTO) was not run: it would add a
dependency for an exploratory method with no test, and notebook 89 already learned and replicated smooth programs.
