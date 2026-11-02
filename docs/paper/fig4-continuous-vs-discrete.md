# Fig. 4: what clustering and bulk miss along the proximal tubule — protocol and results

Fig. 4 compares three summaries of the same control PT tubules: bulk (one PT mean per gene), discrete clustering
(the S1/S2/S3 means) and pseudostructure (a curve along the axis). It shows genes whose pattern the two discrete
summaries cannot express. The figure shows examples and does not headline a count; the count and its null are
reported here.

## Frozen protocol (2026-11-02, before any fit)

**Data.**
- **Fit and units.** The control/AKI package fit (`runs/aki`), only the control sections Ctrl1A2 and Ctrl1A4,
  and only units with a reference label PT-S1, PT-S2 or PT-S3.
- **Counts.** Raw counts from `02_qc/units_qc.h5ad`.
- **Library.** Each unit's library is its count total over non-mitochondrial genes.
- **Genes.** The 13 mitochondrial genes are left out, because they never entered the axis. A gene is screened if
  at least 10% of PT units in each control mouse have a count.

**Each gene gets an axis built without it.**
- **The five gene splits.** The Fig. 3d gene-half reruns split the 19,046 non-mitochondrial genes into halves A
  and B five times, with `default_rng([4, k])`, k = 0–4. Each half is a full rerun of stages 03–04.
- **Per-gene position.** A gene's position for a unit is the mean of the five axes built from the halves that
  exclude that gene. A unit needs a position on at least 3 of those 5 axes.
- **PT coordinate x.** x is the rank of that mean among the control PT units, pooled over both mice, divided by
  their number. It runs from 0 (start of PT) to 1 (end), in equal tubule counts.
- **Why.** No gene can look graded because it helped build its own axis.

**Boundaries and cores.**
- **Boundaries.** b12 is the threshold on x that misclassifies the fewest units between the S1 and S2 labels,
  pooled over both mice; b23 is the same for S2 and S3. They are computed for each axis.
- **Cores.** A segment's core is its x range more than 0.06 from any boundary. The transition window is within
  0.06 of a boundary.

**Models.** All models are NB GLMs with a log(library) offset, fitted separately in each mouse with unit weights:
- **Bulk:** an intercept.
- **Clusters:** the S1/S2/S3 label factor.
- **Clusters + curve:** the label factor plus a cubic B-spline in x with 6 df (knots at x quantiles).
- **Curve:** the spline alone. This is the pseudostructure summary, and it is evaluated on a grid of 201 points
  as log2 counts per 10k.

The dispersion is the ML estimate under clusters + curve, and is held fixed for the nested fit.

**Test beyond clusters.**
- **Statistic.** The likelihood ratio of clusters + curve against clusters.
- **Calibration.** A decoy shuffles x within each label within each mouse (seed 0), then recomputes the
  dispersion and the likelihood ratio.
- **Threshold.** For each mouse, the threshold is the smallest one at which (decoy genes above it) / (real genes
  above it) ≤ 0.05.
- **Pass.** A gene passes when it clears the threshold in both mice.

**Category 1: peak between segments.** The gene passes the test, and in both mice the curve satisfies all of:
- **Location.** Its global maximum over x lies in a transition window, at the same boundary in both mice.
- **Height.** The maximum is at least 0.5 log2 (1.41×) above the mean of the curve over each flanking segment's
  core.

**Category 2: change within a segment.** The gene passes the test, and in both mice:
- **Size.** Inside one segment's core, the curve's mean over the last third differs from its mean over the first
  third by at least 1 log2 (2×).
- **Direction.** The change has the same sign and is in the same segment in both mice.

**Neighbour filter.** A gene is dropped if its mean counts per 10k in any non-PT group of the control sections
(AL, DCT or CD) is more than twice its PT mean. This stops neighbouring-cell signal, such as thick limb beside late
S3, from posing as a PT pattern.

**Examples.**
- **Ranking.** Category 1 is ranked by the smaller of the two mice's peak heights. Category 2 is ranked by the
  smaller of the two mice's within-core changes.
- **Choice.** The figure shows the top three of each. A gene that qualifies for both categories is shown in
  Category 1.

**Figure.**
- **a.** One example drawn as bulk, three cluster steps, and the curve.
- **b.** A heatmap of every gene that passes either category: curves next to their three cluster means and their
  bulk mean, on the same per-gene scale.
- **c.** The Category 1 examples.
- **d.** The Category 2 examples.
- **Example panels.** These show both mice's curves, equal-count bin means of the data, and the cluster means
  drawn as steps.

**What would be a null.**
- **Category 2.** If no gene passes Category 2, the figure shows Category 1 only, and the null is reported here.
- **Category 1.** If no gene passes Category 1 either, the figure is not made in this form.
- **No threshold changes.** Neither threshold is changed after seeing results.

**Known limits.**
- **Position noise.** Positions inside a segment are noisy (notebook 72's C14 bound: 5 of 6 control cells are
  noise-dominated). This noise flattens real within-segment changes, so Category 2 is conservative.
- **No intervals.** The two mice are the replication; there are no per-gene intervals.

## Results

(appended when the screen finishes)
