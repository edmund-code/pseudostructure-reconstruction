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

### Screen (2026-11-02)

`analysis/paper/fig4_screen.py`, output in `results/paper_figures/fig4/screen/`. Inputs: 2,780 (Ctrl1A2) and
2,745 (Ctrl1A4) PT units, 8,039 screened genes, and 32 gene-specific axes. On every axis b12 is 0.29–0.30 and b23
is 0.67–0.68.

| Step | Ctrl1A2 | Ctrl1A4 | Both |
|---|---|---|---|
| Decoy threshold (LR) | 16.0 | 15.8 | |
| Curve beyond clusters | 2,374 | 2,442 | **1,680** (21% of screened genes) |
| Dropped by the neighbour filter | | | 614 of the 1,680 |
| Category 1, peak between segments | | | **5**, all at S2/S3 |
| Category 2, change within a segment | | | **74** (S1 3, S2 60, S3 11) |

- **Category 1** (ranked): Pigr, Cyp4a14, Cyp2d9, Hpd, Hsd17b2. None peaks at S1/S2.
- **Category 2, top ten:** Dpys, Slc6a19, Ugt2b34, Cyp24a1, Slc22a7, Slc22a19, Higd1c, Slc5a12, Slc22a13, Scd1.

**Reading.**
- **Category 1 is real but small.** In each of the 5 genes, the cluster means put S2 and S3 at about the same level.
  The curve instead peaks at the transition and falls through S3. Mixing S2 and S3 tubules cannot produce a peak
  above both.
- **Category 2 is mostly within S2,** the longest segment. Slc6a19 (B0AT1) and Slc5a12 (SMCT2) fall steadily through
  S2, where the cluster means give one value. Neither pair was named in this protocol, so they are not used as
  recovery checks.
- **Many genes beyond the clusters are neighbour signal.** 37% of the genes where the curve beats the clusters are
  more than twice as high in thick limb, DCT or collecting duct. Within-PT structure is partly neighbour RNA; the
  filter keeps it out of both categories.

### Post-hoc check: within a single label

This check is not in the frozen protocol; it was added after the screen, before the figure was shown. Near a
boundary, S1 and S2 units overlap on x, so a pure step between labels can look graded. `fig4_label_check.py`
re-measures each call on units of one label only, from raw counts (`label_check.csv`).
- **Category 2.** 72 of 74 calls change in the same direction within the label in both mice: 68 by at least
  0.5 log2 and 43 by at least 1 log2. Measured within the label, the change is a median 0.83–0.88 of the curve's.
  Label overlap explains little of it.
- **Category 1.** All 5 rise within S2-labelled units toward the boundary (0.63–2.30 log2) and fall within
  S3-labelled units after it (0.39–2.40 log2), in both mice. Hsd17b2's fall in S3 is the weakest (0.39–0.45).

### Figure

`analysis/paper/fig4_panels.py`.
- **Panel a** shows Pigr, the top Category 1 gene.
- **Presentation choice made after the screen:** to avoid showing a gene twice, panel c shows Category 1 ranks 2–4
  (Cyp4a14, Cyp2d9, Hpd). Panel d shows Category 2 ranks 1–3 (Dpys, Slc6a19, Ugt2b34).
- **Panel b** has 79 rows: 5 Category 1 genes, then 74 Category 2 genes.
- **Y-axes** span the bin and cluster means. Where counts run out (Dpys and Ugt2b34 in late S3), the spline
  dives below the data and is clipped.

### Post-hoc rule for what the figure shows (2026-11-02, after the first figure)

**Problem.** Dpys, the top-ranked Category 2 gene, is an artifact. Its curve falls 2.4–5.0 log2 inside the S3 core
because the spline dives where counts run out. On raw counts from S3-labelled units it falls only 0.04 (Ctrl1A2)
and 1.42 (Ctrl1A4). The frozen ranking by the curve's change rewards such dives.

**Rule.** A call is shown only if it holds on raw counts within single labels, in the same direction and by at least
0.5 log2 in both mice (`fig4_label_check.py`).
- **Category 2.** Examples are ranked by that within-label change instead of the curve's.
- **Category 1.** The S2 rise and the S3 fall must each pass.

**Result.**
- **What survives.** 68 of 74 Category 2 calls and 4 of 5 Category 1 calls. Hsd17b2 drops out: its S3 fall is
  0.39–0.45.
- **Panel b** now has 72 rows.
- **Panel d** now shows Slc6a19, Scd1 and Cyp24a1, whose within-label changes are 2.1–3.0 log2. Panels a and c are
  unchanged.
- **The screen's counts above are unchanged.** This rule decides only what the figure shows.

### Panel b redrawn as curves against the cluster staircase (2026-11-02)

The heatmap showed the difference poorly. Its curve block and cluster block share the S1 → S3 trend, so they looked
alike.

**New panel b.** One plot per pattern of the calls that hold within labels:
- peak at S2/S3, 4 genes;
- rises within S2, 32;
- falls within S2, 23;
- falls within S3, 9.

The 4 remaining calls (S1, and one rise within S3) form groups too small to plot. They stay in the counts.

**What each plot shows.**
- **Units:** log2 fold from the gene's bulk PT mean, both mice averaged.
- **Lines:** each gene's curve thin, the median curve thick, and the median of the same genes' S1/S2/S3 cluster
  means as a staircase.
- **Shading:** the gap between the curve and the staircase, which is what the clusters miss.
- **Clipping:** each gene's curve is clipped to its own bin-mean range (±0.25), so spline dives where counts run
  out do not widen the gap.
- **Edge effect:** the upturn in the last ~5% of PT in two groups is a spline edge effect plus the last bin.

### Open issue: the uptick at the PT end (diagnosed 2026-11-08, not yet fixed)

The upturn at the end of PT (for example *Pigr* rising after falling through S3) is not biology. It has two causes.

**1. Mislabelled distal tubules.** The last 2% of control PT units on the gene axes are not late S3:
- median read depth is 1,000–1,800, against about 7,000 for PT;
- 49–89% are tubules the H&E segmentation called distal or collecting;
- they are small (about 380 µm², against about 1,550 for PT and 550–720 for distal tubules labelled distal);
- they sit among PT tubules (median 50–75% of the 8 nearest units are PT);
- they express TAL markers near mTAL levels (Umod about 7 log2 cp10k, Slc12a1 about 4);
- on the package curve they are placed in the TAL–DCT range (t up to 0.57–0.87, against the S3/mTAL cut at 0.28).

Reading: small distal cross-sections inside the cortex pick up PT transcripts from neighbouring bins. Clustering labels them PT, and the curve places them past S3. Rank-scaling PT to (0, 1) squeezes them all into the last 2%.

The same pattern holds on the hvm fit and in the AKI sections. Human seg_class uses other names, so the H&E check does not apply there.

**2. Spline tail.** The last knot is at 0.75, so the end of the curve overshoots the last bin. For *Pigr*: +1.2 at x = 1 against +0.5 for the last bin.

**Screen rerun without these units** (scratch only; 116 + 38 PT-labelled units whose H&E class is not tubule_proximal):

| Quantity | Original | Rerun |
|---|---|---|
| Beyond both mice | 1,680 | 1,599 (1,557 shared) |
| Category 1 | 5 | 6 (adds Cyp4a10) |
| Category 2 | 77 | 89 (71 shared) |

- **LR stability:** Spearman 0.98 and 0.99 by mouse; the median LR ratio is 0.96 and 0.99.
- **Category 2 genes lost (6):**
  - Alpl, Bcat1 and Dpep1 sit at the 1 log2 threshold;
  - Nat8f5, Npl and Osr2 changed sign, so the end units drove their call.
- **End hook** (|curve(1) − curve(0.95)| > 0.5 log2): 39% → 30% of beyond-both genes.
- **Pigr:** the LR is unchanged (938 → 939). Its end value drops from +1.2 to −0.5.

The remaining end units are still half-depth and small, but 49 of 54 are S3. The remaining hook is mostly the spline tail.

**Candidate fix, to freeze before rerunning:**
- drop PT-labelled units whose H&E class is not proximal; this rule uses neither depth nor position;
- clip the drawn curves to the bin-supported range.

### Post-hoc test: does the curve beat straight lines between the clusters? (frozen 2026-11-09, before the fit)

Question (user): would straight lines joining the S1/S2/S3 cluster levels show the same thing as the curve?

- **Null model (`linear`):** NB on raw counts with log(lib) offset. Design [1, x, (x − c₂)₊], with c₂ the S2 centre
  ((b12 + b23)/2 on the gene's own axis). This is the fitted piecewise-linear curve through the three segment centres,
  extended linearly at both ends: 3 parameters.
- **Alternative (`full`):** the null plus the screen's cubic B-spline design (BASIS_DF = 6). The null is nested in it.
- **Dispersion:** per gene, under `full`, as in the screen.
- **Statistics:**
  - NB likelihood ratio;
  - quasi-F, (LR/df₁)/(D_full/df₂).
- **Null calibration:** parametric bootstrap. Each gene is simulated once from its fitted `linear` NB model, with the
  same lib and dispersion. Dispersion is re-estimated and both models refitted. The pooled replicates set a 5% FDR
  threshold per mouse (decoy rule as in the screen). A gene passes if it clears the threshold in both mice.
- **Units:** two runs.
  - The screen's units, so the result matches the current figure.
  - The H&E-concordant units (open issue above).
- **Keeping it:** a gene is a Fig. 4b example only if it passes in both runs.

**Result (2026-11-09, scratch run; screen code unchanged).**
- **Calibration:** bootstrap replicate LRs average 5.0, matching the 5 added degrees of freedom, so the null behaves as
  designed. But the median real gene has LR 7.0, above that null mean, and the bootstrap trusts the NB model, so it is
  permissive.
- **Passing (quasi-F, 5% FDR, both mice):**

  | Units | Genes |
  |---|---|
  | Screen units | 906 |
  | H&E-concordant units | 828 |
  | Both runs | 809 |
  | Both runs, neighbour filter passed | 572 |

  76 of the 79 Fig. 4c calls pass.
- **Significance does not pick examples:**
  - Slc6a19 passes (rank 60) although straight lines reproduce its decline to within the bins.
  - For 319 of the 572, the largest departure lies at x 0.80–0.95, the S3 tail.
  - Adding effect size leaves 35 genes that pass and whose bins depart by more than 1 log2 in both mice. 15 of these
    sit in the interior (x 0.3–0.8).
- **Largest departures:** Cyp4a12a, Pigr, Dpys, Ceacam2, Scd1, Hs6st2, Slc22a13, Slc22a19.
  - Most share one shape: a sharp switch at x 0.40–0.46, partway into S2. Genes with it include Scd1, Slc22a13, Slc22a19,
    Ceacam2, Serpina1f, Higd1c, Mettl7a2, Hs6st2 and Mep1b.
- **Fig. 4b rule:** examples must pass the test in both runs and depart by more than 1 log2.

**Fig. 4b examples (2026-11-09, user: rank by the extra variance the curve explains).** `fig4_linear_gain.py` ranks every
screened gene by (D_lines − D_curve) / D_intercept, taking the smaller value of the two mice. The median gene gains 0.2%.

**Top 12:**

| Rank | Gene | Gain |
|---|---|---|
| 1 | Pigr | 21% |
| 2 | Ugt2b37 | 15% |
| 3 | Serpina1d | 15% |
| 4 | Serpina1f | 14% |
| 5 | Slc22a7 | 14% |
| 6 | Car5a | 14% |
| 7 | Slc22a13 | 13% |
| 8 | Cyp4a12a | 13% |
| 9 | Inmt | 12% |
| 10 | Higd1c | 12% |
| 11 | Apoa4 | 12% |
| 12 | Serpina1a | 12% |

**Panel b shows Pigr, Slc22a7 and Cyp4a12a** (ranks 1, 5 and 8). Each is drawn as one bulk bar, three cluster bars and the curve on top, all averaged over the two mice (user, 2026-11-09). The straight lines are used for ranking only.
- One gene for each shape the straight lines miss: a peak at S2/S3, a switch within S2, a dip just past S2/S3.
- Within each shape, the gene with the cleanest agreement between the two mice. Serpina1d/f and Slc22a13 dive in
  Ctrl1A4 where counts are near zero.
- Curves stop at the last bin.
- Chosen by hand from the top 10. `fig4_panels.py` asserts they stay in the top 10.
