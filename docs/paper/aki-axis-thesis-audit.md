# The AKI and mouse axis: an audit of "continuous beats discrete"

2026-10-21. The species axis round (`hvm-continuous-discovery-round.md`) found that the paper's
second title claim is not supported on the human-versus-mouse axis, because the location statistic
it rested on reads amplitude as position. That verdict was explicitly scoped: the mouse and AKI
evidence for continuity is a different measurement on a different axis, and it was put under its own
audit rather than assumed to survive. This is that audit.

Six published claims were re-derived from the committed tables, each attacked from three directions:
amplitude posing as the effect (the lens that destroyed the species round), coordinate-construction
genes being scored on the axis they built, and whether the blind null can see the claimed effect at
all. A second track rebuilt the species convergence ladder on package-era lines only.

**The headline: the thesis survives in a much narrower form than published.** Two claims fail, one
inverts, three hold with their statistics restated. One result is genuinely continuous-only,
measured, and reproducible — transition width — and it was a secondary bullet in the paper.

## 0 · What I re-derived myself before writing any of this down

Every number in the verdict table below was recomputed from the committed tables in this session, by
me, not taken from an agent's report. Six checks, all exact:

| # | Check | Result | Status |
|---|---|---|---|
| 1 | 28% / 3% from `a_gene_classes.csv.gz`, non-flat and fold-protected | n = 5,259; 1,480 = 0.281422; 142 = 0.027001 | reproduces |
| 2 | The same comparison made symmetric: `min(bic D, bic DS)` vs `min(bic L, bic S, bic B)`, same ±6 margin | 910 = 0.173037 continuous; 348 = 0.066172 steps; ratio **2.615**, not 10.42 | published ratio does not survive |
| 3 | Positional-information ceiling | max set `I` = 1.3669 against log2 3 = 1.5850; **0 of 1,356** sets reach the 3-level ceiling; 69 of 1,356 beat the labels on median error | published statistic is resolution-bound |
| 4 | Change position, from the independent `gene_edges_agreed.csv.gz` (2,966 two-mouse-agreed edges, equal-count bins) | enriched **at** the boundaries: 1.108 / 1.191 / 1.244 / 1.242 at windows 0.10 / 0.15 / 0.20 / 0.30; depleted far from them (0.4295 against 0.5408 uniform) | published claim **inverts** |
| 5 | Planted-truth head-to-head, positional F-test minus segment F-test, paired on the identical planted data | 32 cells; mean margin **exactly 0.0 in 19**; \|mean\| > its own replicate sd in 2; bump +0.125 (15 wins, 0 losses), located −0.031, shape +0.002 | published reading **fails** except for the bump |
| 6 | The interior-dip rule of `pt-positional-decomposition.md` applied to all 16 arms of `i4_program_profiles.csv` | fires in **14 of 16**, including both controls, whose relative excess (0.744, 0.831) **exceeds** the injured mice's (0.566, 0.501) | published claim **fails** |

Also reproduced exactly on the way: `i_survival.csv` 0.522306 / 0.657144 / 0.311265 with relabelling
0.9230–1.2260 (B2 stands); the 34 main rising pathway edges, median width 1.276264, which becomes the
reported 33 edges and median 1.238930 once deduplicated by pathway; and 32 of 33 of those edges wider
than the coordinate's own amplitude-matched response to a true step, median ratio **2.593**.

Nothing in this audit used anything under `data/external/**`.

## 1 · The verdicts

### M3 bullet 1 — gene shape classes, 28% against 3% · HOLDS, NARROWER

The published figures are exactly right and the published *ratio* is not. `dBIC continuous` is
`bic D` minus the best of three continuous models — it lets the continuous family field its best
player and pins the step family to one. The step family has a second member in the same table,
`DS` (steps plus a within-segment slope). Make the comparison symmetric and 28.1% / 2.7% becomes
**17.3% / 6.6%**, a ratio of 2.6 rather than 10.4; 254 of the 1,480 "continuous" genes have `DS` —
a step-family model — as their actual best fit of the six.

What the three attacks found, and this is the part worth keeping:

- **Amplitude does not pose as the effect here.** Spearman(`dBIC continuous`, `range best`) =
  **0.179**, against 0.891/0.956 on the species axis. The discrete-favoured genes have the *larger*
  median amplitude (1.1691) than the continuous ones (0.9151), so "big genes look continuous" is
  false in this table.
- **But decidability is a power readout.** Across deciles of amplitude the continuous-favoured share
  climbs 0.148 → 0.544. So "28%" is substantially the share of genes whose amplitude lets BIC decide
  anything at all.
- **The asymmetry is what is amplitude-free.** Continuous-among-decided is 0.95, 0.93, 0.91, 0.93,
  0.96, 0.93, 0.90, 0.90, 0.85, 0.92 — flat across the entire amplitude range.

**The honest sentence:** of 5,259 non-flat genes that did not build the coordinate, the comparison
decides for 30.8%, and where it decides it chooses a continuous shape **91.2%** of the time, at every
amplitude decile. Published as "28% against 3%" it invites the symmetric recount that halves it.

### M3 bullet 2 — 26 pathways localise a unit better than its segment label · HOLDS, NARROWER

26 and 132 reproduce exactly. The problem is that the claim's verb and its statistic are different
columns. `above labels` is a comparison of **mutual information**, and mutual information about a
25-bin axis is bounded by the predictor's own cardinality: a three-level predictor cannot exceed
log2 3 = 1.585 bits. The labels score 1.123, leaving 0.462 bits of headroom to that ceiling — and
every one of the 26 margins (0.009 to 0.244) fits inside it. No set of 1,356 reaches the ceiling.
So "more bits than the label" is substantially a 25-versus-3 output-resolution statement.

The localisation statistic is the error column, and it is the one to publish: **69 of 1,356** sets
predict a held-out unit's position with lower median error than that unit's own S1/S2/S3 label, and
the 26 span 0.120–0.147 against the label's 0.177.

### M3 bullet 3 — change density concentrates beside, not at, the transitions · FAILS

Every committed number reproduces (excess peaking at 0.8 and 2.1–2.2; 124 and 277 at the
transitions against 1,110 and 1,081 for a step transcriptome). But an independent committed table
disagrees, and it is the better instrument: `pt_transition_sharpness/tables/gene_edges_agreed.csv.gz`
holds 2,966 per-gene edges that agreed between the two mice, fitted on **equal-count bins**, so a
densely populated bin cannot manufacture a peak. On that table change is **enriched at** the label
boundaries at every window I tried — 1.108, 1.191, 1.244, 1.242 at ±0.10, ±0.15, ±0.20, ±0.30 — and
genes further than 0.30 from both boundaries are *depleted* (0.4295 observed against 0.5408 uniform).

The two tables differ in exactly the way bin density predicts, and the equal-count one wins.

**Withdraw "beside, not at".** What survives is the claim that actually bears on the discrete model:
a transcriptome constant within S1/S2/S3 does not describe where the mouse PT changes — 28.7% of the
reproducible excess change sits more than 0.3 units from both boundaries against 1.2% for a step
transcriptome built from the same per-bin label fractions. The audit reports a constant-within-label
control with real count noise at 4.6%, making the honest ratio 6–20× rather than 23×; **I did not
re-derive that control**, so it is carried as the agent's number pending its own check.

### M3 bullet 4 — within-segment positions (84) and transition sharpness (85) · HOLDS, NARROWER — and this is the one real continuous-only result

Split the bullet, because its two halves have opposite outcomes.

**Transition width (85) is askable only on a continuous axis and is answered there.** 33 pathway
edges, median width **1.2389** axis units (IQR 1.03–1.70), and **32 of 33 are wider than this
coordinate's own response to a true step of the same amplitude** — median 2.59×. A step world cannot
produce that, and the comparator is the coordinate's own measured step response, not a model.
Per-edge rank agreement between mice is 0.484 and falls to 0.242 once amplitude is partialled out,
so publish the set-level median and the 32-of-33 count, not a per-pathway width.

**Within-segment position (84) is not answered on this axis.** Its segment window is built from the
mean of all four specimens' transitions rather than the mouse-only mean, so its "S2" spans 48.5% of
the mouse's own S2 and each reported "third of S2" is 0.162 axis units against a measured label-step
width of 0.26–0.33. Its own table reads "no within-segment separation" in both testable mouse
segments. Drop the within-segment-position half of the bullet.

### M4 — the planted benchmark as "the continuous screen detects what the discrete screen misses" · FAILS

The sensitivities are all correctly reported; the comparison they are read as making was never run
in the paper's favour. Paired on the identical planted datasets, the positional F-test minus a
three-segment step comparator gives a mean margin of **exactly 0.0 in 19 of 32 cells**, and exceeds
its own across-replicate spread in 2. At 1 log2 both are at ceiling.

The exception is real and worth keeping precisely because it is specific: the **S2 bump** — the one
planted shape a segment model cannot express — gives **15 positional wins and 0 segment losses**,
mean margin +0.125. That is the benchmark's actual result: the continuous screen's advantage is
confined to interior structure, not to positional structure in general.

What the benchmark does separate cleanly is the two screens the project actually runs: the
positional F-test finds shape and is blind to level; clustering finds level. That is the defensible
reason to run two screens, and it needs no continuous-beats-discrete framing.

### B1–B3 — the AKI geography · B1 and B2 HOLD, "measured continuously" comes out, B3 FAILS

**B1 holds as biology and fails as a methods claim.** The S3 weighting replicates in both mice
(ACUTE-up S3−S1 +0.2106 and +0.1946), 8 of 8 injury and loss programmes are monotone S1 < S2 < S3,
and critically **it passes the amplitude check that destroyed the species axis** (Spearman −0.03 to
+0.11). But a plain segment-label pseudobulk reproduces the continuous drift at **Spearman 1.000000**
(Pearson 0.9939). So "measured continuously" must come out of B1's title: the continuous axis is how
we measured it, not what made it visible.

**B2 holds.** 0.5223 overall and 0.3113 from S2 to S3 against 0.9230–1.2260 under relabelling, all
exact, and recoverable from segment means alone — same caveat as B1, same verdict.

**B3 fails.** The project's own interior-dip rule, applied to all 16 arms rather than the two it was
written for, fires in **14 of 16** — including both control mice, whose relative excess (0.744 and
0.831) is *larger* than the injured mice's (0.566 and 0.501). A rule that fires harder in the
controls cannot establish that late S1 is the least responsive stretch of the injured PT. Notebook
114's separate z-score reading of "S1 is spared" is not touched by this and is where that sentence
should live instead.

## 2 · What the thesis looks like after the audit

Across both axes, the claim "a continuous axis finds biology a segment model misses" now has exactly
three measured supports, and none is the one the paper leads with:

1. **Transition width.** 32 of 33 pathway edges are wider than the coordinate's own step response,
   median 2.59×. Only askable continuously; answered.
2. **Within-segment constancy fails.** 28.7% of reproducible change sits more than 0.3 units from
   both label boundaries, against 1.2% for a step transcriptome. This says the discrete model is
   wrong *inside* segments — which is a different and better claim than saying change avoids the
   boundaries.
3. **Interior structure in the benchmark.** The S2 bump, 15 wins to 0.

Everything else is a **resolution and screening** advantage: the continuous fit supplies the per-gene
quality machinery, the amplitude and separation diagnostics, the within-group nulls, and a 25-bin
readout where the label gives 3. That is a real contribution and it is what the methods section
should claim. It is not a detection advantage, and on the species axis it is not even that
(`hvm-continuous-discovery-round.md`).

## 3 · Where the two axes agree, and where the species round's killers did not transfer

The species round died of one mechanism: a statistic that reads contrast strength as position, acting
on an arm that is systematically flatter. Both halves were tested here and **neither transfers**:

- Amplitude dependence of the shape statistic is 0.179 here against 0.891/0.956 there, and the
  discrete-favoured genes are the *larger* ones.
- The AKI S3-weighting passes the amplitude surrogate at −0.03 to +0.11.

What *does* transfer is the lesson about asymmetric comparisons. The species round's COM compared a
flat arm with a steep one; M3 bullet 1 compares a one-model family with a three-model family;
M3 bullet 2 compares a 25-bin predictor with a 3-level one. In all three the asymmetry, not the
biology, set the size of the effect. **Any statistic that compares the continuous axis with the
segment label must give the label its best available model and state the cardinality of both
predictors.** That is now a standing requirement, beside the amplitude surrogate and the
within-group null.

## 4 · Open questions this leaves

1. The constant-within-label control at 4.6% is the number the surviving bullet 3 claim rests on and
   I have not re-derived it. Do that before it is published.
2. Bullet 1's symmetric recount (17.3% / 6.6%) should be the published pair, with the
   91.2%-of-decidable asymmetry beside it. Decide which is the headline.
3. `fold-protected` has **opposite polarity** between `a_gene_classes.csv.gz` (True for 93.8%) and
   `d_end_contrasts.csv.gz` (unprotected ~5.7%). Nothing in this audit depends on it, but a reader
   comparing the two tables will be misled. Rename one.
4. The stale 141/115 in `pt-continuous-axis-experiments.md` should read **145/118**; the committed
   table and `readout.csv` both say 145.
5. `pt-continuous-axis-experiments.md` names Slc34a1 and Slc5a2 as genes with "continuous shapes
   favoured", but the lowest-BIC model of six for both is a **step** model (`discrete + within`).
   Fix or drop the examples.
