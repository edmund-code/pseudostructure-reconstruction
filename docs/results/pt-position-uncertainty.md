# Uncertainty of reconstructed PT positions (human vs mouse)

Status: protocol frozen 2026-10-21, before any draw was fitted. Results are appended below the protocol.

## Why

Positions on the reconstructed pseudospace are estimates. Human positions are inherited from the mouse through
integration (a human-only reconstruction orders the PT mainly by mitochondrial fraction), and a scratch
sensitivity check across 15 integration recipes moved human positions within a segment about twice as much
as mouse positions (per-unit SD 0.050-0.064 against 0.026-0.036 PT lengths). Downstream claims that depend on
position (boundaries, peak locations, onsets, within-segment trends) need an interval that carries this
uncertainty. The 15 recipes were chosen ad hoc and are not a sample from anything, so they cannot give one.

## Method: spatially blocked delete-half jackknife of the frozen pipeline

- **Recipe.** The frozen human-vs-mouse configuration (`analysis/package_configs/human_vs_mouse`), unchanged:
  stage 03 (two-pass Harmony, Leiden, automatic labels) and stage 04 (elastic curve), starting from the
  stage-02 QC units of the reference run.
- **Resampling unit.** Spatial blocks, not tubules: neighbouring tubules share position (spatial
  autocorrelation of position 0.63-0.87), so tubule-level resampling would understate the variance. Each section
  is cut into a square grid whose cell edge gives a mean of 25 units per occupied cell (edge set per section
  from its own unit density, so pixel scale does not matter).
- **Draw.** In every section independently, keep a uniformly random half of the occupied blocks (floor of
  n_blocks / 2, without replacement), and every unit in them; rerun stages 03 and 04 on the kept units.
- **Draws.** B = 100, seeds 0..99 (`numpy.random.default_rng(seed)`).
- **Threads.** Every draw runs with 4 BLAS/OpenMP threads (`--threads 4`). Added after the two-draw pilot,
  before the full run: the same half-sample gave a different Leiden partition at 8 threads than at 4
  (fingerprints 5f719beb3d1b vs e1e9051a1edd; some positions moved by up to 0.8), while two runs at 4 threads were
  identical. A draw is reproducible only at a fixed thread count. The 8-thread pilot draws were discarded.
- **Scale.** Each draw's positions are put on the PT scale used throughout (PT units, common support across
  specimens, min-max to [0, 1], oriented S1 to S3 by label medians).
- **Variance.** Delete-d jackknife with d = n/2 (Shao & Wu 1989, Ann. Stat. 17:1176): the variance estimate of a
  full-data statistic is (r / (d B)) sum_b (theta_b - theta_bar)^2 with r = n - d; at d = n/2 the factor r/d is 1,
  so it is the mean squared deviation of the half-sample statistics. Intervals are full-data estimate +/- 1.96 SE.
  Half-sampling is used rather than the ordinary bootstrap because duplicated units would create artificial
  neighbours in the kNN graph that Leiden and the curve fit consume.
- **Per-unit positions.** Each unit is present in about half the draws; its SE is computed from the draws that
  contain it.
- **Calls.** A yes/no call (direction flip, boundary order) is reported with the share of draws in which it holds.

Every draw that produces an axis is kept. Draws whose stage-04 order metrics fail (any family out of nephron
order) are counted and reported, and results are also shown without them as a sensitivity check; they are not
silently dropped.

## What it does not cover

- Biological replication. The human arm is one donor; draws resample tubules inside sections and say nothing about
  another person.
- Recipe choice. A pre-declared multiverse (normalisation {total, mitochondria-free} x species offset {none,
  family-weighted} x Harmony theta {2, 6, 10}; 12 fully crossed recipes) is reported separately as a
  sensitivity analysis and is not folded into the intervals.

## Pre-declared first uses

1. Human-minus-mouse S1->S2 and S2->S3 boundaries (balanced logistic), with SE.
2. S3-zone vs S1-zone log2 fold change of Mme, Dcxr and Ctsh per species, and the share of draws in which each
   reverses direction between species.
3. Per-unit position SE by species and zone.

Pass rule for a positional claim: the full-data estimate's 95% interval excludes the null value, and the call
holds in at least 95% of draws.

## Results (2026-10-21, 100 draws, 4 threads each, none failed to produce an axis)

Each draw kept 48.5-51.8% of the QC units (1,084 blocks in total, mean about 25 units per block). 76 draws had a
segment order with no inversions. Of the 24 with inversions, 19 kept the family order PT < AL < DCT < CNT_CD and
the PT label order (the inversions sat elsewhere in the fine-label order). Five draws, all among the 24, are
outright failures of the PT axis: in seed 35 the PT is reversed (mouse Mme comes out S1-high), and in seeds 22, 54
and 94 the human PT is folded (Ctsh flat in human, human S1->S2 boundary about 0.2 too early). The frozen pipeline
therefore breaks on about 5% of half-samples, and those draws dominate the jackknife variance. That rate is for
half the data and need not apply at full n.

Pre-declared statistics (full-data estimate, delete-half jackknife SE, 95% interval, share of draws where the call
holds):

| statistic | all 100 draws (primary) | 76 draws without inversions (declared sensitivity) |
|---|---|---|
| S1->S2 boundary, human - mouse | 0.051 +/- 0.059 (-0.065, 0.166); later in 94% | 0.051 +/- 0.022 (0.007, 0.095); later in 96% |
| S2->S3 boundary, human - mouse | 0.021 +/- 0.072 (-0.121, 0.163); later in 65% | 0.021 +/- 0.025 (-0.029, 0.071); later in 64% |
| Mme S3/S1 log2FC, mouse / human | 3.03 (2.20, 3.85) / -1.34 (-1.80, -0.88) | 3.03 (2.87, 3.18) / -1.34 (-1.47, -1.21) |
| Dcxr S3/S1 log2FC, mouse / human | -0.95 (-1.23, -0.68) / 1.70 (1.11, 2.29) | -0.95 (-1.06, -0.84) / 1.70 (1.56, 1.84) |
| Ctsh S3/S1 log2FC, mouse / human | 0.70 (0.50, 0.91) / -0.90 (-1.25, -0.55) | 0.70 (0.58, 0.83) / -0.90 (-0.99, -0.82) |
| Mme, Dcxr, Ctsh reverse direction between species | 100% of draws each | 100% of draws each |

Against the pass rule:
- **The three direction reversals pass** in the primary analysis: every interval excludes zero, and every draw
  reverses.
- **The human S1->S2 boundary delay does not pass in the primary analysis** (the interval includes zero, and it
  holds in 94% of draws). It passes only in the sensitivity analysis, and only narrowly (lower bound 0.007, 96%).
  It is not reported as established.
- **The S2->S3 boundaries do not differ.**

Median per-unit position SE (PT lengths), all draws: human 0.080 / 0.102 / 0.097, mouse 0.064 / 0.071 / 0.088 in
the S1 / S2 / S3 zones. Not pre-declared: after aligning each draw to the reference by a least-squares affine map
(which removes the draw-wide stretch and shift of the PT scale), human 0.033 / 0.051 / 0.071 and mouse
0.041 / 0.022 / 0.034. Most of a unit's raw spread is the whole scale moving. That component cancels in
within-draw contrasts such as human minus mouse. Human S2 and S3 positions stay roughly 2x less certain than
mouse ones after alignment.

Outputs (private, under the results root): `package_pseudospace/position_half_sampling/` holds the per-draw
tables, `pt_positions_by_draw.parquet`, `pt_position_summary.parquet`, `statistics_by_draw.csv`,
`predeclared_statistics.csv` and `draws.csv`. The pre-declared recipe multiverse is not yet run.
