# Planted truth: power and false calls of the pathway calling methods (notebook 100)

**The question.** Every comparison of the calling methods so far lacked a ground truth. Relabeling bounds false calls,
and external data (notebooks 94, 98) checks generalisation, but neither says how often a method finds a real
difference of a known kind.

This notebook plants known pathway-level differences into our own data, with its real specimens, units, depth and
noise. Each method runs exactly as in notebook 71:
- the positional F-test;
- the segment F-test (its discrete twin);
- clustering.

It then counts how many planted sets each finds, and how many unplanted sets each starts calling.

Protocol frozen on 2026-10-18, before any planted data were made.

## 1 · Labels and inputs

- Notebook `analysis/notebooks/100_pt_planted_benchmark.ipynb`; logic `100.planted.1`; seed 100.
- Results `results/pt_planted_benchmark/`.
- Module `pseudospace/planting.py` (set choice, effect shapes, count planting), tests `tests/test_planting.py`.
- **Base data:** the AKI export through `aki_inputs`, built as notebooks 71 and 93 build it: the universe of 10,198
  genes, 71's sets, the common-support units and the library.
- **Base partition:** notebook 71's relabeling `swap_IR2A2_Ctrl1A4`.
  - Group 1 is IR2A2 + Ctrl1A4, and the true condition is the nuisance shape.
  - Each side holds one IR and one control mouse, so injury cancels.
  - Notebook 95 found this relabeling uncorrelated with IRI (ρ −0.04 at 12 h).
  - It is the relabeling in which the planted effects are the only group difference by design.
- **Cached baseline (read, not refit):**
  - 71's gene model for this partition;
  - its T_spatial and LR_spec;
  - its clustering run, whose fair null is the other relabeling.
- The segment F-test's baseline for this partition is fitted once here with `segment_step_lr`, under 71's dispersion.

## 2 · Planting (fixed now)

- **Sets planted per replicate:** 15 sets drawn from 71's AKI sets with 20–150 members. Each draw is accepted only if
  its Jaccard with every set already chosen is below 0.1. A gene is planted in one set at most.
- **Members planted:** a random 60% of each set's members, at least 10.
- **Direction:** one random sign per set, shared by its planted members.
- **Effect on a unit u of a group-1 specimen,** for a planted gene: multiplier m(u) = exp(ln 2 · A · sign · f(x_u)).
  Group-0 specimens are untouched.
- **Shapes f** (x = position in [0, 1]):

  | Scenario | f(x) | Kind |
  |---|---|---|
  | shape | logistic((x − 0.65) / 0.06) minus its mean over the group-1 units | pure positional |
  | located | logistic((x − 0.65) / 0.06) | an S3-only change (level plus shape) |
  | level | 1 | pure level |
  | bump | exp(−(x − 0.5)² / (2 · 0.1²)) minus its mean over the group-1 units | inside S2 |

  The logistic sits at the controls' S2/S3 cut (0.647).
- **Amplitude** A ∈ {0.5, 1} log2. 8 scenarios × 8 replicates = 64 planted datasets.
- **Counts** (planted genes, group-1 units only):
  - m < 1: binomial thinning, c' ~ Binomial(c, m);
  - m > 1: c' = c + Poisson(c · (m − 1)).
  - *Ponytail:* the added counts are Poisson, so an increase adds a little less overdispersion than the NB would.
    Regenerating planted genes from the fitted NB, as notebook 89 does, is the upgrade.
- **Library:** unchanged. With at most about 600 genes planted, the change in total counts per unit is reported, not
  corrected.

## 3 · The methods (notebook 71's code, refitting only the planted genes)

Planted genes are refitted and spliced into the cached baseline. Every other gene keeps its baseline statistics.

| Method | Refit for planted genes | Then |
|---|---|---|
| Positional F-test | `standard_gene_model` (group = the partition, nuisance = the condition, dispersion fixed to the cached α) and `specimen_shape_lr` | `moderated_f` over all genes (covariate log mean count, robust FALSE, df1 6, df2 6), then `camera_pr` on z. Call: Up and FDR ≤ 0.05. |
| Segment F-test | `segment_step_lr` (same group, nuisance, α) | `moderated_f` over all genes on T_step, then `camera_pr`. Call: Up and FDR ≤ 0.05. |
| Clustering | the gene model's specimen curves for planted genes, spliced into the cached model | `clustering_first_inputs`, `compare_groups` (fair null = the other relabeling), `run_pathways`. Call: matched q < 0.05 (notebook 71's `run_b`). |

## 4 · Metrics and readings

Per method, scenario and replicate:
- **Sensitivity:** planted sets called / 15.
- **New false calls:** unplanted sets called now but not at baseline.
- **Spillover:** the new false calls that share a planted gene.
- **Lost calls:** baseline calls no longer called.
- **Precision:** planted calls / (planted calls + new false calls).

Averages over the 8 replicates, with their range.

**Readings (fixed now).**
1. **A method detects a scenario** when its mean sensitivity at A = 1 is ≥ 0.5, with mean new false calls ≤ 2 per
   replicate.
2. **Which kinds each method sees:** the detection profile over the four scenarios at A = 1 and at A = 0.5.
3. **Best method per scenario:** the highest mean sensitivity among methods that detect it. A tie within 0.1 is
   reported as a tie.

The expectations from the methods' definitions are not hypotheses about biology:
- the positional F-test tests shape only, so it should miss pure level;
- clustering includes level groups;
- the segment test sees steps but may miss a bump inside S2.

The benchmark measures how far each holds on our noise.

## 5 · Outputs and figures

- **Tables:**
  - `replicates.csv` (planted sets, signs, genes);
  - `baseline_calls.csv`;
  - `scores.csv` (method × scenario × A × replicate);
  - `summary.csv`, `readings.csv`, `library_change.csv`;
  - `how_the_protocol_was_read.csv`, `run_record.json`.
- **Figures (plots only):**
  - 100a: the four planted shapes;
  - 100b: sensitivity per method and scenario at A = 0.5 and 1;
  - 100c: new false calls per method and scenario.

## 6 · Caveats

- The truth is planted at the pathway level, with 60% of members moving together. Real programs may be weaker or less
  coherent.
- Splicing keeps every unplanted gene's baseline statistic. That is exact for the F-tests' per-gene step and
  approximate for the moderation prior, which is recomputed over all genes.
- One base partition. The other relabeling carries the A2/A4 injury-like factor, so it is not a null base.

## Notebook 100 results

`analysis/notebooks/100_pt_planted_benchmark.ipynb` (logic `100.planted.1`) writes to `results/pt_planted_benchmark/`.
It followed the protocol frozen in a8f55de, and the module is `pseudospace/planting.py` with `tests/test_planting.py`
(3 tests). Runtime is 36 min for 128 replicates.

**Guards.**
- G1: the universe (10,198 genes) and the sets (1,469) are notebook 71's.
- G2: the cached model reproduces 71's z and its 29 F-test calls on the base partition.
- G3: rerunning 71's clustering reproduces its 0 calls.
- The knots match.
- The segment baseline (df1 2, df2 2) calls 17 sets.

**How the protocol was read.** Planted genes are refitted with the cached dispersion, grid and knots. Clustering's fair
null is the other relabeling. The seeds are (100, scenario, 10·A, replicate). Two labelled post hoc additions follow.

**Post hoc 1: diagnostic of the false calls.** One replicate (shape, A = 0.5) was rebuilt step by step.
- The splice is exact: refitting the original counts reproduces the cached T_spatial at r 0.9997.
- Unplanted genes' z are unchanged (r 0.9999997), and planted genes reach a median z of 2.3.
- The extra calls (29 → 694) are set-overlap spillover. Removing the planted genes from the other sets leaves 10
  calls.
- BH amplifies the spillover when hundreds of sets shift.

**Post hoc 2: the low-dose arm (3 sets) and clean false calls** (new calls on sets sharing no planted gene).

**Sensitivity** (share of planted sets called; A = 1 / 0.5):

| Scenario | Positional F-test | Segment F-test | Clustering |
|---|---|---|---|
| shape, 15 sets | 1.00 / 0.99 | 1.00 / 0.97 | 0.92 / 0.21 |
| located, 15 sets | 1.00 / 1.00 | 1.00 / 1.00 | 0.90 / 0.45 |
| level, 15 sets | 0.06 / 0.02 | 0.03 / 0.02 | **1.00 / 0.92** |
| bump, 15 sets | 1.00 / 1.00 | 1.00 / 0.96 | 0.96 / 0.23 |
| shape, 3 sets | 1.00 / 0.96 | 1.00 / 0.92 | 1.00 / 0.00 |
| located, 3 sets | 1.00 / 0.96 | 1.00 / 0.88 | 0.83 / 0.33 |
| level, 3 sets | 0.00 / 0.04 | 0.00 / 0.04 | **1.00 / 0.92** |
| bump, 3 sets | 1.00 / 0.96 | 1.00 / 0.92 | 1.00 / 0.00 |

**False calls per replicate** (all new calls; clean = sharing no planted gene):

| Arm, A | Positional F-test all (clean) | Segment F-test all (clean) | Clustering all (clean) |
|---|---|---|---|
| 15 sets, 1 | 620–650 (0–0) | 586–610 (0.4–1.4) | 471–689 (0.1–1.0) |
| 15 sets, 0.5 | 545–582 (4.8–5.8) | 458–515 (4.6–6.8) | 51–503 (0.3–3.1) |
| 3 sets, 1 | 164–212 (8.8–9.9) | 138–182 (9.0–10.6) | 124–342 (0.1–4.6) |
| 3 sets, 0.5 | 87–139 (5.9–11.9) | 73–116 (5.0–11.3) | 0–129 (0–3.4) |

The ranges are over the shape, located and bump scenarios (the F-tests) and over all four (clustering). In the level
scenario the F-tests add 0–4 calls.

**Readings.**
- **Frozen rule** (sensitivity ≥ 0.5 with ≤ 2 new false calls): no method detects any scenario, because spillover
  dominates every new call. This is reported, not tuned.
- **Post hoc rule** (≤ 2 *clean* false calls):
  - the F-tests detect shape, located and bump at 15 sets and A = 1;
  - clustering detects all four at 15 sets and shape, located and bump at 3 sets.
  - The F-tests' 5–12 clean false calls at A = 0.5 and at low dose come from BH relaxing when many sets shift: sets
    just below threshold at baseline cross it.

**What this answers.**
1. **Each method sees what it is built to see.**
   - The positional F-test is the most sensitive test of positional change: shape, S3-only and S2 bump changes are
     found at 0.5 log2 at either dose. It is blind to pure level shifts.
   - Clustering finds level shifts at 0.5 log2, and shape changes only at 1 log2.
   - The two complement each other, which confirms the earlier verdict ("better positional test, not better overall";
     "clustering as a level layer") with a known truth.
2. **The segment F-test matches the positional F-test on these planted sets,** S2 bump included. For set-level
   detection, a continuous model buys little sensitivity over steps here. Its gains are in where and how a change
   sits (notebooks 85, 89, 99, 101), not in whether a set is called.
3. **Every set-level method spills.** A planted program makes 100–650 other sets significant through shared genes.
   The number of called sets is not the number of findings. Notebook 102 measures which grouping restores programs.
4. **The F-tests' set calls carry a small BH-driven excess** when real signal is present: 5–12 clean false sets per
   replicate against clustering's 0–5. Their calls should be read as programs (notebook 102), and their marginal calls
   with care.
