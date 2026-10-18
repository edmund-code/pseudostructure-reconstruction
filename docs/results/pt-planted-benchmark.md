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
