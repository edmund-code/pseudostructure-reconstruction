# Planted truth on the human vs mouse data (notebook 103)

**The question.** Notebook 100 measured the calling methods' power and false calls by planting known pathway
differences into our AKI data. Do the same profiles hold under the human vs mouse noise structure: one human donor
with two sections, two mice, and species-specific probe panels?

The protocol is notebook 100's (`pt-planted-benchmark.md`), with its post hoc low-dose arm and its clean-false-call
reading adopted as frozen here. Only the base data change. Frozen on 2026-10-18, before any planted human vs mouse data
were made.

## Changes from notebook 100

- Notebook `analysis/notebooks/103_pt_planted_benchmark_hvm.ipynb`; logic `103.planted_hvm.1`; seed 103.
- Results `results/pt_planted_benchmark_hvm/`.
- **Base data:** the human vs mouse export through `hvm_inputs` (mt excluded), notebook 71's universe and sets.
- **Base partition:** notebook 71's relabeling `swap:HUK1_COR1+Ctrl1A4`.
  - Group 1 is HUK1_COR1 + Ctrl1A4, with species as the nuisance shape.
  - Each side holds one human section and one mouse.
- **Fair null for clustering:** the other relabeling, `swap:HUK1_COR1+Ctrl1A2`, as notebook 71 ran it.
- **Zone cuts:** notebook 71's human vs mouse cuts.
- **Planting:** identical to notebook 100. Shapes are centred over the group-1 units, and the logistic sits at the
  control-derived S2/S3 cut of this dataset's mouse units (`segment_cuts`), not at 0.65.

  *Reading:* notebook 100 used 0.65 for the AKI data, its controls' cut. The equivalent here is the mouse units' cut.
- **Arms:**
  - frozen: 15 sets per replicate;
  - low dose: 3 sets per replicate (notebook 100's post hoc arm, frozen here).

  Each arm runs 4 scenarios × 2 amplitudes × 8 replicates.
- **Readings:**
  - notebook 100's frozen rule (new false calls ≤ 2);
  - notebook 100's post hoc clean-false rule (calls on sets sharing no planted gene ≤ 2), frozen here as the primary
    reading, because notebook 100 showed that all-new-false counts are dominated by spillover.
- **Comparison with notebook 100:** for each method, scenario and amplitude, the sensitivity in human vs mouse
  against AKI. A method's profile **generalises** when its detect or miss reading is the same in both datasets for
  every scenario at A = 1.

## Notebook 103 results

`analysis/notebooks/103_pt_planted_benchmark_hvm.ipynb` (logic `103.planted_hvm.1`) writes to
`results/pt_planted_benchmark_hvm/`. It is notebook 100's code on the human vs mouse base, and ran in 32 min for 128
replicates (20 s from cache).

**Guards and baseline.**
- G1: the universe (11,055 genes) and the sets (1,513) are notebook 71's.
- G2: the cached model reproduces 71's z and its **0** F-test calls on the base partition.
- G3: clustering reproduces its 0 calls.
- The planted logistic sits at the mouse units' S2/S3 cut, 0.605.
- The segment baseline calls 0 sets.

**Sensitivity, human vs mouse against AKI (A = 1 / 0.5):**

| Scenario | Positional F-test | Segment F-test | Clustering |
|---|---|---|---|
| shape, 15 sets | 1.00 / 0.92 (AKI 1.00 / 0.99) | 1.00 / 0.93 | 0.88 / 0.13 (AKI 0.92 / 0.21) |
| located, 15 sets | 1.00 / 0.97 | 1.00 / 0.97 | 0.79 / 0.18 |
| level, 15 sets | 0.00 / 0.00 | 0.00 / 0.00 | 0.93 / 0.69 (AKI 1.00 / 0.92) |
| bump, 15 sets | 1.00 / 0.92 | 1.00 / 0.88 | 0.84 / 0.05 |
| shape, 3 sets | 1.00 / 0.88 | 1.00 / 0.92 | 0.67 / 0.00 |
| located, 3 sets | 0.96 / 0.75 | 1.00 / 0.96 | 0.46 / 0.12 |
| level, 3 sets | 0.00 / 0.00 | 0.00 / 0.00 | 0.96 / 0.88 |
| bump, 3 sets | 1.00 / 0.79 | 1.00 / 0.42 | 0.71 / 0.00 |

**Clean false calls** (sets sharing no planted gene):
- human vs mouse: positional F-test 0.6–1.9 per replicate, segment F-test 0.5–1.4, clustering 0;
- AKI (notebook 100): 5–12 for the F-tests.

**Generalisation** (the same detect reading in both datasets for every scenario at A = 1, with clean false ≤ 2):
- frozen arm (15 sets): all three methods **generalise**;
- low-dose arm (3 sets): none does, but for one reason. On AKI the F-tests' clean false calls (9–10) exceed 2, while on
  human vs mouse they do not (0.9–1.9). Clustering misses "located" at low dose on human vs mouse (0.46 against 0.83).

**What this answers.**
1. **The methods' profiles are properties of the methods, not of one dataset.**
   - The positional F-test finds planted shape, S3-only and S2-bump changes at both doses on both datasets, and never
     finds pure level.
   - Clustering finds level on both, and positional change mainly at 1 log2.
   - Every method is somewhat weaker on human vs mouse: one donor, and the specimen-level noise is larger.
2. **The F-tests' clean false calls in notebook 100 came from the AKI base, not from the test.**
   - That partition already had 29 F-test calls, so a planted program relaxed BH enough to let near-threshold sets
     through.
   - The human vs mouse base had 0 baseline calls, and the excess stays below 2.
   - In practice, a contrast with many marginal calls gains a few more when strong signal is present. Program-level
     reading (notebook 102) handles this.
3. **Spillover is general.** On human vs mouse, 3 planted sets at 1 log2 make 34–129 other sets significant for the
   F-tests and 54–141 for clustering.
