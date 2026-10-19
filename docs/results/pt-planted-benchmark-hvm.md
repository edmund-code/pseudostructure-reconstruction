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
