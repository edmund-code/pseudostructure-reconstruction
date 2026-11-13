# Fig. 2–3 validation experiments: protocol and results

Paper handoff section 6, items 2, 4, 5, 6, 7 and 8, run on the data available locally on 2026-11-02.
The protocol was frozen before any of these fits ran. Results are appended below it and never edited
to fit it.

## Inputs and what they cost

- **Fits.** The package runs under `results/package_pseudospace/runs/{aki,hvm}`, stages 01–04, from the
  v4 segmentations. The AKI fit has 32,939 units on the curve. The paper's reference run
  (v4fix, 32,785) is not readable here. Every number below moves slightly when it becomes available.
- **Labels.** These are the stage-03 segment labels: Leiden clusters of the *same* expression profiles,
  named by marker panels. The curve fit never sees them, except that stage 04 picks its root from
  coarse cluster 0. They are therefore a consistency check, not independent anatomical truth. Figure
  text must not call them "held-out".
- **Order.** Labels map to seven groups in nephron order, PTS1 < PTS2 < PTS3 < AL < DCT < CNT < CD
  (`analysis/paper/fig3_common.py`; mTAL and cTAL are both AL, CCD and IMCD are both CD). The score is
  the unit-level Spearman ρ between a coordinate and the group rank, within each section.
- **Sections.** Ctrl1A2, Ctrl1A4, IR2A2 and IR2A4 are scored in the AKI fit. HUK1_COR1 and HUK1_MED1
  are scored in the human-vs-mouse fit.

## Frozen protocol

**3b. Labels along the axis (item 5).**
- Draw a ridgeline of t per group, one panel per fit.
- Report ρ per section.
- Null: permute t within the section 1,000 times; p = (1 + #null ρ ≥ observed) / 1,001.

**3d. Robustness (item 6).** Each perturbed dataset reruns the frozen stages 03 and 04, with the same
configuration, the R harmony backend and 4 BLAS threads, as `analysis/scripts/position_half_sampling.py`
does.
- **Perturbations.** All three use fractions 10, 20, 30, 50, 70 and 100%, with three seeds each:
  - *tubules*: a random subset of stage-02 QC units, stratified by section;
  - *reads*: binomial thinning of every unit's counts;
  - *genes*: a random subset of the genes offered to variable-gene selection.
- **Score.** ρ between the rerun's t and the reference t on the units both contain, plus the rerun's ρ
  against the *reference* labels.
- **Checks.** The 100% arms are the reproduction check. Gene halves are reported as five
  complementary splits; the ρ between the two halves of each split is reported too.
- **Expected:** ρ ≥ 0.95 down to 30% tubules. Whatever comes out is reported.

**3e. Trajectory baselines (item 7).** Every method runs on each fit's pass-2 Harmony embedding (50
dimensions), and every method gets the same root, the stage-04 root unit.
- **DPT** is the stage-04 DPT baseline, already computed.
- **scFates** uses `tl.tree`, method `ppt`, on the same embedding.
- **Slingshot** uses the pass-2 Leiden clusters (unsupervised), with the root unit's cluster as the start.
  When there are several lineages, its averaged pseudotime is used.
- **GASTON** runs per section on unit centroids plus expression PCs. It is the only spatial method; its
  isodepth is oriented so the root's section end is low.
- **Score:** ρ per section against the label rank.

**3f. Cortical depth (item 8).**
- **Tissue mask.** Threshold the saturation of each H&E slide at about 16 µm/px with Otsu, fill holes,
  and keep the largest component.
- **Capsule.** The capsule is the mask boundary within 100 µm of the mask's convex hull. Boundary
  further inside the hull is the hilum or a cut, and is not capsule.
- **Depth.** Depth is the distance from each unit's centroid to the nearest capsule pixel, in µm.
- **Score.** Spearman ρ between depth and t over PT units per section, plus depth by label.
- **Human sections.** If a human section shows no capsule, it is reported as not computable rather
  than given a substitute edge.
- Depth selects nothing.
- **Expected.** Positive but modest ρ. S1 and S2 both lie throughout the cortical labyrinth, and
  juxtamedullary nephrons start deep.

**Fig. 2c gallery (item 4, partial).**
- **IR2A4 tile.** A tile in the IR2A4 outer stripe, chosen as the 250 µm window with the most
  PT-S3 units in the AKI fit, with the v4 predicted outlines.
- **Human tile.** A tile from HUK1_MED1, the human section with no training blocks, with the v2
  predicted outlines.
- PQ stays unknown until someone draws manual outlines.

**Fig. 2b ablations (item 2, mouse only).**
- **Arms.** Retrain the full model, CNN stem only and ViT only, on the v4 manifest with the v4
  configuration and identical seeds.
- **Selection.** Select each arm by its pooled class-agnostic validation PQ.
- **Score.** Report PQ with D4 TTA on the v4 validation blocks.
- **Reference arm.** The full arm stands in for the locked v4 checkpoint and must land near its
  0.727.
- **Not run.** Human ablations are not run, because the human labels are not readable here.

## Results

(appended as each experiment finishes)

### 3b. Labels along the axis (2026-11-02)

`analysis/paper/fig3_panels.py`, panel b. Unit-level Spearman ρ of t against the seven-group rank:

| Section | Fit | Units | ρ | Largest permuted ρ | p |
|---|---|---|---|---|---|
| Ctrl1A2 | AKI | 7,110 | 0.934 | 0.040 | 0.001 |
| Ctrl1A4 | AKI | 6,703 | 0.946 | 0.044 | 0.001 |
| IR2A2 | AKI | 9,150 | 0.933 | 0.036 | 0.001 |
| IR2A4 | AKI | 9,976 | 0.972 | 0.036 | 0.001 |
| HUK1_COR1 | human vs mouse | 5,489 | 0.949 | 0.051 | 0.001 |
| HUK1_MED1 | human vs mouse | 6,020 | 0.945 | 0.046 | 0.001 |

- **The permutation null is weak by construction.** Shuffling t against the labels gives ρ ≈ 0, so the
  stronger comparison is the baselines in 3e. For reference, the stage-04 DPT baseline scores
  0.88–0.94 on the same sections.
- **The AKI fit has no CNT cluster.** Its labels are DCT2, CCD and IMCD. The human-vs-mouse fit has CNT.
- **The human PTS2 and PTS3 ridges sit earlier than the mouse ones.** The two fits each have their own
  0–1 scale, so the ridges are not on a shared axis.

### 3c. Marker heatmaps (2026-11-02)

- The panel is re-plotted from data, replacing the poster crops, which said "Pseudospace" and came from
  the v4fix run.
- It uses one marker per segment, the same genes in both species: Slc5a2, Slc13a3, Slc7a13, Slc12a1,
  Slc12a3, Slc8a1 and Aqp2. This settles handoff open decision 5 on the human-poster set, because every
  gene is measured in both fits.
- Each row is the mean lognorm expression in 100 equal-width bins of t, z-scored per gene. Bins with
  fewer than 5 units are blank, and the label bar shows the majority label per bin.

### 3f. Cortical depth (2026-11-02)

`analysis/paper/fig3_depth.py`. The capsule and tissue masks were checked by eye on all four mouse
sections; the overlays are in `results/paper_figures/fig3/depth/qc_*.png`. The hilum and cracks are
excluded from the capsule.

**Human sections: not computable.** On both human sections, every outer edge is a cut face. Tubules are
cut open at the edge and there is no fibrous capsule (checked at full resolution).

| Section | PT units | ρ(t, depth) | ρ(label rank, depth) | ρ within PTS1 | ρ within PTS2 | ρ within PTS3 | Median depth S1 / S2 / S3 (µm) |
|---|---|---|---|---|---|---|---|
| Ctrl1A2 | 2,780 | 0.75 | 0.74 | −0.09 | 0.44 | 0.56 | 262 / 459 / 985 |
| Ctrl1A4 | 2,745 | 0.71 | 0.70 | −0.28 | 0.40 | 0.61 | 235 / 407 / 1009 |
| IR2A2 | 3,937 | 0.46 | 0.49 | −0.41 | 0.41 | 0.23 | 429 / 500 / 1297 |
| IR2A4 | 5,551 | 0.27 | 0.29 | −0.46 | 0.29 | 0.23 | 541 / 601 / 1249 |

**Reading.**
- In the controls, t ranks PT cross-sections by depth as well as the segment labels do.
- Within PTS2 and PTS3, t still rises with depth (ρ 0.4–0.6), which is resolution the labels do not carry.
- Within PTS1, depth falls slightly as t rises. This is not interpreted: S1 begins at glomeruli spread
  through the whole cortex.
- In the AKI sections, the correlation is weaker and PTS1-labelled units sit deeper. Whether that comes
  from injury or from the plane of section is not separable here.
- **Not excluded:** a depth-dependent technical gradient, such as capture efficiency, would also produce
  a positive ρ. Depth stays a consistency check, not truth.

### Fig. 2c gallery (2026-11-02)

`analysis/paper/fig2_gallery.py`. Each window was chosen from unit centroids before the image was viewed.

- **IR2A4** has 21 PT-S3 units in the window, at pixel (18033, 4683), with 49 predicted outlines; casts
  are visible.
- **HUK1_MED1** has 18 PT units in the window, at pixel (11912, 4904), with 48 predicted outlines.

The outlines register on the H&E. PQ still needs manual outlines of both tiles.

### 3d. Robustness (2026-11-02)

`analysis/paper/fig3_robustness.py`, 64 reruns per fit (three arms × six fractions × three seeds, plus five gene
splits × two halves). None failed, and none reversed the global orientation.

**Deviations from the protocol.**
- **Thread count.** The 100% reruns at 4 threads do not reproduce the reference exactly: they are deterministic
  among themselves but land on a different Leiden partition. Rerunning the unperturbed data with the thread
  variables unset, as the reference was run, reproduces it exactly (same fingerprints, max |Δt| = 0). So the
  gap comes from the thread count, not the code. The grid keeps 4 threads as frozen; `rho_vs_full_rerun`
  scores each rerun against the 4-thread 100% rerun and barely differs from ρ vs reference.
- **Genes arm.** Only the pool for the integration HVGs is restricted (through `exclude_from_harmony_hvg`). The
  2,000 clustering HVGs that name the labels are not, because the package has no hook for them.
- **Added check.** `pt_limb_reversed` flags reruns whose PT limb is entered from its S3 end (median t of PTS1 above
  PTS3). This partial flip keeps ρ positive, so global orientation alone misses it.

**Median ρ with the reference t** (three seeds; lowest single rerun in brackets where it differs by > 0.05):

| Kept | Mouse fit: tubules | reads | genes | Human + mouse fit: tubules | reads | genes |
|---|---|---|---|---|---|---|
| 10% | 0.983 | 0.458 (0.318) | 0.807 | 0.981 | 0.466 (0.358) | 0.790 |
| 20% | 0.985 (0.664) | 0.772 | 0.891 (0.755) | 0.985 | 0.541 (0.433) | 0.854 (0.773) |
| 30% | 0.988 | 0.788 | 0.947 | 0.982 | 0.755 (0.655) | 0.904 (0.838) |
| 50% | 0.992 | 0.955 | 0.961 (0.549) | 0.987 | 0.930 (0.852) | 0.895 |
| 70% | 0.991 | 0.980 | 0.983 | 0.988 | 0.957 | 0.964 |
| 100% | 0.995 | 0.995 | 0.995 | 0.992 | 0.992 | 0.992 |

- **Expected result met.** ρ ≥ 0.95 down to 30% tubules holds in every 10–70% tubule rerun but one: mouse
  `tubules_f020_s2` (0.664), whose PT limb is reversed. Down to 10% tubules, ρ stays at 0.98.
- **Reads are the weak arm.** Below 50% of reads, ρ falls to 0.46–0.79. Thinned clusters fall below stage 03's
  absolute panel-score threshold and are dropped as Unassigned, so only 13.6k of 32.9k mouse units reach the
  curve at 10% reads. Most limb reversals are in this arm (mouse 6 of 8, human 8 of 9).
- **Against the reference labels,** pooled median ρ stays at 0.95 for every tubule fraction in both fits; the reference
  itself scores 0.952 (mouse) and 0.957 (human + mouse).
- **Gene halves, ρ between half A and half B:** mouse 0.942, 0.908, 0.950, 0.940, 0.952; human + mouse 0.915,
  0.851, 0.874, 0.825, 0.901. Two disjoint gene sets give nearly the same order.

**Reading.** The axis depends little on which tubules are sampled and somewhat more on which genes are
offered. Sequencing depth matters most: half the reads still gives a median ρ ≥ 0.93, but below that, units drop out of
labelling and the PT limb can attach from the wrong end.

### 3e. Trajectory baselines (2026-11-02)

`analysis/paper/fig3_baselines.py` (Slingshot through `fig3_slingshot.R`). Settings, versions and runtimes are in
`results/paper_figures/fig3/baselines/methods_notes.md`. Slingshot and GASTON ran in separate environments,
passed by `--rscript` and `--gaston-python`; `kidney-pseudospace` was not changed.

| Method | Ctrl1A2 | Ctrl1A4 | IR2A2 | IR2A4 | HUK1_COR1 | HUK1_MED1 |
|---|---|---|---|---|---|---|
| pseudostructure t | 0.934 | 0.946 | 0.933 | 0.972 | 0.949 | 0.945 |
| DPT | 0.883 | 0.898 | 0.885 | 0.945 | 0.948 | 0.938 |
| scFates (ppt) | 0.927 | 0.942 | 0.934 | 0.964 | 0.950 | 0.947 |
| Slingshot | 0.837 | 0.863 | 0.824 | 0.874 | 0.908 | 0.862 |
| GASTON | 0.552 | 0.490 | 0.432 | 0.328 | 0.039 | 0.086 |

**Settings.**
- **scFates:** σ was chosen by `explore_sigma` (σ = 1 in both fits), with 50 nodes, λ = 100 and 50 steps.
  It had not converged at 50 steps; 200 and 500 steps gave the same ρ to three decimals. The trees branch (7 tips
  in the mouse fit, 4 in the human-vs-mouse fit), and their pseudotime is used as given.
- **Slingshot:** 5 lineages in the mouse fit and 6 in the human-vs-mouse fit; the score uses
  `slingAvgPseudotime`.
- **GASTON:** the settings of its Visium HD tutorial. It ran on the two mouse sections of the human-vs-mouse
  fit as well, but those are not scored.

**Reading.**
- **On this score, t does not beat scFates.** The two are within 0.01 on every section. t beats DPT by up to
  0.05 in mouse, and ties it in human; it beats Slingshot by 0.04–0.11.
- **The score is close to its ceiling.** The labels are clusters of the same profiles, so any method that orders
  the clusters correctly scores about 0.93–0.97. It cannot separate methods that differ *within* a segment. The
  case for t over scFates has to come from elsewhere: it is one unbranched curve with a fixed root rule, not a
  tree whose pseudotime depends on the branch a unit falls in.
- **GASTON measures something else.** In mouse, its isodepth follows cortex-to-medulla depth: PTS1 and DCT are
  shallow and CD deep, so it partly agrees with segment order. Both human pieces are cortex only, and there it
  finds no axis (ρ < 0.09). Its human orientation rests on a root cluster of 111 and 72 units.
- **Depth is not used to rank the baselines.** Doing so would make depth select a method, which the protocol
  rules out.

### Fig. 2b ablations (2026-11-02)

The architecture switch `model.backbone` (`full` | `cnn_only` | `vit_only`) is in `segmentation/`; `full` is the
default and matches the previous model exactly. `full` is U-Net + OpenMidnight: the CNN stem is the U-Net's encoder
half, and OpenMidnight's tokens join it at stride 16. `cnn_only` is a plain U-Net; `vit_only` drops the U-Net encoder
and its skips, so the decoder upsamples OpenMidnight's tokens alone. Each arm was trained on the v4 manifest with the v4 configuration,
seed 1337, for 54 epochs as v4 was, and selected by pooled validation PQ. Patches were re-extracted from the v4
manifest; all 6,475 match v4 slide by slide and class by class. Software differs from v4: torch 2.10 instead of 2.5.1.

| Arm | Trainable | Epoch | Val PQ, no TTA (1,298) | D4 Nx_4wk_050 (98) | D4 Ctrl_1A2_val (356) | **D4 pooled (454)** |
|---|---|---|---|---|---|---|
| v4 (locked) | 7.31 M | 35 | 0.727 | 0.880 | not run | – |
| U-Net + OpenMidnight (`full`) | 7.31 M | 29 | 0.735 | 0.870 | 0.731 | **0.762** |
| U-Net only (`cnn_only`) | 5.15 M | 39 | 0.729 | 0.867 | 0.767 | **0.789** |
| OpenMidnight only (`vit_only`) | 4.30 M | 27 | 0.726 | 0.852 | 0.735 | **0.760** |

**Reference arm.** The full arm reproduces v4: validation PQ 0.735 vs 0.727, and D4 on Nx_4wk_050 0.870 vs 0.880
(one object of 98).

**Reading.**
- **OpenMidnight does not raise PQ here.** On the frozen score, D4 pooled over the validation blocks, U-Net only
  leads U-Net + OpenMidnight by 0.027 and OpenMidnight only trails it by 0.001.
- **The U-Net-only lead comes from one block.** On the held-out slide U-Net only is 0.004 below the full model. Its lead is on
  Ctrl_1A2_val, mostly the 50-object block at (9709, 10240): 0.71 vs 0.61.
- **The arms are within noise.** There is one seed per arm, and per-region PQ spans 0.61–0.93 (full), 0.71–0.91
  (U-Net only) and 0.59–0.89 (OpenMidnight only). A claim that either branch matters would need several seeds per arm.
- **Scoring window.** That block's 2,048 px window runs 714 px past the annotated validation block. The scorer drops
  predictions on unannotated tissue before scoring, so the overshoot adds no false positives.

**Fig. 2b now shows** the three mouse control bars at D4 pooled. The human bar is still v2's training-time
validation PQ (0.824, no TTA). It is a different score from the mouse bars, and its D4 number needs the human
labels and checkpoint, which are not readable here. Mouse AKI and the human ablations stay empty.

**What the empty bars need.**
- **Human, all three bars.** The HUK1_COR_1 training labels and the human v2 manifest and configuration, which are in
  an account this one cannot read. Then three arms train as the mouse arms did, about 3.5 h each, and score with
  D4 TTA. With the v2 checkpoint as well, the human full bar can be v2's own D4 score instead of a retrained stand-in.
- **Mouse AKI, all three bars.** Dense manual outlines on two or three IR2A2/IR2A4 blocks (cortex and outer stripe,
  about 150 tubules), in the QuPath classes of `Ctrl_1A2_final.geojson`. Each block needs an ROI rectangle that
  marks it fully annotated. The outlines must be drawn from the H&E, not by correcting model outlines, because
  correcting model outlines biases PQ towards the model. No retraining is needed: the three trained mouse arms
  score the new blocks with D4 TTA. AKI sections were never in training, so every block is held out.
- **Fig. 2d from the same work.** If a second person outlines the same AKI or human blocks, the annotator-vs-annotator
  error for Fig. 2d comes from those blocks too.

## COAST fits on the v4fix segmentations (2026-11-12)

The user replaced the mouse segmentations with v4fix and made COAST (github.com/edmundtsou/COAST, commit db5502c) the
method. Both fits were rebuilt with `analysis/paper/coast_fit.py` (results/coast_runs/{aki,hvm}); `fig3_common.RUNS`
points there. The current Fig. 3 has panels a and r only (`make_fig3.py`); the robustness, baseline and depth panels
(3d–f above) are retired and were not rerun. Only the ten gene-half reruns per fit were redone, with COAST, because the
Fig. 4–6 screens place genes on them.

| | control/AKI | human/mouse |
|---|---|---|
| polygons / with bins / after QC | 64,080 / 43,761 / 36,610 | 51,176 / 31,524 / 27,614 |
| units on the curve | 31,928 | 25,352 |
| PT-S1 / S2 / S3 median t | 0.058 / 0.119 / 0.196 | 0.063 / 0.107 / 0.186 |
| Spearman ρ with the package fit (all / PT) | 0.990 / 0.978 | 0.981 / 0.92–0.99 per section |
| label agreement with the package fit | 84% | 92% |

COAST's QC keeps exactly the units the package's stage 02 kept (36,610 and 27,614 vs 27,612). The labels differ: on
the control/AKI fit COAST calls 1,373 package PT-S2 units PT-S1 and has no CNT cluster (package CNT units are DCT2 and
package DCT2 units DCT1), because Leiden runs on COAST's own Harmony neighbour graph. Labels were assigned
automatically and still need the per-cluster review.

**Fig. 3r, ρ per section (p = 0.001 against 1,000 within-section shuffles in every section):** Ctrl1A2 0.950, Ctrl1A4
0.947, IR2A2 0.939, IR2A4 0.971, HUK1_COR1 0.945, HUK1_MED1 0.944 (package on v4: 0.934, 0.946, 0.933, 0.972, 0.949,
0.945).

**Gene halves (COAST):** all ten control/AKI halves keep PT order; one (s1_A) runs the distal limb backwards
(segment-order ρ 0.67) while its PT positions agree with the reference at ρ 0.96. The screens use PT positions only.
