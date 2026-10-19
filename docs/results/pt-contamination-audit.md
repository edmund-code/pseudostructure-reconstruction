# Is the PT signal PT? A bin-level contamination audit (notebook 104)

**The question.** Each PT unit is the sum of the 2 µm Visium HD bins inside one segmented tubule cross-section.
Transcripts that diffuse across the polygon edge, and edge bins that straddle two structures, add the neighbours'
RNA to a unit. Neighbours differ along the PT: S3 runs through the outer stripe beside thick ascending limbs and
collecting ducts, while S1 sits among glomeruli and convoluted tubules. So contamination could make spurious
gradients along the axis. In AKI, interstitial infiltrate and oedema could also make spurious "PT" injury calls.

Earlier notebooks removed a few such genes by hand (Egf, Muc1, Klk1 as neighbouring-segment RNA;
`pt-pathway-literature-novelty.md`). This audit measures contamination for every unit and every gene from the bins
themselves. It then asks whether the paper's positional and AKI findings survive when every unit is rebuilt from its
interior bins only.

Frozen on 2026-10-18, before any bin-level statistic was computed. Structural checks were run before freezing, on
one specimen (IR2A2):
- rebuilding units from the bin-to-unit table reproduces the export counts exactly;
- 41% of PT bins and 40% of PT counts lie ≥ 3 bins (6 µm) from an edge.

## 1 · Labels and inputs

- Notebook `analysis/notebooks/104_pt_contamination_audit.ipynb`; logic `104.contamination.1`; seed 104.
- Results `results/pt_contamination_audit/`.
- Module `pseudospace/contamination.py`, tests `tests/test_contamination.py`.
- **Runs:** both package runs, `aki` (Ctrl1A2, Ctrl1A4, IR2A2, IR2A4) and `hvm` (Ctrl1A2, Ctrl1A4, HUK1_COR1,
  HUK1_MED1). For each run:
  - the bin folders are the inputs listed in its stage-01 manifest;
  - `01_units/<specimen>_bin_to_unit.parquet` maps bins to polygons;
  - the pass-1 labels (`03_integration/pass1_labelled.h5ad`, `coarse_class`) give each polygon's class;
  - the export `pt_units.h5ad` gives the PT units and positions (`support_positions`, rescale `support`; AKI
    14,048 units).
- **Bin classes:**
  - a bin in a polygon takes the polygon's pass-1 `coarse_class` (PT, AL, DCT, CNT_CD, Glomerulus, SmoothMuscle,
    Unassigned, Rejected_ATL);
  - a polygon absent from pass 1 is `unlabelled`;
  - an in-tissue bin in no polygon is `outside`.
- **Gene names:** mouse bins are matched to the export through the stage-01 unit objects' unique var names. Human
  bins are matched through the export's `human_symbol`.

**Guards.**
- G1: for every export PT unit, summing its bins reproduces its export counts exactly.
- G2: the AKI export has 14,048 units in the common support.

## 2 · Depth, core and random aggregates (fixed now)

- **Depth** of a bin: the Euclidean distance, in bins, to the nearest edge bin. An edge bin has a 4-neighbour with
  another polygon label, no polygon, or no tissue (`label_depth`).
- **Core aggregate** of a PT unit: the sum of its bins with depth ≥ 3 (≥ 6 µm inside the edge).
  - Sensitivity ladder: depth ≥ 1, 2, 4 and 5.
- **Random aggregates R1 and R2:** for each unit, as many bins as its core, drawn uniformly from all its bins
  (`size_matched_sample`, seeds 104 and 1104). They carry the core's loss of counts without its loss of edge.

## 3 · References, markers and the non-native fraction (fixed now)

The references are per species and pooled over that species' specimens in the run.
- **PT interior:** core bins of the export PT units.
- **Class k ∈ {AL, DCT, CNT_CD, Glomerulus}:** core bins of class-k polygons.
- **Outside:** all outside bins.

**Markers M_k** (`class_markers`, on the export genes):
- CP10k in reference k ≥ 1;
- reference k ≥ 20 × the PT interior;
- reference k ≥ 2 × every other reference.

**Fraction from class k:** f_k(u) = (s_k(u) − s_k(PT interior)) / (s_k(ref k) − s_k(PT interior)). Here s_k is the
share of counts on M_k. The non-native fraction is F(u) = Σ_k f_k(u).

**Edge enrichment** ε_g, per specimen: log2 CP10k of PT-unit bins with depth ≤ 1 minus that of PT-unit bins with
depth ≥ 3, pooled with 0.5 pseudocounts. A gene is **edge-enriched** when ε ≥ 1 in every specimen of the species.
- G4 checks that edge enrichment tracks non-native RNA: the median ε of the M_k markers is ≥ 1, and that of the 50
  most expressed PT-interior genes is ≤ 0.5. If G4 fails, Q6 is reported as uninterpretable.

## 4 · Questions and readings

**Q1 · Size (descriptive).** F per unit: its median per specimen, and the share of units with F ≥ 0.2. Each f_k is
reported separately.

**Q2 · Mechanism.**
- Depth profile: Σ_k s_k over PT bins in depth layers 0, 1, 2, 3, 4, 5 and ≥ 6, per specimen.
  - **Edge spillover holds** when layer 0 carries ≥ 2 × the share of layers ≥ 3 in every specimen.
- Contact: for each unit, the share of its edge bins touching class k (`contact_share`).
  - **Contact predicts contamination** for class k when Spearman ρ(f_k, contact_k) > 0 in every specimen.

**Q3 · Along the axis (descriptive).**
- The mean f_k in 10 equal-width position bins, per specimen.
- The end difference Δf_k = mean over x ≥ 0.7 minus mean over x ≤ 0.3.
- For AKI, IR − control in mean f_outside.

**Q4 · Do the gradients survive in the core?** (main reading, both species)
- **Graded genes:** notebook 90's `d_end_contrasts`, with |E| ≥ 0.5 and the same sign in both specimens of a species,
  as notebook 101 defined early and late.
- **Contrast per aggregate:** for each specimen and aggregate a ∈ {full, core, R1, R2}, Ê_a = `end_contrast` (pooled
  log2 CP10k of units with x ≥ 0.7 minus x ≤ 0.3, pseudocount 0.5).
- G3 (reported): Spearman ρ between Ê_full and notebook 90's E per specimen. Below 0.7 is recorded as an estimator
  difference.
- **Survives in a:** in both specimens, Ê_a has the sign of E and |Ê_a| ≥ 0.5 · |Ê_full|.
- **Contamination-sensitive:** survives in R1 and R2 but not in the core.
- **The gradients are native** in a species when ≤ 5% of its graded genes are contamination-sensitive.

**Q5 · Do the AKI calls survive in the core?**
- **Calls:** notebook 91's level calls (`level change`: 3,159 genes) and location calls (`location change`: 362).
- **Level:** L̂ for each of the 4 IR–control specimen pairs, as log2 CP10k over all common-support units.
  - Survives when all 4 pairs have the sign of 91's L and mean |L̂_a| ≥ 0.5 · mean |L̂_full|.
- **Location:** P̂ = mean over IR of Ê minus mean over controls of Ê.
  - Survives when its sign is that of 91's `Δ late` and |P̂_a| ≥ 0.5 · |P̂_full|.
- Contamination-sensitive as in Q4. **The AKI calls are native** when ≤ 5% of each call type is
  contamination-sensitive.

**Q6 · Contamination exposure of each calling method** (a new yardstick 8)
- **Calls:** notebook 93's call tables (positional F-test, clustering, row N; AKI and human vs mouse). Each method's
  exposure uses the union of its calls' driver genes.
- **Exposure:** the share of those genes that are edge-enriched. For human vs mouse, a gene counts in either species.
- **Null:** 1,000 draws of the same number of universe genes matched on expression decile (`matched_null`).
- Reported: exposure, null mean, ratio and empirical p.

**Recovery checks (named now).**
- RC1: Egf, Klk1 and Muc1 are edge-enriched in the mouse controls (flagged earlier by hand).
- RC2: Umod ∈ M_AL, Slc12a3 ∈ M_DCT, Aqp2 ∈ M_CNT_CD and Nphs2 ∈ M_Glomerulus, in mouse and human (if on the
  panel).
- RC3 (expected miss): human ALB and TTR, flagged earlier as ambient. Ambient RNA spread evenly over the section is
  not edge-enriched, so the edge test should not flag them.

## 5 · Outputs and figures

- **Tables:**
  - `units.csv` (run, unit, specimen, x, F, f_k, contact_k, counts full and core);
  - `markers.csv`, `depth_profile.csv`, `edge_enrichment.csv.gz`;
  - `gradient_robustness.csv.gz`, `aki_robustness.csv.gz`, `method_exposure.csv`;
  - `recovery_checks.csv`, `how_the_protocol_was_read.csv`, `run_record.json`.
- **Figures (plots only):**
  - 104a: non-native share against depth, per specimen;
  - 104b: f_k along the axis, per specimen;
  - 104c: Ê full against core for graded genes, each species;
  - 104d: method exposure against its matched null.

## 6 · Caveats

- Depth is measured on the bin grid (to within one bin), not on the polygon.
- Marker-based fractions assume the markers are absent from PT cells. The 20 × rule makes the PT interior's share
  small, but injured PT could induce a marker.
- Core bins sit away from the basal edge and the lumen edge. A gene whose mRNA localises to the basal pole is
  "edge-dependent" for a biological reason. A contamination-sensitive gene is therefore edge-dependent, not proven
  contamination.
- Ambient RNA spread evenly is invisible to the core test.
