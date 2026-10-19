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

## Notebook 104 results

`analysis/notebooks/104_pt_contamination_audit.ipynb` (logic `104.contamination.1`) writes to
`results/pt_contamination_audit/`. It followed the protocol frozen in 465e253. Processing the bins took 4 min for 8
specimen runs; the notebook ran from cache in 12 s.

**Guards.**
- G1: summing each export PT unit's bins reproduces its export counts exactly in all 8 specimen runs. This holds for
  HUK1_MED1 too, whose bins exist only as an MTX folder.
- G2: AKI 14,048 units; human vs mouse 11,383.
- G3: the pooled end contrast agrees with notebook 90's curve contrast at ρ 0.81 and 0.83 (mouse) and 0.78 and 0.73
  (human).
- G4: edge enrichment tracks non-native RNA. The median ε of the class markers is 1.0–1.6, and that of the 50 most
  expressed PT genes is −0.25 to −0.02.

**How the protocol was read.**
- The references pool the core bins of each species. In the AKI run they pool all four mice.
- Contact classes are PT, AL, DCT, CNT_CD, Glomerulus, outside and other.
- A scratch run first computed the Q2 depth reading on the fractions f instead of the marker shares s_k. The f values
  subtract the pooled PT-interior share, so they sit near 0 in the core by construction. This was corrected to the
  protocol's s_k before the notebook was executed.
- Post hoc: mean F is added beside the frozen median.
- The random aggregates are size-matched to the depth-≥-3 core. At depth ≥ 4 and ≥ 5, the ladder's sensitive shares
  include count loss.

**Markers** (own data, ≥ 20 × the PT interior and ≥ 2 × every other reference):
- AL: 23–30 genes, among them Umod, Slc12a1, Clcnka and Cldn16.
- DCT: 21–26 genes, among them Slc12a3, Klk1, Calb1 and S100g.
- CNT_CD: 55–71 genes, among them Aqp2, Hsd11b2 and Scnn1g.
- Glomerulus: 198–221 genes, among them Nphs2, Podxl, Ren1 and Emcn.
- Outside: few, because outside bins mix everything.
  - mouse: Tnc and Myh11;
  - human: immunoglobulins (Igkc, Ighg1, Ighm, Jchain), Myh11 and Pln.
- RC2 met: Umod, Slc12a3, Aqp2 and Nphs2 are markers of their class in both species.

| Reading | Result |
|---|---|
| Q1 non-native fraction F per unit | median −0.05 to 0.08 (relative to the pooled interior); mean −0.01 to 0.14; 2–24% of units have F ≥ 0.2 |
| Q2 marker share at the edge (layer 0) against the core (layers ≥ 3) | 1.9–4.5 ×; **not met in every specimen** (Ctrl1A2 1.87 and 1.98, short of 2) |
| Q2 contact with class k predicts f_k (ρ > 0 in every specimen) | AL 0.25–0.53, DCT 0.16–0.48, CNT_CD 0.05–0.33, outside 0.00–0.12: **yes**; Glomerulus: no (one NaN) |
| Q3 change from early to late PT (f, late − early) | AL +0.01 to +0.11; DCT −0.14 to +0.00; outside +0.27 in IR2A2 and +0.07 to +0.09 in human |
| Q3 AKI, IR − control | mean f_outside +0.072, mean F +0.013 |
| **Q4 gradients contamination-sensitive** | **mouse 28 / 1,342 (2.1%); human 22 / 627 (3.5%): native** |
| **Q5 AKI calls contamination-sensitive** | **level 47 / 3,159 (1.5%); location 1 / 362 (0.3%): native** |
| RC1 Egf, Klk1, Muc1 edge-enriched in both control mice | **missed**: ε 0.64 and 1.01, 0.53 and 0.97, 0.28 and 1.10 (edge-enriched in Ctrl1A4 only) |
| RC3 human ALB, TTR not edge-enriched (expected miss) | met: ε −0.98 to −0.26 |

**Q6 · Contamination exposure of the calling methods** (share of the calls' driver genes that are edge-enriched,
against expression-matched random genes):

| Contrast | Method | Driver genes | Exposure | Null | Ratio | p |
|---|---|---|---|---|---|---|
| AKI | positional F-test | 1,685 | 0.015 | 0.013 | 1.18 | 0.21 |
| AKI | clustering | 1,578 | 0.016 | 0.015 | 1.04 | 0.46 |
| AKI | position-split ORA | 1,041 | 0.014 | 0.016 | 0.88 | 0.74 |
| Human vs mouse | **positional F-test** | 1,038 | **0.074** | 0.052 | **1.44** | 0.001 |
| Human vs mouse | clustering | 527 | 0.091 | 0.054 | 1.69 | 0.001 |
| Human vs mouse | position-split ORA | 249 | 0.112 | 0.060 | 1.89 | 0.001 |

The contamination-sensitive genes look like what the test should catch.
- **Human gradients:** Calb1, Klk1 and Slc8a1 (DCT/CNT markers).
- **Mouse gradients:** Plat, Cxcl14 and Clu.
- **AKI level calls:** interstitial and immune genes (Kdr, C1s1, Serping1, Cd74, Ltbp4, Gdf10).

**What this answers.**
1. **Contamination is real, and it follows the anatomy.**
   - Non-native marker RNA is 2–4 times higher in the outermost bin of a unit than in its interior.
   - Each class's signal in a unit rises with how much of the unit's edge touches that class.
   - Along the PT, the signal tracks the neighbours: AL RNA rises toward S3, which lies beside the thick ascending
     limbs, while DCT and glomerular RNA are highest in early S1.
   - In the more injured IR mouse, interstitial RNA rises steeply in late S3.
   - Human S3 units carry immunoglobulin RNA.
2. **The positional findings and the AKI calls are native.**
   - Rebuilt from bins ≥ 6 µm inside the edge, 98% of mouse and 96.5% of human gradients survive, beyond the
     size-matched random bins. So do 98.5% of the AKI level calls and 99.7% of the location calls.
   - Every frozen "native" reading is met. Contamination makes a few dozen genes, not the gradients.
3. **A new yardstick:**
   - The human vs mouse calls of every method are enriched for edge-enriched genes.
   - The positional F-test is the least exposed (1.44×, against 1.69× for clustering and 1.89× for position-split
     ORA).
   - No AKI method is exposed.
   - The flagged genes (`method_exposure.csv`) should be read as neighbouring-structure RNA wherever they drive a call.
4. **The hand-flagged genes are partly confirmed.** Egf, Klk1 and Muc1 are edge-enriched in one control mouse and
   short of the cut in the other (RC1 missed). Klk1 is contamination-sensitive in human.
5. **For the paper:** a methods paragraph and a supplementary figure.
   - The bin-level core rebuild is a general check for any segmentation-aggregated Visium HD analysis.
   - It can say "our gradients are not neighbour RNA" with numbers.

## Literature context (added after the run)

These papers were read in full by the methods scout (`docs/research/methods-scout-contamination-networks-variability-evolution.md`)
after the protocol was frozen. None changed the analysis.
- **Visium HD diffusion is short-range** (Oliveira et al. 2025, Nat Genet, Fig. 2e, colon FFPE):
  - 1.0–1.7% of marker UMIs fall in the adjacent mask;
  - mislocalised UMIs move a median 3–6 µm.
  - Our core rule (≥ 6 µm inside the edge) sits just beyond that median. Our edge-to-core ratio of 2–4× matches a
    signal concentrated in the outer 2 µm.
- **Spillover can create false positional biology** (Mitchel et al. 2026, Nat Genet, preprint read). Goblet cells
  "up-regulate" stem and transit-amplifying markers in crypts and enterocyte markers on villi. This is the
  precedent for the risk tested here.
- **Swapped transcripts decay across boundaries like gradients** (Ni et al. 2022, Nat Commun, SpotClean; Visium v1).
  Decontaminating Visium counts with SoupX or DecontX creates false spatial patterns. This supports diagnosing by
  rebuilding units rather than correcting the counts.
- **Further checks the scout proposes** (not run):
  - a neighbour-expression covariate in the gene model (smiDE, Vasconcelos et al. 2026);
  - an overlap-ratio filter;
  - a planted-contamination arm at 1–4%.
