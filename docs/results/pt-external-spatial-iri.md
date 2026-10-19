# Where along the PT does injury act? Replication in independent spatial IRI data (notebook 105)

**The question.** Our positional AKI claims rest on two IR mice:
- acute injury rises along the PT toward S3 (notebook 96 H1: ρ 0.71 in both IR mice; notebook 99 H1);
- transport loss is S3-weighted (notebook 94);
- the S2→S3 gradient survives injury less than the S1→S2 gradient (notebook 91 I2: pooled survival 0.31 against
  0.66).

External replication so far could not test position. Kirita 2020 is single-nucleus (segment labels only), and human
AKI biopsies carry no positional signal (notebook 98). Two public spatial IRI datasets from labs independent of ours
and of Kirita have biological replicates at our time points:
- **Cheung et al. 2022** (JCI Insight, GSE200115): Visium v1, whole-kidney sections, male C57BL/6J, bilateral IRI.
  3 mice each at quiescent, 12 h and day 1.
- **Gharaie et al. 2023** (Sci Rep, GSE233512): Visium v1, male C57BL/6J, bilateral IRI. 2 control and 2 IRI
  capture areas at 24 h on GEO.

Both papers were read in full by the dataset scout (`docs/research/external-spatial-aki-datasets.md`). Visium v1
spots are 55 µm mixtures of several cells, so a spot's PT position is the position of a PT-dominated mixture.

**A known hazard.** Wu et al. 2024 (Nat Commun, Results) note that injury lowers the PT segment-defining markers. An
axis built from markers could therefore misplace injured S3 spots toward S2. This would bias against the
S3-weighting hypothesis, not toward it. The axis below uses only genes whose gradient survived injury in our own mice.

Frozen on 2026-10-18, before any expression statistic of either dataset was computed. Structural checks only:
- every Cheung barcode maps onto the Visium v1 grid;
- median UMI per spot is 11–27k for the quiescent and 12 h sections, and 1.1–15k on day 1, where two sections are
  shallow.

## 1 · Labels and inputs

- Notebook `analysis/notebooks/105_pt_external_spatial_iri.ipynb`; logic `105.external_spatial.1`; seed 105.
- Results `results/pt_external_spatial_iri/`.
- Module `pseudospace/external_spatial.py`, tests `tests/test_external_spatial.py`.
- **Data:**
  - `data/external/cheung2022_iri_visium/` (9 `.h5`);
  - `data/external/gharaie2023_iri_visium/` (4 MTX sets with positions).
  - Cheung spot positions come from the Visium v1 barcode grid (Gharaie's `tissue_positions_list`, all 4,992
    barcodes).
  - Grid to µm: x = column × 50 µm, y = row × 86.6 µm.
- **Reference for deconvolution:** Kirita 2020 nuclei (control, 12 h and 2 d files), summed per broad type:
  - PTS1, PTS2, PTS3;
  - NewPT (NewPT1 + NewPT2: injured PT);
  - DTL-ATL;
  - TAL (MTAL, CTAL1, CTAL2, MD);
  - DCT (DCT, DCT-CNT);
  - CNT;
  - PC (PC1, PC2);
  - IC (ICA, ICB);
  - Pod, PEC;
  - EC (EC1, EC2);
  - Stroma (Fib, Per);
  - Immune (Mø, Tcell);
  - Uro.
- **Our programs** (notebook 96's): INJ (notebook 76's 18 injury markers), TPT (93 PT transporters) and CYC (297
  G2-M and E2F genes).

## 2 · Spots, PT purity and position (fixed now)

- **Deconvolution:**
  - Markers per broad type: the 50 genes with the highest ratio of that type's CP10k to the maximum over the other
    types, among genes with CP10k ≥ 1 in the type.
  - Each spot's CP10k on the union of markers is fitted by non-negative least squares to the types' CP10k profiles.
  - Weights are scaled to sum to 1.
- **PT-dominant spot:** PTS1 + PTS2 + PTS3 + NewPT weight ≥ 0.6. Sensitivity: ≥ 0.75.
- **Positional genes** (from our data, fixed now):
  - notebook 91's genes with survival `kept` (the gradient survives injury in our IR mice);
  - graded in both our control mice (notebook 90's E ≤ −0.5 for early, ≥ +0.5 for late);
  - PT-enriched in Kirita, with mean CP10k over PTS1/PTS2/PTS3 ≥ 2 × the maximum non-PT broad type (NewPT
    excluded).
- **Position x:**
  - The positional score is the mean z of log1p CP10k over the late genes minus that over the early genes. Each z is
    standardised on the PT-dominant spots of the dataset's control sections.
  - x is the score's within-section rank among PT-dominant spots, scaled to (0, 1).
- **Sensitivity position:** the deconvolution segment position (0 · w_S1 + 0.5 · w_S2 + 1 · w_S3) / (w_S1 + w_S2 +
  w_S3), ranked within section.
- **Cross-check (descriptive):** anatomical depth, the distance in µm from each spot to the nearest grid position
  without tissue. Reported as ρ(x, depth) per section. The hilum also counts as an edge, so medullary depth is
  underestimated.

## 3 · Programs, residuals and readings

- **Expression:** log1p CP10k per spot.
- **Scores:** `program_score` standardised on the control sections' PT-dominant spots (SD floor 0.25).
- **Residual e:** the score minus the control expectation at x (`position_expectation`, width 0.1, cross-fitted for
  control sections).
- **INJ-PT** (a guard against injury genes from neighbouring thick limbs in outer-stripe spots): the INJ markers whose
  Kirita CP10k in NewPT (12 h and 2 d) exceeds their CP10k in TAL at the same time points.

**R1 · Injury rises along the PT (H1).** ρ(e_INJ, x) per section. Replicates in a dataset and time point when:
- every IR section exceeds every control section;
- and mean IR − mean control ≥ 0.2.

The same reading is reported for INJ-PT.

**R2 · Transport loss is S3-weighted.** ρ(e_TPT, x) per section. Replicates when every IR section is below every
control section, with mean IR − mean control ≤ −0.2.

**R3 · Gradient survival (notebook 91's I2).**
- Thirds of x act as S1, S2 and S3.
- A gene is control-graded for a step when its contrast (mean log1p CP10k of the later third minus the earlier
  third) is ≥ 0.5 in absolute value, with the same sign, in every control section.
- **Survival:** the slope through the origin of the IR section's contrasts on the mean control contrasts, over those
  genes.
- **Replicates** when the S2→S3 survival is below the S1→S2 survival in every IR section.

**R4 · Notebook 91's location calls** (`location change`, sign of `Δ late`).
- The IR − control change in the last third of x minus the change over all PT-dominant spots, means over sections.
- **Replicates** when the sign agreement is ≥ 0.6 with a binomial p ≤ 0.05.

**Descriptive.**
- e_INJ, e_TPT and e_CYC along x in 10 bins per section, including the late-S1 dip of notebook 99.
- ρ(x, depth) per section.
- The NewPT weight along x.

**Overall reading.** Our positional injury claim **replicates externally** when R1 replicates at Cheung 12 h and at
least one of Cheung day 1 or Gharaie 24 h.

## 4 · Outputs and figures

- **Tables:**
  - `spots.csv` (dataset, section, condition, time, x, x_seg, depth, weights, scores, residuals);
  - `markers.csv`, `positional_genes.csv`;
  - `r1_r2.csv`, `r3_survival.csv`, `r4_location.csv`, `profiles.csv`;
  - `how_the_protocol_was_read.csv`, `run_record.json`.
- **Figures (plots only):**
  - 105a: a section map coloured by x and by e_INJ, one control and one 12 h section;
  - 105b: e_INJ and e_TPT along x per section;
  - 105c: ρ per section for R1 and R2;
  - 105d: survival per step.

## 5 · Caveats

- 55 µm spots mix cells. PT-dominant spots in the outer stripe still carry thick-limb RNA, which INJ-PT guards
  against.
- One section per mouse. Day 1 has two shallow sections.
- Gharaie's GEO areas are 2 of the paper's 4 per condition. That "1" and "2" are different mice is likely but not
  stated.
- Ranking x within a section assumes each section samples the PT axis similarly. A section missing the outer stripe
  compresses S3.
- The 19 and 30 min clamps differ from ours.
