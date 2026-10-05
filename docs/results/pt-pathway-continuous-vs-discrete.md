# Workstream B · Phase 2 results: pathway approaches and what continuous analysis adds

Written 2026-10-05. Protocol: `plan_phase1.md`, section 8, approved by the coordinator. Notebooks 56–60
are in `analysis/notebooks/`; results are in `results/pt_pathway_continuous/`. Coordinates: SCF13 (the
paper's), DPT13, SCF13-ED (equal depth) and workstream A's LCP (`results/pt_reconstruction_v3/lcp_coordinate.csv`).
Cohort: 2 male control mice and 2 cortex sections from 1 male human donor. All results are
descriptive. Every departure from the frozen protocol is listed in section 7.

## 0 · Summary

1. **"More pathways" fails.** **"More information per program" holds only for within-segment
   gradients.**
   - The order in which genes switch on inside mouse S2 (onset order, (iv)) does not separate from what
     pure segment steps produce. The external check also does not replicate it on SCF13: partial ρ 0.10.
   - Species boundary shifts (iii) do not reach the claim threshold (5 genes against 10).
   - Within-segment slopes (i) do replicate between specimens beyond both step nulls. They do so in more
     cells the less positional noise a coordinate has (LCP: 4 of 6 cells; the others 2). They partly
     agree with external mouse nuclei (partial ρ 0.06–0.19 after removing segment steps).
2. **The beyond-steps pathway signal splits in two.**
   - Beyond **label** steps: label–coordinate discordance alone reproduces 51% of it on SCF13 (62–71%
     on the other coordinates). The pre-registered verdict is discordance.
   - Beyond **coordinate-threshold** steps: 67 calls against null medians of 9–11. That part is not a
     step artefact.
3. **The ORA finding stands, and it changes R4.**
   - Thresholded ORA on segment DE genes calls 118 pathways for species and 0 and 1 under relabeling.
   - The power-matched external nulls also give 0 calls in all 750 same-species splits, with
     cross-species NCDR 0. So ORA is specific by the pre-registered rule.
   - But 60% of its species calls disappear on genes with equal probe counts (118 → 49). It is
     probe-exposed.
4. **Self-contained pathway activity scores tested at specimen level are specific (0 and 0 relabeled)
   but call 85–90% of pathways.** Curve modules fail specificity (NCDR 0.40). Location enrichment
   calls nothing.

## 1 · Notebook 57: is "beyond S1/S2/S3 steps" within-segment biology?

The semi-synthetic arms permute expression rows within specimen × segment, which keeps segment means
and inter-gene correlation and removes all within-segment structure.
- **Arm A** permutes within coordinate-threshold segments (truth = steps on the coordinate; labels are
  noisy).
- **Arm B** permutes within labels (truth = label steps; the coordinate is noisy).

Ten replicates per arm, each with the full joint pathway test.

| Coordinate | Test A observed (beyond label steps) | Arm A median (range) | Ratio (rule ≥ 0.5 → discordance) | Test B observed (beyond coordinate steps) | Arm B median | Arm A median (test B) | Relabeled A; B | Label agreement mouse / human |
|---|---:|---|---:|---:|---:|---:|---|---|
| SCF13 | 79 | 40 (31–51) | **0.51** | 67 | 11 | 9 | 6, 0; 1, 0 | 0.92 / 0.85 |
| DPT13 | 41 | 25.5 | 0.62 | 29 | 2.5 | 7.5 | 7, 6; 0, 0 | 0.87 / 0.85 |
| SCF13-ED | 81 | 53.5 | 0.66 | 84 | 27.5 | 26 | 10, 9; 5, 0 | 0.91 / 0.84 |
| LCP\* | 91 | 64.5 | 0.71 | 84 | 47 | 27 | 1, 0; 0, 0 | 0.96 / 0.87 |

\*On LCP the species are not registered (median human S2 position 0.70 against mouse 0.42), so the
species tests there mostly measure misregistration. Report them, but do not read them.

- **The joint test has a floor in these exact nulls.** Arm B test A and arm A test B should call
  nothing, yet give 3–31 calls. Under relabeling the same tests gave 0–10. About half of the arm-A
  test-A count is therefore this floor, not label noise. The pre-registered ratio does not subtract it;
  the verdict stands as written.
- **Within-segment pathway candidates** (tests A and B both call, on SCF13 and at least one other
  coordinate): 48. Among them are 28 of the robust 36 and 20 of the primary 22.
  - The stricter rule adds that both observed z values exceed every simulated z. On SCF13 it keeps 20
    pathways, 17 of them also cross-coordinate. Six of those 17 are primary.
  - Called by both tests on all three original coordinates: 18 pathways, 8 primary (Bile Acid
    Metabolism, Chemical carcinogenesis, Estrogen Response Early and Late, Fatty Acid Metabolism,
    Metabolism of xenobiotics by P450, Peroxisome, Steroid hormone biosynthesis).
- **Not addressed by these nulls.** A non-anatomical axis shared by the coordinate and the genes (depth,
  co-expression with the construction genes) would also pass test B.

## 2 · Notebook 58: descriptors per specimen, replicated between specimens

Each specimen's segment units come from its own transitions. Genes are fold-protected: genes that
built the coordinate are excluded (611 PCA genes for SCF13 and DPT13, 549 for SCF13-ED, 920 landmarks
for LCP).

### (iv) Onset order inside mouse S2 (primary for F2 against F3)

The statistic is the partial Spearman of onset in Ctrl1A2 against onset in Ctrl1A4, given each
specimen's switch fraction r = (S2 − S1)/(S3 − S1).

| Coordinate | Genes | Partial ρ (95%) | Step-null max (arms A, B) | Post hoc: given r from coordinate segments | Its step-null max |
|---|---:|---|---:|---:|---:|
| SCF13 | 169 | 0.36 (0.16–0.52) | 0.49 | 0.35 | 0.45 |
| DPT13 | 156 | 0.27 (0.08–0.47) | 0.31 | 0.19 | 0.21 |
| SCF13-ED | 38 | 0.64 (0.23–0.86) | 0.24 | 0.55 | 0.22 |
| LCP | 29 | −0.14 (−0.60–0.36) | 0.24 | −0.12 | 0.31 |

- **Pre-registered verdict: inconclusive.** The count-split arm could not be evaluated: notebook 35's
  half-A coordinate spreads mouse S1 over half the axis, leaving only 6 monotone genes per mouse.
- **Added control.** Pure steps plus a spline give onsets that depend on step sizes in a way that r does
  not capture linearly. On SCF13 the null median is 0.20 and the maximum 0.49. The observed 0.36 sits
  inside that range, and so does DPT13.
- Only SCF13-ED exceeds every null, on 38 genes.

### (i) Within-segment gradients

The statistic is the Spearman, across fold-protected genes, of each gene's slope t inside a segment in
one specimen against the other specimen of the same species. Filled = exceeds every step-null
replicate.

| Coordinate | Mouse S1 | Mouse S2 | Mouse S3 | Human S1 | Human S2 | Human S3 |
|---|---|---|---|---|---|---|
| SCF13 | 0.18 | **0.36** (null ≤ 0.15) | 0.30 (null 0.31) | **0.44** (0.38) | 0.18 | 0.12 |
| DPT13 | 0.10 | 0.36 (null 0.48) | **0.43** (0.27) | **0.28** (0.14) | 0.21 | 0.13 |
| SCF13-ED | 0.13 | 0.29 (0.30) | **0.52** (0.12) | **0.21** (0.15) | 0.17 | 0.14 |
| LCP | **0.44** (0.14) | **0.37** (0.10) | **0.39** (0.10) | **0.44** (0.22) | 0.24 | 0.18 |

- **Pre-registered claim on SCF13.** The rule needs ρ ≥ 0.3 and a count-split ρ ≥ 0.2. Mouse S2 meets
  it (0.36; count split 0.20), and so does mouse S3 (0.30; count split 0.37). Of the two, only mouse
  S2 also exceeds the step nulls.
- **Post hoc: slopes inside coordinate-threshold segments.** They exceed both nulls in nearly every
  cell on every coordinate: mouse S2 0.30–0.37, S3 0.33–0.58, human S1 0.19–0.50.
- **Resolution budget.** Over the 24 coordinate × cell combinations, positional noise ÷ within-segment
  spread correlates with slope agreement at Spearman −0.31.

### (ii) Interior extrema (peak or trough inside a segment, replicated in both specimens)

| | SCF13 | DPT13 | SCF13-ED | LCP | Step-null max (SCF13) |
|---|---:|---:|---:|---:|---:|
| Mouse S1 / S2 / S3 | 48 / 34 / 152 | 0 / 25 / 0 | 94 / 10 / 174 | 10 / 6 / 154 | 14 / 37 / 232 |
| Human S1 / S2 / S3 | 61 / 0 / 13 | 11 / 0 / 0 | 60 / 0 / 29 | 73 / 0 / 1 | 4 / 0 / 68 |

- A spline fitted to steps produces most of the mouse counts, mainly at the S3 end of the support.
- **Human S1 is the only cell beyond the nulls on every coordinate**, including equal depth.
- Human S1 is also where the coordinate's order depends most on read depth, and the two human sections
  come from one donor.
- The cross-coordinate rule (at least 2 coordinates including SCF13-ED) gives 44 human and 150 mouse
  genes. Read these counts against the nulls.
- The mouse S2 peak genes (Car4, Cndp2, Cyp2e1, Cyp4b1, Inmt) are notebook 07's "mid-PT" module. They
  do not exceed the S2 step null.

### (iii) Species shifts of a switch relative to each species' own boundaries

- **Rule.** Monotone in all four specimens, sign agreement in all four cross-species pairs, and
  \|shift\| beyond the registration null. The registration thresholds come from notebook 35's refits
  (S1 0.59, S2 1.32, S3 0.84 segment units).
- **Result.** 5 genes pass on SCF13 (Gatm, Idh1, Pdk2, Slc5a1, Slc6a18), 2 on DPT13, 1 on SCF13-ED,
  1 on LCP, and no program.
- Only Slc5a1 keeps its sign on DPT13.
- **Not claimed** (the rule needs ≥ 10 genes or ≥ 3 programs).

## 3 · Notebook 59: external adult male mice (Census, 6 donors, per-donor diffusion pseudotime)

The external pseudotime is built per donor from HVGs that exclude every tested gene. S1 < S2 < S3
holds in all 6 donors.

| Coordinate | Genes with an external onset | Partial ρ, external onset given our r (donor bootstrap) |
|---|---:|---|
| **SCF13 (primary)** | 73 | **0.10 (0.01–0.24)** → not replicated (rule: ≥ 0.2) |
| DPT13 | 62 | 0.29 (0.16–0.63) |
| SCF13-ED, LCP | 16, 13 | too few genes |

The abandon rule for F2(a) needs < 0.1 both here and in notebook 58. Notebook 58 gave 0.36, so F2(a)
is not formally abandoned. It is not supported either.

**Post hoc: within-segment slope agreement with external nuclei** (our mice against the external donor
median).

| Coordinate | S1 | S2 | S3 |
|---|---|---|---|
| SCF13 | 0.29 → 0.12 | 0.36 → 0.17 | 0.20 → 0.12 |
| DPT13 | 0.11 → 0.08 | 0.37 → 0.17 | 0.13 → 0.06 |
| SCF13-ED | 0.18 → 0.10 | 0.34 → 0.18 | 0.14 → 0.10 |
| LCP | 0.40 → 0.17 | 0.33 → 0.18 | 0.31 → 0.19 |

Each cell gives the raw agreement, then the agreement after removing the external S2 − S1 and S3 − S2
steps (label slopes; slopes inside coordinate segments are similar).
- About half of the agreement is shared step structure leaking into within-label slopes. The rest,
  ρ 0.06–0.19, is consistent with a within-segment gradient present in external nuclei too.
- Genes are treated as independent, so read the effect sizes, not significance.

## 4 · Notebook 56: the pathway approaches side by side (SCF13; 1,513 pathways)

Full table: `56_taxonomy/tables/f_taxonomy.csv`. Figure: `56_taxonomy/figures/approaches_side_by_side.pdf`.
Rule = NCDR < 0.25 and each relabeling calls < 5% of tested pathways.

| Approach | Returns | Species | Relabeled | Rule | ∩22 / ∩36 | External |
|---|---|---:|---|:-:|---|---|
| Whole-PT DE + signed GSEA | yes/no, sign | 158 | 193, 255 | fail | 3 / 8 | NCDR 1.34 |
| S1/S2/S3 DE + signed GSEA | sign per segment | 204 | 343, 255 | fail | 7 / 17 | NCDR 0.69 |
| **Whole-PT DE + ORA** | yes/no, sign | 80 | 0, 1 | pass | 6 / 13 | **specific (NCDR 0; 99% of draws)** |
| **S1/S2/S3 DE + ORA** | sign per segment | 118 | 0, 1 | pass | 9 / 18 | **specific; probe-exposed** |
| S1/S2/S3 DE + ORA, equal-probe genes | as above | 49 | 0, 0 | pass | 7 / 13 | – |
| Segment markers + ORA (species pooled) | segment identity | 21 | not run | – | 1 / 1 | – |
| Coordinate-free S1/S2/S3 steps (joint) | step pattern | 50 | 0, 0 | pass | 20 / 28 | NCDR 0.05; 80% |
| Steps on the spline (joint) | step pattern | 16 | 0, 0 | pass | 8 / 9 | – |
| 6 coordinate bins (joint) | 6-bin profile | 49 | 17, 0 | pass | 22 / 29 | – |
| **T_spatial, joint (primary)** | D(s), peak, shape | 60 | 2, 0 | pass | 22 / 36 | – |
| T_spatial, matched / GSEA / size-only | yes/no | 116 / 311 / 344 | 17,13 / 34,24 / 59,55 | pass | 22 / 36 | – |
| Beyond label steps (joint) | within-segment part | 79 | 6, 0 | pass | 22 / 35 | notebook 57 |
| T_level / T_total (joint) | average | 0 / 0 | 46,2 / 72,15 | fail | – | – |
| Pathway-score GAM, structure level (offset / shape) | activity curve | 1,444 / 1,509 | 986,1317 / 67,495 | fail | – | – |
| **Pathway activity (mean-z / AUCell-like), specimen × bin split plot** | activity curve per specimen | 1,355 / 1,291 | 0, 0 | pass | 21–22 / 35–36 | – |
| Gene split plot (8 bins), matched | yes/no | 41 | 0, 0 | pass | 7 / 15 | – |
| Curve modules (k = 8) + ORA | shape families | 44 | 23, 12 | **fail** (NCDR 0.40) | 9 / 20 | LOSO ARI 0.82 |
| Divergence-peak concentration / location | where | 0 / 0 | 0, 0 | – | – | 505 testable |
| tradeSeq equivalents (association, start vs end, early S1) | gene ranks | Spearman with T_spatial 0.48–0.69 | – | – | – | – |

**Overlap** (Jaccard of species calls):
- T_spatial against coordinate-free steps 0.75, against 6 bins 0.76, against beyond-steps 0.62.
- T_spatial against S1/S2/S3 ORA 0.17, against S1/S2/S3 GSEA 0.10, against whole-PT GSEA 0.07.
- Whole-PT and S1/S2/S3 GSEA against each other: 0.74.
- The positional family and the average family are nearly disjoint. Within each family the approaches
  largely agree.

**The ORA finding (for R4).**
- **Our cohort, from saved DESeq2.** S1/S2/S3 ORA calls 118 pathways for species and 0 and 1 under
  relabeling; whole-PT ORA calls 80, also 0 and 1.
- **Power.** The relabeled design has 1 residual df against 2 for species. So the external
  same-species 2 + 2 splits (design ~group, 2 residual df, moderated t) were run.
  - They called 0 pathways in every one of 750 splits: the median DE gene count was 0.
  - Cross-species draws: median 27 species calls against 0 relabeled (NCDR 0; 99% of draws meet the
    rule).
  - By the pre-registered rule, ORA is **specific**.
- **Probe panels.** On genes with equal probe counts, 49 of the species calls remain. 60% are lost,
  so by the rule ORA is **probe-exposed**.
- **What this means.** ORA answers an average-level question with gene-level thresholds. Relabeling
  cannot show any average-level screen false. But its species calls depend heavily on genes whose
  probe counts differ between panels.
- **R4 as written is too general.** The sentence "average-level pathway screens cannot be shown
  specific by relabeling" holds for rank-based screens (GSEA, matched rank-AUC), not for thresholded
  ORA.

## 5 · Notebook 60: the resolution budget

Positional noise is from gene-fold refits: SD in segment units, and noise ÷ within-segment spread.
Mean over the two specimens of each species.

| Coordinate | Mouse S1 | Mouse S2 | Mouse S3 | Human S1 | Human S2 | Human S3 | Cells beyond step nulls (i) | (iv) onset ρ | External (iv) |
|---|---|---|---|---|---|---|---:|---|---|
| SCF13 | 0.56 (1.44) | 0.30 (0.83) | 0.30 (1.02) | 0.41 (1.49) | 0.57 (1.04) | 0.57 (1.42) | 2 | 0.36, null | 0.10 |
| DPT13 | 0.34 (1.04) | 0.20 (0.55) | 0.26 (0.88) | 0.15 (0.76) | 0.34 (0.62) | 0.28 (0.64) | 2 | 0.27, null | 0.29 |
| SCF13-ED | 0.51 (1.33) | 0.40 (0.96) | 0.41 (1.70) | 0.31 (1.17) | 0.83 (1.23) | 0.61 (2.09) | 2 | 0.64 | – |
| LCP | **0.08 (0.27)** | **0.04 (0.16)** | **0.06 (0.25)** | 0.20 (0.65) | 0.26 (0.44) | 0.16 (0.35) | **4** | −0.14 | – |

- LCP's model posterior SD is 0.02–0.03 segment units in mouse and 0.06–0.16 in human. It is narrower
  than its fold noise, so it understates positional error.
- Count-split noise is in `resolution_budget.csv`.
- Figure: `60_resolution_budget/figures/resolution_budget.pdf`.
  - Left panel: noise ÷ spread against between-specimen slope agreement, per cell and coordinate.
  - Right panel: onset order against noise.
- **Reading.** Lower positional noise goes with more replicated within-segment gradient information.
  Onset order does not follow this pattern.
- **Limits.** These are four coordinates and 24 cells; this is not a test. LCP was abandoned by
  workstream A's own rules (held-out prediction), and it does not register the species.

## 6 · What this means for the framing

- **F2 as pre-specified does not survive.** It needed onset order beyond segment means, replicated
  between specimens and externally: SCF13 falls inside the step null, and the external partial ρ is
  0.10. Boundary shifts are 5 genes.
- **What survives.**
  1. Reproducible within-segment gradients. They exceed step nulls, improve as positional noise falls,
     and are partly shared with external nuclei.
  2. A pathway-level positional signal beyond coordinate-threshold steps (test B), which step artefacts
     do not explain.
  3. Steps on the coordinate segment the PT more cleanly than the reviewed labels: a large part of
     "beyond label steps" is that.
- **Recommended frame: F3 plus a narrow F2.**
  - At present resolution, a continuous analysis matches a 3–6-bin analysis in what it calls.
  - Its added information is within-segment gradient direction. That information grows as the
    coordinate's positional noise falls, which is the link to workstream A.
  - Onset order and boundary shifts are not supported.
- **Open risk.** The step nulls cannot rule out a non-anatomical axis that the coordinate shares with
  the genes: read depth, or co-expression with the construction genes. The equal-depth arm and the
  external agreement after removing steps (0.06–0.19) only partly address it.

## 7 · Deviations from the frozen protocol

1. **Notebook 57.** The semi-synthetic exact nulls show a joint-test floor of 3–31 calls. This was not
   anticipated; it is reported beside the pre-registered ratio, which is unchanged.
2. **Notebook 58, step nulls (arms A and B).** Added before parts (i)–(iii) were computed, but after
   (iv) had been seen. The pre-registered verdicts are reported unchanged beside them.
3. **Notebook 58, post hoc variants.** (iv-b), r from coordinate-threshold segments, and (i-b), slopes
   inside coordinate-threshold segments, were added after the first null comparison.
4. **Notebook 58, count split.** The arm was not evaluable for (iv) (6 monotone genes per mouse), so
   the rule returns "inconclusive". For (i) it was evaluated.
5. **Notebook 58 (iii), registration null.** Implemented as the change in the human − mouse segment gap
   across notebook 35's SCF13 refits, divided by the human segment length, 95th percentile.
6. **Notebook 56(d).** Peaks of \|Δ(s)\| replace onsets, because a difference curve has no onset.
7. **Notebook 56(a).** The external pools use notebook 45's moderated t (its declared DESeq2
   substitute); our cohort uses the saved DESeq2.
8. **Notebook 56(b).** The "AUCell-like" score is a top-5% rank-recovery score, not the AUCell package.
9. **Notebook 59.** Human was not tested: no external human dataset resolves S2. The post hoc slope
   section was added after the onset result was seen.
10. **LCP.** Fold protection uses workstream A's 920 landmark genes. Species are not registered on it,
    so notebook 57 there is uninterpretable.

## 8 · Suggested changes to existing documents (for the coordinator)

- **R4 and Supplementary Note 2.** Restrict the "cannot be shown specific" sentence to rank-based
  screens. Add that ORA passes relabeling and external power-matched nulls but loses 60% of its calls
  on equal-probe genes.
- **R5 and Supplementary Note 1.** Replace "35 of 36 robust pathways carry position beyond the steps"
  with the decomposition:
  - beyond label steps: about half reproduced by label–coordinate discordance;
  - beyond coordinate steps: 67 calls against null medians of about 10;
  - 48 cross-coordinate candidates.
- **Outline.** Lead with "continuous ≈ fine discretisation in calls; adds within-segment gradients
  whose replication scales with coordinate resolution". Drop onset ordering and boundary shifts as
  claims.

## Files

**Code** (all new):
- `pseudospace/continuous_descriptors.py` and `pseudospace/coordinate_inputs.py`
- `tests/test_continuous_descriptors.py` (10 tests; full suite 422 passed)

**Notebooks** (all new, executed in place with outputs; stage output-free copies):
- `analysis/notebooks/56_pt_pathway_taxonomy.ipynb`
- `analysis/notebooks/57_pt_beyond_steps_identifiability.ipynb`
- `analysis/notebooks/58_pt_continuous_descriptors.ipynb`
- `analysis/notebooks/59_pt_external_onset_replication.ipynb`
- `analysis/notebooks/60_pt_resolution_budget.ipynb`

**Notebook sources and sensitivity runs** (`workstream-B scratch: `):
- percent-script sources: `nb56.py`–`nb60.py`
- executed copies: `nb57_executed_{dpt13,scf13_ed,lcp}.ipynb` and `nb58_executed_{dpt13,scf13_ed,lcp}.ipynb`
- logs: `run*.log`

**Results** (gitignored): `results/pt_pathway_continuous/`, with the subfolders `56_taxonomy/`,
`57_beyond_steps/`, `58_descriptors/`, `59_external_onsets/`, `60_resolution_budget/` and `inputs/`.

## Addendum (notebook 61, plan addendum 3): order-invariant programs and pre-test deduplication

- **Notebook 37's program grouping is superseded by the order-invariant rule.**
  - Average linkage with tied distances depends on input order. Notebook 37 clustered in library
    order; the new rule sorts by `pathway_id` first.
  - On the 36 robust pathways the count stays at 17, but 7 vitamin pathways change co-members.
    "Metabolism Of Vitamins And Cofactors" moves from the water-soluble/pantothenate program to the
    fat-soluble/retinoid program. Pantothenate and CoA biosynthesis, Vitamin B5 and Water-Soluble
    Vitamins stay together.
  - Notebook 37 itself is unchanged; the moves are in `robust36_program_moves.csv`.
- **Programs under the new rule:**
  - the primary 22 give 12 programs (the same grouping as before);
  - the 60 C1 calls give 25;
  - per method, pathways → programs: segment steps 50 → 25, six bins 49 → 25, whole-PT ORA 80 → 23,
    S1/S2/S3 ORA 118 → 43.
- **Pre-test deduplication** (member Jaccard ≥ 0.7, unions of connected components). The library goes
  from 1,513 sets to 1,160; at ≥ 0.5 it would be 853, and at ≥ 0.9, 1,395.

  | Method | Original library | Collapsed library |
  |---|---|---|
  | T_spatial joint test | 60 calls (2, 0) | 63 (3, 0) |
  | Segment-step screen | 50 calls (0, 0) | 49 (0, 0) |
  | Whole-PT ORA | 80 calls (0, 1) | 44 (0, 1) |

  - All 60 T_spatial calls are kept. There are 6 new calls, all single sets that cross the threshold
    because BH runs over fewer tests.
  - No loss comes from union dilution, and no specificity verdict changes.
  - The frozen rule reports "changes the substantive result", because 4 programs are gained
    (25 → 30). The cause is the number of tests, not the merging.
  - Caveat: chaining merges 67 cell-cycle and proteasome entries into one 429-gene set. Not adopted.
- **Decision.** The pre-test deduplication was tried (1,513 → 1,160 sets; no positional calls lost;
  whole-PT ORA 80 → 44) and dropped by the user's decision in favour of post-test driver-gene programs.
  It has been removed from notebook 61. The order-invariant post-test program rule is kept.
