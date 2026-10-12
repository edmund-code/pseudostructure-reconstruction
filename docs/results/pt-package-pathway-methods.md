

## Notebook 75 results

`analysis/notebooks/75_pt_package_pathway_methods.ipynb` (logic version `75.package.1`, by the coordinator's
instruction) writes to `results/pt_package_pathway_methods/`. It followed the protocol above, and every guard passed:
- G1–G2: export, fingerprints, mt exclusion, and 71's universe, sets and units;
- G3: 71's six gene models were read from its cache and reproduce its z;
- G4: the copied rows reproduce 71's calls and verdicts;
- section 7's guards: the Lake pseudobulk equals 68's cache, L equals 71's S3 level, and the native stage-01 matrices
  hold the export counts.

Its section 9 lists how the protocol was read. The new module is `pseudospace/discrete_comparators.py` (segment
F-test, family ORA, swap rule), tested by `tests/test_discrete_comparators.py`.

**Method table** (`tables/method_table.csv`). Each cell gives calls; relabeled calls; calibration share; verdict under
the one rule.

| Row | Human vs mouse | AKI vs control |
|---|---|---|
| Positional F-test · cameraPR (primary, 71) | 83; 1, 0; 0.071 · pass | 103; 13, 29; 0.120 · fail (calibration) |
| Segment F-test · cameraPR (new; coordinate-free S1/S2/S3 steps) | 76; 0, 0; 0.055 · pass | 38; 24, 17; 0.115 · fail |
| Segment steps, NB joint test (sensitivity) | 34; 0, 0; 0.069 · pass | 0; 0, 0 · fail |
| NB joint test on `T_spatial` (66) | 52; 0, 2; 0.072 · pass | 16; 0, 0; 0.080 · fail (calibration) |
| cameraPR on rank-normal raw `T_spatial` (66) | 115; 53, 50; 0.132 · fail | 153; 106, 107; 0.188 · fail |
| cameraPR use.ranks / GSEA rank-normal / GSEA raw (66) | 97; 28, 19 / 152; 95, 79 / 0; 25, 9 · all fail | all fail |
| GSEA on the positional z (71 S1) | 75; 2, 0; 0.144 · fail (calibration) | 34; 64, 43 · fail |
| Whole-PT DESeq2 + signed GSEA | 146; 189, 261; 0.271 · fail | 622; 287, 58; 0.219 · fail |
| S1/S2/S3 DESeq2 + signed GSEA | 197; 220, 286; 0.266 · fail | 727; 362, 383; 0.349 · fail |
| Whole-PT DESeq2 + ORA | 81; 0, 0; 0.012 · pass | 259; 0, 0; 0.000 · pass |
| S1/S2/S3 DESeq2 + ORA | 164; 0, 1; 0.018 · pass | 411; 0, 0; 0.005 · pass |
| Level / total arms (71 S1; cameraPR and GSEA) | all fail | all fail |
| Level gene list + ORA (71 S3) | 114; 0, 0; 0.000 · pass | — |
| Clustering, default null (71) | 94; 0, 0; 0.003 · pass | 459; 0, 0; 0.017 · pass |
| Curve modules, top 10% genes by z (56(c) twin) | 54; 0, 3; 0.043 · pass | 107; 104, 37 · fail |
| Curve modules, q ≤ 0.05 genes (61 B2 twin) | 59; 0, 0; 0.000 · pass (the relabelings have no significant genes) | 307; 0, 4; 0.018 · pass |
| ssGSEA + limma, shape / level (71 S2) | 866 / 1,208 of 1,513; pass | 1,117 / 1,178 of 1,469; pass |

**Material changes.**

| Rule | Old | New | Met |
|---|---|---|---|
| M1 · cameraPR on raw `T_spatial` passes | fails on both | fails on both | no |
| M2 · NB joint test changes verdict | passes on both | human vs mouse passes; AKI fails calibration (0.080) | **yes** |
| M3 · average-level GSEA passes, or DE + ORA fails (human vs mouse) | fail, fail; pass, pass | fail, fail; pass, pass | no |
| M4 · S1/S2/S3 ORA no longer probe-exposed, or the primary exposed | 0.60 | 0.71; primary 0.17 | no |
| M5 · segment F-test fails on human vs mouse | (Gaussian step screen 50; 0, 0) | 76; 0, 0; 0.055 · pass | no |
| M6 · level reading changes | (i) 17/27 fail; (ii) 90.0%; (iii) 0.60; leak 0.45 | (i) 18/27 fail; (ii) 89.8%; (iii) 0.60; leak 0.46 | no |
| M7 · primary shares > 0.3 of its redundancy clusters with a discrete or average-level row | 0.07–0.17 (sets) | total arm on AKI 0.42 | **yes** |

- **M2.** On AKI the joint test's relabeling-fitted VIF (a = 1.57, median VIF 2.33) leaves a calibration share of 0.080,
  just above the limit. Its cross-fitted share is 0.081 on AKI and 0.074 on human vs mouse.
- **M7 is met through one row only:** the total arm (`T_total` = offset + shape), which was labelled average-level
  before the run. `T_total` contains the positional term by construction. Every DE row and every level row stays at or
  below 0.23 over redundancy clusters.

**What the rows show.**
- **The segment F-test** (the discrete twin of the primary) calls 76 sets on human vs mouse. 71 of them are also
  primary calls (Jaccard 0.81); the spline adds 12 sets and the steps add 5. On AKI it fails both parts of the rule. The
  per-gene tables for notebook 74 are `gene_statistics_segment_{ds}.csv` and `pathway_tests_segment_{ds}.csv`.
- **Packaged tests on raw `T_spatial` are not specific on the package either.** The relabeled cameraPR z are inflated
  by a = 1.65 (human vs mouse) and 2.24 (AKI). Ranking within expression strata brings human vs mouse to 82; 6, 10;
  0.087, which still fails calibration.
- **Average-level screens.**
  - Signed GSEA on pseudobulk DESeq2 calls more sets under the relabelings than for the condition, on both datasets.
  - Thresholded DE + ORA passes the rule but is probe-exposed on human vs mouse: whole PT 0.63, S1/S2/S3 0.71.
  - Its calls barely overlap the primary's (Jaccard 0.09–0.16 over sets). 95–97% of the sets only a DE row calls lean
    to the offset, against 0% of the sets only the primary calls (72% of all sets lean to the offset).
- **Probe exposure, human vs mouse.**

  | Rows | Exposure |
  |---|---|
  | Positional rows (GSEA on the positional z: 0.39) | 0.14–0.25 |
  | Segment F-test | 0.20 |
  | Curve modules | 0.35–0.44 |
  | Clustering | 0.69 |
  | DE rows | 0.63–0.73 |
  | Level arms | 0.55–0.85 |
  | ssGSEA | 0.22–0.28 |

- **Gene-first.**
  - Clustering stays specific; its calls overlap the primary's at Jaccard 0.17 over redundancy clusters (human vs
    mouse).
  - Curve modules of the NB Δ(x) pass on human vs mouse (k = 4 and 3). The top-10% twin fails on AKI.
- **Programs.** The primary's 83 human vs mouse calls form 43 programs, and its 103 AKI calls 56. Figure 4e
  (`figures/fig4e_*_human_vs_mouse`) shows the largest program, "Peroxisomal Lipid Metabolism" (6 calls, 66 genes),
  beside the q ≤ 0.05 module 2 (1,113 genes, 36 enriched sets) and clustering cluster D4 (198 genes, 58 called sets).
- **Level control** (68 on the package; reading unchanged: fail, pass). (i) 18 of 27 literature controls have the
  expected sign (66.7%; Wilson 48–81%), which fails. (ii) 3,575 of 3,982 atlas-strong genes agree (89.8%), which
  passes. (iii) Spearman 0.60 [0.58, 0.61] against an atlas ceiling of 0.97. Families: 196 of 409 (48%). Probe count
  still leaks into the level (Spearman 0.46). Of the 9 failing controls the atlas also disagrees with 6, as in 68.

**Compute.** 605 s wall on 24 workers for the cold development run; the executed notebook reran from the stage cache.
