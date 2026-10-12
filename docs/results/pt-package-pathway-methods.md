# Pathway methods and their specificity on the package pseudospace (notebook 75)

Part of the switch to the package pseudospace; see `pt-package-pseudospace.md`. The protocol was fixed before
any fit it governs.

## Notebook 75 protocol

Frozen on 2026-10-12, before any fit of notebook 75. Nothing in it has been computed on the package inputs. It
applies coordinator decisions D1 (package units and labels), D2 (NB with the mt-free offset; pseudobulk DE keeps its
own model), D5 and D6 (notebook 71's calls and post hoc arms are frozen; the 45-method grid retires), and the
coordinator's decisions O1–O5 on this draft (section 13).

**Why.** The paper's method table needs one row for every pathway method it cites, judged by one rule, on one set of
units. Notebook 71 froze the primary (the positional F-test with cameraPR), the clustering comparator and the
comparison page's post hoc arms. The remaining results of record still rest on SCF13, `mouse_only_v5` and notebook 13's
labels:
- notebook 66: packaged set tests on raw `T_spatial`, and the NB joint test;
- notebook 68: the level positive control;
- notebooks 31, 37, 45, 56 and 61: the conventional comparators, the coordinate-free step screen, probe exposure, the
  gene-first curve modules and the program heatmap (Figure 4e).

**Question.** On the package units, labels and coordinate, which pathway methods are specific against specimen noise
under the swap rule, what does each call, and how far do their calls overlap with the positional F-test? Does the
species level still recover known human–mouse differences? Everything is descriptive: two male control mice, two
cortex sections of one male human donor, two AKI mice. Specimens are the replicates.

**Labels.**
- Notebook `analysis/notebooks/75_pt_package_pathway_methods.ipynb`, logic version `75.methods.1`.
- Results in `results/pt_package_pathway_methods/` (`tables/`, `figures/`, `stage_cache/`, `run_record.json`).
- Datasets: every row runs on human vs mouse (hvm) and AKI vs control (aki) wherever it is defined (O5). hvm only:
  section 4 (probe exposure; AKI has one probe panel), section 7 (the human atlases and literature), 71's S3 row and
  the symmetric-class side bar of section 6.

### 0 · Results of record and where they go

| Source | Result of record (old number) | In 75 |
|---|---|---|
| 66 | cameraPR on rank-normal raw `T_spatial`: hvm 98; 46, 29; 0.121, aki 144; 87, 86; 0.170 (both fail) | **kept** (3.1) |
| 66 | cameraPR `use.ranks`; GSEA on rank-normal and on raw `T_spatial` | **kept** as sensitivities (3.1) |
| 66 | NB joint test, relabeling-fitted VIF: hvm 32; 3, 0; 0.067, aki 12; 0, 0; 0.072; cross-fitted 0.067, 0.075 | **kept** (3.2) |
| 66 | why the rank tests fail: implied inflation a 1.51, 2.13; within-stratum scores 0.082, 0.135 | **kept**, now a pre-specified diagnostic (3.1) |
| 66 | clustering with default, single-split and fair nulls | covered by 71 (`method_table.csv`) |
| 66 | overlap of cameraPR with clustering (Jaccard 0.15, 0.15) | **kept** (section 5) |
| 66 | set-curve calls; "what each finds alone"; worked examples | dropped: never cited; 71 sections 7–8 and the page replace them |
| 68 | (i) literature 17/27 fail; (ii) atlas 90.0% pass; (iii) ρ 0.60; (iv) families 48%; bands; probe leak 0.45; pathways ρ 0.70 | **kept** (section 7) |
| 71 | level/total cameraPR and GSEA, ssGSEA + limma, level list + ORA, atlas replication, winsorizing | covered by 71 (S1–S9); rows copied (3.6) |
| 31/37/56(a) | whole-PT DESeq2 + signed GSEA 158; 193, 255; S1/S2/S3 204; 343, 255 (both fail) | **kept** (3.4) |
| 56(a) | DESeq2 + ORA: whole PT 80; 0, 1; S1/S2/S3 118; 0, 1 (pass) | **kept** (3.4) |
| 56(a) | probe exposure: S1/S2/S3 ORA 118 → 49 on equal-probe genes (60% lost); whole PT → 27 | **kept**, extended to every row (section 4) |
| 56(a) | external power-matched ORA nulls (750 splits: 0 calls; cross-species NCDR 0) | unchanged (atlases only; no unit, label or coordinate of ours) |
| 56(a) | segment markers + ORA (21, no relabeling) | dropped: no swap result; not cited |
| 45/56/61 | coordinate-free step screen 50; 0, 0 (Gaussian joint test) | **kept** as our cohort's row, now the NB segment F-test (3.3; O1, O2) |
| 45 | our-cohort anchor among the external pools ("at the external median") | dropped (O1): the external pools are Gaussian, so an NB row cannot be placed among them |
| 45 | AKI positive control (13 of 14 markers peak in S3; coordinate screen NCDR 0.17) | notebook 76 |
| 56(b) | activity-score split plots: mean-z 1,355, AUCell-like 1,291; 0, 0 (specific, unselective) | covered by 71 S2 (ssGSEA + limma split plot: 866 and 1,208; 0, 0; same family and verdict). Scores are not counts, so D2 defines no NB version |
| 56(c) | curve modules (k = 8, top 10% genes) + ORA: 44; 23, 12; NCDR 0.40 (fail) | **kept**, NB curves (3.5) |
| 61 B2 | modules of significant genes (k = 4) + ORA: 23; 0, 0 (Fig 4e gene-first side) | **kept**, NB curves (3.5) |
| 61 | programs (60 calls → 25) and the program heatmap, pathway-first vs gene-first (Fig 4e) | **kept** on 71's calls, both datasets; 75 owns them (section 6; O4) |
| 61 | quadrants of the continuous method against steps (47/3/13), S1/S2/S3 ORA (26/92/34), whole-PT GSEA (14/144/46); offset-lean share of ORA-only calls (66 of 92) | **kept** (section 5) |
| 37/61 | "continuous adds 13 over steps" (near-miss reading); steps on the coordinate spline (16; 0, 0); six bins | notebook 74 (Fig 3e), which runs after 75 and reads 75's segment F-test (O2) |
| 56(d) | divergence-peak concentration and location (0, 0) | dropped from 75 (a continuous descriptor; 74 if wanted) |

### 1 · Inputs and reuse

**Built** exactly as notebook 71 builds them: `package_inputs.read_pt_units`, `hvm_inputs` and `aki_inputs` on
`results/package_pseudospace/{hvm,aki}/pt_units.h5ad`, with the 13 ortholog-map `mt-` genes excluded, the common
support, the position rule of 71's section 3, the 2% detection rule, the gene sets (10–300) and the 3-bin matching
strata. `partition_groups` gives the condition split and the two balanced relabelings.

**Read, never refit** (notebook 71, `results/pt_pathways_package_pseudospace/`):
- `stage_cache/gene_model_{ds}_{label}`: `T_spatial`, flags, α, Δ(x), the grid and the specimen curves of all three
  partitions;
- `tables/gene_statistics_{ds}.csv` (moderated F, z, q), `pathway_tests_{ds}_cameraPR_{label}.csv`,
  `redundancy_clusters_{ds}.csv`, `method_table.csv`, `s1_pass_table.csv`, `s2_ssgsea_pass_table.csv`,
  `s3_level_gene_table.csv`, `s3_level_variants_summary.csv`, `gene_sets_{ds}.json`, `universe_{ds}.csv`;
- `stage_cache/method_b_{ds}_{contrast}__null_{null}`: the clustering gene groups and clusters (section 4 only).

**Other inputs** (unchanged files): `results/pt_literature_deep/probe_counts.csv`; for section 7, everything notebook
68 reads (the literature controls, sha256 `f1aa7184…6c4f`; the housekeeping list; the Lake v1.5 h5ad; the Census
pseudobulks and family fetch; notebook 62's family catalogue and inventories; notebook 64's dropped pairs), plus
`results/package_pseudospace/runs/hvm/01_units/{sample}_units.h5ad` (native symbols and Ensembl ids, all genes) and,
for the heatmap side bar, notebook 73's `tables/43/gene_classes_symmetric.csv`.

**Guards (not cached).**
- G1 = 71's G1–G3: the export contract, the provenance (notebook 70, package 59da04b, pass-1 fingerprints
  `efafc9ac2e99` and `bb4a80c1b484`) and the mt exclusion.
- G2: the rebuilt inputs are 71's: genes equal `universe_{ds}.csv`, sets equal `gene_sets_{ds}.json`, unit counts equal
  `inputs.csv`.
- G3: each read payload is 71's run: its genes, `T_spatial` and flags equal `gene_statistics_{ds}.csv` (rtol 1e-12),
  and recomputing 71's moderated F from it and 71's `LR_spec` reproduces z. A purged cache stops the notebook; it never
  refits silently.
- G4: the rows copied from 71 reproduce 71's verdicts under section 2's code.

### 2 · One pass rule for every row

Every row, on each dataset where it is defined, is judged by one rule. That includes the rows copied from 71 (G4) and the
conventional DE rows. Each row has three runs: the condition split and the two balanced relabelings
(`partition_groups`). A row **passes** when all of these hold:
- it makes calls (> 0) on the condition split;
- NCDR < 0.25 (mean relabeled calls ÷ condition calls);
- each relabeling calls < 5% of the tested sets;
- calibration share ≤ 0.075: the share of relabeled per-set calibration p ≤ 0.05, pooled over both relabelings.

No row has its own rule or its own threshold. What differs by row is only how each test's per-set p and calls are
read:

- **Calibration p**, by kind of test (71's conventions):
  - one-sided set tests (cameraPR Up calls, GSEA NES > 0 calls, the joint test): `p_one_sided`;
  - signed set tests (calls in either direction): the test's own p (cameraPR two-sided p; GSEA nominal p);
  - several tests per set (S1/S2/S3 segments; ORA lists; modules; clustering lists): the smallest p × the number of
    tests, capped at 1 (notebook 61's list-based convention).
- **Gene universe and sets** are 71's for every row, DE comparators included (O3): `universe_{ds}.csv` and
  `gene_sets_{ds}.json`. Tested sets are counted per row: a set is tested where 10–300 of its members carry that row's
  statistic.
- **Calls in redundancy clusters** use 71's clusters for every row.

### 3 · The rows

**3.1 · Packaged tests on raw `T_spatial`** (66, method A). Read `T_spatial` of each partition from 71's cache.
- Primary: `camera_pr(rank_normal_scores(T_spatial))`, `inter.gene.cor` 0.01, `use.ranks` FALSE; call = Up and FDR ≤ 0.05.
- Sensitivities: `use.ranks` TRUE; `rank_first_gsea` on rank-normal and on raw `T_spatial` (1,000 permutations,
  seed 12, sizes 10–300).
- Diagnostics (66's section 5.2, now pre-specified): the inflation a and ρ implied by the relabeled cameraPR z
  (`decoy_inflation`); cameraPR on scores rank-normalised within the 27 matching strata; the Spearman of `T_spatial`
  with mean expression.

**3.2 · NB joint test** (66, reference row). `matched_pathway_tests` on `T_spatial` (9,999 matched sets, seed 12,
the builder's strata), `normal_matched_q`, the VIF from `decoy_inflation` on the two relabelings' matched z, then
`joint_matched_test`; call = effect > 0 and q ≤ 0.05. Also reported: the cross-fitted calibration share (VIF fitted on
one relabeling, judged on the other). 66's guard a ≈ 1.04 is dropped: it was specific to SCF13.

**3.3 · The coordinate-free step screen: segment F-test** (45/56/61's 50; 0, 0 → NB).
- **Old.** Gaussian `T_steps` (level → steps on `specimen_null.segment_designs`, cell summaries of log1p), matched
  rank-AUC with exact moments and a within-cell VIF. Under NB, notebook 65 found that joint test uncalibrated (A4·JNT:
  15; 1, 1 or 62; 9, 9; 0.079–0.105).
- **New: the segment twin of the primary.** Everything is the positional F-test except the basis, which is the package
  S1/S2/S3 labels instead of the spline. So the row isolates continuous against discrete.
  - Rows: the 71 units; NB GLM on raw counts, offset log(mt-free library), 71's weights n / (2·J·n_j).
  - Designs: `segment_designs` (segment intercepts as the shared shape, within-group specimen contrasts; in a
    relabeling, species × segment as nuisance); `level` adds the group; `steps` adds group × segment; `spec` adds
    within-group specimen contrasts × segment.
  - Numerator `T_step` = LR(level → steps), df₁ 2. Denominator `LR_spec` = LR(steps → spec), df₂ 4 (condition) or 2
    (relabeling), asserted from the design ranks.
  - α: 71's condition-split α, reused for every partition (the NB note's rule for the step screen).
  - `moderated_f` (squeezeVar with log mean count, robust FALSE), z = Φ⁻¹(1 − p), gene q by BH; cameraPR on z (0.01);
    call = Up and FDR ≤ 0.05.
- Per-segment group differences: Δ_S1, Δ_S2 and Δ_S3, the case − reference natural-log rate ratio in each segment at
  equal specimen weights (`standard_gene_model`'s Δ convention), with HC3 SEs, from the steps fit of each partition.
- Sensitivity: `T_step` through 3.2's joint test (the old anchor's set test, on NB).
- Notebook 75 owns this model (O2). It lives in `pseudospace/discrete_comparators.py` (section 12) and its results are
  saved under stable names (section 8). Notebook 74 runs after 75: it imports the module and reads the tables for its
  beyond-steps comparisons, and it does not refit the model.

**3.4 · Conventional pseudobulk DE** (notebook 31's code, unchanged; D2 keeps its model).
- Pseudobulks: raw counts summed per specimen × package segment over the 71 units; whole PT = S1 + S2 + S3. mt genes
  are not in the matrix, so DESeq2's size factors are mt-free.
- Genes: 71's universe (O3; 31's count filter is not applied). In a scope, a gene whose four pseudobulks total zero has
  no Wald stat and is untested there (counted). A gene whose padj PyDESeq2's independent filtering sets to NaN is not
  listed.
- `conventional_pt.deseq2_contrast` (PyDESeq2, Wald, `refit_cooks=False`), per scope (whole PT; S1; S2; S3).
  Designs `~condition` (2 residual df); relabeling `~condition + group`, contrast group (1 residual df).
- **Signed GSEA** (31's `signed_gsea`): `gp.prerank` on the Wald stat, weight 1, `method='multilevel'`, seed 12,
  10–300. BH per partition over sets (whole PT) or set × segment (S1/S2/S3); call = q ≤ 0.05 in either sign.
- **ORA** (56's rule): lists = padj ≤ 0.05 and |log2FC| ≥ 1, one per direction (whole PT 2 lists; S1/S2/S3 6);
  one-sided hypergeometric against the genes with a Wald stat in that scope; BH over set × list (× segment) within
  the partition; call = any test q ≤ 0.05. Implemented in the new module (below), not with `over_representation`,
  whose BH runs only over sets holding a listed gene.

**3.5 · Gene-first curve modules** (56(c) and 61's B2 → NB Δ(x) from 71's cache, 101-point grid).
- Curves: each gene's Δ(x), centred and scaled. k-means (20 starts, seed 56), k from 3–8 by the mean ARI over 20
  resamples of 80% of the genes (56's rule). ORA per module against the partition's tested genes, BH over module × set.
- 56(c) twin: the top 10% of genes by that partition's z (the primary gene statistic, replacing Gaussian `T_spatial`).
- B2 twin: the genes with q ≤ 0.05 on that partition's z (71: 3,596 hvm, 4,358 aki; relabelings 0, 0 and 16, 69, so
  this row is near-trivially specific; reported as such).

**3.6 · Rows copied from 71** (no recomputation; both datasets unless noted): the positional F-test (primary);
clustering, default and single-split nulls; S1 shape GSEA, level and total cameraPR and GSEA; S2 ssGSEA level and
shape; S3 level list + ORA (hvm only, as in 71).

### 4 · Probe exposure (hvm, every row)

Equal-probe genes: `probe_balance == 'balanced'` in `probe_counts.csv` (the same set as 56's equal counts; 9,077 of
71's 11,055 hvm genes). Each row is rerun on them without refitting any
gene model: set tests on the restricted statistic; DE lists with BH recomputed within the restricted genes; the
clustering ORA on 71's cached groups and clusters with the universe restricted (as 71's proteasome sensitivity); sets
re-filtered to 10–300. Reported: calls, calls kept, and **exposure = 1 − |calls ∩ equal-probe calls| / calls**. A row
is probe-exposed when exposure ≥ 0.5 (56's rule).

### 5 · Overlap and what each row calls

- For the primary against every other row (condition split): both / primary only / other only / neither and the Jaccard,
  over sets and over 71's redundancy clusters.
- **Offset lean** of each call, from 71's S1 tables: a set leans to the offset when its level-arm |z| (cameraPR,
  two-sided) exceeds its shape-arm z. Reported per quadrant (61 found 66 of 92 ORA-only calls offset-leaning, 0 of the
  continuous calls).

### 6 · Programs and the Figure 4e heatmap (pathway-first against gene-first)

Notebook 75 owns the programs and the Figure 4e heatmap (O4); later notebooks read `programs_{ds}.csv`.
- **Programs**, per dataset: `pathway_decoys.overlap_programs` on 71's calls of that dataset sorted by `pathway_id`
  (61's order-invariant rule; average linkage on 1 − overlap coefficient, cut 0.5). Contributing genes: the top 10% by
  the condition z (replacing Gaussian `T_spatial`); all members when fewer than 3.
- **Heatmap**, one per dataset: the program with the most calls (ties: lowest best FDR). Rows: its contributing genes, at
  most 40 by z, ordered by reference then case peak bin (hvm mouse then human; aki control then AKI). Columns: specimen
  × 12 equal bins of the 71 position, mean log1p CP10k (mt-free library), z-scored within group. Side bars: z,
  whole-PT DESeq2 log2FC (3.4) and, in hvm only, notebook 73's symmetric class. The hvm heatmap is Figure 4e; the aki
  one is supplementary.
- **Gene-first panel:** the B2-twin module (3.5) with the most enriched sets, same layout; beside it the 71 clustering
  cluster with the most calls.

### 7 · The level positive control on the package (68, hvm)

Notebook 68's protocol unchanged except the inputs:
- **L_g** = `grid_centred_level` on 71's condition-split Δ(x) (reference: converged, not separated), in log2;
  non-converged genes dropped. Guard: equals 71's `s3_level_gene_table.csv` L on shared genes. Pair spread, same-species
  bands and P_g as 68, with package units, labels and the mt-free library. The centring at each x removes the constant
  part of the offset change.
- **A_g** refitted by 68's code (Lake 7 donors against 6 Census adult males; limma-trend robust; centred by the median)
  on the analysed genes: 71's universe present in both atlases. Sensitivities as 68 (12 males, Lake men, Census human).
- **Criteria** (i)–(iv), the 80% thresholds, the housekeeping band, the equal-probe and pairing sensitivities, the
  probe-leak Spearman, and the pathway means over 71's sets: 68's definitions.
- **Families:** sums from the stage-01 native unit matrices for the 71 units (id `<sample>_unit_<k>`), with 68's
  members (notebook 62) and the gene-model library.
- 68's section 12 checks (controls against the atlas sign; pre-named family probes) are re-reported as descriptive.

### 8 · Outputs

All under `results/pt_package_pathway_methods/tables/`, with `{ds}` = `human_vs_mouse` or `aki_vs_control` (71's
names) and column blocks per partition named `condition`, `relabeling 1`, `relabeling 2` (71's layout).

**Segment F-test (stable names; notebook 74 reads these):**
- `gene_statistics_segment_{ds}.csv`: one row per gene of 71's universe; per partition `T_step`, `LR_spec`, `df2`, `s2`,
  `s2_prior`, `d0`, `s2_post`, `F`, `df_total`, `p`, `z`, `q`, `converged`, `separated`, `delta_S1`, `delta_S2`,
  `delta_S3` and their `se_S1`–`se_S3`;
- `pathway_tests_segment_{ds}.csv`: one row per 71 set; per partition cameraPR's `n`, `direction`, `p`, `fdr`,
  `p_one_sided` and `call`;
- `pathway_tests_segment_joint_{ds}.csv`: the joint-test sensitivity, same layout plus `effect`, `z_joint`;
- stage-cache names `segment_ftest_{ds}_{label}`, with 71's partition labels.

**Method table and rows:**
- `method_table.csv`: one row per method × dataset, with family (positional, discrete, average-level, gene-first,
  score), source (71 or 75), gene model, calls (redundancy clusters), relabeled calls, NCDR, tested, calibration share,
  pass, equal-probe calls and exposure (hvm), Jaccard with the primary, the old number, and what the row returns.
- `pathway_tests_{row}_{ds}.csv` (all partitions), with `{row}` one of `raw_camerapr`, `raw_camerapr_ranks`,
  `raw_gsea_ranknormal`, `raw_gsea_raw`, `nb_joint`, `deseq2_wholept_gsea`, `deseq2_segments_gsea`,
  `deseq2_wholept_ora`, `deseq2_segments_ora`, `modules_top10`, `modules_significant`.
- `deseq2_genes_{ds}.csv` (scope × partition), `modules_{ds}.csv`, `diagnostics_raw_tspatial.csv`,
  `probe_exposure.csv` (hvm), `overlap_{ds}.csv`, `programs_{ds}.csv`, `program_heatmap_{ds}.csv` (the plotted
  matrix), `compare_old.csv`, `material_change.csv`, and section 7's tables under 68's names in `tables/level/`.

**Figures:** the method table as a dot plot (66's layout, every row, both datasets); Figure 4e and its aki twin; 68's
main figure.

### 9 · Comparison with the old numbers (descriptive)

`compare_old.csv` sets each row beside its old result: 66's rows (SCF13, v5), 31/56/61's conventional rows and
modules (SCF13), 45's step screen (Gaussian; hvm 50; 0, 0, aki 26; 12, 6 on v5) and 68's criteria. The
unit, label, coordinate and offset changes are 71's (D1–D8 there);
3.3's statistic is new, so its old row is a different test and is labelled so.

### 10 · Material changes (stated in the summary and the ledger)

- **M1.** cameraPR on rank-normal raw `T_spatial` passes the rule on either dataset (old: fails both parts on both).
  It would weaken the case for the specimen-shape denominator.
- **M2.** The NB joint test changes verdict on either dataset (old: passes both).
- **M3.** A signed average-level GSEA row (whole PT or S1/S2/S3) passes on hvm, or a DE + ORA row fails (old: fail,
  pass).
- **M4.** S1/S2/S3 DE + ORA is no longer probe-exposed (exposure < 0.5; old 0.60), or the primary becomes exposed
  (exposure ≥ 0.5).
- **M5.** The segment F-test fails the rule on hvm. R4's "the step screen is specific" rests on the old row.
- **M6.** 68's reading changes: (i) or (ii) changes status (old: fail, pass), (iii) moves by more than 0.10 (old 0.60),
  or the probe-leak Spearman falls below 0.2 (old 0.45).
- **M7.** The primary and any discrete or average-level row share more than 0.3 of their redundancy clusters
  (Jaccard; old over sets: 0.07–0.17, steps 0.75 excepted). The segment F-test is reported but not judged by M7.

Rows with no old number (most aki rows; the probe exposure of rows other than the ORA rows) have no material-change
rule. Everything not listed above is descriptive. Nothing is tuned after a fit; a later check is labelled post hoc.

### 11 · Compute

No 71 model is refitted. On 24 workers, cold:
- inputs about 2 min; 3.1–3.2 about 10 min (GSEA dominates: 12 runs);
- 3.3 about 8 min (three partitions per dataset: nested LR, the denominator and the steps fit for Δ_S; no
  dispersion fit);
- 3.4 about 10 min with scopes in parallel (24 DESeq2 fits of 40–70 s, 24 signed GSEA runs of 20–25 s);
- 3.5, sections 4–6 under 5 min (section 4 reruns 12 GSEA on restricted genes);
- section 7 about 3 min (atlas fits, family sums from four native matrices).

About 40–50 min in all; stage-cached per row.

### 12 · New code

`pseudospace/discrete_comparators.py` with `tests/test_discrete_comparators.py` (synthetic only):
- The segment F-test, the API notebook 74 imports:
  - `segment_step_designs(segment, group, specimen, *, nuisance=None)`: the base, level, steps and spec designs and
    71's weights (built on `specimen_null.segment_designs`);
  - `segment_step_lr(counts, library, segment, group, specimen, dispersion, *, nuisance=None, n_jobs=1)`: `T_step`,
    `LR_spec`, df₁, df₂, the flags (`nb_gam.nested_nb_lr`) and Δ_S1–Δ_S3 with HC3 SEs;
  - `segment_f_test(...)`: the same plus `pathway_pipelines.moderated_f`, returning the gene table of section 8.
  - Tests: nesting; df₁ 2 and df₂ 4 (condition) and 2 (relabeling) from the ranks; `T_step` ≈ 0 when the groups share
    their steps; known per-segment log ratios recovered; a relabeled null gives no excess of small p.
- `family_ora`: 3.4's hypergeometric ORA with BH over the whole family. Tests: against `scipy.stats.hypergeom` and a
  hand BH.
- `swap_rule`: section 2's pass rule. Tests: the 71 rows reproduce 71's verdicts.

### 13 · Coordinator decisions on this draft (2026-10-12)

- **O1.** "Our cohort sits at the external median" is dropped, because the external pools are Gaussian. Our cohort's
  step-screen row stays, as the NB segment F-test (3.3).
- **O2.** Notebook 75 owns the segment F-test, in `pseudospace/discrete_comparators.py` with synthetic tests. Notebook
  74 runs after 75, imports the module and reads the section 8 tables. 74 keeps the coordinate-based steps (on the
  spline, six bins).
- **O3.** Every row uses 71's gene universe and gene sets, DE comparators included.
- **O4.** Notebook 75 owns the programs and the Figure 4e heatmap.
- **O5.** Every row runs on both datasets wherever it is defined.
