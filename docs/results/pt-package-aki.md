# AKI along the PT on the package pseudospace (notebook 76)

Part of the switch to the package pseudospace; see `pt-package-pseudospace.md`.

## Notebook 76 protocol: AKI along the PT on the package pseudospace (Figure 6)

Frozen on 2026-10-12, before any fit on private data. The coordinator's decisions of 2026-10-12 on O1–O6 are
built into the text (section 13). The only private data read so far came from loading checks (sections 2 and 5.4),
which made no contrast and no fit.

- Notebook: `analysis/notebooks/76_pt_package_aki.ipynb`
- Logic version: `76.package_aki.1`
- Results: `results/pt_package_aki/` (`tables/`, `figures/`, `stage_cache/`, `run_record.json`)
- New module: `pseudospace/aki_along_pt.py`, with `tests/test_aki_along_pt.py`. It holds the segment-unit map, the
  program profiles and their summaries, and the kNN placement.

**Why.** The outline's AKI section (user outline, Results 5) claims "progressive, segment-dependent decoupling of
transport functions and cell-cycle activation" along the PT in ischaemic AKI. Figure 6 (`figure_plan_methods.md`) was
blocked on an AKI coordinate; the package now provides one. Notebook 71 tested AKI on a universe restricted to the human
vs mouse genes. That universe lacks Havcr1, Lcn2, Mki67 and Top2a, so Figure 6 needs its own universe.

**Priors.**
- Ischaemia–reperfusion (IRI) injures the outer-stripe S3 most. Models differ in where tubular injury falls
  ([PMID 19759527](https://pubmed.ncbi.nlm.nih.gov/19759527/)).
- Two days after IRI, about half of outer-medullary epithelial cells are Ki67+. Repair comes from surviving
  differentiated cells ([PMID 18371453](https://pubmed.ncbi.nlm.nih.gov/18371453/)).
- After IRI, PT injury starts at the cortico-medullary boundary and spreads to the cortex
  ([PMID 34183416](https://pubmed.ncbi.nlm.nih.gov/34183416/)). Krt20 tracks progressive PT injury
  ([PMID 28931758](https://pubmed.ncbi.nlm.nih.gov/28931758/)).
- Visium IRI maps show region-specific loss of differentiation markers
  ([PMID 34853151](https://pubmed.ncbi.nlm.nih.gov/34853151/)). They resolve regions, not position along the PT.
- G2/M transcripts also mark arrested PT cells, not only dividing ones
  ([PMID 20436483](https://pubmed.ncbi.nlm.nih.gov/20436483/)).

**Cohort.** The cohort is two control and two IRI male mice. The user states that the two IR mice were injured at the
same severity, and that the segmentation training sections are not these sections. Specimens are the replicates, and
all results are descriptive. Two points of wording:
- Never write that the IR mice "differ in severity" (that was notebook 45's and the figure plan's wording).
- Differences between the two IR mice are called "positional differences between the IR mice".

### What Figure 6 must support

| # | Claim (panel) | Needs |
|---|---|---|
| A1 | IR sections hold fewer PT-labelled units in the late PT (6a) | all PT units on one axis; label-loss accounting (section 8) |
| A2 | Positive control: injury markers rise most in S3; the positional screen finds injury Hallmarks (supplement) | mouse-panel universe; notebook 45's rule |
| A3 | Transport loss along the PT, larger toward the late PT (6b) | a fixed transporter program; replication across pairs; relabelings |
| A4 | The cell-cycle program rises in IR, located along the PT (6c; neutral wording, 5.4) | a fixed cell-cycle program; the same rules |
| A5 | Decoupling: activation and loss differ in position (6d) | the precise definition in 5.3 |
| A6 | Specificity under the two relabelings (6e) | notebook 71's pass rule, plus per-claim checks |
| A7 | Label limitation: the results with c19 added (supplement) | placement of off-cohort units (section 8) |

### 1 · Inputs and guards

**Inputs.**
- Notebook 70's AKI export, `results/package_pseudospace/aki/`:
  - `pt_units.h5ad`, read through `package_inputs.read_pt_units` and `aki_inputs`;
  - `all_units_positions.csv`;
  - `tables/pt_stretch.csv`, which gives the pooled PT cut c = 0.285 (raw).
- The package run folder, `results/package_pseudospace/runs/aki/`, read only:
  - `03_integration/pass1_labelled.h5ad`: `leiden_coarse` and `X_harmony` (50 dimensions);
  - `02_qc/units_qc.h5ad`: `layers['counts']`;
  - `01_units/*_units.h5ad`: obs `seg_class` and centroids.
- Gene sets:
  - the ortholog map `results/minimal_pt_scfates/ortholog_map_used.csv`;
  - the three pathway JSONs;
  - `data/external/hgnc/hgnc_group_{752,417}.tsv`, with sha256 asserted as in notebook 74.
- Read only, for comparison: notebook 71's tables and its stage cache `gene_model_aki_vs_control_*`. Notebook 72's
  AKI tables are quoted, never used.

**Guards** (not cached).
- G1: the export contract and provenance, as in notebook 71: package 59da04b and pass-1 fingerprint `bb4a80c1b484`.
- G2: in every specimen, the segment medians increase S1 < S2 < S3.
- G3: the 13 `mt-` genes leave the library and the universe. The library equals notebook 71's exactly.
- G4: df₂ is 12 for the condition and 6 for each relabeling.
- G5 (reuse, decision D2): on notebook 71's 10,198 genes, the counts, library and positions are identical to notebook
  71's (digest). For genes converged in both runs, α, T_spatial, LR_spec and Δ equal notebook 71's cached payloads
  (rtol 1e-6). The refit costs about 4 min. Concatenating notebook 71's payloads instead would be wrong, because the
  median used to fill α for non-converged genes depends on the universe.
- G6: c19 is pass-1 `leiden_coarse` 19, with 80 units labelled cTAL, of which 75 are IR. c15 is Unassigned, with 2,425
  units and no pseudospace. Any difference stops the run.
- G7: the HGNC files are read from `data/external/hgnc/` and must match notebook 74's sha256. If they are missing or
  differ, the T and T_PT sections and the claims that use them (A3, A5) stop and are reported as not run. Every other
  section proceeds.

### 2 · (a) The mouse-panel AKI universe

- **Pool.** Export genes that are mouse symbols of the one-to-one ortholog map: 15,971 of 19,059. Keeping to the pool
  keeps notebook 71's offset, so its fits are reused (G5). The pathways are human libraries mapped through the same map.
  The 3,088 genes without a one-to-one ortholog are not tested; the number of them that pass detection is reported.
- **Universe U76.** Pool genes with equal-specimen detection ≥ 2% over the four AKI mice, minus the 13 mt genes. These
  are AKI's own detection rule and notebook 71's threshold.
  - Built with two `aki_inputs` calls: first with the pool as candidates (mt kept), then with the detected genes and
    `exclude_genes = MT`.
  - The second call builds notebook 37's sets directly on U76: members mapped through the ortholog map, 10–300 tested
    members. Notebook 71's two-step restriction, through the human vs mouse universe, is not used.
- **Loading check** (2026-10-12, no fit). It found:
  - 10,701 genes, a superset of notebook 71's 10,198 (503 new);
  - 1,495 sets: 1,467 of them also tested in notebook 71, 423 of those with the same members;
  - all 18 of notebook 45's injury markers (notebook 71: 14), plus Mki67, Top2a and Pcna;
  - Hallmark G2-M with 176 members and E2F with 186 (notebook 71: 138 and 146).
- **Relation to notebook 71 (decision O3).** Notebook 71's AKI calls stay the record for the method comparison
  (notebook 75's AKI rows), which use the same universe rule as human vs mouse. Notebook 76's calls serve Figure 6
  only. Every notebook 76 table and figure that holds calls is labelled "Figure 6 universe (notebook 76); not the
  method-comparison record". Reported:
  - the Spearman of z on shared genes;
  - calls kept, lost and gained among sets tested in both runs;
  - the Jaccard over redundancy clusters.

### 3 · Gene model and screen (notebook 71's sections 2–3, unchanged, on U76)

- **Gene model.** `standard_gene_model`: an NB GLM on raw counts, offset log(mt-free library), 6-df cubic B-spline, the
  101-point `common_grid`. α comes from the condition's M_full and is reused for the relabelings.
- **Positional F-test.** `specimen_shape_lr`, then `moderated_f` (covariate log mean count, robust FALSE). Gene q is BH.
- **Set test.** `camera_pr` on z (inter-gene correlation 0.01). A call is Up with FDR ≤ 0.05.
- **Partitions** (`partition_groups`): the condition split; relabeling 1, {Ctrl1A2, IR2A2}, which is aligned with the
  A2/A4 factor; relabeling 2, {Ctrl1A4, IR2A2}.
- **Position.** Notebook 71's support min–max rule.
- **Redundancy clusters.** The package's `cluster_gene_sets`.
- **Pass rule (notebook 71).** NCDR < 0.25; each relabeling calls < 5% of tested sets; the calibration share (pooled
  relabeled one-sided p ≤ 0.05) is ≤ 0.075.
- **Gene-level calibration** is reported: the pooled share of relabeled gene p ≤ 0.05.
- The clustering comparator is not run.

### 4 · Programs, fixed now

- **T, the transporters.** Notebook 74's set: the HGNC "Solute carriers" (752) and "ATP binding cassette transporters"
  (417) branches, mapped by `human_symbol` into U76 ([PMID 36243972](https://pubmed.ncbi.nlm.nih.gov/36243972/)).
- **T_PT (primary transport program; decision O2).** T without the SLC25 family (approved symbol prefix `SLC25`; these
  are mitochondrial carriers), restricted to PT-enriched genes in the control mice.
  - A gene is PT-enriched when log₂(rate_PT / rate_nonPT) ≥ 1 in both Ctrl1A2 and Ctrl1A4.
  - rate = Σcounts / Σ(mt-free library) over that specimen's pass-2 units: PT-S1/S2/S3 against all its other pass-2
    units.
  - **The selection is independent of the IR contrast.** It uses control mice only, and only their PT against their own
    non-PT units; no IR unit and no IR − control difference enters it. The one shared quantity is the control PT rate.
    Its regression-to-the-mean bias is negligible, because each control holds about 2,700 PT units.
  - If fewer than 20 genes qualify, T without SLC25 becomes primary.
  - Notebook 74's full T is the sensitivity.
- **C, the cell-cycle program.** The union of MSigDB Hallmark 2020 "G2-M Checkpoint" and "E2F Targets"
  ([PMID 26771021](https://pubmed.ncbi.nlm.nih.gov/26771021/)), mapped into U76: 297 genes.
  - G2-M and E2F are also reported separately.
  - Its interpretation follows the two readings of 5.4.
  - Until the IRI time point is known, all figure text and claims call it the "cell-cycle program", never
    "proliferation" or "arrest" (decision O6).
- **R, injury and repair.** Notebook 45's 18 markers (Havcr1, Lcn2, Krt20, Krt8, Krt18, Spp1, Cd44, Vcam1, Sox9,
  Cdkn1a, Clu, Mt1, Mt2, Atf3, Hmox1, Timp1, Sprr1a, Tnfrsf12a) and its 13 injury Hallmarks, unchanged.
- **Display genes only:** Mki67, Top2a, Havcr1, Lcn2, and Figure 2a's transporters (Slc34a1, Slc5a2, Slc5a12, Slc22a6,
  Slc7a13, Slc22a7).

### 5 · Definitions

**5.1 Segment units.**
- t₁₂ and t₂₃ are the cuts on the analysed position that maximise the balanced accuracy of S1 against S2 and of S2
  against S3. Only control units inside the support are used, so injury cannot move the scale.
- The piecewise-linear map is [0, t₁₂] → [0, 1], [t₁₂, t₂₃] → [1, 2] and [t₂₃, 1] → [2, 3], extended linearly outside
  [0, 1].

**5.2 Program profiles.**
- **Δ_g(x).** The natural-log rate ratio IR/control between the equal-specimen group curves (`GeneModel.delta`).
- **Pair profiles.** Δ_g^{ic}(x) = specimen curve i − specimen curve c, for each of the 4 IR × control pairs.
- **Relabeled profiles.** The relabelings' `delta`.
- **Centring.** m(x) is the median of Δ_g(x) over converged, non-separated genes; each pair and each relabeling uses its
  own median. This is notebook 68's centring, the positional analogue of median-of-ratios normalization: library
  scaling is biased when highly expressed genes change ([PMID 20196867](https://pubmed.ncbi.nlm.nih.gov/20196867/)).
  m(x) is reported, and the uncentred profiles are a sensitivity.
- **Profile.** P(x) = median over member genes (converged, not separated) of [Δ_g(x) − m(x)] / ln 2, in log₂.
- **Loss and activation.** L(x) = −T_PT(x) and A(x) = C(x).
- **Summaries.**
  - level P̄: the mean over the grid;
  - drift d: the mean over the last 21 grid points minus the first 21 (notebook 71's windows);
  - centre of mass μ = Σ s(x)·P⁺(x) / Σ P⁺(x), in segment units s, where P⁺ = max(P, 0);
  - peak: the argmax;
  - within-segment change: P at a segment's end minus its start.
- **Discrete view.** Program medians of the centred log₂ IR/control ratio of specimen pseudobulks (Σcounts/Σlibrary)
  for the whole PT and for S1, S2 and S3. It is drawn beside every profile.

**5.3 Decoupling.**
- **Primary: positional decoupling.** Δμ = μ_A − μ_L, in segment units. Coupled responses (activation proportional to
  loss at every position) give Δμ ≈ 0. Also reported: R² of A(x) on L(x) through the origin across the grid.
- **Secondary, descriptive: unit coupling.** Unit scores are s_P = log₂(10⁴ · Σ_{g∈P} counts / library + 1). Within
  each of 10 control-decile bins and each specimen, report the Spearman ρ(s_{T_PT}, s_C), and the median over bins,
  IR against control. Units are multicellular cross-sections, so this is not a cell-level statement.

**5.4 Two readings of the cell-cycle program (decision O6).** The IRI time point is unknown to us. A rise in C fits
two readings:
- **Division:** repair by surviving epithelial cells, with Ki67 peaking around day 2
  ([PMID 18371453](https://pubmed.ncbi.nlm.nih.gov/18371453/)).
- **Arrest:** PT cells held in G2/M ([PMID 20436483](https://pubmed.ncbi.nlm.nih.gov/20436483/)).

**Readouts, computed only where A4-activation holds:**
- **B1 · phase balance.** G2M′(x) − E2F′(x) at the condition's C peak, in centred log₂. G2M′ is the 111 G2-M-only genes
  and E2F′ the 121 E2F-only genes, without Cdkn1a; the 65 shared genes are left out of both.
  - Balanced is |·| < 0.25. G2/M-dominant is ≥ 0.25 in all 4 pairs.
- **B2 · mitotic genes.** Mki67, Top2a, Cdk1, Ccnb2, Plk1, Aurkb and Cdc20. Their gene Δ is ≥ 0.25 log₂ in all pairs,
  peaking within ±0.5 segment of the E2F′ peak.
- **B3 · arrest genes.** Cdkn1a (p21) and Gadd45a/b/g. Their gene Δ is ≥ 0.25 log₂ in all pairs, peaking within
  ±0.5 segment of the G2M′ peak.
- **B4 · unit coupling** (descriptive). In IR units, per decile bin, ρ(s_G2M′, s_Cdkn1a) against ρ(s_G2M′, s_E2F′).

Ccnb1 is not on the mouse probe panel; Ccnb2 stands in for it. The other genes are all in U76 (loading check).

**Verdict:**
- **Division-compatible:** B1 balanced, B2 holds and B3 fails.
- **Arrest-compatible:** B1 G2/M-dominant and B3 holds.
- **Undetermined:** anything else.

The verdict is reported as compatibility, not mechanism. Transcripts cannot separate a slow cycle from arrest, and
Cdkn1a also rises with injury itself (it is one of the R markers). The figure keeps the neutral "cell-cycle program"
in every case. Once the user gives the time point, it is used only as a prior, never to re-select the verdict.

### 6 · Analyses

**(b) Positive control: notebook 45's D6 on the package axis.**
- **Markers.** For each R marker, the IR − control difference of equal-specimen segment means of Y. Report the share
  that rises most in S3 (notebook 45: 13 of 14, on the mouse-only PT DPT, v5). The continuous version is the share
  whose argmax of Δ_g(x) lies at s ≥ 2.
- **Hallmarks.** The injury Hallmarks called by the screen, and the Mann–Whitney test that they rank above the other
  sets (one-sided).
- **Location.** The share of called sets whose members' median |Δ| peaks in S3.
- **Detection passes** if ≥ 2/3 of the markers rise most in S3, and either ≥ 3 Hallmarks are called or the
  Mann–Whitney p ≤ 0.01.
- **Specificity** is reported under notebook 45's rule (mean relabeled calls ≤ 5% of condition calls) and under the
  R4 rule. Before: notebook 71 had 103 calls and 13 and 29 relabeled, so it passes R4 and fails notebook 45's rule;
  notebook 45 had an NCDR of 0.17.

**(c) Density.**
- **Units.** All 15,013 export PT units, untrimmed, on the support-scaled axis and in segment units.
- **Per specimen.**
  - the shares in 10 bins cut at the equal-specimen deciles of the control units; the tail beyond the control 99th
    percentile;
  - the **late share**, the last two bins: the primary statistic for A1;
  - a KDE, for display only.
- **Context.**
  - PT units per section, and per 1,000 units on the curve;
  - segmentation `tubule_proximal` polygons with bins per section, and the share that stage 03 labels PT;
  - the deep-S3 count, by notebook 45's rule. Depth is the mean distance to the 3 nearest glomerulus polygons of the
    section (stage 01 `seg_class`). A unit is deep S3 when it is deeper than the section's S1 95th percentile
    (notebook 45: 248 and 91 against about 830);
  - in the control S3 units, the Spearman of position against depth. This is the anatomical anchor; x/y is used for
    validation only.
- **Accounting.** The IR late-PT deficit (the control mean late share × IR PT units, minus the observed count) is shown
  beside the IR c19 and c15 units that section 8 places on the late PT.

**(d) Transport, cell cycle and decoupling.**
- **Profiles.** L, A (and G2-M and E2F separately) and R: condition, 4 pairs and 2 relabelings, with the discrete view.
- **Heatmap.** T_PT members, centred log₂ Δ, clipped at ±2, rows ordered by the position of maximal loss.
- **Decoupling.** Δμ, the proportional-coupling R², and the unit ρ, as defined in 5.3.

**Cohort facts** (inventory item 35; descriptive).
- Whole-PT and per-segment pseudobulk log₂ CPM (mt-free library) per specimen, for the R markers, Pcna, Mki67, Top2a
  and the display transporters.
- Over the U76 genes with log₂CPM ≥ 1 in every specimen:
  - the SD of IR2A2 − IR2A4 and of Ctrl1A2 − Ctrl1A4;
  - their Spearman correlation (the A2/A4 factor).

  Both are given for the whole PT and per segment. They describe the data; they do not test equal severity.

### 7 · (e) Specificity, and how Figure 6 handles the failed calibration

- **Expected result.** Notebook 71's AKI screen passes the relabeling rule and fails calibration (0.120). Notebook 76
  is expected to do the same.
- **What Figure 6 shows.**
  - It shows the screen's pass row with "fails calibration" printed on it.
  - It presents no AKI call list as a discovery list. Only the pre-specified programs carry claims.
- **Set p of a program.** Empirically calibrated against the relabelings:
  p_emp = (1 + #{relabeled one-sided cameraPR p ≤ p_S}) / (1 + N), pooled over all U76 sets and both relabelings.
  - A set claim needs p_emp ≤ 0.025 (Bonferroni over T_PT and C) and its own relabeled p > 0.05 in both relabelings.
  - The ceiling: this assumes the relabeled sets are exchangeable nulls. Because they carry real between-mouse
    differences, the test errs conservative.
- **Profile specificity.** For a statistic S in {P̄, d}: |S_relabel| < 0.25·|S_condition| in both relabelings. This is
  the NCDR threshold.
- **Level claims** rest on pair replication and profile specificity only. No calibrated level test exists: notebook
  71's level arms fail.
- **Gene-level claims** (named genes): q ≤ 0.05, but only if the gene-level calibration share is ≤ 0.075. Otherwise
  only the profile is described.

### 8 · (f) The label limitation: c19 and c15

**Loading check, already known.** c19's package positions lie mostly beyond the PT stretch: 7 of 80 at raw ≤ c,
median 0.387. "Added as positioned" is therefore nearly empty after the common-support trim. c15 has no pseudospace,
because stage 03 dropped it before the curve was fitted.

**Placement rule (kNN).**
- Neighbours: k = 30 (the package's neighbour count) in pass-1 `X_harmony`, Euclidean.
- Reference: units with a pseudospace.
- Each placed unit gets:
  - the PT share of its 30 curve neighbours;
  - their median raw pseudospace;
  - its `seg_class`.
- **"On the PT stretch":** a PT share ≥ 0.5 and a median ≤ c. Counts are given per specimen and per placed segment
  (the majority label of its PT neighbours), for c19 and for c15.

**Sensitivity runs** (descriptive, not replacements). Cohort = PT ∪ c19, in two arms:
- (i) at the package positions;
- (ii) at the median position of the 30 nearest PT units, with the majority PT label.

Each arm reruns section 3, the A1–A5 statistics and the positive control. Each change is tabulated against the
primary. Decision O4 approved arm (ii) as the informative arm, because arm (i) is nearly empty.

**c15** gets only the count on the PT stretch (per specimen and placed segment, with its `seg_class` mix). There is no
c15 analysis arm, because the cluster mixes stromal and injured units (decision O4).

### 9 · Other sensitivities (descriptive)

- The package DPT (`dpt_pseudotime`, decision D3): full model refit, with the A2–A5 statistics.
- The full T in place of T_PT.
- Uncentred profiles.

### 10 · Claims, and what makes each not claimable

"All pairs" means all 4 IR × control pairs agree in sign. "Specific" means section 7's profile rule. A claim must also
hold, in direction, in the DPT sensitivity, or the text marks it "package coordinate only".

| # | Claimable when | Not claimable when (then reported as such) |
|---|---|---|
| A1 | both IR late shares are below both control late shares, and the IR mean is ≤ 0.75 × the control mean | either fails. Even when it holds, the wording is "fewer PT-labelled units". It becomes "consistent with relabelling" if the late-PT c19 + c15 units cover the deficit |
| A2 | detection passes (6b) | detection fails. Specificity is reported, not required (as in notebook 45) |
| A3 loss | L̄ > 0 in all pairs; L̄ ≥ 0.25 log₂; specific | any condition fails |
| A3 progressive | d_L > 0 in all pairs; d_L ≥ 0.25 log₂; specific; T_PT set p_emp ≤ 0.025 | any condition fails. Then: "loss without a replicated positional trend" |
| A4 activation | Ā > 0 in all pairs; Ā ≥ 0.25 log₂; specific. Worded "cell-cycle program" | any condition fails |
| A4 reading | the 5.4 verdict, stated as compatibility | it is "undetermined", or A4-activation fails |
| A4 location | each pair's peak within ±0.5 segment of the condition's peak; C set p_emp ≤ 0.025 | either fails. If only specificity fails: "the position of activation differs between the IR mice" |
| A5 | A3-loss and A4-activation both hold; Δμ has the same sign in all pairs; \|Δμ\| ≥ 0.25 segment | any condition fails. Then: "no positional decoupling resolved" |
| Within-segment (A3/A4) | the within-segment change has the same sign in all pairs, is ≥ 0.25 log₂ and is specific | any condition fails, or notebook 72's AKI within-segment order is below the human vs mouse thresholds. Then it is worded "along the axis within the S3 label" |
| A7 | — (descriptive) | a sign or claim verdict that flips in arm (ii) is stated beside the claim |

**Not worth a main figure.** If A3-loss and A4-activation both fail, Figure 6 becomes a supplementary figure (A1, A2),
and the outline's decoupling sentence is withdrawn.

### 11 · Outputs

- **Tables:**
  - `inputs.csv`, `universe.csv` and `gene_sets.json`;
  - `compare71.csv`;
  - `gene_statistics.csv` and `pathway_tests_{partition}.csv`;
  - `pass_rule.csv`;
  - `programs.csv` (members, T_PT enrichment);
  - `segment_units.csv`;
  - `profiles.npz` (condition, pairs, relabelings, m(x)) and `profile_summaries.csv`;
  - `discrete_view.csv`;
  - `positive_control.csv`;
  - `density.csv` and `density_context.csv`;
  - `decoupling.csv` and `unit_coupling.csv`;
  - `cell_cycle_readings.csv` (B1–B4 and the verdict);
  - `cohort_facts.csv`;
  - `placement_c19_c15.csv` and `sensitivity_c19.csv`;
  - `sensitivity_dpt.csv`;
  - `claims.csv`, with each verdict and the rule it used;
  - `run_record.json`.
- **Figure 6, panels a–e,** as in the figure plan:
  - positions in segment units, with the control transitions marked;
  - IR in one colour family and control in another, with pair lines thin.

### 12 · Compute

The gene model on 10,701 genes takes about 4 min (notebook 71: 80 s for the condition split, 35 s per relabeling and
about 15 s per denominator, 24 jobs).

| Item | Time |
|---|---|
| Primary run | about 5 min |
| Two c19 arms | about 10 min |
| DPT refit | about 4 min |
| Pseudobulks, kNN placement (36k × 50), depth, profiles | about 3 min |

A cold run takes about 25 min; a cached rerun takes under 2 min. Memory is below 16 GB. It needs only notebooks 70 and
71.

### 13 · Coordinator decisions (2026-10-12)

- **O1. Decoupling.** Δμ, the centre-of-mass shift in segment units, is the primary statistic. Unit-level coupling is
  descriptive (5.3).
- **O2. Transport program.** T_PT is primary; notebook 74's full T is the sensitivity.
  - T_PT's PT-enrichment uses the control mice only, so it is independent of the IR contrast (section 4).
  - The HGNC files are read from `data/external/hgnc/` with notebook 74's sha256. If they are absent at run time, the
    sections that need them stop and are reported (G7).
- **O3. Record.** Notebook 71's AKI calls stay the record for the method comparison. Notebook 76's calls serve
  Figure 6 only and are labelled so (section 2).
- **O4. c19 and c15.** The c19 neighbour-placement arm (ii) is approved. c15 gets only its count on the PT stretch; there
  is no c15 analysis arm, because it mixes stromal and injured units (section 8).
- **O5. Thresholds.** Approved as pre-specified: 0.25 for the effect sizes and the relabeled ratio, and 0.025 for the
  calibrated cut.
- **O6. Time point.** The IRI time point is unknown; the coordinator will ask the user. The cell-cycle program has two
  pre-specified readings, division and arrest, with the readouts that separate them (5.4). Figure text says
  "cell-cycle program" until the time point is known.
