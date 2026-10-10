# Gene model: negative binomial on raw counts

Decided by the user on 2026-10-09. This note freezes the specification before notebooks 37 onward
are re-run. It supersedes the Gaussian gene models described in SM12 of the supplementary information.

## Why the model changed

Until now every gene model was a weighted Gaussian regression spline on log1p(count / library × 10⁴).
A scratch pilot compared it with a negative-binomial (NB) GLM on raw counts. The pilot ran on the
private data, so this switch is **post hoc** and is recorded as a deviation in Supplementary Note 8.
Notebook 65 reproduces the pilot from tracked code.

**Pilot results (SCF13, species split and the two balanced relabelings):**

| Method | Gaussian on log1p | NB on raw counts |
|---|---|---|
| T_spatial + joint test (A6·JNT) | 60 calls; relabeled 2, 0; calibration 0.060 | 31; 0, 0; 0.067 |
| S1/S2/S3 steps + joint test (A4·JNT) | 50; 0, 0; 0.064 | 15; 1, 1; 0.079 |
| Continuous adds over steps / steps add over continuous | 13 / 3 | 16 / 0 |

**Why NB is better.** On the log1p scale, a constant species ratio is not a constant difference.
It bends where counts are low and where a gene is nearly absent in one species. The Gaussian
T_spatial therefore scores some offsets as shape:
- **Simulation.** Genes given a pure constant ratio reached the Gaussian top 10% at 10.1%, and
  27% in the largest-ratio tertile. Under NB the rate was 0%. The simulation draws counts from the
  NB model, so it favours NB by construction.
- **Real data.** 930 genes are detected in under 1% of one species' structures. 139 of them are in
  the Gaussian top 10% and 6 in the NB top 10%.
- **The relabeling control cannot catch this.** The artefact follows species (depth, probe panels,
  absent genes), and the relabelings mix species.

The voom split-plot model works on specimen × bin pseudobulk and keeps its own mean–variance
weights. Its significant genes are not enriched for nearly-absent genes: 5.1%, against a base rate
of 8.4%; the Gaussian top 10% was 12.5%. It is therefore unchanged.

## Frozen specification

**Model.** Each gene gets an NB GLM:
- log link and offset log(library). The library is the total over the 15,567 orthologs measured in
  all specimens.
- The same designs as before: `contrast_designs`, `segment_designs` and the species trajectory model.
- The same prior weights, wᵢ = n / (2 · J_s · n_j).
- Code: `pseudospace/nb_gam.py`, `GENE_MODEL = 'nb.1'`.

**Dispersion.**
- One α per gene per coordinate run, estimated by weighted maximum likelihood under the
  species-split continuous M_full on all structures. There is no Cox–Reid adjustment, since
  n = 12,272 and p ≤ 22.
- The same α is reused for the nested models, both relabelings, the within-species comparisons,
  every refit (stability, warped coordinate, cortical removal, power draws), the step screen and the
  curves.
- Notebook 45's external pools have no coordinate. Each pool gets its own α under the saturated
  donor × segment model.
- No shrinkage. A DESeq2-style check shrank each α toward a trend α = 0.024 + 0.011/mean, with a
  normal prior on log α. On SCF13 it gave the same 31 joint-test calls (relabeled 0, 0; calibration
  0.068), and Spearman ρ = 0.999 between gene ranks. The rule fixed beforehand (both pass and
  Jaccard ≥ 0.8) keeps the unshrunk estimate.

**Statistic.**
- The likelihood-ratio statistic, LR = max(D_reduced − D_full, 0), replaces the partial F.
- The quasi-likelihood F is computed and reported, and used for nothing.
- Gene statistics remain **ranks**. Their residual degrees of freedom count structures, so they
  are not calibrated gene tests.

**Failures.**
- A gene is non-converged when a fit fails R's glm criterion, the LR is not finite, or α does not
  converge. Its statistics are replaced by the median of the converged genes.
- A gene with zero counts in one species is **separated**. Its fit lies on the boundary and its
  limiting deviances are valid: T_spatial ≈ 0, because there is no positional information in a
  species without counts, and T_level is large. Separated genes are flagged and never filled.

**VIF.** The joint test's VIF comes from the Pearson residuals, (y − μ)/√(μ + αμ²), of the species
M_full, through the unchanged `residual_pathway_correlations`. The Gaussian-residual VIF is reported
in notebook 65 as a sensitivity.

**Species-difference curves.**
- Δ(s) is a natural-log rate ratio, human over mouse.
- Its standard error is an HC3-type sandwich with an expected-information bread; the model-based
  form is also reported.
- Display and descriptor thresholds keep their values, now read on the natural-log rate-ratio scale.
  They are SM17's |Δ| ≥ 0.25 and notebook 58's amplitude, range and prominence cut-offs. For
  moderately expressed genes the two scales coincide.

**Primary method.**
- T_spatial (NB LR) with the joint test is the primary pathway method **by design** (user,
  2026-10-09).
- Notebook 61 compares the methods side by side and no longer selects one by rule. Before this
  decision, its "most calls" step could have picked the split-plot ORA (35 calls) over the pilot
  NB joint test (31).
- Notebook 37's E2 fallback ("keep the earlier list as primary" when fewer than 20 pathways are
  robust) is removed, because the earlier list comes from the superseded model.

**Stop rule.** If the NB joint test fails the relabeling rule or the calibration criterion on SCF13
in notebook 37, the cascade stops and the result is reported. Nothing is tuned.

## Order of work

1. Notebook 65 (new) reproduces the Gaussian results while they are still live, then the NB
   comparison, the gene-level diagnostics and the simulation.
2. Notebook 37 on SCF13, DPT13 and SCF13-ED.
3. Notebooks 39, 43 and 47. Notebook 47 becomes the tracked writer of `coordinate_robustness.csv`,
   its funnel and `d4_coordinate_runs.csv`, which previously came from scratch scripts.
4. Notebook 50 (primary list), then 45 (many-draw and external nulls).
5. Notebooks 56 and 57, then 61 and 63.
6. Notebooks 48, 38, 58–60 and 41; then the draft, the SI and the figure plan.

The Gaussian results, the working notebooks with their outputs and the scratch provenance scripts
are archived in `results/_archive_gaussian_2026-10-09/`, which is outside Git.

## Notebook 65 results

`analysis/notebooks/65_pt_gene_model_comparison.ipynb` (logic version `65.gene_model.1`) writes to
`results/pt_gene_model_comparison/`. It recomputes the Gaussian side from code.
- On its first run every Gaussian statistic, VIF and call matched the live notebooks 37, 45 and 50.
  Notebook 45's pool was rebuilt and matched too.
- The Gaussian reference is saved in `reference/`, and later runs must reproduce it.
- The live-notebook guards switch off once notebook 37's gene table carries `gene_model = nb.1`, so
  notebook 37's NB re-run must write that column.

**Pathway screens on SCF13** (species calls; relabeled calls; share of relabeled p ≤ 0.05, limit 0.075):

| Screen | Gaussian | NB, Gaussian VIF | NB, NB Pearson VIF (frozen) |
|---|---|---|---|
| A6·JNT | 60; 2, 0; 0.060 | 31; 0, 0; 0.067 | 74; 13, 1; **0.085, fails calibration** |
| A4·JNT | 50; 0, 0; 0.064 | 15; 1, 1; 0.079, fails | 62; 9, 9; 0.105, fails |
| Continuous adds / steps add | 13 / 3 | 16 / 0 | 21 / 9 |

- **The pilot is reproduced exactly.** It used the Gaussian VIF, so the pilot table above
  describes the middle column.
- **The frozen primary fails criterion 2.** The NB LR with the Pearson-residual VIF, computed here
  for the first time, passes the relabeling rule (NCDR 0.09) but not the calibration limit.
- **Why.** The median VIF falls from 1.30 to 1.06, and the median within-pathway residual
  correlation from 0.010 to 0.002.
  - Gaussian log1p residuals share a depth component: within each specimen, a structure's mean
    residual correlates with its log library at r = 0.75–0.94.
  - The NB offset absorbs depth, giving r = −0.05 to −0.29.
  - Much of the Gaussian residual correlation was therefore shared depth. It happened to deflate the
    joint test enough to pass the relabelings.
- **Stop rule.** Notebook 37 would compute the same screen. Under the stop rule above, the cascade
  would therefore stop there with this specification. Nothing was tuned, and the choice is open.

**Gene level, as in the pilot:**
- T_spatial Spearman between the models is 0.36.
- The Gaussian-only top-10% genes have median |log ratio| 2.11 (NB-only: 0.83). 25% of them are
  nearly absent in one species (NB-only: 0.2%).
- Of the 930 nearly absent genes, 139 are in the Gaussian top 10% (12.5%) and 6 in the NB top 10%.
- The 4,588 split-plot genes are 5.1% nearly absent, against the 8.4% base rate.
- In the simulation, pure offsets reach the top-10% cut at 10.1% under Gaussian (27% in the
  largest-offset tertile) and at 0% under NB.

**Shrinkage.**
- With the Gaussian VIF it reproduces the pilot check: the same 31 calls (0, 0; 0.068; Jaccard 1.0),
  so the unshrunk α stays primary.
- With the NB Pearson VIF both versions fail calibration (0.085 and 0.087), so the rule cannot be met.

**What size the VIF should be.**
- Notebook 37 validated the Gaussian VIF by fitting E[z²] = a(1 + (m − 1)ρ) to the matched z of all
  1,513 pathways under the two relabelings, using `pathway_decoys.decoy_inflation`.
- The same fit applied to the NB statistic gives a = 1.04, ρ = 0.0096, and an implied median VIF of
  1.35. For the Gaussian statistic it gives 1.03, 0.0086 and 1.30.
- The NB statistic's own null correlation therefore matches the Gaussian-residual estimate (median
  ρ 0.0096). The Pearson-residual estimate (0.0018) is about 5 times too small.
- This uses only the relabelings, never the species split. The choice of VIF source is the user's
  (2026-10-09).

## Notebook 66: packaged set tests and clustering-first on the NB statistic

`analysis/notebooks/66_pt_pathway_pipelines_compared.ipynb` (logic version `66.pipelines.1`) writes to
`results/pt_pathway_pipelines_compared/`. It runs every pathway method on one shared NB gene model
(`pathway_pipelines.standard_gene_model`), for human vs mouse on SCF13 and for AKI vs control mouse
(notebook 45's inputs). Each row gives condition calls; relabeled calls; share of relabeled p ≤ 0.05:

| Method | Human vs mouse | AKI vs control |
|---|---|---|
| limma cameraPR on rank-normal T_spatial | 98; 46, 29; 0.121 (fails both) | 144; 87, 86; 0.170 (fails both) |
| gseapy prerank GSEA on rank-normal T_spatial | 146; 75, 80; 0.147 (fails both) | 203; 120, 146; 0.212 (fails both) |
| NB joint test, VIF fitted on the relabelings | 32; 3, 0; 0.067 | 12; 0, 0; 0.072 |
| Clustering-first ORA (Hallinan's package, single-split null) | 162 and 159; 0, 0 | 469 and 482; 2, 0 |

- **The packaged rank tests are not specific here.** Their relabeled z-scores are inflated by a
  constant factor (1.51 and 2.13), so no setting of `inter.gene.cor` repairs them. Ranking within
  expression × detection × coverage strata (post hoc) nearly passes human vs mouse (share 0.082) but
  not AKI (0.135).
- **The joint test's VIF is fitted on the relabelings that judge it.** Fitted on one and judged on
  the other, its share is 0.067 and 0.075; the AKI value sits at the limit.
- **Clustering-first is specific under a fair null, but returns no p-value per pathway.** It calls
  pathways whose genes share one direction, including constant offsets (oxidative phosphorylation)
  that T_spatial removes by design. It also finds AKI proliferation (G2-M, M phase), which T_spatial
  ranks low. Overlap with cameraPR is Jaccard 0.15 over redundancy clusters on both datasets.

## Notebook 67 protocol: T_spatial normalized by specimen shape variation

Frozen on 2026-10-09, before any fit. The user asked whether `T_spatial` should be normalized.
Notebook 66 suggested two reasons why packaged set tests on `T_spatial` are not specific:
- it grows with counts;
- it treats structures as independent, so shape differences between individual specimens count as
  group differences.

**Statistic.** The NB fits stay the same: the shared gene model, the species-split α, and the
weights and offset of `contrast_designs`.
- **Numerator.** `T_spatial`, unchanged: LR of M_full against M_full without group × basis, with
  df₁ = 6.
- **Denominator.** LR_spec = D(M_full) − D(M_spec), where M_spec adds a specimen-specific shape to
  M_full: specimen-within-group × basis. For the condition split df₂ = 12; for a relabeling, which
  keeps the condition shape as nuisance, df₂ = 6.
- s² = LR_spec / df₂ is moderated by `limma::squeezeVar`, with log mean count per structure as
  covariate. This is the primary. Sensitivities: no covariate, and `robust = TRUE`.
- The moderated F is (T_spatial / df₁) / s²_post, with df (6, df₂ + d₀). It converts to
  z = Φ⁻¹(1 − p).

**Pathway tests on z**, with gene sets as in notebook 66 (sizes 10–300):
- `limma::cameraPR`, with `inter.gene.cor` 0.01 and `use.ranks` FALSE. A call is Up with FDR ≤ 0.05.
  This is the primary test.
- gseapy prerank on z.
- ORA with `gseapy.enrich` on genes with BH ≤ 0.05; the background is all tested genes.
- Sensitivity: `cameraPR` on rank-normal scores of the F.

**Datasets.** Human vs mouse PT on SCF13 and AKI vs control mouse PT, each with the condition split
and both balanced relabelings, exactly as in notebook 66.

**Pass rule.** Notebook 66's rule, unchanged:
- **Relabeling.** NCDR < 0.25, and each relabeling calls < 5% of tested pathways.
- **Calibration.** The pooled share of relabeled one-sided p ≤ 0.05 is ≤ 0.075.

**Gene-level checks**, reported for information and not used for selection:
- under the relabelings, the share of genes with p ≤ 0.01 and p ≤ 0.05, and the median-based λ;
- the Spearman correlation of relabeled z with mean expression.

Genes with zero counts in any specimen are flagged. They stay in the primary and are dropped in a
sensitivity.

**Decision.** If the primary passes on both datasets, it becomes the user's candidate replacement
for the joint test as primary; the user decides. Otherwise the result is reported as is. Nothing is
tuned after the fits, and any later check is labelled post hoc.

## Notebook 67 results

`analysis/notebooks/67_pt_normalized_t_spatial.ipynb` (logic version `67.normalized.1`) writes to
`results/pt_normalized_t_spatial/`. It followed the protocol above; its section 10 lists how the
protocol was read. Cells give condition calls; relabeled calls; NCDR; calibration share:

| Test on the normalized statistic | Human vs mouse | AKI vs control |
|---|---|---|
| cameraPR on z (primary) | 70; 0, 0; 0; 0.044 · pass | 93; 17, 11; 0.15; 0.102 · fails calibration |
| GSEA prerank on z | 77; 0, 0; 0; 0.085 · fails calibration | 51; 34, 31; 0.64; 0.140 · fails both |
| ORA on genes with BH ≤ 0.05 (list-based) | 71; 0, 0; 0; 0.000 · pass | 100; 8, 3; 0.06; 0.029 · pass |

- **The decision rule is not met.** Every cameraPR sensitivity passes on human vs mouse and fails
  calibration on AKI (0.096–0.118).
- **The human pass is against a weak null.** Its replicates are two sections of one donor.
- **Gene level.** Under the relabelings the moderated F is close to calibrated: AKI λ 0.96–1.00;
  human vs mouse 1.35–1.50, against 1.8–2.7 for raw `T_spatial`. Its correlation with expression
  falls from about 0.2 to about 0.
- **Why AKI still fails (post hoc).** The relabeled pathway z are inflated by a constant factor,
  a = 1.49. The relabeled calls are coherent programs that differ between the two AKI mice: DNA
  replication, E2F, EMT and complement. A per-gene denominator cannot see that coherence.
- **The prior's share.** d₀/(df₂ + d₀) is 0.32–0.37 in human vs mouse and 0.55–0.74 in AKI.
- **ORA was a listed secondary test, not the primary.** In human vs mouse no gene reaches BH ≤ 0.05
  under either relabeling, so its relabeling pass there is trivial. In AKI 13 and 57 genes do.
- **Overlap.** The new calls are nearly a subset of notebook 66's raw cameraPR calls (Jaccard 0.70
  and 0.66 over redundancy clusters). They hold 28 of the joint test's 32 calls and 11 of its 12.

## Notebook 68 protocol: does the species level recover known differences?

Frozen on 2026-10-10, before any human-vs-mouse level statistic was computed on our Visium data or
the atlases. The literature controls were compiled by a separate agent that never saw our data
(sha256 below).

**The question.** This is a positive control of the measurement, not a discovery screen. Two things
confound constant (level) human-vs-mouse differences:
- **Probe efficiency.** Each gene's probe efficiency depends on its panel, and the panels differ by
  species (SI, "Probe panels").
- **Replication.** With one human donor and two mice, there is no species-level replication.

The notebook asks whether the level still recovers species differences that are already known. It
reports agreement and effect sizes. **No species-level p-value is computed on our data.**

**Our level statistic, L_g.**
- **Gene model.** The shared NB model of notebooks 66 and 67 (`pathway_pipelines.standard_gene_model`)
  for the human-vs-mouse species split, on notebook 37's inputs (`coordinate_inputs.pt_inputs` on
  SCF13): 12,272 structures and 11,071 genes. It uses offset log(library), basis_df 6, α from M_full
  and the 101-point grid.
  - It is refitted under this notebook's cache.
  - **Guard:** `T_spatial`, the flags, α and Δ must equal notebook 67's cached
    `gene_model_human_vs_mouse_species` (rtol 1e-6).
- **Δ_g(x).** The natural-log rate ratio, human over mouse, between the equal-specimen group curves.
  The grid has 101 equally spaced points on SCF13's common support, 0.053–0.836.
- **Grid weights.** The model has none, so each grid point gets weight 1/101, uniform on the
  coordinate. This is notebook 67's `log_ratio`.
- **Definition.** L_g = (1 / ln 2) · mean over x of [Δ_g(x) − m(x)], where m(x) is the median of Δ at
  x over the reference genes (converged and not separated). L is in log2 units: L = 1 means twice the
  human/mouse ratio of the typical gene.
  - Centring at each x removes three things that change along the PT and between species: depth,
    library composition and the mitochondrial share.
  - It also removes any efficiency difference that is the same for every gene of a panel.
  - **It cannot remove gene-specific probe efficiency.**
- **Pair spread.** The same statistic from the specimen curves, for the four human-section × mouse pairs
  (HUK1_COR1 and HUK1_MED1 × Ctrl1A2 and Ctrl1A4). Each pair is centred by its own median.
  - Reported: the minimum, the maximum, and whether all four pairs share the sign.
  - A pair is skipped where either specimen fit is separated or did not converge.
  - The two human sections come from one donor, so the spread covers sections and mice, not human donors.
- **Same-species reference** (descriptive). Ctrl1A2 − Ctrl1A4 and HUK1_COR1 − HUK1_MED1, by the same
  formula, with their 2.5–97.5% bands over genes.
- **Coordinate-free twin, P_g.**
  - Per specimen × segment: a pseudobulk of raw counts over the 12,272 structures, with
    log-CPM = log2((y + 0.5) / (N + 1) · 10⁶) and N the summed library.
  - Average over S1, S2 and S3 with equal weights.
  - Take the human mean minus the mouse mean (two specimens each), then subtract the median over the
    reference genes.
  - P is built exactly like the atlas statistic. It is the pre-specified sensitivity for L, and it gives
    the direction of separated genes.

**Flags.**
- **Non-converged genes** (21 in notebook 67's cache) are dropped everywhere.
- **Separated genes** (283) are those whose design loses rank on the structures with a count, for
  example a species or a specimen without counts.
  - They are left out of every correlation, band and median.
  - In direction checks they enter with sign(P_g) and with |P_g| in place of |L_g|.
- **Nearly absent genes** (930, detected in < 1% of one species' structures) stay in. Every criterion is
  also reported without them.

**Atlas reference, A_g.** Single-nucleus RNA-seq with poly-A capture and no probes. Donors are the
replicates.
- **Human: Lake et al. 2023 (KPMP v1.5).**
  - Normal cortex (region C, tissue "cortex of kidney"), author labels PT-S1/S2/S3.
  - The 7 donors with ≥ 50 nuclei in every segment: 164-10, 164-6, 18-142, 18-312, 3535, KRP446 and
    KRP460 (4 female, 3 male).
  - Ensembl ids are mapped to our mouse symbols first and symbols second, with duplicates summed
    (notebook 38's rule).
  - The pseudobulk is rebuilt under this notebook's cache at the Ensembl level. **Guards:** it equals
    notebook 38's cache on our genes and notebook 62's cache on all genes.
- **Mouse: CELLxGENE Census 2025-11-08, dataset 25818bf7.**
  - Segments 1/2/3 (`DATA/external/census_pt_segments/`).
  - **Adult males only, 3–6 months and 12 months: 6 donors**, each with ≥ 50 nuclei in every segment.
  - Why not all 12 males:
    - Our controls are adults.
    - The three 21-day males are prepubertal, so androgen-dependent PT genes are not yet induced.
    - Notebook 45 found that age imbalance drives null calls.
    - Notebook 59 used the same 6 donors.
- **Per donor.** log-CPM per segment = log2((y + 0.5) / (N + 1) · 10⁶), with N the donor-segment total
  over the analysed genes. Then the equal-weight mean over S1, S2 and S3. Equal weights mirror P and
  L, which do not weight segments by sampling. Our Visium sampling differs between segments and
  species.
- **Species coefficient.**
  - limma `lmFit` on the donor-level values (genes × 13 donors), design ~ species, then
    `eBayes(trend = TRUE, robust = TRUE)`. The coefficient b_g is in log2.
  - Centring: A_g = b_g − median(b). The median is over the analysed genes: our universe genes present
    in both atlases.
  - Test against the median: t = A_g / (stdev.unscaled · √s²_post), on df.total, with BH over the
    analysed genes. This is an atlas statistic, not one of ours.
- **Atlas-strong:** BH ≤ 0.05 and |A_g| ≥ 1.
- **Sensitivities:**
  - all 12 Census males;
  - Lake's 3 male donors only (all our specimens are male);
  - Census human cortex (dataset 09b518f9) in place of Lake: the 6 donors with ≥ 50 nuclei in both of
    its labels, with convoluted PT weighted 2/3 and S3 1/3 so that S1, S2 and S3 count equally.
- **Atlas ceiling** (descriptive). The Spearman correlation of A_g from Lake with A_g from Census human,
  against the same mice. The mice are shared, so this is an optimistic ceiling for (iii).

**Literature controls.**
- **Source.** `docs/results/pt-level-literature-controls.csv` (34 rows), sha256
  `f1aa71844552327e7f8b11905167fae34c205398c2499053711bcbee2ca06c4f`. It is used as written and never edited
  after the freeze. The search strategy and the excluded candidates are in
  `docs/results/pt-level-literature-controls-notes.md`.
- **What the controls measure.** 29 of the 30 rows with `direction_verified_for_mouse = yes` come from one
  study, Thakur et al. 2024 (Clin Pharmacol Ther, PMID 38711199). It used LC-MS/MS of 78 membrane transporters
  and enzymes in kidney membrane fractions from human (n = 15) and mouse (n = 14).
  - These are **protein** abundances, not mRNA.
  - The rodent samples are whole kidney and the human samples are tissue sections.
  - Values are pooled over sex.
  - So a sign disagreement can come from mRNA–protein discordance or from the tissue, not only from our
    measurement.
- **Primary set for (i).** Gene rows with `direction_verified_for_mouse = yes`. Rows marked `no` (L02, L07,
  L14, L16) are reported separately and do not enter (i). The pathway rows (L16, L17) are checked by the sign
  of their members' mean L. That check is descriptive.
- **Gene rows.** A row is in our universe when its symbol is one of the 11,071 genes: a mouse symbol
  directly, a human symbol through `ortholog_map_used.csv`.
- **Other rows.** Rows that name a non-1:1 gene or a family go to the family analysis. Pathway rows go
  to the pathway analysis. Rows that fit none of these are listed as not evaluable.
- **Correct sign.** A control has the correct sign when sign(L_g) equals its expected direction; for a
  separated gene, sign(P_g). Rows that expect no difference are reported only.

**Housekeeping yardstick.**
- **List.** Eisenberg & Levanon 2013 (`DATA/external/housekeeping/HK_genes.txt`, untracked; 3,804 human
  symbols; sha256 `df7aa10e635efb5b0df7dd607e5c11c91df2bdf3e321fb97c0c088ea96a60756`). Mapped by human symbol through the ortholog map, 3,184 of them are in our universe.
- **Band.** The 2.5–97.5% quantiles of L over the housekeeping genes that are converged and not
  separated.
- **Clearing the band.** A control clears it when its L lies beyond the band on its expected side.
- **Other bands.** The same band is computed for P (used for families) and for A. The band widths, ours
  against the atlas, are reported. The difference in width is what our two probe panels add to the
  spread of level differences.

**Probe and pairing sensitivities.**
- **Equal-probe genes.** (i)–(iii) are repeated on genes with equal included-probe counts in both panels
  (`results/pt_literature_deep/probe_counts.csv`, `probe_balance == balanced`: 9,090 of 11,071). They
  are repeated again on genes that also are probe-clean in both species (notebook 62: 9,019).
- **Probe ratio** (descriptive). The Spearman correlation of L_g − A_g with log2(human probes / mouse
  probes) over the unbalanced genes. A positive value means probe count leaks into the level.
- **Pairing.** (ii) and (iii) are repeated without the pairs that lie in non-1:1 orthogroups (notebook
  64's `dropped_pairs.csv`, 419 in our universe). In those pairs a mis-chosen mouse partner shifts our
  data and the atlas alike.

**Families** (descriptive; the expectation is explicitly weaker).
- **Orthogroups.** Notebook 62's HCOP components (support ≥ 3) of class 1:many, many:1 or many:many,
  with ≥ 1 panel member in each species and unequal HCOP member counts: 720 orthogroups.
- **Family sum.** Per structure, the sum of raw counts over every panel member of the species, from
  `DATA/tubule_by_gene/<sample>_tubule_by_gene_caleb.h5ad` joined on structure id.
  - Summing over all members makes within-family cross-hybridisation cancel, so this is the primary.
  - Probe-clean members only is a sensitivity.
- **Statistic.** F_f is P's statistic on the family sum: equal-segment log-CPM against the gene model's
  library, human minus mouse, minus the median of P over the reference genes.
- **Expectation.** sign(F_f) = sign(n_human − n_mouse), with n the HCOP member counts.
  - Reported as the share that agree, over the eligible families whose summed counts are detected in
    ≥ 2% of the structures of both specimens of at least one species.
  - Reported separately for the pre-named families:
    - expanded in mouse: Cyp4a, Akr1c, Slco1a, Sult2a, Cyp2d, Ces1, Nat8, Cyp2j and Cyp2c;
    - expanded in human: Sult1a;
    - no expectation: Ugt2b, which has 8 and 8 HCOP members.
- **Atlas family sums.** The same statistic in Lake and in the Census adult males. The Census side needs
  the members fetched first; see the inventory.
  - Reported: the Spearman correlation of F with the atlas F, and the sign agreement.
  - If the fetch is not done, the atlas family check is limited to the families already complete. It is
    then labelled as a biased subset.

**Pathways** (descriptive).
- Each of notebook 66/67's 1,513 sets (10–300 tested members) is summarised by the mean L of its members
  that are converged and not separated.
- The atlas mean of A is taken over the same members present in the atlas. The Spearman correlation
  between the two is reported across pathways.
- Literature pathway controls are checked by the sign of their mean L, and their rank among the 1,513 is
  reported.
- No pathway is tested.

**Criteria** (fixed now):
- **(i) Literature.** At least 80% of the primary literature gene controls in our universe have the
  correct sign.
  - This needs ≥ 10 evaluable controls; with fewer, (i) is "not evaluable".
  - Reported beside it: the Wilson 95% interval; how many controls clear the housekeeping band; how many
    have the same sign in all four pairs; the equal-probe subset.
- **(ii) Atlas direction.** Among genes that are atlas-strong and have |L_g| ≥ 1, at least 80% have the
  sign of A_g.
  - This needs ≥ 30 genes; with fewer, (ii) is "not evaluable".
  - Reported beside it: all atlas-strong genes regardless of |L|, and every sensitivity above.
- **(iii) Gene-level agreement** (descriptive). The Spearman correlation of L_g with A_g, with a 95%
  interval from 2,000 gene bootstraps (seed 68).
  - Computed on all genes that are converged, not separated and present in both atlases, and on the
    equal-probe genes.
  - It is set against the atlas ceiling and P's correlation with A.
- **(iv) Families** (descriptive). The share of families in the expected direction, and the agreement of
  F with the atlas.
- **Overall.** "Level recovers known biology" if (i) and (ii) both pass. The primary statistic L
  decides; P is reported beside it.

**Why 80%.**
- Chance is 50%.
- Shape reached 94% direction agreement against the same atlases
  (`docs/results/pt-pathway-literature-novelty.md`). Shape cancels probe efficiency; level does not.
- 80% allows probe-driven sign flips in up to a fifth of strong genes. It still rejects a measurement
  in which probe efficiency rivals biology.
- |L| ≥ 1 and |A| ≥ 1 restrict (ii) to two-fold differences on both sides, where the sign is
  meaningful.

**What counts as failure.**
- A criterion below 80%, or not evaluable, fails. The result is then reported as it is.
- The thresholds, the gene universe, the atlas donors and the statistic are not changed after the fits.
  Any later check is labelled post hoc.
- **A pass licenses no gene-level claim.** It shows that the level agrees with known biology on
  average, not that any one gene's level difference is real. Level differences are then reported with
  the housekeeping band and the pair spread, never with a p-value.

**Readings, fixed now.**

| (i) | (ii) | Reading |
|---|---|---|
| Pass | Pass | Level is a usable descriptive measurement. |
| Fail | Pass | Level agrees with the atlas genome-wide but not with the curated controls, which are protein measurements. List the failing controls, their atlas sign and their probe balance. |
| Pass | Fail | The agreement rests on a few curated genes and does not hold genome-wide. Level stays out of the paper's claims. |
| Fail | Fail | Level differences in our data are not interpretable. The paper keeps the SI's statement. |

**Output.**
- Notebook: `analysis/notebooks/68_pt_species_level_control.ipynb`, logic version `68.level.1`.
- Results: `results/pt_species_level_control/`.
- Reusable logic goes in `pseudospace/species_level.py`, tested by `tests/test_species_level.py`. It
  covers the grid-centred level, the pair levels, the equal-segment pseudobulk log-CPM, the family sums
  and the limma-trend atlas coefficient through rpy2.
- **Assumptions.**
  - The coordinate is treated as exact.
  - Specimens are the replicates on our side and donors in the atlases.
  - The human "replicates" are two sections of one donor.
  - **In the atlases, species is also confounded with dataset:** lab, protocol and donor ages. Agreement
    means the result is consistent across technologies. It does not prove biology.
  - snRNA measures nuclear RNA and Visium HD measures whole-cell RNA. A nuclear-fraction bias that is
    conserved between species cancels in the human-vs-mouse difference.

## Notebook 68 results

`analysis/notebooks/68_pt_species_level_control.ipynb` (logic version `68.level.1`) writes to
`results/pt_species_level_control/`. It followed the protocol above (commit aa15069).

**Stage caches:**
- `gene_model_species__9bfedfdd733d1c53`;
- `lake_cortex_pt__899c17f45c1b5192`;
- `visium_family_sums__981625d9f0e95c86`.

**Guards.**
- The refitted gene model reproduces notebook 67's cache: statistics, flags, α, Δ and curves, rtol 1e-6.
- The Lake pseudobulk equals notebook 62's cache (all genes) and notebook 38's cache (our genes).
- The Census family fetch has the main file's groups and totals.

| Criterion | Result | Status |
|---|---|---|
| (i) primary literature controls with the expected sign | 17 of 27, 63% (Wilson 44–78%) | **fails** |
| (ii) atlas-strong genes with \|L\| ≥ 1 sharing the atlas sign | 3,584 of 3,984, 90.0% (88.9–90.9%) | passes |
| (iii) Spearman L vs A, with a 2,000-gene bootstrap | 0.60 [0.58, 0.61]; atlas ceiling 0.97 | descriptive |
| (iv) families in the member-count direction | 198 of 410, 48% | descriptive |

**Verdict.** The overall criterion is not met. The reading is (fail, pass): level agrees with the atlas
genome-wide but not with the curated controls, which are protein measurements.

**(i) Literature controls.**
- **Coverage.** 27 of the 29 gene rows verified for mouse are in our universe. Slc47a2 and Ephx3 are not.
- **Beside the share:**
  - 8 controls clear the housekeeping band;
  - 14 have the expected sign in all four section × mouse pairs;
  - the equal-probe subset gives 13 of 23;
  - without nearly absent genes, 15 of 25.
- **The failing controls**, with their atlas sign and probe balance (`tables/failing_literature_controls.csv`):
  - All 10 have 3 probes in each panel, and none is nearly absent.
  - **The atlas disagrees with the control too** in 6: Abcc4, Fmo1, Slc16a1, Slc22a6, Slc22a7 and
    Slc22a8. Two examples:
    - Slc22a6 (OAT1): ours +0.2, atlas +2.7.
    - Fmo1: ours −1.9, atlas −2.3.
  - **The atlas agrees with the control** in 4:
    - **Abcc2** is the one clear failure of our measurement: ours +0.6 against an atlas −2.6 (BH < 10⁻³).
    - **Gusb** is small on both sides: ours −0.8, atlas +0.7.
    - **Slc22a2 and Slc22a13** have an atlas level that does not differ from the typical gene (BH 0.31 and
      0.76).
- **The atlas itself** agrees with only 19 of the 27 protein-based controls. Our L agrees with the atlas sign for
  21 of the 27 (post hoc).
- **Rows not verified for mouse** (reported only): Kyat1 has the expected sign and Abcc1 does not; Slc22a3 is
  outside our universe.
- **Pathway rows** (descriptive):
  - L17 (Slc22a1 + Slc22a2) has mean L −2.9. That ranks below all 1,513 pathway means, in the expected
    direction.
  - L16 (13 transporters, not verified) has mean L −0.35, rank 1,131 of 1,513, also in the expected direction.

**(ii) and (iii) sensitivities.**

| Setting | (ii) share | (iii) Spearman |
|---|---|---|
| Primary | 90.0% (n = 3,984) | 0.60 |
| All atlas-strong genes, any \|L\| | 81.0% (n = 6,071) | — |
| P in place of L | 89.8% | 0.60 |
| Equal probes / and probe-clean in both | 91.0% / 91.0% | 0.62 / 0.62 |
| Without nearly absent genes | 88.6% | 0.55 |
| Without separated genes | 89.5% | — |
| Without pairs in non-1:1 orthogroups | 90.0% | 0.60 |
| All 12 Census males | 90.2% | 0.59 |
| Lake's 3 men | 90.9% | — |
| Census human in place of Lake | 91.6% | 0.61 |

- **P is not an independent check.** L and P have Spearman 0.999.
- **Probe count leaks into the level.** Over the 1,918 genes with unequal probe counts, L − A has Spearman 0.45
  [0.42, 0.49] with log2(human / mouse probes).
- **Bands** (2.5–97.5%, log2):

  | Band | Range | Width |
  |---|---|---|
  | Housekeeping, our L | −3.1 to 2.6 | 5.7 |
  | Housekeeping, atlas A | −2.7 to 2.7 | 5.4 |
  | Mouse − mouse | −0.3 to 1.1 | 1.4 |
  | Section − section | −0.5 to 0.4 | 1.0 |

  Housekeeping genes therefore spread almost as widely in the probe-free atlas as in our data, so most of the
  band is between-species or between-dataset spread, not probe noise. The same-species bands are 4–6 times
  narrower.

**(iv) Families** (descriptive; 720 eligible orthogroups plus Ugt2b).
- **The genome-defined expectation fails in both datasets.**
  - Our F has the expected sign in 198 of 410 detected families (48%); with probe-clean members only, 50%.
  - Mouse-expanded families are mouse-higher in 45%, and human-expanded families human-higher in 52%.
  - The atlas F has the expected sign in 52%.
- **Our F tracks the atlas F.** The signs agree in 296 of 407 complete families (73%), with Spearman 0.67
  [0.60, 0.74].
- **Pre-named families** (F, then atlas F):
  - Mouse-higher in both: Slco1a (−7.4, −10.7), Cyp2d (−10.3, −7.0), Cyp2j (−8.0, −9.6), Ces1 (−8.5, −9.8),
    Akr1c (−5.0, −5.4) and Nat8 (−1.5, −3.0).
  - Cyp2c: +1.9 and +2.4, against the expectation in both.
  - Cyp4a: +0.6 against −1.0.
  - Sult2a is not detected.
  - **Sult1a disagrees:** −4.2 against +4.5. The human members SULT1A1 and SULT1A2 have one probe each and are
    detected in under 0.5% of human structures (post hoc). This looks like a probe failure.
  - Ugt2b, which has no expectation: −1.9 against +1.5.

**Pathways** (descriptive). Mean L vs mean A has Spearman 0.70 [0.66, 0.73] across the 1,513 sets. The gene
sets overlap.

**Readings of the protocol** (implementation details not written in it):
- **Census symbol audit.** Census columns whose symbol now belongs to another Ensembl id than our panel's
  measured another gene. They were dropped from that atlas:
  - in the mouse file: Arvcf, Ccdc121, Eppk1, Lin54, Tctn2, Tusc3 and Ugt3a1;
  - in the human file: PAXX.
  - In the existing `census_pt_segments` mouse file, the **Ugt3a1 column is Ugt3a2**, because Census swapped
    the two symbols. Earlier notebooks that read Ugt3a1 from that file (38, 43 and 64) compared the wrong gene.
- **Analysed atlas genes.** Our genes mapped in Lake and listed under our gene in Census: 11,036; 11,048 for
  Census human against mouse.
- **Families.** The member sets are notebook 62's panel members of the orthogroup's non-1:1 classes: 1,652 mouse
  and 1,116 human. The Census mouse members were fetched by Ensembl id (`analysis/scripts/fetch_census_pt_families.py`);
  2,016 of 2,018 were found.
- **Literature pathway rank.** The rank places the members' mean L among the 1,513 pathway means.
- **Post hoc** (section 12 of the notebook): the controls' agreement with the atlas, and the probe counts and
  detection of the pre-named families' members.

**What this cannot tell you.**
- The controls are mostly one protein study of whole mouse kidney, so they test mRNA–protein concordance as much
  as our measurement.
- In the atlases, species is confounded with dataset. Agreement with them shows consistency across technologies,
  not biology.
- With one human donor there is no species-level inference. Every number here is descriptive.

## Notebook 69 protocol: pathway tests without mitochondrial reads in the offset

Frozen on 2026-10-10, before any fit. The user approved the change.

**The failure.**
- **The library holds mitochondrial reads.** Every NB gene model so far uses offset log(library), where the library
  is the total over all orthologs measured in both inputs (`coordinate_inputs.pt_inputs`; notebook 66/67's AKI
  construction). That total includes the mitochondrially encoded genes:
  - human vs mouse: 12 genes (mt-Atp6, mt-Co1, mt-Co2, mt-Co3, mt-Cytb, mt-Nd1, mt-Nd2, mt-Nd3, mt-Nd4, mt-Nd4l,
    mt-Nd5 and mt-Nd6); mt-Atp8 is not measured in human;
  - AKI vs control: the same 12 plus mt-Atp8, which is measured in mouse but is outside the tested universe
    (notebook 12's universe requires both species).
- **Their share changes along the PT, and differently by species.** On SCF13, from the first to the last fifth of
  the coordinate, mt reads fall from 21.5% to 12.4% of the human library and rise from 3.4% to 5.5% of the mouse
  library. Both specimens of each species agree.
- **This is not mitochondrial content.** The 61 nuclear-encoded OXPHOS subunits stay flat within ±0.15 log2 in
  both species. Part of the mouse late rise travels with TAL spillover (Umod, Slc12a1).
- **The effect.**
  - The offset tilts every gene's human/mouse log ratio toward S3. From the shares above, removing the mt reads
    changes [late − early] of every gene's Δ by about −0.13 natural-log units.
  - Genes that are not significant drift too, by about +0.1.
  - The tilt pushes calls toward genes "rising in human toward S3".
- **The mt genes themselves.** In notebook 67 they are among the strongest species-shape genes (z 6.6–7.9), but
  they belong to none of the 1,513 gene sets, because the Enrichr libraries list no mt-encoded gene. Dropping them
  from the cameraPR background alone changes the normalized calls from 70 to 71 (+ KEGG Glycerolipid metabolism).
- **AKI is expected to change little.** mt reads go from 5.4% to 5.8% of the AKI library and from 3.7% to 6.1% of
  the control library, a predicted tilt of about −0.02.
- **Notebook 68 is not rerun.** It centres the level at each grid point on reference genes.

**The single change.** Genes whose mouse symbol starts with `mt-` are removed from two places:
- **The library**, and so from the NB offset and the clustering method's half-count floor:
  - human vs mouse: the total over the orthologs measured in both inputs, without the 12 mt genes
    (`pt_inputs(..., exclude_genes=...)`);
  - AKI: the total over the measured panel genes, without the 13 mt genes.
- **The tested universe:** human vs mouse 11,071 → 11,059 genes; AKI 9,107 → 9,095.

The gene sets are rebuilt on the reduced universe with the same rule (10–300 tested members). No mt gene is in
any set, so the sets are expected to stay 1,513 and 1,429, identical to notebooks 66 and 67.

**Held fixed, as frozen in notebooks 66 and 67.**
- **Inputs.** The same structures, positions, specimens and segments (12,272 and 16,991 structures), and the same
  condition split and two balanced relabelings (`partition_groups`).
- **Gene model.** `pathway_pipelines.standard_gene_model`:
  - NB GLM on raw counts, offset log(library) with the new library;
  - 6-df cubic B-spline with knots at the position quartiles, and the prior weights wᵢ = n / (2 · J · n_j);
  - the 101-point grid.
  - α is re-estimated per gene under the condition split's M_full with the new offset, and reused for both
    relabelings. This is the frozen rule applied to the new offset, not a second change.
- **Normalized method** (notebook 67's primary):
  - the denominator LR_spec (`specimen_shape_lr`), with df₂ = 12 for the condition split and 6 for a relabeling;
  - `limma::squeezeVar` with covariate log mean raw count per structure. Raw counts do not depend on the library,
    so the covariate is unchanged.
  - F and z (`moderated_f`);
  - `cameraPR` with `inter.gene.cor` 0.01 and `use.ranks` FALSE; a call is Up with FDR ≤ 0.05;
  - the gene-level q is the BH of the moderated F's p.
- **Clustering method** (notebook 66's method B), from `pseudospace_reconstruction` at commit 59da04b:
  - stages 06–07 on the shared specimen curves (`clustering_first_inputs`), with the half-count floor over each
    specimen's summed library;
  - `CompareParams(prior_k = 0)`, `PathwayParams(min_genes = 10, max_genes = 300)` and node spacing 1/29;
  - the zone cuts from the package's marker rule on the control mice of the condition split;
  - the condition run with the default null (both relabelings) and with each single-split null; each relabeled
    run with the fair null (the other relabeling).
- **Redundancy clusters.** The package's `cluster_gene_sets` (1 − Jaccard, average linkage, cut at 0.5),
  recomputed on the reduced universe. Expected to equal notebook 66's.
- **Pass rules.**
  - Relabeling: NCDR < 0.25, and each relabeling calls < 5% of tested pathways.
  - Calibration: the pooled share of relabeled one-sided p ≤ 0.05 is ≤ 0.075. For the clustering method, notebook
    61/66's list-based convention: the best matched ORA p × the number of lists.
- **Not changed, because neither method's calls use it.**
  - The log-normalized matrix, normalized inside `rebuild_pt_expression` against the full library. It feeds only
    the zone cuts and descriptive mean expression. The zone cuts are guarded to equal notebook 66's.
  - The covariate strata of the joint test.
- **Not rerun.** Notebook 66's cameraPR and GSEA on raw `T_spatial`, the NB joint test, and notebook 67's
  sensitivities.

**Guards** (not cached).
- Exactly 12 mt genes leave the human vs mouse universe. 12 leave the AKI universe and 13 leave the AKI library.
- **The inputs before the change equal notebook 66's.** With the mt genes kept, the inputs reproduce notebook 66's
  recorded stage keys for its gene models: counts digest, library, position, specimen, genes, group and nuisance.
  The keys are recomputed with the code digest stored in notebook 66's cache sidecars.
- **After the change, only the intended parts differ:**
  - the same structures, positions, specimens and segments;
  - the non-mt count columns are identical;
  - the genes are notebook 66/67's minus the mt genes;
  - the gene sets, redundancy clusters and zone cuts are identical.
- df₂ = 12 and 6, from the design ranks.
- Recomputing everything with the mt genes kept, to reproduce notebook 66/67's calls, is not required.

**Pre-stated checks, each with its expected result.** All are fixed now, before looking.

Definitions:
- d_g, for gene g, is the mean of Δ_g(x) over the grid points in the last 20% of the grid's range, minus the mean over
  the first 20% (21 of the 101 points at each end).
  - Δ is the natural-log rate ratio, case over reference (human/mouse; AKI/control), from the condition split.
  - d_g > 0 means "rising in human (AKI) toward S3".
- **Before** values come from notebook 67's cached `gene_model_<dataset>_species` (Δ, converged, separated), which
  equals notebook 66's. They also use notebook 67's gene-level q (`gene_statistics_<dataset>_primary.csv`, condition
  split). **After** values come from notebook 69.
- Only genes that are converged and not separated enter.

**(a) Null drift.**
- The statistic is the median of d_g over genes with q > 0.5 on the normalized method in both runs. The mt genes are
  absent after the change, so they are excluded from both.
- **Criterion** (human vs mouse): |after| ≤ 0.05 natural-log units, about 5% over the coordinate. A ratio
  criterion (|after| ≤ |before| / 3) was considered and rejected before freezing: the change shifts every d_g by
  about −0.13 deterministically, so with a "before" near +0.1 a ratio rule sits on its own pass/fail line.
- Reported beside it:
  - the predicted "after", the "before" minus the tilt computed from the mt shares;
  - the same median with each run's own q > 0.5 genes;
  - the median Δ curve of those genes along the PT, before and after.
- AKI is reported and not judged, because its predicted tilt (about 0.02) is too small for the ratio to mean
  anything.

**(b) Direction balance.**
- Among genes with q ≤ 0.05 on the normalized method (condition split), count those with d_g > 0 (rising in
  human toward S3) and d_g < 0, before and after, on both datasets.
- **Expected:** the rising share falls after the change on human vs mouse. This check is descriptive and not a pass rule.

**(c) Calls and pass rules**, for both methods on both datasets:
- what is reported: condition calls (pathways, and the redundancy clusters they touch), the calls of each
  relabeling (the fake-group calls), NCDR, the relabeling rule, the calibration share and whether it is within
  0.075;
- the clustering rows: the default null (its calls), each single-split null, and the fair relabeled runs.
- **Expected:**
  - The normalized method still passes on human vs mouse. It still fails calibration on AKI: the change does not
    touch the coherent specimen programs behind that failure (notebook 67, section 6.1).
  - The clustering method stays specific under the fair null. Its calls stay nearly unchanged, because it
    re-references every gap to the typical gene (the median gap at each position), which removes any shift shared
    by all genes.

**(d) Overlap between the methods.**
- **Inputs.** The normalized calls (cameraPR on z) against the clustering calls (its ORA, default null; notebook
  66's "B ORA"), on both datasets.
- **The table:** both, normalized only, clustering only, and neither (tested pathways called by neither). It also
  gives the Jaccard index, over pathways and over redundancy clusters.
- **Sensitivity: without the 8 proteasome subunits.** Psma3, Psma4, Psmb5, Psmb6, Psmc2, Psmc6, Psmd11 and Psmd14
  are removed from both methods' universes, as in the earlier sensitivity:
  - **Normalized.** cameraPR on the same z without the 8 genes, with no squeezeVar refit. Sets below 10 members
    are dropped.
  - **Clustering.** The ORA (`ora_by_group` with its `PathwayParams`) is rerun with the 8 genes outside the
    universe. It keeps the gene groups and clusters of the condition run unchanged.
- **Before values,** human vs mouse:
  - all genes: Jaccard 0.077 over pathways and 0.138 over redundancy clusters;
  - without the 8 subunits: 0.152 and 0.161.
- **Expected:** the overlap stays low. The Jaccard over redundancy clusters stays ≤ 0.2 on both datasets, with and
  without the 8 subunits. Descriptive.

**(e) Lost and gained calls.**
- **What is compared.** For each method and dataset, the calls of notebook 69 against those of notebook 67 for
  the normalized method, and against notebook 66's default-null ORA for the clustering method.
- **Per pathway:** its FDR (normalized) or best matched q (clustering), before and after; and the median d_g over
  its tested members, before and after (its median direction toward S3).
- **Expected:**
  - The normalized method's lost calls are mostly pathways whose genes rise in human toward S3 (median d_g before
    > 0). Its gained calls are mostly falling ones.
  - The clustering method loses and gains few calls.

**Decision rule.**
- **If (a) passes,** notebook 69's calls replace notebook 66/67's calls as the primary for both methods on both
  datasets. The pass rules in (c) are reported for the new calls; this does not reopen notebook 67's decision.
- **If (a) fails,** the result is reported and notebooks 66/67 stay primary. Nothing is tuned.
- Any later check is labelled post hoc.

**Output.**
- Notebook: `analysis/notebooks/69_pt_pathways_without_mitochondrial_reads.ipynb`, logic version `69.no_mito.1`.
- Results: `results/pt_pathways_without_mito/`.
- New option: `pt_inputs(..., exclude_genes=())`, which removes the named genes from the library and from the
  universe. Its default leaves every existing input, and so every existing cache key, unchanged. It is tested in
  `tests/test_coordinate_inputs.py`.

## Notebook 69 results

`analysis/notebooks/69_pt_pathways_without_mitochondrial_reads.ipynb` (logic version `69.no_mito.1`) writes to
`results/pt_pathways_without_mito/`. It followed the protocol above (commit 4fdcdf4). Its section 10 lists how the
protocol was read. Section 8.1 is post hoc.

**Stage caches:**
- `gene_model_human_vs_mouse_species__1c1cadaef98bb99c`;
- `gene_model_aki_vs_control_species__9977e65777d6275d`;
- `specimen_shape_human_vs_mouse_species__486d509e3fe47198`;
- `method_b_human_vs_mouse_species__null_default__4b4253d8cfc3a794`.

**Guards.**
- 12 mt genes left the human vs mouse library and universe (11,071 → 11,059 genes). 13 left the AKI library and 12
  the AKI universe (9,107 → 9,095). Each library fell by exactly those genes' counts.
- With the mt genes kept, the inputs reproduce notebook 66's and 67's recorded gene-model stage keys.
- The structures, positions, non-mt counts, gene sets (1,513 and 1,429), redundancy clusters, zone cuts and grid are
  unchanged.
- df₂ is 12 and 6.
- The clustering ORA recomputed from the gene groups reproduces both this run's calls and notebook 66's.

**Verdicts.**

| Check | Result | Status |
|---|---|---|
| (a) Null drift, human vs mouse: median d_g over the 1,875 genes with q > 0.5 in both runs | +0.079 before, −0.048 after (predicted −0.055) | **passes** (\|after\| ≤ 0.05) |
| (a) AKI, reported only (1,419 genes) | +0.016 → −0.002 (predicted −0.006) | — |
| (b) Significant genes rising in the case group toward S3 | human vs mouse 2,624 of 3,835 (68%) → 1,890 of 3,614 (52%); AKI 55% → 51% | falls, as expected |
| (c) Pass rules | normalized: passes on human vs mouse, still fails calibration on AKI; clustering: specific under the fair null | as expected, except that the clustering calls moved (below) |
| (d) Overlap between the methods | Jaccard over redundancy clusters 0.093–0.198 | ≤ 0.2, as expected |
| (e) Lost and gained calls | normalized: few lost, gains mostly falling; clustering: 36 human vs mouse calls lost | partly as expected |

- **The drift is removed and slightly overshot.** The mt tilt (−0.13) is larger than the drift it removed (+0.08), so
  the null genes now drift slightly the other way. Over each run's own q > 0.5 genes the drift is +0.039 before and
  +0.001 after.
- **The 12 mt genes were all among the significant genes before.**
- **Decision.** Check (a) passes, so notebook 69's calls replace notebooks 66/67's as the primary for both methods on
  both datasets. The normalized method's AKI calibration failure stands.

**(c) Calls and pass rules.** Each cell gives the calls (redundancy clusters); relabeled calls; NCDR; calibration
share (limit 0.075).

| Method | Human vs mouse, 66/67 → 69 | AKI vs control, 66/67 → 69 |
|---|---|---|
| Normalized (cameraPR on z) | 70 (68); 0, 0; 0; 0.044 → **97 (90); 0, 0; 0; 0.044 · passes** | 93 (86); 17, 11; 0.15; 0.102 → **99 (91); 16, 10; 0.13; 0.101 · fails calibration** |
| Clustering, default null | 167 (72); 0, 0 → **132 (51); 0, 0** | 464 (262); 2, 0 → **455 (253); 2, 0** |
| Clustering, single-split nulls | 162, 159 → 162, 46 | 469, 482 → 465, 468 |

The clustering method's list-based calibration share is 0.003 (human vs mouse) and 0.005 (AKI), before and after.

**(d) Overlap between the normalized and clustering calls** (condition split). Each cell gives both / normalized only /
clustering only / neither (Jaccard). Pathways are counted first; redundancy clusters follow in brackets.

| Dataset | Universe | 66/67 | 69 |
|---|---|---|---|
| Human vs mouse | all genes | 17 / 53 / 150 / 1,293 (0.077) [17 / 51 / 55 / 854 (0.138)] | 22 / 75 / 110 / 1,306 (0.106) [21 / 69 / 30 / 857 (0.175)] |
| Human vs mouse | without the 8 proteasome subunits | 14 / 57 / 21 / 1,421 (0.152) [14 / 55 / 18 / 890 (0.161)] | 20 / 77 / 14 / 1,402 (0.180) [20 / 70 / 11 / 876 (0.198)] |
| AKI vs control | all genes | 33 / 60 / 431 / 905 (0.063) [33 / 53 / 229 / 584 (0.105)] | 30 / 69 / 425 / 905 (0.057) [30 / 61 / 223 / 585 (0.096)] |
| AKI vs control | without the 8 proteasome subunits | 34 / 59 / 430 / 906 (0.065) [34 / 52 / 229 / 584 (0.108)] | 29 / 70 / 423 / 907 (0.056) [29 / 62 / 222 / 586 (0.093)] |

The 66/67 rows without the subunits reproduce the earlier sensitivity exactly (71 and 35 calls, 14 shared).

**(e) Lost and gained calls** (`tables/check_e_*.csv`). Median d_g is over each pathway's tested members; > 0 means
rising in the case group toward S3.

| Dataset · method | Kept | Lost (rising before) | Gained (falling after) |
|---|---|---|---|
| Human vs mouse · normalized | 69 | 1 (1) | 28 (20) |
| Human vs mouse · clustering | 131 | 36 (26) | 1 (1) |
| AKI · normalized | 91 | 2 (2) | 8 (8) |
| AKI · clustering | 443 | 21 (14) | 12 (1) |

- **Normalized, human vs mouse.**
  - The one lost call is KRAS Signaling Up (FDR 0.047 → 0.052).
  - 24 of the 28 gains had FDR 0.05–0.10 before. 14 of them already fell in human toward S3 before.
  - Gains include cholesterol and bile acid synthesis, steroid metabolism, mTORC1 signalling, very-long-chain fatty
    acid β-oxidation, glutathione synthesis and vitamin digestion. Among the rising ones are mineral absorption,
    vitamin D metabolism and EMT.
  - The leading calls are unchanged: fatty acid metabolism, inorganic ion and amino acid transport, branched-chain
    amino acid degradation, xenobiotic and arginine/proline metabolism.
  - The Spearman correlation of pathway p before and after is 0.87 (AKI 0.96).
- **Normalized, AKI.**
  - Lost: Myc Targets V1 and Processive Synthesis On Lagging Strand, both rising.
  - Gained, all falling: among them oxidative phosphorylation, Parkinson disease and NAFLD.
- **Clustering, human vs mouse.** The lost calls are mostly Reactome sets that hold proteasome subunits: antigen
  processing, neddylation, deubiquitination, CLEC7A, the ER-phagosome pathway and IL-1 signalling. Cell-cycle sets,
  p53 and EMT are lost too.
- **Clustering, AKI.** 24 of the 33 changes sit at a matched q of 0.05–0.09 on the other side.

**Post hoc (section 8.1): why the clustering calls moved.**
- **Before.** In each human vs mouse condition run, all 8 proteasome subunits sat in one Ward shape cluster. That
  cluster carried 126–130 of the run's calls.
- **After.** The subunits are split between two clusters:
  - with the default null, D4 (126 calls) becomes a 5 / 3 split, and the run falls from 167 to 132 calls;
  - with the relabeling-2 null, D5 (127 calls) becomes a 4 / 4 split, and the run falls from 159 to 46;
  - with the relabeling-1 null, the split is 6 / 2, the larger cluster still carries 129 calls, and the run keeps 162.
- **The stable part.** The gene group "Down · gradient toward S3", which does not depend on the clustering, carries
  108 calls before and 107 after.
- **Without the 8 subunits** the clustering calls are 35 before and 34 after.
- **So the clustering count rests on whether 8 genes land in one cluster.** That is a discrete outcome of the Ward k
  rule; the mt change moved it, but most of the difference is not about mt.
- AKI's clustering calls come from gene groups (Up · uniform, Up · gradient toward S3) and are stable.

**Figure.** `figures/fig1_null_drift_and_overlap.{png,pdf}`:
- a–b: the median Δ along the PT of the null genes, before and after;
- c–d: the 2 × 2 counts of pathways called by the two methods, after the change, with the before counts in brackets.

**What this cannot tell you.**
- Check (a) passes by 0.002.
- The offset now ignores mt reads entirely. The mt share still varies between structures within a specimen.
- The relabelings mix the conditions, so they cannot see an artefact that follows the condition.
- With two specimens per group, and one human donor, every result is descriptive.
