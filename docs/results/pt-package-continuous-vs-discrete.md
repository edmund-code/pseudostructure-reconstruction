# Continuous against discrete on the package pseudospace (notebook 74)

Part of the switch to the package pseudospace; see `pt-package-pseudospace.md`. Frozen on 2026-10-12, before any
fit of notebook 74. Notebook 74 runs after notebooks 71, 72 and 75 and reads their outputs.

## Notebook 74 protocol: what continuous analysis along the PT adds to S1/S2/S3

**Question.** On the package pseudospace, what does a continuous analysis along the PT show that the S1/S2/S3 segment
analysis does not? Which of those additions survive step-only nulls, a placebo coordinate and the coordinate's own
positional noise?

**Why.** Notebooks 57–61 answered this question on SCF13 with Gaussian gene models on log1p. Decision D2 replaces every
gene model with the NB model on raw counts, with the mt-free offset. Decision D1 replaces the units and labels. So 74
is not an input swap: every Gaussian statistic gets its NB counterpart, defined below before any fit. Notebook 72
delegated two items to 74's information side: WSR and the placebo claim gate (its section 8).

**Labels.**
- Notebook: `analysis/notebooks/74_pt_package_continuous_vs_discrete.ipynb`.
- Logic version: `74.continuous_vs_discrete.1`.
- Results: `results/pt_package_continuous_vs_discrete/`, with `tables/`, `figures/`, `stage_cache/`, `refits/placebo/`
  and `run_record.json`.

**Scope.** Human vs mouse, as in 57–61.
- The three pathway-level ladder rows (R2) also run on AKI, descriptively, for 75's method table (75's O5).
- Everything else is human vs mouse only. AKI along the PT belongs to 76.
- The human "replicates" are two sections of one donor: between sections, not between donors.

### 1 · Inputs (read-only)

| Source | What 74 reads |
|---|---|
| 70 | `package_pseudospace/hvm/pt_units.h5ad` through `package_inputs.hvm_inputs` (71's call, mt-free): positions, counts, labels, `dpt_pseudotime`. The base run's integration features (pass 1 ∪ pass 2 `highly_variable_for_harmony`) |
| 71 | stage cache `gene_model_human_vs_mouse_{species,swap_*}`: α, group curves, Δ(x) and HC3 SE, specimen curves (101-point grid). `tables/gene_statistics_human_vs_mouse.csv`: z, q, T_spatial, LR_spec. `pathway_tests_human_vs_mouse_cameraPR_*.csv`, `redundancy_clusters_*.csv`, `gene_sets_*.json`. AKI equivalents for R2 |
| 72 | `tables/noise_budget.csv`, `criteria_hvm.csv` (C14), `registration.csv`. Refit payloads `refit_hvm_{fold0,fold1,fold2,half_A,equal_depth}__*.npz`, read through `package_refits.refit_position`. The P4c cell construction (evaluation genes E, oracle covariates) for WSR. Optional: `unit_set_M.csv` |
| 75 | `discrete_comparators` (import). `gene_statistics_segment_{ds}.csv` and `pathway_tests_segment_{ds}.csv` (the segment F-test). `pathway_tests_deseq2_{wholept,segments}_{ora,gsea}_{ds}.csv` and `programs_{ds}.csv` |
| External | Census mouse PT nuclei (`data/external/census_pt_cells/mouse_pt_cells.h5ad`, 59's rules). HGNC gene groups 752 and 417 in `data/external/hgnc/` (R7). `pt_literature_deep/probe_counts.csv` |

**Guards (not cached).**
- G1: `hvm_inputs` reproduces 71: the counts digest of 71's cache key, positions, 11,055 genes and 1,513 sets.
- G2: `nested_nb_lr` (level → full) with 71's α reproduces 71's cached condition T_spatial to 10⁻⁸.
- G3: 75's segment tables cover exactly 71's universe and sets. Checked at the start of each section that reads them
  (R2, R3, R4 (iii)); such a section stops if they are missing.
- G4: 72's C14 verdict, recomputed from `noise_budget.csv`, equals `criteria_hvm.csv`. Checked at the start of each
  section that reads 72 (R4 count split and (iii), R6, the claim tables); such a section stops if 72 is missing.
- G5: the half-B counts are regenerated with 72's own call: `poisson_count_split` on the base `units_qc` counts, ε 0.5,
  `default_rng(35)`. They are checked against 72's recorded digest where 72 records one.
- G6: every arm, design and refit reports df₁ and df₂ from the design ranks. The expected values are asserted.

### 2 · Shared rules

**Gene model (D2, the frozen NB note).**
- NB GLM on raw counts; offset log(mt-free library); prior weights n / (2·J·n_j) in two-group designs.
- Within one specimen the weights are unit weights. Within one species each specimen gets equal weight.
- **α.** One α per gene per coordinate run.
  - The package coordinate and all its refits (fold, half, equal depth, placebo) reuse 71's condition-split α.
  - Package DPT and each M arm get their own α under their species-split M_full.
  - The half-B reuse is also exact in principle: binomial thinning of an NB count keeps its α.
- **Failures.** In F-tests, non-converged genes are median-filled (71's rule). In descriptors they are NaN. Separated
  genes are flagged and excluded from descriptors.

**The positional F-test pattern** (71) for every new gene test.
- F = (LR_num / df₁) / s²_post, where s² = LR_den / df₂ goes through `moderated_f` (squeezeVar, covariate log mean
  count, robust FALSE).
- z = Φ⁻¹(1 − p); BH over 71's universe.
- Set test: `camera_pr(z, 0.01)`; a call is Up with FDR ≤ 0.05; set z = Φ⁻¹(1 − p_one_sided).
- 71's pass rule: NCDR < 0.25, each relabeling < 5% of tested sets, calibration share ≤ 0.075.
- **The denominator always mirrors the numerator at specimen level.** If the numerator adds group × X beyond group × B,
  the denominator adds within-group specimen × X beyond specimen × B. This is how 71 judges T_spatial against
  LR_spec.

**Coordinates.**

| Code | Coordinate | Role |
|---|---|---|
| P | `pseudospace` on 71's scale | primary |
| D | `dpt_pseudotime` (finite units), `support_positions` | sensitivity (D3) |
| E | 72's equal-depth refit, full counts | the depth leg of the extrema rule only |
| H | 72's half-A refit, half-B counts | the count-split arm |
| PL | placebo refit (R6) | claim gates |
| M | SCF13, DPT13, LCP, LRS-Z on 72's matched set M (separable, O1) | the range of noise in Fig 3d |

- Each coordinate uses its own common support and its own 101-point `common_grid`.
- Genes and sets are always 71's universe.

**Segments.**
- L: the package labels.
- C: coordinate-threshold segments (`coordinate_segments` per specimen and coordinate).
- Segment units: logistic transitions per specimen (58's rule).

**Fold protection.** Excluded: the base integration features and the 15 panel and anchor genes.
- Package limit (72's O4): the labels still use all genes.
- M arms use their own construction genes: SCF13 and DPT13, the 611 PCA genes; LCP, its 920 landmarks; LRS-Z, its
  landmark union.

**Step nulls.** Count rows are permuted together with their library inside specimen × C (arm A) or specimen × L
(arm B). Positions and labels stay with the units.
- 10 replicates per arm.
- Seeds: 57's 5701+ and 5801+ for pathway tests; 58's 5901+ and 6001+ for descriptors.
- α is held fixed. Every numerator and denominator is refitted.

### 3 · Results of record and their NB versions

Old values are SCF13, Gaussian, notebook 13 units (docs `pt-pathway-continuous-vs-discrete.md`).

**R1 · Beyond-steps decomposition** (57; Fig 3e). *Old:*
- Test A (label steps) 79 calls (6, 0); arm A median 40; ratio 0.51, giving "discordance".
- Test B (coordinate steps) 67 (1, 0) against arm medians of 9–11.
- 48 cross-coordinate within-segment candidates.

*NB:*
- **Every step design here includes shared steps** (1 × the same segments, **new**), so that each group carries its
  own full step pattern in the reduced model.
  - 57's `contrast_designs` fitted steps for the case group only. That left the reference group's steps to the shared
    spline, and freeing that spline with group × spline could then score steps as shape.
  - This is a plausible cause of 57's floor: 3–31 calls in exact nulls that should call none.
  - 57's own design runs as a sensitivity on P (observed only).
- Test A numerator: LR(shared spline + shared L-steps + group × L-steps → + group × spline), df₁ 6.
- Denominator: LR(that model + specimen × L-steps → + specimen × spline), df₂ 12 (relabeling 6).
- Test B is the same with C.
- Relabelings carry the species shape as nuisance in the test's own basis: species × spline (57) plus species × the
  same steps (75's convention). This deviates from 57, which had only species × spline; without the steps, a real
  species step would enter the relabeled denominator.
- The df were checked on synthetic designs and are asserted again at run time (G6).
- Run on P and D, all three partitions. Arms on P and D, species split only.
- Unchanged:
  - the primary metric (arm-A median test-A calls ÷ observed);
  - the decision (≥ 0.5 discordance, < 0.2 rejected, otherwise partly explained);
  - the reported floor.
- **Within-segment candidate:** called by A and B on both P and D.
- **Strict candidate:** also, both set z exceed every simulated set z.

**R2 · Continuous against discrete calls** (61 part 2; Fig 3e). *Old:* steps 50, continuous 60, shared 47; continuous
adds 13 (step q 0.058–0.136; 7 called by no specific discrete method); steps add 3.

*NB:*
- A = 71's calls. S = 75's segment F-test calls; no step model is fitted here.
- Reported in sets and in 71's redundancy clusters: |A ∩ S|, |A∖S| and |S∖A|.
- For each set in A∖S:
  - S's FDR ("near miss" if ≤ 0.25, a descriptive label);
  - whether test B calls it;
  - whether 75's segment ORA calls it;
  - probe sensitivity: cameraPR on 71's z restricted to the balanced-probe genes gives FDR > 0.10, or the set's
    mean-z contrast falls below 0.7 of the full one (37's rule).
- **Ladder rows (74's, per 75's O2).** Both have the shared spline and shared steps in every model; relabelings are
  as in R1; the denominators mirror the numerators.
  - "Steps on the spline": group level → + group × L-steps; df₁ 2; df₂ 4 (relabeling 2).
  - "Six coordinate bins": the same with 6 equal-width bins of P; df₁ 5, df₂ 10 (relabeling 5).
  - Old values: 16 (0, 0) and 49 (17, 0).
  - Saved in 75's layout, so the method table can merge them.
- **Categories of A for Fig 3e:**
  - also called by S;
  - continuous-only;
  - beyond coordinate steps (test B);
  - sign change: the members' mean Δ reaches ≥ 0.25·max|Δ| on both signs (61's E4 rule, natural-log scale);
  - location: the segment of max |Δ|.

**R3 · Worked examples** (61 E1–E5; Fig 3a). These are 61's rules with NB statistics.
- Continuous calls are 71's. Step calls are 75's segment F-test. Test B is R1's.
- Whole-PT ORA, segment ORA and segment GSEA are 75's DESeq2 rows. 61 chose its examples by these same rules.
- **Lead gene:** the member with the largest gene z of the rule's statistic:
  - 71's z for E1 and E4; in E4 the member's own Δ must also reverse;
  - test-B z for E2;
  - segment-F z for E5;
  - |DESeq2 Wald| for E3.
- **Each lead gene is drawn three ways:**
  - whole-PT specimen pseudobulk ln rate;
  - specimen × segment pseudobulk;
  - 71's NB species curves, with specimen × 12-bin pseudobulk points.
- The pathway display (members' mean of per-gene centred NB curves) goes to the supplement.
- All of this is display only: no test.

**R4 · Descriptors per specimen** (58; Fig 3c).
- **Curves.** Reuse 71's NB specimen curves on P. D, E, H, PL and M get new unit-weighted NB fits with the same basis.
- **Scale.** Thresholds are read on the natural-log rate scale, as the frozen NB note requires: amplitude 0.5, range
  0.1, prominence 0.2, inner 0.2, |ρ| 0.9.
  - Pre-specified sensitivity: log1p(10⁴ × rate), the old display scale. The ln scale inflates the amplitudes of
    near-absent genes. The nulls use the same scale, so the comparison stays fair.
- **Segment means for r.** Specimen × L pseudobulk ln(max(Σy, 0.5) / Σlibrary).
- **(iv) Onset order** in mouse S2, given r.
  - Unchanged: the gene set; the partial Spearman with 2,000 bootstraps (seed 58); claim if ρ ≥ 0.3 with its interval
    above 0 and the count-split ρ ≥ 0.2; abandon if both are < 0.1; the step-null comparison.
  - *Old:* 0.36 (0.16–0.52); null maximum 0.49; inconclusive.
- **(i) Within-segment slopes.**
  - The NB slope is the signed √LR of (1 → 1 + position) per specimen × segment, with α fixed, at least 30 units, and
    converged, non-separated genes. It replaces the OLS t.
  - Agreement is the Spearman between the two specimens of a species over fold-protected genes; (i-b) uses C.
  - *Old:* mouse S2 0.36 and human S1 0.44 exceed the step nulls; mouse S2 meets the claim rule (count split 0.20).
- **(ii) Interior extrema.** 58's rule on the NB curves.
  - Cross-coordinate rule: ≥ 2 of {P, D, E}, including E.
  - *Old:* human S1 is the only cell beyond the nulls on every coordinate.
- **(iii) Boundary shifts.** 58's rule.
  - Registration threshold: the 95th percentile, over 72's fold and half refits, of |Δgap| ÷ the human segment length.
  - Programs are 75's `programs_human_vs_mouse.csv`.
  - Claim if ≥ 10 genes or ≥ 3 programs, with the same sign on D.
  - *Old:* 5 genes, not claimed.

**R5 · External onsets** (59). *Old:* partial ρ 0.10 (0.01–0.24), not replicated. After the post hoc step removal, the
slope agreement was 0.06–0.19.
- The external pseudotime is unchanged. It is an embedding, not a gene model (72's O3). Its protected genes are R4's
  (iv) genes ∪ the orientation markers.
- The external curves and slopes are NB, offset log(nucleus total without mt). α is estimated per donor under that
  donor's spline design.
- The rules are unchanged. The slope-agreement section is now pre-specified (secondary).

**R6 · Resolution budget, placebo and WSR** (60; Fig 3d). *Old:*
- SCF13 noise ÷ spread 0.83–1.49; 2 cells beyond the nulls (LCP: 4);
- noise against slope agreement, Spearman −0.31 over 24 cells.

**Noise side.** 72's `noise_budget.csv`, context `hvm: all PT units`, arms package and package DPT. The M arms come
from the matched context, if O1 holds.

**Information side.** R1 and R4 for every coordinate. The figure follows 60's layout.

**Placebo (new; 72's section 8).**
1. Rank the eligible genes by the within-species positional z (R7's W1), maximum over species. Eligible means: in 71's
   universe, outside the integration features, the panels and the anchors, and converged.
2. Within each decile of log mean count, take the lowest, as many as the base integration features in that decile
   (53's expression-matched rule).
3. Refit stages 03–04 once, with every other gene flagged `exclude_from_harmony_hvg` (`package_refits`). Guard: the
   refit's integration features ⊂ the placebo set.
4. The placebo's S1→S3 order per specimen is reported and never required.

**WSR.** Notebook 53's statistic, unchanged and log-normalized; held-out residual transfer is in P4c's class (72's O3).
- Folds: 72's fold refits for P and D. For PL, the single refit, with fold-k genes minus the placebo genes.
- 500 within-segment permutations, seed 61.

**Two claims per species × segment cell.** Each component is reported separately.

| Claim | Requires |
|---|---|
| "replicated within-segment gradients" | C14 met; cell not noise-dominated (mean fold noise ÷ spread < 1 over the species' two specimens; **new**); (i) ρ ≥ 0.3 and count-split ρ ≥ 0.2; ρ above every step-null replicate; ρ above the placebo's ρ (phase-3 §5) |
| "within-segment information in held-out genes" | C14 met; cell not noise-dominated; WSR p ≤ 0.01 in both directions; above the placebo in both directions; the placebo itself not significant (53) |

**If C14 is not met (≥ 6 of 12 cells noise-dominated).** Following 72's rule:
- R1 test B, R4 (i), (ii) and (iv), R7, and the Fig 3e "beyond coordinate steps" category become descriptive. No
  within-segment claim is made.
- R2's "continuous adds N" is reported as "calls the segment screen misses", not as within-segment biology.
- Fig 3c and 3d are drawn, captioned "within-segment information not claimed".

**R7 · Solute-transporter micro-gradients** (new; figure plan Fig 3b, gap 2).

*Question.* Do PT solute transporters carry within-segment gradients that replicate between specimens beyond segment
steps?

*Prior.* SLC proteins carry most of the ion and organic-solute transport along the tubule, as the curated SLC map
along the nephron shows ([PMID 34191628](https://pubmed.ncbi.nlm.nih.gov/34191628/)). The reference for
transporter expression by segment is discrete: 14 microdissected rat tubule segments, beginning with the PT
([PMID 25817355](https://pubmed.ncbi.nlm.nih.gov/25817355/)). Live imaging shows axial heterogeneity of function along
the mouse PT ([PMID 41803103](https://pubmed.ncbi.nlm.nih.gov/41803103/)).
*Abandon:* see "not worth showing".

*Gene set T.*
- HGNC gene-group branches "Solute carriers" (752) and "ATP binding cassette transporters" (417), protein-coding,
  downloaded 2026-10-12 (HGNC: [PMID 36243972](https://pubmed.ncbi.nlm.nih.gov/36243972/)).
  - Files `hgnc_group_752.tsv` (474 rows) and `hgnc_group_417.tsv` (53 rows) from
    `genenames.org/cgi-bin/genegroup/download?id=<id>&type=branch`.
  - sha256 752: `131e82911a0840e98f70c6c3075059d0a69a8709c0a3c43db69029aed3a97bee`.
  - sha256 417: `cc6cdac305cd75465de8ae34faac429a0356db916171bfeba4ec8d4ce82cf28a`.
  - Branch membership is used, not a symbol pattern: the SLC branch also holds MFSD, MPC, CTNS and others.
- Mapped by `human_symbol` into 71's universe: 341 genes (309 SLC, of which 41 are SLC25 mitochondrial carriers;
  32 ABC).
- SLC25 (mitochondrial carriers) is its own family in the figure. Every T result is reported with and without it.

*Statistics.* NB on all universe genes, per species, units of that species.

| Test | Numerator | Denominator | df₁ / df₂ |
|---|---|---|---|
| W1 positional | 1 + specimen contrast → + spline | specimen contrast × spline | 6 / 6 |
| W2 beyond L-steps | 1 + specimen contrast + L-steps → + spline | + specimen × steps → + specimen × spline | 6 / 6 |
| W3 beyond C-steps | as W2 with C | as W2 with C | 6 / 6 |

- Moderation and BH run over the universe within a species.
- **Micro-gradient gene:** q ≤ 0.05 in both W2 and W3.
- **Localisation:** a cell where both specimens' slope statistics share a sign with |·| ≥ 2.

*Nulls.*
- W3 under arm A and W2 under arm B are exact nulls (R4's replicates). Calibrated if the pooled share of p ≤ 0.05 is
  ≤ 0.075.
- The T-restricted slope agreement per cell is compared with the same replicates.
- cameraPR of T on W3 z per species is descriptive.

*Figure 3b.*
- Rows: transporters with W1 q ≤ 0.05 in either species; at most 60, chosen by the larger W1 z.
- Values: 71's NB species curves on P, z-scored within species; mouse block, then human block; ordered by mouse peak,
  then human peak.
- Specimen strips; transitions dotted; family annotated.
- Marks show micro-gradient calls and their cells. Noise-dominated cells are hatched.

*Not worth showing* (then 3b becomes a segment-level zonation heatmap without marks, or moves to the supplement) unless
all of these hold:
- C14 met;
- W2 and W3 calibrated;
- ≥ 10 micro-gradient transporters in at least one species;
- the T-restricted slope agreement beyond every step-null replicate in ≥ 1 non-noise-dominated cell.

### 4 · Outputs

Every table in `tables/` has a matching cached stage.

| Result | Tables |
|---|---|
| R1 | `beyond_steps_genes_{P,D}.csv`, `beyond_steps_sets_{P,D}.csv`, `beyond_steps_arms.csv`, `beyond_steps_decision.csv`, `coordinate_thresholds.csv`, `within_segment_candidates.csv` |
| R2 | `continuous_vs_segment.csv`, `continuous_adds.csv`, `segment_adds.csv`, `ladder_rows_{ds}.csv` (75's layout), `fig3e_categories.csv` |
| R3 | `worked_examples.csv` |
| R4 | `descriptors_summary.csv`, `onsets_iv.csv`, `slopes_i.csv`, `extrema_ii.csv`, `shifts_iii.csv`, `step_nulls.csv`, `scale_sensitivity.csv` |
| R5 | `external_donors.csv`, `external_onsets.csv`, `external_slopes.csv` |
| R6 | `placebo_genes.csv`, `placebo_refit.csv`, `wsr.csv`, `resolution_budget.csv`, `cell_claims.csv` |
| R7 | `transporters.csv` (T with W1–W3 and localisation), `transporter_calibration.csv`, `transporter_summary.csv` |
| All | `compare_old.csv`, `material_change.csv` |

Figures 3a–3e and S-panels. Code: a new module, `pseudospace/within_segment_nb.py`, with synthetic tests, holding:
- the beyond-steps, ladder and within-species designs with their mirrored denominators;
- the NB segment slopes;
- the permutation of counts together with their library.

It imports 75's `discrete_comparators` and changes no existing module.

### 5 · Material change from the old record

These rules are fixed now. Old and new differ in units, labels, coordinate and model, so every comparison is
descriptive. A rule that is met means the paper's text must change.

| # | Material change if |
|---|---|
| MC1 | R1's decision category differs from "discordance" (old ratio 0.51) |
| MC2 | test B calls ≤ 2 × the larger arm median, or test B fails the pass rule (old 67 vs 9–11) |
| MC3 | A∖S is 0 redundancy clusters, or \|S∖A\| > \|A∖S\| in clusters, or 75's segment F-test fails its pass rule (old adds 13, steps add 3) |
| MC4 | no cell meets the "replicated within-segment gradients" claim, or a cell other than mouse S2 meets it |
| MC5 | (iv) becomes claim or abandon (old inconclusive) |
| MC6 | (ii): a mouse cell exceeds the nulls on P and E, or human S1 does not |
| MC7 | (iii) or R5 reaches its claim |
| MC8 | the noise against slope-agreement Spearman over cells is ≥ 0 (old −0.31) |
| MC9 | any claim is lost only through C14 or the placebo gate (both new) |

R7 is new and has no old value. MC1 and MC2 are judged on the primary (shared-step) design. 57's own design, run beside
it, shows how much of any change comes from the design and how much from the model.

### 6 · Compute

24 jobs, from 71's timings: a specimen-shape LR takes about 15–20 s and a gene model about 35 s.

| Step | Time |
|---|---|
| R1 observed on P and D (24 nested fits) | about 10 min |
| R1 arms on P and D (2 × 20 replicates × 4 fits) | about 35 min |
| R2 ladder, both datasets | about 5 min |
| R4 on P, D and E, with step nulls on P and D | about 35 min |
| R4 count split | about 3 min |
| R7 W1–W3 and their nulls | about 12 min |
| Placebo refit (one package run; ≥ 40 GB free, about 7 GB transient; needs `PSEUDOSPACE_VISIUM_ROOT`) | about 5 min |
| WSR | about 15 min |
| R5 | about 15 min |
| M arms (O1) | about 45 min |

Total: about 2¼ h, or 3 h with the M arms.

### 7 · Open decisions (coordinator)

- **O1.** Run the M comparators' information side (separable, default on as in 72's O1)?
- **O2.** May 74 run one package refit for the placebo?
- **O3.** The new per-cell noise gate, added to 72's global C14.
- **O4.** The primary threshold scale: ln, as the NB note requires, with log1p as the sensitivity.
- **O5.** T = all of SLC ∪ ABC, or without SLC25.
- **O6.** Store the HGNC files in the data root with their sha256.
- **O7.** AKI: only the R2 ladder rows.
- **O8.** Approve the shared-step correction of 57's designs (R1, R2 ladder), with 57's design kept as a sensitivity.

### 8 · Coordinator decisions (2026-10-12)

- **O1. Yes.** The M comparators (SCF13, DPT13, LCP, LRS-Z on 72's matched set) run as a separable section, reading
  only 72's cached positions. Dropping it changes no other table, figure or rule. The paper's framing of those methods
  is still open with the user.
- **O2. Yes.** One placebo refit: the unmodified package, CLI stages 03–04, written to
  `results/pt_package_continuous_vs_discrete/refits/placebo/`.
- **O3. Approved.** The per-cell noise gate, as pre-specified.
- **O4. Approved.** ln is the primary threshold scale; log1p is the pre-specified sensitivity.
- **O5.** T is the full HGNC SLC branch (no hand pruning) plus the ABC branch. SLC25 is drawn as its own family, and
  every T result is reported with and without it.
- **O6. Approved.** The HGNC files go under `data/external/hgnc/`. The notebook asserts their sha256 and records the
  download date, 2026-10-12.
- **O7. Approved.** AKI runs only the R2 ladder rows.
- **O8. Approved.** Shared steps in every step model. 57's case-only steps run as a sensitivity.
- **Order.** Notebooks 72 and 75 are still running. Phase 2 runs first the sections that need neither: inputs and G1–G2,
  R1, the R2 ladder rows, R4 on P and D with the step nulls, R5, the R7 tests and the placebo refit.
  - Next come the sections that need only 72's refit payloads, once those exist: R4 on E, and WSR, which uses 72's
    fold refits.
  - Last, behind G3 and G4, come the sections that read 75's tables (R2 comparisons, R3, R4 (iii)) or 72's tables (the
    count split, R6 noise and the claim tables).

## Notebook 74 results

**Run.**
- Notebook: `analysis/notebooks/74_pt_package_continuous_vs_discrete.ipynb`, logic version `74.package.1`.
- Output: `results/pt_package_continuous_vs_discrete/`.
- Compute: 77 min uncached at 24 workers. Most of it was the R1 arms, R4 step nulls and comparator nulls; the placebo
  refit took about 2.5 min. A cached rerun takes about 6 min.
- New module: `pseudospace/within_segment_nb.py`, with `tests/test_within_segment_nb.py` (8 tests). The full suite
  passes (520).
- Data: the HGNC files are in `data/external/hgnc/` (sha256 asserted, downloaded 2026-10-12, `SOURCE.txt`).

**Guards.**
- The confirmed pass-1 fingerprints are in place (human vs mouse `efafc9ac2e99`, AKI `bb4a80c1b484`), and the package
  is 0.1.0 at 59da04b.
- G1: the builder reproduces 71's counts digest, universe (11,055) and sets (1,513); the AKI analogues pass too.
- G2: 71's T_spatial is refitted to a relative difference of 7e-6.
- G3: 75's segment tables cover 71's universe and sets.
- G4: 72's C14 is recomputed and agrees.
- G5: 72's half-B counts are reproduced exactly; the stage key of 72's G6 arm b matches.
- df are asserted from the design ranks everywhere: R1 (6, 12/6), ladder (2, 4/2 and 5, 10/5), W tests (6, 6).

**72's verdict sets the frame: C14 is not met.** 10 of 12 cells are noise-dominated. The species-mean fold noise ÷
spread is 0.98–4.25 (old SCF13: 0.83–1.49), and only mouse S3 is below 1. By the protocol, every within-segment
result below is descriptive: no within-segment claim is made. R1 test B, R4 (i), (ii) and (iv), R7, and the Fig 3e
"beyond coordinate steps" category are all reported without a claim.

### How the protocol was read (also in `run_record.json`)

1. **Transitions out of order.** On P, D and E, HUK1_MED1's logistic S1→S2 transition lies after its S2→S3
   transition (0.67 against 0.59 on P). This is 72's C1 failure. On the count-split coordinate both human sections are
   out of order.
   - 58's rule for such a specimen was applied: no segment units, so no onsets, extrema or boundary shifts there. Its
     slopes are kept.
   - Human (ii), human (iv) and all of R4 (iii) are therefore **not evaluable**.
2. **Placebo refit.** The package keeps only variable genes as integration features (Seurat dispersion ≥ 0.5 in each
   species). Non-zonal placebo genes are not variable: only 45 qualified, and stage 03 stopped (50 PCs need ≥ 50
   features).
   - The closest reading of notebook 53's placebo, which used its pipeline with the placebo genes as features, sets
     `integration.hvg_min_disp` to −10⁹ for this refit only. This went through a copy of the tracked config
     (`refits/placebo_config/`). The package code and every other parameter are unchanged.
   - 977 of the 988 placebo genes became features. The refit kept 6,884 of 11,947 PT units in pass 2.
   - The placebo does not order S1→S3 (Spearman with P 0.10). It was flipped, which changes neither slope agreement
     nor WSR.
   - The protocol's decile rule had a shortfall in the lowest-expression decile: 109 features had no eligible
     universe gene.
3. **WSR folds and units.**
   - 72's human vs mouse fold-2 refit failed its orientation guard, so WSR uses folds 0 and 1 for every arm.
   - It runs on the 6,884 units with a position in every arm. The placebo dropped the other units.
4. **R3, E3's lead gene.** It uses 75's segment-F z, because 75 stores DESeq2 Wald statistics per scope.
5. **R4 (iii) programs.** 75's programs carry no driver column, so drivers would be the members in the top 10% of 71's
   z (61's rule). This is moot, because (iii) is not evaluable.
6. **The remaining readings** are minor and listed in `run_record.json`:
   - (iv) nulls use the point estimate;
   - the transporter slope agreement needs > 20 genes;
   - MC8 uses the M cells, where all six arms share one unit set;
   - the within-species tests use genes with counts in both specimens.

### R1 · Beyond-steps decomposition (57)

| | Old (SCF13, Gaussian) | P (package) | D (package DPT) |
|---|---|---|---|
| Test A, beyond label steps (relabeled) | 79 (6, 0) | **64 (0, 0)** | 66 (0, 0) |
| Arm A median test A (range); ratio | 40; **0.51** | 35 (21–46); **0.55** | 44 (30–56); 0.67 |
| Test B, beyond coordinate steps (relabeled) | 67 (1, 0) | **74 (0, 0)** | 65 (0, 0) |
| Arm medians, test B (arm A, arm B) | 9, 11 | 0, 7 | 0, 24.5 |
| Floor (arm B test A, arm A test B) | 3–31 | **0, 0** | 0, 0 |
| 57's case-only design, tests A, B | – | 83, 86 | – |
| Label agreement of coordinate segments, mouse / human | 0.92 / 0.85 | 0.95 / 0.87 | 0.95 / 0.87 |

- **Decision: label–coordinate discordance**, as before. About half of the beyond-label-steps signal is produced by
  coordinate steps plus label noise.
- **Test B** (74 calls) is far above both exact nulls, and the floor is gone with shared steps. 57's case-only design
  gives 83 and 86 calls.
- **All four tests pass 71's rule** on both coordinates (calibration 0.026–0.040).
- **Within-segment candidates** (A and B on P and D): 49. Of them, 28 are strict, with P set z above every simulated z.
  47 of the 49 are 71 calls.
  - They include fatty acid and amino-acid metabolism, SLC-mediated transport, the inorganic cation/anion and organic
    cation/anion transport sets, xenobiotic metabolism and glutathione conjugation.
  - They are descriptive, because C14 is not met.

### R2 · Continuous against discrete calls; ladder rows (61, 37)

- **A (71) against S (75's segment F-test)**, both passing their rules:
  - **71 adds 12 sets (12 redundancy clusters); the steps add 5.** Old: 13 and 3.
  - All 10 of 71's additions with segment FDR ≤ 0.25 are near misses (0.052–0.225). Two are not:
    glycolysis/gluconeogenesis (0.29) and the TCA cycle (0.45).
  - 10 of the 12 are also called by test B.
  - 11 of the 12 are called by no 75 segment ORA. Two are probe-sensitive (VLCFA β-oxidation and organic
    cation/anion/zwitterion transport).
  - The steps add drug metabolism, pantothenate/CoA, cellular hexose transport, Na⁺/Cl⁻ neurotransmitter transporters
    and surfactant metabolism, each with 71 FDR 0.054–0.072.
- **Ladder rows (NB, shared steps, mirrored denominators):**

  | Row | Human vs mouse | AKI |
  |---|---|---|
  | Steps on the spline | 75 (69 clusters); 0, 0; calibration 0.049; pass | 53; 2, 7; 0.088; fail |
  | Six bins | 101 (93); 1, 1; 0.048; pass | 97; 9, 13; 0.103; fail |

  Old (human vs mouse, Gaussian): 16 (0, 0) and 49 (17, 0). The AKI rows fail calibration, as 71's AKI positional
  F-test does.
- **Fig 3e categories of 71's 83 calls** (sets; redundancy clusters):

  | Category | Sets | Clusters |
  |---|---|---|
  | Also called by the segment F-test | 71 | 66 |
  | Continuous-only | 12 | 12 |
  | Beyond coordinate steps (test B; descriptive) | 68 | 65 |
  | Sign change along the PT | 21 | 20 |

  - Location of max |Δ|: S1 33 (22 human-higher), S2 10 (all mouse-higher), S3 40 (31 human-higher).

### R3 · Worked examples (61's rules, NB; Fig 3a)

| Rule | Pathway (sets meeting the rule) | Lead gene |
|---|---|---|
| E1 all agree | Hallmark Fatty Acid Metabolism (21) | Hao2 |
| E2 segment F-test FDR > 0.05, curves structured | Phenylalanine and tyrosine metabolism (8; segment FDR 0.052) | Hpd |
| E3 discrete-only | KEGG Oxidative phosphorylation (133) | Cox7a1 |
| E4 sign change | Valine, leucine and isoleucine degradation (21) | Acadm |
| E5 step-like | Hallmark Bile Acid Metabolism (13) | Gc |

### R4 · Descriptors per specimen (58), on P with 10 + 10 step nulls

- **(iv) Onset order in mouse S2, given r: 0.11 (−0.02 to 0.24)**, 401 genes.
  - Inside the step nulls (maximum 0.34).
  - Count split 0.39 (0.19–0.56), 140 genes.
  - **Decision: inconclusive** (old 0.36, inconclusive).
  - D gives 0.21 and E 0.14, also inside their nulls.
- **(i) Within-segment slope agreement** (fold-protected; descriptive):

  | Cell | Slope ρ (P) | Step-null max | Count split | Placebo |
  |---|---|---|---|---|
  | Mouse S1 | 0.38 | 0.30 | 0.16 | −0.01 |
  | Mouse S2 | 0.37 | 0.21 | 0.16 | −0.03 |
  | Mouse S3 | 0.49 | 0.14 | **0.38** | 0.00 |
  | Human S1 | 0.46 | 0.14 | 0.07 | 0.07 |
  | Human S2 | 0.09 | 0.12 | 0.01 | 0.01 |
  | Human S3 | 0.18 | 0.12 | 0.10 | 0.03 |

  - Five cells exceed every step-null replicate; human S2 does not (old: mouse S2 and human S1).
  - (i-b), inside coordinate segments: mouse S1–S3 and human S1 exceed both nulls.
  - The placebo carries no slope agreement anywhere (−0.03 to 0.07).
  - Only mouse S3 also meets ρ ≥ 0.3 with count split ≥ 0.2.
- **(ii) Interior extrema:** only mouse S3 exceeds the step nulls (333 against a maximum of 309). Mouse S1 and S2 fall
  inside them. The human pair is not evaluable (reading 1). Old: human S1 only.
- **(iii) Boundary shifts:** not evaluable (reading 1). Old: 5 genes, not claimed.
- **Scale sensitivity** (log1p, `scale_sensitivity.csv`):
  - The slopes are scale-free and unchanged.
  - The mouse S3 extrema collapse from 333 to 7 and no longer exceed the nulls, so the (ii) result depends on reading
    amplitudes on the natural-log scale, which inflates near-absent genes.
  - (iv) rises to 0.41 (160 genes) and is still inside its nulls.

### R5 · External mice (59; Census, 6 male adult donors, all ordered S1 < S2 < S3)

- **Onsets: partial ρ 0.30 (donor bootstrap −0.00 to 0.43)**, 61 genes.
  - **Decision: not replicated at the pre-specified level.** The estimate clears 0.2, but the interval touches 0.
  - The raw Spearman of our onsets against external onsets is −0.21.
  - D gives 0.33 (0.21–0.46); E 0.18.
  - log1p sensitivity: 0.02 on P.
  - Old: 0.10 (0.01–0.24).
- **Within-segment slopes against the external donor median** (secondary):
  - raw ρ 0.37–0.41;
  - after removing the external S2−S1 and S3−S2 steps, partial ρ 0.17 (S1), 0.20 (S2) and 0.24 (S3) on P.
  - Old: 0.06–0.19.

### R6 · Resolution budget, placebo, WSR and the claim gates (60)

- **Noise against information** (Fig 3d; 72's fold noise):
  - over the 36 M cells (six arms on one unit set), Spearman of noise ÷ spread against slope agreement **−0.24**
    (old −0.31);
  - LCP and LRS-Z have the lowest noise on M (0.13–0.47).
- **WSR**, P against the placebo:
  - **P:** p ≤ 0.01 in both directions in mouse S1 and S3 and human S1–S3; mouse S2 misses at 0.0499.
  - **Placebo:** never significant (smallest p 0.046), and P exceeds it in every cell and direction.
- **Claims** (`cell_claims.csv`):
  - **No cell meets either claim**, because C14 is not met.
  - Without C14 and the per-cell noise gate, "replicated within-segment gradients" would hold in mouse S3, which is
    also the only cell that is not noise-dominated.
  - Without C14 and the noise gate, the held-out-genes claim would hold in mouse S1 and S3 and human S1–S3.

### R7 · Solute transporters (Fig 3b)

- **T:** 341 genes (44 SLC25, 32 ABC); 297 without SLC25.
- **Zonated (W1 q ≤ 0.05):** 247 in mouse, 147 in human.
- **Micro-gradient genes** (W2 and W3 q ≤ 0.05): mouse 115, human 62 (without SLC25: 103, 59). Across all universe
  genes: 1,410 and 905.
- **Calibration:** both tests are calibrated on their exact nulls (pooled share of p ≤ 0.05: W3 0.035, W2 0.032).
- **Set test:** cameraPR of T on W3 z is Up in both species (mouse p 8e-8, human 1e-5; without SLC25 1e-8 and 1e-6).
- **T-restricted slope agreement:** 0.29–0.75; beyond every step null in five of six cells (all but human S2; without
  SLC25, all but human S3).
- **Gate:**
  - **Not met, because C14 is not met.** Its other three conditions hold: calibration, ≥ 10 genes, and mouse S3 beyond
    the nulls in a non-noise-dominated cell.
  - **Fig 3b is therefore the segment-level zonation heatmap**, with within-segment marks labelled descriptive and
    noise-dominated cells hatched.

### Material change from the old record (protocol section 5)

| # | Old | New | Met |
|---|---|---|---|
| MC1 decision | 0.51, discordance | 0.55, discordance | no |
| MC2 test B | 67 vs 9–11 | 74 vs 7; pass | no |
| MC3 adds | 13 / 3 | 12 / 5 clusters; segment F-test passes | no |
| MC4 gradients claim | mouse S2 | no cell (C14) | **yes** |
| MC5 (iv) | inconclusive | inconclusive (0.11) | no |
| MC6 (ii) | human S1 only | mouse S3 only; human not evaluable | **yes** |
| MC7 (iii), R5 | 5 genes; 0.10 | not evaluable; 0.30 with interval touching 0 | no |
| MC8 noise vs agreement | −0.31 | −0.24 | no |
| MC9 lost only through new gates | – | mouse S3 (gradients); held-out genes in mouse S1, S3 and human S1–S3 | **yes** |

**What changes for the paper text.**
- **Still stands:**
  - "continuous ≈ fine discretisation in calls; adds a few near-miss pathways" (12 against 5);
  - the beyond-label-steps signal is half discordance;
  - test B is real and step-specific.
- **Becomes descriptive:**
  - the within-segment-gradient statements (Fig 3c and 3d);
  - the transporter micro-gradients (Fig 3b).

  They replicate between specimens beyond step nulls and a placebo, and partly in external nuclei. But the package
  coordinate's fold noise exceeds the within-segment spread in 10 of 12 cells, so the frozen rule forbids a
  within-segment claim.
