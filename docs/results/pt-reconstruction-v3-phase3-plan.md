# Workstream A (ws7), phase 3 plan: improving the reconstruction method

Status: **plan for review. No private-data fit has run for phase 3.** The only new computation is a
descriptive re-read of phase-2 posteriors already cached by notebook 51 (§3.2).

**Constraints (user).**
- Kidney data only.
- No x/y anywhere: the method must stay general for any continuous repeated structure. In kidney,
  neighbouring sections are often different nephrons.
- Improve the method.

**Phase-2 starting point** (notebook 51, committed 575f5aa).
- *LCP was not adopted:* abandonment rule (iii) fired, with P4c > 0 in 3 of 6 mouse cells.
- *Landmarks fixed within-segment stability.* Mouse gene-fold agreement: SCF13 0.60, DPT13 0.61,
  LRS 0.91, LCP 0.96.
- *The zero-cost LRS got most of the gain.*
- *No coordinate beat segment + coverage on P4c* (≤ 0 in 44 of 48 cells).
- *Workstream B, notebook 58 (i) on LCP:* within-segment slopes replicate between the two mice
  (ρ 0.37–0.44), against a step-null maximum of 0.10–0.14.

---

## 1. What a model must add beyond LRS, and how each part is tested

LRS is log((late + 0.5)/(early + 0.5)), summed over landmark counts. Its expected value does not
depend on exposure. It is per-structure, needs no fitting, and gives no uncertainty. A model earns its
complexity only through what LRS structurally cannot do. Each candidate addition below gets one
pre-registered ablation on the existing evaluation machinery. Margins are fixed now.

| Component | First-principles reason it could beat LRS | Ablation (arms differ only in this) | Keep if | Abandon (drop the component) if |
|---|---|---|---|---|
| **C-lik:** per-gene NB likelihood with exposure offset | LRS weights genes by raw count and uses only a sign. A likelihood weights each gene by its information (contrast² × rate / dispersion), so genes with high counts but weak contrast no longer dominate | LRS vs M0. M0 is an NB likelihood with log-linear paths (2 parameters per gene, sign from the landmark direction), pooled M-step, same landmarks | Mouse within-segment gene-fold agreement +0.03, **and** P1w (evaluation references) not lower by > 0.02, **and** human S1 gene-fold not lower | M0 ≤ LRS + 0.03 on gene-fold, or worse on P1w by > 0.02 |
| **C-shape:** non-monotone zone shapes | A sign carries no information on genes that peak in S2. LRS cannot use the S2-high genes the data-derived zones provide | M0 vs M1 (spline paths, df 6) | Gene-fold +0.02 **or** WSR (§4) +10%, with nothing else worse by > 0.02 | Neither |
| **C-bal:** shared, specimen-balanced path | LRS has no shared scale. Rank-mapped LRS moves when a specimen's composition changes (for example, cortex-only human). A shared path fixes the gauge | M1 pooled vs M1 balanced. **Composition stress:** delete 70% of one specimen's S3 rows (seed 52), refit, and measure the within-segment median \|Δz\| of untouched structures in other specimens, in segment units. LRS ranks are included | Balanced shift ≤ 50% of the pooled and LRS shifts | Shift not reduced by ≥ 50% |
| **C-unc:** calibrated posterior | LRS needs an extra device for uncertainty. The model's posterior is native, but only useful if calibrated | LRS-U (LRS plus a parametric NB bootstrap of counts plus a landmark-gene bootstrap, B = 200, giving an SD) vs the model posterior. Same calibration tests (§3.2) | Model calibration error (\|SD of standardized differences − 1\|) ≤ LRS-U's, with interval width ≤ LRS-U's | Worse on both |
| **C-reg:** species registration | LRS uses different landmark sets per species, so its scales are not comparable. A model can share shapes across species through conserved landmarks | Separate species fits vs a joint fit (shared shapes for conserved landmarks, species-specific level and amplitude). Test §3.3 | Cross-species agreement of conserved held-out genes' within-zone curves +0.05 | No gain |

**If no component passes,** the honest method is **"LRS with uncertainty" (LRS-U)**: data-derived
landmarks (§2), the ratio, bootstrap SDs, and a zone-unit gauge. That is an acceptable outcome, and
the paper would say so.

---

## 2. General landmark selection from the data's own zones (ZAL, zone-anchored landmarks)

**First principles.** Any tissue with a repeated continuous structure has coarse discrete zones:
clusters or labels. Zones give order at their own resolution. Genes that differ between zones are,
by construction, the genes that change along the axis, which is exactly what a landmark needs. The
same genes carry the within-zone signal if expression changes gradually. External atlases are then
needed only to **evaluate**. This frees all three phase-2 selection references (microdissection,
female mouse snRNA, human half A) as additional evaluation references.

This is also the paper's story: discrete clustering gives coarse order, and the continuous
reconstruction refines it inside zones.

### Procedure (per species)

1. **Zones.** The reviewed S1/S2/S3 labels stand in for any tissue's clusters; no marker list is
   used. Labels came from clustering on all genes, so zone-level information is shared with
   evaluation genes. Within-zone metrics, the target here, are not affected. Stated as a limitation.
2. **Zone order, inferred, not assumed.**
   - Take specimen-balanced zone-mean log-rate profiles on read half A.
   - Choose the shortest Hamiltonian path through the zone centroids (brute force for ≤ 6 zones).
   - Orientation needs one input, the start zone (S1). This is the only domain knowledge allowed.
   - Guard: the inferred path must be S1–S2–S3 in both species. If it is not, report and stop.
3. **Landmark selection on read half A** (Poisson split ε = 0.5, seed 53), excluding the scored gene
   fold.
   - Eligible genes: detection ≥ 10% in the species, plus target-class specificity (no non-PT class
     mean > 2× the PT mean), the generic contamination filter.
   - Statistic: quasi-Poisson zone-versus-constant likelihood-ratio F with a library offset and
     specimen-balanced weights.
   - Take the **top 500** genes by F per species and gene fold. 500 is fixed a priori, close to phase
     2's 420/556.
   - Each gene's shape = its zone means at zone centres.
4. **Gauge.** Zone j occupies [c_{j−1}, c_j] on [0, 1], with widths equal to the specimen-balanced
   occupancy fractions on half A. Positions are reported in **zone units**: S1 → [0, 1), S2 → [1, 2),
   S3 → [2, 3], workstream B's convention.
5. **Positions fitted on read half B** with half-A landmarks, and swapped (half-B landmarks, half-A
   positions). The **final posterior is the normalized product** of the two half-likelihoods times
   the prior. The halves are independent given position, so every read is used once for fitting, and
   never by the selection that chose its gene.
6. **Zone prior (soft).** p(z | label j) = (1 − ε)·U(zone j) + ε·U(other zones), with ε = 0.10 fixed.
   - Ablation M4: the same model with no label prior, shape only.
   - Labels never enter the within-zone order except through which genes were selected.
7. **Placebo control.** The identical pipeline with the 500 *lowest*-F eligible genes, matched in
   expression. Any within-segment metric must beat the placebo coordinate. Otherwise it is measuring
   co-expression programs, not position.

### Comparators

- LCP-ext and LRS-ext, with phase-2 external landmarks;
- LRS-zone, the ratio over zone landmarks;
- SCF13, DPT13, segment labels, and the segment + coverage oracle.

### Rules

- **Keep ZAL** if it is non-inferior to LCP-ext on mouse gene-fold (≥ −0.03) and on P1w measured
  with references neither arm selected on (male mouse snRNA, human half B), within 0.03.
- **Abandon ZAL** if gene-fold agreement drops by > 0.10 or P1w by > 0.05 against LCP-ext.
- If abandoned, external atlases are needed and the method is not general. That would be reported as a
  limitation.

---

## 3. Engineering fixes

### 3.1 Convergence (all 42 phase-2 fits hit the 100-iteration cap)

**Cause.** The likelihood is invariant to monotone re-parameterization of z with matching changes in
the rate paths. The λ = 10 shrinkage is negligible against counts in the thousands, so EM creeps
along a near-flat gauge direction. Ranks were stable (ρ = 1.000 across perturbed starts); the gauge
was not.

**Fixes.**
- Anchor the gauge with the zone prior (§2.6).
- For the shape-only ablation, add occupancy anchoring: after each E-step, monotone-remap the grid so
  that the specimen-balanced cumulative occupancy matches the zone-width reference.
- Use a finer grid, K = 200. Phase-2 posterior SDs in mouse (0.012–0.018) were below the grid step
  (0.020), so posteriors were quantized.
- Convergence: relative objective change < 10⁻⁹ **and** max |Δ posterior mean| < 10⁻⁴ for 5
  consecutive iterations; cap 500.

**Abandon this fix** if more than 10% of fits still hit the cap. Report, and do not raise the cap.

### 3.2 Posterior calibration

**Re-read of phase 2 (cached posteriors, no fit).** Phase 2's "66–83% coverage" was a conservative
check: half B's MAP is itself noisy. The proper statistic is
d = (z_A − z_B)/√(sd_A² + sd_B²), whose SD should be 1 with 80% of |d| ≤ 1.28.

| Pair | Mouse S1 / S2 / S3: SD of d (share \|d\| ≤ 1.28) | Human S1 / S2 / S3: SD of d (share \|d\| ≤ 1.28) |
|---|---|---|
| Read split | 0.96 (0.81) / 0.93 (0.81) / 1.33 (0.80) | 1.30 (0.84) / 0.97 (0.84) / 0.88 (0.86) |
| Gene folds 0–1 and 0–2 | 0.68–1.29 (0.71–0.88) | 0.67–2.34 (0.68–0.93) |

Read-noise calibration is roughly nominal, with heavy tails. Gene-sampling uncertainty is
under-covered, worst in human S1 (SD up to 2.3). Gene-fold pairs share reads, so they understate the
gap.

**Fix (fixed now).** A design-effect temperature T = Var_obs(d_gene)/Var_model(d_gene). It is
estimated from disjoint landmark *gene halves*, which use landmark genes only and never touch
evaluation genes, per species and zone, then frozen. The temperature is applied to the likelihood,
with read-split d reported as a check.

**Calibration targets.**
- SD of d in [0.85, 1.15] for both the read-split and the gene-half pairs;
- share of |d| ≤ 1.28 in [0.75, 0.85];
- reported per species × segment.

**Abandon calibration** if targets fail in more than half the cells. The posterior is then reported
as conditional and not calibrated, and LRS-U's bootstrap is used as the uncertainty.

### 3.3 Species gauge (registration gap 0.28) without x/y

**Cause.** Human external shapes had two zones (early at 1/3, S3 at 5/6) against three in mouse, so the
gauges differed by construction. The gap was constant across refits (0.27–0.28), so this was an
offset, not instability.

**Fix.** Zone-derived shapes give both species three zones, plus the zone-unit gauge (§2.4).
- *Caveat:* in zone units the median-based gap becomes nearly tautological, because zone medians sit
  near zone centres. It is still reported, flagged.
- *Non-tautological test (pre-registered).* Take conserved held-out genes: same sign in the male mouse
  snRNA and human half B evaluation references, and never landmarks. Fit within-zone spline curves
  per species, then report the cross-species Spearman of those curves within each zone, against a
  shift null (human positions shifted ±0.25 zone units).
- *Joint fit (C-reg ablation):* shared shapes for conserved landmarks, species-specific level and
  amplitude.

**Abandon the joint fit** if the curve agreement does not rise by ≥ 0.05 over separate fits.

### 3.4 Human S1 (4.5–7.7 UMI per bin; F5b ceiling 0.20)

No x/y, so the only levers are more informative genes and honest uncertainty.

1. **S1–S2 landmarks.** Human external references had none (convoluted PT against S3 only). Data
   zones provide S1-versus-S2 genes for human, which is the information the external design lacked
   inside early PT.
2. **Information-weighted likelihood** (C-lik) instead of LRS's count weighting.
3. **Calibrated, wide posteriors** (§3.2), so human S1 order is reported with intervals and
   downstream analyses average over them (notebook 27 precedent).

**Abandon the claim of resolved human S1 order** if, under the best surviving arm:
- human S1 gene-fold agreement < 0.50, **or** human S1 P5 < 0.70.

Human S1 is then stated to be unresolved at this depth.

---

## 4. Evaluation (frozen)

**Unchanged.**
- Notebook 51's full benchmark, per species × segment, for every arm, including the ones we lose:
  P1, P1w, gene-fold, P5, LOSO, count split, registration gap, P2, P4, and **P4c, unchanged**.
- P1/P1w gain the freed selection references as extra evaluation references for the ZAL arms, each
  reported separately.

**P4c-dev, repaired a priori.** It failed numerically in phase 2 through separation. The repair,
fixed now: genes with < 10 training counts in a segment are excluded from both models, and spline
coefficients get a ridge of 1. All arms, old and new, are rerun with it.

**New, better-powered test, reported beside P4c, never in its place: WSR (within-segment residual
transfer).**
- **Rationale.** Within-segment position explains a small share of each held-out gene's variance.
  An MSE gain is diluted by count noise, penalized by spline variance, and blurred by
  between-specimen level shifts, so P4c can miss real, consistent direction. A correlation summed
  over thousands of genes, with a permutation null, has power for small consistent effects and
  ignores levels.
- **Procedure.** Per species, gene fold and direction (train specimen → test specimen):
  1. In **each** specimen separately, residualize held-out genes (log-normalized) on that specimen's
     own segment + coverage quadratic oracle. This uses labels and coverage only.
  2. In the training specimen, fit a per-gene, per-segment cubic spline (df 4) of the
     training-range-scaled coordinate to the residuals.
  3. Predict test residuals at test positions, clipped to the training range.
  4. Statistic: the mean over genes of atanh(Pearson(prediction, observed residual)) within segment.
  5. Null: 500 within-segment permutations of the test coordinate (seed 61).
- **A species × segment cell passes if** p ≤ 0.01 in both directions **and** WSR exceeds the
  placebo coordinate (§2.7) in both directions.

**Adoption rule.** Phase 2's (A) and (B) are **unchanged**, including the P4c criteria. WSR, the
placebo comparison and notebook 58(i) are reported beside them. They can support a stated "the
coordinate carries replicated within-segment information" finding, but not adoption.

> **Decision recorded (coordinator, 2026-10-05, before any phase-3 fit): option (b). P4c moves from
> an adoption gate to a claim gate.**
>
> - **Why.** No arm passes P4c: it is ≤ 0 in 44 of 48 cells for SCF13, DPT13, LRS and LCP. P4c
>   therefore measures a limit of the data, not the quality of a method. This was known from phase 2,
>   before any phase-3 fit, and the change is not tuned on phase-3 results.
> - **Adoption** uses (A) and (B) without the P4c items, plus the component ablations against LRS
>   (§1). If no component earns its margin, the adopted method is "LRS with uncertainty".
> - **Claim gates**, reported for every arm whatever is adopted:
>   - "The coordinate improves held-out-gene prediction beyond discrete segments" only if P4c > 0 as
>     originally defined (≥ 5/6 mouse and ≥ 4/6 human cells, mean above SCF13, DPT13 and LRS);
>   - "The coordinate carries replicated within-segment information beyond discrete zones" only if
>     WSR passes (§4) **and** the placebo coordinate fails it.
> - **Deviation table row:** "P4c role changed from adoption gate to claim gate before the phase-3
>   fits, because no existing arm, generic or landmark-based, passes it (data limit, not method
>   quality)."

---

## 5. Downstream (workstream B notebooks, run through their coordinate parameters)

For every surviving arm, and for the placebo coordinate:
- notebook **58 (i)**: within-segment slope agreement between specimens against step nulls;
- notebook **60**: the resolution budget.

`--construction-genes` = that arm's landmark union, for fold protection.

**Pre-registered reading.** A cell counts as continuous information only if it:
- beats every step-null replicate (notebook 58's rule); **and**
- beats the placebo coordinate's slope agreement.

Notebook 57 (beyond-step pathway calls) is optional, on the finally chosen arm only.

---

## 6. Work order, cost and stop points

| Step | Notebook | Content | Cost | Stop if |
|---|---|---|---|---|
| 1 | — | Module: zone order and selection, gauge-anchored EM on a K = 200 grid, design-effect T, product-of-halves posterior, LRS-U, WSR test, repaired P4c-dev; synthetic tests | 1 day | Synthetic tests fail |
| 2 | 52 | Ablation ladder LRS → M0 → M1 → balanced → calibrated, on **external** landmarks, so it is comparable to phase 2 | about 60 fits, ≈ 1 h on 24 workers | Every component fails → LRS-U |
| 3 | 53 | ZAL (data-derived landmarks, cross-fitted) vs external, plus placebo and label-prior ablation | about 60 fits | ZAL abandon rule (§2) |
| 4 | 54 | Full benchmark with WSR for all arms, calibration, species curve test, human S1 rule | ≈ 1 h | — |
| 5 | B's 58(i) and 60 | Surviving arms and placebo | B's runtime | — |

Every protocol detail above goes into the ledger draft before step 2. Nothing is tuned on evaluation
outputs. Null results are reported as they come.

**Files to create in phase 3:**
- `pseudospace/landmark_position.py`: extended. Proposed as additions; existing functions are kept for
  notebook 51's reproducibility.
- `tests/test_landmark_position.py`: extended.
- Notebooks 52–54.
- `workstream-A scratch: ledger_entry_phase3.md`.
