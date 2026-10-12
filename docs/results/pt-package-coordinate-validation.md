# The package coordinate: validation (notebook 72)

Part of the switch to the package pseudospace; see `pt-package-pseudospace.md`. The protocol was fixed before
any fit it governs.

## Notebook 72 protocol: validation of the package coordinate

Draft by ws72, 2026-10-12. It includes the coordinator's decisions on O1–O5 (section 9). No fit on private data has
been run. Phase 2 waits for notebook 70's export and the coordinator's sign-off on its labels. This protocol must be
frozen before notebook 72 runs.

**Question.** Does the package pseudospace (stage 04: an elastic principal curve, nephron-wide, with PT as its first
stretch) measure position along the PT, and only position? The test reruns, on the package's units, the checks that
notebooks 35, 36, 51–53 and 60 ran on SCF13:
- external agreement;
- within-segment order;
- robustness to gene folds, count splitting, equal depth and an injected group effect;
- positional noise.

Only the inputs change. Each refit reruns the package itself, unmodified.

**Why.** The user's decision adopted the package, and notebook 72 neither selects nor un-adopts it. It reports, against
criteria fixed before the earlier fits (notebooks 35, 36 and 51–53), what the package coordinate supports and what the
paper may claim from it. Where 72 needs a rule that did not exist before, it is marked **new**.

**Labels.**
- Notebook: `analysis/notebooks/72_pt_package_coordinate_validation.ipynb`.
- Logic version: `72.package_coordinate_validation.1`.
- Results: `results/pt_package_coordinate_validation/`, with `tables/`, `figures/`, `refits/{hvm,aki}/<refit>/`,
  `stage_cache/` and `run_record.json`.

### 1 · Inputs

- **The export.** `results/package_pseudospace/{hvm,aki}/pt_units.h5ad`, read through `package_inputs.read_pt_units`.
  The analysed position is `package_inputs.support_positions` (notebook 71's rule: min–max over each specimen's
  [1st, 99th] common support).
- **The run folders.** `results/package_pseudospace/runs/{hvm,aki}/` (notebook 70), read-only:
  - `02_qc/units_qc.h5ad`, the base of every refit;
  - `03_integration/pass2.h5ad`, for the embedding, the integration features and the non-PT class means;
  - `04_pseudospace/pseudospace.parquet`.
- **Configurations.** The tracked `analysis/package_configs/{aki,human_vs_mouse}/`, unchanged.
  - `PSEUDOSPACE_VISIUM_ROOT` (or `--visium-root`) is required, as in notebook 70, because the configs expand it.
- **External references**, as in notebook 35: male and female mouse snRNA, human snRNA halves A and B (Census), mouse and
  rat microdissection.
- **Notebook 71, read-only:**
  - `gene_statistics_hvm.csv` (T_spatial for A3);
  - its pass-rule table (the downstream gate);
  - its universe and gene sets, recomputed with `hvm_inputs` and checked equal to 71's.
- **Comparators (section 6 only).**
  - SCF13 and DPT13: notebook 35 and 36 stage caches.
  - LCP: `results/pt_reconstruction_v3/lcp_coordinate.csv` and the nb51 cache.
  - LRS-Z: `results/pt_reconstruction_v3/coordinates/lrs_z.csv` (full and three folds) and nb53's `metrics_long.csv`.
  - AKI: `results/mouse_only_v5/all_mouse_tubules_scanpy_dpt.h5ad` (`pt_subset_scanpy_dpt`).

**Gene model.** D2 applies to every gene model in 72: the canonical-marker curves, A3 and the count-split test (G6). Each
is an NB fit on raw counts with the mt-free offset, using notebook 71's `standard_gene_model` and positional F-test.

D2 covers gene-model fits only (decision O3). P1, P1w and P4/P4c are rank correlations and held-out prediction, not
gene tests, so they keep their inherited definitions on log1p(count / full library × 10⁴). That keeps their inherited
thresholds meaningful.

### 2 · Units, labels and gene sets (shared by every arm)

- **Units.** Each dataset's PT export, all units.
- **Segments.** The package's `segment_class`. The same base labels are used for every refit, as notebook 35 used its
  reviewed labels for every refit.
- **Species groups for P1w.** Human uses early (S1 + S2) and S3, as in notebook 35.
- **Anchors and gene sets, rebuilt on the package units with notebook 35's rules** (`anchor_coordinate`):
  - **Anchors** A_LOOSE (|Δ| ≥ 0.5) and A_STRICT (|Δ| ≥ 1):
    - ≥ 10% detection per species in package PT units;
    - PT-specific by the class means of the package's non-PT families;
    - built from the selection references.
  - **Positive controls** D_OPPOSITE and D_FLAT.
  - **H:** A_STRICT minus the base run's integration features (pass 1 ∪ pass 2).
  - The overlap with notebook 35's saved sets is reported.
  - AKI uses the human vs mouse A_LOOSE for exclusion, because it has no human reference.
- **Evaluation genes E.** Notebook 35's rule: not A_LOOSE; detection ≥ 10% in the evaluated species or group; external
  level ≥ 1. Also excluded are the 15 genes that set the package's labels or orientation: the PT-S1/S2/S3 panels ∪ the
  axis anchors. These cover notebook 35's six markers, plus Slc22a6, Slc13a3, Cyp2e1 and the six late CD anchors.
- **References.**
  - P1 uses evaluation references only, following notebooks 51–53: male mouse snRNA, and human snRNA half B.
  - P1w uses the male-snRNA adjacent contrasts (S2−S1, S3−S1, S3−S2) in mouse, and half B (early, S3) in human.
  - Microdissection is reported for continuity. It is a selection reference for LCP and LRS-Z, so it enters no rule.
- **Folds.** `sha256_folds(genes, 3)`, the same hash as notebooks 35, 36, 51 and 53, so fold k holds the same genes in
  every arm.

### 3 · Refit plan (the package, unmodified, on perturbed inputs)

**Route.** For each refit, write a new run folder `refits/<ds>/<name>/` containing:
1. `02_qc/units_qc.h5ad`: the base QC-retained units (same obs, var and obs order) with the perturbation applied, and
   `layers['lognorm']` and X recomputed by the package's own `qc.pipeline.normalize_log1p` (target 10⁴);
2. a copy of the base `02_qc/manifest.json`, with the perturbation recorded.

Then run `pseudospace-recon run 03 04 --config <tracked config> --data-root <data> --out <refit folder>` as a subprocess.
This is the same CLI route notebook 70 used. DPT is on, as in the config, so every refit also gives a package-DPT refit.

**Stages rerun.** Stages 01–02 are not rerun: the unit set and QC are held fixed. Notebook 35 held its structure set
fixed. Stage 02's split test also re-reads the bin matrices, which a unit-level perturbation does not touch, so a
perturbation cannot pass through stage 02 consistently.

Every perturbation below changes the counts or the integration genes, so it reruns 03 and 04. Only LOSO leaves the
embedding unchanged, and it reruns only the stage-04 curve functions.

| Refit | Perturbation (seeds from notebook 35) | Stages | hvm | aki | Feeds |
|---|---|---|---|---|---|
| identity | none (rewritten copy) | 03–04 | 1 | 1 | reproduction guard |
| floor | lognorm as log1p(x · 10⁴/total), multiply instead of divide; only pass 1 reads it, because stage 03 renormalizes from counts after the gene filter | 03–04 | 1 | 1 | numerical floor |
| fold0–2 | var `exclude_from_harmony_hvg` set on fold-k genes (notebook 35's withholding) | 03–04 | 3 | 3 | fold-protected P1/P1w and positive controls; gene-fold agreement; gap; noise; P4/P4c |
| half A, B | Poisson count split, ε 0.5, seed 35, all units | 03–04 | 2 | 2 | A-vs-B agreement; gap; count-split noise; G6 |
| equal depth | binomial thinning of all units to the pooled PT 20th percentile of the library, seed 352 | 03–04 | 1 | 1 | P5 |
| injection | case-group PT units (human or IR), described below | 03–04 | 1 | 1 | P2 |
| LOSO | none; curve refit on the base embedding without one specimen | 04 functions | 4 | 4 | LOSO |

**The injection.**
- Genes: 10% of the base run's pass-2 integration features (sha256 mod 10 = 0), excluding A_LOOSE, the panels and the
  anchors.
- Effect: × exp(±log 2 · (2s − 1)), where s is the base position on 71's scale clipped to [0, 1]. Signs alternate in
  hash order.
- Applied with `scale_counts`, seed 351.

**Matching back.** Unit ids (`<sample>_<unit>`) are unchanged through stages 03–04.
- A refit position is the refit's `pseudospace.parquet` value for that id. It is NaN when the refit dropped the unit
  from pass 2.
- Scale-dependent metrics (noise, gap) use each coordinate on 71's rule, recomputed on that coordinate.
- Package DPT gets notebook 70's flip rule in every refit.

**Recorded per refit.**
- The share of base PT units missing from pass 2. **New:** more than 5% is flagged as unit-set instability.
- Label agreement with the base (κ over PT-S1/S2/S3, and PT family retention).
- The root, the fingerprints and the integration features.

**Guards (not cached).**
1. The identity refit reproduces notebook 70: Spearman ≥ 0.999 over PT units, identical fingerprints, and max |Δ|
   reported. Otherwise stop.
2. `normalize_log1p` on the base counts reproduces the base lognorm exactly.
3. No fold-k gene is among a fold-k refit's pass-1 or pass-2 integration features.
4. A refit whose PT S1 median is not below its S3 median in some specimen counts as **failed**. It is never repaired, and
   its metrics are NaN, as in notebook 36.
5. Package version 0.1.0 at 59da04b, and R harmony 2.0.5.

**Package limits (accepted, decision O4), stated with every refit result.**
1. **QC is fixed across refits.** Stages 01–02 are not rerun, so the unit set and its QC are the base run's. The
   refits test the coordinate given the units, not the QC.
2. **Fold genes still shape the labels.** The package builds its clustering feature space from all genes, so in a
   fold-k refit the fold-k genes still shape the labels, the pass-1 drops and the root cluster. They do not shape the
   embedding. Withholding them from the labels too would mean modifying the package.

**LOSO** (notebook 36's curve-only design, using public package functions only):
1. Fit `curve_pseudospace` on the base embedding without specimen n, rooted at the base root, or at `choose_root` on
   the remaining units if the root lies in n.
2. Place n's units with `project_onto_curve` and `curve_pseudotime`, signed by the training orientation.
3. Report the Spearman correlation within n's segments against the base.

Harmony saw all specimens, as in notebook 36, so this tests the curve, not the integration.

**Compute** (AKI trial on this machine: stage 03 135 s, stage 04 13 s, DPT about 1 min):
- **Refits:** about 4 min each including input write and extraction. 18 runs take about 75 min sequentially, or about
  40 min with 2 concurrent runs (`PSEUDOSPACE_REFIT_JOBS`, default 2; BLAS threads split; never in a cache key).
- **LOSO:** 8 curve fits, about 3 min.
- **G6:** 3 NB condition-split fits on half B, about 30–60 min with `PSEUDOSPACE_N_JOBS`.
- **Markers, A3, P4/P4c:** about 15 min.
- **Total:** about 2–2.5 h wall.

**Disk.** A refit holds about 6–7 GB while it runs: the perturbed `units_qc`, `pass1_labelled` and `pass2`.
- After extraction, those three files are deleted. The CSVs, manifests and parquet stay (< 50 MB per refit).
- A refit starts only if ≥ 40 GB is free. 124 GB was free on 2026-10-12.
- Each extraction is cached with `cached_payload`. The key holds:
  - the refit spec and base digests (units_qc counts, config sha256);
  - the package commit and the logic version.

### 4 · Validation criteria and what each result lets the paper claim (package arm; package DPT alongside)

**The package is adopted by the user's decision.** Every row below is a validation criterion. It is reported as met or
not met for the package coordinate, with what the paper may claim either way. No criterion can un-adopt the package.
Notebooks 51–53 used groups (A) and (B) as adoption rules for candidate arms; here they are validation criteria only.

- The criteria are defined on human vs mouse and are judged there only. AKI is reported against them descriptively
  (section 5).
- "Relative" criteria compare with SCF13 or DPT13 on the matched set (section 6). If the comparator section is dropped,
  they compare with the published values on the comparators' own units, flagged "different units".

| # | Validation criterion (source) | Rule | Not met: reported as | Met: the paper may claim |
|---|---|---|---|---|
| C1 | Canonical order (33/37 transitions; **new** rule) | logistic S1→S2 < S2→S3 inside the support in every specimen; Slc7a13 and Slc22a7 rise (Spearman > 0) in both mice | "canonical order not recovered in specimen X" | S1→S2→S3 recovered without labels in the curve (labels and orientation share panel genes: say so) |
| C2 | P1, fold-protected (35; 51–53 (A)) | ≥ SCF13 − 0.03 in each species | less external zonation than the previous coordinate | external S3-vs-early zonation recovered as well as SCF13; labels shown as the upper reference |
| C3 | P1w (51–53 (B)) | ≥ DPT13 + 0.05 mouse, ≥ DPT13 human | within-segment order weakly supported (as SCF13: 0.17/0.10) | within-segment order beyond the old coordinates (Fig 2d) |
| C4 | Gene-fold agreement (36; (B)) | ≥ 0.70 mouse, ≥ 0.59 human | within-segment order depends on which genes build it | stable within-segment order |
| C5 | LOSO, curve-only (36; (B)) | ≥ 0.85 | curve depends on one specimen | — (weak test; reported) |
| C6 | P5, equal depth (35; (A)) | ≥ 0.80 overall, ≥ 0.70 human S1; ρ(position, log library) before and after | position partly tracks read depth; depth caveat on every positional result | position is not driven by depth |
| C7 | Human S1 resolved (phase-3 plan) | gene-fold ≥ 0.50 and P5 ≥ 0.70 in human S1 | "human S1 unresolved" | — |
| C8 | Registration gap (36; (A)) | max over fold and half refits of max over segments of \|human − mouse\| ≤ 0.15 | species placement unstable; read cross-species curves in segment units | species segments co-registered on the shared axis |
| C9 | P2 absorption (35) | no material absorption: median ratio ≥ 0.85 **and** human within-segment ρ(refit, base) ≥ 0.90 | the coordinate absorbs group × position effects; 71's shape calls are lower bounds | a species × position effect survives refitting |
| C10 | Downstream gate ((A); read from 71) | NCDR < 0.25 and each relabeling < 5% of tested sets | as 71's M1 | — |
| C11 | Positive controls (35) | sign agreement of r_human − r_mouse with the external interaction, beside labels (35: 91% vs 95%) | descriptive | — |
| C12 | A3 attenuation (35) | flag if ≥ 70% of H genes have human NB amplitude < mouse, or median 71 T_spatial percentile ≥ 0.75 | species amplitude differences partly positional noise | — |
| C13 | G6 count split (35; **new** rule) | arms on half B: (a) refit on A, (b) base coordinate, (c) refit on B; positional F-test + cameraPR; inflation of median F, b/a and c/a; not met if a's calls < 0.5 × b's | calls depend on reusing reads; count-split arm becomes primary | calls are not created by building and testing on the same reads (Fig 2e) |
| C14 | Noise budget (60; **new** reading rule) | per specimen × segment: fold and count-split noise SD in segment units and ÷ within-segment spread. A cell is noise-dominated if fold noise ÷ spread ≥ 1; not met if ≥ 6 of 12 cells | no within-segment resolution claim; 74 restricted to between-segment structure | noise below the within-segment spread (information side in 74) |
| C15 | P4/P4c claim gate (51; 52/53 option b) | P4c > 0 in ≥ 5/6 mouse and ≥ 4/6 human cells (fold × direction) | expected; no held-out-gene prediction claim beyond segment + coverage | continuous position predicts held-out genes beyond segments |

**G6 and noise details.**
- G6 uses 71's universe and gene sets throughout.
- Each arm keeps its own common support, as in notebook 35.
- The offset is the mt-free half-B library.
- Count-split noise comes from half-read refits, so it overstates full-data noise. It is reported, but C14's rule uses
  the fold noise.

**Reported beside the rules.**
- The numerical floor and A-vs-B agreement within segment.
- The order across refits: human S1 < S2 < S3 in each refit.
- The transitions per specimen in segment units.
- The positive-control and A3 references (labels, external).

**Summary rows** (met or not met, never pass or fail). Group (A) without the gate is C2, C6, C8 and C9. Group (B) is
C3, C4 and C5.

**Descriptive panels** (figure plan, Figure 2).
- **(a) Canonical markers** Slc5a2, Slc5a12, Slc34a1, Slc22a6, Slc7a13 and Slc22a7, shown as:
  - NB group curves from `standard_gene_model` on 71's position and grid;
  - specimen curves;
  - specimen means in 20 equal-count bins.

  Only Slc34a1 is outside the panels and anchors; the legend says so.
- **(b) Composition and transitions.** Stacked S1/S2/S3 shares in 20 bins per specimen, with the logistic transitions.
- **(c) P1 gene scatter** against the external contrast.
- **(d) Within-segment order panel.** P1w, gene-fold agreement and LOSO. Read depth is the nuisance reference and labels
  the P1 reference.
- **(e) Robustness.** G6 and P2.

### 5 · AKI (the coordinate behind 71's AKI arm and Figure 6)

**Refits.** The full set of 9 stage-03–04 refits (decision O5), plus 4 LOSO fits, with the metrics of section 4 and
mouse references. Control and IR mice are reported separately.

**P2.** The injection goes into IR PT units. Absorption uses r_IR − r_control; the within-segment ρ is taken in the IR
mice.

**Reporting: descriptive only.** The criteria are defined on human vs mouse, so AKI gets no met / not-met verdict. Each
AKI value is shown beside the matching human vs mouse threshold:
- gene-fold 0.70, P5 0.80, LOSO 0.85, and C9's P2 thresholds;
- for P1 and P1w, the published SCF13 and DPT13 values for the same two control mice on their own units
  (P1 0.745 − 0.03; P1w 0.23 + 0.05).

**Also descriptive:**
- the IR − control registration gap and transition shifts, because injury can move labels.

**Not run:** positive controls, A3 and G6.

**Use in Figure 6.** Figure 6's text quotes these values as descriptions, especially P2, gene-fold and P5 in the IR mice.
It states no AKI verdict.

### 6 · Comparators (separable; runs by default under D4)

This section runs by default (decision O1). It is self-contained, so it can be dropped if the user decides against D4;
sections 4 and 5 do not depend on it, except for the relative criteria (see section 4).

No comparator is rebuilt on the package's units (decision O2). Comparators enter only through their cached positions on
matched units.

**Arms.**
- SCF13 and DPT13: notebook 35 refits (full, fold0–2, equal depth, halves, injection, floor) and notebook 36's SCF13
  LOSO.
- LCP: notebook 51 (full, folds, equal depth, halves, own-landmark injection, LOSO).
- LRS-Z: full and folds only, since notebook 53 did not save its other refits.
- AKI: the v5 PT DPT, full data, descriptive.

**Common unit set M.** Package PT units whose `old_structure_id` (notebook 70's IoU ≥ 0.5 mutual-best match) is an
SCF13 PT structure; for AKI, a v5 structure with a finite `pt_subset_scanpy_dpt`.
- A comparator's position on unit u is its value for old_structure_id(u).
- Every arm on M uses the same package counts, package labels, gene set E, references and folds.
- The package and package DPT are scored on M too. Relative criteria (C2, C3) use M.

**Unmatched units.**
- Package units outside M stay in sections 4–5 and leave only this section. Their share by specimen × segment, and their
  depth against M, are reported.
- Old structures without a package unit are counted, not analysed.

**Recomputed and quoted metrics.**
- **Recomputed on M from cached positions:** P1/P1w (fold-protected), gene-fold, P5, A-vs-B, gap, LOSO, injection
  within-segment ρ, noise budget and P4/P4c.
- **Quoted from the comparator's own units, labelled so** (they need the old injected or perturbed counts):
  - the P2 ratio;
  - the downstream gate (Gaussian, notebooks 37 and 53);
  - LRS-Z's P5, halves and gap.

**Sensitivity.** The within-segment metrics again with the old reviewed labels.

**The matching penalty.** Mouse comparators were built from the old polygons' counts, so on M they are scored against a
different polygon's expression. Each comparator's own-unit value (nb53 `metrics_long.csv`) is shown beside its value on
M; the difference is that penalty.

**Output.** A head-to-head table with every arm × criterion, losers included. No arm is selected here: the package is
the method by the user's decision.
- If LRS-Z or LCP beats the package on C3/C4/C7, the supplement says so. The paper then frames the package as the
  pipeline and the landmark ratio as an alternative ordering step.

### 7 · Outputs

**Tables** in `tables/`:
- refit metadata (`refit_record_{ds}.csv`);
- gene sets and counts;
- `p1`, `p1w` and `positive_controls`;
- `gene_fold`, `loso`, `p5`, `count_split_agreement` and `registration`;
- `p2_absorption`, `g6_count_split`, `a3_attenuation`, `noise_budget` and `p4_p4c`;
- `transitions` and `markers_curves`;
- `criteria_{ds}.csv`, one row per criterion with value, rule and source. Human vs mouse rows also hold met / not met;
  AKI rows hold the human vs mouse threshold beside the value;
- `comparators_matched.csv` (section 6) and `unit_set_M.csv`.

**Other outputs.** Figures 2a–e, the S4 head-to-head panel, and `run_record.json`.

**Code.** New module `pseudospace/package_refits.py`:
- write a perturbed stage-02 folder;
- run stages 03–04 through the CLI;
- extract, match back, guard orientation and clean up.

Its tests (`tests/test_package_refits.py`) are synthetic:
- the perturbed writer reproduces `normalize_log1p`;
- ids survive;
- fold flags reach the HVG step;
- the orientation-failure path is exercised.

No existing module changes.

### 8 · Not carried over, with the reason

- **n_map.** The package has no posterior mappings.
- **scFates seed stability.** The package's curve seed changes only round-off; the floor refit covers it.
- **Retired arms:** SCF13-ED, CAC and SCF-PT.
- **A3b noise matching.** It ran on the retired 66-pathway list.
- **WSR and the placebo claim gate.** They need a placebo coordinate and belong to notebook 74's information side.

### 9 · Decisions (coordinator, 2026-10-12)

- **Framing.** The package is adopted. Section 4's rows are validation criteria, reported as met or not met with the
  claim each result allows; none can un-adopt the package.
- **O1.** Section 6 runs by default under D4, and stays separable.
- **O2. Declined.** No LRS-Z rebuild on package units; comparators enter only through their cached positions on matched
  units.
- **O3.** D2 applies only to gene-model fits. P1, P1w, P4 and P4c keep their inherited log-normalized definitions, so
  the inherited thresholds stay meaningful.
- **O4. Accepted.** Fixed QC across refits, and fold genes still shaping the labels, are recorded as package limits
  (section 3).
- **O5.** AKI gets the full refit set (9 runs). The criteria apply to human vs mouse only; AKI is reported descriptively
  against the same thresholds (section 5).

## Notebook 72 results

**Run.**
- Notebook: `analysis/notebooks/72_pt_package_coordinate_validation.ipynb`, logic version
  `72.package_coordinate_validation.1` (refits `72.refits.1`).
- Output: `results/pt_package_coordinate_validation/`.
- **Refits:** 19 package runs of stages 03–04 on perturbed stage-02 inputs: 9 per dataset, plus one 12-thread identity
  run for human vs mouse. They ran at most 2 at a time, about 3–6 min each. Their large h5ad files were deleted after
  extraction.
- **Other fits:** 8 curve-only LOSO fits, and the three NB count-split arms on notebook 71's universe and gene sets.
- New module `pseudospace/package_refits.py`, with `tests/test_package_refits.py` (6 tests). The full suite passes
  (513).

**Guards.**
- Notebook 70's confirmed pass-1 fingerprints are in place: human vs mouse `efafc9ac2e99`, AKI `bb4a80c1b484`.
- Both identity refits reproduce notebook 70 bit for bit (maximum |Δ| 0).
- No withheld gene leaked into a fold refit's integration features.
- **One refit failed the orientation guard:** human vs mouse fold 2. The curve runs through the mouse PT from S3 to S1
  and through the human PT from S1 to S3. Its positions are NaN, by notebook 36's rule.
- Unit-set instability (> 5% of base PT units missing from pass 2) was flagged in four refits: human vs mouse half B,
  and AKI fold 1, half A and half B.

**Package sensitivity found while building the refits (beyond the protocol; reported, not tuned).**
- The package's stage 03 is chaotic. Two changes alone altered the pass-1 clustering and moved PT units:
  - an equal-valued lognorm layer with unsorted column indices;
  - 12 instead of the default BLAS/OpenMP threads.

  In the shipped notebook, the 12-thread identity has within-segment ρ 0.98 against the base and a different
  fingerprint.
- The numerical floor (1e-16) gives within-segment ρ 0.96.
- So the refits keep the package's own file layout and its default threads, as notebook 70 ran. Results of the
  package depend on the machine's thread count.

**How the protocol was read** (listed in `run_record.json`).
1. **Library.** The library for P1/P1w/P4 and for equal depth is the total over every export gene. In human vs mouse
   these are all shared orthologs, which is notebook 35's measured-in-both total.
2. **Excluded genes.** Evaluation genes and positive controls exclude the 15 label/orientation genes (the PT panels
   plus the axis anchors).
3. **Label references.** P1's label reference is notebook 51's (labels coded 0/1/2 as a coordinate). The positive
   controls' reference is notebook 35's S3 − early contrast.
4. **Failed refits** are NaN. Medians and maxima ignore them (notebooks 36 and 51); their P4c cells count as not
   positive.
5. **Gap and noise** rescale each coordinate by notebook 71's rule within the context's units. The gap is 0.369 on all
   PT units and 0.217 on M. C8 uses all PT units.
6. **C12** counts as not met when the flag is raised. The label-based and external fractions are shown beside it.
7. **LOSO** is computed for the package only (no DPT analogue). LRS-Z's gap uses its folds only. LCP's LOSO, and every
   comparator's P2 ratio and downstream gate, are quoted from their own units.
8. **Notebook 71's tables keep notebook 69's dataset names**: hvm = `human_vs_mouse`, aki = `aki_vs_control`. That
   covers `gene_statistics_human_vs_mouse.csv`, `universe_{ds}.csv`, `gene_sets_{ds}.json`, `method_table.csv` and the
   cached stage `gene_model_human_vs_mouse_species`.
   - C10 reads `method_table.csv` (normalized · cameraPR, human vs mouse).
   - C13 runs on 71's universe (11,055 genes) and gene sets (1,513), which equal the loader's. Its tested sets equal 71's.
   - C12's percentile is 71's condition-split T_spatial.
9. **A3 label amplitudes** are specimen pseudobulk log rates, the count analogue under D2.
10. **P4c sensitivity.** A few PT units lie far along the nephron-wide curve, and the within-segment spline of P4c
    extrapolates there: one mouse cell reaches −1.43. A common-support sensitivity is reported beside C15 and not used
    for it.

**Validation criteria, human vs mouse** (met or not met; none un-adopts the package).

| # | Value | Verdict |
|---|---|---|
| C1 canonical order | HUK1_MED1 transitions cross (S1→S2 0.666 > S2→S3 0.590; human PT-S2 has 259 units); Slc7a13 and Slc22a7 rise in both mice (ρ 0.69–0.85) | not met |
| C2 P1 on M | mouse 0.746 vs SCF13 0.746; human 0.295 vs 0.246 (labels 0.749 / 0.314) | met |
| C3 P1w on M | mouse 0.367 vs DPT13 0.288; human 0.156 vs 0.130 | met |
| C4 gene-fold | mouse 0.751; human 0.488 | not met |
| C5 LOSO | 0.959 | met |
| C6 P5 | 0.731; human S1 0.552 | not met |
| C7 human S1 | gene-fold 0.516; P5 0.552 | not met (unresolved) |
| C8 gap | 0.369, from half A (M 0.217; full fit 0.076) | not met |
| C9 P2 | ratio 0.985; human ρ 0.922 | met |
| C10 gate (71) | 83 calls; relabeled 1, 0; NCDR 0.006 | met |
| C11 positive controls | 0.899 vs labels 0.959 | descriptive |
| C12 A3 | flag raised: 0.829 of H genes are flatter in human (labels 0.829, external 0.771); 71 T_spatial percentile 0.943 | not met |
| C13 count split | calls: a 66, b 81, c 66; median F b/a 1.22, c/a 1.06; z of a vs b ρ 0.52; 71's 83 calls kept by a / b / c: 59 / 73 / 61 | met |
| C14 noise | 10 of 12 cells noise-dominated (median ratio 1.73) | not met |
| C15 P4c | mouse 1/6, human 0/6 (common support 2/6, 0/6) | not met |

The package DPT is reported alongside. It is weaker on most criteria (gene-fold 0.54 / 0.29; P5 0.70), and its gap
reaches 0.96 in the half-A refit.

**AKI (descriptive, against the same thresholds).**
- Transitions are ordered in all four specimens.
- P1: control 0.745, IR 0.363. P1w: control 0.300, IR 0.181.
- Gene-fold: control 0.812, IR 0.705. LOSO 0.976. P5 0.890.
- P2: ratio 0.984, IR ρ 0.961. IR − control gap 0.086.
- Noise-dominated cells: 11 of 12. P4c > 0 cells: control 3/6, IR 3/6.

**Comparators on M.**
- **M (human vs mouse):** 11,488 of 11,947 package PT units (94.5–96.7% per specimen). 1,378 SCF13 PT structures have
  no package PT unit.
- **M (AKI):** 14,383 of 15,013 package PT units.

| Metric | package | package DPT | SCF13 | DPT13 | LCP | LRS-Z |
|---|---|---|---|---|---|---|
| P1 mouse / human | 0.746 / 0.295 | 0.751 / 0.215 | 0.746 / 0.246 | 0.741 / 0.290 | 0.745 / 0.320 | 0.746 / 0.316 |
| P1w mouse / human | 0.367 / 0.156 | 0.313 / 0.112 | 0.240 / 0.107 | 0.288 / 0.130 | 0.428 / 0.202 | 0.439 / 0.195 |
| gene-fold mouse / human / human S1 | 0.760 / 0.487 / 0.516 | 0.540 / 0.287 / 0.300 | 0.537 / 0.229 / 0.195 | 0.634 / 0.590 / 0.601 | 0.962 / 0.820 / 0.647 | 0.933 / 0.965 / 0.979 |
| P5 / human S1 | 0.739 / 0.556 | 0.696 / 0.506 | 0.658 / 0.550 | 0.842 / 0.749 | 0.982 / 0.652 | own units 0.976 / 0.980 |
| count split A vs B | 0.255 | 0.205 | 0.071 | 0.394 | 0.887 | – |
| gap (max over refits) | 0.217 | 0.984 | 0.352 | 0.334 | 0.292 | 0.072 (folds) |
| LOSO | 0.963 | – | 0.935 | – | own units 0.99 | – |
| noise-dominated cells (of 12) | 8 | 4 | 7 | 2 | 0 | 0 |
| P4c > 0 cells mouse / human | 2 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 2 / 0 | 0 / 0 |

- The quoted own-unit values (P2 ratio, downstream gate, notebook 53's summary) are in
  `comparators_quoted_own_units.csv`.
- The old-label sensitivity is in `comparators_matched_old_labels.csv`.
- AKI on M (full data, not fold-protected): P1 control 0.747 / 0.750 / 0.752 for the package, package DPT and v5 DPT;
  P1w control 0.326 / 0.296 / 0.266.

**What the paper may claim.**
- The package coordinate recovers external S3-versus-early zonation as well as the labels and SCF13 do.
- Its within-segment slopes agree with external adjacent contrasts better than SCF13 and DPT13.
- Notebook 71's positional calls pass its gate, and the count split does not create them.
- An injected group-by-position effect is not absorbed, in either dataset.

**What it may not claim.**
- Stable within-segment order. Gene folds disagree (human 0.49), and one fold refit reverses the mouse PT.
- Depth-invariant positions (P5 0.73; human S1 0.55, unresolved).
- A stable species registration (gap 0.37).
- Positional noise below the within-segment spread (10 of 12 cells).
- Held-out-gene gains beyond segments plus coverage.

The landmark coordinates (LCP, LRS-Z) remain more stable within segments on the same units. Continuous claims should
rest on between-segment structure and specimen replication (notebook 74).
