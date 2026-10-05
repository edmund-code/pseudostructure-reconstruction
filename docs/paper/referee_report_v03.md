# Referee report on draft v0.3: "Where human and mouse proximal tubule differ along the tubule axis"

Reviewed: `docs/paper/draft.md` (v0.3), `docs/paper/outline.md`, `docs/paper/referee_report_v01.md`,
`docs/paper/referee_response_v01.md`, the saved tables under `results/` (notebooks 35–46 and their
upstreams) and the figures in `results/paper_figures/` (Figures 1–5, S9, S10 inspected).

Reviewer checks are marked **[reviewer check]**. One of them re-ran notebook 43's own index code
(`pseudospace.zonation_reliability.conservation_index`) in memory with different winsorisation
settings. It reproduces the published index exactly at the published setting. The script and its
output are in `SCRATCH/review/nb43_check/` (`rerun_index.py`, `index_winsor_sensitivity.csv`). Nothing
in the repository or `results/` was modified.

---

## Verdict

**Major revision.**

The central biological result is supported and robust. Relative to same-species reproducibility,
human–mouse agreement of PT gene gradients is low (index about 0.1–0.35 across contrasts, labels,
winsorisation settings and data sources, and in public atlases alone). Several statements built on it
are overstated or misstated, however:
- the "conserved core, mostly transporters" descriptor;
- the "not skewed toward mouse" asymmetry claim;
- the AKI specificity verdict;
- the glutathione and Cyp2e1 placements;
- a selective use of failed-validation Visium donors.

The coordinate that the paper presents as its premise does not deliver the main claim. The text is
about twice journal length and still narrates its internal revision history. Ethics and specimen
records are also missing, and those block submission by themselves.

---

## 1. Numbers checked against source tables

**Verified: about 90 quantitative statements match their tables** (rounding aside). By section:

- **Abstract.**
  - 12,866 structures; 100 cross-species draws.
  - Index 0.17 (0.12–0.22) and 0.28 (0.22–0.33); external 0.27–0.28.
  - 388 = 196 + 75 + 41 + 76; spread 0.58; 26 of 36 (`conservation_index.csv`,
    `class_counts_by_setting_and_contrast.csv`, `posthoc_reference_free_amplitude.csv`,
    `decisions.csv`).
- **R1 and Methods: counts and partitions.**
  - 26,839 (= 7,231 + 6,676 + 6,098 + 6,834) of 37,559 structures.
  - Per-specimen S1/S2/S3 counts: 2,852, 2,808, 3,438, 3,768; S1 5,777, S2 3,806, S3 3,283;
    24,335 nephron structures; 12,272 in common support.
  - Transitions 0.26–0.28, 0.35–0.38 and 0.63–0.64 (`segment_transitions_by_specimen.csv`).
  - Partition ARI 0.979 and 0.981; 98.5%; 12,871 and 12,864.
- **R1: coordinate validation.**
  - P1 0.74, 0.73, 0.53, 0.51, 0.25 and 0.27; unprotected 0.75 and 0.31.
  - Rat 0.10; absorption 1.01; 91% of 244; 59 of 66 and 39 of 42; 145 against 130; 57 against 59.
  - P1w 0.17 and 0.10; gene-fold agreement 0.34–0.70; P5 0.36–0.46; library correlations
    −0.31/−0.38 → −0.04/−0.13; registration gap 0.29 with the fold-2 order break; numerical shift 0.44.
  - DPT13 0.80, 0.72, 0.14, 0.29, 0.59, 0.49 and 0.29 (`nb35/*`, `nb36/selection_rule.csv`,
    `arm_metrics.csv`).
- **R1: label checks.**
  - Glomerulus distances 56/93/297 and 0.65; human 149/168/259 uncorrected, 131/148/228 corrected,
    0.23.
  - Attachment 9 polygons (OR 1.0), 46 (3.0), 237 (2.0) and 543 (4.5).
  - Transfer 0.71/0.54, 0.87/0.80, 0.92/0.87; 81%, 99% and 65%; human S1 fraction 0.44–0.49 against
    0.25–0.34.
  - Integration genes 16% against 2%; spillover 0.07/0.10/0.45% and 0.14/0.13/0.29%; zero-bin
    polygons inside the hull 0.2% and 1.2%; coverage 97.4% and 76.5% (`pt_revision_labels/*`).
- **R1: coordinate choice and the attenuation check.**
  - 22 of 36; 32 of 36; 40 of 40; ρ 0.71; 13/3/0/1/5; 37/5 against 16/0
    (`d4_coordinate_runs.csv`, `coordinate_robustness*.csv`).
  - σ 0.089 → 11.9/24.1, 15% and 7%; σ 0.10 → 14% and 8% (`final_checks/attenuation_36_summary.csv`).
- **R2: our cohort.**
  - The full R2 table (`table_S1_strategy_specificity.csv`).
  - Minimum FDP 1.00/0.75/0.95; 1,087 decoy calls; DESeq2 9,010/2/8 (`conventional_deseq2_all_partitions.csv`).
  - Decoy a and ρ 1.40/0.019 and 1.03/0.009; 59.5 against 15; within-species 48/16 and 94/0.
- **R2: many-draw nulls and AKI.**
  - Medians 22/40/34/13 and 11/5/0/0; NCDR 1.34, 0.69 and 0.05; shares 19%, 21% and 80%.
  - Precision 71% against a 3% base rate; calibration 1.0%, 0.6%, 3.8–4.8% and 2.4–7.5%.
  - Coherence 0.26–0.50, 0.03–0.20 and 0.013; age 22, 62 and 1.
  - AKI 26/12/6, 50/10/7, 0.35, 0.17 and 0.42 (`pt_revision_method/tables/*`).
- **R3.**
  - Ceiling table ranges (0.70; 0.70–0.93; 0.63–0.88; 0.73–0.93; 0.09–0.20; 0.13–0.24).
  - Index table, including S3c 0.29 (0.23–0.35).
  - Cross-fit 0.31, 0.38, 0.33 and 0.43; human ceiling 0.64–0.67 → 0.53–0.55; human amplitude about
    0.50.
  - Spread 0.58, 0.47 and 0.84; projection ratios 0.15–0.23, 1.96 (1.41–3.05), 0.14 and 1.66;
    Deming and orthogonal 0.16–0.17; tertile 0.18; probe-balanced 0.16.
  - AKI 0.74 (r 0.83) and 0.41–0.44 (r 0.59); Lake 0.69/0.50/0.35; donor −0.45, 0.50 and +0.33.
  - Time of day 0.99–1.35 and 73%; PDK4 94% against 7%.
  - All class counts (31/58, 13/28, 1/110, 38/71, 31/78; 247/66, 59/141); 33 of 43 transporters
    (77%) and 51%.
  - Funnel 60/58/51/36/26 and 17 programs; categories 26/3/2/5/0; the 8 probe-sensitive pathways;
    26 graded, 4 localised, 6 reversing; peaks 0.16 and 0.63; 27 of 36 and 17 of 26.
- **R4.**
  - 28 of 36 and 32 of 36; split-plot 41/0/0, 15 of 36, 2.6 against −0.5, a 1.10 against 2.07,
    ρ 0.21.
  - Gene table 0.67/0.57 (0.62)/0.64 and 0.73/0.60 (0.66)/0.72; 17 and 6 of 36; p 0.016, 0.033, 0.19
    and 0.21.
  - Slope −0.81 to −0.90, R² 0.60–0.82; 1–3 beyond the global component; 184 and 6.
  - Lake genome-wide 0.43, 6,787, 93.1%, 335, 0.70, 96.7%, 276, 0.67 and 93.7%.
- **R5.**
  - Every gene value in `story_gene_evidence.csv` and `D_rat_segment_proteome.csv` quoted for Gatm,
    Gamt, Acadm, Acaa2, Gclc, Gclm, Gss, Cyp2e1, UGT1A9, Dcxr, Ugt3a1, RBP4, AOX1 and OAT1.
  - Every symmetric class stated for 50 named genes (`gene_classes_symmetric.csv`).
  - 39/30/5/4.
- **Figure legends.**
  - Fig 1a; Fig 2f–h; Fig 3a [7,391; 109/37/31/58]; Fig 3b [0.65, 0.87]; Fig 3d; Fig 4a [22, 40];
    Fig S1e [76/71/79%]; Fig S9 [75–76%]; Fig S10 [71%, κ 0.54; 87%, κ 0.80].

**Mismatches and misstatements (12).** None changes a headline number. Items 1, 4 and 8 change how a
claim must be worded.

| # | Location | Draft says | Table says | File |
|---|---|---|---|---|
| 1 | Introduction ("This study", last paragraph); R3 "conserved core…"; Discussion 2 | conserved core "mostly apical transporters"; "concentrated in strong transport and segment-identity genes" | 33 of 196 conserved genes are Slc/Abc/Aqp (17%). Transporters are *enriched* in the conserved class, but they are not most of it | `pt_revision_classes/reviewed/tables/gene_classes_symmetric.csv` |
| 2 | Abstract | index "at most 0.31 with … transferred" labels | S3c − early, declared co-primary in the R3 table, reaches 0.33 (upper bound 0.43) | `crossfit/tables/crossfit_conservation_index.csv` |
| 3 | R3, cross-fit paragraph | "did not outnumber … in 11 of 12 checks" | the pre-specified decision "species-only counts unskewed" evaluates **False**; the text does not say that a pre-registered criterion failed | `crossfit/tables/decisions.csv` |
| 4 | R2 "Injury severity"; Fig 2h; Discussion | T_spatial "specificity failed … against the 5% condition" (17%) | Under the paper's own rule (Methods 6.14: NCDR < 0.25 **and** relabeled calls < 5% of *tested* pathways), T_spatial passes in AKI: NCDR 0.17, with 10 and 7 relabeled calls of about 1,500. It fails only notebook 45's stricter D6 criterion (relabeled ≤ 5% of *condition* calls). The offset joint screen and the whole-PT joint rank test have the lowest ratios (0.09) | `table_aki_strategy_calls.csv`; notebook 45 protocol |
| 5 | R4 "Six checks keep the effect" | six checks | five are listed and Figure 4a shows eight. Under "detected in both" and "sex-unbiased", 34 and 35 of 36 pathways keep ≥ 0.7 of the effect, not all | `table_S2_pathways.csv` |
| 6 | Methods 6.25; Fig S10 | distances "corrected by this factor (−13.5%)" | the factor is 0.44068/0.5 = 0.881, i.e. −11.9%; 13.5% is the overstatement. The quoted 131/148/228 correctly use −11.9%. Figure S10a plots the *uncorrected* 149/168/259 while R1 quotes the corrected values | `pt_revision_labels/summary.json` |
| 7 | R5 RBP4 | "7 of 8 informative Visium samples" | the pre-specified table has 7 of 7 informative; 7 of 8 comes from the post-hoc centred table, which is not labelled post hoc | `C_visium_lead_genes.csv`, `C_posthoc_centred_lead_genes.csv` |
| 8 | R5 glutathione "Reading"; Discussion 3 and "Implications" | human side "weaker or no gradients, not a different placement" | Lake GCLM S3 − S1 −0.76 (male −1.05) against mouse +2.1 to +2.4, which is opposite. GSS is *conserved* (Lake +0.67) | `story_gene_evidence.csv`; `gene_classes_symmetric.csv` |
| 9 | R1 coordinate choice against Fig S2b | DPT13 gives 37 mouse-vs-mouse calls | Fig S2b says 36. These are different runs (notebook 36 diagnostic against the notebook 37 rerun), which the text does not explain | `d4_coordinate_runs.csv`; `dpt_diagnostic/summary.json` |
| 10 | R2 "Our cohort" percentile; Fig 2f stars | our relabelings at the 70th–100th percentile | the percentiles use notebook 45's re-implementation of our cohort (whole-PT GSEA relabeled 89 and 266; step-screen *species* calls 50), not Table S1/notebook 37 (193 and 255; 16). The 50 against 16 discrepancy is not reported | `anchor_our_cohort.csv` |
| 11 | Figure 4d | prints "robust > uncalled, Census: p = 2 × 10⁻¹³" | the response to M7 says the nominal, overlap-ignoring p-values are no longer quoted | `external_replication_summary.csv` |
| 12 | R3 classes | "7,407 probe-balanced genes … summed over S2 − S1 and S3c − early" | 7,526 probe-balanced genes, of which 7,407 carry a union class (119 untested). The class is a union, not a sum | `gene_classes_symmetric.csv` |

Further smaller wording problems:
- R1 says that reusing reads inflated statistics by "about 30%"; the table gives 34%.
- Figure 2d is still titled "Specimen differences are offsets" (v0.1 minor 10).
- R1's 91% sign agreement omits that segment labels reach 95% on the same 244 genes (Figure 1d).

---

## 2. Does each claim follow from its evidence?

### (a) "Largely species-specific beyond a conserved core" against the class counts

The two lines of evidence are consistent, but only because they measure different things, and the
paper never shows the link. Of the classifiable zonated genes, 51% are conserved. Among genes zonated
in **both** species, 72% (196 of 271) agree in direction, against 50% expected without association.
The low index is therefore not "most zonated genes disagree". It means that **agreement is concentrated
in the strongest gradients and absent from the large mass of weakly zonated genes.**

**[Reviewer check]** I re-ran notebook 43's gradients and binned genes by |mouse gradient|.

| Contrast | Bottom 80% of genes: raw r (sign agreement) | Top 1% of genes: raw r (sign agreement) |
|---|---|---|
| S2 − S1 | −0.02 to −0.03 (0.47–0.51) | 0.61 (0.74) |
| S3 − early | about 0 (0.57–0.64) | 0.66 (0.82) |
| S3c − early | about 0 | 0.55 (0.77) |

This curve is the direct evidence for the title. It belongs in Figure 2 and should replace the
unsupported "mostly apical transporters" (mismatch 1).

The transporter claim also conflates power with conservation. Transporters are strongly expressed and
strongly zonated, so 40 of 43 are called zonated in both species (93%), against 271 of 388 overall.
Among genes zonated in both species, direction agreement is 33 of 40 (83%) for transporters against
196 of 271 (72%) overall. Report this comparison, not "77% against 51%".

"Species-only and reversal genes are few" is not right for reversals: 75 reversals are 28% of the
genes zonated in both species. They are also label-sensitive in S2 − S1. With Lake-transferred labels,
S2 − S1 reversals fall from 37 to 24 and 15, and halve in the cross-fit (21 → 11/9 and 15 → 8/5).
S3c − early reversals are stable. State this.

**The species-only asymmetry is a threshold choice, and the paper uses the same quantity in two
incompatible ways.**

| Rule | Mouse-only : human-only |
|---|---|
| Notebook 40 (v1) | 247 : 66 |
| A8 alone (external flat confirmation) | 187 : 31 |
| A3b alone (amplitude-matched thresholds) | 59 : 141 |
| Symmetric (A3b + A8) | 41 : 76 |

The "no skew" result exists only because A3b scales the human floor and flat margin by the human/mouse
amplitude ratio (0.58). That is correct if the lower human amplitude is technical: capture, diffusion,
the 76% bin coverage or label impurity. It is wrong if the lower amplitude is biological, and the paper
reports it as biology ("Human gradients are smaller as well as different"). The answer to "which genes
are zonated only in mouse" therefore depends on an unresolved question.

The pre-specified cross-fit criterion for "unskewed" also failed (mismatch 3). The defensible statement
is: "the balance of species-only genes depends on whether zonation is judged on an absolute scale (more
mouse-only genes) or relative to each species' amplitude (no excess); either way confidently
species-only genes are few." Drop "not skewed toward mouse" from the Introduction and the Abstract, or
qualify it there.

### (b) Construction of the conservation index

- **Winsorisation moves the estimate. [Reviewer check: notebook 43's code, q varied.]** The index
  winsorises gradients at 1/99%, which clips precisely the strong conserved markers.

  | Contrast | No winsorisation | Published (1/99%) | q = 0.05 |
  |---|---|---|---|
  | S2 − S1 | 0.25 (0.19–0.32) | 0.17 | 0.10 |
  | S3 − early | 0.35 (0.27–0.42) | 0.28 | 0.19 |
  | S3c − early | 0.34 (0.27–0.40) | 0.29 | 0.21 |

  External-only values span 0.25–0.31. Every upper bound stays below 0.5, so the pre-registered wording
  survives, but the point estimate shifts by up to 50%. Report the range.

  The disattenuation is also mismatched: r is computed on winsorised values, but reliability on
  unwinsorised ones (`zonation_reliability.disattenuated`).
- **The ceilings.**
  - The S2 − S1 human ceiling rests on a single dataset pair (ours against Lake). Census human has no
    S2, and Lake S2 comes from few nuclei. The "bars are medians" in Fig 3b therefore summarise one
    point.
  - The ceilings compare across platforms (Visium HD against snRNA), while the numerator compares
    within one platform. This probably biases the index *upward*, which is conservative for the claim.
    Say so.
  - A triad check on S3 − early, which gives each dataset's validity from three same-species datasets,
    gives a cross-species correlation of about 0.24–0.31 **[reviewer check]**, consistent with the
    index.
- **Intervals.** The 95% intervals are gene bootstraps. They treat genes as independent, although
  gradients are module-correlated, and they condition on the specimens.
  - Methods 6.23 states the second point. The first should be stated too.
  - A module- or pathway-block bootstrap, and a donor bootstrap of Lake and mouse snRNA for the
    external index, are cheap and would give honest intervals.
- **No evolutionary benchmark.** Rat failed the reliability rule, so nothing says whether 0.2–0.35 is
  unusual for two mammals about 90 My apart. v0.1's M2 alternative ("PT zonation is evolutionarily
  labile in general") remains open and belongs in the Discussion.

### (c) R2: conditional specificity and what it licenses for the smooth coordinate

- **External calibration covers the segment-step screen only.** The specificity of the smooth
  coordinate rests on our cohort's two relabelings and on AKI, where under the paper's own rule
  T_spatial *passes* (mismatch 4) despite a known positional nuisance (coherence 0.42).
  - The honest conclusion is that the rule NCDR < 0.25 cannot detect positional non-specificity of
    this size.
  - Our cohort's 0.02 is reassuring because of its margin, not because it passes a rule.
- **The robust list is driven by within-segment position.** Of the 36 robust pathways, **35 are among
  the 79 calls of the "smooth position beyond S1/S2/S3 steps" joint test** (6 and 0 relabeled). Only
  **9 of 36 are among the 16 step-screen calls** (`matched_and_joint_tests_all_partitions.csv`).
  - R1 says within-segment order is weakly supported: P1w 0.17 and 0.10, gene-fold agreement 0.34–0.70,
    human S1 order depth-dependent. R1 also says the coordinate adds "not resolution finer than
    segments".
  - Both cannot stand unqualified. The equal-depth rerun (32 of 36 robust) mitigates the depth concern
    but does not remove the tension.
  - Report the decomposition. Present the 9 step-screen pathways as the externally calibrated top
    tier.
- **Our cohort's within-species positional coherence is atypically low** (0.013, against 0.14–0.20 in
  the human pools). The human pair is two sections of one donor, so the relabelings under-represent
  donor-level positional variation. The external human nulls compensate for the step screen, not for
  T_spatial.
- **"Not specific" for average-level screens means vulnerable to relabeling, not wrong.** In
  cross-species draws, 71% of whole-PT GSEA species calls recur in the all-donor reference (base rate
  3%), against 81% for the step screen. R2 and Methods 6.14 say this, but the Abstract and the R2
  heading do not.

### (d) R5 biology against the symmetric classes

| Claim | Problem | Fix |
|---|---|---|
| Glutathione synthesis is late in mouse, with "weaker or no" human gradients and "not a different placement" | Gss is conserved and Lake GCLM is opposite (mismatch 8). The pathway is "shared" and probe-sensitive | Remove glutathione from the Discussion's short list and from "Implications: metabolism" |
| Cyp2e1 is a "rodent placement with weaker or no human gradient" | Human CYP2E1 is not expressed (1.6% detection; Lake ≤ 3 CPM), so this is an expression difference, not a positional one | Re-label it as mouse-specific expression |
| Acaa2 reversal (Discussion 3) | With reviewed labels it is a reversal only in S2 − S1, the contrast in which 65% of human "S2" transfers to Lake S3; our human S3 − early is +0.04. With full-gene transferred labels it becomes an S3c − early reversal, which is reassuring (`human_lake_centred__*/gene_classes_symmetric.csv`). The Visium donors show ACAA2 falling (0 of 7 agree with Lake) | Keep the claim with the label caveat, and report the Visium value or drop Visium altogether (see F2) |
| "Weak or no counterpart" for β-oxidation and sterol synthesis (Implications) | Acadm and Hmgcr are conserved | Use "weaker" |
| Hsd17b4 is "the only confidently mouse-only" peroxisomal gene | It is indeterminate under transferred labels | Add the caveat |
| Gatm: human "extends along the PT" | Acceptable as an amplitude claim, but Census human (−0.10, v0.1 minor 22) and the Visium donors (−0.95 in 10 of 10 samples) are omitted | Add both |
| Dcxr: "Visium donors cannot test a rising gene" | RBP4's *rise* is cited from the same donors | See F2 |

### (e) Does the coordinate choice affect a headline claim?

**R3 and R5 do not depend on it**: they use labels. These do:
- the robust and core lists: 36/26 under scFates against 26/15 under DPT13; 12 of 26 core pathways
  shared;
- "26 of 36 summarise mouse zonation", which becomes 15 of 26 under DPT13;
- the within-mouse positional calls: 16 against 37;
- **peak positions**: the human-high median peak is 0.16 under scFates and 0.44 under DPT13 (the
  mouse-high median is 0.63 and 0.68).
  - In segment terms, human-high peaks still fall in human S1 (DPT13 S1→S2 at 0.48–0.50).
  - Under SCF13-ED, peaks correlate with scFates at only ρ 0.54 (unreported), and 6 of 32 move by more
    than 0.25. Several jump from about 0.6–0.75 to 0.01 (Vitamin digestion, Fat-soluble vitamins,
    Nitrogen metabolism).
  - Report peaks in segment units only.

The disclosure is now honest. But the pre-registered rule selected DPT13, and on the primary joint test
DPT13 passes the gate. Keeping scFates is therefore a forking path that a referee will not accept as
"the coordinate on which the analysis was designed". The simplest fix: make the **coordinate-robust
intersection** (22 robust; 13 mouse summaries, 3 shared, 0 human-led, 1 mixed, 5 unresolved) the
primary reporting set, and report each coordinate separately as a sensitivity. The qualitative message
does not change.

---

## 3. Novelty and framing

- **What the paper is now.** It is a biology paper (R3, R5) with a methods component (R2) and a
  reconstruction (R1). The order puts the reconstruction first and the methods second, and the biology
  uses neither.
- **The outline's central claim is not delivered.** "Reconstructing a continuous PT coordinate … lets a
  human–mouse comparison ask *where*" is contradicted by R1 itself: no finer-than-segment resolution,
  and the main claim is coordinate-free.
- **Novelty risk.** The main claim is reproducible from public atlases alone (external index
  0.27–0.28). The data's added value is real but modest: intact tissue on one platform, region
  matching, and replication within specimens.
- **Single most important message:** *measured against same-species reproducibility, human and mouse
  PT share the axial zonation of their most strongly zonated genes and little else, and this holds in
  public atlases alone.* Hence mouse PT is an incomplete guide to human PT outside the canonical
  segment markers.
- **The Abstract does not deliver it efficiently.**
  - It opens with 12,866 structures on a coordinate that the main result does not use.
  - It spends its second paragraph on R2.
  - It quotes a threshold-dependent asymmetry (41 against 76) as a finding.
  - It never shows that the conserved part is the strong-gradient part.
- **Recommendation.** Frame the paper as a biology and analysis paper. Lead with the ceiling-calibrated
  comparison and the strength-dependent agreement curve. Present the few-specimen specificity result
  as a secondary, general methodological finding (one paragraph and one figure). Move the coordinate to
  a Supplementary Note, as a tool that adds within-specimen replication. Title, for example:
  *"Human and mouse proximal tubule share the zonation of their strongest segment markers and little
  else."*

---

## 4. Length and structure: concrete proposal (Results about 4,500 words, 6 main figures)

Current length: Results 7,569 words (R1 1,599; R2 1,326; R3 1,937; R4 725; R5 1,979); Discussion
1,551; Methods 9,230.

**Results**

| New section | Content (source) | Words | Main figure |
|---|---|---|---|
| R1 Data and segment labels | Cohort, segmentation and matrices; joint labels with distance, attachment and reference-transfer validation; spillover in one sentence (current R1 cohort + label block) | 600 | Fig 1: design (1a), label validation (S10a–c condensed), one coordinate overview panel (1b right) |
| R2 Conservation against same-species ceilings | Ceilings, index with the winsorisation range, agreement by gradient strength (new panel), cross-fit labels, external-only index, amplitude spread (current R3 first half, S9, S4h) | 1,000 | Fig 2: 3b, strength curve, S9 summary, spread |
| R3 Which genes are shared | Symmetric classes; direction agreement among genes zonated in both species; transporter enrichment stated correctly; asymmetry depends on the threshold; tissue and physiological state in one paragraph (R3 second half, notebook 44) | 900 | Fig 3: 3a, 3c, class by gradient strength |
| R4 Pathway screens with two specimens per species | Our cohort's table (condensed to 4 strategies), many-draw nulls, precision, AKI with the corrected rule reading, mechanism in two sentences (R2) | 800 | Fig 4: 2a, condensed 2b, 2f, 2g, 2h |
| R5 Position-dependent pathways summarise mouse zonation | Coordinate-robust intersection as the primary list; step-screen tier; categories; robustness and external agreement in one paragraph (R3 pathway part + R4) | 550 | Fig 5: 3d + condensed 4a + 4c |
| R6 Lead genes | Gatm, Acaa2, Dcxr, Ugt3a1, RBP4, with the Slc7a13 reference. Glutathione, Cyp2e1, sterol, OAT and "what did not hold" go to a supplementary table (R5) | 650 | Fig 6: Fig 5 cards, Gamt replaced by Rbp4 or UGT1A9, with a rat-protein column |

**Move to Supplementary Notes.**
- SN1: coordinate construction, validation, limitations, choice and the DPT13 disclosure table (R1
  "Coordinate", "does not resolve", "Choice", "adds"; Methods 6.9–6.10; Figs S1, S2).
- SN2: decoy calibration, the pathway layer, probe panels and the split-plot test (R2 decoy and
  pathway-layer paragraphs; R4 split-plot; Figs 2c–e, 4e).
- SN3: history of the class rules (v1 → symmetric) in one table (replaces Fig S8 and all "withdrawn"
  paragraphs).
- SN4: physiological and tissue state (notebook 44; Fig S11).
- SN5: external agreement details and genome-wide Lake (R4 external; Fig S7).
- SN6: Visium donors (the full lead-gene table) and rat proteome.
- SN7: segmentation model (Methods 6.3).
- SN8: **one table of deviations from protocol**, which replaces the scattered "(disclosed)" text and
  the seven-item list.

**Cut.**
- The earlier 82 → 66 → 42 list and notebook 12's historical DPT row (Figure 4a).
- "Earlier drafts", "v0.2 call withdrawn" and every other narration of internal history.
- Notebook, workstream and internal arm names (SCF13-ED, DPT13, CAC, P1w, A3b, A8) from the main text.
  Add a short glossary for NCDR, T_spatial, joint test and S3c.
- Figure S2 (the matched-only DPT diagnostic, superseded by the D4 reruns).

**Merge.**
- R4 robustness into R5.
- The tissue-state checks into R3.
- The two label-validation passages (R1 and Methods 6.25) into one.

**Methods.** About 3,000 words in the main text: cohort and ethics, Visium HD, segmentation summary,
matrices, labels, gradients and index, classes, pathway screens and specificity, external data,
statistics. The rest goes to Supplementary Methods. The figure legends in the draft (about 4,100
words) also need about a 50% cut, with source-file provenance moved to a manifest.

---

## 5. Ranked issues

### Fatal (fix before submission)

**F1. Mandatory records are missing, and one of them bears on a result.**
- *What is missing.*
  - Ethics approvals for the animal and human work.
  - Human donor source, consent, age, kidney function, procurement, and whether the two sections are
    serial.
  - Mouse strain, age and supplier.
  - The H&E–Visium registration method; provenance of the human `v2` segmentation; the pixel scale.
  - Data and code deposition.
- *Why it matters beyond paperwork.* Notebook 45 shows that positional screens fail when mice differ in
  (pubertal) age. R2's claim that "our two control mice" share positional biology therefore depends on
  an unrecorded age.
- *Fix.* Lab records; no analysis. If the ages differ, R2's cohort conclusion must be revisited.

**F2. Failed-validation Visium donors are used selectively.**
- Abedini et al. (2024) failed the pre-specified S3 label check. The draft still cites those samples as
  replication for UGT1A9, UGT3A1, ACOX2, SLC6A19 and the *rise* of RBP4.
- It declares them unable to test a rising gene for DCXR. It omits that GATM, ACAA2 and DCXR fall there,
  in the mouse direction (`C_visium_lead_genes.csv`, `C_posthoc_centred_lead_genes.csv`).
- *Fix.* Either drop every Visium lead-gene statement, or report the full lead-gene table in SN6 with
  the label-failure caveat, applied identically to every gene. Label 7 of 8 as post hoc. No new
  analysis.

**F3. The title-level claim is not shown, and its descriptor is false.**
- "Beyond a conserved core of strong markers, mostly apical transporters" is never demonstrated. The
  strength-dependent agreement is not in the paper, and transporters are 17% of the conserved class.
- *Fix.* Add the agreement-by-gradient-strength analysis (a small analysis on existing objects; my
  version is in `SCRATCH/review/nb43_check/`). Reword to "a conserved component concentrated among the
  most strongly zonated genes, in which transporters are over-represented".

### Major

**M1. The species-only asymmetry depends on thresholds.**
- Mouse-only to human-only runs 247:66, 187:31, 59:141 and 41:76 across rules. The pre-specified
  cross-fit criterion failed, and A3b presupposes that the human amplitude loss is technical.
- *Fix.* Remove "not skewed toward mouse" from the Abstract and Introduction, or state the dependence.
  Add a paragraph distinguishing absolute from relative zonation. Report the failed criterion. No
  analysis.

**M2. Sensitivity of the conservation index.**
- *Fix.*
  - Report the index under q ∈ {0, 0.01, 0.05} (0.10–0.35) and fix the reliability mismatch.
  - Add a module-block bootstrap, and a donor bootstrap for the external index.
  - State that the S2 − S1 human ceiling is a single pair, and give the likely direction of the
    cross-platform ceiling bias.
- *New analysis:* small. It would change the reported point estimates, probably not the wording.

**M3. The pathway signal is within-segment, and the coordinate is not validated there.**
- *Fix.*
  - Report that 35 of 36 robust pathways are called by the beyond-steps statistic and 9 of 36 by steps.
  - Present the 9 step-screen pathways as the top tier.
  - Express peaks in segment units and report the SCF13-ED peak correlation (0.54).
- No new analysis.

**M4. The coordinate decision is a forking path.**
- *Fix.* Make the scFates ∩ DPT13 robust set (22 pathways) primary, and report all three coordinates
  as sensitivities. Abstract and R3 numbers become 22 pathways (13/3/0/1/5). No new analysis.

**M5. AKI specificity is misread (mismatch 4).**
- *Fix.* State that T_spatial passes the paper's R2 rule and fails only notebook 45's stricter
  criterion. Conclude that the R2 rule is insensitive to a known positional nuisance. Note that average-level joint screens
  had the lowest AKI ratios (0.09). No analysis.

**M6. The R2 headline overstates "not specific".**
- *Fix.* Rephrase the heading and the Abstract as "cannot be shown specific by relabeling". Give the
  71% against 81% precision in the main text. No analysis.

**M7. R5 claims and the Discussion's short list.**
- *Fix.* Remove glutathione and Cyp2e1 from the positional list; use "weaker", not "no", where
  Acadm, Hmgcr and Gss are conserved; add the S2 − S1 label caveat for Acaa2 and the label dependence
  of Hsd17b4; add Census and Visium values for GATM. Section 2(d) gives the details. No analysis.

**M8. Framing and structure (section 3).** Rewrite the Abstract around the single message; restructure
as in section 4; update the outline's central claim. No analysis.

**M9. Sufficiency of the two running addenda.**
- *Addendum 1.*
  - **Index without the 611 integration genes.** Useful, but it is weaker evidence than the cross-fit
    already done, and it will probably *lower* the index (the conserved class is enriched for
    integration genes), so it supports the claim.
  - **Integration-gene enrichment.** Must be restated for the **symmetric** classes; the current
    "16% against 2%" uses superseded notebook 40 classes.
  - **Transporter class shares.** Should be reported as direction agreement among genes zonated in
    both species, as in 2(a).
  - **Table S2 coordinate columns.** Needed.
- *Addendum 2.*
  - **Overlap-respecting replication test.** Fine, but it changes no conclusion unless it is positive,
    because pathway-level replication is already disclaimed.
  - **Symmetric-class enrichments.** Needed before any species-only biology is named. With 41 and 76
    genes, a null result should be pre-declared as "no pathway naming", and Fig S8d should then go.
- *What the addenda miss.* Neither addresses M2's winsorisation and bootstrap issue, M3's
  decomposition or F3's strength curve. These three are the analyses that bear on conclusions.

### Minor

1. Mismatches 2, 5, 6, 7, 9, 10, 11 and 12 (table above).
2. Unresolved minor items from v0.1:
   - item 17: mouse S3c calls are confirmed against references whose S3 is mostly outer stripe;
   - item 20: conventional-only pathways replicate better as a group (Mann–Whitney p = 7 × 10⁻⁴);
   - item 22: Census human GATM −0.10;
   - item 37: jargon density, which remains high (NCDR appears on 34 lines; DPT13 and SCF13-ED on another 34).
3. Female-mouse references replicate only 8 and 4 of the 36 robust pathways. Report this beside 17 and
   6.
4. Figure 4a uses fill for both core status and DPT13 robustness, and still shows the historical DPT
   row.
5. Figure 5: replace Gamt, which is not replicated, indeterminate and expressed at log-norm < 0.7, with
   RBP4 or UGT1A9. Consider Gclm in place of Gclc only if glutathione stays in the paper.
6. Figure 1b: the coordinate trace loops back on itself in the Harmony 1–2 projection. Say that the
   curve lives in five dimensions, or show a projection in which it is monotone.
7. Figure 2b has 13 rows; keep 4–5 in the main text.
8. Figure 3b: "Bars are medians" summarises a single within-human S2 − S1 point.
9. The [VERIFY]/[CITE] items:
   - UGT1A9 exon-1 probes;
   - the AI-agent wording for the literature search, with archived queries;
   - OpenMidnight (model card only);
   - citations for GSE267280 and GSE277302;
   - whether to report the segmentation PQ;
   - the θ = 6 rationale.
10. The Kim et al. (2011) citation for S3 vulnerability is still unchanged (v0.1 item 35).
11. The Discussion's "first comparison" claim needs a final check against Klötzer et al. (2025), whose
    PT depth was not verified, and against human spatial atlases with PT subsegment labels. A Europe
    PMC search (generic terms) found no position-by-position human–mouse PT comparison.

---

## 6. Status of the v0.1 issues

| v0.1 item | Status in v0.3 |
|---|---|
| M1 (labels) | Largely resolved for the index (cross-fit). Not resolved for S2 − S1-specific gene claims (reversals halve with transferred labels; Acaa2, Slc9a3 and Hsd17b4 change class). Medullary-ray location of human S3 is untested |
| M2 (ceiling, rat) | Resolved in design. New sensitivities in 2(b). The evolutionary-lability alternative is still open |
| M3 (flat calls) | Mechanics resolved (A8). The A3b dependence is new (M1 above) |
| M4 (confounders) | Resolved as far as the data allow |
| M5 (coordinate) | Disclosure resolved. The decision itself is not (M4 above) |
| M6 (two draws) | Resolved for the step screen. Not resolved for T_spatial; the AKI reading is wrong (M5) |
| M7 (pathway layer) | Reframing resolved. The overlap-respecting test is pending (addendum 2); Figure 4d still prints the nominal p |
| M8 (citations) | Resolved. A new overstatement replaces the old one (mismatch 1) |
| Minor 36 (structure) and 37 (jargon) | Not resolved (sections 3 and 4) |
| Minor 17, 20, 22, 25, 34 and 35 | Open |
