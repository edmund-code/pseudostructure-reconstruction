# Notebook 48: overlap-respecting replication (M7 fix 4) and symmetric-class enrichments (D3)

Notebook `analysis/notebooks/48_pt_pathway_replication_overlap.ipynb`, logic version
`48.pathway_addenda.1`; results in `results/pt_revision_addenda_pathways/`. Protocol in
the workstream-5 plan (round 2), written before any statistic was computed.

## (a) Replication against overlap-preserving nulls

**Test.** The replication score is notebook 37's per-pathway joint z. Recomputing it with exact
matched-null moments reproduces the saved z (median |Δ| 0.009). A null draw relabels genes within the
matching strata and rescores every set at once, so set sizes, strata and every pairwise overlap are
kept. Each draw recomputes the matched z and the VIF. There were 9,999 draws per reference and
statistic.

| Census cortex vs male mice | Mean z | p vs matched random | Notes |
|---|---:|---:|---|
| 36 robust | 1.53 | 1 × 10⁻⁴ | null median 0.00 |
| 26 core | | 1 × 10⁻⁴ | |
| 17 robust programs | | 1 × 10⁻⁴ | |
| 21 relabel-called (not robust) | 0.59 | 0.027 | |
| Robust − relabel-called, pathways | 0.94 | **0.005** | 28-pathway comparator: 0.010 |
| Robust − relabel-called, programs (17 vs 13) | 0.53 | 0.09 (gene relabeling); **0.059 (label permutation)** | |

- **Lake vs male mice.** Robust 1.24 (p 1 × 10⁻⁴); D p = 0.015; program label permutation p = 0.17.
- **Selection-matched null** (strata also split by our `T_spatial` quintile): p = 0.004 (Census) and
  0.035 (Lake); 0.06 and 0.19 with female mice; D p = 0.04–0.32.
- **Slope-free score** (uniform flattening removed): robust beats random (p 0.001 Census) but not
  relabel-called (D 0.07, p 0.42).

**Pre-specified decision: better than random only**, for both Census and Lake. The pathway-level
comparator condition passes. The program-level condition fails (p = 0.059 against a 0.05 threshold,
with no margin either way).

**Proposed R4 text.**

> Against 9,999 random gene-set collections that keep the 36 pathways' sizes, expression strata
> and pairwise overlaps, the robust pathways replicated in both human references (mean joint z 1.53
> and 1.24, null median 0; p = 10⁻⁴). They also replicated more than the 21 pathways called only under
> specimen relabeling (mean z 0.59; p = 0.005 and 0.015). That margin did not hold for the 17 programs
> as units (p = 0.06 and 0.17), nor once the uniform flattening of human gradients was removed
> (p = 0.42). It weakened against random genes as strongly selected on our `T_spatial`
> (p = 0.004 and 0.035). By the pre-specified rule this is replication better than chance, not
> pathway-specific replication: the robust pathways replicate mainly because they are where mouse
> zonation is strongest.

This replaces the Mann–Whitney p-values (10⁻¹³, 10⁻²⁰, 0.016, 0.033).

## (b) Class claims under the symmetric classes (notebook 40's code)

The guard passed: rerunning with notebook 40's classes reproduced all 5,680 rows of its saved table.
Symmetric classes have 222 conserved, 92 reversal, 51 mouse-only and 102 human-only genes (notebook
40: 175, 65, 285, 81). Enriched pathways: 23 conserved, 4 reversal, 7 mouse-only, 20 human-only (notebook
40: 30, 2, 17, 21).

| Draft claim (notebook 40) | Old | New (symmetric) | Decision |
|---|---|---|---|
| Lines 311–312: integration genes are 16% of conserved and 2% of mouse-only genes (611 PCA genes) | 16% / 2% | 13% / 0% (752 HVGs: 16% / 0%) | change numbers |
| Lines 613–614: mouse-only enriched for sterol synthesis | q 5 × 10⁻⁵ (Cyp51, Ebp, Lipa, Lss) | q 0.006, **1 member (Ebp)**; Cyp51, Lss, Lipa now indeterminate | keep by rule, but a single-gene result: reword |
| Mouse-only peroxisomal β-oxidation | q 3 × 10⁻⁴ (4 genes) | q 0.31; only Hsd17b4 remains mouse-only | change (drop) |
| Mouse-only tyrosine catabolism | q 0.002 (5 genes) | q 2 × 10⁻⁵ (Comt, Fah, Maoa); robust to expression | keep |
| Mouse-only pyrimidine catabolism (notebook 40 text) | q 0.007 | q 0.74 | change |
| New mouse-only enrichments | none | steroid-hormone metabolism (Hsd17b11, Stard3nl), histidine, arginine biosynthesis (2 genes each) | small; report only if needed |
| Conserved: amino-acid transport, bile salts/organic acids, bile secretion, renin–angiotensin, apical surface, semaphorin | q 10⁻¹⁶ to 0.002 | q 4 × 10⁻¹⁷ to 0.013; first four robust to expression | keep |
| Reversal: haematopoietic lineage (Cd55, Mme), viral myocarditis | q 2 × 10⁻⁴, 0.002 | q 0.011, 0.07; not robust to expression | keep (weak) |
| Human-only programs (Figure S8d) | 6 programs, 1 robust to expression | 3 of 6 kept, each carried by 1–2 genes; DAG/IP3, phospholipids and sphingolipid lost; new hits are mostly single-gene DNA-replication sets; 1 of 20 robust to expression | change: no robust human-only program |

**Suggested replacement for draft lines 613–614.**

> Under the symmetric classes, mouse-only genes are enriched for tyrosine catabolism (Comt, Fah,
> Maoa; q = 2 × 10⁻⁵, robust to expression matching). Sterol synthesis rests on Ebp alone and
> peroxisomal β-oxidation on Hsd17b4 alone, so neither is a program-level result. The conserved class
> remains enriched for amino-acid and organic-anion transport, bile secretion and the
> renin–angiotensin genes.

**Lines 311–312.** Replace with "13% of conserved genes against none of the 51 mouse-only genes".

Figure S8d/e remains notebook 40's by design (superseded figure).

## Round 2b (referee report v0.3: 2(c), M5, M6, M9)

**Declared naming rule.** The rule was added to the plan after the all-gene run had been seen, and
before the probe-balanced run. An enrichment counts only if q ≤ 0.10, it is not expression-driven, and
≥ 3 class genes carry it.

Probe-balanced classes (41 mouse-only, 76 human-only):
- mouse-only: **Tyrosine metabolism (Comt, Fah, Maoa)**, q 2 × 10⁻⁶;
- human-only: Cardiac Conduction (Ahcyl1, Hipk2, Wwtr1) and TP53 phosphorylation (Hipk2, Prkaa2, Taf4).
  These are generic labels on overlapping 3-gene sets, so name the genes, not the pathways.

The result is not null, so the rule names only these. Sterol synthesis (Ebp alone; q 0.08) and
peroxisomal β-oxidation (Hsd17b4 alone) are not named. Figure S8d shows notebook 40's v1 programs and
should go; if a panel is wanted, it shows only the counting pathways.

**AKI under both rules** (`aki_strategies_under_both_rules.csv`; 1,429 tested pathways):

| Strategy | Condition | Relabeled | NCDR | Relabeled ÷ tested | R2 rule | D6 rule |
|---|---:|---|---:|---:|:-:|:-:|
| Whole-PT GSEA | 649 | 208, 20 | 0.18 | 8.0% | fail | fail |
| Segment GSEA | 790 | 292, 207 | 0.32 | 17.5% | fail | fail |
| Offset, joint | 56 | 10, 0 | 0.09 | 0.3% | pass | fail |
| Steps, joint | 26 | 12, 6 | 0.35 | 0.6% | fail | fail |
| Whole-PT, joint rank | 394 | 70, 0 | 0.09 | 2.4% | pass | fail |
| Segment, joint rank | 447 | 81, 11 | 0.10 | 3.2% | pass | fail |
| `T_spatial`, joint | 50 | 10, 7 | 0.17 | 0.6% | **pass** | fail |

**Recommended sentence (checked against these numbers).**

> In the AKI design, the two AKI mice differ in positional injury (within-group specimen × segment
> coherence 0.42). Even so, `T_spatial` passed the paper's R2 rule (NCDR 0.17; 10 and 7 relabeled
> calls of 1,429 tested), as did the average-level joint tests (0.09–0.10); it failed only notebook 45's
> stricter criterion. The rule therefore cannot detect a positional nuisance of this size. Our cohort's
> specificity rests on its margin (`T_spatial` NCDR 0.02; within-species coherence 0.013), not on
> passing the rule.

One correction to the requested wording: the rule is not blind to everything here. The step screen
failed it (NCDR 0.35).

**50 against 16.** Both are joint tests of species × S1/S2/S3 steps beyond a species offset, on
identical structures, genes and pathways.
- Notebook 37's version (16 calls; 0 and 0 relabeled) adds the steps to a shared 6-df spline in the
  scFates coordinate, so it is not coordinate-free.
- Notebook 45's version (50; 0 and 0) uses segment intercepts as the shared shape, with no coordinate,
  exact matched moments and a segment-centred VIF.
- The 16 are almost a subset of the 50 (15 of 16). The 50 contain 28 of the 36 robust pathways
  (against 9 for the 16) and 47 of the 60 `T_spatial` candidates. The referee's "robust list is driven
  by within-segment position" therefore reflects the spline-based step design.
- Figure 2f–g's stars use notebook 45's version.
- R2 should read "the coordinate-free segment-step screen called 50 pathways against 0 and 0 under
  relabeling (notebook 45; 16 with notebook 37's spline-based design)".
