# R2 revision: many-draw specimen nulls and an AKI positive control (notebook 45)

Workstream 5 draft for `docs/results/pt-pathway-method-selection.md` (new section) and the R2 text.
Notebook `analysis/notebooks/45_pt_specimen_null_many_draws.ipynb`, logic version
`45.specimen_null.1`; results in `results/pt_revision_method/`. Protocol, priors and decision rules
were fixed before any statistic was computed (the workstream-5 plan). One post-hoc check is labelled as such.

## What was run

- **Pools.** Census mouse PT nuclei (dataset 25818bf7; 12 male and 12 female mice, 3 per age at 3 weeks,
  3–6 months, 12 months, 21 months); Lake/KPMP healthy cortex (7 donors, S1/S2/S3); Census human cortex
  (6 donors, convoluted PT and S3); Lake human + Census male mouse (cross-species); our cohort
  (notebook 37 structures); our AKI design (mouse_only_v5: Ctrl1A2, Ctrl1A4, IR2A2, IR2A4).
- **Draws.** 300 + 300 random same-species 2 + 2 splits (male, female mice), all 105 (Lake) and all 45
  (Census human) splits; 100 random 2 human + 2 mouse draws, each with the species split and its two
  relabelings; an all-donor reference (7 human against 12 mouse donors).
- **Strategies at segment resolution** (the atlases have no coordinate): O1 whole-PT pseudobulk + signed
  GSEA; O2 S1/S2/S3 pseudobulk + signed GSEA; O3 unit-level offset `T_level`, joint test; P1 unit-level
  group × segment steps `T_steps`, joint test. Same weights, specimen contrasts, matched rank-AUC,
  VIF and BH rule as notebook 37. Two declared substitutions, both checked on our cohort: moderated t
  (voom-style) for DESeq2, and the exact matched-null mean/SD for 9,999 random sets (recomputing
  notebook 37's species `T_spatial` test gives 59 joint calls against its 60; median |Δz| 0.010).

## Results

**Same-species nulls** (pathways called per split, BH q ≤ 0.05; ~1,300–1,420 tested; share of splits
with ≥ 1 call in brackets):

| Pool | O1 whole-PT + GSEA | O2 segments + GSEA | O3 offset, joint | P1 steps, joint |
|---|---|---|---|---|
| Census mouse, male (300) | 22, p90 66 (94%) | 38, p90 92 (98%) | 7, p90 28 (77%) | 11, p90 43 (80%) |
| Census mouse, female (300) | 40, p90 137 (96%) | 62, p90 165 (99%) | 1, p90 26 (55%) | 5, p90 21 (82%) |
| Lake human cortex (105) | 34, p90 125 (95%) | 41, p90 155 (97%) | 3, p90 83 (64%) | **0, max 9 (6%)** |
| Census human cortex (45) | 13, p90 39 (91%) | 12, p90 45 (93%) | 0, p90 17 (49%) | **0, max 0 (0%)** |

Values are medians. Per-pathway calibration (share of pathways at p ≤ 0.01; 1% expected): P1 1.0%
and 0.6% in the human pools, 3.8–4.8% in the mouse pools; average-level screens 2.4–7.5% everywhere.

**Cross-species draws** (100; medians):

| Strategy | Species calls | Relabeled calls | NCDR | Draws meeting nb31 rule | Precision vs all-donor reference (base rate) |
|---|---|---|---|---|---|
| O1 whole-PT + GSEA | 21 | 25 | 1.34 | 19% | 0.72 (0.030) |
| O2 segments + GSEA | 58.5 | 40 | 0.69 | 21% | 0.73 (0.058) |
| O3 offset, joint | 0 (calls in 4 draws) | 1 | – | 3% | – |
| P1 steps, joint | 15 | 0.5 | **0.05** | **80%** | 0.81 (0.019) |

**Pre-specified decision: R2 is general at segment resolution** (P1 NCDR 0.05 < 0.25 with 80% of
draws meeting the rule; O1 NCDR 1.34 ≥ 0.5; P1 90th percentile ≤ 5% of tested in 4 of 4 pools). The
80% condition is met exactly, with no margin.

**Our cohort with the same machinery** (species; relabel Ctrl1A2, Ctrl1A4): O1 158; 89, 266 (notebook 37
DESeq2: 158; 193, 255). O2 196; 273, 243. O3 0; 49, 5. P1 50; 0, 0 (notebook 37's spline-based
steps: 16; 0, 0). Our P1 relabelings sit at the median of the N5 relabel distribution; our
average-level relabelings at its 70th–100th percentile.

**Mechanism.** Median within-pathway correlation of member genes across specimens:

| Pool | Specimen offsets | Specimen × segment | Within-unit (VIF source) |
|---|---|---|---|
| Census mouse male / female | 0.26 / 0.40 | 0.06 / 0.03 | 0.006 / 0.006 |
| Lake / Census human | 0.32 / 0.50 | 0.14 / 0.20 | 0.014 / 0.014 |
| Our cohort (within species) | 0.17 | **0.013** | 0.010 |
| Our AKI design (within condition) | 0.07 | **0.42** | 0.021 |

Specimen offsets move pathways together far more than the unit-level residual correlations the VIF
uses, so no competitive test with a unit-level VIF can absorb them. Positional deviations are less
coherent, except where specimens differ in positional biology. Per pathway, null-call frequency
tracks coherence only weakly (Spearman −0.22 to 0.33).

**Where positional screens are not specific.**
- *Age in male mice.* Age-balanced male splits: P1 median 0 (26 splits). Post hoc (after the pooled
  results): P1 median 1, 22 and 62 null calls when the two pairs differ by 0, 1 or 2 three-week-old
  mice; O1 does not depend on it (24, 23.5, 7). Female mice: 2, 7, 15. Sex-biased PT programs emerge
  after 3 weeks of age, with outer-stripe patterns ([PMID 40259083](https://pubmed.ncbi.nlm.nih.gov/40259083/)).
  The step screen is detecting a real pre- versus post-pubertal positional difference between
  specimens; a 2 + 2 design cannot separate it from the group factor.
- *AKI (D6).* See below.

**D6 · AKI against control (2 + 2).** Calls for condition; relabelings (Ctrl1A2 + IR2A2, Ctrl1A2 + IR2A4):

| Strategy | Condition | Relabelings | Relabel ÷ condition |
|---|---|---|---|
| O1 whole-PT + GSEA | 649 | 208, 20 | 0.18 |
| O2 segments + GSEA | 790 | 292, 207 | 0.32 |
| O3 offset, joint | 56 | 10, 0 | 0.09 |
| P1 steps, joint | 26 | 12, 6 | 0.35 |
| `T_spatial`, joint (mouse-only PT DPT) | 50 | 10, 7 | 0.17 |

- **Detects known positional biology.** 13 of 14 injury markers rise most in S3 (median AKI − control
  0.45 in S1, 1.25 in S3, log scale). P1 and `T_spatial` call TNF-α/NF-κB, Hypoxia and EMT (3 injury
  Hallmarks; Mann–Whitney p 3 × 10⁻⁴ and 8 × 10⁻⁵), and 77–80% of their called pathways change most
  in S3.
- **Outer stripe not resolved.** Deep (outer-stripe-like) and cortical S3 rise equally (7 of 14
  markers higher in deep S3). Deep S3 structures are depleted in AKI sections (248 and 91 against 833
  and 826 in controls), so outer-stripe injury shows as loss of labelled S3, not as change within it.
- **Fails the pre-specified 5% relabeling condition** (0.35 and 0.17). Within the AKI group the two
  mice differ in injury severity (deep S3 structures: IR2A4 91, IR2A2 248), and that difference
  is positional (specimen × segment coherence 0.42). **D6 decision: fail** on specificity; pass on
  detection.

## What this supports, and what the split-plot argument can claim

Supported, at segment resolution: across 750 same-species and 100 cross-species 2 + 2 designs,
average-level competitive screens report pathways in nearly every null split and as many under
relabeling as for species; a step screen with the joint test reports almost none when specimens of a
group share their positional biology (adult human cortex; age-matched mice). Our cohort's relabelings
behave like the external null.

Not supported: that position-dependent screens are specific by design. They are specific when
specimen × position deviations are small and incoherent; they are not when specimens differ in
positional biology (puberty, injury severity). NCDR is also not a false-discovery proportion: in
cross-species draws, 72% of whole-PT GSEA species calls were also called with all 19 donors (base
rate 3%). A strong true effect dominates a competitive ranking, so relabeled counts overstate how many
of the species calls are specimen noise. NCDR measures vulnerability to specimen noise.

The split-plot argument can claim: (1) relabelings are exchangeable draws of the section-level null;
(2) a specimen-level offset test has 2 df, so per-pathway offset inference is impossible in this
design, whereas positional terms are replicated within specimens. It cannot predict the specificity of
unit-level competitive screens, which depends on two empirical quantities the df count does not
contain: the size and the within-pathway coherence of specimen × position deviations. It should be
stated as motivation, with the coherence check (within-group positional deviations) as the evidence
for our cohort: within-species coherence 0.013, within-mouse 16 and within-human 0 `T_spatial` calls
(notebook 31, matched test).

What transfers to `T_spatial`: the external test is about the statistic class (position-dependent,
competitive, joint test) at segment resolution. It does not test the coordinate itself (registration
noise; DPT's NCDR 0.40), smoothing, or anything shared by both human sections of one donor. In our
cohort steps and `T_spatial` behave alike (16; 0, 0 and 60; 2, 0), and in the AKI design (26; 12, 6 and
50; 10, 7), which is the only coordinate-based evidence beyond the original two draws.

## Proposed R2 text (≤ 300 words)

> **A specificity problem in cross-species pathway analysis.** With two specimens per species, a
> pathway screen can only be judged against specimen variation. We measured how often each strategy
> reports pathways when no group difference exists: in our cohort's two balanced relabelings, and in
> 750 same-species and 100 cross-species 2 + 2 designs drawn from public snRNA atlases with segment
> labels (24 mice; 13 human donors). The atlases have no coordinate, so this test compares average-level
> screens with segment-step position-dependent screens, not the smooth coordinate.
>
> Average-level competitive screens reported pathways in almost every null split (whole-PT pseudobulk
> GSEA: median 13–40 per split in four pools) and, in cross-species draws, as many under relabeling as
> for species (median NCDR 1.34). The segment-step screen with the joint test called none in adult human
> cortex (median 0 in 150 splits), few in female mice (median 5), and in cross-species draws a median of
> 15 species against 0.5 relabeled calls (NCDR 0.05); our cohort's step-screen relabelings (0 and 0)
> sit at the median of that distribution. NCDR measures vulnerability to specimen noise, not the false-discovery proportion: most
> average-level species calls in the draws (72%) recurred with all donors.
>
> This specificity is empirical, not a consequence of the design. Specimen offsets were far more
> coherent within pathways than specimen-by-segment deviations (median member correlation 0.26–0.50
> against 0.03–0.20), which unit-level variance inflation cannot absorb. Where specimens differ in
> positional biology, positional screens also respond to relabeling: male mice split unevenly by
> prepubertal age (median 22–62 calls), and AKI mice of unequal injury severity (relabeled calls 17–35%
> of the AKI-versus-control calls, although the screen placed injury programs in S3). The split-plot
> argument motivates the approach; in our cohort, specificity rests on small within-species positional
> deviations (coherence 0.013), which we report beside every call.

(294 words.)

## Methods paragraph (for Methods 6.14)

Many-draw specimen nulls (notebook 45). PT nuclei with author segment labels were taken from CELLxGENE
Census 2025-11-08 (mouse dataset 25818bf7, embryonic and newborn stages removed; human cortex dataset
09b518f9) and from Lake et al. 2023 (KPMP v1.5; healthy cortex), keeping donors with ≥ 50 nuclei in
every segment. Counts over the accepted ortholog panel were normalized as for our structures. Within
each pool, tested genes (≥ 2% detection, equal-donor mean), matching strata, pathway sets (notebook
37's 1,513, 10–300 tested members) and a CAMERA variance inflation from within-cell residual
correlations were fixed once. For each 2 + 2 partition, unit-level weighted partial F statistics for a
group offset and for group × segment steps (notebook 37's designs with segment intercepts as the
shared shape, computed exactly from cell means and within-cell sums of squares) were tested with the
covariate-matched rank-AUC (exact null moments) and the joint test; whole-PT and per-segment
pseudobulk moderated t statistics (voom-style weights, empirical-Bayes variances) were tested with
signed preranked GSEA. Calls are BH q ≤ 0.05 within partition and strategy. Same-species splits are
nulls; in cross-species draws the species split is compared with its two balanced relabelings (species
segment shape as nuisance) and with an all-donor reference. The AKI design (mouse_only_v5) adds
notebook 37's coordinate `T_spatial` on the mouse-only PT DPT.

## Changes proposed to existing files (not made)

- `docs/results/pt-pathway-method-selection.md`: add this section after "R2: specificity"; replace
  "Under independent specimen deviations, each relabeling is one draw" paragraph's implication of
  specificity-by-design with the conditional statement above; add NCDR ≠ FDP.
- `docs/paper/outline.md` R2 row: evidence "Notebooks 31, 37, 45"; claim as in the proposed text.
- Draft/abstract: remove "as the split-plot design predicts" (referee 38); Discussion point 1 to say
  specificity is empirical and conditional.
- Commit `ws5/fetch_census_pt_cells.py` as `analysis/scripts/fetch_census_pt_cells.py` (writes
  `data/external/census_pt_cells/`, input to notebook 45).
