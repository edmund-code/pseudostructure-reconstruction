# Workstream B · Phase 1: pathway analysis along the PT coordinate, continuous against discrete

Written 2026-10-05 for the coordinating session. Phase 1 uses saved tables and stage-cache payloads
only. No model was fitted to private data. The scripts that produced every new number are in
`workstream-B scratch: q1_beyond_steps.py`, `q2_segment_ora.py` and `q4_curve_sizing.py`, with their outputs
saved as `.out` files beside them. The gene-score, overlap and resolution summaries (`q3_gene_scores.out`,
`q5_overlaps.out`, `q6_resolution.out`) came from inline commands on the same cached tables.
Coordinate: notebook 13 scFates (SCF13), with DPT13 and SCF13-ED as sensitivities, all from the
notebook-37 runs (`results/pt_pathway_final/{scfates,dpt13,scf13_ed}`).

## 0 · Summary

1. **Pathway calls do not separate continuous from discrete analysis.**
   - A coordinate-free S1/S2/S3 step screen (notebook 45) calls 20 of the 22 primary pathways. Six
     equal coordinate bins call all 22.
   - Eight of the 36 robust scFates pathways are missed by the step screen. All eight are near
     misses (step q 0.057–0.136).
   - The two primary pathways it misses (Cornified Envelope; xenobiotics by P450) are the two
     probe-sensitive ones.
   - "More pathways" cannot be the claim.
2. **New in phase 1: thresholded ORA passes the relabeling rule.** ORA on S1/S2/S3 DESeq2 genes calls
   118 pathways for species and 0 and 1 under relabeling (whole PT: 80; 0, 1).
   - This contradicts the unqualified R4 sentence that average-level screens "cannot be shown
     specific by relabeling".
   - The pass is partly built in. The relabeled DESeq2 design (`~species + group`) has 1 residual df
     against 2 for the species design, so the relabelings find only 0.5–4 DEGs per contrast.
   - ORA's species calls are offset calls, which probe-panel differences can create.
   - R4 needs narrowing to rank-based competitive screens, or an external, power-matched test of ORA
     (phase 2, notebook 56).
3. **The "beyond steps" signal is real at the scFates pathway level but coordinate-fragile.**
   - The joint test of smooth position beyond S1/S2/S3 steps calls 79 pathways on scFates (6, 0
     relabeled), 41 on DPT13 (7, 6) and 81 on SCF13-ED (10, 9).
   - Of the robust 36, 35 are called beyond steps on scFates, 18 on DPT13 and 18 on all three
     coordinates. Of the primary 22, the counts are 22 and 12.
   - At gene level the beyond-step score agrees across coordinates (Spearman 0.84 with DPT13, 0.79
     with SCF13-ED). It carries about 35–39% of the positional species-difference sum of squares in the
     top 10% of genes.
   - Its signal-to-noise against relabeling is the worst of the positional statistics. The 99th
     percentile of species scores over relabeled scores is 8× for beyond steps, 16× for steps and 43×
     for T_spatial.
   - The "beyond label steps" term cannot be separated from label noise with the existing tables:
     6–15% of structures per specimen sit out of order relative to the best monotone thresholds on the
     coordinate.
4. **Resolution budget.** The gene-fold positional noise of the coordinate (a lower bound) is about the
   same size as the within-segment spread of positions. The ratio is 0.74–0.77 in mouse S2/S3 and
   1.0–1.4 elsewhere. For the full coordinate, the noise SD is about 0.17 segment units in mouse
   S2/S3 and 0.24–0.33 elsewhere.
5. **The strongest "continuous-only" candidate is the ordering of switch-on positions inside mouse S2.**
   - 202 mouse genes switch on inside S2 under both coordinates. Their order agrees between scFates
     and DPT13 at Spearman 0.90 (0.69 for the 91 genes away from the boundaries).
   - The discrete comparator recovers part of it. The segment "switch fraction" r = (S2 − S1)/(S3 − S1)
     predicts the continuous order at ρ 0.69, and r itself replicates between the two mice at 0.95.
   - The continuous gain is therefore the part of the onset order that r does not explain. Whether it
     replicates across specimens and in external snRNA is untested; it is the first test of phase 2.
   - Interior peaks are the weakest category. In mouse, 360 genes have an interior extremum on scFates
     against 51 on DPT13.
6. **Recommendation.**
   - Keep the T_spatial joint test as the call-making screen.
   - Report every called pathway with its decomposition: steps, from the coordinate-free and externally
     calibrated screen; beyond-steps; and location descriptors in segment units, each with a
     replication flag.
   - Frame continuous analysis as "more information per program at a measured resolution", not as more
     pathways (framing F2 in section 6). Workstream A's coordinate enters through the resolution
     budget.

## 1 · What was read and what exists

Read in full: RULES.md, AGENTS.md, CLAUDE.md, `docs/paper/outline.md`,
`docs/results/pt-pathway-method-selection.md`, draft v0.6 R4, R5 and the Methods on pathway screens,
and Supplementary Notes 1, 2, 5 and 8. I also read the notebook-37, 45, 47, 48 and 50 tables, the
notebook-12 conventional tables, the notebook-07 module tables (built on 03's DPT, not scFates) and the
notebook-10 summary (on 03's DPT).

Prior "continuous against conventional" attempts, and why they do not settle the question:
- **Notebook 10** (03's DPT): 2 strict trajectory-only genes; 0 within-bin peak shifts and 0
  near-boundary cases among strict gene candidates. This is a prior negative at gene level.
- **Notebook 14**: the headline "66 discrete → 89 continuous" uses `T_total`, which is not specific
  (NCDR 1.72). It must not carry the claim; the method-selection doc already says so.
- **Notebook 07** (03's DPT): curve-module clustering. Response modules are unstable when one specimen
  is left out (ARI 0.39–0.53 against the primary solution). Specimen pairs reproduce a gene's
  cross-species state transition in 0.49 of pairs on average. No relabeling control was run.
- **Notebook 12** (scFates): the "information classes" framework and peak, width and sign-change
  descriptors. These are descriptive and unreplicated.

## 2 · Goal 1: taxonomy of pathway approaches and comparison

### 2.1 What each approach asks and returns

| Family | Approach | Question it answers | Returns |
|---|---|---|---|
| Discrete, average | Whole-PT pseudobulk DE + GSEA or ORA | Does the program differ on average between species? | yes/no, one sign |
| Discrete, segment | S1/S2/S3 pseudobulk DE + GSEA or ORA | Does it differ within S1, S2 or S3? | yes/no and sign per segment (3-bin location) |
| Discrete, segment | Segment markers (segment vs rest) + ORA | What characterises each segment, species pooled? | segment-specific programs (a zonation question, not a species one) |
| Discrete, positional | Species × S1/S2/S3 steps beyond an offset, joint test | Does the segment profile differ in shape between species? | yes/no; step pattern (2 df) |
| Discretised coordinate | Species × 6 equal coordinate bins | The same question, with 6 bins | yes/no; 5-df profile |
| Continuous | Gene GAM species × smooth position (T_spatial) + matched rank-AUC (joint test) | Is the positional shape of member genes more species-dependent than matched non-members? | yes/no; D(s) and S(s) curves; peak; shape class |
| Continuous | T_spatial + preranked GSEA or size-only sets | The same question, with a different set null | yes/no |
| Continuous | Smooth position beyond S1/S2/S3 steps | Does shape differ beyond what segment steps explain? | yes/no; within-segment component |
| Continuous, average | Gene GAM offset `T_level` or `T_total` | Average difference (plus shape) | yes/no, sign |
| Continuous, self-contained | Pathway score (module mean-z, AUCell- or decoupleR-like) GAM along the coordinate | Does the program's activity curve differ? | activity curve; yes/no |
| Continuous, specimen-level | Specimen × 8-bin pseudobulk split-plot (voom, moderated F) | Positional difference judged against specimen × position variation | yes/no, calibrated at gene level |
| Continuous | tradeSeq-style association, pattern/condition, startVsEnd and earlyDE tests | Position dependence per species; shape difference; end-point contrast; local difference | gene calls, then sets |
| Continuous | Curve-module clustering + ORA (notebook 07) | Which shapes recur, and which pathways fill them? | modules (shape families) with enriched pathways |
| Continuous | Switch-point, onset and peak-position enrichment | Where does each program switch or peak, and in which order? | location, sharpness, order |
| Continuous | Trajectory alignment (Genes2Genes; notebooks 09/10 on 03's DPT) | Where do species curves match or mismatch? | alignment path; mismatch regions |

### 2.2 Comparison table filled from existing results (species split of 2 + 2; 1,508–1,513 pathways)

NCDR = mean relabeled calls ÷ species calls. Rule = NCDR < 0.25 and each relabeling calls < 5% of
tested pathways. "∩22 / ∩36" are overlaps with the primary and robust lists. External = notebook-45
cross-species draws (median species / relabeled calls, median NCDR, share of draws meeting the rule),
where they exist.

| # | Approach | Species | Relabeled | NCDR | Rule | ∩22 / ∩36 | External | Source |
|---|---|---:|---|---:|:-:|---|---|---|
| D1 | Whole-PT DESeq2 + signed GSEA | 158 | 193, 255 | 1.42 | fail | 3 / 8 | 21 / 25, NCDR 1.34, 19% | nb37 S1; this phase (overlap) |
| D2 | S1/S2/S3 DESeq2 + signed GSEA | 204 | 343, 255 | 1.47 | fail | 7 / 17 | 58.5 / 40, 0.69, 21% | nb37 S1 |
| D3 | **S1/S2/S3 DESeq2 + ORA** (padj ≤ 0.05, \|LFC\| ≥ 1) | **118** | **0, 1** | **0.00** | **pass\*** | 9 / 18 | not run | **this phase** (q2) |
| D3b | Whole-PT DESeq2 + ORA | 80 | 0, 1 | 0.01 | pass\* | 6 / 13 | not run | this phase |
| D3c | S1/S2/S3 ORA without the LFC floor | 191 | 0, 0 | 0 | pass\* | 7 / 13 | not run | this phase |
| D4 | Segment markers + ORA (species pooled) | 9/6, 0/1, 11/10 (S1, S2, S3; up/down) | – | – | n/a | – | – | nb12 |
| D5 | Coordinate-free step screen, joint test | 50 | 0, 0 | 0 | pass | 20 / 28 | 15 / 0.5, 0.05, 80% | nb45 |
| D6 | Spline-based step screen, joint test | 16 | 0, 0 | 0 | pass | 8 / 9 | – | nb37 |
| D7 | 6 equal coordinate bins, joint test | 49 | 17, 0 | 0.17 | pass | 22 / 29 | – | nb37 |
| D8 | S1/S2/S3 pseudobulk joint rank test (notebook-45 unit-level) | 0 | 27, 3 | – | fail | – | 61 / 6.75, 0.09, 51% | nb45 |
| C1 | T_spatial, matched rank-AUC only | 116 | 17, 13 | 0.13 | pass | 22 / 36 | – | nb37 |
| **C2** | **T_spatial, joint test (primary)** | **60** | **2, 0** | **0.02** | **pass** | 22 / 36 | none (no external coordinate) | nb37 |
| C3 | T_spatial, preranked GSEA | 311 | 34, 24 | 0.09 | pass | – | – | nb37 |
| C4 | T_spatial, size-only random sets | 344 | 59, 55 | 0.17 | pass | – | – | nb37 |
| C5 | Beyond S1/S2/S3 steps, joint test (scFates / DPT13 / SCF13-ED) | 79 / 41 / 81 | 6, 0 / 7, 6 / 10, 9 | 0.04 / 0.16 / 0.12 | pass | 22 / 35 (scFates) | – | nb37; this phase |
| C6 | `T_level`, joint test | 0 | 46, 2 | – | fail | – | 0 / 1 | nb37 |
| C7 | `T_total`, joint test | 0 | 72, 15 | – | fail | – | – | nb37 |
| C8 | Pathway-score GAM, offset | 1,444 | 986, 1,317 | 0.80 | fail | – | – | nb31/37 |
| C9 | Pathway-score GAM, shape | 1,509 | 67, 495 | 0.19 | fail (33% of pathways) | – | – | nb31/37 |
| C10 | Specimen × 8-bin split plot, matched only (joint: 0) | 41 | 0, 0 | 0 | pass | 7 / 15 | – | nb37 |
| C11 | tradeSeq-style tests | not run | | | | | | needs fits |
| C12 | Curve modules + ORA | not on scFates; 03-DPT version unstable | | | | | | needs fits |
| C13 | Onset or peak enrichment | descriptive only (peaks: ρ 0.71 scFates vs DPT13; 7 of 22 change segment) | | | | | | needs fits |

\*D3 passes our relabeling, but the comparison is not power-matched: 1 residual df (relabeled design)
against 2 (species design). Its species calls are offset calls, exposed to probe-panel efficiency, which
relabeling cannot see. External same-species 2 + 2 draws (2 residual df) are the fair control (phase 2).

### 2.3 Overlaps between approaches (scFates)

- T_spatial joint (60) against the coordinate-free steps (50): 47 shared.
- Beyond-steps (79) against T_spatial: 53 shared. Against the coordinate-free steps: 43.
- Six bins (49) against T_spatial: 47 shared.
- S1/S2/S3 GSEA (204) shares 24 with T_spatial; segment ORA (118) shares 26.
- The robust 36:
  - 28 are called by the coordinate-free steps.
  - The 8 that are not are Metabolism of xenobiotics by P450, Valine/leucine/isoleucine
    degradation, Adipogenesis, Biological Oxidations, Diseases of Metabolism, Formation of Cornified
    Envelope, Water-soluble vitamins and Mitochondrial FA β-oxidation.
  - All 8 have step q between 0.057 and 0.136. Five of them are also beyond-step calls on both scFates
    and DPT13. Only two are in the primary 22, and both are probe-sensitive.

### 2.4 What exists and what needs new fits

- **Exists:** D1–D3 (D3 computed here from saved DESeq2), D5–D8, C1–C10, peak positions, shape
  classes, gene-level beyond-step scores on three coordinates, and LOSO stability for T_spatial.
- **Needs new fits (protocol in section 8):**
  - D3 under external, power-matched nulls and on equal-probe genes;
  - specimen-level pathway-score split plot (the self-contained alternative to C8/C9);
  - curve modules + ORA on scFates under relabeling;
  - onset and peak enrichment;
  - tradeSeq-equivalent tests inside the existing nested GAM (association within species,
    startVsEnd, earlyDE), reported as gene-rank concordance only;
  - per-specimen descriptors for the information categories.

### 2.5 Recommendation (Goal 1)

Calls come from **C2 (T_spatial joint test)**. It is the most specific competitive screen (2 and 0
relabeled), answers the positional question, and is protected against probe-panel offsets, which cancel
within species.

Each call is reported with three readouts:
- the **coordinate-free step screen (D5)** as its externally calibrated anchor, which says whether a
  segment-resolution difference exists without any coordinate;
- the **beyond-step component (C5)**, flagged when it holds on at least two coordinates;
- **location and order descriptors** in segment units (onset, peak), each with a specimen-replication
  flag.

D2/D3 are reported as the conventional comparators, with the ORA caveat. Self-contained scores (C8/C9)
are dropped unless the specimen-level version passes. Curve modules are kept only if they pass
relabeling and LOSO stability.

Why the T_spatial joint test and not the step screen as the primary: in this cohort they are close.
T_spatial is preferred because it is the statistic that returns location, but the paper should say
plainly that the step screen gives almost the same list.

## 3 · Goal 2: what continuous analysis can show that segment means cannot

### 3.0 Formal frame

Write a gene's positional profile in species k as f_k(s). The S1/S2/S3 pseudobulk is the projection
P f_k onto the 3-dimensional space of segment-wise constants, weighted by structure density. Two things
follow.

1. **Invisible.** Every within-segment function g with zero weighted mean in each segment satisfies
   P g = 0. Within-segment gradients, interior peaks and the location of a switch inside a segment
   produce no segment-mean signal unless they change a segment's mean.
2. **Blurred.** For a monotone gene, the three means carry one shape number, r = (m2 − m1)/(m3 − m1). A
   switch model has at least two shape numbers, location and width, and r confounds them: an early
   sharp switch and a gradual linear rise can give the same r. A transient of width w inside a segment
   of length L moves that segment mean by about a·w/L.

The species comparison follows the same algebra applied to Δ(s) = f_h − f_m. The beyond-steps
statistic tests (I − P)Δ ≠ 0 with label-defined segments.

Two artefacts can create (I − P) signal without biology:
- **(a) Label noise.** Structures labelled S1 that are really S2 carry the S2 level and sit at
  S2-like coordinate positions. 6–15% of structures per specimen are out of order relative to the best
  monotone thresholds.
- **(b) Coordinate noise.** A true step is smeared into a ramp of width about the positional noise.

Test A (beyond label steps) is protected against (b) but not (a). Test B (beyond coordinate-threshold
steps) is protected against (a) but not (b). A within-segment claim therefore needs both tests plus
replication (section 8, notebook 57).

### 3.1 Resolution budget of the coordinate (new, from notebook-35 fold noise)

| Species | Segment | Per-fold noise SD (segment units) | Noise ÷ within-segment position SD | Full-coordinate noise SD (÷√3, optimistic) |
|---|---|---:|---:|---:|
| mouse | S1 | 0.56 | 1.27 | 0.32 |
| mouse | S2 | 0.30 | 0.77 | 0.17 |
| mouse | S3 | 0.30 | 0.74 | 0.17 |
| human | S1 | 0.41 | 1.44 | 0.24 |
| human | S2 | 0.57 | 1.02 | 0.33 |
| human | S3 | 0.57 | 1.27 | 0.33 |

- Gene-panel noise only, so this is a lower bound on total positional error. Segment lengths are
  taken within notebook 37's common support (0.018–0.864).
- Only mouse S2 and S3 have noise clearly below the spread of positions.
- What noise does to each feature:
  - It cannot create a consistent ordering of genes; to first order it shifts all of them together.
  - It inflates every transition width.
  - It erases features narrower than about 0.2–0.3 segment units.
- This orders the categories by robustness: ordering > gradient sign > switch location relative to a
  boundary > interior peaks > sharpness.

### 3.2 Categories

**(i) Within-segment gradients**

- **Test:**
  - Per species and specimen: a within-segment linear term (degree-1 Legendre per segment) beyond
    segment steps, fitted under both test A and test B.
  - For species: the beyond-steps term C5, plus a test-B version.
- **Discrete comparator:** segment means. Gradients are invisible to them (3.0 (1)).
- **Existing numbers:**
  - Beyond-step pathway calls 79 / 41 / 81 on the three coordinates. Gene scores agree across
    coordinates at ρ 0.84 and 0.79.
  - About 35–39% of the positional species-difference SS in the top 10% of genes lies beyond steps.
    In relabelings the share is 36–81%, so specimen noise is disproportionately beyond-step.
  - Within-segment slopes agree with external adjacent-segment contrasts at only 0.17 (mouse) and
    0.10 (human) (P1w, fold-protected).
  - Caveat on P1w: it is low both when the coordinate is poor and when the biology is step-like, so it
    is not a clean coordinate metric. Gene-fold within-segment agreement (0.49; human S1 0.03–0.46) is
    the cleaner one.
- **Validation:**
  - Between-specimen correlation of slopes on fold-protected genes;
  - the count-split coordinate;
  - DPT13, SCF13-ED and later workstream A;
  - an external within-segment pseudotime in Census adult male mouse snRNA (59,518 PT nuclei with
    S1/S2/S3 labels are on disk).
- **Status:** unproven. With existing tables, label noise cannot be excluded.

**(ii) Non-monotone or transient programs that peak inside a segment**

- **Test:** the extremum lies in the inner 60% of a segment and exceeds both segment-end values by at
  least 0.2 log units, in both specimens and on at least 2 of 3 coordinates.
- **Discrete comparator:** attenuated by w/L. "S2-high" programs are visible to segment means and do
  not count.
- **Existing numbers** (saved species curves, amplitude ≥ 0.5, descriptive):

  | | Mouse | Human |
  |---|---|---|
  | Interior-extremum genes, scFates | 360 (116 S1, 28 S2, 216 S3) | 132 |
  | Interior-extremum genes, DPT13 | 51 (47 in S2) | 32 |
  | Shared between the two coordinates | 44 | all 32, all in S1 |

  - Human S1 order is depth-dependent: position correlates with library size at −0.31 and −0.38. The
    equal-depth refit agrees with the original order in human S1 at only 0.36–0.46.
  - Pathway level: 4 of 36 robust pathways are "localized", and their shape classes agree between
    scFates and DPT13 for 17 of 22.
- **Status:** the weakest category. It is coordinate-dependent and attenuated by noise. Do not lead
  with it.

**(iii) Location and sharpness of transitions; species shifts relative to segment boundaries**

- **Test:** per gene and specimen, a switch fit (logistic: location and width) in that species' own
  segment units. The species shift is onset_h − onset_m. Claim it only when it exceeds the
  registration noise and has the same sign in all 4 cross-species specimen pairs.
- **Discrete comparator:**
  - The segment switch fraction r confounds location with width, so a shift is blurred.
  - "Opposite signs in adjacent segments" cannot separate a reversal from a boundary shift. Of the 70
    reversal genes, 54 are visible to segment DESeq2 with opposite signs. Only 21 survive a per-gene
    shift; the method-selection doc describes the rest as species shifts of the expression boundary.
    This changes interpretation, not detection.
- **Existing numbers:**
  - 41 genes are monotone in both species on scFates, with a median |onset shift| of 0.42 segment
    units. On DPT13 there are 11 such genes, with 0.37.
  - Registration noise: human segment transitions move by up to 0.29 coordinate units across refits,
    and human S1 < S2 < S3 broke in one gene-fold refit.
  - Sharpness is bounded below by the noise (3.1). No "gradual vs sharp" claim is allowed below about
    0.3 segment units.
- **Status:** plausible for a handful of genes; the bar is high.

**(iv) Order of program onsets along the axis (lead candidate)**

- **Test:** half-amplitude onsets per gene and specimen for monotone genes. The claim is the partial
  correlation of onsets between specimens given r. The same test applies at pathway level with member
  medians.
- **Discrete comparator:** r gives a coarse order, which is blurred, not absent.
- **Existing numbers:**
  - 202 mouse genes switch on inside S2 under both coordinates, with onset ρ 0.90 between scFates and
    DPT13 and a median |difference| of 0.13 segment units. Their onsets spread over 1.16–1.92 (IQR).
  - Away from both boundaries (91 genes), ρ is 0.69.
  - r predicts the continuous onset at ρ 0.69, and r replicates between the mice at 0.95.
  - Pathway-level median onsets agree only at ρ 0.53 across coordinates (13 pathways). Gene-level
    order is the robust unit.
- **Precedents:** liver bile-acid enzymes ordered along the lobule as in their cascade (Halpern 2017);
  enterocyte absorption programs in sequence along the villus (Moor 2018).
- **Status:** the best candidate. It needs specimen-level and external replication beyond r.

**(v) Power and specificity at a fixed specimen count**

- **Existing numbers:** T_spatial 60 (2, 0) against coordinate-free steps 50 (0, 0) and 6 bins 49
  (17, 0).
- Only near-threshold pathways are coordinate-only, and the two primary ones are probe-sensitive.
- The dichotomisation literature predicts a power loss only for the within-segment share, about a third
  of the positional SS here (Royston 2006).
- **Status:** not a selling point. Report it as "the continuous screen is as specific as the step screen
  and slightly more powerful".

## 4 · How much of the 35/36 beyond-step signal replicates (existing tables only)

| Check | Result |
|---|---|
| Robust 36 called beyond steps: scFates / DPT13 / SCF13-ED / all three | 35 / 18 / 34 / 18 |
| Primary 22: scFates / DPT13 / SCF13-ED / all three | 22 / 12 / 21 / 12 |
| Core 26: scFates / DPT13 / all three | 25 / 14 / 14 |
| Beyond-step pathway calls shared across all three coordinates (of 79) | 33; 10 of them are not coordinate-free step calls |
| Beyond-step z_joint, robust 36: median species / median max-relabel | 3.56 / 2.28 |
| Gene-level beyond-step score, Spearman across coordinates; overlap of top 10% | 0.84 (DPT13), 0.79 (SCF13-ED); 0.88, 0.86 |
| Genes in the top 10% beyond steps but below the median for steps (scFates) | 56; 79% still top 10% beyond on DPT13, 61% on SCF13-ED |
| Specimens: LOSO retention of the 36 (T_spatial joint; no beyond-step LOSO saved) | 35, 36, 36, 32 of 36 (omitting Ctrl1A2, Ctrl1A4, COR1, MED1) |
| Within-species positional calls (joint, mouse vs mouse, human vs human) | 0 and 0 on every coordinate |
| Count split (old 66 list, T_spatial) | 59 of 66 robust kept with an independent coordinate |

What this means:
- About half the robust list (18 of 36; 12 of 22) carries a beyond-step signal that survives a change
  of trajectory method on the same embedding.
- None of these checks is independent of the shared embedding or of the segment labels. Specimen-level
  and count-split replication of the beyond-step term specifically does not exist yet.

## 5 · Weaknesses that cut against the framing

1. Within-segment order is weak:
   - P1w is 0.17 (mouse) and 0.10 (human);
   - gene-fold within-segment agreement is 0.49, and 0.03–0.46 in human S1;
   - fold noise is at least the size of the within-segment spread outside mouse S2/S3;
   - the coordinate predicts held-out genes no better than segment means (P4 gain ≈ 0).
2. Human S1 order is depth-dependent. Position correlates with library size at −0.31 and −0.38, and the
   equal-depth refit agrees with the original order at only 0.36–0.46.
3. Peaks move between coordinates: ρ 0.71, and 7 of 22 change segment. Three pathways move from about
   0.6–0.75 to near 0 under SCF13-ED.
4. Beyond-step calls halve on DPT13 (35 → 18 of the robust 36).
5. The coordinate-free step screen recovers 20 of 22, and six coordinate bins recover 22 of 22.
6. Notebook 10 (03's DPT) found 2 strict trajectory-only genes. Notebook 07's modules were unstable
   under LOSO.
7. Lake 2023 kept PT-S1 and PT-S2 merged because they "could not be accurately resolved" in integrated
   data. Human S1/S2 is weakly separable even discretely. I found no source describing healthy human PT
   as a continuum: Lake's PT trajectory is a repair trajectory. The brief's "Lake 2023 continuous PT
   states" is not supported by what I read.
8. ORA passing relabeling (2.2) weakens the clean "average-level fails, positional passes" contrast in
   R4.

## 6 · Framing options

| Option | Claim | Result it needs (pre-registered in section 8) | Current status |
|---|---|---|---|
| F1 "More pathways / more power" | Continuous analysis finds programs that segment analysis misses | Coordinate-only calls that replicate externally beyond the step screen | Fails now: 2 of 22 not called by steps, both probe-sensitive; 8 of 36, all near misses. **Not recommended** |
| F2 "More information per program, at a measured resolution" (recommended) | For each program the coordinate adds where it switches, in what order, and how it shifts between species. Claims are limited to scales above the coordinate's positional noise | (a) Within-S2 mouse onset order replicates between specimens beyond r (partial ρ ≥ 0.3), and in external mouse snRNA beyond r (partial ρ ≥ 0.2). (b) At least 10 genes, or at least 3 programs, show a species boundary shift above registration noise in all 4 specimen pairs. (c) Resolution budget reported per coordinate | (a) Cross-coordinate ρ 0.90; r explains ρ 0.69. Untested at specimen and external level. (b) Untested. (c) Done for SCF13 (3.1) |
| F3 "Continuous ≈ fine discretisation, and that is the finding" (fallback) | At Visium-HD structure resolution with these coordinates, a continuous analysis equals a 3–6-bin analysis in calls. Its value is within-specimen replication and location read-out at segment resolution | F2(a) fails; the decomposition shows that the beyond-step share does not replicate | Consistent with every existing number |

Workstream A fits into F2 directly. A better reconstruction lowers the positional noise and should
raise the replicated continuous-only information. The paper can then show that information gain
tracks reconstruction quality across coordinates (DPT13, SCF13, SCF13-ED and A's method). That link
is the strongest version of the re-centred paper. If the information gain does not track the
reconstruction metric, the method's value has to be argued on other grounds.

## 7 · Precedents (verified in Europe PMC; WebSearch was capped)

- **Liver.** Reconstructing lobule coordinates showed about 50% of genes zonated, with abundant
  non-monotonic profiles peaking mid-lobule. Bile-acid enzymes are ordered as in their cascade:
  Halpern et al. 2017 ([PMID 28166538](https://pubmed.ncbi.nlm.nih.gov/28166538/)); review by
  Ben-Moshe & Itzkovitz 2019 ([PMID 30936469](https://pubmed.ncbi.nlm.nih.gov/30936469/)).
  Mixed-effect space × time models: Droin et al. 2021 ([PMID 33432202](https://pubmed.ncbi.nlm.nih.gov/33432202/)).
  All are precedents for categories (ii) and (iv).
- **Intestine.** Landmark-based villus reconstruction found sequential carbohydrate, peptide and fat
  absorption programs in distinct villus compartments: Moor et al. 2018
  ([PMID 30270040](https://pubmed.ncbi.nlm.nih.gov/30270040/)). ClumpSeq extended this to rare cells:
  Manco et al. 2021 ([PMID 34031373](https://pubmed.ncbi.nlm.nih.gov/34031373/)).
- **Kidney, discrete and continuous.**
  - Ransick et al. 2019, anatomy-guided scRNA-seq with sex and regional diversity
    ([PMID 31689386](https://pubmed.ncbi.nlm.nih.gov/31689386/); abstract only).
  - Chen, Chou & Knepper 2021, all 14 microdissected tubule segments, 3,709 segment-specific
    transcripts ([PMID 33769951](https://pubmed.ncbi.nlm.nih.gov/33769951/)).
  - Lake et al. 2023, PT-S1/S2 merged and a PT repair trajectory
    ([PMID 37468583](https://pubmed.ncbi.nlm.nih.gov/37468583/), OA full text).
  - Hinze et al. 2021, cells ordered along the corticomedullary axis; "spatial gene expression in the
    kidney changes gradually" ([PMID 33239393](https://pubmed.ncbi.nlm.nih.gov/33239393/)).
- **Continuous variation within discrete types:** Cembrowski & Menon 2018
  ([PMID 29576429](https://pubmed.ncbi.nlm.nih.gov/29576429/)); Trapnell 2015
  ([PMID 26430159](https://pubmed.ncbi.nlm.nih.gov/26430159/)).
- **Trajectory DE.**
  - tradeSeq (Van den Berge et al. 2020; [PMID 32139671](https://pubmed.ncbi.nlm.nih.gov/32139671/)).
    It argues that cluster-based procedures "fail to exploit the continuous resolution".
  - condiments (Roux de Bézieux et al. 2024; [PMID 38280860](https://pubmed.ncbi.nlm.nih.gov/38280860/)).
  - Lamian (Hou et al. 2023; [PMID 37949861](https://pubmed.ncbi.nlm.nih.gov/37949861/)), whose
    sample-variability argument matches ours.
  - PseudotimeDE (Song & Li 2021; [PMID 33926517](https://pubmed.ncbi.nlm.nih.gov/33926517/)) on
    pseudotime uncertainty.
  - switchde (Campbell & Yau 2017; [PMID 28011787](https://pubmed.ncbi.nlm.nih.gov/28011787/)) on switch
    location and rate.
  - GeneSwitches (Cao et al. 2020; [PMID 32058565](https://pubmed.ncbi.nlm.nih.gov/32058565/)) on
    switch ordering.
  - Genes2Genes (Sumanaweera et al. 2025; [PMID 39300283](https://pubmed.ncbi.nlm.nih.gov/39300283/)).
- **Activity scores:** AUCell in SCENIC (Aibar et al. 2017;
  [PMID 28991892](https://pubmed.ncbi.nlm.nih.gov/28991892/)); decoupleR (Badia-i-Mompel et al. 2022;
  [PMID 36699385](https://pubmed.ncbi.nlm.nih.gov/36699385/)).
- **Dichotomisation:** Royston, Altman & Sauerbrei 2006, on loss of power and on cutpoint bias
  ([PMID 16217841](https://pubmed.ncbi.nlm.nih.gov/16217841/)); Altman & Royston 2006 BMJ
  ([PMID 16675816](https://pubmed.ncbi.nlm.nih.gov/16675816/); title only seen).
- **Post-selection inference:** count splitting (Neufeld et al. 2023;
  [PMID 36511385](https://pubmed.ncbi.nlm.nih.gov/36511385/)).

How they bear on the framing:
- The liver and intestine precedents rest on coordinates validated by smFISH or LCM landmarks. Ours is
  validated only weakly within segments, so an ordering claim must carry its own replication.
- In kidney, the precedent for continuity is corticomedullary (Hinze), not within-PT.

## 8 · Phase-2 protocol (pre-registered; notebooks 56–60, `results/pt_pathway_continuous/`)

These rules hold throughout:
- Coordinates are SCF13 (primary), DPT13 and SCF13-ED, and later workstream A's, all run unchanged.
- Every thresholded statistic is frozen here.
- Nothing is tuned on the species contrast.
- A result that fails its rule is reported as failed.

**Notebook 56 · `56_pt_pathway_taxonomy_completion`**
- **(a) ORA control.**
  - Rerun D3/D3b on genes with equal probe counts.
  - Run ORA in notebook 45's same-species pools: 2 + 2 donor splits with design `~group` (2 residual
    df), plus the 100 cross-species draws, using the same moderated-t pseudobulk statistics.
  - ORA is called **specific** if all three hold:
    - the 90th percentile of null calls is ≤ 5% of tested pathways in all four pools;
    - the median NCDR across cross-species draws is < 0.25;
    - at least 80% of draws meet the rule.
  - Otherwise R4 is reworded: "passes our relabeling through power loss".
  - If ≥ 50% of the 118 calls are lost on equal-probe genes, ORA is labelled probe-exposed.
- **(b) Specimen-level pathway-score split plot.**
  - Mean-z and AUCell-like rank scores per structure; specimen × 8-bin means.
  - Moderated F for species × bin against specimen × bin (14 df), under the species split and both
    relabelings.
  - Kept as a self-contained comparator if the rule passes.
- **(c) Curve modules.**
  - Top 10% of T_spatial genes; Δ(s) curves standardised.
  - k chosen by the stability criterion on the species run, then frozen. Hypergeometric ORA per
    module.
  - The identical pipeline is run on each relabeling's top 10%.
  - Kept if NCDR < 0.25, each relabeling calls < 5%, and median LOSO ARI ≥ 0.6.
- **(d) Onset and peak enrichment.**
  - Matched random sets, as in notebook 37 strata. Statistic: the dispersion of member onsets, and
    separately the member-median onset.
  - Same specificity rule.
- **(e) tradeSeq equivalents** (association per species, startVsEnd, earlyDE in the first 0.33 of S1)
  inside the existing nested design. Report gene-rank Spearman with T_spatial; no new pathway calls.

**Notebook 57 · `57_pt_beyond_steps_identifiability`**
- Two semi-synthetic arms on the real structures, coordinates and labels:
  - **(A) truth = coordinate-threshold steps** (per-specimen thresholds maximising agreement with
    labels, frozen);
  - **(B) truth = label steps.**
- In each arm, residuals are resampled within specimen × segment and genes keep their fitted step
  amplitudes. Ten replicates are run, each with the full joint pathway test.
- **Primary metric:** beyond-label-step pathway calls in arm A, divided by the observed 79.
- **Decision:**
  - If the arm-A median is ≥ 40 (≥ 50% of 79), "beyond steps" is reported as label–coordinate
    discordance, and no within-segment pathway claim is made from it.
  - If < 16 (20%), the label-noise explanation is rejected.
- **Also:** the test-B (beyond coordinate steps) statistic on real data, with arm B as its null. A
  pathway is a within-segment candidate only if it passes both A and B on SCF13 and on one other
  coordinate.

**Notebook 58 · `58_pt_continuous_descriptors_replication`** (categories i–iv)
- **Data and fits:**
  - Per specimen: smooth curves (6-df spline) within the common support.
  - Monotone genes: amplitude ≥ 0.5 log units and |ρ(position)| ≥ 0.9 in that specimen.
  - Per gene and specimen: onset (half-amplitude), width (10–90%), within-segment slopes and interior
    extrema, all in that specimen's own segment units. Also r from that specimen's segment means.
  - Gene sets: all genes, and fold-protected genes (not among the 611 coordinate-construction genes).
  - Coordinate arms: SCF13, the count-split coordinate (nb35 `refit_split_A`, tested on the B counts),
    the gene-fold refits, DPT13 and SCF13-ED.
- **Primary (iv):**
  - Within mouse S2: the partial Spearman of onset(Ctrl1A2) with onset(Ctrl1A4), given each
    specimen's r, over genes monotone in both. Gene-bootstrap 95% CI, on fold-protected genes, on SCF13.
  - **Claim** if the partial ρ is ≥ 0.3 with the CI excluding 0 and is ≥ 0.2 on the count-split arm.
  - **Abandon** the onset claim if the partial ρ is < 0.1 on both arms.
  - Human (COR1 against MED1, within S2 and S3) is reported, but cannot support a donor-level claim.
- **(i):**
  - Per species × segment cell: the between-specimen Spearman of within-segment slopes over
    fold-protected genes with |slope| SE-standardised.
  - Claimed for a cell if ρ ≥ 0.3 on SCF13 and ≥ 0.2 under count split. Otherwise the cell is "not
    resolved".
- **(ii):**
  - An interior extremum is claimed for a gene only if it appears in the same segment in both
    specimens and on ≥ 2 of 3 coordinates, and survives the equal-depth arm (human S1 especially).
  - Expected: few. If more than 20, inspect for depth.
- **(iii):**
  - Species boundary shift = onset_h − onset_m in own-species segment units, with all 4 cross-species
    pairs having the same sign.
  - |shift| must exceed the 95th percentile of a registration null: the same statistic with the human
    transitions moved to each refit's registration (nb35 `species_registration_across_refits`).
  - Claimed for the paper if ≥ 10 genes or ≥ 3 programs pass.

**Notebook 59 · `59_pt_external_onset_replication`**
- **External coordinate.**
  - Census mouse snRNA PT nuclei (`data/external/census_pt_cells`), adult males only (3–6 and
    12 months; 21-day and Theiler stages excluded).
  - A PT pseudotime is built per donor from a gene fold disjoint from the test genes: diffusion
    pseudotime on PCA, rooted by the same early-marker rule.
  - Onsets and r are computed per donor and averaged over donors.
- **Primary metric:** partial Spearman of our mouse onsets (fold-protected, SCF13) with the external
  onsets, given our r, within S2.
- **Decision:**
  - **Claim** if partial ρ ≥ 0.2 with a donor-bootstrap CI excluding 0.
  - **Abandon** F2(a) if partial ρ < 0.1 here and in notebook 58.
- **Human:** Lake/Census cortex PT, with S1/S2 merged, so only early → S3 onsets. Descriptive only.

**Notebook 60 · `60_pt_resolution_budget_summary`**
- One table per coordinate:
  - positional noise per species × segment (gene folds and count split);
  - the replicated information metrics from notebooks 58/59;
  - beyond-step identifiability from notebook 57;
  - the call overlap with D5.
- Run on workstream A's coordinate when it exists.
- **Paper rule.** A continuous-only statement goes into the paper only if it passes on the paper's
  coordinate and is not reversed on DPT13. The relation between noise and replicated information
  across coordinates is reported as a figure, with no test (k ≤ 4 coordinates).

**Order and cost.**
- Run 58 (iv) and 57 first: they decide F2 against F3.
- 56(a) is cheap and decides the R4 wording.
- 59 depends on 58.
- All of these reuse cached coordinates and the existing nested-design code.
- New reusable code goes into `pseudospace/continuous_descriptors.py` with tests on synthetic sigmoids,
  bumps and label-noise fixtures.

## 9 · Changes suggested to existing files (for the coordinator)

- `docs/paper/draft.md` R4 conclusion and Supplementary Note 2: ORA on segment DESeq2 passes our
  relabeling (118; 0, 1). Narrow the sentence to rank-based competitive screens, or add the notebook-56
  external control before claiming generality.
- `docs/paper/outline.md`: replace "more pathways" with the resolution-budget framing (F2/F3); note that
  the beyond-step share of the robust list is 18 of 36 when it must hold on all three coordinates.

## Files created in this phase

- `workstream-B scratch: plan_phase1.md` (this file)
- `workstream-B scratch: q1_beyond_steps.py`, `q1_beyond_steps.out`
- `workstream-B scratch: q2_segment_ora.py`, `q2_segment_ora.out`
- `workstream-B scratch: q3_gene_scores.out` (inline analysis of cached gene scores; the commands are recorded in
  this session)
- `workstream-B scratch: q4_curve_sizing.py`, `q4_curve_sizing.out`
- `workstream-B scratch: q5_overlaps.out`, `q6_resolution.out`
