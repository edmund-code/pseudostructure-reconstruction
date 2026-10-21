# The human-vs-mouse axis: what is paper-ready after the package switch

Written 2026-10-20. The species-axis companion to `paper-ready-results.md`, which ranks both axes together.

> **Corrected 2026-10-20, after an adversarial audit of every claim on this page.** The user doubted the claims and
> commissioned a parallel audit: 11 claim groups checked against the committed tables, plus 4 skeptics attacking the
> interpretations. Of 71 numeric claims, 56 reproduced, 9 were misquoted and 6 were contradicted. **Two of the three
> structural claims this page opened with do not survive**, and the audit found a provenance problem that is worse
> than the one this page named as blocking. The headline sections below are rewritten; the record of what was wrong
> is in §7 so the mistakes are not silently absorbed.

Sources: `docs/results/pt-package-issues.md`, `pt-package-species-calls.md` (77),
`pt-package-continuous-vs-discrete.md` (74), `pt-coexpression-conservation.md` (97), `pt-non-one-to-one.md` (86),
`pt-positional-decomposition.md` part II (91), `pt-gene-atlas.md` and `pt-evidence-model.md` (112, 113),
`pt-published-findings.md` (81, 82), `pt-disease-genetics.md` (108), `pt-transporter-handoffs.md` (101),
`pt-package-label-comparison.md` (73).

## 0 · The structural conclusion, corrected

The earlier version of this section claimed "the species result is gene-level, not pathway-level" on three legs.
Two of the three legs do not support it.

- **Leg 1 was mis-stated.** Notebook 77's calls do lose specificity against a selection-matched null (p 0.35, 0.64),
  but that null matches random sets on the *z-quintile of the statistic that defined the calls*, and it drives every
  statistic in the table to zero — including the comparison against noise calls that this page cited as a positive.
  The gene arm was never given a comparably strict null, so "genes, not sets" follows from the choice of nulls
  rather than from a like-for-like comparison. Worse, the page omitted notebook 77's own per-call table:
  **13 of 83 sets replicate at BH ≤ 0.10** against Census male and 13 against Lake male, 11 in at least two of four
  references, and three in all four — Reactome Biological Oxidations, Reactome Phase II Conjugation, Hallmark
  Hypoxia. Set-level replication exists; this page dropped it.
- **Leg 2 is honest but is a power statement.** Of 381 sets none is called at q ≤ 0.05 (verified; min q 0.789).
  The raw p distribution really is close to null (16 of 381 at p ≤ 0.05 against 19.1 expected, KS p 0.258), so
  nothing is hiding behind the BH burden. But 159 of 381 sets meet every *structural* conserved criterion except q,
  and the top of the raw ranking is the same biology that replicates in 77 — SLC-mediated Transmembrane Transport
  (z 2.47), Biological Oxidations (the global minimum p). "No set survives BH over 381 sets with two specimens per
  species" is not "no set-level species program exists".
- **Leg 3 is irrelevant to the claim.** Notebook 113's species model contains no set-level predictor at all; all
  seven species lines are per-gene calls, and the model's *outcome* is membership of a curated gene set. It cannot
  rank a gene-level method above a set-level one.

**The defensible version, which is what the paper should say:** *set calls add nothing beyond their member genes'
positional strength, and the 381-set screen has no per-set power in this cohort.* That is narrower, it is true, and
it still points the paper at genes — but it is not evidence that pathway-level species biology is absent.

## 1 · The strongest species results, ranked

### S1 · The reversed gradients — one real line of evidence, not three
The earlier version claimed three independent corroborations. **It was circular, and I verified the circularity
directly rather than taking the auditor's word for it.**

Notebook 112's species atlas has **seven** lines (not nine — nine is the AKI count), and two of them are
`S5 rewired` (notebook 97's own statistic) and `S6 reversed gradient` (the reversal call itself). So the convergence
core counts the claim as evidence for the claim. Measured on `results/pt_gene_atlas/tables/species_atlas.csv.gz`:

| core definition | genes | reversed | share | base rate |
|---|---|---|---|---|
| k ≥ 5 of 7 (as published) | 24 | 15 | 0.625 | 0.132 |
| drop S6, k ≥ 5 of 6 | 14 | 5 | 0.357 | 0.132 |
| drop S5 and S6, k ≥ 5 of 5 | 7 | 1 | 0.143 | 0.132 |
| drop S3, S5, S6, k ≥ 4 of 4 | 19 | 1 | **0.053** | 0.132 |

S3 has to go too: `P(S3 zonation-class call | zonation class == 'reversal')` is **1.000** over 65 genes, so S3 is a
structural superset of the reversal call, not an independent vote. With the three reversal-informative lines
removed the core is *below* the base rate for reversal — and what it actually is is a **transporter list**:
Slc23a1, Slc23a3, Slc5a12, Slc27a2, Slc9a8, Abcd3, Slc46a3, Slc6a6, Slc6a13, Slc6a18, Slc7a8, Slc7a13, Slc16a9,
Slc16a12, Slc39a5, Slc22a8, Slc25a17, Abcc4, Sfxn4.

The other two "independent" strands also fail:
- Notebook 97's rewiring D is quoted from the **total** layer, where a gene's top co-expression partners are largely
  the genes with the same positional profile — so a reversed gradient mechanically acquires different partners. On
  the **residual** layer, with position removed, the separation falls from 0.761 / 0.213 to **0.386 / 0.295**. And
  notebook 97 failed its own reproducibility gate (median W_m 0.670 against a 0.7 threshold; verdict "not shown
  useful"), which this page did not mention.
- The "replicates in independent atlases (ρ 0.30)" number is a genome-wide correlation validating the D measure,
  not a per-gene replication of the reversal list.
- Notebook 106's sequence negative is computed on the gene table the reversal list is *read from*, so it is not an
  outside check; and at n ≈ 70 reversed genes it can only exclude an AUC above ~0.55.
- The list overlap was also miscounted: **11** of 97's 25 rewired-in-both genes are on the reversal list (Mme,
  Sh3gl2, Lap3, Mgll, Hpgd, Mfap3l, Slc16a9, Slc41a2, Slc17a3, Acot12, Ugt3a1), not 7. An 11-of-25 overlap against
  a 0.8% base rate is the best available *proof that S5 and S6 are one line*, not two.

**What survives:** one level-free method — co-expression neighbourhood rewiring — separates reversed from conserved
genes in our own specimens (0.761 / 0.213 on the total layer, 0.386 / 0.295 with position removed). There is also a
per-gene external corroboration the auditor found while trying to refute the claim, which the original page did not
cite and which is stronger than what it did cite; it is recorded in `pt-coexpression-conservation.md` and is set
aside here because the user has scoped the current work to our own data.

**And the repo's own calibration ranks this last.** In notebook 113 the two species lines carrying this headline are
the only two whose coefficients are indistinguishable from zero: rewired 0.188 [−0.486, 0.887] and reversed
gradient **−0.196** [−0.760, 0.343]. The reversal list should not lead the paper's biology.

### S2 · Conservation by gradient strength — the re-run already exists
The earlier version called this "the single blocking item". **That was wrong.** Notebook 73 re-runs notebook 49's
part A on the package units and labels (`part49(...)` on `label_units(...)`, which reads notebook 70's export) and
its outputs are committed at `results/pt_package_label_comparison/tables/49/`. `pt-package-label-comparison.md`
records it as rule M5, "holds (gene and module); same verdict".

| | retired coordinate (49) | package (73) |
|---|---|---|
| index, top 1% (S2−S1) | 0.569 | 0.614 |
| index, top 1% (S3−early) | 0.576 | 0.564 |
| index, all genes | 0.165 / 0.274 / 0.287 | 0.189 / 0.267 / 0.288 |
| share conserved, strongest bin | 0.889 | up to 1.000 (n 53) |

Two further misquotes to fix in the paper text: the top-1% point estimates are **0.569 and 0.576**, not "≈ 0.55"
(0.54 was a confidence bound quoted as a point estimate); and the 0–50% bin index is **not estimable** (NaN in every
bootstrap draw), not "near zero". "Direction agreement rising 0.5 → 0.85–0.89" conflates two columns: raw sign
agreement over all genes runs 0.482 → 0.667, while *share conserved* among genes zonated in both species runs
0.571 → 0.889.

### S3 · The conserved core is transport, and this is now the best-supported species claim
With the reversal story demoted, transport is what is left standing — and it is standing on more than before,
because the de-circularised convergence core in S1 is itself a transporter list.

- 16 genes conserved in both species' co-expression and in gradient (97): Slc4a4, Slc6a19, Slc7a8, Dpep1, Gpx3,
  Apoe, Iyd, Lgmn, Tpp1, Cotl1, Sptan1, Itpr2, Fam151a, Tbx2, Atp6v0c, Gnptg.
- **33 of 43 zonated transporters conserved (0.767) against 196 of 388 zonated genes (0.505).** This page and three
  others called it a hand count; it is not. Notebook 47 writes it to
  `results/pt_revision_addenda/tables/transporter_class_shares.csv`, and it is reproducible from notebook 43's
  `gene_classes_symmetric.csv`. Open item 18 is closed; `docs/paper/open_items.md`, `draft_v03_changelog.md` and
  `referee_response_v01.md` all still call it a hand count and are stale.
- **11 conserved transporter handoffs, 0 reversed, against 0.549 expected by chance** — amino acids, organic anions,
  SGLT2→SGLT5. Two corrections: external confirmation is 11/11 in mouse microdissection but **9 of 11** in the human
  atlas (both Slc22a13 pairs fail); and "inside S2" holds only loosely — by the notebook's own margin rule just
  3 of 11 human crossings sit inside a segment, with 7 within 0.05 of the S2/S3 cut.
- Mendelian proximal-tubulopathy genes sit toward S1 in both species (mean E −1.47 mouse, −0.56 human, beyond
  matched nulls); 25 of 27 literature positions reproduced; OCT2 rises toward S3 in mouse and is flat in human (108).

### S4 · Non-1:1 orthology — the section that survived the audit intact
Every claim in this section reproduced from the tables (86):
- Strong-zonation shares by orthology class: 1:1 3.8%, 1:many 9.9%, many:many 14.2%, many:1 **1.7%** — the excess is
  in mouse genes of mouse-expanded families, not in mouse single-copy genes whose human side expanded. MH OR
  1.84–4.05 under every sensitivity. 41 strong probe-clean non-1:1 mouse genes, 25 outside the shared space.
- Sex-biased mouse PT genes are enriched in non-1:1 families: any bias OR 1.50 (p 0.001), **male bias 1.90**
  (p 10⁻⁴), female 1.16 (p 0.40) — an open question in the gene-age literature, never tested in kidney, and the sex
  bias comes from data independent of our all-male cohort.
- Tandem paralogs under 1 Mb share position: median r 0.31 against 0.07 dispersed and 0.04 non-paralog neighbours
  (human 0.48 / 0.07 / 0.04); tandem status predicts shared position and sequence identity does not (t 3.7 vs 0.5).
- Counterweights all verified: human counterparts zonated 20% against 22%; family sums fail calibration at 0.089;
  paralogs agree with microdissection 65% against 70%.

This is the most robust biological section on the species axis, and it should be promoted in the Results order.

### S5 · Where the species differ
36 human-high and 47 mouse-high calls; peak segments human-high S1 14 / S2 4 / S3 18 and mouse-high 17 / 4 / 26
(verified). One wording fix: human-high S3 is 18 of 36, exactly 50% — a plurality, not a majority — so "most calls
differ most in S3 in both directions" holds for mouse-high (0.55) and should be softened for human-high. 89% of
resolved calls summarise mouse zonation the human sections do not share (73 of 82, verified), but note that
`category_counts.csv` carries an alternate mt-in-library split of 68/6/1/7/1 which would change that fraction.

The 83-call list is also weaker than the headline count: 0 of 83 survive either relabeling, only 60 are "core",
only 32 pass the NB joint segment-step test, and 8 are probe-sensitive. Gene overlap between calls is substantial —
single-linkage at Jaccard > 0.25 gives 53 components, not the 77 redundancy clusters — so the independent-statement
count is 43–53.

## 2 · What the switch broke

1. **"Continuous beats discrete" is weak at set level.** Segment F-test 76 calls, 71 shared with the positional
   F-test; the spline adds 12 clusters (10 near misses), the steps add 5. The supported message is narrower:
   testing *shape against specimen shape noise* finds species differences that level screens and conventional DE
   miss or confuse with probe-panel differences.
2. **Within-segment claims are closed off by the frozen noise gate:** human gene-fold 0.49, equal depth 0.73 (human
   S1 0.55), species registration gap 0.37, 10 of 12 cells noise-dominated, one fold refit reversing the mouse PT.
3. **The landmark comparators:** 0 of 12 noise-dominated cells for both LCP and LRS-Z, but the gene-fold band
   0.93–0.97 describes **LRS-Z only** — LCP is 0.962 mouse / 0.820 human and 0.647 in human S1. Corrected.
4. Species F-test calibration passes narrowly at 0.071 against the 0.075 limit.
5. Human PT-S2 shrinks to 632 units from 1,813; S2−S1 classified genes 9,013 → 6,269; human units are effectively
   undefined for within-species statements.
6. Human zonation is the weakest measurement: human S3−S1 agrees with microdissected human S1/S3 at ρ 0.14 against
   0.46 for mouse.

## 3 · Three confounds specific to this axis

1. **Probe panels differ, so levels are never readable.** Probe exposure: F-test 0.17, clustering 0.69. Structural
   zeros include mouse Gsta1 and human UGT2B15, AKR1C1.
2. **Neighbour signal is asymmetric between species**, which is worse than its being present. Verified per
   specimen: mouse centred S3−S1 is S100g −5.00/−4.70, Slc12a3 −4.21/−3.40, Calb1 −4.49/−4.56, Pvalb −4.98/−3.31,
   Slc8a1 −4.25/−3.85, Nphs2 −3.70/−4.33, Ren1 −2.34/−2.01. The human bound is **−1.4 per specimen** (Ren1 −1.36 in
   HUK1_MED1), not −1.2; −1.2 is only reachable by averaging the two human sections, which would compare a human
   mean against mouse per-specimen values. The qualitative contrast (mouse 3–5× stronger) is unaffected. Species
   calls are 1.4–1.9× enriched for edge-enriched genes, the F-test least at 1.44×.
3. **The human tissue is one donor, two cortex sections**, with no counterpart to mouse outer-stripe S3.

## 4 · Recommended Results structure, revised

| # | Section | Rests on | Changed? |
|---|---|---|---|
| R1 | Cohort, structures, the reconstructed axis and its measured limits | 70, 72 | — |
| R2 | A gene-level species comparison: what shape-against-noise testing finds that level screens and conventional DE do not | 71, 74, 75 | — |
| R3 | **The conserved transport core**: the de-circularised convergence core, the 11 handoffs, 33 of 43 transporters, disease genes toward S1 | 112 (de-circularised), 97, 101, 108, 47 | **promoted** from R4 |
| R4 | **Non-1:1 orthology and the sex program inside it** | 86 | **promoted** from R5 |
| R5 | Conservation rises with gradient strength | 49 **and 73's package re-run** | no longer blocked |
| R6 | Reversed gradients, on one level-free line of evidence, with the circularity stated | 90 D, 91 S, 97 residual layer | **demoted** from R3 |
| R7 | Negatives: no per-set power in the 381-set screen; no sequence prior; co-expression below its own reproducibility bar; paralog substitution fails all three tests | 91 S2, 106, 97, 88 | — |

## 5 · What to package, revised

1. **The positional F-test with its relabeling calibration**, shipped with the probe-exposure and neighbour-marker
   filters as a required pre-step. *The justification has changed:* the old version cited notebook 113's 0.88
   coefficient, which does not support it (see §7). The justification is notebook 71's calibration (0.071 against
   0.075) and notebook 100/103's planted-truth profile.
2. **Orthology handling** — the 1:1 space, family sums with their calibration failure wired in as a refusal, and
   the non-1:1 class comparisons. Promoted to second because S4 is the section that survived the audit intact.
3. **The level-free comparison layer**, shipped with its scope stated: notebook 97 failed its own genome-wide gate
   (median W_m 0.670 < 0.7), the headline total layer is position-contaminated, and the residual layer is the one to
   use. It is a corroboration layer on an existing list, not a discovery screen.

**Do not ship** set-level species calls, or any cross-species level comparison.

## 6 · The real blocking item (it is not notebook 49)

**Four of notebook 112's seven species lines were computed on the retired SCF13 coordinate**, from the notebook's
own source map: `S1 level` and `S4 atlas level` ← notebook 68; `S2 F-test` ← notebook 67; `S3 zonation class` ←
notebook 40. Notebooks 67 and 68 read `results/pt_reconstruction_v2/primary_pt_coordinate.csv` (12,272 structures);
notebook 40 reads `results/minimal_pt_scfates/`. The package export is 11,947 units, 11,383 in common support.

Consequences:
- Any gene reaching k ≥ 5 must be called by at least two SCF13 lines. All 23 filtered core genes carry 3 or 4 SCF13
  calls, and 20 of 23 carry only one package-era line. The gene list would probably survive — all 23 are also
  called by notebook 71's package F-test — but **the convergence statistic itself (the enrichment ladder and the k
  threshold) is not a package-era number.**
- The four top-ranked lines in notebook 113's species ranking are all SCF13-era; only handoffs (0.397, the lowest)
  is package-based. A package F-test exists and was not used: notebook 71's gives 3,584 calls at q < 0.05 against
  notebook 67's 3,889, Jaccard 0.62, with 1,020 SCF13-only calls. The lines are not interchangeable.
- There is **no package replacement at all** for the notebook 68 level/atlas lines or the notebook 40
  zonation-class line.

So the blocking item is: rebuild the species atlas on package-era lines, or state the convergence result as
SCF13-era. Notebook 49 is not blocking and never was.

## 7 · Audit record: what this page got wrong

Kept so the corrections are traceable.

| Claim as published | What the tables say |
|---|---|
| "nine species lines" | seven (S1–S7); nine is the AKI count |
| reversal list corroborated three independent ways | ~one; S5/S6/S3 are reversal-informative lines of the same core |
| 7 of 97's rewired-in-both genes on the reversal list | 11 |
| 97's rewiring separates reversed from conserved 0.761/0.213 | true on the position-contaminated total layer; 0.386/0.295 on the residual layer |
| F-test "most informative species screen", 0.88 | intervals overlap zonation class (0.746–1.019 vs 0.388–0.900); the species model adds nothing over expression (+0.010 / −0.031); no relabeling control exists for any species line; and the model's own fitted weighting is rejected by its own notebook |
| notebook 49 is the only pre-package number and is blocking | notebook 73 already re-ran it; four of 112's seven species lines are the real pre-package problem |
| 33 of 43 transporters is a hand count | notebook 47 writes it |
| "replication supports genes, not gene sets" | 13 of 83 sets replicate at BH ≤ 0.10; three in all four references |
| species anchored share rises to 0.870 at the highest k | 0.870 is k ≥ 5; it falls to 0.800 at k ≥ 6 |
| landmark gene-fold 0.93–0.97 | LRS-Z only; LCP 0.820 human, 0.647 human S1 |
| human neighbour bound "at most −1.2" | −1.4 per specimen |
| 11 handoffs "confirmed in microdissection and in the atlases" | 11/11 mouse, 9/11 human |
| conservation index "≈ 0.55" top 1% | 0.569 / 0.576; 0.54 was a CI bound |
| index "near zero" below the moderate bins | not estimable (NaN in every draw) |
