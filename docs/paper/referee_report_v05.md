# Referee report on draft v0.5 (second round, same referee)

**What I reviewed.**
- The manuscript: `docs/paper/draft.md` (v0.5) and `docs/paper/supplementary_information.md`.
- The changelogs for v0.4 and v0.5 and `docs/paper/outline.md`.
- The tables of notebooks 47–50:
  - `results/pt_revision_addenda/`;
  - `results/pt_revision_addenda_pathways/`;
  - `results/pt_conservation_strength/`;
  - `results/pt_primary_pathways/`.
- The notebook 49 protocol (`SCRATCH/ws6/plan.md`).

Figures were not reviewed, because they are being rebuilt.

**Reviewer checks.** Two new checks are marked **[reviewer check]**. They use saved tables only and
modify nothing:
- `SCRATCH/review/v05_checks/reviewer_checks_v05.py`
- `SCRATCH/review/v05_checks/reviewer_checks_v05.csv`

---

## Verdict: minor revision

The revision answers the v0.3 report thoroughly:
- Every mismatch is fixed.
- The analyses I asked for were done and pre-registered: the strength curve, index sensitivity, the
  overlap-preserving replication test and the coordinate-robust primary list.
- Negative outcomes are reported plainly: "largely different genes" fails under module blocks, rule (ii)
  fails for S2 − S1, and the conventional-only group claim was dropped.
- The title is supported by the pre-specified concentration rule.

Three things stop me from saying "ready once lab records are added":
1. **The abstract's asymmetry claim is reproduced by ranking precision alone.** My check shows that
   unequal precision of the human and mouse strength references fully reproduces the 0.00 against 0.18
   asymmetry.
2. **The top-bin index is reported without the module-block intervals.** These intervals include 0, yet
   the same bootstrap was used to soften the all-gene wording.
3. **The weak-gene verdict is summarised incompletely.** For the S3 contrasts, the pre-specified rule
   "reproducible but not shared" actually holds.

All three need text changes. Only the first needs a cheap control, and only if the authors want to keep
the claim. With these fixes, the paper is ready once the lab records (F1) are added.

---

## 1. Status of the v0.3 issues

| v0.3 item | Status | Evidence in v0.5 |
|---|---|---|
| **F1** lab records | Not resolved (needs the user) | [VERIFY] kept for ethics, donor, strain, age, registration, pixel scale and deposition. R4 and Limitations now state that R4 depends on the unrecorded ages of the mice |
| **F2** selective Visium use | Resolved | No Visium statement in the main text (R6 says only that the labels failed). SN6 gives all 16 genes, pre-specified and post-hoc columns, one caveat for every row, and says that GATM, ACADM, ACAA2, DCXR, AOX1 and SLC22A6 fall there |
| **F3** title not shown; "mostly transporters" | Resolved, with new issues (section 3) | Notebook 49's pre-registered strength curve. Transporters are given as 33 of 196 (17%) |
| **M1** asymmetry depends on thresholds | Resolved | R3 table of four settings (247:66, 187:31, 59:141, 41:76); the failed cross-fit criterion is reported (R3; SN8 row 6); absolute against relative zonation is explained. But a new asymmetry claim replaces the old one (section 4, N1) |
| **M2** index sensitivity | Resolved | Winsorisation range 0.09–0.35; reliabilities matched; module-block bounds up to 0.54 lead to dropping "largely different genes" (SN8 row 18); donor-bootstrap bias acknowledged and leave-one-donor-out ceilings 0.63–0.68; single-pair ceiling and platform bias stated |
| **M3** pathway signal within segments | Resolved | Decomposition: beyond steps 22 of 22, spline-based steps 8. The coordinate-free step screen calls 20 of 22 (all 12 core); the two exceptions are probe-sensitive. Peaks in segment units; ρ 0.54 for the equal-depth refit (SN1) |
| **M4** forking path | Resolved | The primary list is the 22 robust under both coordinates (13/3/0/1/5; 12 core). Each coordinate is a sensitivity; SN8 rows 1 and 13 |
| **M5** AKI misread | Resolved | R4: the coordinate screen passes the R4 rule (0.17), fails the positive-control rule, and "the margin, not the rule" |
| **M6** "not specific" | Resolved | R4 heading "cannot be shown specific"; precision 71% against 81% |
| **M7** R5 claims | Resolved | Glutathione and Cyp2e1 moved to SN9 (GCLM opposite; Gss conserved; Cyp2e1 an expression difference); "weaker"; Acaa2 label caveat; Hsd17b4 and Slc9a3 label dependence; GATM −0.10 |
| **M8** framing and structure | Resolved | Biology-led abstract; the coordinate is a tool (R1, SN1); six sections, Box 1, SN1–10, a single deviation table. The R-section table in `outline.md` is stale (minor) |
| **M9** addenda | Resolved | Integration genes on symmetric classes (12% against 0%); transporter direction agreement (83% against 72%); naming rule for enrichments; overlap-preserving replication on the 22 |
| Mismatch 1 | Resolved | "Over-represented … 33 of 196 (17%)" |
| Mismatch 2 | Resolved | R2 gives 0.31 (0.38) and 0.33 (0.43); the abstract no longer quotes it |
| Mismatch 3 | Resolved | R3 "Failed pre-specified check"; SN8 row 6 |
| Mismatch 4 | Resolved | R4 AKI paragraph; SN2 |
| Mismatch 5 | Resolved | Four all-pass checks plus three counted (20/14, 21/20, 22 with 7 losing q); all match `table_S2` subset to the 22 |
| Mismatch 6 | Resolved | Methods: factor 0.881; SM25 "−11.9%; uncorrected values 13.5% too large" |
| Mismatch 7 | Resolved | SN6 labels "7 of 8" as post hoc |
| Mismatch 8 | Resolved | SN9 |
| Mismatch 9 | Resolved | SN1 "37 against 36" |
| Mismatch 10 | Resolved | R4 table lists both step designs (50 and 16; 15 of 16 within the 50) |
| Mismatch 11 | Resolved in text | No overlap-naive p-value remains; Figure 4d is pending the figures round |
| Mismatch 12 | Resolved | R3 gives "7,526 … 7,407 could be tested"; the class is a union |
| Smaller wording items | Resolved | 34% (SN1); 91% against 95% (SN1) |
| 2(d) glutathione | Resolved | Removed from the Discussion's list; GCLM opposite stated (SN9) |
| 2(d) Cyp2e1 | Resolved | "An expression difference, not a positional one" (R6, Discussion) |
| 2(d) Acaa2 | Resolved | S2 − S1 caveat, and the class under transferred labels |
| 2(d) "weaker", Hsd17b4, Gatm | Resolved | As above |
| 2(d) Dcxr and Visium | Resolved | Visium out of the main text |
| v0.1 minor 17 (S3c confirmed against outer-stripe references) | Resolved | R3 caveat; Limitations |
| v0.1 minor 20 (conventional-only group replication) | Resolved | Tested with overlap kept, p 0.41 and 0.24; claim dropped (SN5, SN8 row 20) |
| v0.1 minor 22 (GATM −0.10) | Resolved | R6 |
| v0.1 minor 37 (jargon) | Largely resolved | Box 1 glossary; internal arm names confined to SN1 |
| Female-mouse replication | Resolved | 5 and 3 of 22 (R5); 8 and 4 of 36 (SN5) |
| Figure minors (4–8) | Deferred | Figures round |
| [VERIFY]/[CITE] items | Partly resolved | Still open: OpenMidnight [CITE]; GSE277302's accession not confirmed in the full text; the θ = 6 rationale (Methods line 908) remains a [VERIFY] that is not a lab record, contrary to the changelog. The AI-search wording is stated, but the queries are not archived. UGT1A9 is stated as unchecked, which is acceptable |
| Kim et al. 2011; "first comparison" claim | Resolved | Reference removed; claim softened |

---

## 2. Numbers checked

**About 95 statements match their tables** (rounding aside). By section:

- **Abstract.**
  - Index 0.09–0.35; 0.17 and 0.28; module-block bound 0.54; public atlases 0.25–0.31
    (`index_sensitivity.csv`).
  - Top-1% index 0.55–0.58 and 0.09–0.37 (`strength_curves.csv`).
  - Direction agreement about 0.5 rising to 0.84–0.89 (`direction_agreement_by_strength.csv`).
  - Asymmetry 0.00 and 0.18 (S3 − early, single-species schemes); 100 draws; 22 pathways; 13 of 22
    mouse summaries.
- **R1.** Every count, label-check and coordinate value; unchanged and correct.
- **R2.**
  - Ceiling table; index table, including the module-block intervals (−0.07 to 0.40; 0.01–0.51;
    0.01–0.51), winsorisation ranges and public-atlas ranges.
  - Leave-one-donor-out 0.63–0.68; ceiling 0.13 (`human_ceiling_dependence.csv`).
  - All 18 cells of the index-by-strength table.
  - Rule (i): 0.46 (module-block interval 0.11–0.67) and 0.44 (0.10–0.63).
  - 1,308/3,021/2,343 and 5,509/3,358/2,391 genes; cross-species r* −0.13 to 0.01.
  - Post hoc 0.68 and −0.11.
  - Atlas pair 0.56/0.55 and 0.08/0.15 (SN10i).
  - Direction 0.57/0.48/0.58 → 0.89/0.84/0.89; 196 of 271.
  - Single-species 0.00 (−0.06 to 0.06; ceilings 0.70 and 0.71) and 0.18 (0.12–0.25).
  - Transferred labels 0.31/0.38 and 0.33/0.43; without integration genes 0.15 (0.20), 0.22 (0.27),
    0.25 (0.31) and at most 0.24 (0.33).
  - Spread 0.58/0.47/0.84; projection ratios 0.15–0.23 and 1.96; tertile 0.18.
- **R3.**
  - Class table; 7,526 and 7,407; reversals 28% of 271 (37 → 24/15, 5–11 against 15/21; 41 and 11–26
    against 16/24).
  - Transporters 33/196, 40/43 against 271/388, 33/40 against 196/271.
  - Asymmetry table; tyrosine q 2 × 10⁻⁶; human-only genes.
  - Integration-gene shares 24/196, 4/75, 1/76, 0/41, 359/6,307; PDK4 94% against 7%.
- **R4.**
  - Strategy table (158/193,255; 204/343,255; 0/46,2; 50/0,0; 16/0,0; 60/2,0).
  - 15 of 16; many-draw medians and NCDR; calibration; precision 71% and 81% (base rates 3% and 2%).
  - Coherence; age 22/62/1; AKI 13 of 14, 77–80%, 0.42, 0.17 with 10 and 7 of 1,429, 0.09–0.10, and
    26/12/6 with 0.35 (`aki_strategies_under_both_rules.csv`).
- **R5.**
  - Funnel 60/58/51/36; 26 under the second coordinate; 40.
  - Categories 13/3/1/5/0; 12 core; 10 against 12; peaks 0.50/0.95/0.54 and 1.77–1.95.
  - ρ 0.71 with 7 segment changes; coordinate-free 20 of 22 (12 of 12), the two exceptions, 28 of 36;
    spline 8; beyond steps 22.
  - Robustness 22/22/22/22, 20/14, 21/20, 22 with 7 (list correct).
  - Genes 0.67/0.73, 0.57/0.60 and 0.64/0.72; pathways 12/5/5/3.
  - Overlap null 1.62/1.44 with p 10⁻⁴; 0.004/0.007; 0.039/0.065; 0.26/0.13
    (`table_replication_primary_overlap_null.csv`).
- **R6.** Unchanged gene values, plus GATM −0.10.
- **SN10.** Tables a–i spot-checked row by row against `index_sensitivity.csv`, `strength_curves.csv`,
  `strength_bin_classes.csv`, `direction_agreement_by_strength.csv`,
  `posthoc_mouse_strength_full_lake.csv` and `human_ceiling_dependence.csv`. All correct.

**Mismatches and misstatements (9).** Items 1, 2 and 5 change wording in the abstract.

| # | Location | Draft says | Table says | File |
|---|---|---|---|---|
| 1 | Abstract; Discussion "Weaker genes" | "For the weakest half of genes, human gradients could not be shown to reproduce" | That holds for the 0–50% bin only. For S2 − S1 the human ceiling is not estimable over the weakest **80%** (5,509 genes). For S3 − early and S3c − early, the weakest 80% *does* reproduce in human (ceiling 0.61 [0.56–0.66] and 0.60 [0.55–0.66]) and is not shared (index 0.03 [−0.03 to 0.08]), so rule (ii) **holds** for these contrasts | `strength_decision_rules.csv`; SN10d |
| 2 | Abstract ("0.09–0.37 among … genes whose gradients reproduced"); R2 "In bins where both species' gradients reproduce (ceiling lower bounds ≥ 0.3), the index stayed between 0.09 and 0.37" | reproducing bins have an index of 0.09–0.37 | Bins whose gradients reproduce but whose index upper bound is ≥ 0.5 are classed "uncertain" and excluded. Examples: S3 − early 95–99% (ceilings 0.84/0.90, index 0.40) and every 99–100% bin. The range is partly produced by the classification rule | `strength_bin_classes.csv`; `strength_curves.csv` |
| 3 | R2 "Pre-specified verdict"; SN8 row 16 | verdict "agreement is confined to genes that are reliably zonated" | The pre-specified string is "…, and few genes are". The clause was dropped, and the title uses a different wording again, but SN8 says the verdict "was applied" | `decisions.csv`; WS6 protocol |
| 4 | R2 and SN10 note on dashes | "the half-reference's reliability was below 0.2" | In two dashed bins the fold-averaged human-reference reliability is 0.24 (S2 − S1 50–80%) and 0.23 (S3 − early 0–50%). The rule is "any dataset in a fold", so say "in at least one fold" | `strength_curves.csv` |
| 5 | Abstract "Asymmetry" | general statement, 0.00 against 0.18 | S3 − early only. For S2 − S1 the weakest-80% index is not estimable under either single-species scheme. These schemes were pre-specified as **secondary, "not used for the decision"** | SN10g; WS6 protocol |
| 6 | R5 "12 and 5 … 5 and 3 (BH ≤ 0.10)" | BH over the 22 is implied | BH was applied within the 36-pathway robust list (`q_within_robust`) | `table_S3_external_replication.csv` |
| 7 | R3 "tyrosine catabolism (Comt, Fah, Maoa)" | tyrosine catabolism | The enriched set is KEGG "Tyrosine metabolism". Comt and Maoa are catecholamine-degradation enzymes (the set also drives "Dopaminergic synapse"); only Fah is in tyrosine catabolism. Name it "tyrosine and catecholamine metabolism" and note that it rests on 3 genes | `class_pathway_enrichment_symmetric_probe_balanced.csv` |
| 8 | SN1 table "Median peak … 0.50 / 1.99; 0.95 / 1.90; 0.54 / 1.88" against R5 "mouse-high 1.77–1.95" | – | SN1 uses each coordinate's own robust list (36, 26, 51); R5 uses the 22. Label the SN1 row | `primary_pathway_list.csv` |
| 9 | Methods, θ = 6 | "[VERIFY: authors to state its original rationale]" | The v0.5 changelog says this marker was resolved and that only lab-record [VERIFY]s remain | `draft_v05_changelog.md` §5 |

---

## 3. Does the strength curve follow from notebook 49?

**Protocol.** It was pre-registered before any new result (`ws6/plan.md`, 2026-10-04 23:20).
- It reproduces notebook 43 exactly.
- Strength is cross-fitted: donor folds by hash, defined on fold A of each atlas and evaluated on fold
  B plus our data, with the folds swapped. Our data never enter strength.
- Strength is the *maximum* of the two species' ranks, so the selection does not favour shared genes.
- The design avoids the regression-to-the-mean problem of my v0.3 same-data check, which the authors
  report as a scheme in its own right.
- The decision rules and verdict wordings were fixed in advance.

This is well done.

**Per-bin ceilings.** The within-bin winsorisation and matched reliabilities are correct. Two points
need stating:
1. **The human half-reference is very noisy.** It is 3–4 Lake donors, some with 59–88 S2 nuclei. In
   the mid bins its reliability is 0.41–0.68, so the disattenuated ceilings are imprecise.
   **[Reviewer check]** On the full atlas, Lake's |S3 − S1| ranks agree with our human |S3 − S1| at
   Spearman only 0.09, against 0.34–0.46 for the mouse references with our mouse.

   Consequences:
   - The "max rank" strength axis is in practice driven by mouse strength. Genes zonated strongly in
     human only cannot be ranked reliably.
   - The title's "most strongly zonated genes" therefore effectively means "genes most strongly zonated
     in mouse (and in human where measurable)". This belongs in R2 or Limitations.
2. **"Reproducible" is not the same as "biological zonation."** Per-bin human ceilings can be raised by
   gradients that both human datasets share for non-zonation reasons: procurement, injury state, ambient
   signal (the paper notes TAL pickup in late PT). The paper's own "shared confounders" paragraph
   applies here. One sentence in R2 is enough.

**"Not estimable" bins.**
- The rule (reliability < 0.2 in any dataset) is pre-specified and the dashes are honest.
- The summary is incomplete (mismatch 1). At this half-atlas resolution only S2 − S1 cannot decide
  between "not reliably zonated" and "not shared". For both S3 contrasts the pre-specified rule (ii)
  holds: the weakest 80% reproduce within each species (0.61/0.78, lower bounds ≥ 0.55) and share
  nothing (index 0.03).
- The post-hoc full-atlas scheme (properly labelled post hoc in R2 and SN8 row 17) points the same way
  for S2 − S1 (ceiling 0.68, index −0.11).
- The applied global verdict ("confined to genes that are reliably zonated") is the most conservative
  reading. The abstract should say what the data show per contrast: for S3 contrasts, weak gradients
  are reproducible and not shared; for S2 − S1 the half-atlas is too small to tell, and post hoc the
  same holds.

**Cross-fitting of strength.**
- Correct.
- An approximate check of label-construction circularity is reassuring. Segment labels were built from
  the 611 integration genes and named with marker panels, and both are concentrated among the strongest
  genes. **[Reviewer check]** Among about 2% of the strongest genes (S3 − S1, full references), 29 of
  147 are integration genes and 8 are label markers. Removing them changes the raw cross-species r only
  from 0.585 to 0.568 (top ~5%: from 0.443 to 0.388).
- The authors should add the formal version, the strength curve without integration and marker-panel
  genes, as one row of SN10. I expect it to hold.

**The post-hoc label.** It is correctly applied to the mouse-strength full-atlas scheme. It is *not* applied
to the single-species schemes, which the protocol calls secondary ("to show the asymmetry"; "not
used for the decision"). These schemes are now in the abstract and Discussion as a finding (section 4,
N1).

**Is the title supported?** "Human and mouse proximal tubule agree mainly on their most strongly zonated
genes" — **yes, as a concentration claim**:
- Rule (i) holds in both primary contrasts under gene and module-block bootstraps (TOP − BOTTOM 0.46
  [module-block interval 0.11–0.67]; 0.44 [0.10–0.63]).
- Public atlases alone show the same pattern (0.55–0.56 against 0.08–0.15).
- Direction agreement rises to 0.84–0.89.
- Selection noise in the human rank dilutes, rather than creates, the concentration.

Two qualifications must accompany it:
- **(a)** The *level* of agreement in the top bin is poorly determined. Module-block intervals for the
  top 1% are −0.04 to 0.82 (S2 − S1), −0.02 to 0.74 (S3 − early) and 0.06 to 0.79 (S3c − early), and
  for the top 5% 0.03–0.63, 0.10–0.66 and 0.09–0.64. The abstract's 0.55–0.58 needs these intervals,
  or should be replaced by the top-5% contrast with its module-block interval. The paper used module
  blocks to retract "largely different genes"; it must apply them symmetrically.
- **(b)** Even in the top 1%, agreement is about half of same-species agreement. "Agree mainly" is
  right; "agree" alone would not be.

Title option 3 ("confined to strongly zonated genes") is too strong: bins at 80–99% keep an index of
0.2–0.4.

---

## 4. New problems introduced by the revision

**N1 (major; abstract, R2, R3, Discussion). The asymmetry is reproduced by unequal precision of the
strength references.**
- *The claim.* Genes weakly zonated in mouse share nothing with human (0.00), whereas genes weakly
  zonated in human "still follow part of the mouse pattern" (0.18).
- *The precision difference.* The human ranking comes from 3–4 Lake donors per fold, and its |gradient|
  ranks barely agree with our human data (ρ 0.09). The mouse ranking agrees at ρ 0.46. A noisy human
  ranking leaks genes that are truly strongly zonated, mostly conserved strong markers, into its
  "weakest 80%", and the bin then shows agreement.
- **[Reviewer check]** Degrading the mouse ranking with noise until its rank agreement matches the human
  atlas (ρ ≈ 0.10) raises the weakest-80% cross-species r from 0.03 to 0.15–0.19. The human-ranked value
  is 0.17. The whole asymmetry is thus reproduced without any biological difference (approximation:
  S3 − S1, raw r, full references).
- *Other problems.* The result exists only for S3 − early and S3c − early, and the schemes were
  pre-specified as secondary.
- *Fix.*
  - Remove the "Asymmetry" bullet from the abstract, the R3 bullet and the Discussion's "Binning by
    species" reading of the 0.18.
  - Keep the defensible half: genes weakly zonated in mouse but reproducible in both species are not
    shared (0.00, ceilings 0.70/0.71). It supports "human zonation outside the mouse pattern is
    reproducible".
  - Alternatively, run a matched-precision control in notebook 49 (rank mouse on a donor subsample, or a
    noise-degraded reference, matched to the Lake fold's reliability) and keep the claim only if 0.18
    survives.

**N2 (major; abstract, R2, Discussion). Asymmetric use of the module-block bootstrap.** See 3(a). Report
the module-block intervals for the strength bins in R2 and SN10d, and quote the top-bin value with its
interval.

**N3 (moderate; abstract, R2). The range of the weaker bins and the per-contrast verdict.** See
mismatches 1–3. Rewrite as follows: "In the 50–95th strength percentiles, where gradients reproduced
in both species, the index was 0.09–0.30. For S3 contrasts, the weakest 80% reproduced in both species
but shared nothing (0.03); for S2 − S1 the half-atlas could not tell." Avoid "not shared" for bins with
an index of 0.3–0.4; "weakly shared" is accurate. List the dropped clause of the pre-specified verdict
in SN8.

**N4 (minor). The strength axis is effectively mouse-defined** (section 3). Say so beside the title
claim.

**N5 (minor). Selection-matched replication.** Under the selection-matched null, the replication z of
the 22 halves (0.82/0.69; p 0.001/0.006). Against relabel-called pathways, program-level p becomes 0.27
and 0.37. R5's closing sentence ("replicate mainly because they are where mouse zonation is strongest")
is right. Quote the selection-matched values in SN5 so that the reader can see it.

**N6 (minor). Stale outline.** The R-section table of `outline.md` still says "scFates primary" and
"36 robust". It is not part of the manuscript, but it will mislead later workstreams.

**N7 (minor).** Mismatches 4 and 6–9.

**No regressions.**
- The Visium, glutathione, Cyp2e1 and AKI fixes hold throughout the main text and the SNs.
- No overlap-naive p-value reappears.
- Results stay at about 4,700 words.

---

## 5. Required before submission

1. Records (F1).
2. N1: remove the asymmetry claim, or support it with a matched-precision control.
3. N2: module-block intervals for the strength bins in the abstract, R2 and SN10.
4. N3 and mismatches 1–3: rewrite the weak-gene sentences and the verdict per contrast; add the SN8
   entry.
5. Add one SN10 row: the strength curve without integration and marker-panel genes. My approximate
   check suggests it holds.
6. Minor items N4–N7 and mismatches 4 and 6–9.

**Verdict: minor revision.** The central claims are now supported, pre-registered and honestly qualified.
What remains is one secondary claim in the abstract that ranking precision alone reproduces, the
symmetric reporting of an uncertainty that is already computed, and lab records.
