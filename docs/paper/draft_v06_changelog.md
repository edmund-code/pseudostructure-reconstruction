# Draft v0.6 change log

**Files.**
- *Base:* the committed `docs/paper/draft.md` and `docs/paper/supplementary_information.md` (v0.5,
  3f4df61, with the coordinator's θ = 6 edit).
- *Outputs:* `draft_v06.md` and `supplementary_notes_v06.md`, built by `apply_v06.py` (exact-match edits
  only). The repository files were not edited.

**Review.** `docs/paper/referee_report_v05.md` (minor revision). All changes are text-only. New numbers
were read from `results/pt_conservation_strength/tables/` and
`results/pt_primary_pathways/tables/`.

**Lengths.**
- *Results:* 4,719 words (tables excluded), against 4,696 in v0.5. R2 grew by about 100 words for the
  per-contrast verdict, the top-bin intervals and the strength-axis caveat; R1, R3 and R6 were
  tightened.
- *Abstract:* about 315 words.

---

## 1. Coordinator decisions

| # | Decision | Change | Source |
|---|---|---|---|
| 1 | Remove the asymmetry claim (0.00 against 0.18) | Removed from the abstract, the R2 block, the R3 bullet and the Discussion "Binning by species" bullet. SN10g now has one sentence saying the contrast is not interpretable as an asymmetry because, in a reviewer check, unequal precision of the strength references reproduced it. The S12 legend says the same. The Methods call the single-species schemes pre-specified secondary schemes. No control was run | referee N1; `strength_curves.csv` |
| 2 | Module-block intervals wherever the top-bin index appears | Abstract: "0.55–0.58 in the top 1% (module-block 95% intervals between −0.04 and 0.82)". R2: top 1% −0.04 to 0.82, −0.02 to 0.74, 0.06 to 0.79; top 5% 0.03–0.63, 0.10–0.66, 0.09–0.64; "about half of same-species agreement". Discussion: "intervals as wide as −0.04 to 0.82". SN10d now has a module-block column for every bin. Legend 2b shows both intervals | `strength_curves.csv` (scheme max rank, q 0.01, bootstrap module) |
| 3 | Weak-gene verdict in full | R2 "The weakest 80%, per contrast": for S3 − early and S3c − early, rule (ii) holds (human ceilings 0.61 and 0.60 with lower bounds 0.56 and 0.55; mouse 0.78 and 0.80; index 0.03 [−0.03 to 0.08; −0.04 to 0.10]). For S2 − S1 it cannot be evaluated; the post hoc full-atlas check gives 0.68 and −0.11. The verdict is quoted in full ("…, and few genes are"), with the note that it holds as stated only for S2 − S1 at half-atlas resolution and that the title states rule (i) instead. SN8 row 16 records the dropped clause and the title wording. Abstract and Discussion give the per-contrast reading | `strength_decision_rules.csv`; `decisions.csv` |
| 4 | Replace the 0.09–0.37 range | Per-bin values with one explanatory clause: 50th–95th percentiles where gradients reproduce 0.09–0.30 (S2 − S1 from the 80th); 95th–99th 0.33–0.40, "weakly shared"; top 1% 0.55–0.58. Abstract and Discussion follow | `strength_curves.csv` |
| 5a | Replication BH family | R5: "BH ≤ 0.10 within the 36 principal-curve robust pathways" | `table_S3_external_replication.csv` (`q_within_robust`) |
| 5b | Tyrosine naming | R3: "the KEGG tyrosine-metabolism set (q = 2 × 10⁻⁶), but only through three genes: the catecholamine-degrading enzymes Comt and Maoa and the tyrosine-catabolic enzyme Fah". SN3 is changed to match | `class_pathway_enrichment_symmetric_probe_balanced.csv` |
| 5c | The 9 mismatches | Section 2 below | – |
| 5d | Label-circularity check | Left out. No committed table holds the strength curve without integration and marker-panel genes; the reviewer's numbers come from a scratch script. Notebook 47's all-gene index without integration genes is already in R2 | – |
| 6 | Changelog markers | Section 4 below | – |

## 2. The referee's mismatch list

| # | Fix |
|---|---|
| 1 | "Weakest half" is replaced by the per-contrast statement (abstract, R2, Discussion) |
| 2 | The 0.09–0.37 range is replaced by per-bin values (decision 4) |
| 3 | The verdict is quoted in full; SN8 row 16 records the change |
| 4 | The dash note now reads "reliability below 0.2 in at least one fold" in the R2 table note, SN10 and legend 2b. SN10d adds the fold-averaged half-reference reliabilities per bin |
| 5 | The asymmetry is removed from the abstract and main text (decision 1) |
| 6 | BH family stated as the 36 (decision 5a) |
| 7 | Tyrosine naming (decision 5b) |
| 8 | The SN1 peak row is labelled "each coordinate's own robust list: 36, 26, 51" |
| 9 | θ = 6: section 4 |

## 3. Other referee items

| Item | Change |
|---|---|
| Section 3 point 1 and N4: the strength axis is mouse-driven | R2 bullet *What "strongest" means*: the human half-reference's reliability was 0.10–0.78 below the 95th percentile, against 0.67–0.98 for mouse (`strength_curves.csv`, fold-averaged reliabilities). Also added to the Discussion main finding, a Limitations bullet ("Strength axis") and SN10. The reviewer's Spearman 0.09 against 0.34–0.46 is not quoted, because it comes from a reviewer check, not a committed table |
| Section 3 point 2: reproducible is not the same as zonated | R2 bullet pointing to the shared-confounders paragraph |
| Section 3(a) and (b) | "About half of same-species agreement" (R2, Discussion) |
| Title option 3 too strong | Replaced by "Human–mouse agreement of proximal tubule zonation is concentrated in strongly zonated genes" |
| N5: selection-matched replication | SN5 table: 22 mean z 0.82 and 0.69 (p 0.001, 0.006); program-level p against relabel-called 0.27 and 0.37 (`table_replication_primary_overlap_null.csv`). R5's closing sentence is unchanged |
| N6: stale outline | Not a manuscript file. For the coordinator: the R-section table in `docs/paper/outline.md` still says "scFates primary" and "36 robust" |
| Section 5 item 5: SN10 row without integration and marker genes | Not added; it needs an analysis (decision 5d) |

## 4. Remaining markers (corrects the v0.5 changelog)

The v0.5 changelog said that only lab-record [VERIFY] items remained. That was wrong.

**Open items that are not lab records:**
- **θ = 6 rationale**: `[VERIFY: authors to state its original rationale]` (Methods; SM7), as edited by
  the coordinator.
- **OpenMidnight citation**: [CITE] (Kaplan et al. 2025, model card only). The citation form is still to
  be chosen.
- **GSE277302**: cited as Nguyen et al. 2025 from Europe PMC's accession link and a matching abstract.
  The full text is not open access, so the accession is not confirmed in the paper itself. It needs
  confirmation.

**Lab-record [VERIFY] items** (unchanged from v0.5):
1. Mouse strain, ages (R4 depends on them), supplier, housing and euthanasia.
2. Human donor: source, age, kidney function, procurement, and whether the sections are serial.
3. Sex from records.
4. Ethics approvals.
5. Tissue processing and imaging.
6. Library preparation and sequencing depth.
7. Upstream mouse quality-control rules.
8. Human segmentation provenance and pixel scale (scanner metadata).
9. H&E–Visium registration method.
10. Segmentation inference software versions.
11. Data and code deposition.

**Counts.** No [PENDING] remains. [VERIFY]: 17 in the draft and 18 in the supplement, all in the
categories above. [CITE]: the OpenMidnight reference, plus the marker legends.

## 5. Numbers checked this round

| Number | Table value |
|---|---|
| Top 1% module-block intervals | −0.042 to 0.816; −0.022 to 0.739; 0.057 to 0.790 |
| Top 5% module-block intervals | 0.030–0.628; 0.097–0.660; 0.088–0.638 |
| 50th–95th bins | 0.239, 0.216 (S2 − S1, 80–95%); 0.092, 0.256, 0.274; 0.105, 0.215, 0.299, giving 0.09–0.30 |
| 95th–99th bins | 0.328, 0.396, 0.373, giving 0.33–0.40 |
| Weakest 80%, S3 contrasts | human ceilings 0.61 (0.56–0.66), 0.60 (0.55–0.66); mouse 0.78, 0.80; index 0.028 (−0.029 to 0.085), 0.029 (−0.040 to 0.101) |
| Half-reference reliabilities below the 95th percentile | human 0.10–0.78; mouse 0.67–0.98 |
| Selection-matched replication | 0.8237 / 0.6866; p 0.0013 / 0.0057; programs 0.2734 / 0.3659 |
| Tyrosine metabolism | q_joint 2 × 10⁻⁶; class genes Comt, Fah, Maoa |
