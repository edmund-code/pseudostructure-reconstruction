# Draft v0.4 change log

**Files.**
- `draft_v04.md`: main text, Methods, legends and references.
- `supplementary_notes_v04.md`: Supplementary Notes 1–9, supplementary figure legends and
  Supplementary Methods SM1–SM33.
- `figure_plan_v04.md`: the figure plan.
- `draft_v04_changelog.md`: this file.

**Sources.**
- The referee report is `docs/paper/referee_report_v03.md`.
- The decisions are in `docs/paper/outline.md` ("Decisions after the independent review").
- New results come from notebook 47 (`results/pt_revision_addenda/`) and notebook 48
  (`results/pt_revision_addenda_pathways/`; `docs/paper/r4_replication_and_classes.md`).
- Table paths below are relative to `results/`.

**Length.**
- *Results.* 4,780 words without tables (4,690 without the WS6 slots), against the target of about
  4,500. By section: R1 596, R2 707, R3 952, R4 833, R5 786, R6 900. R2 will grow when the WS6 results
  arrive.
- *Methods.* About 2,800 words in the main text; everything else is in SM1–SM33.
- *Discussion.* About 1,600 words.
- *Abstract.* 273 words without slots.

**Structure** (referee section 4).
- *Results.* R1 data and labels; R2 conservation against ceilings; R3 which genes are shared; R4
  pathway screens with two specimens per species; R5 position-dependent pathways; R6 lead genes.
- *Main figures.* Six.
- *Supplementary.* Supplementary Notes 1–9 (SN1–SN9 below) replace the coordinate, decoy, class-rule,
  state, external, Visium/rat, segmentation, deviation and gene-panel passages.
- *Glossary.* Box 1.

---

## 1. Coordinator instructions

| # | Instruction | Done |
|---|---|---|
| 1 | Mismatches 1–12 | Section 2 below |
| 1 | Section 2(d), the R5 table | R6 and SN9 (section 3 below) |
| 1 | F2, Visium | Every Visium lead-gene statement is out of the main text. R6 keeps one sentence saying the donors could not be used because the labels failed. SN6 has the full 16-gene table (pre-specified and post hoc columns), with the same caveat for every row; "7 of 8" is labelled post hoc |
| 1 | M1 | R3 table of four rule settings (247:66, 187:31, 59:141, 41:76). The absolute and relative scales are explained. The failed pre-specified cross-fit criterion is reported (17 against 4; other half 7 against 22). "Not skewed toward mouse" is removed from the abstract and Introduction |
| 1 | M5, AKI | R4 and SN2: the coordinate screen passes the R4 rule (NCDR 0.17; 10 and 7 of 1,429) and fails the positive-control rule; the average-level joint tests pass (0.09–0.10); the rule catches the step screen (0.35). Conclusion: the margin, not the rule. Source: `pt_revision_addenda_pathways/tables/aki_strategies_under_both_rules.csv` |
| 1 | M6 | The R4 heading, abstract and Figure 4 title say "cannot be shown specific by relabeling". Precision 71% against 81% (base rates 3% and 2%) is in the main text. Source: `pt_revision_method/tables/table_cross_species_draws_summary.csv` (0.7149; 0.8125) |
| 1 | M7 | R6 and Discussion: glutathione and Cyp2e1 are out of the positional list (SN9). "Weaker", not "no" counterpart. Acaa2 label caveat; Hsd17b4 and Slc9a3 label dependence; GATM second-atlas value −0.10 |
| 1 | Minor 2–11 | Section 4 below |
| 2 | Remove internal history | Main text has no "earlier drafts", "withdrawn", notebook, workstream or arm names (SCF13-ED, DPT13, CAC, P1w, A3b, A8). They remain only in the editorial note and the code-availability sentence. Box 1 is the glossary. SN8 is a single deviation table that replaces the scattered "(disclosed)" text and the seven-item list. The DPT13 disclosure is in SN8 (row 1) and SN1 |
| 3 | Abstract and titles | Rewritten around the single message. Pathway specificity is the second finding. The title options note that the WS6 verdict may change the wording |
| 4 | Slots | WS6: index sensitivity (R2, Methods, Table 1, Figure 2c), agreement by strength (R2, Figure 2b, 3d, Discussion, abstract). NB47 and NB48: filled (see below) |
| 4 | 22-pathway numbers | Section 5 below |
| 5 | Figure plan | `figure_plan_v04.md`: 6 main and 11 supplementary figures, one line per panel with its source table, and a list of removed panels |

## 2. Mismatches (referee section 1)

| # | Fix in v0.4 | Source |
|---|---|---|
| 1 | "Conserved core, mostly transporters" removed everywhere. Now: transporters are over-represented among conserved genes but make up 33 of 196 (17%); zonated in both 93% against 70%; direction agreement 83% against 72% | `pt_revision_addenda/tables/transporter_direction_agreement.csv` |
| 2 | The abstract no longer says "at most 0.31". R2 gives 0.31 (upper 0.38) for S2 − S1 and S3 − early, and 0.33 (0.43) for S3c − early. The Discussion says "at most 0.33" | `crossfit/tables/crossfit_conservation_index.csv` |
| 3 | The failed pre-specified "unskewed" criterion is reported in R3 and in SN8 row 6 | `crossfit/tables/decisions.csv` |
| 4 | The AKI reading is corrected under both rules (M5 above) | `aki_strategies_under_both_rules.csv` |
| 5 | "Six checks" is replaced by four all-pass checks plus three listed checks, with counts for the 22: detected in both 20 (14 keep q), sex-unbiased 21 (20), equal probes 22 (15 keep q). SN1 and Table S2 hold the 36-pathway values (34 and 35 keep the effect) | `table_S2_pathways.csv` subset to the 22 |
| 6 | Human distances are converted at 0.44068 µm/px, × 0.881 (−11.9%); the uncorrected values are 13.5% too large (Methods, SM25). Figure S8a shows both scales | `pt_revision_labels/summary.json` |
| 7 | RBP4 "7 of 8" is out of the main text. SN6 gives 7 of 7 (pre-specified) and 7 of 8 (post hoc, labelled) | `C_visium_lead_genes.csv`, `C_posthoc_centred_lead_genes.csv` |
| 8 | Glutathione is out of the Discussion short list and the metabolism implications. SN9 states that GCLM is opposite (atlas −0.76; men −1.05) and Gss conserved (+0.67) | `story_gene_evidence.csv`, `gene_classes_symmetric.csv` |
| 9 | SN1 explains DPT13 within-mouse 37 (rerun) against 36 (diagnostic) | `d4_coordinate_runs.csv` |
| 10 | Figure 4c and SN2 state that the percentiles use the re-implementation. "50 against 0 and 0 (16 with the spline-based step design)" is quoted, with the design difference | `step_screen_reconciliation.csv`, `anchor_our_cohort.csv`, `r4_replication_and_classes.md` |
| 11 | No nominal overlap-ignoring p-value appears anywhere. The Figure 4d (old) p is removed in the figure plan. The overlap-preserving test replaces the Mann–Whitney values | `table_replication_overlap_null.csv` |
| 12 | "7,526 genes have equal probe counts; 7,407 could be tested"; the class is the union, not a sum | `gene_classes_symmetric.csv` |
| – | Count-split inflation is 34%, not "about 30%" (SN1) | `nb35/g6_count_split.csv` (4.324 ÷ 3.234) |
| – | "Specimen differences are offsets" is retitled "Within-species controls" (SN2; Figure S3c) | – |
| – | 91% sign agreement is given with the labels' 95% (SN1; Figure S1c) | `nb35/positive_control_recovery.csv` |

## 3. Claims against evidence (referee section 2)

- **(a) Strength and the asymmetry.**
  - 72% direction agreement among genes zonated in both species is now in R2 and the abstract.
  - The agreement-by-strength curve is [PENDING WS6].
  - Reversals are stated as 28% of the genes zonated in both species. Their S2 − S1 label dependence
    is given from notebook 47: 37 → 24 and 15 with transferred labels; 5–11 in cross-fit halves
    against 15 and 21 reviewed. S3c − early reversals are stable (41; 11–26 against 16 and 24).
  - The asymmetry is handled as M1 (above).
  - Source: `pt_revision_addenda/tables/reversal_label_sensitivity.csv`.
- **(b) Index construction.**
  - R2 states the single-pair S2 − S1 human ceiling, the probable upward bias of cross-platform
    ceilings, and that gene-bootstrap intervals ignore module correlation.
  - The winsorisation range, the reliability mismatch and block or donor bootstraps are [PENDING
    WS6].
  - The missing evolutionary benchmark is in R2 and the Discussion.
- **(c) R4.**
  - The external calibration covers the step screen; this is stated.
  - Our cohort's coherence (0.013) is atypically low against human pools (0.14–0.20); this is stated.
  - The robust-list decomposition is reported (R5; SN1; Figure S2d):
    - spline-based steps 8 of 22 (9 of 36);
    - beyond steps 22 of 22 (35 of 36);
    - coordinate-free steps 20 of 22 (28 of 36), the new top-tier statement.
  - Precision 71% against 81% is in R4.
- **(d) R5 table** (now R6 and SN9).

  | Gene | Fix |
  |---|---|
  | Glutathione | Out of the short list and implications (SN9) |
  | Cyp2e1 | An expression difference, not a positional one (R6 summary; SN9; Discussion) |
  | Acaa2 | Reversal rests on S2 − S1 (S3c − early indeterminate); with full-gene transferred labels, a reversal in S3c − early. Visium in SN6 only |
  | β-oxidation and sterol | "weaker" (Discussion) |
  | Hsd17b4 | Indeterminate with transferred labels (R3; SN9) |
  | Gatm | The second atlas value −0.10 is added (R6). Visium −0.95 (10 of 10) is in SN6 only, per F2 |
  | Dcxr | The Visium statement is removed from the main text |

  Sources: `gene_classes_symmetric.csv` for the reviewed and `human_lake_centred__*` settings,
  and `story_gene_evidence.csv`.
- **(e) Coordinate.**
  - The primary list is the 22 robust under both coordinates (decision).
  - Peaks are given in segment units only: human-high median 0.50 / 0.95 / 0.54, mouse-high
    1.77–1.95 for the 22.
  - The equal-depth refit's peak ρ is 0.54 (SN1).
  - Sources: `robust_peaks_in_segment_units.csv`, `peak_agreement_between_coordinates.csv`.

## 4. Ranked issues (referee section 5)

| Item | Fix | Residual |
|---|---|---|
| **F1** records | Methods, Limitations and `open_items` keep [VERIFY] for ethics, donor, mouse strain, age and supplier, registration, human segmentation, pixel scale and deposition. New: the R4 conclusion and Limitations state that the mice's ages are unrecorded and that unequal ages would reopen R4 | Needs the user |
| **F2** Visium | As instruction 1 | – |
| **F3** title descriptor | "Mostly apical transporters" is removed. The new descriptor is "a small set of strongly zonated genes, in which transporters are over-represented" | The strength claim waits for WS6 |
| **M1** asymmetry | As instruction 1 | – |
| **M2** index sensitivity | Caveats added | [PENDING WS6] |
| **M3** within-segment drive | Coordinate-free step screen 20 of 22; spline steps 8; beyond-steps 22; peaks in segment units; SCF13-ED ρ 0.54 | – |
| **M4** forking path | 22-pathway intersection is primary (13/3/0/1/5; 12 core); each coordinate is a sensitivity (SN1 table) | – |
| **M5** AKI | As instruction 1 | – |
| **M6** "not specific" | As instruction 1 | – |
| **M7** R5 claims | As instruction 1 and section 3(d) | – |
| **M8** framing | Abstract, titles, Introduction "This study", R1 (coordinate as a tool) and the Discussion are restructured. The outline's central claim must still be updated in `outline.md` (coordinator) | – |
| **M9** addenda | Notebook 47 (integration genes, transporters, decomposition, peaks, reversals, Table S2) and notebook 48 (replication, enrichments, AKI, 50 against 16) are integrated. WS6 covers M2 and F3 | – |
| Minor 1 | Mismatches 2, 5, 6, 7, 9, 10, 11 and 12 (section 2) | – |
| Minor 2 | Item 17 (S3c confirmation against outer-stripe references) is a caveat in R3 and Limitations. Item 20 (conventional-only group replication) is stated in SN5 [VERIFY: no saved group test]. Item 22 (Census GATM −0.10) is in R6. Item 37 (jargon) is handled by the glossary, arm names removed and the coordinate moved to SN1 | – |
| Minor 3 | Female-mouse replication is given beside the male values: 8 and 4 of 36 (SN5); 5 and 3 of 22 (R5) | – |
| Minor 4 | Figure 5b uses a separate marker for coordinate robustness; the DPT row is removed (figure plan) | Figures round |
| Minor 5 | Figure 6 replaces Gamt with Rbp4. Gclc, Gamt and Cyp2e1 move to Figure S11 | Figures round |
| Minor 6 | Five-dimension note on the coordinate trace (Figure 1e legend; Figure S1a) | – |
| Minor 7 | Figure 4b has six rows; all strategies are in Figure S3a | – |
| Minor 8 | Figure 2a shows points, with no median bar on the single within-human pair | – |
| Minor 9 | [VERIFY] and [CITE] items kept: UGT1A9 probes (SN9), AI-agent wording (SM22), OpenMidnight, GSE267280 and GSE277302, PQ reporting (SN7), θ = 6 (Methods) | Needs the user |
| Minor 10 | The Kim et al. 2011 sentence is deleted from the Introduction and the reference removed. A Europe PMC search (generic terms) found no better dedicated source; Heyman et al. 2010's abstract does not name S3 | – |
| Minor 11 | "First comparison" is softened to "a structured literature search found no position-by-position comparison", with the Klötzer et al. 2025 caveat (Discussion) | – |

## 5. Numbers computed or filled this round (with sources)

**The 22-pathway primary list**, from `pt_revision_labels/coordinate_robustness.csv` and
`coordinate_robustness_funnel.csv`:
- 22 pathways (13 / 3 / 0 / 1 / 5); 12 core under both coordinates; 10 human-high and 12 mouse-high.

**Computed by workstream 4 by subsetting saved tables to the 22.** Notebook 41 or the coordinator should
confirm these:
- *Robustness:* `table_S2_pathways.csv` ratio and q columns (20/14, 21/20, 22/15, 7 probe-sensitive).
- *Pathway replication:* `table_S3_external_replication.csv` (12, 5, 5, 3 at BH ≤ 0.10 within the
  36-list).
- *Split-plot matched calls:* 7 of 22 (SN2).
- *Decomposition:* `pt_pathway_final/scfates/tables/matched_and_joint_tests_all_partitions.csv`
  (spline steps 8, beyond steps 22). The method reproduces notebook 47's 9 and 35 for the 36.
- *Coordinate-free step screen 20 of 22 (all 12 core):* reproduced notebook 48's mapping of the
  `pt_revision_method/stage_cache/anchor_our_cohort__*.csv` "called" indices to pathway ids.
  - The mapping reproduces notebook 48's 28 of 36, and DPT13's robust list gives 24 of 26.
  - The two exceptions are Formation Of Cornified Envelope and Metabolism of xenobiotics by CYP450.
  - The per-pathway list is not saved by any notebook.
- *Peaks in segment units:* `pt_revision_addenda/tables/robust_peaks_in_segment_units.csv`.

**Notebook 47 values quoted:**
- index without integration genes 0.15 (0.20), 0.22 (0.27), 0.25 (0.31); cross-fit ≤ 0.24 (0.33) and
  S3c ≤ 0.29 (0.38);
- transporters 40/43 (93%), 271/388 (70%), 33/40 (83%), 196/271 (72%), 33/196 (17%);
- integration-gene shares;
- reversals;
- peak correlations 0.71 (7 change segment) and 0.54 (10 change segment, 6 move > 0.25).

**Notebook 48 values quoted:**
- replication 1.53 / 1.24, p 10⁻⁴; 0.59; 0.005 / 0.015; programs 0.06 / 0.17; slope-free 0.42;
  selection-matched 0.004 / 0.035;
- AKI table;
- naming rule and enrichments (tyrosine q 2 × 10⁻⁶; human-only genes Ahcyl1, Hipk2, Wwtr1, Prkaa2,
  Taf4; conserved transport, bile, apical-surface and renin-secretion genes);
- 50 against 16 (15 of 16 within the 50; 28 of 36; 47 of 60).

## 6. Numbers not verified, or not verifiable from saved tables

1. **[PENDING WS6]**: every slot (index sensitivity, block and donor bootstraps, agreement by
   gradient strength, class by strength).
2. **Integration-gene shares.** Notebook 47 reports two scopes. The coordinator's summary used "all
   classified genes" (conserved 12.6%, reversal 6.5%, human-only 1.0%, mouse-only 0%, indeterminate
   5.5%). The draft quotes the probe-balanced scope, which matches the primary class table (24/196
   = 12%, 4/75 = 5%, 1/76 = 1%, 0/41, 359/6,307 = 6%).
3. **Overlap-preserving replication for the 22.** It cannot be recomputed, because the cached nulls hold
   collection-level draws only. The text says the test used the 36.
4. **Conventional-only pathways replicate better as a group** (v0.1 minor 20, p = 7 × 10⁻⁴). No saved
   table holds this test; it is marked [VERIFY] in SN5.
5. **Pixel scale.** The 131 / 148 / 228 µm values depend on the unconfirmed 0.44068 µm/px human scale.
6. **"3-week-old" pre-pubertal mice** (R4). Carried from `docs/paper/r2_revision_results.md`; the ages
   behind `posthoc_null_calls_by_juvenile_imbalance.csv` were not re-read.
7. **UGT1A9 exon-1 probes**: [VERIFY].
8. **Mouse ages in our cohort** (F1): unrecorded; R4's cohort conclusion depends on them.

## 7. Other changes

- **References.**
  - Kim et al. 2011 is removed (minor 10).
  - No new references were added. Le Hir et al. 1982, Bastin et al. 1990, Feola et al. 2026,
    Kaminska et al. 2026 and the glutathione and Cyp2e1 references are now cited from SN4 and SN9.
  - A script checked the list against the citations in the main text and supplement: 98 entries, all
    cited.
- **Supplementary Methods.**
  - Built from the v0.3 Methods by `build_sm_v04.py` (exact-match edits) and renumbered SM1–SM33.
  - SM31 adds notebooks 47 and 48, and SM32 is the analysis-to-notebook map.
  - Internal-role wording ("workstream", "coordinator", "(disclosed)", "v0.2") was removed.
  - The AKI rule and the −11.9% factor were corrected.
- **Figure-legend changes against round 3**, all recorded in `figure_plan_v04.md`:
  - Figure 2 → Figure 4, with panels 2f–h becoming 4c–e;
  - Figure 3d → Figure 5a on the 22;
  - Figure 5 → Figure 6, with the new gene set and a rat-protein column;
  - Figures S2 and S8 dropped;
  - new Figure S2 (coordinate sensitivity), S6 (state) and S11 (supplementary genes).
- **Still to update outside this draft** (coordinator):
  - the outline's central claim and R-section table;
  - `docs/paper/r2_revision_results.md` 72% → 71%;
  - `revision_state_citations.md` OAT1 −6.6 → −6.5 and the "first protein-level support" wording;
  - `r3_revision_wording.md` Pkhd1.
