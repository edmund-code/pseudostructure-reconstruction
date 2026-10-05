# Draft v0.3 change log

`draft_v03.md` replaces `docs/paper/draft.md` (v0.2). It integrates `v02_changes.md` items 1–15, the
outline's revised claim, notebooks 41–46 and the referee report. Every number below was re-read from
the saved table named beside it during this round. Paths are relative to the results root
(`results/`) unless they start with `docs/`.

Short names for the main sources:
- **nb42**: `pt_revision_labels/` (`summary.json`, `transfer_agreement.csv`, `transfer_confusion_*.csv`,
  `glomerulus_attached_pt.csv`, `segment_fractions.csv`, `spillover_index*.csv`, `registration_qc.csv`,
  `integration_genes_by_class.csv`);
- **D4**: `pt_revision_labels/d4_coordinate_runs.csv`, `coordinate_robustness.csv`,
  `coordinate_robustness_funnel.csv`;
- **nb43**: `pt_revision_classes/reviewed/tables/` (`conservation_index.csv`,
  `pairwise_gradient_correlations.csv`, `class_counts_by_setting_and_contrast.csv`,
  `gene_classes_symmetric.csv`, `robust_pathways_reframed.csv`, `decisions.csv`);
- **nb44**: `pt_revision_state/tables/` (`A1_*`, `A2_*`, `B_*`, `C_*`, `D_rat_segment_proteome.csv`);
- **nb45**: `pt_revision_method/tables/` (`table_null_call_distribution.csv`,
  `table_cross_species_draws_summary.csv`, `anchor_our_cohort.csv`, `table_mechanism_coherence.csv`,
  `posthoc_null_calls_by_juvenile_imbalance.csv`, `table_aki_strategy_calls.csv`,
  `aki_injury_pathways.csv`, `aki_injury_markers_by_segment_and_depth.csv`, `decisions.csv`);
- **nb46**: `pt_revision_classes/crossfit/tables/` (`crossfit_conservation_index.csv`,
  `crossfit_class_counts.csv`, `crossfit_label_sets.csv`, `decisions.csv`);
- **S1/S2**: `pt_pathway_final/scfates/tables/table_S1_strategy_specificity.csv`,
  `table_S2_pathways.csv`.

---

## 1. Front matter

| Passage | Change | Source |
|---|---|---|
| Header | v0.3 note, marker legend and standing caveat | – |
| Titles | Three options on the revised claim; option 4 withdrawn as an overclaim (referee B39) | outline.md |
| Abstract | Rewritten (≈ 270 words) around the outline claim. R2 is conditional ("specific only when the specimens within a group shared their positional biology"). "As the split-plot design predicts" removed. Conservation index defined in words; values 0.17 (0.12–0.22), 0.28 (0.22–0.33), external 0.27–0.28, cross-fitted ≤ 0.31; classes 388 zonated = 196 / 75 / 41 / 76; spread 0.58; 26 of 36 pathways | nb43 `conservation_index.csv`, `class_counts_by_setting_and_contrast.csv` (symmetric, union); nb46 `decisions.csv`; nb45 `table_cross_species_draws_summary.csv` (100 draws); notebook 39 posthoc amplitude; nb43 `robust_pathways_reframed.csv` |
| Introduction, protein maps | Breljak et al. 2016 cited as a discordance for OAT1/OAT2 (citation sentence 2; M8) | `docs/paper/revision_state_citations.md` |
| Introduction, methods context | Added Polański et al. 2024 (Visium HD aggregation, B32), Van den Berge et al. 2020 and Roux de Bézieux et al. 2024 (B30), Murphy & Skene 2022 (B31; author list confirmed in Europe PMC) | Europe PMC by PMID |
| Introduction, "This study" | "Every zonated call was confirmed"; species-only calls need the other species flat in its own atlas (M3 fix 5); central claim restated | nb43 protocol |

## 2. Results

### R1

| Passage | Change | Source |
|---|---|---|
| Cohort and labels | Labels described as joint Leiden cluster memberships, not anatomy (M1 fix 5) | Methods 6.7 |
| Rat | Rat excluded: reliability 0.03 (S2 − S1) and 0.45 (S3 − early) < 0.5 | nb43 `decisions.csv` |
| Choice of coordinate | Rewritten per item 12: rule selected DPT13 (P5 0.80 vs 0.72; gap 0.14 vs 0.29); gate applied after DPT13 failed it on the matched-only test (NCDR 0.40); joint test DPT13 45 vs 3, 0 (NCDR 0.03), 26 robust, 15 core; authors kept scFates. "The coordinator applied" → "we applied" (referee B27) | D4 `d4_coordinate_runs.csv`; notebook 36 `selection_rule.csv` |
| How the choice matters | 22 / 36 robust on DPT13, 32 / 36 on SCF13-ED; direction agrees in 40 / 40; peak Spearman 0.71 (recomputed: 0.714 on 22); categories 13 / 3 / 0 / 1 / 5; within-species 37, 5 vs 16, 0 | D4 `coordinate_robustness.csv`, `coordinate_robustness_funnel.csv` |
| Segment labels (new) | Distances mouse 56 / 93 / 297 µm (ρ 0.65); human 131 / 148 / 228 µm corrected, 149 / 168 / 259 uncorrected (ρ 0.23); attachment at 2 µm (human 9 polygons, OR 1.0; mouse 46, OR 3.0) and 10 µm (OR 2.0 vs 4.5); transfer 0.71 / κ 0.54, 0.87 / 0.80, 0.92 / 0.87; S1 81%, S3 99%, S2 → Lake S3 65%; transferred human S1 fraction 0.44–0.49 vs mouse 0.25–0.34 | nb42 `summary.json`, `glomerulus_attached_pt.csv`, `transfer_agreement.csv`, `transfer_confusion_human_lake_primary.csv`, `segment_fractions.csv` |
| Integration genes (new) | 16% of conserved vs 2% of mouse-only (notebook 40 classes) | nb42 `integration_genes_by_class.csv` |
| Spillover and capture (new) | Spillover medians human 0.07 / 0.10 / 0.45%, mouse 0.14 / 0.13 / 0.29%; zero-bin polygons outside the bin-covered hull (0.2% and 1.2% of polygons inside the hull lack bins) | nb42 `spillover_index.csv`, `registration_qc.csv` |
| What the coordinate adds | Power, not specificity: 60 vs 16 joint calls (M5 fix 3); attenuation 15% / 7% (σ 0.10: 14% / 8%) | S1; workstream 1 final checks |

### R2 (rewritten; item 14)

| Passage | Change | Source |
|---|---|---|
| Heading and framing | Average-level screens not specific; positional screens specific only conditionally. Split-plot argument described as motivation, not proof (M6 fix 1) | `docs/paper/r2_revision_results.md` |
| Cohort table | Unchanged values, matched and joint tests labelled per row (B8) | S1 |
| Joint test across coordinates | scFates 60 / 2, 0; DPT13 45 / 3, 0; SCF13-ED 66 / 1, 0 | D4 |
| Many-draw nulls (new) | Median whole-PT GSEA calls 13–40 per pool; cross-species NCDR 1.34; step screen median 0 (human), 5 and 11 (mice), cross-species 15 vs 0.5 (NCDR 0.05), 80% of draws meeting the rule; our relabelings at the 54th percentile (steps) and 70th–100th (average-level); calibration 1.0% / 0.6% and 3.8–4.8% | nb45 `table_null_call_distribution.csv`, `table_cross_species_draws_summary.csv`, `anchor_our_cohort.csv` |
| NCDR is not FDP (new) | 71% of whole-PT GSEA species calls recur with all donors (base rate 3%) | nb45 `table_cross_species_draws_summary.csv` (precision 0.715) |
| Mechanism (new) | Offset coherence 0.26–0.50 vs positional 0.03–0.20; our cohort 0.013 | nb45 `table_mechanism_coherence.csv` |
| Where positional screens fail (new) | Age: median 22 and 62 null calls vs 1 when balanced (post hoc). AKI: 13 of 14 markers highest in S3; 3 injury Hallmarks; 77–80% of called pathways change most in S3; relabeled ÷ condition 0.35 (steps, 26 vs 12, 6) and 0.17 (T_spatial, 50 vs 10, 7); AKI coherence 0.42 | nb45 `posthoc_null_calls_by_juvenile_imbalance.csv`, `decisions.csv`, `aki_injury_pathways.csv`, `table_aki_strategy_calls.csv`, `table_mechanism_coherence.csv` |
| Disclosures | Many-draw nulls added after review, with a rule fixed before they ran | `r2_revision_results.md` |

### R3 (rewritten; items 9, 10, 13 and the outline)

| Passage | Change | Source |
|---|---|---|
| Ceilings table (new) | Within human 0.70 / 0.70–0.93; within mouse 0.63–0.88 / 0.73–0.93; human–mouse 0.09–0.20 / 0.13–0.24 | nb43 `pairwise_gradient_correlations.csv` (r_corrected) |
| Conservation index (new) | 0.17 (0.12–0.22), 0.28 (0.22–0.33), S3c 0.29 (0.23–0.35); external 0.27 / 0.28 / 0.27; rat withdrawn | nb43 `conservation_index.csv` |
| Cross-fit (new) | ≤ 0.31 (upper 0.38) for S2 − S1 and S3 − early; S3c − early up to 0.33 (upper 0.43); human ceiling 0.64–0.67 → 0.53–0.55; human amplitude ≈ 0.50 | nb46 `crossfit_conservation_index.csv`, `crossfit_class_counts.csv`, `decisions.csv` |
| Amplitude | SD ratios 0.58 / 0.47 / 0.84; axis ratios 0.15 (0.15–0.23), 1.96 (1.41–3.05), 0.14 / 1.66; Deming 0.16–0.17; projection caveat (B33) | notebook 39 tables |
| Tissue state | Heading replaced by the M4 wording (item 9). Donor at the median of 14 healthy KPMP donors (calibrated −0.45; predicted altered fraction 0.50); AKI mice score 0.33 above controls; Lake 0.69 / 0.50 / 0.35, ρ −0.41 / −0.24; mouse AKI 0.74 (r 0.83) and S3 0.41–0.44 (r 0.59); ZT4/ZT16 ratio 0.99–1.35, sterol 73%; fasting point estimate; PDK4 94% vs 7%; HMGCS2 low | nb44 `B_our_donor_placement.csv`, `B_calibration_scores.csv`, `A1_*`, `A2_*`; notebook 39 `lake_donor_scores.csv`, AKI ratios |
| Classes (new) | 7,407 probe-balanced; 196 / 75 / 41 / 76 / 712 / 6,307; per contrast 31 / 58, 13 / 28, 1 / 110; female 38 / 71, 31 / 78; strong mouse-only 1; v1 247 / 66 → A3b 59 / 141 → 41 / 76; cross-fit 11 of 12 unskewed (exception 17 / 4, other half 7 / 22) | nb43 `class_counts_by_setting_and_contrast.csv`, `decisions.csv`; nb46 `crossfit_class_counts.csv` |
| Shared and not shared | Conserved examples; transporters 33 / 43 (77%) vs 196 / 388 (51%) under symmetric classes, 30 / 54 (56%) vs 29% under notebook 40; OAT1/OAT2 discordance with Breljak (M8); SLC34A3 usage (Fernandes et al. 2026); species-only and reversal examples; class enrichments flagged as notebook 40 only | nb43 `gene_classes_symmetric.csv`; item 9 |
| Pathways | Funnel 60 / 58 / 51 / 36 / 26; categories 26 / 3 / 0 / 2 / 5 with names; "cannot recover a conserved core" (M7 fixes 1, 3); 8 probe-sensitive pathways named (M7 fix 5); beyond-flattening reduced to Estrogen Early/Late | nb43 `robust_pathways_reframed.csv`; S2 |

### R4

| Passage | Change | Source |
|---|---|---|
| Robustness | Pass counts 28 / 36 and 32 / 36 instead of 94% / 97% (C5); 8 lose q with equal probes (C6) | S2 |
| Coordinate (new) | 22 / 36 on DPT13, 32 / 36 on SCF13-ED, directions never disagree; notebook 12 DPT run historical only (M5 fix 4) | D4 |
| External | Relabel-called comparison p 0.016 / 0.033 / 0.19 / 0.21 (M7 fix 2); pathway p-values nominal; the 7 × 10⁻²⁰ and 2 × 10⁻¹³ values no longer quoted (M7 fix 4, partial) | notebook 37 `external_replication_summary.csv` |
| Additional donors (new) | D9 sentence (item 9) | nb44 `C_visium_*` |

### R5 (rewritten; item 10)

Every gene now carries its notebook 43 union class, and T_spatial rank is stated separately.

| Gene | Class stated | Other numbers | Source |
|---|---|---|---|
| Gatm | indeterminate | −0.78 ours, +0.51 Lake; 6.8 / 6.3 / 12.1 / 11.0; 23.5 vs 4,270 TPM; rat S3/S2 −3.7; top 0.1% T_spatial | nb43; `story_gene_evidence.csv`; nb44 `D_rat_segment_proteome.csv` (−3.70); notebook 12 `gene_statistics.csv` |
| Gamt | indeterminate (probe-imbalanced) | rat S3/S2 +1.6 (1.649) | nb43; nb44 D |
| Acadm | conserved | +2.1 / +2.3 vs +0.46; rat +0.7 (0.660) | nb43; `story_gene_evidence.csv`; nb44 D |
| Acaa2 | reversal (S2 − S1 and S3 − early) | −4.0 / −4.2; +0.58, +0.95 men; rat −2.7 (−2.676) | same |
| Acsm3, Crot, Nudt19 | indeterminate | – | nb43 |
| Hsd17b4 | mouse-only | – | nb43 |
| Gclc, Gclm | indeterminate (Gclc probe-imbalanced) | rat +2.8 / +1.9; GCLC +0.61 | nb43; nb44 D |
| Gss | conserved | rat +2.4 | nb43; nb44 D |
| Cyp2e1 | not classifiable | rat −5.5 | nb44 D |
| Dcxr | reversal (S3c − early, S3 − early) | +0.97 / +1.53 (+1.59 / +1.48); −0.46 / −0.54; −1.94 / −2.03; rat −1.1 | nb43; `story_gene_evidence.csv`; nb44 D |
| Ugt3a1 | reversal (S3c − early) | +1.0 / +1.5 vs −0.62 / −0.98 | same |
| Pah, Cyp24a1 / Igfbp4 / Acox2, Glyat, Nt5e, Cyp7b1 | conserved / mouse-only / indeterminate | – | nb43 |
| Rbp4, Aox1 | outside the classes (expressed in human only) | +2.46 / +3.46; 7 of 8 Visium; AOX1 +1.13 | `pt_zonation_classes/tables/one_species_expressed_genes.csv`; nb44 C |
| Slc22a6 | indeterminate | 970 TPM; rat −6.5 (−6.550) | nb43; nb44 D |
| Slc13a3 | conserved (mouse-only in S2 − S1) | – | nb43 |
| Sterol: Hmgcr / Ebp / Hmgcs1, Cyp51, Lss | conserved / mouse-only / indeterminate | time-of-day flag; rat-mRNA argument dropped (item 9) | nb43; nb44 A2 |
| Slc9a3 / Slc4a4 / Slc6a19 | human-only / conserved / indeterminate | – | nb43 |
| Psat1 | conserved (probe-imbalanced) | – | nb43 |

Wording changes in R5:
- Creatine is a large amplitude difference, not a human reversal (item 9).
- "Flat in human" withdrawn for Acadm.
- The glutathione story is demoted to supporting.
- UGT1A9 and Ephx1 are not in Figure 5. The Margaillan et al. 2015 sentence on UGT1A9 is kept.
- RBP4: bin diffusion and tubular uptake of filtered protein added as alternatives (referee B26).
- The Gamt rat-protein sentence no longer says "first protein-level support". It now reads "supports, at protein level".

### Discussion

| Passage | Change | Source |
|---|---|---|
| What is new | Three results restated: conditional specificity; conservation index; short list at stated evidence levels | items 12–14 |
| What the spatial data add (new) | Answers referee B36 | – |
| Prior work | Limbutara et al. 2020 and Song et al. 2023 added | Europe PMC |
| Flattening | Mouse-centric view; reference-free spread | notebook 39 |
| Implications | Fernandes et al. 2026, Rakhshandehroo et al. 2009, Aomura et al. 2026, Bignon et al. 2023, Heyman et al. 2010, Layton et al. 2019 | `revision_state_citations.md` |
| Limitations | Labels, classification, measurement and coordinate dependence. "Shared confounders" paragraph pasted verbatim from `revision_state_citations.md`, with the two PMID links rendered as author–year citations (Xiong et al. 2023; Rakhshandehroo et al. 2009). The PMIDs are in section 8 | item 9 |
| Disclosed decisions | Seven items, including the coordinate rule and the notebook 45 post hoc age analysis | items 12, 14 |

## 3. Methods

Built by `apply_v03_methods.py` from the v0.2 Methods (exact-match edits) plus new subsections.

| Section | Change | Source |
|---|---|---|
| Header | Lists the v0.3 additions | – |
| 6.1 | AKI kidneys also used in notebooks 44 and 45 | nb44, nb45 |
| 6.3 | Human `v2` µm fields assume 0.5 µm/px [VERIFY] | nb42 `summary.json` |
| 6.4 | Registration [VERIFY] kept. Notebook 42 QC added: centroid offset 0 px; coverage 97.4% vs 76.5%. Zero-bin polygons resolved: 16 and 90 inside the hull (0.2% and 1.2% of polygons there) | nb42 `summary.json`, `registration_qc.csv` |
| 6.7 | "What the labels are" paragraph; Song et al. 2023; θ = 6 rationale [VERIFY] | referee M1, B29 |
| 6.10 | Count splitting cites Neufeld et al. 2023 (B28). Override text rewritten. [PENDING D4] resolved with a three-coordinate table | D4 `d4_coordinate_runs.csv` |
| 6.12 | Why tradeSeq/condiments were not used (B30) | – |
| 6.14 | "Why relabeling is a valid null" → "Motivation (not a guarantee)". NCDR ≠ FDP. Notebook 45 paragraph added (from `r2_revision_results.md`, Methods paragraph, reformatted) | nb45 |
| 6.16 | [PENDING D4] resolved: DPT historical; coordinate robustness in 6.29 / Table S2 | D4 |
| 6.19 | "Dilution-free ratio" → "reference-axis (projection) ratio" (B33). Injury-marker [VERIFY] resolved with notebook 39's lists | notebook 39 source (`INJURY_MARKERS`, `STRESS_MARKERS`) |
| 6.20 | Marked superseded as primary by 6.26 | – |
| 6.21 | Rat excluded; pointer to revision datasets | nb43 |
| 6.23 | Notes on two-draw NCDR, NCDR ≠ FDP, gene-bootstrap intervals, sensitivity analyses | – |
| 6.24 | `fetch_census_pt_cells.py` added; "Notebook 40 not committed" → notebooks 35–46 committed (git ls-files) | git |
| 6.25–6.30 (new) | Label validation (nb42), ceilings and symmetric classes (nb43), state / donors / protein (nb44), cross-fit (nb46), coordinate robustness (D4), figure assembly (nb41) | notebook protocols |
| 6.31 | Data availability renumbered (was 6.25) | – |
| Wording | "added by the coordinator" and "(coordinator correction)" removed from 6.20 and 6.27 (B27) | – |

[PENDING D8] (v0.2 R3 and Discussion) is resolved by the M4 wording and notebook 44 results in R3.

## 4. Figures, tables and references

- **Legends.** Section 7 pastes the round 3 legends (`SCRATCH/figures/figure_legends.md`, commit e23a987),
  changed only where they conflict with tables or text:
  - Figure 2b, c: "Table S2" → Table S1, the strategy table (`table_S1_strategy_specificity.csv`).
  - Library tags: "(W) WikiPathways" dropped, because Table S2 contains no WikiPathways set.
  - Figure 4a: the "DPT" row is named as notebook 12 on notebook 03's DPT, historical (item 15).
  - Figure S1e: a note giving the final-list attenuation.
  - Figure S2: a note that DPT13 is specific on the joint test (D4).
  - Figure S4a: a projection-ratio note.
  - Figure S5: a note that verdicts predate the symmetric classes. The examples are headline genes
    only (Gatm, Gclm, Acsm3, Acox2); Gclc is "supporting" in notebook 38.
  - Figure S9a: S3c − early reaches 0.33 (upper 0.43), which the bracketed values omit.
  - Figure S10a: mouse glomeruli fail the size rule; the pixel-scale [VERIFY] is added.
- **Proposed Figure S11** (notebook 44 state figure), not yet in notebook 41.
- **Figure references in the text** follow round 3 numbering: label checks are Figure S10, cross-fit Figure S9, many-draw nulls Figure 2f–h, AKI Figure 2h, state Figures S4 and S11.
- **Table 1** is new.
- **Tables S1–S11** are renumbered to match the notebook 37 file names. S2 gains the coordinate
  ("DPT13-robust") columns, and S9–S11 are new.
- **References.** 14 added, each checked by PMID in Europe PMC:
  - Abedini 2024; Aomura 2026; Bignon 2023; Feola 2026; Fernandes 2026;
  - Limbutara 2020; Murphy & Skene 2022; Neufeld 2023; Polański 2024;
  - Rakhshandehroo 2009; Roux de Bézieux 2024; Song 2023; Van den Berge 2020; Wigger 2026.

  Every listed reference is cited in the text, checked by script.

---

## 5. Discrepancies found between sources (for the coordinator)

1. **Transporter share.** Item 9's 30 of 54 (56%) against 29% uses the superseded notebook 40 classes.
   Under the symmetric classes it is 33 of 43 (77%) against 196 of 388 (51%). R3 gives both, primary
   first. Only Slc17a1 is a mouse-only transporter under symmetric classes. OAT1 is indeterminate and
   NaDC3 conserved, so the referee's "18 mouse-only transporters including OAT1, NaDC3, SMCT2" no
   longer holds.
2. **"72%" (item 14, `r2_revision_results.md`)**: the saved precision is 0.7149, so the draft says 71%.
3. **Rat OAT1**: `revision_state_citations.md` gives S3/S2 −6.6; the table gives −6.5495, so the draft
   says −6.5.
4. **Gatm wording (item 9)**: "declines modestly at most" conflicts with Lake +0.51 (a rise). The draft
   says "changes little and not consistently: −0.78 in our sections and +0.51 in Lake".
5. **Gamt protein**: `revision_state_citations.md` calls rat Gamt "the first protein-level support".
   No literature check backs the priority claim, so it is softened.
6. **Cross-fit bound (item 13)**: "≤ 0.31 (upper 0.38)" holds for S2 − S1 and S3 − early, the decision
   contrasts. S3c − early reaches 0.33 (upper 0.43) in setting (ii). Stated in R3, the Discussion
   and Figure S9.
7. **Slc13a3**: earlier drafts and the v0.2 R5 called it mouse-only. The symmetric union class is
   conserved (S3c − early down in both); it is mouse-only in S2 − S1 only.
8. **Pkhd1**: `docs/paper/r3_revision_wording.md` lists it as human-only. Under the symmetric union it
   is indeterminate (the mouse atlas is not flat-compatible). The draft uses Sel1l3, Slc9a3 and Ace2.
9. **Gclc vs Gclm in Figure 5**: Gclc is a notebook 38 "supporting" gene with unequal probes (3 mouse
   against 6 human). Gclm is "headline" and probe-balanced. Both are indeterminate.
10. **Notebook 45 anchor**: notebook 45's re-implementation of the step screen gives 50 species calls
    in our cohort, against 16 in notebook 37. Relabelings are 0 and 0 in both
    (`anchor_our_cohort.csv`). R2 quotes notebook 37's 16. The difference comes from the
    implementation (unit-level exact F, not structure-level) and should be stated if the anchor is
    shown in Figure 2g.
11. **Abedini PMID**: `data/external/abedini2024_visium/SOURCE.txt` cites PMID 38514613, which is a
    different paper. The correct PMID is 39048792 (Nat Genet 56:1712–1724).
12. **θ = 6**: `human_vs_healthy_mouse/diagnostics/pass1_clustering_sensitivity.csv` (2026-09-08) is no
    longer written by notebook 03. Its reference partition matches θ = 2 (ARI 1.0), not θ = 6 (ARI
    0.42–0.46), so it cannot support the choice. θ = 6 is left [VERIFY].
13. **Figure legends (round 3)** refer to the strategy table as Table S2; corrected (section 4).

## 6. Numbers not verified against a saved table

1. **33 of 43 transporters (77%) and 196 of 388 (51%).** Counted by workstream 4 from
   `gene_classes_symmetric.csv` (probe-balanced, union class; Slc*, Abc*, Aqp* symbols). Reproducible,
   but no notebook saves it.
2. **Corrected human distances 131 / 148 / 228 µm** and "−13.5%". Arithmetic on notebook 42's values
   under the unconfirmed 0.44068 µm/px human scale. The uncorrected values (149 / 168 / 259) are
   saved.
3. **Human `v2` segmentation pixel size.** Whether it ran at 0.44068 µm/px is not recorded.
4. **UGT1A9 probe positions** relative to exon 1 (R5 [VERIFY]).
5. **"Human S1 ordering depends on read depth"** and the other coordinate limits are carried from v0.2
   (notebooks 35–36), not re-read this round.
6. **Workstream 1 attenuation (15% / 7–8%)** is carried from `SCRATCH/ws1/final_checks.md`, not a
   results table.
