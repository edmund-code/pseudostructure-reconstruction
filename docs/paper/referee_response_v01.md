# Response to the internal referee report (v0.1 → draft v0.3)

The report is `docs/paper/referee_report_v01.md` and the revised manuscript `draft_v03.md`. Table paths
are relative to `results/`. "D4" means `pt_revision_labels/d4_coordinate_runs.csv` and
`coordinate_robustness*.csv`.

## A. Major issues

| Item | Change made | Evidence | Residual risk |
|---|---|---|---|
| **M1** Human labels not validated; asymmetry where they are weakest | R1 says labels are joint cluster memberships. A new label-validation block covers glomerular distance, attachment and reference transfer. S3c − early is co-primary (index 0.29). R3 is repeated with cross-fitted transferred labels. Integration-gene enrichment in the conserved class is reported (16% vs 2%). The S2 − S1 mouse-only excess is gone under symmetric classes (31 vs 58) | Notebook 42 (`summary.json`, `transfer_agreement.csv`, `glomerulus_attached_pt.csv`, `integration_genes_by_class.csv`); notebook 46 (`crossfit_conservation_index.csv`, `crossfit_class_counts.csv`); notebook 43 `conservation_index.csv` | **High.** Human labels have weak support: transfer 0.71 (κ 0.54); 65% of human "S2" → Lake S3; distance ρ 0.23; only 9 human polygons attached at 2 µm. Fixes 2 (species-separate clustering) and 3 (R3 without integration genes) were not run; held-out-gene cross-fitting replaces them. Medullary-ray location of human S3 not tested. Human pixel scale unconfirmed |
| **M2** No same-species ceiling; rat | Disattenuated within-human, within-mouse and cross-species correlations for every dataset pair. Conservation index = cross-species ÷ geometric-mean ceiling: 0.17 / 0.28 / 0.29 ours, 0.27–0.28 external. Rat fails the pre-registered reliability rule (0.03, 0.45), and every rat-mRNA argument is withdrawn | Notebook 43 `pairwise_gradient_correlations.csv`, `conservation_index.csv`, `decisions.csv` | **Medium.** Human ceiling 0.65 rests on our two sections against Lake (7 donors). Intervals are gene bootstraps, not donor bootstraps. The ceiling is lower with transferred labels (0.53–0.55) |
| **M3** Species-only classes depend on unconfirmed, asymmetric flat calls | Symmetric classes: A3b (human floor and flat margin scaled by amplitude) + A8 (species-only needs the other atlas flat-compatible). A6/A7 female confirmation and all contrasts reported. Strong mouse-only genes: 1. Introduction corrected | Notebook 43 `class_counts_by_setting_and_contrast.csv`, `gene_classes_symmetric.csv`, `mouse_only_gradient_distribution.csv` | **Medium.** 6,307 of 7,407 genes are not classifiable. Lake has 7 healthy donors against 12 male mice. Flat calls remain section-level for human |
| **M4** Shared confounders; "tissue state ruled out" overclaims | Heading replaced with the recommended M4 wording. S3 numbers reported (AKI 0.41–0.44, r 0.59). Donor placed on the Lake injury axis (median of 14 healthy donors). Time-of-day zonation check (ZT4/ZT16). Fasting point estimate. PDK4/HMGCS2 as suggestive. "Shared confounders" Limitations paragraph added verbatim | Notebook 44 `B_our_donor_placement.csv`, `A1_*`, `A2_*`; notebook 39 `lake_donor_scores.csv`; `docs/paper/revision_state_citations.md` | **Medium–high.** Fasting is untestable with segment-resolved public data (GSE267280 is one pooled column per condition). Sterol zonation is not time-of-day stable. The platform calibration offset is estimated in mouse and assumed for human. Donor source, age and procurement are unrecorded |
| **M5** Coordinate selection outcome-driven; specificity coordinate-dependent | R1 and Methods 6.10 state that the pre-registered rule selected DPT13, that the gate was applied afterwards on the matched-only test, and that the authors kept scFates. Notebook 37 rerun on DPT13 and SCF13-ED: both specific on the joint test (45 / 3, 0 and 66 / 1, 0). The coordinate-free step screen is the primary specificity demonstration (16 / 0, 0), and the coordinate adds power. The notebook 12 DPT run is historical only (R4, Figure 4a note). Per-pathway coordinate robustness added to Table S2 | D4; S1 `table_S1_strategy_specificity.csv`; notebook 36 `selection_rule.csv` | **Medium.** The coordinate is a decision, not the rule's outcome; a reader may prefer DPT13. 14 of 36 robust pathways are not robust on DPT13 (direction never disagrees). DPT13 gives 37 mouse-vs-mouse positional calls (scFates 16) |
| **M6** Two null draws; split-plot does not predict specificity | Split-plot argument recast as motivation (R2, Methods 6.14). Many-draw nulls: 750 same-species and 100 cross-species 2 + 2 designs from Census and Lake. Calibration at p ≤ 0.01 reported. Mechanism (pathway coherence of offsets vs positional deviations). AKI positive control (detection passes, specificity fails). Specificity stated as empirical and conditional | Notebook 45 `table_null_call_distribution.csv`, `table_cross_species_draws_summary.csv`, `anchor_our_cohort.csv`, `table_mechanism_coherence.csv`, `table_aki_strategy_calls.csv` | **Medium.** External nulls test the step screen, not the smooth coordinate. Our cohort's NCDR still rests on two draws. 80% of cross-species draws meet the rule, exactly at the pre-specified 80%. Age-imbalanced male splits break the step screen (post hoc). NCDR is not an FDP |
| **M7** Pathway layer restates mouse zonation; replication not specific | Robust list reframed as "programs whose mouse zonation the human sections do not share" (26 of 36 summary of mouse zonation, 0 human-led). "Conserved core" reading dropped. Relabel-called comparison reported (p 0.016–0.21). Probe-sensitive pathways named (8). Two "beyond flattening" pathways dropped | Notebook 43 `robust_pathways_reframed.csv`; S2; notebook 37 `external_replication_summary.csv` | **Medium.** Fix 4 is partial: pathway-level Mann–Whitney p-values are now called nominal and the extreme values are not quoted, but no program-level or permutation comparison was run. Pathway-level replication beyond the global component (1–3 pathways) is not shown |
| **M8** Citations contradict claims | Breljak cited as a discordance for OAT1/OAT2, with mRNA–protein, antibody and interindividual explanations. AQP1 and OAT1–3 removed from "agrees with". Transport claim quantified (33 of 43 conserved vs 51%; v1 30 of 54 vs 29%). Conserved zonation separated from conserved usage (Fernandes et al. 2026) | Notebook 43 `gene_classes_symmetric.csv`; `revision_state_citations.md` | **Low–medium.** The transporter share is a workstream-4 count, not a notebook output. Human protein evidence (HPA) is presence-only |

## D. Missing analyses

| Item | Change made | Evidence | Residual risk |
|---|---|---|---|
| **D1** Anatomical validation of human labels | Done in part: glomerular distance, attachment, reference transfer, cross-fitted R3 | Notebooks 42, 46 | Species-separate clustering and R3 without integration genes not run. Mouse glomerulus polygons fail the size rule (descriptive). Human pixel scale unconfirmed |
| **D2** Same-species ceiling; rat decision | Done. Rat excluded | Notebook 43 | Ceilings rest on few human donors |
| **D3** Symmetric classification | Done (A3b + A8 primary; A6/A7 reported) | Notebook 43 | Most genes unclassifiable; class enrichments not recomputed under symmetric classes (the cited enrichments are notebook 40's) |
| **D4** Notebook 37 on DPT13, SCF13-ED and steps | Done. Joint test specific on all three; per-pathway robustness in Table S2 | D4; S1 | The "DPT13-robust" columns still have to be merged into `table_S2_pathways.csv`. Figure 4a's historical DPT row remains (labelled) |
| **D5** Many-draw null | Done | Notebook 45 | Segment resolution only; no coordinate in external atlases |
| **D6** Positive control (AKI) | Done. Detection passes (13 of 14 markers highest in S3; injury Hallmarks called; 77–80% of calls change most in S3). Specificity fails (relabeled ÷ AKI 0.35 steps, 0.17 T_spatial) because the two AKI mice differ in severity | Notebook 45 `table_aki_strategy_calls.csv`, `aki_injury_pathways.csv`, `decisions.csv` | Specificity under unequal positional biology is not guaranteed. Reported as a limit, not resolved |
| **D7** Tissue state and measurement quality | Done: donor injury placement, spillover index (small and similar), registration QC (zero-bin polygons lie outside the capture area) | Notebook 44 `B_*`; notebook 42 `spillover_index*.csv`, `registration_qc.csv` | Calibration offset transferred from mouse. H&E–Visium registration method unrecorded |
| **D8** Programs exposed to physiological state | Partly done. Time of day: FAO, glutathione, tyrosine and lead genes stable; sterol not. Fasting: whole-kidney point estimate only | Notebook 44 `A1_*`, `A2_*` | Fasted segment-resolved data do not exist, so a state contribution to sterol and PPARα programs is not excluded (Limitations) |
| **D9** More human donors with spatial data | Done; uninformative. Abedini Visium labels failed the S3 marker check. S1-high arms of UGT1A9, UGT3A1, ACOX2, SLC6A19 and the RBP4 rise replicate descriptively | Notebook 44 `C_visium_*` | No additional human donor validates the S3 side or any rising gene |
| **D10** Protein or in situ validation | Rodent side done with rat segment proteomics (S3/S2 ratios). Supports the rodent arm of creatine, glutathione, β-oxidation, Dcxr, OAT1 and Cyp2e1 | Notebook 44 `D_rat_segment_proteome.csv` | Human protein is presence-only (HPA). No human in situ validation; listed as future work (GATM, ACAA2, DCXR, UGT3A1, RBP4). Rat transcript data excluded (D2) |

## C. Mismatches (section C table)

| Item | Draft v0.1 said | Change made | Evidence | Residual risk |
|---|---|---|---|---|
| **C1** | Mouse snRNA 0.75 | 0.74 (fold-protected 0.745) | notebook 35 `p1_external_concordance.csv` | None |
| **C2** | Count split "without changing which pathways were called" | 145 vs 130 calls; 57 vs 59 of 66 | notebook 35 `g6_count_split.csv` | None |
| **C3** | "DPT was more stable" | Rule selected DPT13, disclosed (M5) | notebook 36 `selection_rule.csv`; D4 | See M5 |
| **C4** | No average-level screen reached FDP 10% | Pathway-score offset model reaches it (1,087 calls) and is explained | S1 | None |
| **C5** | 94% and 97% kept | 28 / 36 and 32 / 36 pass | S2 | None |
| **C6** | ≥ 70% effect under equal probes | Effect kept, but 8 / 36 lose q; named in R3 | S2 | Probe sensitivity affects xenobiotic and glutathione pathways |
| **C7** | Tissue state "ruled out" | M4 wording; S3 and Lake ρ reported; donor placed | notebook 39; notebook 44 | See M4 |
| **C8** | Correlation 0.14 (S2 − S1) only | S3c − early reported as co-primary (index 0.29; SD ratio 0.84) | notebook 43; notebook 39 | None |
| **C9** | "Every alternative rule" | All settings and contrasts reported; under symmetric rules S3 − early is 1 vs 110 | notebook 43 `class_counts_by_setting_and_contrast.csv` | None |
| **C10** | ACAA2 +0.58 "nearly flat" | Acaa2 is a reversal (+0.58; +0.95 in men) | notebook 43; `story_gene_evidence.csv` | None |
| **C11** | "Checked under pre-specified rules" for sterol, glutathione, peroxisomal, AOX1 | Each claim carries its symmetric class and its notebook 38 verdict. Sterol and glutathione are supporting, and AOX1 is supporting | notebook 43; `story_gene_evidence.csv` | Figure S5 verdicts predate the symmetric classes (noted in legend) |
| **C12** | Human "about half as large overall" | 0.58 (S2 − S1), 0.84 (S3c − early) | notebook 39 posthoc amplitude | None |
| **C13** | "Agrees with AQP1 and OAT1–3" | Removed; OAT1/OAT2 discordance stated; AQP1 indeterminate | notebook 43; Breljak et al. 2016 | None |
| **C14** | Rat omitted from R1 table | Rat stated and excluded by the reliability rule, with its concordance 0.10–0.13 | notebook 43 `decisions.csv` | None |
| **C15** | Robust ≫ uncalled | Relabel-called comparison added (p 0.016–0.21); p-values nominal | notebook 37 `external_replication_summary.csv` | See M7 |

## B. Minor issues

Minor issues 1–24 and 37–39 were addressed in v0.2 (`docs/paper/draft_changelog.md`) or by the R3/R5
rewrites above. Those addressed in v0.3:
- 27: coordinator language removed;
- 28: Neufeld et al. 2023;
- 29: integration caveat and Song et al. 2023 (the θ = 6 rationale is still open; see open items);
- 30: tradeSeq/condiments, in the Introduction and Methods 6.12;
- 31: Murphy & Skene 2022;
- 32: Polański et al. 2024 and registration QC;
- 33: projection ratio;
- 36: "What the spatial data add";
- 26: RBP4 alternatives (ambient RNA, bin diffusion, tubular uptake of filtered protein) in R5.

Still open:
- 25: UGT1A9 exon-1 probe positions [VERIFY];
- 34: AI-agent wording and the query archive (Methods 6.22 [VERIFY]);
- 35: Kim et al. 2011 remains the S3-vulnerability source.

## E. Central claim

The claim is now the outline's: "beyond a conserved core of strong markers, PT zonation is largely
species-specific". It rests on:
- a conservation index well below the pre-registered 0.5 bound against same-species ceilings, in our
  data and in public atlases alone;
- symmetric classes that show no excess of mouse-only genes;
- cross-fitted label checks.

Of the referee's four alternative explanations, three are addressed: label definition (M1, partly),
dataset heterogeneity (M2) and asymmetric flat calls (M3). Physiological state (M4) is bounded but
not excluded and is stated as a limitation. The single human donor remains the main limit on any
population statement.
