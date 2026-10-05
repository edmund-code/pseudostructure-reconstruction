# Figure plan for draft v0.4

> **Superseded** by `figure_plan_methods.md` (2026-10-05), which follows the authors' methods-paper outline.

There is one line per panel, with the source table relative to `results/`. "nb41 §3" means notebook
41's descriptive aggregates (no model, no test). Panels marked [WS6] wait for the running analysis. The
legends are in `draft_v04.md` (main) and `supplementary_notes_v04.md` (supplementary).

**General rules.**
- Every figure keeps the standing caveat in its legend.
- Pathway labels carry library tags (H, K, R) only.
- Peaks are shown in segment units only (0–1 S1, 1–2 S2, 2–3 S3).
- Robustness to a coordinate gets its own marker; fill stays reserved for core status (referee
  minor 4).
- Nominal pathway-level p-values that ignore overlap are not printed anywhere (referee mismatch 11).

## Main figures

**Figure 1 · Data and segment labels**
- 1a Design schematic with counts: `minimal_pt_scfates/pass2_nephron_filter_counts.csv` and `cross_species_pt_scfates.h5ad` (26,839 structures; 12,866 PT; per-specimen PT counts).
- 1b Glomerular distance by label and species (median, IQR; ρ 0.65 and 0.23): `pt_revision_labels/glomerulus_distance_by_label.csv` and `summary.json`. Plot human distances at 0.44 µm/px, with the uncorrected values in Figure S8a [VERIFY pixel scale].
- 1c S1 share among PT within 10 µm of a glomerulus against all PT, per species, with odds ratios: `pt_revision_labels/glomerulus_attached_pt.csv` (morphological, 10 µm).
- 1d Transfer confusion as row shares, human and mouse microdissection, with agreement and κ: `pt_revision_labels/transfer_confusion_human_lake_primary.csv`, `transfer_confusion_mouse_microdissection_primary.csv` and `transfer_agreement.csv`.
- 1e Coordinate overview: per-specimen segment medians and transitions: `pt_pathway_final/scfates/tables/segment_transitions_by_specimen.csv`, `minimal_pt_scfates/pt_segment_order_by_specimen.csv` and `pt_reconstruction_v2/primary_pt_coordinate.csv`.

**Figure 2 · Conservation against same-species ceilings**
- 2a Left: r* for every dataset pair, as points, with no median bar for the single within-human S2 − S1 pair (referee minor 8). Right: index and 95% interval, ours and external, with dashed lines at 0.5 and 0.8. Sources: `pt_revision_classes/reviewed/tables/pairwise_gradient_correlations.csv` and `conservation_index.csv`.
- 2b [WS6] Agreement by gradient strength: r and sign agreement in bins of |mouse gradient|, with per-bin ceilings, for three contrasts. The source is WS6's output table.
- 2c Index per label setting and gene half, with and without the 611 integration genes: `pt_revision_classes/crossfit/tables/crossfit_conservation_index.csv` and `pt_revision_addenda/tables/no_integration_genes_conservation_index.csv`. [WS6] Add the winsorisation and bootstrap variants when ready.
- 2d Noise-corrected human ÷ mouse spread per contrast (ours, external, transferred labels): `pt_zonation_amplitude/tables/posthoc_reference_free_amplitude.csv` and `pt_revision_classes/crossfit/tables/crossfit_class_counts.csv` (amplitude column).

**Figure 3 · Which genes are shared**
- 3a S2 − S1 scatter coloured by symmetric class (7,391 genes): `pt_revision_classes/reviewed/tables/gene_classes_symmetric.csv` with nb41 §3 means (as in the current Figure 3a).
- 3b Zonated-in-both rate (93% against 70%) and direction agreement (33/40 against 196/271), transporters against all genes: `pt_revision_addenda/tables/transporter_direction_agreement.csv`.
- 3c Mouse-only : human-only under four rule settings (absolute/scaled × with/without flat confirmation) and cross-fit: `pt_revision_classes/reviewed/tables/class_counts_by_setting_and_contrast.csv` and `crossfit/tables/crossfit_class_counts.csv`.
- 3d [WS6] Class composition by |mouse gradient| bin: `gene_classes_symmetric.csv` joined with WS6's per-gene gradients.

**Figure 4 · Pathway screens with two specimens per species** (replaces Figure 2)
- 4a Design schematic of the three balanced contrasts (current 2a), with the caption "motivation, not guarantee".
- 4b Six strategies: species and relabeled calls, with NCDR. Rows come from `pt_pathway_final/scfates/tables/table_S1_strategy_specificity.csv`; the coordinate-free step screen row comes from `pt_revision_method/tables/anchor_our_cohort.csv` (50; 0, 0). Referee minor 7: four to six rows only.
- 4c Same-species null calls per split, four pools; stars for our relabelings on the same rows as their strategies, in black (item 15): `pt_revision_method/tables/same_species_null_calls.csv`, `table_null_call_distribution.csv` and `anchor_our_cohort.csv`.
- 4d Cross-species NCDR per strategy (medians; 0.25 line; share meeting the rule), plus a precision inset (71% against 81%; base rates 3% and 2%): `pt_revision_method/tables/cross_species_draws.csv` and `table_cross_species_draws_summary.csv`.
- 4e AKI: condition and relabeled calls per strategy, marked pass or fail under the R4 rule and the positive-control rule: `pt_revision_addenda_pathways/tables/aki_strategies_under_both_rules.csv`.

**Figure 5 · Position-dependent pathways (primary list of 22)**
- 5a Member categories for the 22, with ★ for core under both coordinates and ◆ for the 20 also called by the coordinate-free step screen. Sources: `pt_revision_classes/reviewed/tables/robust_pathways_reframed.csv`, filtered by `pt_revision_labels/coordinate_robustness.csv` (robust under SCF13 and DPT13).
  - The ◆ list is not saved per pathway. Notebook 41 must reproduce notebook 48's mapping of the `pt_revision_method/stage_cache/anchor_our_cohort__*.csv` "called" indices to pathway ids. The two pathways not called are Formation Of Cornified Envelope and Metabolism of xenobiotics by CYP450.
- 5b Effect retained under seven checks, for the 22 (no historical DPT row; referee minor 4 and item 15), with markers for robustness under DPT13 and SCF13-ED: `pt_revision_addenda/table_S2_pathways_with_coordinate.csv`.
- 5c Gene-level agreement with the cortex atlas interaction, top 10% by T_spatial in black: `pt_pathway_final/scfates/tables/gene_level_statistics.csv` and `external_replication_summary.csv` (ρ only; no pathway p-values).
- 5d Peaks in segment units under the principal curve and DPT13, for human-high and mouse-high pathways: `pt_revision_addenda/tables/robust_peaks_in_segment_units.csv` (pooled transitions).

**Figure 6 · Lead genes** (replaces Figure 5)
- Cards: Gatm, Acadm, Acaa2, Dcxr, Ugt3a1, Rbp4 and Slc7a13 (reference). Gamt, Gclc and Cyp2e1 move to Figure S11 (referee section 4 R6; minor 5).
- Left of each card: notebook 37's fitted curves with specimen 12-bin means (as in the current Figure 5, `pt_pathway_final/scfates` cache, nb41 §3).
- Right of each card:
  - external S1/S2/S3 values: `pt_literature_deep/lake_cortex_segment_means.csv` and `story_gene_evidence.csv` (microdissection TPM; snRNA contrasts);
  - a new rat-protein S3/S2 column: `pt_revision_state/tables/D_rat_segment_proteome.csv` (Ugt3a1 not quantified; Rbp4 detected in S1 only).
- Titles:
  - symmetric class: `gene_classes_symmetric.csv` (Rbp4: "expressed in human only");
  - probes: `pt_literature_deep/probe_counts.csv`;
  - T_spatial rank: `pt_pathway_remodeling_scfates/gene_statistics.csv`.
- Rbp4's mouse curve is flat at zero, so state "not expressed" on the card.

**Table 1** (main): `pt_revision_classes/reviewed/tables/conservation_index.csv`, `crossfit/tables/crossfit_conservation_index.csv`, `pt_revision_addenda/tables/no_integration_genes_conservation_index.csv` and `class_counts_by_setting_and_contrast.csv`. [WS6] Add the sensitivity rows.

## Supplementary figures

**S1 · Coordinate construction and validation**
- S1a Embedding and trace (current 1b left), with the five-dimension note (referee minor 6): notebook 13 pass-2 embedding cache.
- S1b P1 bars (current 1c): `pt_reconstruction_v2/nb35/p1_external_concordance.csv` and `nb36/arm_metrics.csv`.
- S1c Identifiability (current 1d), with the label value of 95% shown beside the 91% (referee): `nb35/p2_absorption.csv`, `positive_control_recovery.csv` and `g6_count_split.csv`.
- S1d–f Equal depth, gene folds and numerical floor, and held-out prediction (current S1a–d): `nb35/p5_exposure_invariance.csv`, `p3a_gene_fold_within_segment_agreement.csv`, `numerical_floor.csv` and `nb36/p4_heldout_prediction.csv`.

**S2 · Pathway calls under three coordinates** (new; replaces the matched-only DPT diagnostic)
- S2a Robust and core counts and overlaps: `pt_revision_labels/d4_coordinate_runs.csv` and `coordinate_robustness_funnel.csv`.
- S2b Peaks in segment units under three coordinates: `pt_revision_addenda/tables/robust_peaks_in_segment_units.csv` and `peak_agreement_between_coordinates.csv`.
- S2c Within-species positional calls per coordinate: `d4_coordinate_runs.csv`.
- S2d Decomposition by screen: `pt_revision_addenda/tables/robust_list_decomposition.csv` and `pt_revision_addenda_pathways/tables/step_screen_reconciliation.csv`.

**S3 · Specificity details**
- S3a All strategies (current 2b): `table_S1_strategy_specificity.csv`.
- S3b Decoy FDP curves (current 2c): `matched_and_joint_tests_all_partitions.csv` and `pt_pathway_method_selection/conventional_signed_gsea_all_partitions.csv`.
- S3c Within-species controls (current 2d), retitled "Within-species controls" (referee): `within_species_controls.csv`.
- S3d Decoy inflation (current 2e): `decoy_inflation_by_statistic.csv`.
- S3e Split-plot Q–Q (current 4e): notebook 37 split-plot cache and `split_plot_summary.csv`.

**S4 · The 36 principal-curve robust pathways** (current S3b–e): `table_S2_pathways.csv`, `table_S4_programs.csv` and the fitted-curve cache. Mark the 22 primary pathways and give peaks in segment units (`robust_peaks_in_segment_units.csv`).

**S5 · Zonation amplitude** (current S4a–f, h–i; drop S4g): `pt_zonation_amplitude/tables/*.csv`.

**S6 · Physiological and tissue state** (new in notebook 41): `pt_revision_state/figures/fig_revision_state.pdf` panels a–c, plus `pt_revision_state/tables/B_lake_donor_injury_scores.csv` for the donor panel. Drop the Visium panel (F2).

**S7 · Genome-wide agreement with the cortex atlas** (current S7): `pt_literature_deep/gene_level_lake_contrasts.csv` and `genome_wide_lake_check.csv`.

**S8 · Segment-label checks** (current S10 plus additions)
- S8a Glomerular distance, uncorrected and pixel-corrected: `glomerulus_distance_by_label.csv` and `summary.json`.
- S8b Attachment at 2 and 10 µm: `glomerulus_attached_pt.csv`.
- S8c, d Confusion matrices for human, mouse microdissection and mouse snRNA: `transfer_confusion_*_primary.csv`.
- S8e Spillover index: `spillover_index.csv`.

**S9 · Index under transferred labels and without integration genes** (current S9a plus notebook 47): `crossfit/tables/crossfit_conservation_index.csv` and `pt_revision_addenda/tables/no_integration_genes_conservation_index.csv`. Show S3c − early with its 0.33 (0.43) maximum.

**S10 · Lead-gene verification matrix** (current S5), with added columns for symmetric class and rat S3/S2: `story_gene_evidence.csv`, `gene_classes_symmetric.csv` and `D_rat_segment_proteome.csv`.

**S11 · Supplementary gene panels**
- S11a–f Cards for Gamt, Gclc, Gclm, Gss, Cyp2e1 and UGT1A9, built as in Figure 6.
- S11g, h Whole-PT averaging for Slc22a6, Slc13a3 and Cyp24a1 (current S6).

## Removed from the figure set

| Panel | Reason |
|---|---|
| Current Figure 2f–h | Moved to Figure 4c–e |
| Current Figure 3d | 36-pathway categories; now Figure 5a on the 22 |
| Current Figure 4a "DPT (notebooks 03, 12)" row | Historical; not support (referee minor 4; item 15) |
| Current Figure 4d nominal p (2 × 10⁻¹³) | Overlap-ignoring p (referee mismatch 11); replaced by the overlap-preserving test in text |
| Current Figure S2 (matched-only DPT diagnostic) | Superseded by the reruns (referee section 4) |
| Current Figure S8 (notebook 40 classes, including S8d programs) | Replaced by Supplementary Note 3 tables; v1 enrichments do not hold (notebook 48) |
| Current S4g (lead genes against flattening) | Relied on the projection ratio |
| Current S10a uncorrected-only distances | Both scales shown in S8a |
| Visium donor panel (notebook 44 figure d) | F2: Visium donors only as a table in Supplementary Note 6 |
