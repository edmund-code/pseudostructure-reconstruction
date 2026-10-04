# Change log: draft v0.1 → v0.2 (partial)

Base: `docs/paper/draft.md` (commit d47342e, which includes the notebook 40 corrections). Output:
`SCRATCH/ws4/draft_v02_partial.md`. The edits were applied by `SCRATCH/ws4/apply_v02.py`, which
checks that each target passage occurs exactly once; 29 edits were applied. R2, R3 and the abstract
were not rewritten; only the listed passages changed.

## Referee section C (each number checked against its source table)

| # | Location | Change | Source checked |
|---|---|---|---|
| 1 | R1 table | Mouse snRNA (coordinate column) 0.75 → 0.74, the fold-protected value; 0.75 is now given as the unprotected value | `nb35/p1_external_concordance.csv` (0.745 and 0.750) |
| 2 | R1 count split | "without changing which pathways were called" → circular arm 145 against 130 calls; robust 57 against 59 of 66 | `nb35/g6_count_split.csv` |
| 3 | R1 alternatives | "DPT was more stable" → the pre-registered rule selected DPT, with a metric table (see M5) | `nb36/selection_rule.csv`, `nb35/summary.json` |
| 4 | R2 | "No average-level screen" → "No competitive average-level screen"; the pathway-score offset model reached FDP 10% with 1,087 calls | `table_S1_strategy_specificity.csv` |
| 5 | R4 | 94% / 97% → pass counts 28/36 (detected in both) and 32/36 (sex-unbiased), which give 26 core | `table_S2_pathways.csv` |
| 6 | R4 | Equal probes: effect kept for 36/36, but 8 of 36 lose q ≤ 0.10 | `table_S2_pathways.csv` (`q_probe_equal`) |
| 7 | R3 | "Ruled out ... tissue state" softened; added mouse AKI S3 − early 0.41–0.44 (r 0.59), human/mouse S3c 0.28–0.33, Lake medians 0.69/0.50/0.35, WS3 ρ −0.42 (−0.66 to −0.10); marked [PENDING D8] | `aki_amplitude_ratios.csv`, `amplitude_ratios.csv`, `posthoc_reference_free_amplitude.csv`, `lake_donor_scores.csv`, `zonation_check/summary.json` |
| 8 | R3 | Added S3c − early disattenuated r 0.29 (0.23–0.34) | `posthoc_reference_free_amplitude.csv` |
| 9 | R3 | "every alternative rule" → A1–A7 and both region-matched contrasts; all-S3 contrast 77 against 111 | `class_sizes_by_setting.csv`, `class_sizes_by_contrast.csv` |
| 10 | R5 | ACAA2 no longer "nearly flat": +0.58, male donors +0.95, above the 0.5 floor; possible reversal (notebook 39) | `story_gene_evidence.csv`, `lead_gene_flattening_classes.csv` |
| 11 | Discussion | "Checked list ... under pre-specified rules" rewritten by evidence level: headline / supporting / class enrichment only / indeterminate in our data | `story_verdicts.csv`, `lead_gene_classes.csv` |
| 12 | Discussion (two places) | "About half as large overall" qualified: S2 − S1 0.47–0.58; S3c − early 0.84 | `posthoc_reference_free_amplitude.csv` |
| 13 | Introduction, R3, Discussion | Conserved-transport claim limited to SGLT1/2 (Vrhovac et al. 2015). Aqp1 indeterminate, OAT1 mouse-only, OAT2 conserved in our transcript data but conflicting with Breljak et al. 2016 protein data | `gene_zonation_classes.csv` |
| 14 | R1 table | Added rat row (0.11 / 0.10) and the notebook 34 label comparison (0.13; 55% direction agreement) | `nb35/p1_external_concordance.csv`; `pt_external_segment_validation/coordinate_check.csv` |
| 15 | R4 external | Added robust against relabel-called p = 0.016 / 0.19 (Census, male / female mice) and 0.033 / 0.21 (Lake) | `external_replication_summary.csv` |

## M5: coordinate selection disclosure

- R1, Methods 6.10, Discussion limitation 1 and open item 14 now state:
  - the pre-registered rule selected DPT13 (P5 0.80, P1 0.74 / 0.29, gap 0.14);
  - SCF13-ED failed it;
  - scFates itself would have failed it (P5 0.72, gap 0.29);
  - the coordinator applied C6 criterion 5 to every candidate after DPT13 was seen to fail it
    (NCDR 0.40), and scFates was kept as the fallback.
  Sources: `nb36/selection_rule.csv`; `SCRATCH/ws1/plan.md` (coordinator decision, 18:54).
- DPT was removed from the R4 robustness list. These places now carry [PENDING D4]: R1, R4, Methods
  6.10 and 6.16, and the Figure 4a legend.

## M4: tissue state

The text now says that injury and mouse AKI do not reproduce the S2 − S1 pattern. It no longer says
that tissue state is ruled out. Unmeasured differences that the external atlases share (procurement,
age, fasting, medication) are named. Marked [PENDING D8] in R3 and the Discussion.

## Queued v0.2 items

- **Item 1.** R1 attenuation: 10 / 26 pathways (median peaks 0.16 / 0.63); about 15% of the peak
  gap and 7% of the mouse-high count at σ = 0.089, and 14% / 8% at σ = 0.10 (`SCRATCH/ws1/final_checks.md`).
- **Item 2.** The partition-concordance sentences (ARI 0.979 / 0.981, 99.1%, 98.5%, ARI 0.955, 125)
  were added to Methods 6.7, and the glomerulus-source sentence (22 of 1,581; Jaccard 0.97) to
  Methods 6.16. Open items 6 and 10 are marked resolved.
- **Item 4.** The draft already used notebook 40 classes for reversals. Open item 20 now records
  that `r5_draft.md` and the novelty doc still carry the old labels; those files lie outside this
  task.

## Deliberately not changed (wait for the revision analyses, or outside this task)

- The abstract still says "human gradients were about half as large" (C12). Fix it when the abstract
  is rewritten.
- R2 table layout (minor 8, 2 × 2 of matched and joint), R3 Deming/orthogonal note (minor 11),
  female-biased background (minor 13), Dcxr "five references" (minor 24), and other section B items.
- R2 and R3 structure and the M1–M3, M6–M8 analyses.
