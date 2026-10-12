# The human vs mouse calls of notebook 71, downstream (notebook 77)

Part of the switch to the package pseudospace; see `pt-package-pseudospace.md`. Frozen on 2026-10-12, with the
coordinator's decisions of the same day (section 11), before any private-data run of notebook 77.

## Notebook 77 protocol: what the positional F-test calls contain, where they differ, and whether they replicate

**Why.** Decision D5 makes notebook 71's positional F-test calls (83 human vs mouse sets) the paper's species list. They
replace the retired 22-pathway list. The analyses that described the 22 must therefore rerun on the 83:
- what the calls summarise, by symmetric zonation class (43 §7; 47 e);
- where along the PT each call differs, and in which direction (47 e; 50; Figure 5d);
- how far a coordinate-free screen and the robustness checks support them (50 part 1; 47 c/e; R5);
- whether they replicate in external atlases beyond random sets (48 A; 50 part 2; 37 E3);
- the lead genes: ranks, segment values, verdicts and curves (37; 38; 41; R6; Figure 5e).

Notebook 72's verdicts bound the wording. C8 is not met (species registration unstable: read cross-species curves in
segment units), and C14 is not met (10 of 12 specimen × segment cells noise-dominated: no within-segment claim). Every
location below is therefore reported at segment resolution only.

**Question.** On the package units, labels and coordinate, what do the 83 calls summarise, in which segment does each
differ, which direction is higher there, and do they replicate beyond random gene sets and beyond a non-specific
comparator? Everything is descriptive: two male control mice and two cortex sections of one male human donor.
Specimens are the replicates. Human vs mouse only; AKI calls belong to notebook 76.

**Labels.**
- Notebook: `analysis/notebooks/77_pt_package_species_calls.ipynb`.
- Logic version: `77.species_calls.1`.
- Results: `results/pt_package_species_calls/`, with `tables/`, `figures/`, `stage_cache/` and `run_record.json`.
- Order: runs after 71–75. It reads their outputs and refits no gene model.

### 0 · Results of record and where they go

| Source | Result of record (old number; SCF13, Gaussian, notebook 13 units) | In 77 |
|---|---|---|
| 47 c, 50 | primary list: 22 robust under SCF13 and DPT13; 12 core | the 83 calls (D5); "core" kept as a gene-subset property (3.4) |
| 43 §7, 50 | categories of the 22: 13 summary of mouse zonation, 3 shared, 1 mixed, 5 unresolved, 0 human-led (36: 26/3/0/7) | 3.2 |
| 47 e, 50 | direction at peak: 10 human-high, 12 mouse-high | 3.3 |
| 47 e, 50 | peaks in segment units: human-high median 0.50 / 0.95 / 0.54 (SCF13 / DPT13 / ED), mouse-high 1.77–1.95; SCF13 vs DPT13 Spearman 0.71, 7 of 22 change segment | 3.3, at segment resolution |
| 50 | coordinate-free step screen 20 of 22 (all 12 core); both exceptions probe-sensitive; spline steps 8 of 22 | 3.4 (75's segment F-test, 74's ladder) |
| 47 e | beyond S1/S2/S3 steps: 35 of 36 robust (22 of 22) | 3.4 (74's test A) |
| R5 | robustness: four checks 22/22; detected in both 20/14; sex-unbiased 21/20; equal probes 22 effect, 7 lose q | 3.4 |
| 48 A, 50 | overlap-preserving replication, mean z 1.62 / 1.44, p_random 10⁻⁴; vs relabel-called p 0.004 / 0.007; programs 0.039 / 0.065; decision Census "specific", Lake "better than random only"; slope-free vs relabel 0.26 / 0.13 | 3.5 |
| 37 E3 | per pathway (BH ≤ 0.10 within list): 12 / 5 of 22 (male mice), 5 / 3 (female); genes ρ 0.68 / 0.73 vs matched random 0.57 / 0.60 (97.5th 0.62 / 0.66), top-expressed 0.64 / 0.72; global slope ours −0.81, Census −0.86, Lake −0.88 | 3.5 |
| 50 part 3 | conventional-only 184; Mann–Whitney 5.9 × 10⁻⁴; overlap-respecting D p 0.31 / 0.19: "not supported" | 3.6 (secondary) |
| 37, R6 | T_spatial percentile: Gatm 0.14%, Ugt3a1 0.52%, Acadm 0.83%, Acaa2 1.04%, Dcxr 1.05% | 3.7 |
| 38, R6 | ours human S3 − S1: Gatm −0.78, Dcxr +0.97, Rbp4 +2.46; stories creatine lead, fao lead, detox supporting, reversal lead | 3.7 |
| 41, Fig 5e/6 | lead-gene cards (curves, specimen bins, external profiles) | 3.7, 3.8 |
| 61, Fig 5d | programs by peak and direction, one program heatmap | 3.8 |
| 39 §9, 62, 48 8b/8c | beyond-flattening classes; non-1:1 orthologs; AKI rules; step-screen reconciliation | not carried over (section 4) |

### 1 · Inputs and guards

**Read, never refit.**
- 70: `package_pseudospace/hvm/pt_units.h5ad` through `package_inputs.hvm_inputs` (71's call: 13 ortholog-map `mt-`
  genes excluded). Run folder `runs/hvm/` for the integration features (73's `harmony_features`).
- 71 (`pt_pathways_package_pseudospace/`): `pathway_tests_human_vs_mouse_cameraPR_{species,swap_*}.csv`,
  `gene_statistics_human_vs_mouse.csv`, `redundancy_clusters_human_vs_mouse.csv`, `gene_sets_human_vs_mouse.json`,
  `universe_human_vs_mouse.csv`, `s9_winsor_lost_new.csv`; stage cache `gene_model_human_vs_mouse_species` (Δ(x), SE,
  group and specimen curves on the 101-point grid).
- 72: `tables/transitions.csv`; stage cache `g6_{a,c}` (per-set calls of the count-split arms).
- 73: `tables/43/gene_classes_symmetric.csv` and `gene_gradients.csv` (primary arm, mt-free);
  `tables/mt_in_library/43/gene_classes_symmetric.csv` (sensitivity).
- 74: `beyond_steps_sets_P.csv` and `pathway_tests_ladder_steps_on_spline_human_vs_mouse.csv`.
- 75: `gene_statistics_segment_human_vs_mouse.csv`, `pathway_tests_segment{,_joint}_human_vs_mouse.csv`,
  `pathway_tests_raw_tspatial_human_vs_mouse.csv`, `pathway_tests_deseq2_{wholept,segments}_gsea_human_vs_mouse.csv`,
  `deseq2_genes_human_vs_mouse.csv`, `programs_human_vs_mouse.csv`.
- Unchanged external inputs, as 37/38/50: Census PT segment pseudobulks (human, mouse), the Lake v1.5 h5ad,
  `pt_literature_deep/{probe_counts,story_gene_evidence,gene_level_lake_contrasts,lake_cortex_segment_means}.csv`,
  38's microdissection and snRNA aggregates, `xiong2023_sex_bias/PT-S{1,2,3}_fWTvsmWT.csv`.
- Read only, for the comparison: `pt_primary_pathways/`, `pt_revision_addenda{,_pathways}/`,
  `pt_revision_classes/reviewed/`, `pt_pathway_final/scfates/`, `pt_literature_deep/`.

**Guards (not cached).**
- G1: export contract and pass-1 fingerprint `efafc9ac2e99`; `hvm_inputs` gives 71's universe (11,055 genes), sets
  (1,513) and units (11,383).
- G2: 71's cached condition model reproduces its gene table (T_spatial to 10⁻¹²; z through `moderated_f`, 75's G3).
- G3: the calls are 71's (condition `call` true: 83); `overlap_programs` on them, with contributing genes = top 10% of
  the condition z, reproduces 75's `programs_human_vs_mouse.csv`.
- G4: `continuous_descriptors.logistic_transitions` on 71's positions and package labels reproduces 72's transitions
  (10⁻⁶).
- G5: 75's segment tables cover 71's universe and sets. 74's tables, if absent, leave their columns NaN with a note.
- G6: 73's classes come from its run with the confirmed fingerprint.
- G7: the Lake pseudobulk rebuilt over 71's universe equals 50's cached pseudobulk on the genes both share.

### 2 · Shared rules and the changes from the old definitions

**Shared.**
- Units: 71's units (package PT, common support). Labels: the package `segment_class`.
- Calls: 71's condition cameraPR calls (Up, FDR ≤ 0.05). Every count is given over sets and over 71's redundancy
  clusters; a cluster takes the label of its lowest-FDR called member.
- Contributing members: a set's members in the top 10% of 71's condition z; all members when fewer than 3 (75's rule).
- Programs: 75's, on the 83 calls.
- Descriptive per-unit expression: log1p(count / mt-free library × 10⁴), the paper-wide library (73's primary arm,
  75's heatmaps).
- Transitions: per specimen from 72; pooled = mean over the four specimens (47 e's rule); species = mean over its two.
- Segment units: `continuous_descriptors.segment_units` on [0, 1]. A mapping whose transitions are not increasing is
  undefined and reported so (72 C1: HUK1_MED1's transitions cross).

**Changes, each fixed now (K1–K11).**
- **K1 · The list (D5).** The 83 calls replace the 22. "Robust under both coordinates" has no counterpart: 71 fitted no
  second-coordinate gene model, and 77 fits none (D3; O2).
- **K2 · Gene statistic.** 71's condition z replaces Gaussian T_spatial wherever a gene ranking is needed: contributing
  members, top-10% genes, the selection-matched strata, lead-gene ranks. The NB T_spatial percentile is reported beside
  it (75's precedent).
- **K3 · Location (C8, C14; confirmed by O6).** Primary: the segment where the species differ most, from 75's
  coordinate-free segment model (3.3 b). The coordinate peak in segment units (47 e's rule) becomes a sensitivity, with
  pooled and species-specific transitions. Only the segment index is reported. Within-segment fractions and wording such as "late
  S2 near the S2→S3 transition" are dropped.
- **K4 · Second ruler.** The agreement of the coordinate peak with the coordinate-free peak replaces SCF13 vs DPT13.
- **K5 · Coordinate-free support.** 75's NB segment F-test replaces 45's Gaussian step screen. 74's steps on the spline
  and test A (beyond label steps) replace 37's `T_discrete_spatial` and `T_position_given_segments`.
- **K6 · Robustness.** Gene-removal checks run as cameraPR on 71's z without refitting (71's proteasome rule), with
  37's thresholds. 72's count split and 71's winsorizing are added per call. Refit-based checks (anatomical
  registration; cortical-like mouse S3) are not run (O1). Coordinate and registration robustness are notebook 72's
  results, and the text cites them:
  - C8, the registration gap, not met: species placement is unstable, so cross-species curves are read in segment
    units;
  - C14, the noise budget, not met: there is no within-segment claim.
- **K7 · Replication comparators (O3).** The old comparator structure is kept, with two comparators:
  - **Matched random sets.** Overlap-preserving within-strata gene relabelings keep set sizes, expression strata and
    overlaps, exactly as 48/50. They give p_random.
  - **Noise calls.** Old: the sets called under relabeling by 37's matched test without the VIF (21 not robust). 71
    itself yields only 1 relabeled call (relabelings 1 and 0), which cannot serve as a comparator. The analog is
    therefore the relabeling calls of the same positional statistic without its specificity correction: 75's
    cameraPR on rank-normal raw T_spatial. It calls 53 and 50 sets under the two relabelings; the comparator is
    their union minus the 83.
- **K8 · Classes.** 73's symmetric classes (primary arm) replace 43's reviewed ones; 73's secondary arm (mt in library)
  is the sensitivity, because 73 found class counts depend on the library.
- **K9 · Lead-gene display (C8, C14).** Curves are drawn on each species' own segment units. Claims are between-segment
  contrasts only.
- **K10 · Circularity flags.** Gatm, Cyp2e1 and Slc7a13 are package label-panel genes, and Gatm and Cyp2e1 are axis
  anchors (73 U2). Lead genes and contributing members carry label-panel, axis-anchor and integration-feature flags.
- **K11 · Reading of mouse-led calls (C12).** 72 found human amplitudes partly attenuated by positional noise. A call
  is called "mouse-led" only through the label-based classes, never through its coordinate curve.

### 3 · Analyses

**3.1 · The call table** (47 c's Table S2, rebuilt). One row per tested set: 71's n, direction, FDR and call for the
condition and both relabelings; redundancy cluster; program; every column of 3.2–3.5; and membership of the old 22, 36
and 60 and of 69's 97.

**3.2 · Categories** (43 §7, unchanged rules). Per call, each contributing member is classed over S2 − S1 and
S3c − early:
- *shared*: conserved or reversal in some contrast;
- *mouse-led*: zonated in mouse only, and never zonated in both;
- *human-led*: the converse;
- *mixed*: both mouse-led and human-led.

A call with fewer than 3 zonated contributing members is *unresolved*. Otherwise it is a *summary of mouse zonation*
(≥ 50% mouse-led), *shared* (≥ 50% shared), *human-led* (≥ 50% human-led) or *mixed*. Counts are given per call,
redundancy cluster and program, with the K8 sensitivity, and with the label-panel members counted.

**3.3 · Location and direction.** For 3.3 a and b, genes that are non-converged, separated or non-finite leave the
ranking; sets keep their remaining members.
- **(b) Primary, coordinate-free.** z_S = Δ_S / SE_S per gene and segment (75's condition steps fit). On these three
  points, `local_pathway_curves` gives each set's rank divergence D(S) against all genes. The peak segment is the argmax
  of D(S). The direction is the sign of the median member z_S there (> 0 is human-high).
- **(a) Sensitivity, the coordinate (47 e unchanged).** z(x) = Δ(x) / SE(x) from 71's cache; `spatial_descriptions`
  gives the peak x and the median member z there. The peak is converted to segment units with pooled transitions, and
  again with each species' transitions; the peak segment is the integer part.
- **Reported.**
  - Per call: segment and direction under b, a pooled, a human and a mouse.
  - Agreement of b with a pooled, in segment and in direction; this is the K4 replacement.
  - The segment × direction table under b.
  - Per program: the majority direction (37's rule) and the median member peak segment.
  - The median a-pooled segment units by direction, labelled "within-segment fraction not interpreted (C14)", for
    the comparison with 0.50 and 1.77–1.95.

**3.4 · Coordinate-free support and robustness.** Per call:
- **Specimens:** 71's relabeling FDRs.
- **Coordinate-free and step designs:** the 75 segment F-test call, its NB joint-test sensitivity, and 74's steps on
  the spline and test A.
- **Refits and caps:** retained in 72's count-split arms a and c; kept under all four of 71's winsor caps.
- **Gene-removal checks:**
  - the subsets: detected in ≥ 5% of units in each species; mouse sex-unbiased (Xiong 2023: no segment with
    |log2FC| ≥ 1 and padj ≤ 0.05); 37's 16 flagged genes removed; equal probe counts; the 8 proteasome subunits removed;
  - the test: `camera_pr` (0.01) on 71's condition z restricted to the subset, sets with 10–300 remaining members;
  - effect = mean z of members − mean z of the other subset genes;
  - retained: Up, FDR ≤ 0.10 over the subset's tested sets, and effect ratio ≥ 0.7 (37's thresholds).

**Core** = called and retained under detected-in-both and sex-unbiased (37's definition, with "robust" read as
"called"). **Probe-sensitive** = not retained under equal probes. For calls the segment F-test misses, the probe status
is listed (old: both exceptions probe-sensitive).

**3.5 · External replication** (48 A, 50 part 2, 37 E3). Definitions unchanged except K2, K7 and the units.
- **Ours.** Per specimen, S3 − early (early = S1 + S2) of the mean per-unit expression, by package labels. It is
  median-centred over informative genes; interaction = human − mouse. It is coordinate-free.
- **References.**
  - Human: Census cortex (6 donors) and Lake healthy cortex (donors with ≥ 50 nuclei per segment).
  - Mouse: Census male (primary) and female (sensitivity).
  - Informative genes: detected in ≥ 10% of units in each species, and mean log2CPM ≥ 1 in both external species.
- **Per-set score:** `JointReplicationZ` on a_g = sign(ours) × Welch t (human vs mouse donors). The slope-free version
  removes the uniform flattening from both sides (37's β).
- **Collections:** the 83 calls; core calls; call programs; the K7 noise-call comparator; comparator programs (same
  program rule).
- **Matched random sets (K7)** come from 9,999 within-strata gene relabelings (seed 50), applied to every set at once.
  Sizes, strata and overlaps are kept. Three nulls:
  - matched: the strata are external mean, our detection and our coverage (primary);
  - selection-matched: the strata × the quintile of z;
  - slope-free matched.

  Female references get the signed matched null only.
- **Decision (unchanged, Census male, signed, matched, 83 calls):**
  - *specific*: p_random ≤ 0.01, p vs the noise-call comparator ≤ 0.05, and program label permutation ≤ 0.05;
  - *better than random only*: p_random ≤ 0.01 alone;
  - *not replicated*: otherwise.
- **Per call (37 E3):** replicated = BH ≤ 0.10 within the 83 on p_joint, with effect > 0, for each of the four
  references.
- **Genes (37 E3):**
  - Spearman of ours against the external interaction on the top 10% of informative genes by z;
  - against 1,000 expression-matched random draws (median, 97.5th percentile; seed 12);
  - against the top 10% by mean expression;
  - sign agreement.
- **Global component:** the slope and R² of the interaction on the mouse gradient, for ours and for each reference.
- **Context, not rerun:** 71 S7 (atlas replication by method category).

**3.6 · Conventional-only sets** (50 part 3; secondary, kept by O4). Conventional calls are 75's whole-PT or S1/S2/S3 DESeq2
with signed GSEA (condition). Conventional-only means not among the 83. The test is unchanged:
- per-donor whole-PT log2CPM, median-centred;
- Welch t of human against male mouse donors, signed by 75's whole-PT log2FC;
- joint z;
- the collection null and D under joint relabeling, plus the overlap-naive Mann–Whitney.

The decision rule is 50's.

**3.7 · Lead genes** (37, 38, 41, R6). The genes are 41's nine (Gatm, Gamt, Acadm, Acaa2, Dcxr, Ugt3a1, Gclc, Cyp2e1,
Slc7a13) plus Rbp4. No gene is added or dropped after the run. Per gene:
- **Rank and calls:** percentile of 71's condition z among the universe, q, the NB T_spatial percentile, and the calls
  in which it contributes.
- **Segment values (38's "ours"):**
  - per-specimen S3 − S1 on 71's units and package labels, the species means and the interaction;
  - library: mt-free; 71's full-library Y as the sensitivity.
- **Verdicts:**
  - 38's `ours` and `ours_attenuation` tests are rerun with 38's claims and rules;
  - the attenuation a is recomputed by 38's rule. External contrasts come from `gene_level_lake_contrasts.csv`, so
    genes absent from it cannot be controls;
  - every external test is read from `story_gene_evidence.csv`; detection is on package units;
  - gene and story verdicts come from `replication_rules`.
- **Class and flags:** 73's union class (K8), and the K10 flags.
- **Program drivers** (descriptive, 37's rule on 75's programs): contributing genes in ≥ half of a program's calls, top
  8 by z.

**3.8 · Figures.**
- **Figure 5d.** Left: programs by peak segment (3.3 b) and direction, sized by calls. Right: one program heatmap,
  75's layout, rows = contributing genes (at most 40 by z).
  - Program: the largest program whose direction differs from that of 75's Figure 4e program (ties: lowest best FDR).
  - Columns: per species, 4 bins per segment on its own segment units (K3, K9), with specimen strips.
  - Side bars: z, 75's whole-PT log2FC and the 73 class.
- **Figure 5e.** Lead-gene cards, laid out as 41's:
  - left: 71's cached group and specimen curves (log CP10k), with 20 equal-count specimen bins (72's panel a rule),
    on species segment units;
  - right: the external S1/S2/S3 profiles, unchanged;
  - title: story, new verdict, class, probes, z percentile and flags.

  A gene outside 71's universe shows bins only.
- **Supplementary.** Location agreement (a against b); collection-null histograms (48's layout); the lead genes on 71's
  scale with transitions marked.

**3.9 · Export for the literature pass (O5).** `calls_for_literature.csv` has one row per call (83). It holds names and
labels only, no test statistics beyond the FDR, so that queries can be written without numbers from our data. Columns:
- `pathway_id`, `name`, `library`, the 71 FDR;
- redundancy cluster, and program with its name;
- category (3.2);
- direction and peak segment under 3.3 b, with the 3.3 a pooled segment beside them;
- lead genes: the call's contributing members ordered by z, at most 8, each with its 73 union class and K10 flags;
- which of 3.7's lead genes the call contains;
- whether it was in the old 22, 36 or 82 confident sets, and whether the old literature check (2026-10-04) covered it.

### 4 · Not carried over, with the reason

- **39 §9, beyond-flattening classes** (and 43 §7's column built from them). Notebook 39 is retired (73, decision 5).
  Its r = 0.152 is a mouse-referenced ratio, and 39's own post hoc diagnostics showed it depends on the reference species
  (0.15 against 2–3). "Beyond flattening" by that r is therefore not identified. The outline's need, pathway
  differences beyond a shared flattening, is met by 3.5's slope-free statistic, which estimates the flattening from
  external data on both sides.
- **62, non-1:1 orthologs along the coordinate.** Retired, for three reasons:
  - No outline claim, figure or draft sentence cites its coordinate results.
  - Its stage 1 and 2(b) rest on within-axis peak separations of 0.2 coordinate units and per-specimen spline zonation.
    Both lie below the package coordinate's noise: C14 (10 of 12 cells noise-dominated) and C8 (gap 0.37). 62 itself
    found that segment means explain 16 of its 17 separations.
  - It works outside the 1:1 space, so 71's gene models do not cover it, and a rerun needs new native-panel fits.

  Its coordinate-free parts (orthogroups, probe audit) stay unchanged inputs to 64/73 and 68/75.
- **48 8b** (AKI under both rules): 75 judges every AKI row by one rule; AKI belongs to 76. **48 8c** (50 against
  16): replaced by 75's segment F-test and 74's ladder.
- **37's funnel** (specimen-sensitive, stability runs, registration, region): the list is 71's by D5. Specimen
  sensitivity is 71's relabeling (1, 0). Registration and region are not refitted (O1); 72's C8 and C14 cover the
  coordinate.
- **47 c's DPT13 / SCF13-ED columns**: no second-coordinate gene model (O2).
- **Shape class, averaging loss, onsets and within-segment location**: 74's R4.
- **Literature status**: not computed here. 3.9 exports the calls for the coordinator's separate literature agent (O5).

### 5 · Outputs (`tables/`)

- `table_S2_calls.csv` (3.1).
- `call_categories.csv` and `category_counts.csv` (3.2).
- `call_locations.csv` and `location_summary.csv` (3.3).
- `call_robustness.csv` and `robustness_summary.csv` (3.4).
- `replication_overlap_null.csv`, `replication_decision.csv`, `replication_per_call.csv`, `replication_genes.csv` and
  `global_component.csv` (3.5).
- `conventional_only_overlap_null.csv` (3.6).
- `program_summary.csv` (direction, peak segment, drivers).
- `lead_genes.csv` and `lead_stories.csv` (3.7).
- `calls_for_literature.csv` (3.9).
- `compare_old.csv`, `material_change.csv` and `run_record.json`.

Figures as in 3.8.

### 6 · Comparison with the old results (descriptive; `compare_old.csv`)

Every row of section 0, old beside new. Also:
- the overlap of the 83 with the old 22, 36 and 60 and with 69's 97, over sets and over 71's redundancy clusters;
- the old 22 not called, and new calls outside the old 36;
- the lead genes' old and new ranks, values, classes and verdicts.

The changes come from units, labels, coordinate, gene model and the list together (71's D1–D8); 77 does not
apportion them.

### 7 · Material changes (pre-specified; stated in the summary and the ledger)

- **M1 · Categories.** Among resolved calls, "summary of mouse zonation" falls below 50% (old 13 of 17), or any call is
  human-led (old 0). R5's headline would change.
- **M2 · Direction by segment.** Under 3.3 b, the median peak segment of human-high calls is not S1, or that of
  mouse-high calls is S1. Old: human-high in S1 under every coordinate; mouse-high in S2–S3.
- **M3 · Location reliability.** 3.3 a (pooled) and b agree on the segment for fewer than half of the calls (old SCF13
  vs DPT13: 15 of 22). Then no per-call segment is reported; direction and programs only.
- **M4 · Replication.** The Census decision is not "specific", or p_random > 0.01 in either atlas.
- **M5 · Beyond flattening.** The slope-free comparison with the noise-call comparator gives p ≤ 0.05 in either atlas
  (old 0.26, 0.13). The calls would then carry replicated, set-specific differences beyond flattening.
- **M6 · Lead genes.** Any of these:
  - a lead story's verdict changes;
  - Gatm is outside the top 0.5% by z, or Acadm, Acaa2, Dcxr or Ugt3a1 is outside the top 2%;
  - Acaa2, Dcxr or Ugt3a1 is no longer a reversal in 73's union class.
- **M7 · Conventional-only.** 3.6's decision changes (old: "not supported once overlap is respected").

Everything else is descriptive. Nothing is tuned after a fit. A later check is labelled post hoc.

### 8 · Code-equivalence check (phase 2, untracked scratch script, before the package run)

The notebook's functions are fed the old inputs:
- notebook 12's T_spatial (43 §7's contributing genes) and 37's fit cache;
- 43's reviewed classes and 37's transitions;
- the old 22/36 and 50's comparator. They must reproduce:
- 43 §7's 26/3/0/7 (36) and 50's 13/3/1/5 (22);
- 47 e's scFates segment units;
- 50's `table_replication_primary_overlap_null.csv` and `conventional_only_overlap_null.csv` (same seeds);
- 38's story verdicts.

The notebook does not repeat the check, because it needs the old inputs.

### 9 · New code

`pseudospace/species_calls.py`, with `tests/test_species_calls.py` (synthetic only). No existing module changes.
- `member_category` and `call_categories`: 43 §7's rules. Tests cover each class rule and the fewer-than-3 case.
- `subset_retention`: 3.4's check. Tests: the full gene set gives a ratio of 1; a set with no remaining members is
  untestable.
- `segment_peaks`: 3.3 b, around `local_pathway_curves`. Test: a planted segment is recovered.
- Segment-unit mapping with the undefined rule, around `continuous_descriptors.segment_units`. Tests: monotone, and
  crossing transitions are refused.

### 10 · Compute (cold; no NB fit, no refit)

| Step | Time |
|---|---|
| inputs and guards | about 3 min |
| locations, categories, gene-removal cameraPR | about 3 min |
| Lake rebuild | about 5 min |
| 10 collection nulls of 9,999 draws | 10–20 min |
| lead genes and figures | about 3 min |
| **total** | **about 25–35 min** |

- The scratch equivalence check takes about 20 min, mostly cache hits.
- No refit or package-DPT model is run (O1, O2).

### 11 · Coordinator decisions (2026-10-12)

- **O1 · No refit-based robustness in 77.** Notebook 72 already reports coordinate and registration robustness, and the
  text cites it:
  - C8 not met: registration gap 0.37, so cross-species curves are read in segment units;
  - C14 not met: 10 of 12 cells noise-dominated, so there is no within-segment claim.

  77 keeps only the no-refit checks of 3.4.
- **O2 · No package-DPT gene models.**
- **O3 · The old comparator structure is kept (K7).** Matched random sets (sizes and expression strata, overlaps kept;
  48/50's null) give p_random. 75's raw-T_spatial relabeling calls are the noise-call comparator, because 71 itself
  yields only 1 relabeled call.
- **O4 · The conventional-only section (3.6) is kept**, as secondary.
- **O5 · No literature pass inside 77.** The coordinator runs a separate literature agent on the new calls afterwards.
  3.9 exports a clean table of the 83 calls with category, direction, segment and lead genes for it.
- **O6 · Confirmed.** The coordinate-free peak (75's per-segment Δ/SE, 3.3 b) is primary. The coordinate peak in
  segment units (3.3 a) is the sensitivity.

## Notebook 77 results

`analysis/notebooks/77_pt_package_species_calls.ipynb` (logic version `77.package_species_calls.1`, by the coordinator's
instruction) writes to `results/pt_package_species_calls/`. It followed the protocol above and refits no gene model.
- New module `pseudospace/species_calls.py` holds the rules of notebooks 37 E3, 38, 43 §7, 47 (e) and 50.
- Its tests are in `tests/test_species_calls.py` (9 tests). The full suite passes (529).
- Cold run 155 s; the executed notebook reran from the stage cache in 81 s.

**Guards.** G1–G7 passed:
- the export contract and the pass-1 fingerprint `efafc9ac2e99`;
- 71's universe (11,055 genes), sets (1,513) and units (11,383);
- 71's cached model reproduces its T_spatial and z;
- the 83 calls, and 75's 43 programs rebuilt exactly;
- 72's transitions;
- 73's confirmed run;
- the Lake pseudobulk, which equals notebook 50's on the shared genes.

**Code equivalence** (scratch, before the package run). Fed the old inputs (notebook 13 structures, reviewed labels,
notebook 12/37 statistics, 43's classes, 37's transitions, 50's comparator), the module reproduces:
- 43 §7's categories per pathway (26/3/0/7 of 36; 13/3/1/5 of the 22);
- 47 (e)'s peaks and segment units (|Δ| ≤ 10⁻¹⁵);
- 50's replication table to 2 × 10⁻¹⁶, with its decisions;
- 50's conventional-only table;
- 38's attenuation (0.859) and all 42 gene verdicts.

**Material changes.**

| Rule | Old | New | Material |
|---|---|---|---|
| M1 · categories | 13 of 17 resolved mouse summaries; 0 human-led | 73 of 82 (0.89); 1 human-led (KEGG Citrate cycle) | **yes** (the human-led clause only) |
| M2 · direction by segment | human-high in S1; mouse-high S2–S3 | median segment human-high S2 (index 1.5: S1 14, S2 4, S3 18); mouse-high S3 (S1 17, S2 4, S3 26) | **yes** |
| M3 · location reliability | 15 of 22 | 64 of 83 (0.77) | no |
| M4 · replication | Census specific; p_random 10⁻⁴ both | Census and Lake specific; p_random 10⁻⁴ both | no |
| M5 · beyond flattening | slope-free p vs comparator 0.26 / 0.13 | 0.35 / 0.63 | no |
| M6 · lead genes | stories lead/lead/supporting/lead; Gatm top 0.14%; Acaa2, Dcxr, Ugt3a1 reversals | stories unchanged; Gatm top 0.11% (z); Acaa2 and Dcxr no longer reversals in 73's classes | **yes** (classes only) |
| M7 · conventional-only | not supported once overlap is respected | same | no |

**The list.** The 83 calls fall in 77 redundancy clusters and 43 programs. 12 of the old 22 are among them (Jaccard
0.13), as are 24 of the old 36, 41 of the old 60 and 82 of notebook 69's 97. The 10 of the old 22 that are not called:
Drug ADME, O-linked Glycosylation Of Mucins, Formation Of Cornified Envelope, Aspirin ADME, Chemical carcinogenesis,
Retinol metabolism, Metabolism of xenobiotics by cytochrome P450, Steroid hormone biosynthesis, Vitamin digestion and
absorption, and Nitrogen metabolism.

**Categories (3.2).**

| Count | Mouse summary / shared / mixed / unresolved / human-led |
|---|---|
| Calls | 73 / 4 / 4 / 1 / 1 |
| Redundancy clusters | 67 / 4 / 4 / 1 / 1 |
| Calls, 73's mt-in-library classes | 68 / 6 / 7 / 1 / 1 |

- The calls still mostly summarise mouse zonation that the human sections do not share. The share is higher than
  before: 0.89 of resolved calls, against 0.76.
- 16 calls have a label-panel gene among their contributing members.

**Where the species differ (3.3).** Coordinate-free (b, primary):

| Peak segment | Human-high | Mouse-high |
|---|---|---|
| S1 | 14 | 17 |
| S2 | 4 | 4 |
| S3 | 18 | 26 |

- Direction: 36 human-high, 47 mouse-high (old 10 / 12). The coordinate peak (a) agrees on the direction for 82 of 83
  calls and on the segment for 64 (pooled transitions).
- The old pattern (human-high in S1, mouse-high in late S2) does not hold. Most calls differ most in S3, in both
  directions, and human-high calls split between S1 and S3.
- Two standing limits bear on this:
  - S3 compares human cortical S3 with mouse S3 that includes the outer stripe.
  - Human PT-S2 is small (632 units), so S2 is rarely the peak.
- Within-segment fractions are not read (C14). The human species mapping is undefined (section 4 of the reading).

**Support and robustness (3.4; no refit).**

| Check | Calls retained or called (of 83) |
|---|---|
| Segment F-test (75) | 71; 2 of the 12 misses are probe-sensitive |
| NB joint step test (75) | 32 |
| Steps on the spline (74) | 71 |
| Beyond label steps, test A (74) | 61 |
| Beyond coordinate steps, test B (74) | 68 |
| Count split, arms a / c (72) | 59 / 61 |
| All four winsor caps (71) | 79 |
| Detected in both species | 75 |
| Sex-unbiased | 60 |
| Flagged genes removed | 82 |
| Equal probe counts | 75 (8 probe-sensitive) |
| Proteasome subunits removed | 83 |
| **Core** (detected in both and sex-unbiased) | **60** |

Old: 20 of 22 coordinate-free, 8 of 22 spline steps, 22 of 22 beyond steps; sex-unbiased 21 / 20.

**External replication (3.5).** Matched random sets give p_random; the noise-call comparator has 26 sets (67
relabeling calls of raw-T_spatial cameraPR, 41 of which are among the 83).

| Reference, signed, matched | Mean z | p_random | p vs noise calls | Programs (label permutation) | Decision |
|---|---|---|---|---|---|
| Census − male mice | 1.27 | 0.0001 | 0.0008 | 0.0001 | specific |
| Lake − male mice | 1.08 | 0.0001 | 0.0038 | 0.0004 | specific |
| Census − female mice | 1.12 | 0.0001 | 0.021 | 0.0033 | specific |
| Lake − female mice | 0.95 | 0.0001 | 0.045 | 0.0094 | specific (core: better than random only) |

- **Selection-matched null (strata × z quintile):** not replicated in either atlas (p_random 0.35 and 0.64; old 0.0013
  and 0.0057). The calls replicate no better than random genes with the same positional strength. Replication is a
  property of strongly position-dependent genes, not of the call sets as such.
- **Slope-free:** better than random (0.0015, 0.0055), but not better than the noise calls (0.35, 0.63; old 0.26,
  0.13). Beyond the shared flattening, the calls show no set-specific replication.
- **Per call (BH ≤ 0.10 within the list):** 13 and 13 of 83 with male mice, 4 and 9 with female mice (old 12 and 5 of
  22; 5 and 3).
- **Genes:** the top 10% by z agree with the atlas interaction at ρ 0.78 and 0.78. That is above expression-matched
  random genes (median 0.56 / 0.58; 97.5th percentile 0.61 / 0.63) and now also above the most highly expressed genes
  (0.70 / 0.70; old 0.64 / 0.72).
- **Global slope:** ours −0.79, Census −0.85, Lake −0.88 (old −0.81, −0.86, −0.88).

**Conventional-only (3.6).**
- 199 conventional calls; 168 conventional-only sets tested.
- Mann–Whitney p 1.7 × 10⁻³ (Census) and 1.8 × 10⁻⁸ (Lake).
- Overlap-preserving: p 0.43 / 0.29 against random collections; D p 0.34 / 0.24.
- Decision: not supported once overlap is respected (unchanged).

**Lead genes (3.7).**

| Gene | z rank, top % (old T_spatial %) | Ours human S3 − S1, mt-free (old) | 73 class (old) | Flags |
|---|---|---|---|---|
| Gatm | 0.11 (0.14) | −0.74 (−0.78) | indeterminate (same) | label panel, axis anchor |
| Acadm | 0.33 (0.83) | −0.03 (−0.01) | indeterminate (conserved) | — |
| Acaa2 | 0.32 (1.04) | −0.02 (+0.04) | indeterminate (reversal) | — |
| Dcxr | 0.07 (1.05) | +0.86 (+0.97) | indeterminate (reversal) | — |
| Ugt3a1 | 0.05 (0.52) | −0.60 (−0.62) | reversal (same) | — |
| Rbp4 | 50.7 (0.08) | +2.41 (+2.46) | unclassified (same) | — |

- Every gene verdict and every story verdict is unchanged: creatine lead, fao lead, detox supporting, reversal lead.
  The attenuation is 0.854 (old 0.859).
- Gclc ranks in the top 1.5%. Cyp2e1 (a label-panel gene and axis anchor) is at 36.7%, Slc7a13 at 7.1%.
- Rbp4 is not expressed in mouse, so its positional F z is mid-ranked. R6 cites no rank for it.
- The full-library sensitivity moves the human values by ≤ 0.1 (Gatm −0.65, Dcxr +0.89).

**Programs and figures.**
- **Largest programs:**
  - Peroxisomal Lipid Metabolism (6 calls; mouse-high; S3; Figure 4e);
  - Glycosaminoglycan degradation (6; human-high; S1; the Figure 5d heatmap);
  - Bile Acid Metabolism (5; mouse-high; S3);
  - Transport Of Inorganic Cations/Anions And Amino Acids/Oligopeptides (4; human-high; S3);
  - Glutathione metabolism (4; human-high; S3).
- **Figure files:** `fig5d_programs_by_segment`, `fig5d_program_heatmap`, `fig5e_lead_genes`, and supplementary
  `supp_location_agreement_and_null` and `supp_lead_genes_coordinate`.
- **Literature export:** `calls_for_literature.csv` has the 83 calls with category, direction, segment, up to 8 lead
  genes with classes and flags, and old-list membership.

**How the protocol was read** (also in `tables/how_the_protocol_was_read.csv`).
1. The logic version is `77.package_species_calls.1`, by the coordinator's instruction; the protocol named
   `77.species_calls.1`.
2. **Human transitions.** The mean of the two sections crosses (HUK1_MED1, 72 C1), so the human segment-unit mapping
   is undefined. Where the protocol asks for each species' own segment units, pooled transitions are used for both
   species, and each species' label spans are marked. The mouse-transition sensitivity is reported (65 of 83 agree).
3. **Categories.** They use contributing members without notebook 75's fewer-than-3 fallback, as 43 §7 did.
   Contributing genes are 75's `nlargest` top 10% of z.
4. **Program categories** take the most common category of their calls.
5. **Gene-removal checks** recompute cameraPR's FDR over the subset's tested sets without a squeezeVar refit. The
   effect is the members' mean z minus the mean of every other subset gene.
6. **Winsorizing:** a call is kept when no cap lists it as lost.
7. **Per-call replication** uses the exact-moment joint z from the collection run, not 37's Monte Carlo test.
8. **Gene-level comparison scores:** the top 10% by mean expression and by NB T_level (71's cache).
9. **Noise-call programs** use the call programs' rule.
10. **Attenuation controls** need external contrasts from 38's old-universe table. Our contrasts use the mt-free
    library.
11. **The ours arm** uses every export gene on 71's units.
12. **The literature export** flags membership of the old 60/36 lists, which the 2026-10-04 check covered.
13. **Figure 5d heatmap program.** Figure 4e's program is Peroxisomal Lipid Metabolism (mouse-high), so the heatmap
    shows Glycosaminoglycan degradation.
14. **M2** reads the median of the segment index; a median between segments counts as "not S1".
15. **M6** reads ranks as percentiles of 71's z.
16. **Seeds:** 50 for the nulls, 12 for the gene draws, 38 for the attenuation bootstrap.
