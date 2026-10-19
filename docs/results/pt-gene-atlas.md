# One row per gene: joining every per-gene result, and which method yields the most biology (notebook 112)

**The question.** Twenty-odd notebooks have each produced their own per-gene table (108 gene-keyed tables under
`results/*/tables/`). No analysis has ever put them side by side. Two things follow from doing so, and neither has
been measured:

1. **Convergence.** Which genes are called by several *independent* lines of evidence, and is that more than the
   overlap expected from each line's size and the shared expression advantage? A gene carried by four independent
   measurements is a paper claim; a gene carried by one is a hypothesis.
2. **Method yield.** The scoreboard in `docs/results/pt-discovery-methods.md` scores methods on calibration and
   external replication. It has never scored them on the thing the paper actually needs: **how many clean,
   externally supported, interpretable genes each one delivers**, and how many of those no other method found.

This notebook is a synthesis, not a new fit. It computes no new model on the expression data. Every statistic it
uses was computed and committed by an earlier notebook. What is new is the join, the convergence null, and the
yield criteria — all fixed below before any of them is evaluated.

**Why this is worth doing now.** The method line has reached the point where the remaining question is not "does
method X calibrate?" but "which combination of methods should the paper run, and what does it find that the others
miss?". That is a question about the union and the intersection of the gene lists, which is exactly what an atlas
answers.

Frozen on 2026-10-19, before any join was executed. Before freezing I listed the gene-keyed tables and their
columns (shapes and column names only; no values were read other than the first two rows of each table, which were
inspected to confirm the key column).

## 1 · Inputs and keys

- Notebook `analysis/notebooks/112_pt_gene_atlas.ipynb`; logic `112.gene_atlas.1`; seed 112.
- Results `results/pt_gene_atlas/`.
- New helper module `pseudospace/gene_atlas.py` with `join_lines`, `convergence_counts`, `stratified_line_null`
  and `yield_profile`, each tested on synthetic fixtures.
- **Two atlases, two keys.**
  - **AKI atlas**, keyed on mouse gene symbol. Universe: the genes in notebook 91's shared AKI universe
    (`results/pt_positional_decomposition/tables/i_genes.csv.gz`).
  - **Species atlas**, keyed on the shared-space gene symbol. Universe: the eligible genes of
    `results/pt_pathways_package_pseudospace/tables/universe_human_vs_mouse.csv`.
- A gene is **testable on a line** when it is in that line's own universe. Lines are only compared on the genes
  they could both have called; the convergence count is a count over testable lines, and the null preserves this.

## 2 · The lines of evidence (fixed now)

Each line is one binary call per gene, with a sign where the source defines one. The source table and column are
fixed here so the join cannot be tuned.

**AKI lines (mouse).**

| Line | Meaning | Source | Call |
|---|---|---|---|
| A1 | whole-PT level change | 91 `i_genes` | `level change` ≠ "none" |
| A2 | location change beyond level | 91 `i_genes` | `location change` true |
| A3 | positional F-test (shape) | `pt_normalized_t_spatial/gene_change_aki_vs_control` | `q` < 0.05 and `converged` |
| A4 | position-split selection | `pt_pathway_networks/row_n_genes_aki_vs_control` | `selected (row N)` |
| A5 | external mouse IRI (Kirita 12 h) | `pt_external_iri/external_genes` | same sign as our level and \|t\| ≥ 2 |
| A6 | external human AKI (KPMP) | `pt_human_aki_benchmark/external_genes` | same sign as our level and \|t\| ≥ 2 |
| A7 | injury phase membership | `pt_injury_phases/phase_programs` | in the acute or sustained program |
| A8 | tubule-to-tubule variability | 107 `genes` | `beyond both relabelings` |
| A9 | androgen instrument | 111 `sets` | AR-maintained (any segment) and `L` < 0 |

**Species lines (human vs mouse).**

| Line | Meaning | Source | Call |
|---|---|---|---|
| S1 | level difference | `pt_species_level_control/gene_level` | `direction` ≠ "none" |
| S2 | positional F-test (shape) | `pt_normalized_t_spatial/gene_change_human_vs_mouse` | `q` < 0.05 and `converged` |
| S3 | zonation-class disagreement | `pt_zonation_classes/gene_zonation_classes` | mouse and human calls differ on any contrast |
| S4 | atlas replication of the level call | `pt_species_level_control/gene_level` | `atlas_strong` |
| S5 | co-expression rewiring | 97 `genes` | `category` names a rewired class |
| S6 | gradient reversal | 106 `genes` | `category` = reversed |
| S7 | confirmed transporter handoff | 101 `pairs` | member of a `confirmed` pair |

**Filters (not lines).** These remove calls rather than support them, and are applied after the counts so their
effect is measurable.

- **F1 core rebuild** (104): `survives core3` for the matching kind (AKI call or species gradient).
- **F2 neighbour marker** (109): not in the post-hoc suspect set.
- **F3 edge enrichment** (104): the gene's edge enrichment below the 95th percentile of expression-matched genes in
  the relevant specimens.
- **F4 probe** (species only): `probe_clean_both`.

**Annotation (not lines, not filters).** Named membership used for the "anchored" criterion: the literature
programs of notebooks 79–82 (`pt_known_injury_programs`, `pt_known_species_programs` members), the transporter
family list (101 eligible), the disease sets of 108, the AR sets of 111, the injury-stable positional genes of 105,
and the hand-named genes of 108's `h4_named`.

## 3 · Tests (fixed now)

**Q1 · Does evidence converge beyond chance?**
- Statistic: the number of genes called by ≥ k testable lines, k = 2…6, separately for the AKI and species atlases.
- **Null:** within each line, permute the call labels among testable genes **inside expression deciles** (deciles
  from the control mean of `universe_*`), keeping the number of calls per decile. 2,000 draws. This destroys the
  gene-specific agreement while keeping each line's size and its expression advantage.
- **Holds** when the observed count at k = 3 exceeds the 99th percentile of the null in both atlases.
- **Relabeled control (C1):** rebuild A1, A2, A8 and A9 from their stored R1 and R2 relabeled columns, leaving the
  other lines as they are, and recompute. The excess at k = 3 should shrink towards the null.

**Q2 · How independent are the lines?**
- Statistic: the pairwise Jaccard and the tetrachoric-free phi coefficient between every pair of lines, on the
  genes testable on both; then a correlation-matrix clustering (average linkage) into line groups.
- **Reported, not passed or failed.** The effective number of independent lines is the number of groups at phi 0.5.
- The convergence count is recomputed with one line per group (the line with the largest universe) as a
  sensitivity.

**Q3 · Which method yields the most interpretable biology?**
- For each **calling** method (A1, A2, A3, A4 on AKI; S1, S2, S3 on species — the lines that are our own screens
  rather than external or instrumental evidence), on its own called gene list:
  - **clean share**: passes F1, F2, F3 (and F4 for species);
  - **supported share**: among clean genes, the share with external support (A5 or A6 for AKI; S4 for species);
  - **anchored share**: among clean and supported genes, the share in a named annotation set;
  - **unique share**: the share called by no other calling method on the same axis;
  - **yield**: the number of genes that are clean **and** supported;
  - **novel yield**: clean and supported and **not** anchored and called by ≥ 2 lines overall.
- **Reported as a profile.** The ranking by yield is the answer to the question; no threshold passes or fails.
- **Sensitivity:** yield recomputed with the external threshold at \|t\| ≥ 1 and ≥ 3.

**Q4 · The gene list the paper should carry.**
- The **core list**: genes clean, supported and called by ≥ 3 lines. Reported with every line's call and sign.
- The **novel tail**: core-list genes in no named annotation set. These are read against the literature after the
  run and labelled as a post-hoc reading.
- **Prediction fixed now:** the core list will be dominated by S3 identity and transport genes already named in
  notebooks 91, 94 and 99, and the novel tail will be small (I predict under 40 genes on the AKI axis). A large
  novel tail would mean the earlier notebooks under-reported; an empty one would mean the atlas adds nothing beyond
  what the single best screen already gave.

**Q5 · Does the atlas reproduce its sources?**
- Recovery checks against numbers already committed in the result documents: 91's level and location call counts,
  104's survival shares, 109's suspect counts, 101's confirmed pair count, 111's AR-maintained set sizes.
- Any mismatch is a join error and must be fixed before Q1–Q4 are read.

## 4 · Outputs and figures

- **Tables:** `aki_atlas.csv.gz`, `species_atlas.csv.gz` (one row per gene, every line, filter and annotation);
  `convergence.csv`, `line_similarity.csv`, `method_yield.csv`, `core_list.csv`, `novel_tail.csv`,
  `recovery_checks.csv`, `how_the_protocol_was_read.csv`, `run_record.json`.
- **Figures (plots only):**
  - 112a: observed genes at ≥ k lines against the null band, both atlases;
  - 112b: the line-similarity matrix;
  - 112c: method yield — clean, supported, anchored, unique — as stacked bars per method;
  - 112d: the core list as a line-by-gene call matrix for the top genes by convergence.

## 5 · Caveats

- The lines are not equally powered, and three of them (A5, A6, S4) are external replication, which is evidence of a
  different kind from a screen call. They are kept as lines for Q1 and excluded from the methods compared in Q3.
- Convergence rewards genes that are easy to measure. The expression-decile null removes the marginal advantage but
  not the fact that an abundant gene is called more confidently by every method.
- A gene absent from the core list is not a negative result; it may be testable on only one line.
- The atlas inherits every caveat of its sources, in particular the 2-vs-2 design.

## Protocol amendment 1 (2026-10-19, before any join was run)

Loading the source tables showed that three of the frozen rules name a column that does not hold a call. Only column
names, dtypes and the value sets of categorical columns were read; no statistic, count or gene list was evaluated.

- **S1 (species level).** `direction` in `pt_species_level_control/gene_level` is a numeric statistic, not a
  category. The operational level call in notebook 68 is **\|L\| ≥ 1**, the threshold of its criterion (ii). S1
  becomes \|L\| ≥ 1 with `converged` and not `separated`.
- **S4 (atlas replication).** Read as notebook 68's criterion (ii) per gene: `atlas_strong` **and** the atlas
  coefficient `A` sharing the sign of `L`.
- **S5 (co-expression rewiring).** `category` in notebook 97's `genes.csv` is the gradient category, not a rewiring
  class. S5 becomes membership of notebook 97's doubly rewired set: `post_hoc_gene_lists.csv` with
  `group == "rewired in both"`.
- **Readings (not changes).** A1 and A2 are the boolean columns `level change` and `location change`. A7 is a
  `phase_programs` program beginning with ACUTE or SUSTAINED. F3 flags a gene whose mean edge enrichment over the
  relevant specimens is above the 95th percentile **of its own expression decile**.
