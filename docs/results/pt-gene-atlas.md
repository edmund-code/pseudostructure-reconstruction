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

## Notebook 112 results

`analysis/notebooks/112_pt_gene_atlas.ipynb` (logic `112.gene_atlas.1`) writes to `results/pt_gene_atlas/`. It
followed the protocol frozen in 14e685d with amendment 1 (d4c863f) and ran in 12 s. 10,361 genes on the AKI axis
and 11,055 on the species axis.

**Q5 · Recovery.** All nine checks reproduce the numbers already committed in the source documents: 3,159 level and
362 location calls (91); 126 AKI calls and 171 species gradients failing the core rebuild (104); 84 suspect genes
(109); 66 confirmed handoff pairs (101); 186 AR-maintained S3 genes (111); 43 injury-stable positional genes (105);
237 eligible transporters.

**Line sizes.** AKI: level 3,159, location 362, F-test 4,244, position-split 2,719, Kirita 12 h 3,899, KPMP 2,461,
phase 1,010, variability 1,217, androgen 195. Species: level 5,775, F-test 3,889, zonation class 614, atlas level
4,905, rewired 25, reversed gradient 74, handoff 58.

### Q1 · Evidence converges, and the more it converges the more it is known biology

| k lines | AKI observed | null mean | enrichment | species observed | null mean | enrichment |
|---|---|---|---|---|---|---|
| 1 | 7,909 | 9,029 | **0.88** | 8,614 | 9,248 | **0.93** |
| 2 | 5,312 | 5,970 | 0.89 | 4,936 | 4,833 | 1.02 |
| 3 | 3,158 | 2,880 | 1.10 | 1,522 | 1,165 | 1.31 |
| 4 | 1,712 | 1,036 | 1.65 | 238 | 90 | 2.65 |
| 5 | 795 | 285 | **2.79** | 24 | 3.3 | **7.36** |
| 6 | 295 | 58 | **5.08** | 5 | 0.09 | **58** |

- **Q1 holds** (k = 3 beyond the 99th percentile in both atlases), but the interesting part is the shape. At k = 1
  and 2 the observed count is *below* the null: the lines concentrate on fewer genes than chance would spread them
  over. The excess only becomes large at k ≥ 4.
- **C1, the relabeled control, collapses.** Rebuilding the level, location, variability and androgen lines with a
  balanced relabeling in the condition's place drops the k = 5 count from 795 to 150 and 162, and k = 6 from 295 to
  11 and 9. The reconstruction is exact: it reproduces notebook 91's own published relabeled counts (level 40 and 3
  against 3,159; location 0 and 1 against 362), so the control is that notebook's I1 check carried through the join.
- **Sensitivity.** Keeping one line per similarity group (8 of 9 on AKI, 5 of 7 on species) leaves k = 3 beyond the
  99th percentile in both.

**The validity curve (post hoc).** The annotation sets were never used to build a line, so the share of a core list
that falls in a named literature, disease, transporter or androgen set is an independent check of the convergence
score. It rises monotonically with k:

| k | AKI anchored share | species anchored share |
|---|---|---|
| background | 0.209 | 0.154 |
| ≥ 2 | 0.318 | 0.200 |
| ≥ 3 | 0.401 | 0.340 |
| ≥ 4 | 0.490 | 0.576 |
| ≥ 5 | **0.622** | **0.870** |
| ≥ 6 | **0.775** | 0.800 |

Counting independent lines therefore does measure something real: a gene carried by five lines is three to six times
more likely to be a gene the field has already named than a gene drawn at random from the universe.

### Q2 · The lines are close to independent

- AKI: 8 groups among 9 lines at phi 0.5. The only pair above 0.5 is the positional F-test and the position-split
  selection (0.54), which are two readings of the same positional statistic. Everything else sits at 0.02–0.30; the
  level screen and Kirita's 12 h level agree most among the rest (0.30).
- Species: 5 groups among 7. The F-test pairs with the zonation-class disagreement (0.58) and rewiring with the
  reversed gradients (0.57). The level screen and the atlas level agree at only 0.28, and the handoff line is
  *negatively* related to reversed gradients (−0.22), as a confirmed handoff requires a conserved order.
- So the convergence count is not a count of restatements. The two positional readings should be treated as one.

### Q3 · Which screen yields the most interpretable biology

| Atlas | Method | called | clean share | supported share | anchored share | unique share | **yield** | **novel yield** |
|---|---|---|---|---|---|---|---|---|
| AKI | level (91) | 3,159 | 0.902 | **0.708** | 0.372 | 0.380 | 2,017 | 1,266 |
| AKI | location (91) | 362 | 0.843 | **0.711** | **0.724** | 0.019 | 217 | 60 |
| AKI | positional F-test | 4,244 | **0.934** | 0.586 | 0.367 | 0.329 | **2,323** | **1,470** |
| AKI | position-split | 2,719 | 0.917 | 0.611 | 0.417 | 0.093 | 1,523 | 888 |
| species | level | 5,775 | 0.911 | **0.568** | 0.192 | **0.632** | **2,991** | 2,417 |
| species | positional F-test | 3,889 | 0.913 | 0.466 | 0.312 | 0.403 | 1,655 | 1,139 |
| species | zonation class | 614 | 0.904 | 0.470 | **0.452** | 0.049 | 261 | 143 |

- **On volume the positional F-test wins the AKI axis** (2,323 supported genes, 1,470 of them unanchored) and the
  level screen wins the species axis (2,991).
- **On precision the level and location screens win.** 71% of their clean calls are supported externally against
  59% for the F-test. The ordering is the same at external \|t\| ≥ 1 and ≥ 3 (sensitivity table), so it is not a
  threshold effect.
- **The location screen is the most literature-recognised** (72% of its supported calls are in a named set, against
  37% for the level and F-test screens) and the least redundant in the other direction: only 2% of its calls are
  unique, because nearly every location call is also a level or shape call.
- **The species level screen is the most private**: 63% of its calls are made by no other species screen, and only
  19% land in a named set — the lowest anchoring of any method here.
- **Contamination costs little.** The clean share is 0.84–0.93 everywhere, so the 104 and 109 filters remove under a
  sixth of any list. They are cheap to run and should be run, but they do not separate the methods.

### Q4 · The core list, and the prediction that failed

- **The frozen prediction was wrong.** I predicted under 40 unanchored genes in the AKI core list. At the frozen
  k ≥ 3 there are 2,450 core genes, 1,468 of them unanchored; on the species axis 1,230 and 812. The reason is
  visible in the Q1 table: at k = 3 the convergence enrichment is only 1.10×, so the frozen bar sits where
  convergence is barely beyond chance.
- **Post-hoc strict core (k ≥ 5, chosen after seeing the enrichment column).** AKI 699 genes, 264 unanchored;
  species 23 genes, 3 unanchored.
- **The AKI strict core is one coherent program:** peroxisomal and mitochondrial fatty-acid and amino-acid
  catabolism falling together with S2/S3 identity — Csad, Pank1, Fmo5, Amacr, Dao, Crot, Agps, Scp2, Acadm, Abat,
  Hgd, Haao, Kyat3, Mccc1, Acy3, Glyat, Aldh1l1, Pipox, Cth, Grhpr, Aspa — against a small induced set (Anxa2,
  Anxa3, Krt8, Tgm2, Zwint, Stmn1, Cdk2).
- **The species strict core is a reversal list:** Slc16a9, Lap3, Mme, Hpgd, Slc17a3, Slc5a12, Slc27a2, Vnn1,
  Acot12, Gclm, Napsa, Slc9a8, Slc38a10, Slc6a6, Abcc4, Pdk4, Sfxn4, Spats2l, Pitpnc1, Larp4b, Ninl, Sulf2,
  Tmem243 — 12 of the 23 carry a mouse gradient that the human does not share (zonation class "reversal" or
  "mouse-only"), and 20 of 23 are in a named set.
- **Unanchored genes worth a look** (AKI, ≥ 5 lines, in no named set): Rorc, Dhrs3, Aldh8a1, Susd3, Chpt1, Aamdc,
  Esrrg, Sirt3, Cdo1, Acss2, Pipox, Echdc2, Grtp1, Thnsl2, Aspa, Cth, Grhpr, Bphl. Several are the expected
  company of the catabolic program (Sirt3, Esrrg and Rorc are its regulators; Pipox, Aspa, Cth and Grhpr are
  peroxisomal or sulfur-amino-acid enzymes), which is a hint that "unanchored" here means under-curated rather
  than novel.

**Where the best-supported genes change (post hoc, from `aki_atlas.csv.gz`).** The share of each injury location
among all genes carrying one, and among the core lists:

| injury location | background | k ≥ 3 | k ≥ 5 | k ≥ 6 | k ≥ 7 | k ≥ 7 vs background |
|---|---|---|---|---|---|---|
| inside S1 | 30.4% | 23.1% | 18.6% | 17.0% | 12.5% | **0.41×** |
| at S1/S2 | 24.5% | 21.4% | 14.8% | 14.2% | 6.2% | **0.25×** |
| inside S2 | 26.0% | 28.5% | 29.3% | 27.0% | 15.6% | 0.60× |
| at S2/S3 | 10.8% | 14.4% | 18.6% | 18.4% | 37.5% | **3.48×** |
| inside S3 | 8.2% | 12.6% | 18.8% | 23.4% | 28.1% | **3.41×** |

The convergence ranking was built from gene-level calls that carry no position of their own, so this is an
independent recovery of the late-PT weighting of injury.

### What this answers

1. **The convergence count is a usable quality score, and it should be read at k ≥ 5, not k ≥ 3.** It survives the
   expression-stratified null, collapses under balanced relabeling, and is validated by an annotation set it never
   saw.
2. **There is no single best screen; there are two.** The positional F-test finds the most, the level and location
   screens are right most often, and they overlap at phi 0.17–0.28. Running both and intersecting is what the atlas
   actually rewards — which is what the scoreboard already recommends, now with a gene count attached.
3. **Convergence recovers the late-PT weighting of injury without a new fit.** Among all genes with a location,
   10.8% are at the S2/S3 boundary and 8.2% inside S3; among the k ≥ 5 core these rise to 18.6% and 18.8%, and
   among k ≥ 7 to 37.5% and 28.1% (3.5× and 3.4× the background), while "inside S1" and "at S1/S2" fall to 0.25–0.41×.
   The ranking that the atlas produces agrees with notebook 99's acute S3 weighting and notebook 105's external
   sections, by a route that shares none of their statistics.
4. **The 104 and 109 filters are cheap and nearly method-neutral.** They remove 7–16% of every list and do not
   change the ranking of the methods.
5. **For the paper:** the k ≥ 5 strict core is the gene-level table to carry, with the line-by-line calls beside
   each gene; `aki_atlas.csv.gz` and `species_atlas.csv.gz` are the supplementary tables that make every other
   notebook's gene list reproducible from one file.
