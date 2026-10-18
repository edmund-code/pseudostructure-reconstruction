# Pathway networks along the PT: scoranet's ideas against the positional F-test and clustering (notebook 93)

## Notebook 93 protocol

Frozen on 2026-10-18, before any fit of notebook 93. Nothing below has been computed on our data.

**Why.** On 2026-10-18 the user asked whether scoranet's single-cell over-representation networks could help interpret
our results. The follow-up request was to borrow its ideas, add what it cannot do, and test on our data whether it is
a better way than the positional F-test or the clustering method. The trigger was notebook 91's injury screen: 538
level sets and 236 location sets. 71's redundancy step (whole-set Jaccard ≥ 0.5) still leaves 261 and 121 groups, too
many to read.

**scoranet**, as read from its code (GitHub mikemartinez99/scoranet, v1.0.0, released 2026-10-17; R/cluster_de.R,
cluster_ora.R, network_utils.R, network_themes.R, DIAGNOSTICS.md). For each cell type it:
1. runs Seurat `FindMarkers` (Wilcoxon, cells as replicates) and keeps genes with padj ≤ 0.05 and |avg_log2FC| > 0.5;
2. runs GO over-representation (clusterProfiler `enrichGO`, BH), optionally on the tested genes as the universe and
   split into up and down genes;
3. builds a network of the enriched terms (p.adjust < 0.05). Each term is represented by its DE genes, an edge joins
   two terms whose DE genes have Jaccard ≥ 0.1, and terms without an edge are dropped;
4. groups the terms with Leiden (igraph `cluster_leiden`, modularity, resolution 1, edge weights, iterated until stable,
   seed 1234);
5. gives each term a direction score z = (n_up − n_down)/√n over its DE genes: up, down or mixed at |z| ≥ 1;
6. names each group by up to two GO ancestors, picked greedily. Each pick maximises (terms it newly covers) ×
   information content log(ontology size / (descendants + 1)), must newly cover ≥ 15% of the group's terms, and may
   have at most 1,000 descendants. If no ancestor qualifies it falls back to keywords or the hub term;
7. plots diagnostics: overlap threshold, resolution sweep, and driver genes (the genes in most terms of a group).

The idea is the same as an enrichment map (Cytoscape EnrichmentMap; clusterProfiler `emapplot`). Those papers were
not read in full for this protocol.

**What we take, change and add.**
- **Not taken: its DE step.** It treats units as replicates, which with 2 vs 2 specimens produces false hits. Every
  gene statistic here is 71's NB gene model, read and never refit.
- **Taken as a calling row (row N):** the gene-first ORA on the tested universe, with scoranet's up/down split extended
  to up/down × segment (A1 below).
- **Taken as an organising layer, for every method's calls:** the driver-gene network, Leiden, ancestor naming and
  driver genes.
- **Changed:** terms without an edge are kept as one-call programs, so every call sits in exactly one program.
- **Added, because scoranet has no axis:**
  - A1, a position split of the gene lists;
  - A2, positional node attributes (profile, peak, trough, a position score);
  - A3, program profiles and their coherence, against a null;
  - A4, replication of each program's profile across the two case specimens;
  - A5, a second network that links calls by where along the PT their difference lies, not by shared genes (NET-pos);
  - A6, specificity at the level of programs (relabeled programs);
  - A7, program stability across reconstructions (notebook 69 → 71).

**Questions.**
- **Q1 · Calling.** Is position-split ORA (row N) specific, stable, unexposed to probe differences and robust to hub
  genes, compared with the positional F-test (F) and clustering (C)?
- **Q2 · Organising.** On each method's calls, does the network group them into fewer programs that are more coherent
  along the PT, replicate across specimens, survive relabeling and survive a change of reconstruction better than the
  current grouping (75's overlap programs, OVL) and 71's redundancy clusters (RED)?
- **Q3 · Discovery.** What programs do the methods find, where along the PT, and do the methods find different programs
  or the same programs through different sets?

Everything is descriptive: two control mice, two AKI mice, and two cortex sections of one human donor. Specimens are
the replicates.

### 0 · What was seen before freezing

- scoranet's code, documentation and vignette (read in full).
- **Hierarchy coverage of 71's sets** (human vs mouse, a count of ids and names only):
  - 1,169 of the 1,181 Reactome sets have their R-HSA id in the current human release;
  - 262 of the 282 KEGG sets match a br08901 map by exact normalised name.
- limma is reachable through rpy2 in the analysis environment.
- The structure of 71's caches: gene models with specimen curves, and clustering runs with their gene groups and the
  ORA table including the overlap genes.
- 71's unit-matching table holds notebook 69's positions and labels for the units in both runs.

No gene statistic, call, program or network of ours was computed.

### 1 · Inputs (read, never refit)

**Notebook 71** (`results/pt_pathways_package_pseudospace/`) for both datasets (human vs mouse, AKI vs control):
- the universe and gene sets;
- `gene_statistics_{ds}.csv` (z, q, `converged`, `separated` per partition);
- the cameraPR tables per partition;
- `redundancy_clusters_{ds}.csv`;
- the stage-cached gene models (Δ(x) on the 101-point grid and the four specimen curves);
- the clustering runs (`method_b_*`: gene groups, clusters and ORA, the condition run with the default null and both
  relabelings with the fair null);
- the S1 tables (offset lean);
- the unit-matching tables.

Inputs are rebuilt as 75 builds them (`package_inputs`, mt genes out) for the units, positions and segment labels, the
package gene-set object, and the probe exposure.

**Notebook 69** (`results/pt_pathways_without_mito/`), for A7 and the reconstruction check:
- the gene models, gene statistics and cameraPR condition tables;
- the default-null clustering run;
- its redundancy clusters;
- its gene sets, rebuilt by 71's code (`pt_inputs` on SCF13 for human vs mouse; notebook 12's universe for AKI).

**Notebook 75** (`results/pt_package_pathway_methods/tables/`): `programs_{ds}.csv` (guard), `method_table.csv`
(reference rows) and `probe_exposure.csv`.

**Notebook 91** (`results/pt_positional_decomposition/tables/`): `i5_library_sets.csv` and `i_genes.csv.gz`, for the
application in section 7.

**External** (`data/external/pathway_hierarchies/`, downloaded 2026-10-18):

| File | Source | sha256 |
|---|---|---|
| `ReactomePathwaysRelation.txt` | reactome.org/download/current, release 97 | `fd49a624…3f1e` |
| `ReactomePathways.txt` | same | `f6d7a2bf…df3a` |
| `br08901.keg` | KEGG REST `get/br:br08901` | `86f3c003…a406` |

`probe_counts.csv` (notebook 38) gives the equal-probe genes.

**Guards (not cached):**
- G1: the rebuilt universe, sets and units equal 71's, and the cached gene models reproduce 71's T_spatial and z
  (75's G1–G3).
- G2: F's and C's calls and relabeled calls reproduce 71's method table.
- G3: OVL on F's calls reproduces 75's `programs_{ds}.csv`.
- G4: 69's sets equal its saved tables (71's guard), and F's and C's 69 → 71 Jaccard values reproduce 71's
  `compare69_pathways_summary.csv`.
- G5: the clustering ORA recomputed from a run's gene groups reproduces that run's calls (71's G7).
- G6: the hierarchy files match the sha256 above.
- G7: Leiden run twice with the same seed gives the same membership.

### 2 · Definitions

- **Grid and segments.** Positions are 71's 101-point grid. The segment cuts c12 and c23 are computed on the reference
  specimens' units. Each cut is the position that minimises misassigned units: S1 against S2 + S3 for c12, and S1 + S2
  against S3 for c23. Ties go to the midpoint of the tied range. Grid points are S1 below c12, S3 at or above c23, and
  S2 between.
  - For notebook 69 (A7), the cuts use 69's positions and labels of the reference units present in both runs (71's
    unit-matching table).
- **Usable gene** in a partition: converged, not separated, finite Δ(x).
- **Gene shape:** δ_g(x) = (Δ_g(x) − mean over the grid of Δ_g) / ln 2, the positional part of the case − reference
  difference in log2.
  - Half-range a_g = (max δ_g − min δ_g)/2.
  - **Extreme:** "up" at the argmax when max δ_g ≥ −min δ_g, otherwise "down" at the argmin. The extreme's segment is
    the segment of that grid point.
  - **Late or early:** a gene is *late* when it is up in S3 or down in S1 (the case is relatively higher late), and
    *early* when it is up in S1 or down in S3.
- **Driver genes** of a call, by method (relabelings use their own run's statistics):
  - **F** (71's cameraPR: Up and FDR ≤ 0.05): usable members with gene q ≤ 0.05. If fewer than 3, usable members in
    the top 10% of z (75's contributing genes); if still fewer than 3, all usable members.
  - **C** (71's clustering, default null; fair null in the relabelings): the union of the overlap genes over the lists
    where the set has matched q < 0.05.
  - **N** (A1): the union of its members in the lists where it has q ≤ 0.05.

### 3 · Q1: the calling row (row N) against F and C

**Row N (A1).**
- Selected genes per partition: usable, gene q ≤ 0.05 (that partition's positional F-test q), and a_g ≥ 0.5 log2
  (scoranet's |log2FC| > 0.5).
- Six lists: up or down × the extreme's segment (S1, S2, S3).
- ORA with `discrete_comparators.family_ora` over the six lists:
  - background = the partition's usable genes (scoranet's tested universe);
  - 71's sets as members over that background;
  - BH over set × list; call = any list q ≤ 0.05;
  - calibration p = min(1, 6 × best p).
- Tested sets: those with 10–300 members in the background.

**Criteria**, computed the same way for F, C and N on each dataset:
- **V1 · swap rule:** 75's rule (calls > 0, NCDR < 0.25, each relabeling below 5% of tested sets, calibration share
  ≤ 0.075).
  - F and C are read from 71's runs; their calibration p follow 75's conventions.
  - In a relabeling with no selected genes, N calls nothing and has calibration p 1.
- **V2 · reconstruction stability:** the Jaccard of the condition calls in notebooks 69 and 71, over the sets tested in
  both. N is rerun on 69's gene models, gene q and sets with 69's cuts.
- **V3 · probe exposure** (human vs mouse only): 1 − |calls ∩ equal-probe calls| / calls, following 75's section 4.
  - The equal-probe genes are `probe_balance == 'balanced'`.
  - N is rerun with lists and background restricted to them; F's and C's values are read from 75's
    `probe_exposure.csv`.
- **V4 · hub-gene robustness:** the share of a method's calls kept after removing its 10 hub genes, i.e. the driver
  genes of the most calls (ties broken by higher condition z).
  - F: cameraPR on the same z without them, sets below 10 members dropped (71's proteasome convention).
  - C: the package ORA (`ora_by_group`) rerun on the condition run's groups and clusters, with them outside the
    universe (71's convention).
  - N: its ORA rerun without them in the lists and background.
- **Descriptive:**
  - calls and calls in 71's redundancy clusters;
  - the offset lean of calls (75's section 5 rule);
  - the overlap of N with F and with C (Jaccard over sets and over redundancy clusters).

**Verdict per dataset, N against F and N against C** (dominance):
- A criterion is equal when the difference is within its tolerance: Jaccard and shares 0.05; V1 pass equals pass.
- N is **better** when it is at least as good on every criterion defined on that dataset and strictly better on one.
- N is **worse** in the reverse case. Otherwise it is **mixed**, and the profile is listed.
- Direction: V1 pass > fail; V2 and V4 higher is better; V3 lower is better.

### 4 · The network and the positional additions

**NET (scoranet's network; isolates kept):**
- Nodes are a method's calls; the edge weight is the Jaccard of their driver genes; an edge exists where it is ≥ 0.1.
- Leiden: python-igraph `community_leiden`, modularity, edge weights, resolution 1, β 0.01, n_iterations −1, with
  Python's `random` seeded at 1234 before each call.
  - This is the same igraph C core as R's `cluster_leiden`, but the random streams differ from scoranet's.
- Isolates are one-call programs.
- Programs are numbered by number of calls, ties by the smallest best p.

**A2 · Call attributes:**
- the call profile, the mean δ_g over its usable driver genes;
- the peak (grid position of the profile maximum: where the case is relatively highest) and the trough;
- the profile half-range;
- the position score (n_late − n_early)/√n over its usable driver genes: toward S3 at ≥ 1, toward S1 at ≤ −1,
  otherwise mixed (scoranet's threshold).

**A3 · Program profile and coherence:**
- The program profile is the mean δ_g over the union of its calls' usable driver genes.
- Coherence R² = 1 − Σ_g Σ_x (δ_g − profile)² / Σ_g Σ_x δ_g², over that union.

**A4 · Replication:**
- For each case specimen k: δ_g,k(x) = (curve_k(x) − mean of the two reference specimens' curves(x)), centred over the
  grid, in log2.
- Only genes converged and not separated in all four specimen fits are used.
- Program profile_k is the mean over the program's driver genes; r = Pearson(profile_1, profile_2) over the grid.
- Human vs mouse: the two case specimens are two sections of one donor, so r is section to section. AKI: two mice.

**A5 · NET-pos (a network scoranet cannot build):**
- An edge joins two calls whose call profiles correlate at Pearson r ≥ 0.8, weighted by r.
- Same Leiden settings; isolates are one-call programs.
- It groups calls by where along the PT their difference lies, whatever genes they share.

**Driver-gene panel:** per program, the genes that are driver genes of the most of its calls (scoranet's), with their
δ_g along the PT.

**Naming** (scoranet's ancestor rule):
- **Reactome:** the R-HSA id at the end of a set's name, and ancestors through `ReactomePathwaysRelation.txt` (human
  rows, transitive).
- **KEGG:** match the map name in `br08901` (exact normalised; then without parenthetical text and without "- animal").
  Its ancestors are its B (subcategory) and A (category) entries.
- **Hallmark:** no ancestors.
- **Candidates:**
  - each call covers itself and its ancestors;
  - a candidate may have at most 4% of its hierarchy's pathways as descendants (scoranet's 1,000 of about 26,000 GO
    BP terms);
  - IC = log(hierarchy size / (descendants + 1)).
- **Greedy picks:** up to 2. Each maximises (calls newly covered) × IC and must newly cover ≥ 15% of the program's
  calls. Ties go to the lower id.
- **Name:** the picks joined by " · ". If none qualifies, the hub call (highest degree, then smallest p).
- Unnamed calls are counted.

### 5 · Q2: groupings on the same calls

**Groupings,** for each method's condition calls (F, C, N) on each dataset:
- NET;
- NET-pos;
- OVL: `pathway_decoys.overlap_programs` with 75's rule (contributing genes = top 10% condition z, cut 0.5);
- RED: 71's redundancy clusters restricted to the calls.

**Criteria:**
- **W1 · programs per call** (lower = more compressed).
- **W2 · coherence z:**
  - the call-weighted mean program R² (A3);
  - against 1,000 random partitions of the same calls with the grouping's own program sizes (seed 93);
  - z = (observed − null mean) / null SD; undefined when every program is a single call.
  - Not judged for NET-pos, which groups by profile similarity by construction (circular); it is reported.
- **W3 · replication z:** the call-weighted mean r (A4), against the same random partitions.
- **W4 · reconstruction stability:**
  - the adjusted Rand index between the grouping built on 69's calls and statistics and the one built on 71's, over
    the calls in both runs;
  - RED uses 69's and 71's own clusters.
- **W5 · program NCDR:**
  - the mean number of programs over the two relabelings divided by the condition programs;
  - every grouping is rebuilt on each relabeling's calls and driver genes.
- **Descriptive:**
  - settings sensitivity: NET over threshold {0.05, 0.1, 0.2, 0.3} × resolution {0.5, 1, 2}, and OVL over cut {0.3,
    0.4, 0.5, 0.6, 0.7}, as ARI against the default (median and range);
  - the share of programs named by an ancestor;
  - for V4, the programs (NET) that keep at least one call.

**Verdict:** NET against OVL and NET against RED, and NET-pos against OVL, per method and dataset, by the dominance
rule of section 3. Directions: W1 lower; W2, W3 and W4 higher; W5 lower. Tolerances: W1 0.02; W2 and W3 0.5; W4 and W5
0.05. A criterion undefined for either side is skipped.

### 6 · Q3: discovery

- **Program tables** for each method × dataset (NET and NET-pos). Each row gives:
  - the name, calls and driver genes (top 5);
  - the peak and trough segments and the position score class;
  - coherence and replication r;
  - the member calls.
- **Union network.**
  - Calls are F ∪ C ∪ N (condition). Driver genes follow F's rule for every set, so the edges rest on one gene
    evidence.
  - Settings are NET's.
  - For each program: which methods call at least one of its calls.
  - Methods are compared at the program level: programs found by each, by pairs and by all three, and the program
    Jaccard for each pair.
  - This separates "different biology" from "the same biology through different sets".

### 7 · Application: notebook 91's injury screen (descriptive)

NET is run separately on 91's level-up, level-down and location sets.
- Driver genes:
  - level sets: members with `level change` True and the sign of the set's L;
  - location sets: members with `location change` True.
- Programs are named as in section 4.
- The peak of each call is 91's set `peak` (in 91's segment units).
- Reported: the programs, their sizes and names, and the share of calls in programs of ≥ 3 calls.

The question is whether this makes 91's 774 calls readable. No verdict.

### 8 · Outputs

Under `results/pt_pathway_networks/`, `tables/` and `figures/`:
- **Tables:**
  - `inputs.csv`, `segment_cuts.csv`;
  - `row_n_{ds}.csv` (per set: p and q per list, call, calibration p, per partition), `row_n_genes_{ds}.csv`;
  - `calling_criteria.csv` (V1–V4 and the descriptive columns, F, C and N, both datasets), `calling_verdicts.csv`;
  - `calls_{ds}.csv` (every call of every method, with its driver genes and A2 attributes);
  - `programs_{method}_{grouping}_{ds}.csv`;
  - `grouping_criteria.csv`, `grouping_verdicts.csv`, `grouping_settings.csv`;
  - `union_programs_{ds}.csv`, `union_overlap.csv`;
  - `screen91_programs.csv`;
  - `hub_genes.csv`;
  - `run_record.json`.
- **Figures** (plots only):
  - 93a: the NET network of F's calls on each dataset, nodes coloured by peak position;
  - 93b: program profiles with the two case specimens;
  - 93c: the grouping criteria;
  - 93d: the calling criteria;
  - 93e: the union network coloured by method;
  - 93f: 91's screen networks.

### 9 · What is fixed now, and caveats

Everything above, including the thresholds (gene q 0.05, a_g 0.5, Jaccard 0.1, resolution 1, r 0.8, 10 hub genes),
seeds (1234, 93) and tolerances. A later check is labelled post hoc.

**Caveats:**
- **Row N in the relabelings.** It is near-trivially specific when they have no selected genes; 75 reports the same
  for its q ≤ 0.05 module row.
- **Replication in human vs mouse** is section to section within one donor.
- **The grouping criteria** measure organisation, not truth.
- **No network program is a test.** Program NCDR reuses the set-level relabelings.
- **Reactome and KEGG** releases differ from the 2022/2019 libraries; unmatched sets are counted.

## Notebook 93 results

Not yet run.
