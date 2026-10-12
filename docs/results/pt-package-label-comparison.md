# The species comparison by label on the package units (notebook 73)

Part of the switch to the package pseudospace; see `pt-package-pseudospace.md`. The protocol was fixed before
any fit it governs.

## Notebook 73 protocol: the label-based species comparison on the package units

Draft by ws73, 2026-10-12, before any private-data run of notebook 73. It includes the coordinator's decisions of
2026-10-12 (sign-off constant, integration-gene arms, 64 on the package units, mt-free library, 44 B in and 39
retired, no change to notebook 70). It is frozen once the coordinator supplies the label fingerprint.

**Why.** Under coordinator decision D1 the paper uses one unit set. The label-based, coordinate-free species results
(notebooks 42, 43 and 49, with 46, 47, 48 part B and 64 as sensitivities) must therefore move to the package's units,
counts and automatic S1/S2/S3 labels. Notebook 42's reference transfer is the first independent test of the package's
human labels.

**Question.** Do the label-based results hold when only the inputs change? The inputs are the units, their counts,
the S1/S2/S3 labels, the expression glomeruli and the integration genes. No coordinate is used anywhere. Everything
is descriptive: two male control mice and two cortex sections from one male human donor, with specimens as the
replicates.

**Labels.**
- Notebook: `analysis/notebooks/73_pt_package_label_comparison.ipynb`.
- Logic version: `73.package_labels.1`.
- Results: `results/pt_package_label_comparison/`, with `tables/<part>/`, `figures/`, `stage_cache/` and
  `run_record.json`.
- Dataset: hvm. The aki export enters only the 44 B part.

**Two library arms** (coordinator, decision 4).
- **Primary: the mt-free library, the paper-wide rule since notebook 69.** The mt genes are the ortholog map's `mt-`
  symbols that are in the export. As in notebook 69, they leave the library and the tested genes. Here the library
  means:
  - every pseudobulk (log2 CPM, log ratios, levels);
  - every per-unit normalisation (log1p CP10k);
  - the spillover share.

  Every part runs on this arm, and the M rules apply to it.
- **Secondary: mt genes in the library and in the gene space** (43's rule), which is the pure input swap.
  - It runs 42 (c) and the spillover, 43 in full, and 49 part A with the gene bootstrap.
  - Old → secondary isolates the unit and label change; secondary → primary isolates the library change.
  - The scratch equivalence check (section 2) uses this arm.

### 0 · Gene model (D2)

None of these analyses fits a gene GAM, so D2 changes nothing here.

| Notebook | Statistic | D2 |
|---|---|---|
| 42 | Spearman or centred Pearson transfer on log1p CP10k profiles; count shares; polygon geometry | none |
| 43, 46, 47 a/b/b′/d/f, 64 | specimen × segment pseudobulk log2 ratios (prior 0.5), squeezeVar-moderated per-gene contrasts; log1p means only as a sensitivity | none |
| 49 | as 43, plus an SVD of log1p for the co-expression blocks | none |
| 48 B | matched set tests on class indicators; residual correlations of log1p | none |
| 44 B | percentile ranks of pseudobulks | none |

The parts that rested on Gaussian gene fits or on the retired 22-pathway list are not carried over (section 5).
These are 43 §7 (notebook 12's T_spatial), 47 c/e and 48 A (notebook 37 runs). Notebook 39 is retired (section 5).
The mt-free library is D2's offset rule applied to pseudobulks: the primary arm above.

### 1 · Inputs: what each old input becomes

| Old input | Package replacement |
|---|---|
| Notebook 13 PT structures plus `tubule_by_gene/*_caleb.h5ad`, through `rebuild_pt_expression` | `package_pseudospace/hvm/pt_units.h5ad`: every PT unit, no support trim. Gene space: the export's genes, all measured in both species unless a var `measured_in_both` says otherwise; the mt genes are dropped in the primary arm. Library: the total over the gene space (primary without mt; secondary with mt, 43's rule) |
| Reviewed labels (notebook 13's `REVIEWED_13`) | the export's `segment_class`: automatic pass-1 labels after the coordinator's per-cluster sign-off (G2) |
| Notebook 13 cohort of every class (centroids, `feature_index`, the Glomerulus cluster) | `hvm/all_units_positions.csv` (label, family) joined to `runs/hvm/02_qc/units_qc.h5ad` obs and counts |
| Depth anchors: notebook 03 pass-1 `Glomerulus` structures | QC units of family `Glomerulus` (label `Podocyte`). Stage 03 drops them before pass 2, but they are in the pass-1 table |
| Fine GeoJSON plus `_read_and_validate_segmentation`, which reads `*_kept_tubules_labeled_fine` only | the stage-01 segmentations (`Ctrl_1A*_v4`, `HUK1_*_v2`), read with the package's `load_geometries` and `polygons_for_obs` |
| Morphological glomeruli (class `glomerulus` in the GeoJSONs) | unchanged. Mouse geometry now comes from the units' own file, so 42's frame-offset check is the identity |
| Zero-bin polygons (`tubule_by_gene` `n_spots`) | `runs/hvm/01_units/{sample}_units.h5ad`, one row per polygon, and `binning_summary.csv` |
| Notebook 13's 611 PCA genes and 752 Harmony HVGs | `runs/hvm/03_integration/pass1_harmony_features.csv`: the embedding whose clusters became the labels |
| Native per-specimen matrices for the alternative pairings (64) | `runs/hvm/01_units/{sample}_units.h5ad`: native symbols, all genes, the same unit ids without the sample prefix. They go through `rebuild_pt_expression` unchanged |
| `tubule_by_gene` var `gene_ids` (Lake and probe mapping) | the stage-01 var `gene_ids`, from the same 10x feature table |
| Notebook 62's `pt_expressed` and `log2_cp10k` (64's choice of main paralog) | recomputed on the package PT units by notebook 62's rule (section 3) |

**Unchanged inputs.**
- The ortholog map.
- Notebook 38's `probe_counts.csv` and `lake_cortex_segment_means.csv`.
- The Lake v1.5 h5ad.
- The Census pseudobulks and the paralog extension.
- GSE150338, GSE212213 and GSE56743.
- Notebook 62's orthogroups and probe status.
- The pathway libraries.

**Read only, for the comparison.** `pt_revision_labels/`, `pt_revision_classes/{reviewed,crossfit}/`,
`pt_revision_addenda/`, `pt_revision_addenda_pathways/`, `pt_conservation_strength/`, `pt_paralog_pairing_check/`,
`pt_revision_state/`, and 70's `hvm/tables/` (`stage_counts.csv`, the label crosstab).

**Loader** (new `pseudospace/package_labels.py` and `tests/test_package_labels.py`, synthetic fixtures).
`package_inputs.py` needs no change: it is coordinate-specific (support trim, universe) and is used only through
`read_pt_units`.
- `label_units(results_root, dataset, exclude_genes=())`:
  - takes `read_pt_units` and returns 43's objects: obs with sample, species, segment S1/S2/S3, `feature_index`,
    centroids, `n_spots`, `total_counts` and `old_structure_id`;
  - also the gene space, CSR counts, library, log1p CP10k, and the pass-1 fingerprint from `uns['provenance']`;
  - `exclude_genes` leave the gene space and the library: the mt genes in the primary arm, none in the secondary.
- `all_units(results_root, dataset, genes=(), exclude_genes=())`:
  - every QC unit, with label, family, `feature_index`, centroids, `n_spots`, `total_counts` and `seg_class`;
  - the label and family come from `all_units_positions.csv`; the other columns are read from
    `runs/{dataset}/02_qc/units_qc.h5ad`, since notebook 70 is not changed;
  - raw counts for `genes` (the spillover markers), and the library without `exclude_genes`;
  - it raises if a unit is missing from `units_qc`.
- `cortical_depth(units, anchors, k=3)`: 43's mean centroid distance to the k nearest anchors of the same sample.
- `segmentation_glomeruli(path)`: the `glomerulus` polygons and µm per px (42's median √(Area µm² / Area px)).
- `native_sources(results_root, dataset)`: the stage-01 paths, plus a `pt_obs` in the format `rebuild_pt_expression`
  expects.
- `harmony_features(results_root, dataset, which)` and `label_panel_genes(provenance)`: the PT-S1/S2/S3 panels, 9
  genes.

### 2 · Guards (not cached)

- **G1.** `read_pt_units(…, 'hvm')`, the export contract.
- **G2 · The label sign-off.**
  - The config cell holds `CONFIRMED_PASS1_FINGERPRINT = None`, a placeholder. The coordinator fills in the
    membership sha1 after reviewing notebook 70's `cluster_labels_for_review.csv`.
  - The notebook stops while the constant is the placeholder, or if it differs from the pass-1 fingerprint in
    `uns['provenance']`.
- **G3.** Each specimen has PT-S1, PT-S2 and PT-S3. Each mouse specimen has at least 3 `Glomerulus` units.
- **G4.** Every unit polygon joined by `feature_index` has its centroid within 2 px (42's rule), for PT units and for
  `Glomerulus` units.
- **G5.** `rebuild_pt_expression` on the stage-01 sources with the paper map reproduces the export exactly: the same
  gene space and identical counts. 64's native route is therefore the export.
- **G6.** Internal reproduction to 1e-9 against 73's own 43 part, in each arm it runs: the index and symmetric class
  counts from 46's `evaluate` (labels, all genes), from 49's estimator, and from 64's `primary_analysis` (baseline).
  These replace the old guards against the published tables.
- **G8.** In the primary arm, each library falls by exactly the counts of the mt genes it drops. The number of mt
  genes leaving the library and the tested genes is reported (notebook 71's G3).
- **G7.** 64's Census extension reproduces its control genes (64's guard).

**Equivalence check before the package run** (phase 2, untracked scratch script, as for notebook 71). The notebook's
code is fed the *old* inputs, in the secondary arm: notebook 13's PT structures, the reviewed labels, notebook 03's glomeruli, and the 611
and 752 genes. It must reproduce the published tables exactly:
- 42: `transfer_agreement`, `glomerulus_distance_by_label`, `glomerulus_attached_pt`;
- 43: `conservation_index`, `gene_classes_symmetric`, the class counts;
- 46: the cross-fit index and counts;
- 47: `no_integration_genes_*`;
- 48: `class_pathway_enrichment_symmetric*`;
- 49: `strength_curves`, `decisions`;
- 64: the baseline.

The notebook does not repeat this check, because it needs the old inputs.

### 3 · Inherited definitions (unchanged)

**42.**
- Polygon join: by `feature_index`, with the centroid within 2 px.
- Distance: boundary to boundary, to the nearest glomerulus polygon (STRtree), in µm.
- Glomeruli: morphological glomeruli are primary; expression glomeruli (the Glomerulus family) are the sensitivity.
  Glomerulus-sized means a median area ≥ 3 × the median PT area.
- Attached: within 2 µm (primary) or 10 µm (sensitivity). The S1 odds ratio uses +0.5 cells.
- Transfer genes:
  - measured in both species;
  - detected in ≥ 10% of the species' PT units;
  - reference maximum ≥ 1 log2;
  - a largest segment difference ≥ 1 log2.
- Transfer scores:
  - primary: Spearman argmax;
  - sensitivity: the centred score (per-specimen PT mean against the reference centred over its segments).
- References: Lake/KPMP cortex for human; GSE150338 microdissection (primary) and Census male snRNA for mouse.
- Reading: the human labels are **supported** if human primary agreement is ≥ 0.70 and differs from mouse
  microdissection agreement by ≤ 0.10. Otherwise the support is **weak**.
- Spillover markers: 42's list. A marker is used if it is detected in ≥ 1% of non-PT units. The index is its share
  of the unit's counts.
- Registration QC: zero-bin polygons inside the hull of the covered centroids; coverage = bins × 4 µm² / area.

**43.**
- Units: every PT unit, with no coordinate filter.
- Cortical-like mouse S3: depth at or below the 95th percentile of the same specimen's S1 depth.
- Contrasts: S2 − S1, S3c − early and S3 − early, where early pools S1 + S2 counts.
- Gene filters:
  - ≥ 20 counts in every pseudobulk;
  - external datasets: ≥ 50 nuclei and a level ≥ 3 log2.
- Disattenuated r: winsorized at 1/99%, over ≥ 200 genes.
- Conservation index: primary, external, and alternative ceilings, with 1,000 gene bootstraps (seed 43).
- Wording: 0.5 and 0.8.
- Rat: used only if reliability is ≥ 0.5 and ≥ 70% of genes map.
- Classes: α 0.05, effect 0.5, margin 0.3, confirmation α 0.05.
  - The 8 settings; the symmetric setting is A3b + A8.
  - Probe-balanced genes (notebook 38).
  - D3 asymmetry rule: a ratio ≥ 2.
  - Mouse-only distribution: |mouse| ≥ 1.5 with |human| ≤ 0.3.

**49.**
- Donor folds: sha256 round-robin (Lake healthy 4/3, Census male 6/6).
- Strength: the maximum percentile rank of |fold-A gradient|.
- Bins: 0, 50, 80, 90, 95, 99 and 100%. BOTTOM is below 80%; TOP is 95% and above.
- Settings: q 0.01; minimum reliability 0.2; at most 10% of draws missing.
- Bootstraps: 1,000 gene and 1,000 module-block draws (seed 49).
- Modules: SVD with 50 components, k-means k = 100 (50 and 200 as sensitivity).
- Rules:
  - (i) TOP − BOTTOM ≥ 0.25, with the interval excluding 0;
  - (ii) BOTTOM ceilings with lower bounds ≥ 0.3 and an index upper bound < 0.5.
- The four verdicts.
- Part B: q ∈ {0, 0.01, 0.05} × reliability {43, matched}, with gene, module (k 50/100/200), donor and donor + gene
  bootstraps, and the leave-one-donor-out human ceiling.
- Also run: the secondary schemes, direction agreement by strength, and the external pair binned on our strength.
- The mouse-strength / full-Lake diagnostic is now **pre-specified**, but stays secondary.

**46.**
- Gene halves: sha256 by gene name, so each gene keeps its old half.
- Labels: transferred by the centred score on one half's reference-zonated genes; evaluated on the other half.
- Settings:
  - (i) human Lake, mouse package;
  - (ii) human Lake, mouse microdissection;
  - (iii) package labels.
- Claim stands: index ≤ 0.5 with an upper bound < 0.6, for S2 − S1 and S3 − early, in every setting and half.
- Unskewed: ratio < 2 in S2 − S1 and S3c − early.

**47.**
- **(a)** Evaluation genes exclude the integration genes; 46's evaluation and decision apply (coordinator,
  decision 2).
  - Primary gene set: the pass-1 Harmony features.
  - Secondary gene set: those plus the 9 PT panel genes. These are new circularity sources: the panels named the
    labels.
  - The pass-2 features are counted descriptively.
- **(b, b′)** Transporters, `^(Slc|Abc|Aqp)`, probe-balanced, on the union of the S2 − S1 and S3c − early classes.
- **(d)** Integration genes by symmetric class.
- **(f)** Reversal counts under each label setting:
  - package labels;
  - the full-gene centred transfers;
  - the cross-fit halves;
  - the integration-gene-free sets.

**48 B.**
- Universe: the classified genes. Sets: 10–300 classified members.
- Strata: 3 (5 as sensitivity) on mean, detection and minimum species detection over all PT units.
- Residual correlations: after removing the specimen × segment-region means.
- Tests: 9,999 matched draws (seed 40), then the joint test with VIF.
- Calls:
  - enriched: q ≤ 0.10 with effect > 0;
  - not expression-driven: also q ≤ 0.10 in the 5-bin and upper-expression runs;
  - programs: overlap cut 0.5.
- Also run on the probe-balanced genes.
- Naming rule: q ≤ 0.10, not expression-driven, and ≥ 3 class genes.

**64.**
- Families: notebook 62's orthogroups (HCOP ≥ 3; 1:many, many:1, many:many) and probe status.
- Main paralog: the highest pooled CP10k among the members that are PT-expressed (detected in ≥ 2% of units in both
  specimens of the species) and probe-clean.
  - The expression is recomputed on all package PT units, from native stage-01 counts (coordinator, decision 3).
  - The library of the CP10k is the species panel without the mt genes.
  - The choices that differ from notebook 62's are counted.
- Pairings: the swap and the drop (`ortholog_pairing`); probe balance is recomputed for every pair.
- Swapped-in genes missing from the Census extension are excluded and counted.
- Unaffected: every index changes by < 0.03, and each zonated class count by < 10% (≤ 2 genes when the count is < 20).

**44 B.**
- Signature: Lake aPT against healthy PT nuclei, top 30 up and top 30 down. Donors need ≥ 30 aPT and ≥ 100 healthy
  PT nuclei, and a gene needs mean CPM ≥ 20.
- Score: the mean within-sample percentile of the up genes minus that of the down genes.
- Platform offset: our control mice against the healthy Census mice.
- Universe: hvm PT ∩ aki PT ∩ Lake ∩ Census genes, without the mt genes in the primary arm.
- Whole-PT pseudobulks per specimen come from both package exports. The aki export supplies the AKI − control score.

### 4 · What changes unavoidably

- **U1 · Units and QC.**
  - New: package stage 01 (mouse: every v4 polygon; human: v2) and stage-02 QC.
  - Old: notebook 01 on `kept_tubules_labeled_fine`, plus the 03/13 structure filters.
  - Unit numbers and depths differ (old: 12,866 PT; 2,852 / 2,808 / 3,438 / 3,768).
- **U2 · Labels.**
  - Package: per-cluster panel argmax. Each PT segment is named by 3 mouse panel genes: S1 Slc5a2, Slc5a12, Gatm;
    S2 Slc22a6, Slc13a3, Cyp2e1; S3 Slc7a13, Slc22a7, Cyp7b1. Hand review is replaced by the coordinator's
    per-cluster sign-off.
  - Lead genes Gatm and Slc7a13 are label-defining. Gene tables carry a `label_panel_gene` flag.
- **U3 · Gene space.** The export's shared 1:1 orthologs (about 15,567), against notebook 13's genes measured in all
  four specimens. The Lake pseudobulk is recomputed once over all genes and mapped to the gene space (64's route).
- **U4 · Expression glomeruli and depth anchors.** The package's Glomerulus family.
- **U5 · Integration genes.** The package's per-group HVG intersection (pass 1), plus the panel-gene arm. Notebook
  13's 611 have no exact analogue.
- **U6 · Guards.** The guards against the published tables become internal guards (G6). Code equivalence moves to
  the scratch check.
- **U7 · 64's main-paralog choice** is made from package PT expression.
- **U8 · The library is mt-free in the primary arm** (notebook 69's rule). Unlike U1–U7, this is a chosen change.
  The secondary arm keeps the old library so that its effect can be read separately.

### 5 · Not carried over

These parts rest on the retired 22-pathway list, the Gaussian runs, or notebook 40's classes. Notebook 77 reruns them
on 71's calls where the outline needs them:
- 43 §7 (M7, robust-pathway categories);
- 47 (c) and (e);
- 48 A, 8b and 8c;
- 48 B's guard against notebook 40;
- 42's integration genes by notebook-40 class (47 d covers it).

**Notebook 39 is retired** (coordinator, decision 5), for three reasons:
- Its headline, the human/mouse spread (S2 − S1 0.58, S3c − early 0.84), is 43's noise-corrected amplitude, and 43
  runs in full here.
- Its projection ratios (0.15–0.23 on mouse references, 1.96 on Lake) serve only the biology-led R2 paragraph.
- Its AKI amplitude arm and §9 rest on `mouse_only_v5` labels and on the retired 22-pathway list.

### 6 · Outputs

**Tables.** Each part keeps its old table names under `tables/{42,43,46,47,48B,49,64,44B}/`, so readers change only
the folder. The secondary arm writes to `tables/mt_in_library/{42,43,49}/`.

**Exports for 74/77.**
- `gene_classes_symmetric.csv`, with `probe_balanced` and `label_panel_gene`;
- `gene_gradients.csv`: per specimen, contrast and gene, plus the moderated species means;
- `segment_labels.csv`;
- `cohort.csv`, for items 23–24:
  - 70's stage counts, beside the old ones;
  - the quantiles of bins and UMI per unit, by specimen and species, for PT units and for all units.

**Comparison and record.** `compare_old.csv`, figures as in 42, 43 and 49, and `run_record.json` with the
provenance digest and the confirmed fingerprint.

### 7 · Comparison with the old results (descriptive)

- **Headline values: old, then secondary (mt in), then primary (mt-free).**
  - 42: transfer agreement and κ (human 0.71 / κ 0.54; mouse 0.87 / κ 0.80 and 0.92 / κ 0.87); distance Spearman
    (0.65 / 0.23); S1 odds ratio at 10 µm (4.5 / 2.0); human S1 share (0.55); spillover range.
  - 43: index 0.17 / 0.28 / 0.29; union classes 196 / 75 / 41 / 76; amplitude ratio (S2 − S1 0.58); the wording and
    the D3 decision.
  - 49: top-1% index 0.57 / 0.58 / 0.55; rules and verdict; direction agreement from 0.57 / 0.48 / 0.58 to 0.89 /
    0.84 / 0.89.
  - 46: maximum index 0.31 (upper bound 0.38).
  - 47: without integration genes, at most 0.24 (0.33); transporters 33 / 196 conserved and 33 / 40 against 196 / 271.
  - 48 B: the naming-rule table.
  - 64: Δ ≤ 0.002 and 0.007.
  - 44 B: donor percentile 0.5 among Lake normal donors.
  - Diagnostic: the median gradient of each contrast, per specimen. It shows the constant shift that the library
    change induces.
- **Genes.** Spearman of per-gene gradients, old against new, per species and contrast, on genes in both runs. Also a
  crosstab of the old symmetric union class against the new one.
- **Labels.** Notebook 70's crosstab is cited, not recomputed. A new matched-unit comparison runs on the units that
  carry an `old_structure_id`: the reference transfer's agreement with the old reviewed labels and with the package
  labels, on the same units and the same transferred labels.

### 8 · Material changes from the old results (pre-specified)

Each of these is stated in the notebook summary and in the ledger. Nothing is tuned.
- **Which arm.** Every rule applies to the **primary arm** (mt-free library), against the old values.
- **Attribution.** When a rule fires, the secondary arm says whether the units and labels or the library move the
  result. That reading is descriptive.
- **M1.** 42's reading of the human labels changes (old: weak support, 0.71 against 0.87), or agreement falls by
  ≥ 0.10 in either species.
- **M2.** In either species, the label-order Spearman with morphological glomerulus distance is ≤ 0, or the S1 odds
  ratio at 10 µm is ≤ 1.
- **M3.** 43's D2 wording changes (old: "largely different genes (stands)"), or any primary index moves by ≥ 0.10.
- **M4.** 43's D3 decision changes (old: no consistent asymmetry).
- **M5.** 49's rule (i) fails in S2 − S1 or S3 − early under either bootstrap, or its verdict changes (old: "agreement
  is confined to genes that are reliably zonated").
- **M6.** The claim "largely species-specific" fails in any 46 or 47 (a) evaluation.
- **M7.** Any of these:
  - direction agreement among the genes zonated in both species moves by ≥ 0.10 (old 0.72);
  - transporters no longer agree in direction more often than all genes (old 0.83 against 0.72);
  - the probe-balanced species-only naming changes (old: mouse-only Tyrosine metabolism through Comt, Fah and Maoa).
- **M8.** 64's swap verdict becomes "affected".
- **M9.** In 44 B, either human section leaves the Lake normal range (old: percentile 0.5, within range).

### 9 · Effort

- **Size.** One notebook, about 1,900 lines. The analysis code is copied verbatim from 42–64 so that the swap can be
  audited by diff.
- **Run time.** About 2 h cold:
  - one pass over Lake;
  - the 49 bootstraps (about 25 min, plus the secondary arm's gene bootstrap);
  - 48 B's 9,999-draw runs.
- **Phase 2 waits** for notebook 70's export and the coordinator's fingerprint.
