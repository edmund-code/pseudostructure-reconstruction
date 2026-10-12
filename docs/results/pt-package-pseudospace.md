# The package pseudospace: notebooks 70 and 71

On 2026-10-11 the user froze Caleb's package `pseudospace_reconstruction` (stages 01–04: bins → units, QC, two-pass
Harmony with labels, elastic principal curve) as the method that builds the pseudospace for both datasets. On
2026-10-12 the user asked to "switch all results to use calebs package now", without waiting for the reference
`v4fix` segmentations. The package replaces the PT scFates coordinate (SCF13, human vs mouse) and the mouse-only PT
DPT (AKI vs control), together with their unit sets and reviewed labels.

- Notebook 70 builds the package pseudospace for both datasets and exports the PT units.
- Notebook 71 reruns notebook 69's analysis on it, unchanged, and freezes the comparison page's post hoc analyses as
  secondary sections.
- Later notebooks rerun every other coordinate- or label-dependent result through the same loader
  (`pseudospace/package_inputs.py`). Every later analysis uses the NB gene model with the mt-free offset. The
  package's own DPT, on the same embedding, is the sensitivity coordinate.

Both protocols below were fixed before any fit they govern. Notebook 72 (validation of the package coordinate) and
notebook 73 (the label-based species comparison on the package units) have their own documents:
`pt-package-coordinate-validation.md` and `pt-package-label-comparison.md`. Notebook 74 (continuous against discrete) is in
`pt-package-continuous-vs-discrete.md`, notebook 75 (pathway methods and their specificity) is in
`pt-package-pathway-methods.md`, notebook 76 (AKI along the PT, Figure 6) is in `pt-package-aki.md`, and notebook 77 (the human vs mouse calls
downstream) is in `pt-package-species-calls.md`.

## Notebook 70 protocol: the package pseudospace for both datasets

Frozen on 2026-10-12, before any private-data run of notebook 70. The user froze Caleb's package
`pseudospace_reconstruction` 0.1.0 (commit 59da04b) as the method that builds the pseudospace and asked to switch
all results to it now, on the local segmentations.

**Question.** What do stages 01–04 of the frozen package give on our inputs, for the AKI dataset and for the human
vs mouse dataset? Notebook 70 builds both, exports the PT units for the downstream loader, and describes how the
new coordinate and labels agree with the old ones (SCF13; the mouse-only PT DPT). Nothing here selects or tunes the
method. The agreement numbers are descriptive.

**Inputs.**

| Dataset | Sample | Group | Segmentation (polygons) | Bins (2 µm) |
|---|---|---|---|---|
| aki, hvm | Ctrl1A2 | Control / Mouse | `Ctrl_1A2_v4.geojson` (13,630) | `$PSEUDOSPACE_VISIUM_ROOT/Ctrl1A2/square_002um` |
| aki, hvm | Ctrl1A4 | Control / Mouse | `Ctrl_1A4_v4.geojson` (14,592) | `$PSEUDOSPACE_VISIUM_ROOT/Ctrl1A4` |
| aki | IR2A2 | IR | `IR_2A2_v4.geojson` (18,359) | `$PSEUDOSPACE_VISIUM_ROOT/IR_2A2/square_002um` |
| aki | IR2A4 | IR | `IR_2A4_v4.geojson` (17,691) | `$PSEUDOSPACE_VISIUM_ROOT/IR_2A4/square_002um` |
| hvm | HUK1_COR1 | Human | `HUK1_COR1_v2.geojson` (8,747) | `HUK_COR_1/outs/binned_outputs/square_002um` |
| hvm | HUK1_MED1 | Human | `HUK1_MED1_v2.geojson` (14,186) | `HUK_MED_1/outs/binned_outputs/square_002um` |

- Segmentation paths are relative to the data root. Both human slices are healthy cortex; "MED" is a misnomer.
- The human v2 files have the same feature schema as the mouse v4 files: `classification.name` (`tubule`,
  `glomerulus`), the `Score` measurement and `instance_id`. Neither species has `properties.id`. The package then
  names units `unit_<feature_index>`, as in the AKI trial.
- The ortholog map is `results/minimal_pt_scfates/ortholog_map_used.csv`: 17,449 HCOP one-to-one pairs, read as in
  `coordinate_inputs.pt_inputs`.
- The old coordinates are compared only, never used as inputs:
  - hvm: SCF13, `results/pt_reconstruction_v2/primary_pt_coordinate.csv`. Its 12,866 structures are the obs index
    of `cross_species_pt_scfates.h5ad`. The mouse structures come from `*_kept_tubules_labeled_fine.geojson` and
    the human structures from `HUK1_*_v2.geojson`.
  - aki: the mouse-only PT DPT, `results/mouse_only_v5/all_mouse_tubules_scanpy_dpt.h5ad` (`pt_subset_scanpy_dpt`).
    Its structures come from the same v4 files as the new run.

**Deviations from Caleb's reference run.**
1. **Segmentations.** The mouse segmentations are the local v4 files: 64,272 polygons, against 64,080 in his
   `v4fix_*`. The human segmentations are the raw panoptic output.
2. **AKI config.** It is his shipped `mouse_kidney_aki` config with the trial's samples, and `gene_sets.dir` is
   relative to the data root (our Enrichr JSONs). Stage 07 alone reads that entry.
3. **Two species.** The guide requires shared 1:1 names before stage 02, and we implement the mapping. Stage 02 of
   hvm runs through the package's Python route (`stage02_inputs` → `qc.run_qc` → `save_stage02`), the same calls
   as its notebook 02. It gets the renamed units, plus `bin_dirs` that point to ortholog-space copies of the bin
   matrices.
   - Why the copies are needed: the split test (doublet tier 2) re-reads each sample's bin matrix, and
     `build_bin_table` requires that bin totals sum to the unit totals.
   - Each copy holds the shared genes only, under mouse symbols. It keeps the same barcodes and links to the
     original `spatial/`.
   - The CLI then resumes at stage 03, because stage 02 has its manifest.
4. **hvm `qc.bin_area_px`.** The package default is derived per sample from `spot_diameter_fullres` (about
   20.58 px²). AKI keeps 20.6.
5. **Human segmentation classes.** They are `tubule` and `glomerulus`, not the mouse
   proximal/distal/collecting classes, so the human MAD strata are coarser.
6. **Nephron-wide curve.** The package fits one curve through every kept nephron unit, so the PT is a stretch of
   it. The old coordinates (SCF13, the mouse-only PT DPT) were fitted on PT units only. Compare positions by rank
   within the PT, not by absolute value.
7. **Addition to Caleb's shipped config (coordinator, 2026-10-12).** `parameters: dpt: run: true` in both
   datasets. Stage 04 then also writes the package's diffusion-pseudotime baseline on the same embedding and
   root, with defaults: max(30, √n) neighbours, 15 diffusion components, 10 DCs, clipping at the 5th–95th
   percentile, then [0, 1].
   - The curve does not change: DPT is computed after it.
   - It is exported as `dpt_pseudotime`, oriented S1 → S3 (flipped and recorded if needed, NaN for units
     disconnected from the root). It is a sensitivity coordinate only; `pseudospace` stays primary.
8. **Build.** Linux x86-64: Python 3.11.16, numpy 1.26.4, scikit-learn 1.9.0, scanpy 1.11.5, R 4.5.3, harmony
   2.0.5. The guide says stage 03 onward can differ slightly between builds.

**Parameters.** These are the AKI config verbatim. hvm differs only where noted.
- **Segmentation.** `coordinate_scale` 1, `id_property` id, `class_property` classification, `score_measurement`
  Score.
- **Gene patterns.** mt `^mt-`; ribo `^(Rpl|Rps|Mrpl|Mrps)`. In hvm these apply to mouse symbols after the rename.
- **qc.**
  - `bin_area_px` 20.6 (hvm: derived) and `doublet_exempt_classes` [glomerulus].
  - Defaults: `min_total_counts` 500, `min_bin_coverage` 0.90, `mad_nmads` 3 on [log total, log per bin,
    complexity residual] within sample × seg_class (`mad_min_group_size` 50).
  - Split test: 99 shuffles, α 0.05. Discard tiers: split test, morphology (solidity < 0.75), neck.
  - Normalisation to 10⁴; seed 0.
- **integration.**
  - `backend` r and `label_source` pass1.
  - Defaults: HVG mean 0.0125–3, dispersion 0.5, groupwise intersection; PCA 50, Harmony on 50 PCs with theta 6,
    λ 1, 30 iterations.
  - Clustering: 2,000 HVGs, 30 neighbours, Leiden resolution 2 in both passes; `min_label_score` 0.5,
    `min_label_margin` 0.1.
  - Pass-2 gene filter: 5% of units and 20 counts.
- **dpt** (addition, see deviation 7). `run` true; other settings are the package defaults.
- **Panels.** Caleb's shipped mouse kidney panels, unchanged. That includes `rejected_labels` ATL → Rejected_ATL
  and `drop_classes` [Glomerulus, Vessel, Stroma, SmoothMuscle, Immune]. There are no `label_overrides`.
- **curve.** `n_dims` 50; defaults `n_nodes` 50, `lam` 0.01, `mu` 0.1, `end_stretch` 2, `seed` 0. The root is the
  medoid of the lowest-marker-axis PT cluster. The axis anchors are early [Cyp2e1, Gatm, Slc13a3, Slc22a6, Slc5a12,
  Slc5a2] and late [Aqp2, Aqp4, Atp6v0d2, Rhcg, Slc14a2, Wnt7b]. (Caleb ships DPT off; see the dpt line).
- **Groups.** AKI: reference Control, case IR. hvm: reference Mouse, case Human. Stages 05–07 are not run.

**Ortholog step** (`pseudospace/package_pseudospace.py`).
1. Stage 01 runs per sample on the original bins, so its outputs keep native names.
2. The map is validated as in `rebuild_pt_expression`: nonempty, no blank symbols, and no duplicated symbol in
   either column, compared case-insensitively.
3. A pair is shared when its mouse symbol is a stage-01 var name of every mouse sample (exact match) and its human
   symbol is a var name of every human sample (case-insensitive, as `map_human_to_mouse_space`).
4. Human columns are renamed to mouse symbols. Every non-shared gene is dropped in both species, and the gene order
   follows mouse stage 01. Expected: about 15,567 genes.
5. Panel genes without a shared 1:1 ortholog: Akr1b3 (ATL panel). Nphs1 is not in the mouse probe set, and
   Cox2/Fxyd4/Nephrin/Podocin match no gene, as in the reference. Every PT, axis-anchor and QC-region gene is
   shared.

**PT extraction and export.**
- The PT units are the stage-03 `pass2.h5ad` units with `broad_tubule_marker_call == 'PT'`, from pass-1 labels.
- X holds raw integer counts from `02_qc/units_qc.h5ad` `layers['counts']`: every gene kept through QC.
- obs, var and `uns['provenance']` follow the export contract.
- `all_units_positions.csv` lists every QC-retained unit. Units that stage 03 dropped have no pseudospace.

**Label review file.** `{aki,hvm}/cluster_labels_for_review.csv` has one row per cluster of each pass. Each row
holds the label, family, n units, units per sample, the top three panel scores with the package's margin and
flags, and the top 10 markers.
- The markers are the package's own ranking, rebuilt with its functions: Wilcoxon `rank_genes_groups` on that
  pass's clustering feature space. Pass 1 uses `expression_gene_filter` → `clustering_feature_space` on
  `pass1_labelled.h5ad`; pass 2 uses `clustering_feature_space` on `pass2.h5ad`.
- The first five must equal the package's `pass*_cluster_top_markers.csv`.
- No override is applied. The coordinator confirms or flags each cluster (AGENTS.md coarse-label checkpoint).

**Old-structure matching.**
1. Each old structure gets its source polygon through its feature_index, verified by centroid within 1 px.
2. Old and new polygons are paired when IoU ≥ 0.5 and each is the other's best match.
3. Where old and new come from the same file (human; all of AKI), this is expected to reduce to feature_index
   identity with IoU 1. The identity share is reported.
4. Match rates are reported both ways per sample: new PT units with an old structure, and old structures with a
   new PT unit.

**Descriptive checks** (per dataset).
- **Stage counts.** Polygons, polygons with bins, units after QC, units in pass 2, clusters per pass, and PT-S1/S2/S3
  units per specimen.
- **Label checkpoint.** For every pass-1 and pass-2 cluster: label, family, size, the top three panel scores,
  margin, top markers, and in hvm the species share. Flags:
  - weak (best score < 0.5);
  - the package's cross-family ambiguity;
  - a near-tie at any margin < 0.1;
  - Unassigned or Rejected;
  - a cluster that is more than 95% one species (hvm).
- **PT order.**
  - Adjacent-segment AUCs P(S2 > S1) and P(S3 > S2) along `pseudospace`, pooled and per specimen.
  - The orientation and the root unit, with its label and sample.
- **PT on the curve.**
  - The 1–99% range of PT positions.
  - The PT/non-PT boundary c that maximises balanced accuracy of "PT iff pseudospace ≤ c".
  - The PT units beyond c and the non-PT units before c, pooled and per specimen.
- **Agreement with the old coordinate** (matched PT units). Spearman of `pseudospace` against SCF13 (hvm) or the
  mouse-only PT DPT (aki), pooled and per specimen; also on notebook 69's common-support subset.
- **Label agreement.** New and old S1/S2/S3 for matched units: % identical, Cohen's κ and the crosstab.
- **AKI determinism.** Comparison with the scratch trial: stage summaries, fingerprints, labels and the maximum
  |Δ pseudospace|. A difference is reported, not fixed.

**Stop conditions.** The run stops and reports, writing no export, if any of these holds:
- the package is not 0.1.0 at 59da04b, or R harmony is not 2.0.5;
- a human segmentation lacks `classification` or `Score`;
- `validate` fails;
- PT-S1, PT-S2 or PT-S3 is missing in any specimen;
- the PT segment medians are out of order (S1 < S2 < S3), pooled or in any specimen;
- the curve does not run S1 → S3.

No `label_overrides` are made without the coordinator's sign-off; doubtful labels are reported.

## Notebook 70 results

`analysis/notebooks/70_pt_package_pseudospace.ipynb` (logic version `70.package_pseudospace.1`) writes to
`results/package_pseudospace/`. It followed the protocol above. No stop condition triggered: the package is 0.1.0 at
59da04b with R harmony 2.0.5, PT-S1/S2/S3 are present and in order in every specimen, and both curves run S1 → S3 from
a PT-S1 root in Ctrl1A4. The AKI run is identical to the scratch trial: same fingerprints and labels, maximum
|Δ pseudospace| 0.

**Implementation, beyond the protocol.**
- Human vs mouse stage 02 ran through the package's Python route on ortholog-space copies of the bin matrices
  (deviation 3); 15,567 shared one-to-one genes.
- No segmentation file, mouse v4 included, has `properties.id`, so units are named `unit_<feature_index>`.
- The mouse bins live outside the data root; the notebook takes `--visium-root` (`PSEUDOSPACE_VISIUM_ROOT`).
- `pt_units.h5ad` holds every gene kept by QC; `var.in_pass2` marks the pass-2 genes.

**Stage counts.**

| Dataset | Polygons | With bins | After QC | On the curve | Clusters (pass 1 / 2) | PT-S1 / S2 / S3 |
|---|---|---|---|---|---|---|
| AKI | 64,272 | 43,984 | 36,445 | 32,939 | 20 / 18 | 4,965 / 6,069 / 3,979 |
| Human vs mouse | 51,155 | 31,533 | 27,395 | 25,048 | 18 / 15 | 6,338 / 2,361 / 3,248 |

Human PT-S2 has only 373 (HUK1_COR1) and 259 (HUK1_MED1) units.

**Order and the PT stretch.**
- Adjacent-segment AUCs along the curve, pooled (S1→S2, S2→S3): AKI 0.976, 0.985; human vs mouse 0.893, 0.975.
  Per specimen: mice 0.95–0.99; HUK1_COR1 0.89, 0.94; HUK1_MED1 0.85, 0.96.
- The PT occupies the start of the nephron-wide curve: raw 0–0.285 (AKI) and 0–0.202 (human vs mouse); 1.2% and 2.6%
  of PT units lie beyond the PT stretch.

**Agreement with the old coordinates** (matched PT units; descriptive).

| | AKI vs mouse-only PT DPT | Human vs mouse vs SCF13 |
|---|---|---|
| Spearman, pooled (per specimen) | 0.888 (0.850–0.930) | 0.908 (mice 0.913, 0.957; human 0.870, 0.876) |
| Spearman, package DPT baseline | 0.946 | 0.871 |
| Labels identical / Cohen's κ | 89.7% / 0.84 | 87.4% / 0.80 |
| Main label shifts | 15% of old S1 → S2; 11% of old S3 → S2 | 47% (COR) and 60% (MED) of old human S2 → S1; mice ≤ 6% |
| Old PT structures on a new PT unit | 84.7% (13.6% removed by package QC) | 89.3% |
| New PT units with an old structure | 96.5% | 96.2% |

Mouse polygons in human vs mouse were matched from the kept-tubule files to the v4 files by mutual-best IoU ≥ 0.5:
97.6% and 97.8% matched, median IoU 0.99. The human and AKI units come from the same files as the old structures, so
they match by feature index with IoU 1.

### Label checkpoint (coordinator, 2026-10-12)

Every pass-1 cluster of both datasets was reviewed against its panel scores and top markers
(`cluster_labels_for_review.csv`). All labels are confirmed as the package assigned them; no `label_overrides`. The
confirmed pass-1 fingerprints are human vs mouse `efafc9ac2e99` (18 clusters, 27,395 units) and AKI `bb4a80c1b484`
(20 clusters, 36,445 units). Flags kept on record:
- **Human vs mouse c2 and c3, PT-S1, 94% human.** They score clearly S1 on the panel (margins 1.26 and 1.32; PT-S2
  never above 0.36). Consequence: human PT-S2 is small, and about half of the old human S2 units are now S1.
- **Human vs mouse c17, PT-S1, 481 units, 99% HUK1_MED1.** Heat-shock genes (Hspa1a/b, Dnajb1) lead its markers,
  followed by PT-S1 markers (Aldob, Bhmt, Gatm). Kept as a section-specific state.
- **AKI c19, "cTAL", 80 units, 94% IR.** Its markers (Havcr1, Spp1, Krt20, Lcn2, Slc6a18) are those of injured PT, so
  these units leave the PT cohort. **AKI c15, Unassigned, 2,425 units, 72% IR** (Cryab, Tnc, Vim, Thbs1), is dropped.
  Under the package's rule the most injured units therefore leave the AKI PT cohort; this is recorded as a limitation
  of every AKI result.
- **cTAL/mTAL near-ties** (human vs mouse c6 and c15, AKI c8) stay within one family and do not touch the PT.

## Notebook 71 protocol: notebook 69's PT pathway analysis, unchanged, on the package pseudospace

Frozen on 2026-10-12, before any fit on the package pseudospace.

**Why.** The user froze Caleb's package `pseudospace_reconstruction` as the method that creates the pseudospace for
both datasets. Its stages are 01–04: bins → units, QC, two-pass Harmony with labels, and an elastic principal curve.
The user asked to "switch all results to use calebs package now". The package replaces two coordinates:
- SCF13, the notebook 13 scFates coordinate on notebook 03's structures (human vs mouse);
- the mouse-only PT DPT of `mouse_only_v5` (AKI vs control).

By the user's decision, notebook 71's results become primary. This protocol does not decide that. It fixes the
analysis before any fit, and it says beforehand which results would count as a material change from notebook 69.

**The single change.** The units, their counts, their PT labels and their position now come from the package export.
Notebook 69's analysis is otherwise unchanged:
- the gene model and the specimen-shape denominator;
- the positional F-test (primary);
- the clustering comparator, run through our pipeline as notebook 69 runs it;
- the swaps and their pass rules;
- the gene sets and redundancy clusters;
- the mt exclusion.

Notebook 71's labels:
- notebook `analysis/notebooks/71_pt_pathways_package_pseudospace.ipynb`;
- logic version `71.package.1`;
- results in `results/pt_pathways_package_pseudospace/` (`tables/`, `figures/`, `stage_cache/`, `run_record.json`).

### 1 · Inputs

**Files.**
- `results/package_pseudospace/{hvm,aki}/pt_units.h5ad` (ws70's export; contract below);
- `results/minimal_pt_scfates/ortholog_map_used.csv`, the map notebook 69 uses;
- `data/mouse_vs_human/pathway_gene_sets/{Reactome_2022,MSigDB_Hallmark_2020,KEGG_2019_Mouse}.json`;
- read only, for the comparison with notebook 69: `results/pt_pathways_without_mito/` (tables and stage cache),
  `results/pt_reconstruction_v2/primary_pt_coordinate.csv` and `results/mouse_only_v5/all_mouse_tubules_scanpy_dpt.h5ad`
  (69's positions, for the unit matching);
- for the secondary sections: the files listed in section 6.

**Export contract,** checked by `package_inputs.read_pt_units`; violations raise:
- PT units only, after the package's QC;
- X holds raw nonnegative integer counts; var index is unique mouse symbols;
- human vs mouse: the shared one-to-one orthologs, renamed to mouse symbols, with a `human_symbol` column that must
  equal the ortholog map's pair;
- obs index is a unique unit id. Columns:
  - `sample`, one of the six specimens;
  - `group`: Mouse/Human or Control/IR, one per specimen;
  - `segment_class`: PT-S1, PT-S2 or PT-S3, the package's labels;
  - `pseudospace`: the finite raw stage-04 value;
  - `polygon_id` and `old_structure_id`;
- `uns['provenance']`, which must not be the ws71 development object.

**Builders** (`pseudospace/package_inputs.py`, tests in `tests/test_package_inputs.py`). They apply notebook 69's
rules to the export. `exclude_genes` is the 13 ortholog-map mouse symbols that start with `mt-`.

- **Shared rules.**
  - Each specimen's [1st, 99th] percentile common support (`common_support_mask`).
  - The position is the stage-04 value min–max scaled over the units inside the support (section 3).
  - The S1 median must lie below the S3 median.
- **`hvm_inputs`** returns `pt_inputs`' dict (notebook 37's construction):
  - library: the total over the measured shared orthologs without the mt genes (the NB offset);
  - `Y`: log1p(count / full library × 10⁴), with mt reads in the library;
  - universe: measured in both inputs, not mt, equal-specimen detection ≥ 2%;
  - pathways: members mapped through the ortholog map, 10–300 tested members;
  - covariates and strata (3 bins), and the 61-point specimen-overlap grid.
- **`aki_inputs`** returns notebook 69's AKI objects (its cell 4, notebook 45's construction):
  - panel: var genes in the ortholog map;
  - measured: genes with any count in the kept units;
  - library: the total over the measured panel without the mt genes; `Y` uses the full library;
  - pool: the human vs mouse universe ∩ measured. Notebook 71 passes `hvm_inputs`' universe with the mt genes kept,
    as notebook 69 passed notebook 12's universe;
  - universe: the pool with equal-specimen detection ≥ 2%, minus mt;
  - pathways: notebook 37's on (human vs mouse universe minus mt), then restricted to the AKI universe with the same
    10–300 rule;
  - covariates: notebook 45's, meaning mean, detection and the share of segments with a detection.

**Checked before any fit (phase 1, loading only).** Fed pt_units-like exports of the *old* inputs, the builders
reproduce notebook 69's inputs exactly.
- Human vs mouse: SCF13 with the tubule_by_gene counts, under `rescale='unit'`. AKI: mouse_only_v5 with
  `pt_subset_scanpy_dpt`.
- Reproduced, with no difference: positions, specimens, segments, the case flag, counts (same digest), library, the
  library with mt, `Y`, covariates, strata, grid, genes (11,059 and 9,095), and the gene sets (1,513 and 1,429; same
  ids, order and members). The mt exclusions are also reproduced: 12 and 13 out of the libraries, 12 and 12 out of
  the universes.
- The human vs mouse builder with mt kept gives notebook 12's universe (11,071 genes). Fed that universe, the AKI
  builder reproduces notebook 69's AKI genes and sets.
- The check ran as an untracked scratch script. The notebook does not repeat it: it needs the old inputs.

**Where the new inputs necessarily differ from notebook 69's.**
- **D1 · Units and counts.** Package units come from stage 01 (2 µm bins aggregated into segmentation polygons) and
  are filtered by stage 02 QC. Notebook 69 used the `tubule_by_gene` structures and mouse_only_v5. Unit numbers,
  depths and counts all differ.
  - On the AKI trial run: 15,013 PT units, 14,048 after the support, against 16,991; median library 5,281.
- **D2 · The PT cohort and the S1/S2/S3 labels** are the package's stage-03 labels. They replace notebook 03/13
  (human vs mouse) and mouse_only_v5 (AKI). In notebook 69 the labels enter only:
  - the cohort;
  - the orientation guard;
  - AKI's coverage covariate, which no call uses.

  The zone cuts come from the marker rule, not from the labels, so that rule is unchanged.
- **D3 · The coordinate** is the stage-04 elastic principal curve, fitted on the whole nephron and restricted to PT,
  then min–max scaled over the common support. Notebook 69's human vs mouse coordinate was SCF13, used unscaled.
- **D4 · The human vs mouse panel** is the export's shared orthologs, all treated as measured in both inputs. A var
  column `measured_in_both`, if present, overrides this. Notebook 69 instead took per-specimen panel membership from
  the source files. The universe changes through detection on the new units.
- **D5 · The AKI universe** is restricted to the *new* human vs mouse universe, not notebook 12's. It therefore moves
  with human vs mouse detection.
  - On the trial run with the old universe, 10,193 AKI genes pass, against 9,095.
- **D6 · The log-normalized matrix** is rebuilt from the export's raw counts by notebook 69's rule. The package's own
  `lognorm` layer is not used. The matrix feeds only the zone cuts and the descriptive covariates.
- **D7 · No stage-key guard against notebooks 66/67.** The equivalence check above replaces notebook 69's KEPT guard.
- **D8 · The package version** of stages 01–04 is ws70's and is recorded from `uns['provenance']`. The comparator
  functions keep commit 59da04b, which is asserted as in notebook 69: `compare_groups`, `run_pathways`,
  `cluster_gene_sets`, `zone_cuts_from_markers` and `build_gene_sets`.

**Guards (not cached).**
- G1: the export contract above.
- G2: in every specimen, the segment medians increase S1 < S2 < S3 along the analysed position.
- G3: the mt genes are exactly the 13 `mt-` symbols of the ortholog map. The guard reports how many leave each library
  (expected 12 human vs mouse, 13 AKI) and each universe (at most 12). Each library falls by exactly their counts.
- G4: df₂ = 12 for the condition split and 6 for each relabeling, from the design ranks. The specimen curves are
  identical across the three partitions, and the grid is `common_grid`.
- G5: the zone cuts are increasing. Clean marker runs are reported; if a run is not clean, the notebook says so and
  proceeds.
- G6: the package's gene-set builder and ours agree on the tested sets and members (notebook 69's guard).
- G7: the clustering ORA recomputed from the gene groups reproduces the run's own calls (notebook 69's guard).
- G8, reported only: whether the AKI universe admits Havcr1, Lcn2, Mki67 and Top2a.
  - On the trial run all four are in the package AKI matrix and in the one-to-one ortholog map.
  - Their equal-specimen AKI detection is 46%, 25%, 8.7% and 8.7%, so the AKI detection rule alone would admit them.
  - They are excluded only because they are absent from notebook 12's human vs mouse universe. They are not detected
    in ≥ 2% of the healthy units.
  - Notebook 71 reports whether the *new* human vs mouse universe admits them. The rule is not changed.

### 2 · Inherited parameters (all from notebook 69, unchanged)

| Item | Value |
|---|---|
| mt genes | ortholog-map mouse symbols with prefix `mt-` (13); out of library and universe; `Y` keeps the full library |
| Detection rule | equal-specimen mean detection ≥ 0.02 |
| Gene sets | Reactome 2022, MSigDB Hallmark 2020, KEGG 2019 Mouse; 10–300 tested members (≤ universe − 1) |
| Redundancy clusters | package `cluster_gene_sets`: 1 − Jaccard, average linkage, cut at Jaccard 0.5, medoid names |
| Partitions | `partition_groups`: condition split plus two balanced relabelings (species as nuisance shape) |
| Gene model | `standard_gene_model`: NB GLM on raw counts, offset log(library without mt), cubic B-spline with 6 df (internal knots at position quartiles, boundary knots 0 and 1), prior weights n / (2·J·n_j), ML α per gene under the condition split's M_full, reused for the relabelings; non-converged filled by median, separated flagged |
| Grid | `common_grid`: 101 points over the intersection of the specimens' [1st, 99th] percentile ranges |
| Denominator | `specimen_shape_lr`: LR_spec = D(M_full) − D(M_spec); df₂ 12 / 6 |
| Positional F-test (primary) | s² = LR_spec / df₂; `squeezeVar` with covariate log mean raw count per unit, robust FALSE (notebook 69's call); F = (T_spatial / 6) / s²_post on (6, df₂ + d₀); z = Φ⁻¹(1 − p); gene q = BH over all genes |
| Set test (primary) | `cameraPR` on z, inter.gene.cor 0.01, use.ranks FALSE; call = Up and FDR ≤ 0.05; one-sided p = p/2 if Up else 1 − p/2 |
| Clustering comparator | package stages 06–07 (`compare_groups`, `run_pathways`) on our specimen curves (`clustering_first_inputs`); curves floored at 0.5 count over each specimen's summed library (mt-free); `CompareParams(prior_k=0)`, other defaults (window_half_width 0.02, peak_window_half_width 0.02); `PathwayParams(min_genes=10, max_genes=300)`; node spacing 1/29; zone cuts from `zone_cuts_from_markers` (Slc5a12, Slc13a3, Slc16a9) on the reference specimens' `Y` (control mice in both datasets); runs: condition × {default null, each single-split null}, each relabeling × the fair null (the other relabeling); a call is matched ORA q < 0.05 in any gene group or cluster |
| Pass rules | relabeling: NCDR < 0.25 and each relabeling calls < 5% of tested sets; calibration: pooled share of relabeled one-sided p ≤ 0.05 is ≤ 0.075 (clustering: list-based convention, best matched p × number of lists) |
| Proteasome sensitivity | the 8 subunits Psma3, Psma4, Psmb5, Psmb6, Psmc2, Psmc6, Psmd11, Psmd14 removed: cameraPR on the same z without them (no squeezeVar refit, sets < 10 dropped); clustering ORA rerun with them outside the universe on the condition run's groups and clusters |
| Drift statistic | d_g = mean Δ over the last 21 of 101 grid points minus the first 21; converged, non-separated genes |
| Parallel workers | `PSEUDOSPACE_N_JOBS` (never in a cache key) |

### 3 · The position scale

**The question.** Notebook 69's positions lay in [0, 1]. Is every downstream step invariant to an affine change of
position after trimming?

**Invariant, or equivariant:**
- **The NB gene model and its denominator.** `make_gam_design` fixes the boundary knots at 0 and 1, but the internal
  knots are position quartiles, which are equivariant. The basis therefore spans the same cubic-spline space on the
  data for any affine map that keeps the positions in [0, 1], and `nb_gam` fits on an orthonormal basis of that space.
  - Synthetic check (untracked scratch script): positions in [0, 1] against the same positions compressed to
    [0.02, 0.18].
  - T_level, T_spatial, T_total, α, LR_spec and Δ agree to 10⁻¹¹. The convergence and separation flags are identical.
  - So the F, z, q, cameraPR calls and d_g are unaffected by the scale.
- The 1–99% trimming and `common_grid`, which are quantile-based.
- The drift windows, which are grid indices.
- The zone cuts: 120 equal-width bins between the trimmed reference units' extremes.
- The overlap grid weights: σ in grid bins.
- The gene-covariate coverage bins: equal width over the range.

**Not invariant: the clustering comparator.** Three of its parameters are in axis units:
- `CompareParams.window_half_width` = 0.02, the local-window selection rule;
- `peak_window_half_width` = 0.02, the pair peak search;
- the node spacing of 1/29: `locate` calls a location a border within 1/58 of a cut.

These change the selected genes, the shape-class locations, the group names, and so the clustering calls. On a
compressed scale the windows cover a larger share of the PT.

**The rule.** Min–max over the units inside the common support, giving [0, 1]. This is notebook 69's AKI rule. It
keeps those three parameters as fixed fractions of the trimmed PT. Rescaling does not change the trim, because the
mask is affine-invariant.

In grid-span terms, the 0.02 window was:
- 2.6% in notebook 69's human vs mouse run (grid 0.053–0.836);
- 2.0% in its AKI run (grid 0–1);
- 2.3% on the AKI trial run under the new rule (grid 0.053–0.939).

**The alternatives would compress the PT on the trial run.**
- **`pt_inputs`' rule** (rescale only outside [0, 1]) would leave the raw values as they are. The kept PT spans only
  0.031–0.249, so the windows would cover about 9% of the PT.
- **Min–max over all PT units:** 58 PT units lie beyond 0.35 raw, up to 0.87. The kept core would map to
  [0.035, 0.287], a 3.5-fold compression.

**What the rule changes, beyond the coordinate itself:** only the clustering comparator's numbers. These are its
selected genes, groups, clusters and calls, the NCDR and calibration values of its rows, and the overlap and Jaccard
with the normalized method.
- Positional F-test results do not depend on this choice.
- The default `CompareParams` and 1/29 are inherited unchanged. The package curve's own node spacing (0.018 raw over
  the whole nephron on the trial run, about 0.08 on the rescaled PT) is not substituted.

### 4 · Primary analysis and outputs (notebook 69's sections, unchanged)

Sections 1–5, 7 and 9 of notebook 69 run as they are. The run keeps 69's table and figure names and its stage names
(`gene_model_{ds}_{label}`, `specimen_shape_{ds}_{label}`, `method_b_{ds}_{contrast}__null_{null}`), so readers that
glob them still work.

**Tables, as in notebook 69:**
- `inputs.csv`;
- `redundancy_clusters_{ds}.csv`;
- `gene_model_summary.csv`;
- `gene_statistics_{ds}.csv`, with the condition split and both relabelings;
- `pathway_tests_{ds}_cameraPR_{label}.csv`;
- `moderation.csv`;
- `clustering_runs.csv`;
- `method_table.csv`, with notebook 69's rows marked "69 ·" as the before (replacing 66/67's);
- `check_d_overlap.csv`, with `run` = "notebook 71" and "notebook 69";
- `posthoc_clustering_lists.csv`, now pre-specified: the list that carries most calls, where the 8 subunits sit, and
  the calls through "Down · gradient toward S3";
- `run_record.json`.

**Figure 1** keeps notebook 69's layout:
- a–b: the median Δ of the null genes along the PT, for notebook 71 with notebook 69 dashed;
- c–d: the 2 × 2 overlap of the methods, with notebook 69's counts in brackets.

**Changed meaning.** Notebook 69's checks (a), (b) and (e) were before/after tests of the mt change. In notebook 71 the
"before" is notebook 69, so they become the descriptive comparison of section 5. Check (a)'s null-drift statistic is
also reported on notebook 71's own q > 0.5 genes as a diagnostic; see M7.

**Level-dominance of the clustering method** (new; computed from the condition run's gene table and ORA):
- **L1:** the share of selected genes without a replicated shape (`shape_class` "uniform").
  - Notebook 69: human vs mouse 6,100 / 8,659 = 0.704; AKI 2,864 / 5,217 = 0.549.
- **L2:** the share of clustering calls that are significant only through level-only lists, meaning the "· uniform"
  gene groups or clusters U0/D0.
  - Notebook 69: 11 / 132 = 0.083 and 196 / 455 = 0.431.
- **The package's own `pathway_call` classes** among the calls are tabulated beside them.

### 5 · Comparison with notebook 69 (descriptive; tables `compare69_*.csv`)

**Units: the old → new matching (exported for later notebooks).** `unit_matching_{ds}.csv` holds one row per
package unit:
- `unit_id`, `sample`, `polygon_id`, `old_structure_id`;
- `in_support_71`, `position_71`, `pseudospace_raw`, `segment_71`;
- `old_position_69` (notebook 69's analysed position, NaN when the old structure is outside 69's trimmed set) and
  `old_segment_69`;
- `in_set_69`.

Old ids map to:
- human vs mouse: `primary_pt_coordinate.csv` `structure_id` (positions rescaled exactly as `pt_inputs` did);
- AKI: mouse_only_v5 obs names (positions min–max scaled over 69's support).

Reported per specimen:
- the share matched;
- the Spearman of position_71 against old_position_69 on units in both sets;
- the 3 × 3 agreement of the segment labels;
- the units in only one set.

**Genes** (genes tested in both runs, condition split and each relabeling):
- the Spearman of z, T_spatial and d_g;
- the Jaccard of the gene q ≤ 0.05 sets;
- the counts of genes tested only in 71 or only in 69.

**Pathways** (sets tested in both runs; sets tested in only one run are listed separately). For each method:
- kept, lost and gained calls, each with FDR or best matched q, before and after, and the median d_g of its members,
  before and after (notebook 69's check (e) layout);
- the Jaccard over pathways and over redundancy clusters, using 71's clusters, with 69's clusters as a check;
- the Spearman of cameraPR p.

**Pass rules and level-dominance, 69 → 71,** side by side for both methods and datasets.

### 6 · Secondary sections, pre-specified (the comparison page's unfrozen analyses)

Each section is defined exactly as the page's post hoc script defined it. The definitions below are complete; the cited
scripts (`SC/...`) are untracked session scratch, named for provenance only. They
are all secondary: the primary stays the positional F-test with cameraPR, and level stays descriptive.

Inputs are notebook 71's own objects throughout: inputs, the α of the condition split, the GeneModels, the moderated
F, `log_mean_count`, and the clustering runs. Each section saves the tables the page reads, so the page can be rebuilt
from notebook 71's tables alone.

**Shared pass rule for every swap arm** (`SC/li2022/tlevel/run_ds.py:95-97`, `ssgsea.py:74-77`). It needs:
- real calls > 0;
- NCDR < 0.25;
- each swap calls < 5% of the tested sets;
- a calibration share ≤ 0.075, where the share is the fraction of p ≤ 0.05 pooled over both swaps.

The swaps are the two `partition_groups` relabelings.

**S1 · Level and total gene z** (`SC/li2022/tlevel/run_ds.py`, which uses the same code for both datasets).
- **Offset:** log(library) + log(sf_specimen).
  - The specimen pseudobulk rate is summed counts over summed library.
  - Only genes with a pseudobulk > 0 in all 4 specimens are used.
  - sf = exp(median_g[log rate_sg − mean_s log rate_sg]).
- **Designs:** `contrast_designs(position, group, specimen, knots, nuisance_shape)`. `nospec` is M_full without the
  within-group specimen contrasts. A relabeling adds the 0/1 species column.
- **α:** notebook 71's condition-split α, reused.
- **`nested_nb_lr`** gives:
  - T_level = base → level;
  - T_total = base → full;
  - LR_spec_level = nospec → full, with df 2 for the condition split and 1 for a relabeling.
- **Sign:** the group coefficient of an M_level `nb_fit` with the same offset.
- **Level z:** `moderated_f(T_level, 1, LR_spec_level, df_lvl, covariate=log_mean_count)`. The signed z is
  sign · Φ⁻¹(1 − p/2).
- **Total z:** `moderated_f(T_total, 7, LR_spec_level + s²·df₂, df_lvl + df₂, covariate=log_mean_count)`. s² is
  notebook 71's own LR_spec / df₂. The z is the upper tail.
- **Set tests** on the level, total and shape z:
  - **cameraPR (0.01).** Level calls are FDR ≤ 0.05 in either direction, with the two-sided p as calibration p. Shape
    and total calls are Up and FDR ≤ 0.05, with `p_one_sided` as calibration p.
  - **GSEA prerank** (`rank_first_gsea`, raw metric, 1,000 permutations, seed 12, sizes 10–300). Level calls are
    q ≤ 0.05. Shape and total calls are NES > 0 and q ≤ 0.05. The direction is the sign of NES.
- **Saved:** gene z per arm and partition; per set p, FDR/q, direction and call; `pass_table.csv`.
  - FDR and direction are saved explicitly, because the page reran cameraPR to recover them.

**S2 · ssGSEA + limma** (`SC/li2022/tlevel/ssgsea.py`).
- **Bins:** 6, at the quantiles 1/6…5/6 of the pooled position (`searchsorted`, right).
- **Pseudobulk:** specimen × bin log2(Σcounts / Σlibrary · 10⁶ + 1).
- **Scoring:** `gseapy.ssgsea(sample_norm_method='rank', min_size=10, max_size=300)` with weight 0.25. The score is
  the NES.
- **Level test:** limma `eBayes` on the specimen-mean score, design [1, group (+ nuisance)], coefficient group.
- **Shape test:**
  - design: [specimen indicators, bin dummies 1–5, group × bin, (+ nuisance × bin)];
  - the moderated F on the group × bin columns;
  - a call is BH ≤ 0.05.
- **Saved:** `ssgsea_pass_table.csv` and the per-set results.

**S3 · Level gene list + ORA** (`SC/li2022/level/level_variants.py:56-83`; human vs mouse only, as before).
- **L, dh and dm** follow notebook 68's definitions, applied to notebook 71's condition-split GeneModel:
  - L is `species_level.grid_centred_level(delta, reference=converged & ~separated)`;
  - dh and dm are `pair_levels` on the same-species pairs.

  Notebook 68 had mt reads in the offset; notebook 71 does not. The median-gene centring removes the constant part.
- **Kept genes:** converged, not separated, with dh and dm present, and present in notebook 68's atlas coefficients.
  The atlas coefficients are coordinate-free and read from `results/pt_species_level_control/tables/gene_level.csv`,
  columns `A` and `A_census_human`. Universe genes absent from that table are dropped and counted.
- **Effect and variance:**
  - species split: eff = L, s² = (dh² + dm²)/4;
  - relabelings: as in the script.
- **Moderation:** `squeezeVar(s², 2, covariate = log mean count)`. t = eff / √s²_post on 2 + d₀ df, then BH.
- **Lists:** Up and Down at q ≤ 0.05.
- **ORA:** hypergeometric. The background is the kept genes, BH runs over set × direction, and the best row per set
  is kept.

**S4 · p-value histograms** (`SC/pages/work_cmp/build_compare.py:276-287`). The cameraPR calibration p of each arm,
for the real split and each swap, in 20 equal bins on [0, 1], with the share ≤ 0.05.

**S5 · Implied inter-gene correlation** (`build_compare.py:289-300`). Per arm and partition:
- μ and σ² (ddof 1) are taken over all genes;
- for each set, y_s = n_s(z̄_s − μ)²/σ² − 1;
- ρ̂ = Σ(n_s − 1)y_s / Σ(n_s − 1)².

**S6 · Overlap of every arm with the clustering calls** (`build_compare.py:333-369`).
- Counts: both, arm only, clustering only and neither.
- The Jaccard over pathways and over redundancy clusters.
- **Direction agreement** for the signed arms. It uses clustering's best list (minimum q_matched, then p_matched,
  then larger overlap). The direction is read from the list name with the corrected rule ("Up"/"U…" is up,
  "Down"/"D…" is down). It is not read with `run_ds.py:110-114`'s rule.

**S7 · Atlas positional replication by category** (`SC/li2022/population/01_build.py`, `02_population_tests.py`,
`04_our_calls.py`).
- **Atlas pseudobulks:**
  - Census human and mouse: `data/external/census_pt_segments/*`, with the script's early/late and stage rules;
  - Lake: the notebook 37 cache, checked against the h5ad.

  They are restricted to notebook 71's human vs mouse universe. The atlas statistics are coordinate-free but depend on
  the universe, so they are recomputed.
- **Gene test:** the voom paired late × human test. Unsigned gene z = Φ⁻¹(1 − p), without Ugt3a1/2.
- **Pathway test:** cameraPR (0.01) on sets of ≥ 5.
- **Categories** come from notebook 71's calls. Per category, the table gives tested, called, share and median z,
  with a one-sided Mann-Whitney test against "neither".
- **Gene agreement:** among notebook 71's q ≤ 0.05 genes, the share with voom BH ≤ 0.05 and the same sign (sign =
  71's d_g). This is compared with uncalled genes reweighted over 5 × 5 strata of log mean count × atlas level.

**S8 · Atlas level replication** (`SC/li2022/level/atlas_level_camera.py`).
- **Per-donor level:** 2/3 · log2CPM(early) + 1/3 · log2CPM(late).
- **Donors:** Lake 7 or Census human 6, against the 6 Census males aged 3–6 or 12 months.
- **Fit:** limma `~species`, trend and robust.
- **Set tests:** camera with inter.gene.cor NA, camera with 0.01, and ORA by direction.
- **Groups:** notebook 71's call categories plus S3's ORA calls.
- **Null:** 40 random 6 v 6 splits of the Census males (`default_rng(11)`), now saved as a table.

**S9 · Winsorized sensitivity** (`SC/li2022/audit/winsor.py`, `compare_b.py`).
- **Cap:** y' = min(y, μ + c·√(μ + αμ²)) from notebook 71's condition M_full and α, with c ∈ {8, 8f, 4, 4f}. In "f"
  variants the cap is never below 1 count.
- **Rerun:** the full positional F-test with α re-estimated, using the covariate log of the winsorized mean.
- The c values were ad hoc in the script and are fixed here before the fit. The clustering method is not refit.

**S10 · Per-gene mechanics tables** (`SC/pages/work_mech/reproduce.py`, `build.py:93-160`). Notebook 71 saves them
directly, so no refit is needed afterwards.
- **From the gene model and the F-test:** log mean count, T_spatial/6, LR_spec/12, s²_post, F, z, q and fit_ok.
- **level_lfc:** the group coefficient of an all-gene M_level fit, divided by ln 2.
- **level_z** from S1, and the size-factor shift.
- **From `compare_groups`:**
  - selected, direction, shape_class, location, cluster;
  - rms, level, shape_rms;
  - null_rms, min_pair_reproduction, shape_null_ratio, min_pair_shape_reproduction;
  - the expression stratum and the selection, shape and local-window floors;
  - the local-window statistics.
- **Curves:** D, P, N, S and the weights; `specimen_filled` and the specimen curves on every second grid point; the
  typical gap; the 20-bin observed rates and segment shares; and the cluster mean with p10/p90 of S.

**Not carried over.**
- The nb34 external-segment table (`pt_external_segment_validation/pathway_level_interactions.csv`) is on an older
  statistic and an older coordinate. The page drops it, or a separate protocol reruns it.
- Notebook 66's raw-`T_spatial` GSEA, the joint test and the set-curve calls stay unrun, as in notebook 69.

### 7 · What would be reported as a material change from notebook 69

Every check runs and is reported plainly; nothing is tuned. Any of the following is stated as a material change in
the notebook summary and in the ledger:
- **M1. The positional F-test fails its pass rule on human vs mouse.** That means no calls, NCDR ≥ 0.25, a
  relabeling calling ≥ 5% of the tested sets, or a calibration share > 0.075.
  - Notebook 69: 97 calls; 0 and 0 relabeled; calibration 0.044.
- **M2. The AKI verdict of the positional F-test changes.** Notebook 69 passed the relabeling rule (NCDR 0.13) and
  failed calibration (0.101).
- **M3. The clustering method loses fair-null specificity** on either dataset (relabeling rule under the fair null).
- **M4. The clustering method stops being level-dominated:** L1 < 0.5 on human vs mouse (69: 0.704). AKI's L1 and both
  datasets' L2 are reported.
- **M5. The methods start to agree:** Jaccard over redundancy clusters > 0.2 on either dataset (69: 0.096–0.198).
- **M6. The human vs mouse normalized call list moves.** Either of these counts:
  - fewer than half of notebook 69's called redundancy clusters, among sets tested in both runs, are called again;
  - any of 69's five strongest calls is lost: Hallmark Fatty Acid Metabolism; Reactome Transport Of Inorganic
    Cations/Anions And Amino Acids/Oligopeptides; KEGG Valine, leucine and isoleucine degradation; Hallmark
    Xenobiotic Metabolism; KEGG Arginine and proline metabolism.
- **M7. The offset tilt returns:** the median d_g over notebook 71's own q > 0.5 genes (converged, not separated) has
  |value| > 0.05 on human vs mouse. Notebook 69: +0.001.

Everything else is descriptive: the gene z correlations, the unit matching, the proteasome sensitivity, and the
secondary arms.

### 8 · What is fixed now

Everything above is fixed, including the position rule (section 3), the AKI universe rule (restricted to the new
human vs mouse universe) and the winsorizing c values. The two items the draft left open were settled by the
coordinator before any fit:
1. **Old ids.** `old_structure_id` holds notebook 70's match: for human vs mouse, the `structure_id` of
   `primary_pt_coordinate.csv` (the obs index of `cross_species_pt_scfates.h5ad`); for AKI, the obs name in
   `mouse_only_v5/all_mouse_tubules_scanpy_dpt.h5ad`.
2. **The comparison page's example genes and pathways.** Each notebook 69 example is kept only if it still meets its
   role's screen; otherwise the screen's top gene replaces it. The screens' thresholds (for example q < 0.01,
   |level| > 1.5) were set after seeing notebook 69, so the page labels the examples illustrative.


## Notebook 71 results

`analysis/notebooks/71_pt_pathways_package_pseudospace.ipynb` (logic version `71.package.1`) writes to
`results/pt_pathways_package_pseudospace/`. It followed the protocol above. Every guard passed, including the
confirmed pass-1 fingerprints. Its section 13 lists how the protocol was read:
- notebook 69's human vs mouse sets were rebuilt with `pt_inputs`, because a rebuild from its gene list alone mismaps
  symbols;
- the Lake atlas uses notebook 37's cache, so new-universe genes are absent from that pair;
- T_spatial and LR_spec were added to the gene tables.

**Inputs.** Human vs mouse: 11,383 units, 11,055 genes, 1,513 sets. AKI: 14,048 units, 10,198 genes, 1,469 sets.

**Material changes.** None of M1–M7 is met. The new calls are primary for both methods on both datasets.

| Rule | Notebook 69 | Notebook 71 |
|---|---|---|
| M1 · positional F-test, human vs mouse | 97 calls; relabeled 0, 0; calibration 0.044 | 83 calls; relabeled 1, 0; NCDR 0.006; calibration 0.071 · **passes narrowly** |
| M2 · positional F-test, AKI | passes relabeling, fails calibration (NCDR 0.13, 0.101) | 103 calls; relabeled 13, 29; NCDR 0.20; calibration 0.120 · same verdict |
| M3 · clustering, fair null | relabeled 0, 0 and 2, 0 | 0, 0 on both datasets |
| M4 · clustering level-dominance L1 (AKI); L2 | 0.704 (0.549); 0.083 (0.431) | 0.708 (0.563); 0.138 (0.349) |
| M5 · Jaccard between the methods (redundancy clusters) | 0.175, 0.096 | 0.171, 0.121 |
| M6 · 69's called clusters called again; five strongest calls | — | 76 of 90; all five kept |
| M7 · null drift, human vs mouse | +0.001 | +0.006 |

The clustering method calls 94 sets on human vs mouse (69: 132) and 459 on AKI (69: 455). Its list-based
calibration share is 8 of 3,026 on human vs mouse in both runs (a coincidence of counts from different swap runs) and
0.017 on AKI (69: 0.005).

**Notebook 69 → 71** (descriptive).
- Matched units: position Spearman 0.87–0.96 per specimen; segment labels agree for 80–82% of human units.
- Gene z Spearman 0.83 (human vs mouse) and 0.88 (AKI).
- Positional F-test calls are stable. Human vs mouse keeps 82, loses 15 (all with 69 FDR 0.013–0.048) and gains 1
  (Jaccard 0.84); AKI Jaccard 0.72.
- Clustering calls move more: Jaccard 0.42 (human vs mouse) and 0.66 (AKI). The list "Down · gradient toward S3",
  which carried 107 of 69's human vs mouse calls, now carries none. Without the 8 proteasome subunits the clustering
  calls fall from 94 to 43.

**Secondary sections.**
- **Swap rule** (S1): only the primary arm, shape with cameraPR, passes on human vs mouse. Every GSEA arm and every
  level and total arm fails on both datasets; the level arms fail by wide margins (human vs mouse cameraPR: 20 calls,
  relabeled 102 and 149).
- **ssGSEA + limma** (S2) passes the swap rule on both datasets but calls most sets: 866 (shape) and 1,208 (level) of
  1,513 on human vs mouse.
- **Level gene list + ORA** (S3): 114 calls, 0 and 0 under the swaps.
- **Atlas replication** (S7): 53–56% of our significant genes against 20–22% of matched uncalled genes; 0 of 75
  clustering-only sets replicate.
- **Winsorizing** (S9) keeps 79–83 of the 83 human vs mouse calls; no cap changes a verdict.
